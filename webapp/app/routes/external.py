"""
External annotation proxy.

Fetches UniProt + AlphaFold metadata for a sequence and serves it through
the CarboDB API. Results are cached in SQLite for 30 days. The proxy validates
that the uniprot_id belongs to a sequence in our database before hitting any
external service (so /api/v1/external/AAAA returns 404 immediately without
touching UniProt). Per-IP rate limiting prevents abuse.

Endpoints
---------
GET /api/v1/external/{uniprot_id}
    Returns combined UniProt summary + AlphaFold metadata.
    Lazy: only call when the user opens the "Extended details" section.

GET /api/v1/external/{uniprot_id}/structure
    Proxies the AlphaFold PDB file (same-origin so the browser's NGL viewer
    can fetch it without CORS gymnastics).

Cache table is created on first import of this module if it doesn't exist.
"""
from __future__ import annotations

import json
import os
import sqlite3
import time
from collections import defaultdict, deque
from typing import Optional

import httpx
from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import Response

router = APIRouter(prefix="/external", tags=["external"])

# ─── Configuration ────────────────────────────────────────────────────────────
CACHE_TTL_SECONDS = 30 * 24 * 3600              # 30 days fresh
HTTP_HEADERS = {"User-Agent": "CarboDB/1.0 (https://github.com/johannes-balkenhol/CarboDB)"}
STALE_TTL_SECONDS = 365 * 24 * 3600             # serve stale up to 1 year if upstream is down
HTTP_TIMEOUT_SECONDS = 8.0
DB_PATH_ENV = "DB_PATH"
DEFAULT_DB = "data/primary/carbodb.sqlite"

# Per-IP rate limit: 60 external annotation requests per minute per IP.
# This is generous for a single user opening many panels but cuts off scrapers.
RATE_LIMIT_REQUESTS = 60
RATE_LIMIT_WINDOW_SECONDS = 60
_RATE_LIMIT_BUCKETS: dict[str, deque[float]] = defaultdict(deque)

# In-process semaphore to prevent thundering herd on the same uniprot_id.
# Two simultaneous requests for AAAA shouldn't both call UniProt.
_INFLIGHT: dict[str, float] = {}

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{uid}.json"
ALPHAFOLD_META_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uid}"
ALPHAFOLD_PDB_URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v{v}.pdb"
# AlphaFold model versions seen so far: 4, 5, 6. Try newest first.
ALPHAFOLD_VERSIONS_TO_TRY = [6, 5, 4]


# ─── DB helpers ───────────────────────────────────────────────────────────────
def _db() -> sqlite3.Connection:
    path = os.environ.get(DB_PATH_ENV, DEFAULT_DB)
    if not os.path.exists(path):
        raise HTTPException(503, f"Database not found at {path}")
    conn = sqlite3.connect(path, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA cache_size = -524288")
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def _ensure_cache_table() -> None:
    """Idempotently create the cache table on import."""
    conn = _db()
    try:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS external_annotations_cache (
                uniprot_id  TEXT PRIMARY KEY,
                data_json   TEXT NOT NULL,
                cached_at   INTEGER NOT NULL,
                fetch_status TEXT NOT NULL  -- 'ok' | 'partial' | 'error'
            )
        """)
        conn.execute("""
            CREATE INDEX IF NOT EXISTS idx_extcache_cached_at
                ON external_annotations_cache(cached_at)
        """)
        conn.commit()
    finally:
        conn.close()


_ensure_cache_table()


# ─── Validation: uniprot_id must exist in our sequences table ─────────────────
def _validate_uniprot_id(uniprot_id: str) -> None:
    """Reject IDs that aren't in our sequences table. Prevents the proxy
    from being used as a free UniProt scraping endpoint."""
    if not uniprot_id or len(uniprot_id) > 20:
        raise HTTPException(400, "Invalid uniprot_id format")
    # UniProt accessions match this regex; reject anything else without a DB hit
    import re
    if not re.match(r"^[A-Z][A-Z0-9]{5,9}(-\d+)?$", uniprot_id):
        raise HTTPException(400, "Invalid uniprot_id format")
    conn = _db()
    try:
        row = conn.execute(
            "SELECT 1 FROM sequences WHERE uniprot_id = ? LIMIT 1",
            (uniprot_id,)
        ).fetchone()
        if row is None:
            raise HTTPException(404, f"uniprot_id {uniprot_id} not in CarboDB")
    finally:
        conn.close()


# ─── Rate limiting ────────────────────────────────────────────────────────────
def _check_rate_limit(client_ip: str) -> None:
    """Sliding-window rate limit. Raises 429 if exceeded."""
    now = time.time()
    bucket = _RATE_LIMIT_BUCKETS[client_ip]
    # Drop timestamps older than the window
    while bucket and bucket[0] < now - RATE_LIMIT_WINDOW_SECONDS:
        bucket.popleft()
    if len(bucket) >= RATE_LIMIT_REQUESTS:
        raise HTTPException(
            429,
            f"Rate limit: {RATE_LIMIT_REQUESTS} requests per "
            f"{RATE_LIMIT_WINDOW_SECONDS}s. Try again shortly."
        )
    bucket.append(now)


def _get_client_ip(request: Request) -> str:
    """Extract client IP, respecting X-Forwarded-For if behind a reverse proxy."""
    xff = request.headers.get("x-forwarded-for")
    if xff:
        return xff.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


# ─── Cache layer ──────────────────────────────────────────────────────────────
def _read_cache(uniprot_id: str) -> Optional[dict]:
    """Return cached data dict + cached_at + freshness, or None if missing."""
    conn = _db()
    try:
        row = conn.execute(
            "SELECT data_json, cached_at, fetch_status "
            "FROM external_annotations_cache WHERE uniprot_id = ?",
            (uniprot_id,)
        ).fetchone()
        if row is None:
            return None
        return {
            "data": json.loads(row["data_json"]),
            "cached_at": row["cached_at"],
            "fetch_status": row["fetch_status"],
        }
    finally:
        conn.close()


def _write_cache(uniprot_id: str, data: dict, status: str = "ok") -> None:
    conn = _db()
    try:
        conn.execute(
            "INSERT OR REPLACE INTO external_annotations_cache "
            "(uniprot_id, data_json, cached_at, fetch_status) "
            "VALUES (?, ?, ?, ?)",
            (uniprot_id, json.dumps(data), int(time.time()), status)
        )
        conn.commit()
    finally:
        conn.close()


# ─── External API fetchers ────────────────────────────────────────────────────
def _fetch_uniprot(uniprot_id: str) -> Optional[dict]:
    """Fetch UniProt JSON and reduce to the fields we care about. Returns
    None on hard failure (network, 5xx). Returns {} on 404 (entry doesn't
    exist on UniProt — distinct from network failure)."""
    url = UNIPROT_URL.format(uid=uniprot_id)
    try:
        with httpx.Client(timeout=HTTP_TIMEOUT_SECONDS, headers=HTTP_HEADERS) as client:
            r = client.get(url, headers={"Accept": "application/json"})
        if r.status_code == 404:
            return {}
        r.raise_for_status()
        raw = r.json()
    except (httpx.HTTPError, ValueError):
        return None

    return _summarize_uniprot(raw)


def _summarize_uniprot(raw: dict) -> dict:
    """Reduce the giant UniProt JSON to a slim summary for the frontend."""
    out: dict = {}

    # Names
    pname = raw.get("proteinDescription", {})
    rec_name = pname.get("recommendedName", {}).get("fullName", {}).get("value")
    out["protein_name"] = rec_name or ""
    out["gene_name"] = ""
    genes = raw.get("genes", [])
    if genes:
        gn = genes[0].get("geneName", {}).get("value")
        if gn:
            out["gene_name"] = gn

    # Organism + lineage
    org = raw.get("organism", {})
    out["organism"] = org.get("scientificName", "")
    out["lineage"] = org.get("lineage", []) or []
    out["taxon_id"] = org.get("taxonId")

    # Function (free text)
    function_texts = []
    for c in raw.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            for t in c.get("texts", []):
                v = t.get("value")
                if v:
                    function_texts.append(v)
    out["function_text"] = " ".join(function_texts)

    # Subcellular location
    locations: list[str] = []
    for c in raw.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                loc = sl.get("location", {}).get("value")
                if loc:
                    locations.append(loc)
    out["subcellular_location"] = locations

    # GO terms grouped by aspect
    go: dict[str, list] = {"molecular_function": [], "biological_process": [], "cellular_component": []}
    aspect_map = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}
    for xref in raw.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        go_id = xref.get("id")
        term = ""
        evidence = ""
        for prop in xref.get("properties", []):
            k = prop.get("key", "")
            v = prop.get("value", "")
            if k == "GoTerm":
                # Format is "F:RuBisCO activity"
                if ":" in v:
                    aspect_letter, term = v.split(":", 1)
                else:
                    term = v
            elif k == "GoEvidenceType":
                # Format is "IDA:UniProtKB"
                evidence = v.split(":")[0] if v else ""
        if go_id and term:
            aspect = aspect_map.get((v.split(":", 1)[0] if v and ":" in v else "F"), "molecular_function")
            # Recompute aspect from the actual GoTerm prefix, not from local var
            for prop in xref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    pv = prop.get("value", "")
                    if pv and ":" in pv:
                        aspect = aspect_map.get(pv.split(":", 1)[0], aspect)
            go.setdefault(aspect, []).append({
                "id": go_id,
                "name": term,
                "evidence": evidence,
            })
    out["go_terms"] = go

    # Active sites, binding sites, domains from the "features" array
    active_sites = []
    binding_sites = []
    domains = []
    ptms = []
    for f in raw.get("features", []):
        ftype = f.get("type", "")
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        desc = f.get("description", "")
        if ftype == "Active site" and start:
            active_sites.append({"position": start, "description": desc})
        elif ftype == "Binding site" and start:
            ligand = ""
            for li in f.get("ligand", {}).items() if isinstance(f.get("ligand"), dict) else []:
                if li[0] == "name":
                    ligand = li[1]
            binding_sites.append({
                "position_start": start,
                "position_end": end or start,
                "ligand": ligand,
                "description": desc,
            })
        elif ftype == "Domain" and start and end:
            domains.append({
                "start": start, "end": end, "name": desc, "source": "UniProt"
            })
        elif ftype in ("Modified residue", "Glycosylation", "Disulfide bond") and start:
            ptms.append({
                "position": start, "type": ftype, "description": desc
            })
    out["active_sites"] = active_sites
    out["binding_sites"] = binding_sites
    out["domains_uniprot"] = domains
    out["ptms"] = ptms

    # Cross-references we care about
    cross_refs: dict[str, list[str]] = defaultdict(list)
    interesting = {"PDB", "KEGG", "Reactome", "BRENDA", "InterPro", "Pfam", "PANTHER", "AlphaFoldDB"}
    for xref in raw.get("uniProtKBCrossReferences", []):
        db = xref.get("database")
        if db in interesting:
            xid = xref.get("id")
            if xid:
                cross_refs[db].append(xid)
    # Deduplicate while keeping order
    out["cross_refs"] = {k: list(dict.fromkeys(v)) for k, v in cross_refs.items()}

    # Sequence length & checksum (handy for sanity checks)
    seq = raw.get("sequence", {})
    out["sequence_length"] = seq.get("length")
    out["sequence_checksum"] = seq.get("crc64")

    return out


def _fetch_alphafold_meta(uniprot_id: str) -> dict:
    """Returns {available, mean_plddt, fragment_count, version}.
    Never raises — returns {available: False} on any failure or 404."""
    url = ALPHAFOLD_META_URL.format(uid=uniprot_id)
    try:
        with httpx.Client(timeout=HTTP_TIMEOUT_SECONDS, headers=HTTP_HEADERS) as client:
            r = client.get(url)
        if r.status_code == 404:
            return {"available": False}
        r.raise_for_status()
        entries = r.json()
        if not entries:
            return {"available": False}
        first = entries[0]
        return {
            "available": True,
            "version": first.get("latestVersion", "v4"),
            "mean_plddt": first.get("globalMetricValue"),
            "fragment_count": len(entries),
            "model_id": first.get("entryId"),
        }
    except (httpx.HTTPError, ValueError):
        return {"available": False}


# ─── Main combined fetch with stale-cache fallback ────────────────────────────
def _fetch_combined(uniprot_id: str) -> dict:
    """Fetch UniProt + AlphaFold and combine. Handles partial failures."""
    uniprot_data = _fetch_uniprot(uniprot_id)
    alphafold_data = _fetch_alphafold_meta(uniprot_id)

    status = "ok"
    if uniprot_data is None:
        status = "error"
        uniprot_data = {}
    elif uniprot_data == {}:
        status = "partial"  # entry doesn't exist on UniProt

    return {
        "uniprot_id": uniprot_id,
        "uniprot": uniprot_data,
        "alphafold": alphafold_data,
        "_fetch_status": status,
    }


# ─── Routes ───────────────────────────────────────────────────────────────────
@router.get("/{uniprot_id}")
def get_external_annotation(
    uniprot_id: str,
    request: Request,
    refresh: bool = Query(False, description="Bypass cache and refetch"),
):
    """Return combined UniProt + AlphaFold annotation for a CarboDB sequence.

    Checks cache first (fresh = under 30 days). On cache miss or `refresh=true`,
    fetches upstream and caches the result. If upstream fails and there is
    stale cache available, returns the stale data with `stale: true`.
    """
    _check_rate_limit(_get_client_ip(request))
    _validate_uniprot_id(uniprot_id)

    cached = _read_cache(uniprot_id)
    now = int(time.time())

    # Cache hit, fresh
    if cached and not refresh and (now - cached["cached_at"]) < CACHE_TTL_SECONDS:
        return {
            **cached["data"],
            "from_cache": True,
            "cached_at": cached["cached_at"],
            "stale": False,
        }

    # Try upstream
    fresh = _fetch_combined(uniprot_id)
    if fresh["_fetch_status"] != "error":
        _write_cache(uniprot_id, fresh, status=fresh["_fetch_status"])
        return {
            **fresh,
            "from_cache": False,
            "cached_at": now,
            "stale": False,
        }

    # Upstream failed — fall back to stale cache if available
    if cached and (now - cached["cached_at"]) < STALE_TTL_SECONDS:
        return {
            **cached["data"],
            "from_cache": True,
            "cached_at": cached["cached_at"],
            "stale": True,
            "_warning": "Upstream APIs unavailable, returning cached data.",
        }

    # No cache, upstream down
    raise HTTPException(
        503,
        "External annotation services (UniProt/AlphaFold) are currently "
        "unavailable and no cached data exists for this sequence."
    )


@router.api_route("/{uniprot_id}/structure", methods=["GET", "HEAD"])
def get_alphafold_structure(uniprot_id: str, request: Request):
    """Proxy the AlphaFold PDB file. Same-origin for the browser.

    NOTE: PDB files are not cached on disk. They're typically 100-500 KB and
    AlphaFold's own CDN handles caching; we just proxy through to avoid CORS.
    """
    _check_rate_limit(_get_client_ip(request))
    _validate_uniprot_id(uniprot_id)

    last_error = None
    with httpx.Client(timeout=HTTP_TIMEOUT_SECONDS * 2, headers=HTTP_HEADERS) as client:
        for v in ALPHAFOLD_VERSIONS_TO_TRY:
            url = ALPHAFOLD_PDB_URL_TEMPLATE.format(uid=uniprot_id, v=v)
            try:
                r = client.get(url)
            except httpx.HTTPError as e:
                last_error = str(e)
                continue
            if r.status_code == 200:
                return Response(
                    content=r.content,
                    media_type="chemical/x-pdb",
                    headers={
                        "Content-Disposition": f'inline; filename="{uniprot_id}.pdb"',
                        "Cache-Control": "public, max-age=604800",  # 7 days
                        "X-AlphaFold-Version": f"v{v}",
                    }
                )
            elif r.status_code == 404:
                last_error = f"v{v} not available"
                continue
            else:
                last_error = f"HTTP {r.status_code}"
    raise HTTPException(404, f"AlphaFold structure not available for {uniprot_id} ({last_error})")


@router.get("/_admin/cache_stats")
def cache_stats():
    """Quick diagnostic endpoint for cache state."""
    conn = _db()
    try:
        row = conn.execute(
            "SELECT COUNT(*) AS total, "
            "       SUM(CASE WHEN fetch_status='ok' THEN 1 ELSE 0 END) AS ok_count, "
            "       SUM(CASE WHEN fetch_status='partial' THEN 1 ELSE 0 END) AS partial_count, "
            "       MIN(cached_at) AS oldest, "
            "       MAX(cached_at) AS newest "
            "FROM external_annotations_cache"
        ).fetchone()
        return dict(row)
    finally:
        conn.close()
