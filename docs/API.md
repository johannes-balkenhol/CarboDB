# CarboDB Annotation API

**Script:** `scripts/11_annotate_sequence.py`  
**Purpose:** Annotate a single protein sequence or batch of sequences with carboxylase predictions.  
**Runtime:** ~30s with ESM-2 (full accuracy) | ~3s without ESM-2 (fast mode)

---

## Usage

### Command line

```bash
# Single sequence from FASTA file (full pipeline with ESM-2)
python scripts/11_annotate_sequence.py --fasta my_protein.fasta --pretty

# Inline sequence string
python scripts/11_annotate_sequence.py --sequence MSPQTETKASVEFK... --pretty

# Fast mode (no ESM-2, ~3s, EC accuracy reduced)
python scripts/11_annotate_sequence.py --fasta my_protein.fasta --no-esm2

# Specify kingdom for better Km prediction
python scripts/11_annotate_sequence.py --fasta my_protein.fasta --kingdom plant

# Save output to file
python scripts/11_annotate_sequence.py --fasta my_protein.fasta --out result.json --pretty

# Multi-sequence FASTA (returns JSON array)
python scripts/11_annotate_sequence.py --fasta batch.fasta --out results.json
```

### Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `--fasta` | PATH | — | Input FASTA file (single or multi-sequence) |
| `--sequence` | STRING | — | Raw amino acid sequence (alternative to --fasta) |
| `--out` | PATH | stdout | Output JSON file path |
| `--kingdom` | STRING | `bacteria` | Organism kingdom for Km prediction: `bacteria`, `plant`, `archaea`, `fungi` |
| `--no-esm2` | flag | off | Skip ESM-2 embedding (faster, less accurate EC prediction) |
| `--pretty` | flag | off | Pretty-print JSON output |

---

## Output JSON Schema

Single sequence returns an object. Multi-sequence FASTA returns an array of objects.

```json
{
  "cdb_query_id": "sp|P00875|RBL_SPIOL",
  "sequence_length": 475,
  
  "is_carboxylase": true,
  "carboxylase_probability": 1.0,
  "confidence": "high",
  
  "ec_predicted": "4.1.1.39",
  "ec_name": "ribulose-bisphosphate carboxylase (RuBisCO)",
  "ec_probabilities": {
    "4.1.1.39": 1.0,
    "4.1.1.49": 0.0,
    "4.2.1.1": 0.0,
    "6.3.4.14": 0.0,
    "6.3.5.5": 0.0
  },
  
  "km_predicted_mM": 0.0101,
  "km_predicted_log10": -1.9957,
  "km_ec_used": "4.1.1.39",
  
  "features_used": ["composition", "pfam", "interpro", "esm2"],
  "pfam_hits": ["PF00016", "PF02788"],
  
  "warnings": [],
  "runtime_seconds": 28.4
}
```

### Field Descriptions

| Field | Type | Description |
|---|---|---|
| `cdb_query_id` | string | Sequence ID from FASTA header |
| `sequence_length` | int | Length in amino acids (after cleaning invalid chars) |
| `is_carboxylase` | bool | True if predicted CO₂-fixing carboxylase |
| `carboxylase_probability` | float | Model confidence [0–1]; threshold = 0.5 |
| `confidence` | string | `high` (≥0.90), `medium` (0.70–0.90), `low` (0.50–0.70), `non_carboxylase` (<0.50) |
| `ec_predicted` | string | Top predicted EC class (e.g. `4.1.1.39`) |
| `ec_name` | string | Human-readable EC name |
| `ec_probabilities` | object | Top-5 EC classes with probabilities |
| `km_predicted_mM` | float\|null | Predicted CO₂ Km in mM; null if EC not in trainable set |
| `km_predicted_log10` | float\|null | log₁₀(km_predicted_mM) |
| `km_ec_used` | string\|null | EC class used for Km prediction |
| `features_used` | array | Feature groups computed: composition, pfam, interpro, esm2 |
| `pfam_hits` | array | Pfam domain accessions found by HMMER |
| `warnings` | array | Non-fatal issues (e.g. ESM-2 skipped, short sequence) |
| `runtime_seconds` | float | Total annotation time |

### Confidence levels

| Label | co2_prob | Interpretation |
|---|---|---|
| `high` | ≥ 0.90 | Strong prediction — reliable for downstream use |
| `medium` | 0.70–0.90 | Likely carboxylase — verify with domain knowledge |
| `low` | 0.50–0.70 | Borderline — treat as candidate only |
| `non_carboxylase` | < 0.50 | Not predicted as carboxylase |

### EC classes covered (26 classes)

The model was trained on 26 EC classes. Top classes by database frequency:

| EC | Name | n in DB |
|---|---|---|
| 4.1.1.39 | RuBisCO | 155,198 |
| 4.2.1.1 | Carbonic anhydrase | 69,517 |
| 6.3.4.16 | ACC biotin carboxylase | 47,185 |
| 6.3.4.14 | Pyruvate carboxylase | 44,496 |
| 6.3.5.5 | Carbamoyl-phosphate synthase | 33,639 |
| 6.3.4.18 | 3-Methylcrotonyl-CoA carboxylase | 20,984 |
| 4.1.1.49 | PEP carboxylase | 17,742 |
| 4.1.1.31 | PEP carboxykinase | 17,536 |

### Km-trainable EC classes (10 classes)

Km predictions are only available for sequences predicted in one of these 10 EC classes:
`4.2.1.1`, `4.1.1.39`, `4.1.1.31`, `4.1.1.49`, `6.3.4.14`, `4.1.1.32`, `6.4.1.1`, `6.4.1.2`, `6.4.1.3`, `6.4.1.4`

For other EC classes, `km_predicted_mM` will be `null`.

---

## Flask Integration (for CarboDB-App backend)

### Synchronous endpoint (fast mode only)

```python
from flask import Flask, request, jsonify
import subprocess, json, tempfile, os

app = Flask(__name__)

@app.route('/api/annotate', methods=['POST'])
def annotate():
    data = request.get_json()
    sequence = data.get('sequence', '')
    kingdom  = data.get('kingdom', 'bacteria')
    fast     = data.get('fast_mode', False)
    
    if not sequence:
        return jsonify({'error': 'No sequence provided'}), 400
    
    cmd = [
        'python', 'scripts/11_annotate_sequence.py',
        '--sequence', sequence,
        '--kingdom', kingdom,
    ]
    if fast:
        cmd.append('--no-esm2')
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    
    if result.returncode != 0:
        return jsonify({'error': result.stderr}), 500
    
    # Parse JSON from stdout (last line)
    output = result.stdout.strip().split('\n')
    json_line = next((l for l in reversed(output) if l.startswith('{')), None)
    
    return jsonify(json.loads(json_line))
```

### Asynchronous endpoint (recommended for ESM-2 full pipeline)

ESM-2 takes ~30s — use a task queue (Celery/RQ) with polling:

```python
from celery import Celery
import subprocess, json

celery = Celery('carbodb', broker='redis://localhost:6379/0')

@celery.task
def annotate_task(sequence, kingdom='bacteria', fast=False):
    cmd = ['python', 'scripts/11_annotate_sequence.py',
           '--sequence', sequence, '--kingdom', kingdom]
    if fast:
        cmd.append('--no-esm2')
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    output = result.stdout.strip().split('\n')
    json_line = next((l for l in reversed(output) if l.startswith('{')), None)
    return json.loads(json_line)

# Flask routes
@app.route('/api/annotate/submit', methods=['POST'])
def submit():
    data = request.get_json()
    task = annotate_task.delay(data['sequence'], data.get('kingdom','bacteria'))
    return jsonify({'task_id': task.id})

@app.route('/api/annotate/result/<task_id>', methods=['GET'])
def result(task_id):
    task = annotate_task.AsyncResult(task_id)
    if task.state == 'PENDING':
        return jsonify({'status': 'pending'})
    elif task.state == 'SUCCESS':
        return jsonify({'status': 'done', 'result': task.result})
    else:
        return jsonify({'status': 'error', 'error': str(task.info)})
```

---

## Performance

| Mode | Runtime | EC accuracy | Km accuracy | Use case |
|---|---|---|---|---|
| Full (with ESM-2) | ~30s | 99.8% | R²=0.95 | Production annotation |
| Fast (--no-esm2) | ~3s | ~60–70% | Not reliable | Quick screening |

**Note:** ESM-2 downloads ~2.5 GB on first run. Set `TORCH_HOME` to a persistent cache directory:

```bash
export TORCH_HOME=/path/to/esm2_cache
python scripts/11_annotate_sequence.py --fasta seq.fasta
```

---

## Dependencies

```
# conda env: carboxylase (see environment.yml)
xgboost>=1.7
fair-esm              # ESM-2 model
torch>=2.0
pandas
numpy
hmmer                 # hmmscan must be in PATH
```

HMMER database required at `data/dbs/pfam/Pfam-A.hmm` (configured in `config.py`).

---

## Error Handling

| Warning | Meaning | Impact |
|---|---|---|
| `ESM-2 skipped` | --no-esm2 flag used | EC prediction less accurate |
| `ESM-2 failed: ...` | GPU/memory error | EC prediction falls back to composition+Pfam only |
| `Sequence too short (n < 50)` | Input below minimum length | Prediction may be unreliable |
| `Non-standard amino acids removed` | Input contains X, B, Z, etc. | Cleaned before prediction |
| `Km prediction not available for EC X` | EC not in trainable set | km_predicted_mM = null |
| `Pfam HMM not found` | Missing Pfam-A.hmm | Pfam features set to zero |

---

## Example curl (once Flask API is running)

```bash
# Submit annotation job
curl -X POST http://localhost:5000/api/annotate/submit \
  -H "Content-Type: application/json" \
  -d '{"sequence": "MSPQTETKASVEFKAGVK...", "kingdom": "plant"}'

# {"task_id": "abc123"}

# Poll for result
curl http://localhost:5000/api/annotate/result/abc123
```
