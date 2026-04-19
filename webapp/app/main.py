import os, time, logging
from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from .startup import ModelStore, load_all_models
from .routes.predict import router as predict_router
from .routes.browse import router as browse_router
from .routes.batch import router as batch_router

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger(__name__)

@asynccontextmanager
async def lifespan(app: FastAPI):
    log.info("CarboDB webapp starting...")
    t = time.time()
    load_all_models()
    log.info(f"Models loaded in {time.time()-t:.1f}s")
    yield
    log.info("CarboDB webapp shutting down")

app = FastAPI(title="CarboDB v5",
              description="Carboxylase sequence · function · Km database",
              version="5.0.0", lifespan=lifespan)

app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])

# API routes first
app.include_router(predict_router, prefix="/api/v1")
app.include_router(browse_router, prefix="/api/v1")
app.include_router(batch_router, prefix="/api/v1")

@app.get("/api/v1/health")
def health():
    return {"status": "ok", "models_loaded": ModelStore.ready,
            "db_path": os.environ.get("DB_PATH", "not set")}

# Static frontend LAST — catches everything else
static_dir = os.path.join(os.path.dirname(__file__), "..", "static")
if os.path.exists(os.path.join(static_dir, "index.html")):
    app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")
