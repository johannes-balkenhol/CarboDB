"""
startup.py — Load all models ONCE at startup, keep in memory.
CRITICAL: Never reload models per request — costs 30s overhead each time.
"""
import os
import json
import logging
import numpy as np
import xgboost as xgb

log = logging.getLogger(__name__)

class ModelStore:
    """Singleton holding all loaded models and metadata."""
    ready = False
    
    # XGBoost models
    xgb_binary = None
    xgb_ec = None
    xgb_km = None
    
    # ESM-2
    esm_model = None
    esm_alphabet = None
    esm_device = "cpu"
    
    # Feature metadata
    feature_names = None       # list of 1793 feature names in order
    ec_label_map = None        # {'4.1.1.39': 8, ...}
    ec_inv_map = None          # {8: '4.1.1.39', ...}
    pfam_features = None       # set of pfam feature names
    n_features = 1793


EC_NAMES = {
    '1.1.1.39': 'Malate dehydrogenase',
    '1.2.7.7': 'Pyruvate:ferredoxin oxidoreductase',
    '4.1.1.112': 'beta-Carbonic anhydrase',
    '4.1.1.21': 'Phosphogluconate decarboxylase',
    '4.1.1.31': 'Phosphoenolpyruvate carboxylase (PEPC)',
    '4.1.1.32': 'Phosphoenolpyruvate carboxykinase',
    '4.1.1.38': 'Pyruvate carboxylase (ATP)',
    '4.1.1.39': 'Ribulose bisphosphate carboxylase (RuBisCO)',
    '4.1.1.49': 'Phosphoenolpyruvate carboxylase',
    '4.1.1.61': 'Methylmalonyl-CoA decarboxylase',
    '4.1.1.90': 'Dimethylmalonyl-CoA decarboxylase',
    '4.1.1.93': 'Phosphoenolpyruvate carboxykinase (GTP)',
    '4.2.1.1': 'Carbonic anhydrase',
    '6.3.3.3': 'Phosphoribosylformylglycinamidine synthase (PFAS)',
    '6.3.4.14': 'Biotin carboxylase (ACC)',
    '6.3.4.16': 'Carbamoyl-phosphate synthase (small subunit)',
    '6.3.4.18': 'Phosphoribosylamine-glycine ligase',
    '6.3.4.6': 'Urea carboxylase',
    '6.3.5.5': 'Carbamoyl-phosphate synthase (large subunit)',
    '6.4.1.1': 'Pyruvate carboxylase',
    '6.4.1.2': 'Acetyl-CoA carboxylase',
    '6.4.1.3': 'Propionyl-CoA carboxylase',
    '6.4.1.4': 'Methylcrotonyl-CoA carboxylase',
    '6.4.1.5': '3-Methylglutaconyl-CoA carboxylase',
    '6.4.1.6': 'Geranoyl-CoA carboxylase',
    '6.4.1.8': 'Acetopropionyl-CoA carboxylase',
}

# EC classes that have trained Km models
KM_EC_CLASSES = {
    '4.1.1.39', '4.2.1.1', '4.1.1.49', '4.1.1.31',
    '6.3.4.14', '6.4.1.1', '6.4.1.2', '6.4.1.3', '6.4.1.4',
    '6.3.5.5', '6.3.3.3',
}


def load_all_models():
    """Load all models into ModelStore. Call once at startup."""
    models_dir = os.environ.get("MODELS_DIR", "models")
    
    # Load feature names
    feat_path = os.path.join(models_dir, "feature_names_binary.json")
    ModelStore.feature_names = json.load(open(feat_path))
    ModelStore.n_features = len(ModelStore.feature_names)
    ModelStore.pfam_features = {f for f in ModelStore.feature_names if f.startswith("pfam_")}
    log.info(f"Feature names loaded: {ModelStore.n_features} features")

    # Load EC label map
    ec_path = os.path.join(models_dir, "ec_label_map.json")
    ModelStore.ec_label_map = json.load(open(ec_path))
    ModelStore.ec_inv_map = {v: k for k, v in ModelStore.ec_label_map.items()}
    log.info(f"EC label map: {len(ModelStore.ec_label_map)} classes")

    # Load XGBoost models
    ModelStore.xgb_binary = xgb.Booster()
    ModelStore.xgb_binary.load_model(os.path.join(models_dir, "binary_v5.json"))
    log.info("Binary model loaded")

    ModelStore.xgb_ec = xgb.Booster()
    ModelStore.xgb_ec.load_model(os.path.join(models_dir, "ec_v5.json"))
    log.info("EC model loaded")

    ModelStore.xgb_km = xgb.Booster()
    ModelStore.xgb_km.load_model(os.path.join(models_dir, "km_v5_weighted.json"))
    log.info("Km model loaded")

    # Load ESM-2 (optional — skip if not available)
    esm_device = os.environ.get("ESM2_DEVICE", "cpu")
    try:
        import esm, torch
        log.info(f"Loading ESM-2 on {esm_device}...")
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model = model.to(esm_device).eval()
        ModelStore.esm_model = model
        ModelStore.esm_alphabet = alphabet
        ModelStore.esm_device = esm_device
        log.info("ESM-2 loaded successfully")
    except Exception as e:
        log.warning(f"ESM-2 not available: {e} — fast/pfam modes only")
        ModelStore.esm_model = None

    ModelStore.ready = True
    log.info("All models ready")
