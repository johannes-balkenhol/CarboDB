"""
CarboxyPred API Routes
REST API endpoints for predictions
"""

from flask import Blueprint, request, jsonify, current_app
from werkzeug.utils import secure_filename
import uuid
import os
from datetime import datetime

from ..pipeline.prediction import get_predictor, predict_sequence
from ..pipeline.feature_extraction import extract_features
from ..database.models import db, Sequence, Prediction, Feature, BatchJob

# Create blueprint
api = Blueprint('api', __name__, url_prefix='/api')


# =============================================================================
# HEALTH CHECK
# =============================================================================

@api.route('/health', methods=['GET'])
def health_check():
    """API health check"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.utcnow().isoformat(),
        'version': '2.0.0'
    })


# =============================================================================
# SINGLE SEQUENCE PREDICTION
# =============================================================================

@api.route('/predict', methods=['POST'])
def predict_single():
    """Predict for a single sequence
    
    Request body:
        - sequence: str (required) - amino acid sequence
        - uid: str (optional) - sequence identifier
        - save: bool (optional) - save to database
    
    Returns:
        - prediction results
    """
    data = request.get_json()
    
    if not data or 'sequence' not in data:
        return jsonify({'error': 'Missing sequence'}), 400
    
    sequence = data['sequence']
    uid = data.get('uid', f'user_{uuid.uuid4().hex[:8]}')
    save_to_db = data.get('save', False)
    
    try:
        # Get predictor
        model_dir = current_app.config.get('MODEL_DIR', 'ml_models')
        predictor = get_predictor(model_dir)
        
        # Make prediction
        result = predictor.predict(sequence, uid)
        
        # Optionally save to database
        if save_to_db:
            _save_prediction(uid, sequence, result)
        
        return jsonify({
            'success': True,
            'prediction': result.to_dict()
        })
    
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        return jsonify({'error': f'Prediction failed: {str(e)}'}), 500


@api.route('/predict/features', methods=['POST'])
def predict_with_features():
    """Predict and return detailed features
    
    Request body:
        - sequence: str (required)
        - uid: str (optional)
    
    Returns:
        - prediction + feature breakdown
    """
    data = request.get_json()
    
    if not data or 'sequence' not in data:
        return jsonify({'error': 'Missing sequence'}), 400
    
    sequence = data['sequence']
    uid = data.get('uid', 'unknown')
    
    try:
        # Extract features
        features = extract_features(sequence, uid)
        
        # Make prediction
        model_dir = current_app.config.get('MODEL_DIR', 'ml_models')
        predictor = get_predictor(model_dir)
        result = predictor.predict(sequence, uid)
        
        # Group features by type
        feature_groups = {
            'amino_acids': {k: v for k, v in features.items() if k.startswith('aa_')},
            'dipeptides': {k: v for k, v in features.items() if k.startswith('dp_')},
            'physicochemical': {k: v for k, v in features.items() 
                               if k in ['hydrophobic', 'charged', 'polar', 'aromatic', 
                                       'small', 'glycine', 'net_charge']},
            'motifs': {k: v for k, v in features.items() if k.startswith('motif_')},
            'invariant': {k: v for k, v in features.items() if k.startswith('inv_')},
        }
        
        return jsonify({
            'success': True,
            'prediction': result.to_dict(),
            'features': feature_groups,
            'length': features.get('length', len(sequence))
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


# =============================================================================
# BATCH PREDICTION
# =============================================================================

@api.route('/predict/batch', methods=['POST'])
def predict_batch():
    """Submit batch prediction job
    
    Request:
        - file: FASTA file upload
        - email: notification email (optional)
    
    Returns:
        - job_id for status tracking
    """
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
    
    email = request.form.get('email')
    
    try:
        # Generate job ID
        job_id = str(uuid.uuid4())
        
        # Save uploaded file
        filename = secure_filename(file.filename)
        upload_dir = current_app.config.get('UPLOAD_DIR', 'uploads')
        os.makedirs(upload_dir, exist_ok=True)
        filepath = os.path.join(upload_dir, f'{job_id}_{filename}')
        file.save(filepath)
        
        # Count sequences
        from Bio import SeqIO
        sequences = list(SeqIO.parse(filepath, 'fasta'))
        n_sequences = len(sequences)
        
        # Create job record
        job = BatchJob(
            job_id=job_id,
            user_email=email,
            status='pending',
            total_sequences=n_sequences,
            input_file=filepath
        )
        db.session.add(job)
        db.session.commit()
        
        # Queue background task (Celery)
        from ..tasks import process_batch_job
        process_batch_job.delay(job_id)
        
        return jsonify({
            'success': True,
            'job_id': job_id,
            'total_sequences': n_sequences,
            'status': 'pending',
            'message': f'Job submitted. {n_sequences} sequences will be processed.'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@api.route('/jobs/<job_id>', methods=['GET'])
def get_job_status(job_id):
    """Get batch job status"""
    job = BatchJob.query.filter_by(job_id=job_id).first()
    
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    
    return jsonify({
        'job_id': job.job_id,
        'status': job.status,
        'total_sequences': job.total_sequences,
        'processed': job.processed,
        'progress': job.processed / job.total_sequences * 100 if job.total_sequences > 0 else 0,
        'created_at': job.created_at.isoformat() if job.created_at else None,
        'completed_at': job.completed_at.isoformat() if job.completed_at else None,
        'error': job.error_message
    })


@api.route('/jobs/<job_id>/results', methods=['GET'])
def get_job_results(job_id):
    """Download batch job results"""
    job = BatchJob.query.filter_by(job_id=job_id).first()
    
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    
    if job.status != 'completed':
        return jsonify({'error': f'Job not completed. Status: {job.status}'}), 400
    
    if not job.result_file or not os.path.exists(job.result_file):
        return jsonify({'error': 'Results file not found'}), 404
    
    from flask import send_file
    return send_file(job.result_file, as_attachment=True)


# =============================================================================
# DATABASE QUERIES
# =============================================================================

@api.route('/sequences/<uid>', methods=['GET'])
def get_sequence(uid):
    """Get sequence and predictions by UID"""
    seq = Sequence.query.filter_by(uid=uid).first()
    
    if not seq:
        return jsonify({'error': 'Sequence not found'}), 404
    
    predictions = [p.to_dict() for p in seq.predictions]
    domains = [d.to_dict() for d in seq.domains]
    
    return jsonify({
        'sequence': seq.to_dict(),
        'predictions': predictions,
        'domains': domains
    })


@api.route('/sequences/search', methods=['GET'])
def search_sequences():
    """Search sequences by criteria
    
    Query params:
        - organism: organism name filter
        - ec: EC number filter
        - gene: gene name filter
        - min_km: minimum Km value
        - max_km: maximum Km value
        - is_co2: only CO2 enzymes (bool)
        - limit: max results (default 100)
    """
    query = Sequence.query
    
    # Apply filters
    if request.args.get('organism'):
        query = query.filter(Sequence.organism.ilike(f"%{request.args['organism']}%"))
    
    if request.args.get('ec'):
        query = query.filter(Sequence.ec_number.like(f"{request.args['ec']}%"))
    
    if request.args.get('gene'):
        query = query.filter(Sequence.gene_name.ilike(f"%{request.args['gene']}%"))
    
    # Join with predictions for Km filtering
    if request.args.get('min_km') or request.args.get('max_km') or request.args.get('is_co2'):
        query = query.join(Prediction)
        
        if request.args.get('is_co2'):
            query = query.filter(Prediction.is_co2_enzyme == True)
        
        if request.args.get('min_km'):
            query = query.filter(Prediction.km_predicted >= float(request.args['min_km']))
        
        if request.args.get('max_km'):
            query = query.filter(Prediction.km_predicted <= float(request.args['max_km']))
    
    limit = int(request.args.get('limit', 100))
    results = query.limit(limit).all()
    
    return jsonify({
        'count': len(results),
        'sequences': [s.to_dict() for s in results]
    })


@api.route('/stats', methods=['GET'])
def get_stats():
    """Get database statistics"""
    total_sequences = Sequence.query.count()
    total_predictions = Prediction.query.count()
    co2_enzymes = Prediction.query.filter_by(is_co2_enzyme=True).count()
    
    # EC class distribution
    ec_dist = db.session.query(
        Prediction.ec_class, 
        db.func.count(Prediction.id)
    ).filter(
        Prediction.is_co2_enzyme == True
    ).group_by(Prediction.ec_class).all()
    
    return jsonify({
        'total_sequences': total_sequences,
        'total_predictions': total_predictions,
        'co2_enzymes': co2_enzymes,
        'ec_distribution': {ec: count for ec, count in ec_dist if ec}
    })


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _save_prediction(uid, sequence, result):
    """Save sequence and prediction to database"""
    # Check if sequence exists
    seq = Sequence.query.filter_by(uid=uid).first()
    
    if not seq:
        seq = Sequence(
            uid=uid,
            sequence=sequence,
            length=len(sequence),
            source='api'
        )
        db.session.add(seq)
        db.session.flush()
    
    # Add prediction
    pred = Prediction(
        sequence_id=seq.id,
        model_version='v3',
        is_co2_enzyme=result.is_co2_enzyme,
        co2_probability=result.co2_probability,
        ec_class=result.ec_class,
        ec_probability=result.ec_probability,
        km_predicted=result.km_predicted,
        km_log=result.km_log,
        km_confidence=result.km_confidence
    )
    db.session.add(pred)
    db.session.commit()
