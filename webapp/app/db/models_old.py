"""
CarboxyPred Database Models
SQLAlchemy ORM models for the prediction database
"""

from datetime import datetime
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Index

db = SQLAlchemy()

class Sequence(db.Model):
    """Protein sequence with metadata"""
    __tablename__ = 'sequences'
    
    id = db.Column(db.Integer, primary_key=True)
    uid = db.Column(db.String(20), unique=True, nullable=False, index=True)
    accession = db.Column(db.String(20))
    sequence = db.Column(db.Text, nullable=False)
    length = db.Column(db.Integer, nullable=False)
    description = db.Column(db.Text)
    organism = db.Column(db.String(255), index=True)
    organism_id = db.Column(db.Integer)
    gene_name = db.Column(db.String(50), index=True)
    ec_number = db.Column(db.String(20), index=True)
    source = db.Column(db.String(50), default='user')
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    updated_at = db.Column(db.DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    predictions = db.relationship('Prediction', backref='sequence', lazy='dynamic', cascade='all, delete-orphan')
    features = db.relationship('Feature', backref='sequence', lazy='dynamic', cascade='all, delete-orphan')
    domains = db.relationship('Domain', backref='sequence', lazy='dynamic', cascade='all, delete-orphan')
    known_enzyme = db.relationship('KnownEnzyme', backref='sequence', uselist=False, cascade='all, delete-orphan')
    
    def __repr__(self):
        return f'<Sequence {self.uid}>'
    
    def to_dict(self):
        return {
            'uid': self.uid,
            'accession': self.accession,
            'length': self.length,
            'organism': self.organism,
            'gene_name': self.gene_name,
            'ec_number': self.ec_number,
            'description': self.description,
            'source': self.source,
        }


class Prediction(db.Model):
    """ML model predictions for a sequence"""
    __tablename__ = 'predictions'
    
    id = db.Column(db.Integer, primary_key=True)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequences.id', ondelete='CASCADE'), nullable=False, index=True)
    model_version = db.Column(db.String(20), nullable=False)
    is_co2_enzyme = db.Column(db.Boolean, index=True)
    co2_probability = db.Column(db.Float)
    ec_class = db.Column(db.String(50))
    ec_probability = db.Column(db.Float)
    km_predicted = db.Column(db.Float, index=True)
    km_log = db.Column(db.Float)
    km_confidence = db.Column(db.String(20))
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def __repr__(self):
        return f'<Prediction {self.sequence_id} v{self.model_version}>'
    
    def to_dict(self):
        return {
            'model_version': self.model_version,
            'is_co2_enzyme': self.is_co2_enzyme,
            'co2_probability': self.co2_probability,
            'ec_class': self.ec_class,
            'ec_probability': self.ec_probability,
            'km_predicted': self.km_predicted,
            'km_confidence': self.km_confidence,
        }


class Feature(db.Model):
    """Calculated sequence features"""
    __tablename__ = 'features'
    
    id = db.Column(db.Integer, primary_key=True)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequences.id', ondelete='CASCADE'), nullable=False, index=True)
    feature_type = db.Column(db.String(50), nullable=False, index=True)  # 'amino_acid', 'dipeptide', 'invariant'
    feature_name = db.Column(db.String(50), nullable=False, index=True)  # 'aa_G', 'dp_GK', etc.
    feature_value = db.Column(db.Float, nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def __repr__(self):
        return f'<Feature {self.feature_name}={self.feature_value}>'


class Domain(db.Model):
    """Protein domain annotations"""
    __tablename__ = 'domains'
    
    id = db.Column(db.Integer, primary_key=True)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequences.id', ondelete='CASCADE'), nullable=False, index=True)
    domain_id = db.Column(db.String(20), nullable=False, index=True)  # PF00016
    domain_name = db.Column(db.String(100))
    source = db.Column(db.String(20), nullable=False, index=True)  # 'pfam', 'prosite', 'interpro'
    start_pos = db.Column(db.Integer)
    end_pos = db.Column(db.Integer)
    e_value = db.Column(db.Float)
    score = db.Column(db.Float)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def __repr__(self):
        return f'<Domain {self.domain_id}>'
    
    def to_dict(self):
        return {
            'domain_id': self.domain_id,
            'domain_name': self.domain_name,
            'source': self.source,
            'start': self.start_pos,
            'end': self.end_pos,
            'e_value': self.e_value,
        }


class BatchJob(db.Model):
    """Batch processing jobs"""
    __tablename__ = 'batch_jobs'
    
    id = db.Column(db.Integer, primary_key=True)
    job_id = db.Column(db.String(36), unique=True, nullable=False, index=True)
    user_email = db.Column(db.String(255))
    status = db.Column(db.String(20), default='pending', index=True)  # pending, running, completed, failed
    total_sequences = db.Column(db.Integer, default=0)
    processed = db.Column(db.Integer, default=0)
    input_file = db.Column(db.String(255))
    result_file = db.Column(db.String(255))
    error_message = db.Column(db.Text)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    started_at = db.Column(db.DateTime)
    completed_at = db.Column(db.DateTime)
    
    def __repr__(self):
        return f'<BatchJob {self.job_id}>'
    
    def to_dict(self):
        return {
            'job_id': self.job_id,
            'status': self.status,
            'total_sequences': self.total_sequences,
            'processed': self.processed,
            'progress': self.processed / self.total_sequences * 100 if self.total_sequences > 0 else 0,
            'created_at': self.created_at.isoformat() if self.created_at else None,
            'completed_at': self.completed_at.isoformat() if self.completed_at else None,
        }


class KnownEnzyme(db.Model):
    """Known CO2-fixing enzymes with experimental data"""
    __tablename__ = 'known_enzymes'
    
    id = db.Column(db.Integer, primary_key=True)
    sequence_id = db.Column(db.Integer, db.ForeignKey('sequences.id', ondelete='CASCADE'), nullable=False, index=True)
    km_experimental = db.Column(db.Float, index=True)
    km_source = db.Column(db.String(255))  # BRENDA, paper DOI
    ec_verified = db.Column(db.String(20), index=True)
    substrate = db.Column(db.String(100))  # CO2, HCO3-
    organism = db.Column(db.String(255))
    reference = db.Column(db.Text)
    notes = db.Column(db.Text)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def __repr__(self):
        return f'<KnownEnzyme {self.sequence_id}>'
    
    def to_dict(self):
        return {
            'km_experimental': self.km_experimental,
            'km_source': self.km_source,
            'ec_verified': self.ec_verified,
            'substrate': self.substrate,
            'reference': self.reference,
        }
