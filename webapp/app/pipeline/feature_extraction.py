"""
CarboxyPred Feature Extraction Pipeline
Calculates all features needed for ML predictions
"""

import numpy as np
from typing import Dict, List, Optional
from dataclasses import dataclass

AA_LIST = 'ACDEFGHIKLMNPQRSTVWY'
AA_SET = set(AA_LIST)

@dataclass
class SequenceFeatures:
    """Container for all extracted features"""
    uid: str
    length: int
    amino_acids: Dict[str, float]
    dipeptides: Dict[str, float]
    physicochemical: Dict[str, float]
    invariant: Dict[str, float]
    motifs: Dict[str, float]
    
    def to_dict(self) -> Dict[str, float]:
        """Flatten all features into single dict"""
        result = {}
        result.update({f'aa_{k}': v for k, v in self.amino_acids.items()})
        result.update({f'dp_{k}': v for k, v in self.dipeptides.items()})
        result.update(self.physicochemical)
        result.update(self.invariant)
        result.update(self.motifs)
        result['length'] = self.length
        return result


class FeatureExtractor:
    """Extract sequence features for CO2 enzyme prediction"""
    
    def __init__(self):
        self.aa_list = AA_LIST
        
    def clean_sequence(self, seq: str) -> str:
        """Clean and validate sequence"""
        seq = ''.join(c for c in str(seq).upper() if c in AA_SET)
        return seq
    
    def extract_amino_acid_composition(self, seq: str) -> Dict[str, float]:
        """Calculate amino acid frequencies"""
        n = len(seq)
        return {aa: seq.count(aa) / n for aa in self.aa_list}
    
    def extract_dipeptide_composition(self, seq: str) -> Dict[str, float]:
        """Calculate dipeptide frequencies"""
        n = len(seq) - 1
        if n <= 0:
            return {f'{a1}{a2}': 0.0 for a1 in self.aa_list for a2 in self.aa_list}
        
        counts = {}
        for a1 in self.aa_list:
            for a2 in self.aa_list:
                counts[f'{a1}{a2}'] = 0
        
        for i in range(len(seq) - 1):
            dp = seq[i:i+2]
            if dp in counts:
                counts[dp] += 1
        
        return {k: v / n for k, v in counts.items()}
    
    def extract_physicochemical(self, seq: str) -> Dict[str, float]:
        """Calculate physicochemical properties"""
        n = len(seq)
        
        return {
            'hydrophobic': sum(seq.count(a) for a in 'AILMFVWY') / n,
            'charged': sum(seq.count(a) for a in 'DEKR') / n,
            'polar': sum(seq.count(a) for a in 'STNQ') / n,
            'aromatic': sum(seq.count(a) for a in 'FWY') / n,
            'small': sum(seq.count(a) for a in 'AGST') / n,
            'glycine': seq.count('G') / n,
            'proline': seq.count('P') / n,
            'cysteine': seq.count('C') / n,
            'net_charge': (seq.count('K') + seq.count('R') - seq.count('D') - seq.count('E')) / n,
        }
    
    def extract_invariant_features(self, seq: str) -> Dict[str, float]:
        """Extract catalytic/invariant features"""
        n = len(seq)
        
        # Catalytic core region (middle 50%)
        start, end = n // 4, 3 * n // 4
        core = seq[start:end]
        core_len = len(core) if core else 1
        
        features = {}
        
        # Catalytic residue content in core
        for aa in 'DEHKCST':
            features[f'inv_cat_{aa}'] = core.count(aa) / core_len
        
        # Catalytic residue clustering
        cat_positions = [i for i, a in enumerate(seq) if a in 'DEHKC']
        if len(cat_positions) > 1:
            distances = [cat_positions[i+1] - cat_positions[i] 
                        for i in range(len(cat_positions) - 1)]
            features['inv_cat_mean_dist'] = np.mean(distances)
            features['inv_cat_std_dist'] = np.std(distances)
            features['inv_cat_clustering'] = 1 / (np.mean(distances) + 1)
        else:
            features['inv_cat_mean_dist'] = 0
            features['inv_cat_std_dist'] = 0
            features['inv_cat_clustering'] = 0
        
        return features
    
    def extract_motifs(self, seq: str) -> Dict[str, float]:
        """Extract enzyme-specific motifs"""
        n = len(seq)
        
        return {
            # RuBisCO motifs
            'motif_rubisco_kk': seq.count('KK') / n * 100,
            'motif_rubisco_gk': seq.count('GK') / n * 100,
            'motif_gxg': sum(1 for i in range(len(seq)-2) 
                           if seq[i] == 'G' and seq[i+2] == 'G') / n * 100,
            # CA motifs
            'motif_ca_hh': seq.count('HH') / n * 100,
            'motif_ca_his_cluster': sum(1 for i in range(len(seq)-9) 
                                       if seq[i:i+10].count('H') >= 3) / n * 100,
            # PEPC motifs
            'motif_pepc_rr': seq.count('RR') / n * 100,
            # Biotin carboxylase
            'motif_biotin_amk': seq.count('AMK') / n * 100,
            # Key dipeptides for Km
            'motif_dp_nc': seq.count('NC') / (n-1) * 100 if n > 1 else 0,
            'motif_dp_he': seq.count('HE') / (n-1) * 100 if n > 1 else 0,
            'motif_dp_fd': seq.count('FD') / (n-1) * 100 if n > 1 else 0,
        }
    
    def extract_all(self, seq: str, uid: str = 'unknown') -> SequenceFeatures:
        """Extract all features from a sequence"""
        seq = self.clean_sequence(seq)
        
        if len(seq) < 50:
            raise ValueError(f"Sequence too short: {len(seq)} aa (minimum 50)")
        
        return SequenceFeatures(
            uid=uid,
            length=len(seq),
            amino_acids=self.extract_amino_acid_composition(seq),
            dipeptides=self.extract_dipeptide_composition(seq),
            physicochemical=self.extract_physicochemical(seq),
            invariant=self.extract_invariant_features(seq),
            motifs=self.extract_motifs(seq),
        )
    
    def extract_for_model(self, seq: str, feature_cols: List[str]) -> np.ndarray:
        """Extract features in the order expected by ML model"""
        features = self.extract_all(seq)
        feature_dict = features.to_dict()
        
        # Return features in correct order
        return np.array([feature_dict.get(col, 0.0) for col in feature_cols])


# Convenience function
def extract_features(seq: str, uid: str = 'unknown') -> Dict[str, float]:
    """Quick feature extraction"""
    extractor = FeatureExtractor()
    return extractor.extract_all(seq, uid).to_dict()
