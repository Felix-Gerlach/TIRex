"""
TIR-Tuner engine — predicts nucleotide edits that increase translation
initiation at a target start codon, or decrease it at a downstream start
codon while preserving the encoded protein. Pure logic + OSTIR scoring.
"""

from .scoring import OstirScorer
from .candidate import Candidate
from .increase_engine import enumerate_edits, greedy_optimize
from .decrease_engine import suppress_downstream

__all__ = [
    'OstirScorer', 'Candidate',
    'enumerate_edits', 'greedy_optimize', 'suppress_downstream',
]
