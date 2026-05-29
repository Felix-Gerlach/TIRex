"""Candidate edit dataclass shared by both engines and the UI table."""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class Candidate:
    edit_type: str                 # 'substitution' | 'deletion' | 'insertion'
                                   #  | 'dinucleotide' | 'codon swap' | 'RBS synonymous'
    position: int                  # 1-based position of the edit (anchor)
    change: str                    # human-readable, e.g. 'A->G', 'ins T', 'del C', 'AT->GC'
    new_sequence: str              # full mutated sequence
    tir: float                     # predicted TIR at the relevant start codon
    baseline_tir: float            # wild-type TIR at that start codon
    target_start: int              # 1-based start position scored (after indels)
    # decrease-mode extras
    main_tir: Optional[float] = None        # main start TIR after edit
    main_baseline_tir: Optional[float] = None
    notes: str = ''
    history: List[str] = field(default_factory=list)   # for greedy stacks
    warnings: List[str] = field(default_factory=list)  # constraint soft-flags

    @property
    def warning_text(self) -> str:
        return '; '.join(self.warnings)

    @property
    def delta(self) -> float:
        return self.tir - self.baseline_tir

    @property
    def fold(self) -> float:
        if self.baseline_tir <= 0:
            return float('inf') if self.tir > 0 else 1.0
        return self.tir / self.baseline_tir

    @property
    def main_delta(self) -> Optional[float]:
        if self.main_tir is None or self.main_baseline_tir is None:
            return None
        return self.main_tir - self.main_baseline_tir
