"""
RBS library / TIR ramp design.

Given a target start codon, enumerate single + dinucleotide edits in the
upstream window, then pick a spread of variants whose predicted TIRs are spaced
(log-scale) across the achievable range. The result is an expression "ramp" —
a panel of constructs covering weak → strong initiation — handy for tuning
experiments.

Returns a list of Candidate (one per ramp rung, ascending TIR). The wild-type
is included as the natural anchor.
"""

import math
from typing import Callable, List, Optional

from .candidate import Candidate
from .increase_engine import enumerate_edits
from .scoring import OstirScorer
from .constraints import ConstraintConfig


def design_library(seq: str, target_start: int, scorer: OstirScorer,
                   n_variants: int = 8, upstream_len: int = 20,
                   include_dinucleotide: bool = True,
                   cfg: Optional[ConstraintConfig] = None,
                   progress: Optional[Callable[[int, int, str], None]] = None
                   ) -> List[Candidate]:
    cfg = cfg or ConstraintConfig()
    cands = enumerate_edits(seq, target_start, scorer, upstream_len,
                            include_dinucleotide, cfg=cfg, progress=progress)
    baseline = scorer.score_start(seq, target_start)

    # Candidate pool = wild-type + all single/dinucleotide edits.
    wt = Candidate(edit_type='wild-type', position=target_start,
                   change='(none)', new_sequence=seq, tir=baseline,
                   baseline_tir=baseline, target_start=target_start,
                   notes='unmodified reference')
    pool = [wt] + cands
    pool = [c for c in pool if c.tir > 0]
    if not pool:
        return [wt]

    lo = min(c.tir for c in pool)
    hi = max(c.tir for c in pool)
    if hi <= lo:
        return [wt]

    # Log-spaced TIR targets across [lo, hi].
    log_lo, log_hi = math.log10(lo), math.log10(hi)
    targets = [10 ** (log_lo + (log_hi - log_lo) * i / (n_variants - 1))
               for i in range(max(2, n_variants))]

    chosen = []
    used = set()
    for t in targets:
        # nearest unused candidate (by log distance) to this rung
        best = None
        for idx, c in enumerate(pool):
            if idx in used:
                continue
            d = abs(math.log10(c.tir) - math.log10(t))
            if best is None or d < best[0]:
                best = (d, idx)
        if best is not None:
            used.add(best[1])
            chosen.append(pool[best[1]])

    # De-dup by sequence, sort ascending TIR, annotate rung.
    seen = set()
    ramp = []
    for c in sorted(chosen, key=lambda c: c.tir):
        if c.new_sequence in seen:
            continue
        seen.add(c.new_sequence)
        ramp.append(c)
    for i, c in enumerate(ramp, 1):
        fold = (c.tir / baseline) if baseline else float('inf')
        c.notes = (f'rung {i}/{len(ramp)} · {fold:.2f}× WT'
                   + (f' · {c.notes}' if c.notes and c.edit_type == 'wild-type'
                      else ''))
    return ramp
