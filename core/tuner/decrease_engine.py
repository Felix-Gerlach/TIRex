"""
Decrease mode: suppress translation initiation at a *downstream* internal
start codon WITHOUT changing the main protein's amino-acid sequence.

Frame is anchored at the main start codon. Every candidate edit is synonymous
in that frame, so the encoded protein is preserved. Two strategies:

  (a) 'codon swap'      : mutate the main-frame codon(s) that overlap the
                          downstream start triplet so it is no longer ATG/GTG/TTG.
  (b) 'RBS synonymous'  : mutate codons overlapping the downstream start's SD /
                          spacer window to weaken its ribosome binding site.

Both report the downstream TIR change *and* verify the main start TIR.
"""

from typing import Callable, List, Optional

from .candidate import Candidate
from .genetic_code import is_start_codon
from .mutations import synonymous_variants
from .scoring import OstirScorer
from .constraints import analyze, ConstraintConfig


def _fmt_edits(edits) -> str:
    return ', '.join(f'{p}:{o}->{n}' for p, o, n in edits)


def _build_candidates(variants, seq, main_start, down_start, scorer,
                      down_baseline, main_baseline, edit_type, cfg, note_fn):
    """Batch-score a set of synonymous variants and build Candidates."""
    d0 = down_start - 1
    down_tirs = scorer.score_many([(v[0], down_start) for v in variants])
    main_tirs = scorer.score_many([(v[0], main_start) for v in variants])
    out: List[Candidate] = []
    for (new_seq, edits), dtir, mtir in zip(variants, down_tirs, main_tirs):
        out.append(Candidate(
            edit_type=edit_type, position=down_start,
            change=_fmt_edits(edits), new_sequence=new_seq,
            tir=dtir, baseline_tir=down_baseline, target_start=down_start,
            main_tir=mtir, main_baseline_tir=main_baseline,
            notes=note_fn(new_seq, edits),
            warnings=analyze(new_seq, seq, frame_start=main_start, cfg=cfg),
        ))
    return out


def codon_swap(seq: str, main_start: int, down_start: int, scorer: OstirScorer,
               main_baseline: Optional[float] = None,
               max_variants: int = 200, cfg: Optional[ConstraintConfig] = None,
               progress: Optional[Callable[[int, int, str], None]] = None
               ) -> List[Candidate]:
    """Strategy (a): destroy the downstream start triplet synonymously."""
    cfg = cfg or ConstraintConfig()
    frame0 = main_start - 1
    d0 = down_start - 1                      # 0-based index of downstream start
    down_baseline = scorer.score_start(seq, down_start)
    if main_baseline is None:
        main_baseline = scorer.score_start(seq, main_start)

    if progress:
        progress(0, 1, 'codon swap')
    variants = list(synonymous_variants(
        seq, frame0, d0, d0 + 3, max_variants=max_variants,
        only_changing=(d0, d0 + 3),
    ))

    def note_fn(new_seq, edits):
        triplet = new_seq[d0:d0 + 3]
        return ('start codon removed' if not is_start_codon(triplet)
                else f'still {triplet}')

    out = _build_candidates(variants, seq, main_start, down_start, scorer,
                            down_baseline, main_baseline, 'codon swap',
                            cfg, note_fn)
    out.sort(key=lambda c: (c.tir, abs(c.main_delta or 0.0)))
    return out


def rbs_synonymous(seq: str, main_start: int, down_start: int,
                   scorer: OstirScorer, main_baseline: Optional[float] = None,
                   rbs_window: int = 20, max_variants: int = 400,
                   cfg: Optional[ConstraintConfig] = None,
                   progress: Optional[Callable[[int, int, str], None]] = None
                   ) -> List[Candidate]:
    """Strategy (b): weaken the downstream SD/spacer synonymously."""
    cfg = cfg or ConstraintConfig()
    frame0 = main_start - 1
    d0 = down_start - 1
    down_baseline = scorer.score_start(seq, down_start)
    if main_baseline is None:
        main_baseline = scorer.score_start(seq, main_start)

    idx_lo = max(frame0, d0 - rbs_window)
    idx_hi = d0
    if idx_hi <= idx_lo:
        return []

    if progress:
        progress(0, 1, 'RBS synonymous')
    variants = list(synonymous_variants(
        seq, frame0, idx_lo, idx_hi, max_variants=max_variants,
    ))
    out = _build_candidates(
        variants, seq, main_start, down_start, scorer,
        down_baseline, main_baseline, 'RBS synonymous', cfg,
        lambda ns, edits: f'{len(edits)} synonymous change(s)')
    out.sort(key=lambda c: (c.tir, abs(c.main_delta or 0.0)))
    return out


def suppress_downstream(seq: str, main_start: int, down_start: int,
                        scorer: OstirScorer, strategy: str = 'both',
                        rbs_window: int = 20,
                        cfg: Optional[ConstraintConfig] = None,
                        progress: Optional[Callable[[int, int, str], None]] = None
                        ) -> List[Candidate]:
    """Run the requested decrease strategy/strategies and return a combined,
    ranked candidate list (lowest downstream TIR first)."""
    cfg = cfg or ConstraintConfig()
    main_baseline = scorer.score_start(seq, main_start)
    results: List[Candidate] = []
    if strategy in ('codon swap', 'both'):
        results += codon_swap(seq, main_start, down_start, scorer,
                              main_baseline, cfg=cfg, progress=progress)
    if strategy in ('RBS synonymous', 'both'):
        results += rbs_synonymous(seq, main_start, down_start, scorer,
                                 main_baseline, rbs_window=rbs_window,
                                 cfg=cfg, progress=progress)
    results.sort(key=lambda c: (c.tir, abs(c.main_delta or 0.0)))
    return results
