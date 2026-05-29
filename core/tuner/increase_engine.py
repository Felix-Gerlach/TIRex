"""
Increase mode: raise the TIR at the target start codon by editing the
upstream window between the SD/RBS and the start codon.

Two strategies, user-selectable:
  * 'enumerate'  : independent single subs/dels/ins + dinucleotide subs (depth 1)
  * 'greedy'     : repeatedly apply the best edit until a target fold-change
                   is reached or no further improvement is found.
"""

from typing import Callable, List, Optional

from .candidate import Candidate
from .mutations import (
    gen_substitutions, gen_deletions, gen_insertions, gen_dinucleotides,
)
from .scoring import OstirScorer
from .constraints import analyze, ConstraintConfig


def _window(target_start: int, upstream_len: int, seq_len: int):
    """0-based half-open window [w_lo, w_hi) ending just before the start
    codon (target_start is 1-based)."""
    start0 = target_start - 1
    w_hi = start0
    w_lo = max(0, start0 - upstream_len)
    return w_lo, w_hi


def enumerate_edits(seq: str, target_start: int, scorer: OstirScorer,
                    upstream_len: int = 20,
                    include_dinucleotide: bool = True,
                    cfg: Optional[ConstraintConfig] = None,
                    progress: Optional[Callable[[int, int, str], None]] = None
                    ) -> List[Candidate]:
    """Independent depth-1 edits across the upstream window, ranked by TIR."""
    cfg = cfg or ConstraintConfig()
    baseline = scorer.score_start(seq, target_start)
    w_lo, w_hi = _window(target_start, upstream_len, len(seq))

    generators = [
        gen_substitutions(seq, w_lo, w_hi, target_start),
        gen_deletions(seq, w_lo, w_hi, target_start),
        gen_insertions(seq, w_lo, w_hi, target_start),
    ]
    if include_dinucleotide:
        generators.append(gen_dinucleotides(seq, w_lo, w_hi, target_start))

    edits = [e for g in generators for e in g]
    # Batch-score all candidate sequences (parallel when worthwhile).
    tirs = scorer.score_many([(e[3], e[4]) for e in edits], progress=progress)

    candidates: List[Candidate] = []
    for (etype, pos, change, new_seq, new_start), tir in zip(edits, tirs):
        candidates.append(Candidate(
            edit_type=etype, position=pos, change=change,
            new_sequence=new_seq, tir=tir, baseline_tir=baseline,
            target_start=new_start,
            warnings=analyze(new_seq, seq, frame_start=new_start, cfg=cfg),
        ))
    candidates.sort(key=lambda c: c.tir, reverse=True)
    return candidates


def greedy_optimize(seq: str, target_start: int, scorer: OstirScorer,
                    upstream_len: int = 20,
                    include_dinucleotide: bool = True,
                    target_fold: float = float('inf'),
                    max_rounds: int = 6,
                    cfg: Optional[ConstraintConfig] = None,
                    progress: Optional[Callable[[int, int, str], None]] = None
                    ) -> List[Candidate]:
    """Hill-climb: each round apply the single best-improving edit. Returns the
    cumulative trajectory (one Candidate per accepted round)."""
    cfg = cfg or ConstraintConfig()
    cur_seq = seq
    cur_start = target_start
    base0 = scorer.score_start(seq, target_start)
    history: List[str] = []
    trajectory: List[Candidate] = []

    for rnd in range(max_rounds):
        if progress:
            progress(rnd, max_rounds, f'greedy round {rnd + 1}')
        round_cands = enumerate_edits(
            cur_seq, cur_start, scorer, upstream_len,
            include_dinucleotide, cfg=cfg, progress=None,
        )
        if not round_cands:
            break
        best = round_cands[0]
        cur_round_base = scorer.score_start(cur_seq, cur_start)
        if best.tir <= cur_round_base + 1e-9:
            break   # no improvement
        history = history + [f'{best.edit_type} {best.change}@{best.position}']
        accepted = Candidate(
            edit_type=best.edit_type, position=best.position,
            change=best.change, new_sequence=best.new_sequence,
            tir=best.tir, baseline_tir=base0, target_start=best.target_start,
            history=list(history),
            notes=f'round {rnd + 1}', warnings=list(best.warnings),
        )
        trajectory.append(accepted)
        cur_seq = best.new_sequence
        cur_start = best.target_start
        if base0 > 0 and best.tir / base0 >= target_fold:
            break
    return trajectory


def beam_optimize(seq: str, target_start: int, scorer: OstirScorer,
                  upstream_len: int = 20, include_dinucleotide: bool = True,
                  beam_width: int = 5, max_rounds: int = 4,
                  target_fold: float = float('inf'),
                  cfg: Optional[ConstraintConfig] = None,
                  progress: Optional[Callable[[int, int, str], None]] = None
                  ) -> List[Candidate]:
    """Beam search: keep the best `beam_width` partial solutions each round and
    expand each by one more edit. Returns the best trajectory found (one
    Candidate per round, cumulative)."""
    cfg = cfg or ConstraintConfig()
    base0 = scorer.score_start(seq, target_start)

    # Each beam item: (seq, start, tir, [Candidate history])
    beam = [(seq, target_start, base0, [])]
    best_overall = None

    for rnd in range(max_rounds):
        if progress:
            progress(rnd, max_rounds, f'beam round {rnd + 1}')
        expansions = []
        for (bseq, bstart, btir, bhist) in beam:
            cands = enumerate_edits(bseq, bstart, scorer, upstream_len,
                                    include_dinucleotide, cfg=cfg)
            for c in cands[:beam_width]:
                if c.tir <= btir + 1e-9:
                    continue
                step = Candidate(
                    edit_type=c.edit_type, position=c.position,
                    change=c.change, new_sequence=c.new_sequence,
                    tir=c.tir, baseline_tir=base0,
                    target_start=c.target_start,
                    history=[*[h.change for h in bhist], c.change],
                    notes=f'round {rnd + 1}', warnings=list(c.warnings))
                expansions.append((c.new_sequence, c.target_start, c.tir,
                                   bhist + [step]))
        if not expansions:
            break
        # Dedup by sequence, keep best tir, take top beam_width.
        seen = {}
        for item in expansions:
            s = item[0]
            if s not in seen or item[2] > seen[s][2]:
                seen[s] = item
        beam = sorted(seen.values(), key=lambda x: x[2], reverse=True)[:beam_width]
        top = beam[0]
        if best_overall is None or top[2] > best_overall[2]:
            best_overall = top
        if base0 > 0 and top[2] / base0 >= target_fold:
            break

    return best_overall[3] if best_overall else []
