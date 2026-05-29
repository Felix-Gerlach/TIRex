"""
Q5 / NEBaseChanger-style mutagenesis primer design.

Given an original sequence and an edited sequence, design a back-to-back,
non-overlapping primer pair that introduces the edit by whole-plasmid PCR +
KLD (kinase/ligase/DpnI). Convention:

  * The 5' ends of the two primers meet at the edit site.
  * Forward primer extends 3' into the sequence downstream of the edit and
    carries the substitution/insertion at/near its 5' end.
  * Reverse primer is the reverse complement of the region immediately 5'
    (upstream) of the edit.

Tm is an estimate (NEB recommends their Tm calculator for Q5); annealing length
is grown to reach a target Tm on the homology portion.
"""

from dataclasses import dataclass
from typing import Optional

_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def revcomp(seq: str) -> str:
    return ''.join(_COMP.get(b, 'N') for b in reversed(seq.upper()))


def _gc(seq: str) -> float:
    return 100.0 * sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0.0


def tm_estimate(seq: str) -> float:
    """Rough primer Tm (°C). Wallace rule for <14 nt, else a GC%-based formula."""
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    if n < 14:
        a = seq.count('A') + seq.count('T')
        g = seq.count('G') + seq.count('C')
        return 2 * a + 4 * g
    gc = _gc(seq)
    return 64.9 + 41 * (gc / 100 * n - 16.4) / n


def _grow_to_tm(seq_from_5p: str, target_tm: float = 62.0,
                min_len: int = 16, max_len: int = 36) -> str:
    """Return the shortest prefix (>=min_len) of seq_from_5p whose Tm >= target,
    capped at max_len."""
    best = seq_from_5p[:min_len]
    for L in range(min_len, min(max_len, len(seq_from_5p)) + 1):
        cand = seq_from_5p[:L]
        best = cand
        if tm_estimate(cand) >= target_tm:
            break
    return best


def _edit_span(ref: str, new: str):
    """(lo, hi_new, hi_ref): common-prefix/suffix diff. Edit occupies
    new[lo:hi_new] and ref[lo:hi_ref]."""
    p = 0
    while p < min(len(ref), len(new)) and ref[p] == new[p]:
        p += 1
    s = 0
    while (s < min(len(ref), len(new)) - p
           and ref[len(ref) - 1 - s] == new[len(new) - 1 - s]):
        s += 1
    return p, len(new) - s, len(ref) - s


@dataclass
class PrimerPair:
    forward: str
    reverse: str
    fwd_tm: float
    rev_tm: float
    fwd_len: int
    rev_len: int
    kind: str
    note: str = ''

    def as_dict(self):
        return {
            'forward (5\'->3\')': self.forward,
            'reverse (5\'->3\')': self.reverse,
            'fwd_len': self.fwd_len, 'rev_len': self.rev_len,
            'fwd_Tm': round(self.fwd_tm, 1), 'rev_Tm': round(self.rev_tm, 1),
            'kind': self.kind, 'note': self.note,
        }


def design_primers(original: str, mutated: str, target_tm: float = 62.0,
                   circular: bool = True) -> Optional[PrimerPair]:
    """Design a Q5-style back-to-back primer pair realising original->mutated."""
    original = original.upper()
    mutated = mutated.upper()
    if original == mutated:
        return None

    lo, hi_new, hi_ref = _edit_span(original, mutated)
    ins = mutated[lo:hi_new]          # changed/inserted bases (in mutant)
    del_len = hi_ref - lo             # bases removed from original
    sub_len = hi_new - lo             # bases present at edit site in mutant

    up = mutated[:lo]                 # unchanged region 5' of the edit
    down = mutated[hi_new:]           # unchanged region 3' of the edit
    if not up or not down:
        return None                   # edit at a terminus; needs circular context

    # Reverse primer: revcomp of upstream homology, 5' end at the edit junction.
    reverse = _grow_to_tm(revcomp(up), target_tm)

    # Classify the edit and build the forward primer (5' end at the edit).
    classify = ('substitution' if sub_len == del_len and sub_len > 0 else
                'insertion' if del_len == 0 else
                'deletion' if sub_len == 0 else 'replacement')
    if classify == 'deletion':
        forward = _grow_to_tm(down, target_tm)
    else:
        # changed/inserted bases form the 5' end, then anneal into downstream
        anneal = _grow_to_tm(down, target_tm, min_len=max(8, 16 - sub_len))
        forward = ins + anneal

    return PrimerPair(
        forward=forward, reverse=reverse,
        fwd_tm=tm_estimate(forward[len(ins):] or forward),
        rev_tm=tm_estimate(reverse),
        fwd_len=len(forward), rev_len=len(reverse),
        kind=classify,
        note='Whole-plasmid PCR with Q5 + KLD (NEBaseChanger style). '
             'Tm is estimated — verify with the NEB Tm calculator.')
