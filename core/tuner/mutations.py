"""
Edit generators.

Increase-mode edits operate on an upstream window [w_lo, w_hi) (0-based,
half-open) that lies strictly *before* the start codon. Each generator yields
tuples:  (edit_type, position_1based, change_str, new_seq, new_start_1based)

Decrease-mode synonymous enumeration is in `synonymous_variants`.
"""

import itertools
from typing import Iterator, List, Tuple

from .genetic_code import CODON_TABLE, synonymous_codons

BASES = ('A', 'C', 'G', 'T')

Edit = Tuple[str, int, str, str, int]


# ====================================================================== #
#  Increase-mode generators
# ====================================================================== #
def gen_substitutions(seq: str, w_lo: int, w_hi: int,
                      target_start: int) -> Iterator[Edit]:
    for i in range(w_lo, w_hi):
        old = seq[i]
        for b in BASES:
            if b == old:
                continue
            new_seq = seq[:i] + b + seq[i + 1:]
            yield ('substitution', i + 1, f'{old}->{b}', new_seq, target_start)


def gen_deletions(seq: str, w_lo: int, w_hi: int,
                  target_start: int) -> Iterator[Edit]:
    for i in range(w_lo, w_hi):
        old = seq[i]
        new_seq = seq[:i] + seq[i + 1:]
        # everything downstream (incl. start codon) shifts left by 1
        yield ('deletion', i + 1, f'del {old}', new_seq, target_start - 1)


def gen_insertions(seq: str, w_lo: int, w_hi: int,
                   target_start: int) -> Iterator[Edit]:
    # insertion sites are *between* bases: i in [w_lo, w_hi]
    for i in range(w_lo, w_hi + 1):
        for b in BASES:
            new_seq = seq[:i] + b + seq[i:]
            # start codon shifts right by 1
            yield ('insertion', i + 1, f'ins {b}', new_seq, target_start + 1)


def gen_dinucleotides(seq: str, w_lo: int, w_hi: int,
                      target_start: int) -> Iterator[Edit]:
    for i in range(w_lo, w_hi - 1):
        old = seq[i:i + 2]
        for b1, b2 in itertools.product(BASES, BASES):
            new = b1 + b2
            if new == old:
                continue
            new_seq = seq[:i] + new + seq[i + 2:]
            yield ('dinucleotide', i + 1, f'{old}->{new}', new_seq, target_start)


# ====================================================================== #
#  Decrease-mode: synonymous enumeration
# ====================================================================== #
def _codon_span(idx_lo: int, idx_hi: int, frame_start0: int) -> Tuple[int, int]:
    """Return [c_lo, c_hi) codon-start indices (0-based, frame-aligned) that
    fully cover the nucleotide window [idx_lo, idx_hi)."""
    rel_lo = idx_lo - frame_start0
    rel_hi = idx_hi - frame_start0
    c_lo = frame_start0 + (rel_lo // 3) * 3
    c_hi = frame_start0 + ((rel_hi + 2) // 3) * 3
    return c_lo, c_hi


def synonymous_variants(seq: str, frame_start0: int, idx_lo: int, idx_hi: int,
                        max_variants: int = 400,
                        only_changing: Tuple[int, int] = None
                        ) -> Iterator[Tuple[str, List[Tuple[int, str, str]]]]:
    """
    Yield (new_seq, edits) where only nucleotides in [idx_lo, idx_hi) differ
    from `seq`, every codon in the affected frame keeps its amino acid, and
    nucleotides outside [idx_lo, idx_hi) are held fixed.

    `edits` is a list of (position_1based, old_base, new_base).
    If `only_changing` (a,b) is given, at least one changed base must fall in
    [a, b) — used to require the start-codon triplet itself to mutate.
    """
    frame_start0 = frame_start0 % 3 if frame_start0 < 0 else frame_start0
    c_lo, c_hi = _codon_span(idx_lo, idx_hi, frame_start0)
    c_lo = max(c_lo, frame_start0)

    # For each affected codon, list synonymous codons compatible with the
    # fixed (out-of-window) positions of that codon.
    per_codon_options: List[List[str]] = []
    codon_starts: List[int] = []
    for cs in range(c_lo, c_hi, 3):
        if cs + 3 > len(seq):
            break
        orig = seq[cs:cs + 3]
        if orig not in CODON_TABLE:
            per_codon_options.append([orig])
            codon_starts.append(cs)
            continue
        opts = []
        for cand in synonymous_codons(orig):
            ok = True
            for k in range(3):
                pos = cs + k
                # positions outside the editable window must stay as original
                if not (idx_lo <= pos < idx_hi) and cand[k] != orig[k]:
                    ok = False
                    break
            if ok:
                opts.append(cand)
        if orig not in opts:
            opts.append(orig)
        per_codon_options.append(opts)
        codon_starts.append(cs)

    if not codon_starts:
        return

    count = 0
    for combo in itertools.product(*per_codon_options):
        if count >= max_variants:
            break
        new_seq = list(seq)
        edits = []
        changed_in_target = False
        for cs, cand in zip(codon_starts, combo):
            for k in range(3):
                pos = cs + k
                if cand[k] != seq[pos]:
                    new_seq[pos] = cand[k]
                    edits.append((pos + 1, seq[pos], cand[k]))
                    if only_changing and only_changing[0] <= pos < only_changing[1]:
                        changed_in_target = True
        if not edits:
            continue
        if only_changing and not changed_in_target:
            continue
        count += 1
        yield (''.join(new_seq), edits)
