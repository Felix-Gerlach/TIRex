"""
Codon-optimization layer (host-aware, protein-preserving).

Computes CAI (Codon Adaptation Index) and GC content against a host codon-usage
table, and synonymously rewrites a CDS to raise CAI without changing the
encoded protein. The first N codons can be preserved (they overlap the RBS /
translation-initiation region and affect TIR).

Ships an E. coli K-12 table; more hosts can be added to CODON_USAGE.
"""

import math
from typing import Dict, List, Tuple

from .tuner.genetic_code import CODON_TABLE, SYNONYMOUS

# ---------------------------------------------------------------------------
# Host codon usage: fraction of each codon WITHIN its amino-acid group.
# Values are approximate published E. coli K-12 usage (Kazusa-style).
# ---------------------------------------------------------------------------
_ECOLI_FRACTION: Dict[str, float] = {
    'TTT': 0.58, 'TTC': 0.42,
    'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47,
    'ATT': 0.51, 'ATC': 0.39, 'ATA': 0.10,
    'ATG': 1.00,
    'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35,
    'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.24,
    'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49,
    'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.24,
    'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33,
    'TAT': 0.59, 'TAC': 0.41,
    'CAT': 0.57, 'CAC': 0.43,
    'CAA': 0.34, 'CAG': 0.66,
    'AAT': 0.49, 'AAC': 0.51,
    'AAA': 0.74, 'AAG': 0.26,
    'GAT': 0.63, 'GAC': 0.37,
    'GAA': 0.68, 'GAG': 0.32,
    'TGT': 0.46, 'TGC': 0.54,
    'TGG': 1.00,
    'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.10, 'AGA': 0.07, 'AGG': 0.04,
    'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15,
    'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30,
}

CODON_USAGE: Dict[str, Dict[str, float]] = {
    'E. coli K-12': _ECOLI_FRACTION,
}


def _weights(fractions: Dict[str, float]) -> Dict[str, float]:
    """Relative adaptiveness w = fraction / max(fraction in aa group)."""
    by_aa: Dict[str, float] = {}
    for codon, frac in fractions.items():
        aa = CODON_TABLE[codon]
        by_aa[aa] = max(by_aa.get(aa, 0.0), frac)
    return {codon: (frac / by_aa[CODON_TABLE[codon]]) if by_aa[CODON_TABLE[codon]] else 1.0
            for codon, frac in fractions.items()}


def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return 100.0 * sum(1 for b in seq if b in 'GC') / len(seq)


def cai(seq: str, host: str = 'E. coli K-12') -> float:
    """Codon Adaptation Index of a CDS (frame 0). Met/Trp/stop excluded."""
    w = _weights(CODON_USAGE[host])
    logs = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        aa = CODON_TABLE.get(codon)
        if aa in (None, '*', 'M', 'W'):
            continue
        wi = w.get(codon)
        if wi and wi > 0:
            logs.append(math.log(wi))
    if not logs:
        return 1.0
    return math.exp(sum(logs) / len(logs))


def best_codon(aa: str, host: str = 'E. coli K-12') -> str:
    """Highest-usage codon for an amino acid."""
    fr = CODON_USAGE[host]
    cands = SYNONYMOUS.get(aa, [])
    return max(cands, key=lambda c: fr.get(c, 0.0)) if cands else ''


def optimize_cds(seq: str, host: str = 'E. coli K-12',
                 preserve_first: int = 0) -> Tuple[str, dict]:
    """Synonymously rewrite a CDS to maximise CAI (greedy max-weight codon per
    position). `preserve_first` codons are left untouched. Returns
    (optimized_seq, stats)."""
    seq = seq.upper()
    out = []
    changed = 0
    n_codons = len(seq) // 3
    for ci in range(n_codons):
        codon = seq[ci * 3:ci * 3 + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None or aa == '*' or ci < preserve_first:
            out.append(codon)
            continue
        bc = best_codon(aa, host)
        if bc and bc != codon:
            changed += 1
            out.append(bc)
        else:
            out.append(codon)
    # keep any trailing partial bases unchanged
    tail = seq[n_codons * 3:]
    new_seq = ''.join(out) + tail
    stats = {
        'host': host,
        'cai_before': cai(seq, host),
        'cai_after': cai(new_seq, host),
        'gc_before': gc_content(seq),
        'gc_after': gc_content(new_seq),
        'codons_changed': changed,
        'codons_total': n_codons,
        'preserved': preserve_first,
    }
    return new_seq, stats


def translate_equal(a: str, b: str) -> bool:
    from .tuner.genetic_code import translate
    return translate(a) == translate(b)
