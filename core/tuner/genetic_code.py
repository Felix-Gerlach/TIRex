"""
Standard genetic code helpers for TIR-Tuner.

Pure data + functions, no external deps. Used by the mutation engines to
enumerate *synonymous* edits (decrease mode) and to recognise start codons.
"""

from typing import Dict, List

# Bacterial start codons recognised by OSTIR's scan.
START_CODONS = ('ATG', 'GTG', 'TTG')

STOP_CODONS = ('TAA', 'TAG', 'TGA')

# Standard genetic code (DNA, uppercase). '*' = stop.
CODON_TABLE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Reverse map: amino acid -> list of synonymous codons.
SYNONYMOUS: Dict[str, List[str]] = {}
for _codon, _aa in CODON_TABLE.items():
    SYNONYMOUS.setdefault(_aa, []).append(_codon)


def translate(seq: str) -> str:
    """Translate a DNA string (frame 0) to amino acids; trailing partial
    codon ignored. Unknown codons -> 'X'."""
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa.append(CODON_TABLE.get(seq[i:i + 3], 'X'))
    return ''.join(aa)


def synonymous_codons(codon: str) -> List[str]:
    """All codons (including the input) encoding the same amino acid."""
    aa = CODON_TABLE.get(codon.upper())
    if aa is None:
        return [codon.upper()]
    return list(SYNONYMOUS[aa])


def is_start_codon(triplet: str) -> bool:
    return triplet.upper() in START_CODONS


def is_stop_codon(triplet: str) -> bool:
    return triplet.upper() in STOP_CODONS
