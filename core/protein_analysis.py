"""
Protein property calculations using BioPython's ProtParam module.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis


def analyze_protein(aa_sequence: str) -> dict:
    """
    Calculate biophysical properties for an amino acid sequence.

    Parameters
    ----------
    aa_sequence : str
        Single-letter amino acid sequence. Stop codon symbols ('*', '?')
        are stripped before analysis.

    Returns
    -------
    dict with keys:
        length           : int   – number of amino acids
        molecular_weight : float – MW in kDa (None on error)
        isoelectric_point: float – pI (None on error)
        gravy            : float – GRAVY hydrophobicity score (None on error)
        aa_composition   : dict  – {one_letter_code: fraction} (None on error)
        error            : str   – error message if calculation failed
    """
    clean_seq = aa_sequence.replace('*', '').replace('?', '').replace('-', '')

    result = {
        'length': len(clean_seq),
        'molecular_weight': None,
        'isoelectric_point': None,
        'gravy': None,
        'aa_composition': {},
        'error': None,
    }

    if not clean_seq:
        result['error'] = 'Empty sequence'
        return result

    try:
        pa = ProteinAnalysis(clean_seq)
        result['molecular_weight'] = round(pa.molecular_weight() / 1000, 3)
        result['isoelectric_point'] = round(pa.isoelectric_point(), 2)
        result['gravy'] = round(pa.gravy(), 3)
        result['aa_composition'] = {
            aa: round(frac * 100, 1)
            for aa, frac in pa.get_amino_acids_percent().items()
            if frac > 0
        }
    except Exception as exc:
        result['error'] = str(exc)

    return result
