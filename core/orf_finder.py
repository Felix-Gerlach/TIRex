"""
Find ORF extents (start codon to next in-frame stop codon) from a DNA sequence.
Supports linear and circular sequences.
"""

CODON_TABLE = {
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

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
START_CODONS = {'ATG', 'GTG', 'TTG'}


def find_orf(sequence: str, start_pos_1indexed: int, circular: bool = False):
    """
    Scan from start_pos_1indexed in-frame until first stop codon.

    Parameters
    ----------
    sequence         : DNA sequence string (T-based)
    start_pos_1indexed: 1-indexed position of the start codon
    circular         : if True, wraps around the end of the sequence

    Returns
    -------
    dna_seq  : str  – full ORF DNA including stop codon (or up to seq end)
    aa_seq   : str  – translated amino acids (stop codon excluded)
    end_pos  : int  – 1-indexed position of first base AFTER the stop codon
                      (or len(sequence) + 1 if no stop found)
    has_stop : bool – whether a stop codon was found
    """
    seq = sequence.upper().replace('U', 'T')
    seq_len = len(seq)
    start_0 = start_pos_1indexed - 1  # convert to 0-indexed

    codons_dna = []
    aa_list = []
    has_stop = False
    end_pos = seq_len + 1  # default: runs to end

    pos = start_0
    max_codons = seq_len  # safety limit to prevent infinite loops on circular

    while max_codons > 0:
        max_codons -= 1

        if circular:
            idx = pos % seq_len
            # make sure we can read 3 bases wrapping around
            if seq_len < 3:
                break
            codon = (seq[idx % seq_len] +
                     seq[(idx + 1) % seq_len] +
                     seq[(idx + 2) % seq_len])
        else:
            if pos + 3 > seq_len:
                break
            codon = seq[pos:pos + 3]

        codons_dna.append(codon)
        aa = CODON_TABLE.get(codon, '?')

        if aa == '*':
            has_stop = True
            end_pos = (pos % seq_len if circular else pos) + 3 + 1  # 1-indexed base after stop
            break
        else:
            aa_list.append(aa)

        pos += 3

        # For linear: if we've passed the end, stop
        if not circular and pos >= seq_len:
            break

        # For circular: stop when we've looped back past start
        if circular and pos >= start_0 + seq_len:
            break

    return ''.join(codons_dna), ''.join(aa_list), end_pos, has_stop
