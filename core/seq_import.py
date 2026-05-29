"""
Sequence file import.

Supported now:
  * FASTA            (.fasta .fa .fna .txt)
  * GenBank          (.gb .gbk .genbank)  — sequence + CDS feature coordinates
  * SnapGene         (.dna)               — via BioPython 'snapgene' (if available)

Pending a sample file:
  * VectorBee        (.vbee)              — format unknown; raises a clear error

Returns an ImportedSeq with the cleaned nucleotide sequence, a display name,
and any CDS features found (1-based start positions on the forward strand),
which the UI can use to pre-seed ORFs/targets.
"""

import os
import re
from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class ImportedSeq:
    sequence: str                       # cleaned A/C/G/T
    name: str = 'sequence'
    cds: List[Tuple[int, int, str]] = field(default_factory=list)
    # each CDS: (start_1based, end_1based, label) on the forward strand
    source_format: str = ''
    note: str = ''


def _clean(seq: str) -> str:
    return re.sub(r'[^ACGTU]', '', (seq or '').upper().replace('U', 'T'))


def _from_fasta_text(text: str, default_name: str) -> ImportedSeq:
    lines = text.strip().splitlines()
    name = default_name
    body_lines = []
    for ln in lines:
        if ln.startswith('>'):
            if name == default_name:
                name = ln[1:].strip().split()[0] or default_name
        else:
            body_lines.append(ln)
    return ImportedSeq(sequence=_clean(''.join(body_lines)),
                       name=name, source_format='fasta')


def _from_biopython(path: str, fmt: str, default_name: str) -> ImportedSeq:
    from Bio import SeqIO
    rec = next(SeqIO.parse(path, fmt))
    seq = _clean(str(rec.seq))
    cds = []
    for feat in getattr(rec, 'features', []):
        if feat.type == 'CDS':
            try:
                start = int(feat.location.start) + 1      # 1-based
                end = int(feat.location.end)
                strand = feat.location.strand
                label = (feat.qualifiers.get('gene')
                         or feat.qualifiers.get('product')
                         or feat.qualifiers.get('label')
                         or ['CDS'])[0]
                if strand == 1:                           # forward only
                    cds.append((start, end, str(label)))
            except Exception:
                continue
    name = (rec.name or rec.id or default_name)
    return ImportedSeq(sequence=seq, name=str(name), cds=cds,
                       source_format=fmt)


def load_sequence(path: str) -> ImportedSeq:
    """Load a sequence file, dispatching on extension."""
    ext = os.path.splitext(path)[1].lower()
    default_name = os.path.splitext(os.path.basename(path))[0]

    if ext in ('.fasta', '.fa', '.fna', '.txt', ''):
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            text = fh.read()
        # GenBank disguised as .txt?  (starts with LOCUS)
        if text.lstrip().startswith('LOCUS'):
            return _from_biopython(path, 'genbank', default_name)
        return _from_fasta_text(text, default_name)

    if ext in ('.gb', '.gbk', '.genbank'):
        return _from_biopython(path, 'genbank', default_name)

    if ext == '.dna':
        try:
            return _from_biopython(path, 'snapgene', default_name)
        except Exception as exc:
            raise ValueError(
                'Could not read SnapGene .dna file. Your BioPython build may '
                f'lack the snapgene parser ({exc}).')

    if ext == '.vbee':
        raise ValueError(
            "VectorBee (.vbee) import isn't implemented yet — the format isn't "
            "documented here. Please share a sample .vbee file so it can be "
            "supported; meanwhile export to GenBank or FASTA.")

    # Fallback: try FASTA text.
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        return _from_fasta_text(fh.read(), default_name)


SUPPORTED_FILTER = (
    'Sequence files (*.fasta *.fa *.fna *.gb *.gbk *.genbank *.dna *.txt);;'
    'FASTA (*.fasta *.fa *.fna);;GenBank (*.gb *.gbk *.genbank);;'
    'SnapGene (*.dna);;All files (*)'
)
