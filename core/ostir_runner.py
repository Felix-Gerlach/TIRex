"""
QThread wrapper around OSTIR's run_ostir() function.
Emits enriched results that include ORF extents and protein properties.
"""

import re

from PyQt6.QtCore import QThread, pyqtSignal

from ostir import run_ostir
from .orf_finder import find_orf
from .protein_analysis import analyze_protein


class OSTIRRunner(QThread):
    """
    Runs OSTIR analysis in a background thread so the GUI stays responsive.

    Signals
    -------
    results_ready(list)  : emitted with list of enriched result dicts
    error_occurred(str)  : emitted if an exception is raised
    progress_update(str) : status messages during processing
    """

    results_ready = pyqtSignal(list)
    error_occurred = pyqtSignal(str)
    progress_update = pyqtSignal(str)

    def __init__(self, sequence: str, params: dict, parent=None):
        super().__init__(parent)
        self.sequence = sequence
        self.params = params

    @staticmethod
    def _sanitize_seq(seq: str) -> str:
        return re.sub(r'[^ACGTU]', '', (seq or '').upper().replace('U', 'T'))

    @staticmethod
    def _valid_asd(asd):
        """Return a valid 9-nt anti-SD or None. Guards against OSTIR's broken
        error handler, which crashes on invalid aSD/sequence input."""
        if not asd:
            return None
        a = re.sub(r'[^ACGT]', '', asd.upper().replace('U', 'T'))
        return a if len(a) == 9 else None

    def run(self):
        try:
            self.progress_update.emit('Running OSTIR TIR prediction…')

            # Defensive sanitisation: OSTIR's own input-validation path is
            # buggy (it calls rich.print(style=...) which raises), so we must
            # never hand it a non-A/C/G/T sequence or a non-9-nt anti-SD.
            self.sequence = self._sanitize_seq(self.sequence)
            if not self.sequence:
                self.error_occurred.emit(
                    'No valid A/C/G/T nucleotides in the input sequence.')
                return
            primary_asd = self._valid_asd(self.params.get('aSD')) or None
            if self.params.get('aSD') and primary_asd is None:
                self.error_occurred.emit(
                    'The anti-SD must be exactly 9 bases (A/C/G/T).')
                return

            raw_results = run_ostir(
                self.sequence,
                start=self.params.get('start'),
                end=self.params.get('end'),
                name=self.params.get('name', 'sequence'),
                aSD=primary_asd,
                threads=self.params.get('threads', 1),
                decimal_places=self.params.get('decimal_places', 4),
                circular=self.params.get('circular', False),
                constraints=self.params.get('constraints') or None,
                verbosity=0,
            )

            total = len(raw_results)
            self.progress_update.emit(
                f'Found {total} start codon(s). Analysing ORFs…'
            )

            circular = self.params.get('circular', False)
            enriched = []

            for idx, result in enumerate(raw_results):
                self.progress_update.emit(
                    f'Analysing ORF {idx + 1}/{total}…'
                )
                start_pos = result['start_position']
                dna_seq, aa_seq, end_pos, has_stop = find_orf(
                    self.sequence, start_pos, circular
                )

                props = analyze_protein(aa_seq) if aa_seq else {
                    'length': 0, 'molecular_weight': None,
                    'isoelectric_point': None, 'gravy': None,
                    'aa_composition': {}, 'error': 'No AA sequence',
                }

                enriched.append({
                    **result,
                    'orf_index': idx,
                    'visible': True,
                    'dna_sequence': dna_seq,
                    'aa_sequence': aa_seq,
                    'end_position': end_pos,
                    'has_stop': has_stop,
                    **{f'protein_{k}': v for k, v in props.items()},
                })

            # ---- Optional second pass with a secondary anti-SD ---------
            secondary = self._valid_asd(self.params.get('aSD_secondary'))
            primary = primary_asd
            if secondary and secondary != primary:
                self.progress_update.emit(
                    'Scoring against secondary anti-SD…'
                )
                try:
                    sec_results = run_ostir(
                        self.sequence,
                        start=self.params.get('start'),
                        end=self.params.get('end'),
                        name=self.params.get('name', 'sequence'),
                        aSD=secondary,
                        threads=self.params.get('threads', 1),
                        decimal_places=self.params.get('decimal_places', 4),
                        circular=circular,
                        constraints=self.params.get('constraints') or None,
                        verbosity=0,
                    ) or []
                    # Match by start position.
                    by_pos = {int(r['start_position']): r for r in sec_results}
                    for e in enriched:
                        sp = int(e.get('start_position', -1))
                        sr = by_pos.get(sp)
                        if sr is not None:
                            e['expression_secondary'] = sr.get('expression')
                            e['dG_total_secondary'] = sr.get('dG_total')
                        else:
                            e['expression_secondary'] = None
                            e['dG_total_secondary'] = None
                    # Annotate which anti-SDs were used (handy for export).
                    for e in enriched:
                        e['aSD_primary'] = primary
                        e['aSD_secondary'] = secondary
                except Exception:
                    # Secondary scoring is best-effort; never fail the run.
                    for e in enriched:
                        e.setdefault('expression_secondary', None)

            self.progress_update.emit(
                f'Done — {total} ORF(s) analysed.'
            )
            self.results_ready.emit(enriched)

        except Exception as exc:
            self.error_occurred.emit(str(exc))
