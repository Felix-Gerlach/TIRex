"""
QThread wrapper around OSTIR's run_ostir() function.
Emits enriched results that include ORF extents and protein properties.
"""

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

    def run(self):
        try:
            self.progress_update.emit('Running OSTIR TIR prediction…')

            raw_results = run_ostir(
                self.sequence,
                start=self.params.get('start'),
                end=self.params.get('end'),
                name=self.params.get('name', 'sequence'),
                aSD=self.params.get('aSD') or None,
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

            self.progress_update.emit(
                f'Done — {total} ORF(s) analysed.'
            )
            self.results_ready.emit(enriched)

        except Exception as exc:
            self.error_occurred.emit(str(exc))
