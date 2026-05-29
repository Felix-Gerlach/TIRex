"""Background worker for the TIR-Tuner dialog (keeps the UI responsive)."""

from PyQt6.QtCore import QThread, pyqtSignal

from core.tuner.scoring import OstirScorer
from core.tuner.increase_engine import (
    enumerate_edits, greedy_optimize, beam_optimize,
)
from core.tuner.decrease_engine import suppress_downstream
from core.tuner.library import design_library


class TunerWorker(QThread):
    progress = pyqtSignal(int, int, str)     # done, total, message
    finished_ok = pyqtSignal(list, dict)     # candidates, stats
    failed = pyqtSignal(str)

    def __init__(self, seq, mode, params, parent=None):
        super().__init__(parent)
        self.seq = seq
        self.mode = mode               # 'increase' | 'decrease' | 'library'
        self.params = params

    def _emit(self, done, total, msg):
        self.progress.emit(done, total, msg)

    def run(self):
        try:
            p = self.params
            scorer = OstirScorer(aSD=p.get('aSD'), threads=p.get('threads', 1))
            if self.mode == 'increase':
                target = int(p['target_start'])
                upstream = int(p.get('upstream_len', 20))
                dinuc = bool(p.get('include_dinucleotide', True))
                depth = p.get('depth')
                if depth == 'greedy':
                    cands = greedy_optimize(
                        self.seq, target, scorer,
                        upstream_len=upstream, include_dinucleotide=dinuc,
                        target_fold=p.get('target_fold', float('inf')),
                        max_rounds=int(p.get('max_rounds', 6)),
                        progress=self._emit,
                    )
                elif depth == 'beam':
                    cands = beam_optimize(
                        self.seq, target, scorer,
                        upstream_len=upstream, include_dinucleotide=dinuc,
                        beam_width=int(p.get('beam_width', 5)),
                        max_rounds=int(p.get('max_rounds', 4)),
                        target_fold=p.get('target_fold', float('inf')),
                        progress=self._emit,
                    )
                else:
                    cands = enumerate_edits(
                        self.seq, target, scorer,
                        upstream_len=upstream, include_dinucleotide=dinuc,
                        progress=self._emit,
                    )
            elif self.mode == 'library':
                cands = design_library(
                    self.seq, int(p['target_start']), scorer,
                    n_variants=int(p.get('n_variants', 8)),
                    upstream_len=int(p.get('upstream_len', 20)),
                    include_dinucleotide=bool(p.get('include_dinucleotide', True)),
                    progress=self._emit,
                )
            else:  # decrease
                cands = suppress_downstream(
                    self.seq,
                    int(p['main_start']), int(p['down_start']), scorer,
                    strategy=p.get('strategy', 'both'),
                    rbs_window=int(p.get('rbs_window', 20)),
                    progress=self._emit,
                )
            stats = {'ostir_calls': scorer.calls, 'cache_hits': scorer.hits,
                     'candidates': len(cands)}
            self.finished_ok.emit(cands, stats)
        except Exception as exc:
            import traceback
            self.failed.emit(f'{exc}\n{traceback.format_exc()}')
