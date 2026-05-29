"""
Main application window.
Layout:
  ┌──────────────┬────────────────────────────────────┐
  │  InputPanel  │  VisualizationWidget (top)          │
  │  (left dock) ├────────────────────────────────────┤
  │              │  ORFTableWidget      (bottom)       │
  └──────────────┴────────────────────────────────────┘
"""

import os
import csv
import pandas as pd

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QSplitter, QStatusBar, QMenuBar, QMenu,
    QMessageBox, QFileDialog, QProgressBar, QLabel,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QFont

from .input_panel import InputPanel
from .visualization_widget import VisualizationWidget
from .orf_table_widget import ORFTableWidget
from .sds_page_widget import SDSPageDialog
from .tuner_dialog import TunerDialog
from .codon_optimizer_dialog import CodonOptimizerDialog
from .batch_dialog import BatchDialog
from core.ostir_runner import OSTIRRunner
from core import session as _session


class MainWindow(QMainWindow):

    # ------------------------------------------------------------------ #
    def __init__(self):
        super().__init__()
        self.setWindowTitle('TIRex — Translation Initiation Rate Explorer')
        self.resize(1400, 860)

        self._sequence = ''
        self._results = []
        self._runner = None
        self._target_orf_index = None

        self._build_ui()
        self._build_menu()
        self._build_statusbar()

    # ------------------------------------------------------------------ #
    #  UI                                                                  #
    # ------------------------------------------------------------------ #
    def _build_ui(self):
        # Root splitter: input panel left, main area right
        root_splitter = QSplitter(Qt.Orientation.Horizontal)
        root_splitter.setChildrenCollapsible(False)

        # Left: input panel (fixed reasonable width)
        self.input_panel = InputPanel()
        self.input_panel.setMinimumWidth(290)
        self.input_panel.setMaximumWidth(360)
        self.input_panel.run_requested.connect(self._on_run_requested)
        root_splitter.addWidget(self.input_panel)

        # Right: vertical splitter (viz on top, table on bottom)
        right_splitter = QSplitter(Qt.Orientation.Vertical)
        right_splitter.setChildrenCollapsible(False)

        self.viz_widget = VisualizationWidget()
        right_splitter.addWidget(self.viz_widget)

        self.orf_table = ORFTableWidget()
        self.orf_table.setMinimumHeight(180)
        self.orf_table.visibility_changed.connect(self._on_visibility_changed)
        self.orf_table.selection_changed.connect(self._on_orf_selected)
        self.orf_table.aa_copied.connect(self._on_aa_copied)
        self.orf_table.target_requested.connect(self._on_target_requested)
        self.orf_table.target_cleared.connect(self._on_target_cleared)
        self.orf_table.open_gel_requested.connect(self._open_gel_dialog)
        self.orf_table.tune_requested.connect(self._open_tuner_dialog)
        self.viz_widget.tir_range_changed.connect(self.orf_table.set_tir_range)
        right_splitter.addWidget(self.orf_table)

        right_splitter.setSizes([500, 300])
        root_splitter.addWidget(right_splitter)
        root_splitter.setSizes([290, 1110])

        self.setCentralWidget(root_splitter)

    # ------------------------------------------------------------------ #
    def _build_menu(self):
        mb = self.menuBar()

        # File menu
        file_menu = mb.addMenu('File')

        load_action = QAction('Load sequence (FASTA/GenBank/.dna)…', self)
        load_action.setShortcut('Ctrl+O')
        load_action.triggered.connect(self.input_panel._load_fasta)
        file_menu.addAction(load_action)

        file_menu.addSeparator()

        save_session_action = QAction('Save session…', self)
        save_session_action.setShortcut('Ctrl+S')
        save_session_action.triggered.connect(self._save_session)
        file_menu.addAction(save_session_action)

        open_session_action = QAction('Open session…', self)
        open_session_action.triggered.connect(self._open_session)
        file_menu.addAction(open_session_action)

        file_menu.addSeparator()

        quit_action = QAction('Quit', self)
        quit_action.setShortcut('Ctrl+Q')
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        # Export menu
        export_menu = mb.addMenu('Export')

        csv_action = QAction('Results as CSV…', self)
        csv_action.setShortcut('Ctrl+Shift+S')
        csv_action.triggered.connect(self._export_csv)
        export_menu.addAction(csv_action)

        fasta_action = QAction('Protein sequences as FASTA…', self)
        fasta_action.triggered.connect(self._export_fasta)
        export_menu.addAction(fasta_action)

        export_menu.addSeparator()

        png_action = QAction('Visualization as PNG…', self)
        png_action.triggered.connect(lambda: self._export_figure('png'))
        export_menu.addAction(png_action)

        svg_action = QAction('Visualization as SVG…', self)
        svg_action.triggered.connect(lambda: self._export_figure('svg'))
        export_menu.addAction(svg_action)

        # Tools menu
        tools_menu = mb.addMenu('Tools')
        tuner_action = QAction('Tune translation rate (TIR-Tuner)…', self)
        tuner_action.setShortcut('Ctrl+T')
        tuner_action.triggered.connect(lambda: self._open_tuner_dialog(None))
        tools_menu.addAction(tuner_action)
        codon_action = QAction('Codon optimizer (CAI)…', self)
        codon_action.triggered.connect(self._open_codon_optimizer)
        tools_menu.addAction(codon_action)
        batch_action = QAction('Batch OSTIR (multi-FASTA)…', self)
        batch_action.triggered.connect(self._open_batch_dialog)
        tools_menu.addAction(batch_action)
        tools_menu.addSeparator()
        gel_action = QAction('Simulate SDS-PAGE gel…', self)
        gel_action.triggered.connect(self._open_gel_dialog)
        tools_menu.addAction(gel_action)
        clear_target_action = QAction('Clear target protein', self)
        clear_target_action.triggered.connect(self._on_target_cleared)
        tools_menu.addAction(clear_target_action)

        # Help menu
        help_menu = mb.addMenu('Help')
        about_action = QAction('About TIRex', self)
        about_action.triggered.connect(self._show_about)
        help_menu.addAction(about_action)

    # ------------------------------------------------------------------ #
    def _build_statusbar(self):
        sb = self.statusBar()

        self.status_label = QLabel('Ready')
        self.status_label.setFont(QFont('Segoe UI', 8))
        sb.addWidget(self.status_label, 1)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)   # indeterminate
        self.progress_bar.setFixedWidth(180)
        self.progress_bar.setFixedHeight(14)
        self.progress_bar.setVisible(False)
        sb.addPermanentWidget(self.progress_bar)

    # ------------------------------------------------------------------ #
    #  OSTIR run                                                           #
    # ------------------------------------------------------------------ #
    def _on_run_requested(self, sequence: str, params: dict):
        if self._runner and self._runner.isRunning():
            QMessageBox.information(self, 'Running',
                                    'OSTIR is already running. Please wait.')
            return

        self._sequence = sequence
        self._results = []
        self._target_orf_index = None
        self.input_panel.set_running(True)
        self.progress_bar.setVisible(True)

        self._runner = OSTIRRunner(sequence, params, parent=self)
        self._runner.results_ready.connect(self._on_results_ready)
        self._runner.error_occurred.connect(self._on_error)
        self._runner.progress_update.connect(self._on_progress)
        self._runner.start()

    def _on_progress(self, message: str):
        self.status_label.setText(message)

    def _on_results_ready(self, results: list):
        self._results = results
        self.input_panel.set_running(False)
        self.progress_bar.setVisible(False)

        if not results:
            self.status_label.setText('No start codons found.')
            QMessageBox.information(self, 'No results',
                                    'OSTIR found no start codons in the '
                                    'given range.')
            return

        self.status_label.setText(
            f'Found {len(results)} start codon(s) in {len(self._sequence)} bp sequence.'
        )
        # Populate the table first so it has data when the viz broadcasts
        # its reset TIR range.
        self.orf_table.update_results(results)
        self.viz_widget.set_data(self._sequence, results)

    def _on_error(self, message: str):
        self.input_panel.set_running(False)
        self.progress_bar.setVisible(False)
        self.status_label.setText(f'Error: {message}')
        QMessageBox.critical(self, 'OSTIR Error', message)

    # ------------------------------------------------------------------ #
    #  Visibility & selection                                              #
    # ------------------------------------------------------------------ #
    def _on_visibility_changed(self, orf_index: int, visible: bool):
        self.viz_widget.update_visibility(orf_index, visible)

    def _on_orf_selected(self, orf_index: int):
        self.viz_widget.highlight_orf(orf_index)

    def _on_aa_copied(self, orf_index: int, aa_length: int):
        self.status_label.setText(
            f'Copied {aa_length} aa of ORF #{orf_index + 1} to clipboard.'
        )

    # ------------------------------------------------------------------ #
    #  Target protein                                                      #
    # ------------------------------------------------------------------ #
    def _on_target_requested(self, orf_index: int):
        target = next((r for r in self._results
                       if r.get('orf_index') == orf_index), None)
        if target is None:
            return
        target_start = target.get('start_position')
        target_end = target.get('end_position')
        self._target_orf_index = orf_index
        # Hide every fragment that is NOT in the same reading frame as the
        # target start codon (frame = start position modulo 3).
        for r in self._results:
            sp = r.get('start_position')
            r['visible'] = (sp is not None
                            and (sp - target_start) % 3 == 0)
        # Remove out-of-frame rows from the table entirely (not just uncheck).
        self.orf_table.set_frame_filter(target_start)
        # Push visibility into viz for every ORF
        for r in self._results:
            self.viz_widget.update_visibility(
                r.get('orf_index'), r.get('visible', True)
            )
        n_vis = sum(1 for r in self._results if r.get('visible'))
        self.status_label.setText(
            f'Target set: ORF at position {target_start} (stop {target_end}). '
            f'Showing {n_vis} in-frame fragment(s); out-of-frame hidden.'
        )

    def _on_target_cleared(self):
        self._target_orf_index = None
        for r in self._results:
            r['visible'] = True
        self.orf_table.clear_frame_filter()
        for r in self._results:
            self.viz_widget.update_visibility(
                r.get('orf_index'), True
            )
        self.status_label.setText('Target cleared. Showing all ORFs.')

    def _open_tuner_dialog(self, orf_index=None):
        if not self._results:
            QMessageBox.information(self, 'No data',
                                    'Run OSTIR first to populate ORF results.')
            return
        # Offer both configured anti-SDs to the tuner so the user can choose
        # which ribosome the optimizer scores against.
        ip = self.input_panel
        asd_options = [(f'Primary · {ip.primary_cmb.currentText()}',
                        ip.primary_asd())]
        sec = ip.secondary_asd()
        if sec and sec != ip.primary_asd():
            asd_options.append((f'Secondary · {ip.secondary_cmb.currentText()}',
                                sec))
        dlg = TunerDialog(parent=self, sequence=self._sequence,
                          results=self._results,
                          selected_orf_index=orf_index,
                          asd_options=asd_options)
        dlg.apply_sequence.connect(self._apply_tuned_sequence)
        dlg.exec()

    def _open_codon_optimizer(self):
        dlg = CodonOptimizerDialog(parent=self, results=self._results)
        dlg.exec()

    def _open_batch_dialog(self):
        dlg = BatchDialog(parent=self, default_asd=self.input_panel.primary_asd())
        dlg.exec()

    # ------------------------------------------------------------------ #
    #  Session save / load                                                 #
    # ------------------------------------------------------------------ #
    def _save_session(self):
        if not self._results and not self._sequence:
            QMessageBox.information(self, 'Nothing to save',
                                    'Run OSTIR first, then save the session.')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save session', 'session.tirex',
            'TIRex session (*.tirex);;All files (*)')
        if not path:
            return
        try:
            ip = self.input_panel
            params = {
                'name': ip.seq_name_edit.text(),
                'aSD': ip.primary_asd(),
                'aSD_secondary': ip.secondary_asd(),
            }
            _session.save_session(
                path, sequence=self._sequence, params=params,
                results=self._results, target_orf_index=self._target_orf_index,
                tir_min=getattr(self.viz_widget, '_tir_min', 0.0),
                tir_max=(0.0 if getattr(self.viz_widget, '_tir_max', float('inf'))
                         == float('inf') else self.viz_widget._tir_max))
            self.status_label.setText(f'Session saved to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Save error', str(exc))

    def _open_session(self):
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open session', '',
            'TIRex session (*.tirex);;All files (*)')
        if not path:
            return
        try:
            data = _session.load_session(path)
        except Exception as exc:
            QMessageBox.critical(self, 'Open error', str(exc))
            return
        self._sequence = data.get('sequence', '')
        self._results = data.get('results', []) or []
        self._target_orf_index = data.get('target_orf_index')
        self.input_panel.seq_edit.setPlainText(self._sequence)
        name = (data.get('params') or {}).get('name')
        if name:
            self.input_panel.seq_name_edit.setText(str(name))
        if self._results:
            self.orf_table.update_results(self._results)
            self.viz_widget.set_data(self._sequence, self._results)
        self.status_label.setText(
            f'Loaded session: {len(self._results)} ORF(s), '
            f'{len(self._sequence)} bp.')

    def _apply_tuned_sequence(self, new_seq: str):
        """Load a tuner-proposed edit into the input panel and re-analyze,
        enabling rapid iterative editing."""
        import re
        new_seq = (new_seq or '').upper()
        self.input_panel.seq_edit.setPlainText(new_seq)

        # Bump a version suffix on the sequence name so iterations are distinct.
        name = self.input_panel.seq_name_edit.text().strip() or 'sequence'
        base = re.sub(r'_v\d+$', '', name)
        self._tuner_iter = getattr(self, '_tuner_iter', 0) + 1
        self.input_panel.seq_name_edit.setText(f'{base}_v{self._tuner_iter}')

        self.status_label.setText(
            f'Applied tuner edit → re-analyzing {len(new_seq)} bp sequence…')
        # Re-run OSTIR via the normal input-panel path (reuses current options).
        self.input_panel._on_run()

    def _open_gel_dialog(self):
        if not self._results:
            QMessageBox.information(self, 'No data',
                                    'Run OSTIR first to populate ORF results.')
            return
        # Pass the current TIR range so the gel respects the same filter
        # as the sequence view / table.
        tir_min = getattr(self.viz_widget, '_tir_min', 0.0)
        tir_max = getattr(self.viz_widget, '_tir_max', float('inf'))
        dlg = SDSPageDialog(parent=self, results=self._results,
                            target_orf_index=self._target_orf_index,
                            tir_min=tir_min, tir_max=tir_max)
        dlg.exec()

    # ------------------------------------------------------------------ #
    #  Export                                                              #
    # ------------------------------------------------------------------ #
    def _export_csv(self):
        if not self._results:
            QMessageBox.information(self, 'No data', 'Run OSTIR first.')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save CSV', 'ostir_results.csv',
            'CSV files (*.csv);;All files (*)'
        )
        if not path:
            return
        try:
            # Flatten results (exclude internal keys)
            skip = {'dna_sequence', 'aa_composition', 'visible',
                    'protein_aa_composition'}
            rows = []
            for r in self._results:
                row = {k: v for k, v in r.items()
                       if k not in skip and not isinstance(v, (dict, list))}
                # Add composition as separate columns
                comp = r.get('protein_aa_composition') or {}
                for aa, pct in comp.items():
                    row[f'aa_{aa}_%'] = pct
                rows.append(row)
            df = pd.DataFrame(rows)
            df.to_csv(path, index=False)
            self.status_label.setText(f'CSV saved to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Export error', str(exc))

    def _export_fasta(self):
        if not self._results:
            QMessageBox.information(self, 'No data', 'Run OSTIR first.')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save FASTA', 'ostir_proteins.fasta',
            'FASTA files (*.fasta *.fa);;All files (*)'
        )
        if not path:
            return
        try:
            lines = []
            for r in self._results:
                if not r.get('visible', True):
                    continue
                idx = r.get('orf_index', 0) + 1
                pos = r.get('start_position', '?')
                codon = r.get('start_codon', '?')
                tir = r.get('expression', '?')
                mw = r.get('protein_molecular_weight', '?')
                aa = r.get('aa_sequence', '')
                header = (f'>ORF_{idx}_pos{pos}_{codon}_TIR{tir}_MW{mw}kDa')
                lines.append(header)
                # Wrap at 60 chars
                for i in range(0, len(aa), 60):
                    lines.append(aa[i:i + 60])
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write('\n'.join(lines) + '\n')
            self.status_label.setText(f'FASTA saved to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Export error', str(exc))

    def _export_figure(self, fmt: str):
        path, _ = QFileDialog.getSaveFileName(
            self, f'Save {fmt.upper()}',
            f'ostir_visualization.{fmt}',
            f'{fmt.upper()} files (*.{fmt});;All files (*)'
        )
        if not path:
            return
        try:
            self.viz_widget.save_figure(path)
            self.status_label.setText(f'Figure saved to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Export error', str(exc))

    # ------------------------------------------------------------------ #
    #  About                                                               #
    # ------------------------------------------------------------------ #
    def _show_about(self):
        QMessageBox.about(
            self, 'About TIRex',
            '<b>TIRex</b> — Translation Initiation Rate Explorer<br><br>'
            'A graphical interface for the <b>Open Source Translation '
            'Initiation Rate</b> (OSTIR) predictor.<br><br>'
            'Prediction engine: <a href="https://github.com/barricklab/ostir">'
            'barricklab/ostir</a><br>'
            'Based on the Salis Lab RBS Calculator.<br><br>'
            'Requires ViennaRNA 2.6.4+',
        )
