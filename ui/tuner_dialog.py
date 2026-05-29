"""
TIR-Tuner dialog for TIRex.

Seeded directly from the ORFs OSTIR already found, this dialog proposes
nucleotide edits that either
  * INCREASE translation initiation at a chosen target start codon, or
  * DECREASE it at a chosen downstream start codon while keeping the main
    protein's amino-acid sequence unchanged.

Results appear as a ranked table and a ΔTIR substitution heatmap (with
user-defined region labels).
"""

from PyQt6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
    QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox, QRadioButton,
    QButtonGroup, QGroupBox, QPushButton, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QProgressBar, QMessageBox, QSplitter,
    QStackedWidget, QTabWidget, QFrame, QListWidget, QTextEdit, QApplication,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont, QColor, QBrush, QGuiApplication

from .tuner_worker import TunerWorker
from .result_plot import ResultPlotWidget
from core.tuner.scoring import OstirScorer


class TunerDialog(QDialog):

    # Emitted when the user checks "Apply" on a candidate edit. Carries the
    # full mutated sequence so the main window can load + re-analyze it.
    apply_sequence = pyqtSignal(str)

    def __init__(self, parent=None, sequence='', results=None,
                 selected_orf_index=None, aSD=None, asd_options=None):
        super().__init__(parent)
        self.setWindowTitle('TIR-Tuner — translation-rate edit predictor')
        self.resize(1240, 820)

        self._seq = (sequence or '').upper().replace('U', 'T')
        self._original_seq = self._seq          # never changes (for diff)
        self._results = list(results or [])
        # Anti-SD options the optimizer can score against.
        if asd_options:
            self._asd_options = [(str(n), s) for n, s in asd_options if s]
        else:
            self._asd_options = [('anti-SD', aSD)] if aSD else [
                ('E. coli (native)', 'ACCTCCTTA')]
        self._aSD = self._asd_options[0][1]
        self._worker = None
        self._candidates = []
        self._mode = 'increase'
        self._populating = False
        self._applied = []                      # [{desc, before, after}]

        # Build the (position, codon, TIR) list from existing ORFs.
        self._starts = sorted(
            [{'pos': int(r.get('start_position', 0)),
              'codon': r.get('start_codon', '?'),
              'tir': float(r.get('expression') or 0.0),
              'orf_index': r.get('orf_index')}
             for r in self._results if r.get('start_position')],
            key=lambda s: s['pos'],
        )
        self._selected_pos = None
        if selected_orf_index is not None:
            for s in self._starts:
                if s['orf_index'] == selected_orf_index:
                    self._selected_pos = s['pos']
                    break

        self._build_ui()
        self._populate_start_combos()

    # ------------------------------------------------------------------ #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        # Header strip
        head = QFrame()
        head.setObjectName('TunerHead')
        head.setStyleSheet(
            '#TunerHead { background: qlineargradient(x1:0,y1:0,x2:1,y2:0,'
            'stop:0 #0f172a, stop:1 #0d9488); border-radius: 10px; }')
        hl = QVBoxLayout(head)
        hl.setContentsMargins(14, 9, 14, 9)
        t = QLabel('TIR-Tuner')
        t.setFont(QFont('Segoe UI', 14, QFont.Weight.Bold))
        t.setStyleSheet('color:white;')
        hl.addWidget(t)
        self._sub_lbl = QLabel()
        self._sub_lbl.setStyleSheet('color:#99f6e4; font-size:8pt;')
        hl.addWidget(self._sub_lbl)
        root.addWidget(head)
        self._update_header()

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        root.addWidget(split, 1)

        # ---------- Left: controls ----------
        left = QWidget()
        left.setMaximumWidth(390)
        left.setMinimumWidth(330)
        ll = QVBoxLayout(left)
        ll.setSpacing(10)

        # Scoring anti-SD selector
        asd_box = QGroupBox('Scoring anti-SD')
        al = QVBoxLayout(asd_box)
        self.asd_cmb = QComboBox()
        for name, seq in self._asd_options:
            self.asd_cmb.addItem(f'{name}  ({seq})', seq)
        self.asd_cmb.setToolTip(
            'Which anti-SD (ribosome) the optimizer scores edits against. '
            'Switching re-scans the current sequence.')
        self.asd_cmb.currentIndexChanged.connect(self._on_asd_changed)
        al.addWidget(self.asd_cmb)
        ll.addWidget(asd_box)

        mode_box = QGroupBox('Mode')
        ml = QVBoxLayout(mode_box)
        self.rb_increase = QRadioButton('Increase initiation at a target start')
        self.rb_decrease = QRadioButton('Decrease a downstream start (protein-preserving)')
        self.rb_library = QRadioButton('Build RBS library / TIR ramp')
        self.rb_increase.setChecked(True)
        g = QButtonGroup(self)
        g.addButton(self.rb_increase)
        g.addButton(self.rb_decrease)
        g.addButton(self.rb_library)
        self.rb_increase.toggled.connect(self._on_mode_changed)
        self.rb_decrease.toggled.connect(self._on_mode_changed)
        self.rb_library.toggled.connect(self._on_mode_changed)
        ml.addWidget(self.rb_increase)
        ml.addWidget(self.rb_decrease)
        ml.addWidget(self.rb_library)
        ll.addWidget(mode_box)

        self.opt_stack = QStackedWidget()
        self.opt_stack.addWidget(self._build_increase_opts())
        self.opt_stack.addWidget(self._build_decrease_opts())
        self.opt_stack.addWidget(self._build_library_opts())
        ll.addWidget(self.opt_stack)

        self.predict_btn = QPushButton('Predict edits  ▶')
        self.predict_btn.setObjectName('AccentButton')
        self.predict_btn.setMinimumHeight(40)
        self.predict_btn.clicked.connect(self._predict)
        ll.addWidget(self.predict_btn)

        # Applied-edits history (this tuning session)
        hist_box = QGroupBox('Applied edits (this session)')
        hv = QVBoxLayout(hist_box)
        self.applied_list = QListWidget()
        self.applied_list.setFixedHeight(96)
        hv.addWidget(self.applied_list)
        hrow = QHBoxLayout()
        self.undo_btn = QPushButton('↩ Undo last')
        self.undo_btn.setObjectName('GhostButton')
        self.undo_btn.clicked.connect(self._undo_last)
        self.compare_btn = QPushButton('Compare…')
        self.compare_btn.setObjectName('GhostButton')
        self.compare_btn.clicked.connect(self._show_diff)
        self.report_btn = QPushButton('Export report…')
        self.report_btn.setObjectName('GhostButton')
        self.report_btn.clicked.connect(self._export_report)
        hrow.addWidget(self.undo_btn)
        hrow.addWidget(self.compare_btn)
        hrow.addWidget(self.report_btn)
        hv.addLayout(hrow)
        ll.addWidget(hist_box)

        ll.addStretch()
        split.addWidget(left)

        # ---------- Right: results ----------
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        self.result_title = QLabel('Results')
        self.result_title.setFont(QFont('Segoe UI', 11, QFont.Weight.Bold))
        rl.addWidget(self.result_title)

        self.tabs = QTabWidget()
        tbl_tab = QWidget()
        tl = QVBoxLayout(tbl_tab)
        tl.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.cellClicked.connect(self._on_cell_clicked)
        self.table.itemChanged.connect(self._on_item_changed)
        tl.addWidget(self.table)
        hint = QLabel('Check “Apply” to load that edit back into TIRex and re-analyze · '
                      'click “Copy seq” to copy the full mutated sequence.')
        hint.setProperty('muted', True)
        hint.setWordWrap(True)
        tl.addWidget(hint)
        self.tabs.addTab(tbl_tab, 'Table')

        self.plot_widget = ResultPlotWidget()
        self.tabs.addTab(self.plot_widget, 'Graph')
        rl.addWidget(self.tabs)
        split.addWidget(right)
        split.setSizes([360, 860])

        # ---------- Bottom: progress + close ----------
        bottom = QHBoxLayout()
        self.status_lbl = QLabel('Pick a mode and predict edits.')
        self.status_lbl.setProperty('muted', True)
        bottom.addWidget(self.status_lbl, 1)
        self.progress = QProgressBar()
        self.progress.setFixedWidth(220)
        self.progress.setVisible(False)
        bottom.addWidget(self.progress)
        self.apply_tirex_btn = QPushButton('✓ Apply to TIRex & close')
        self.apply_tirex_btn.setObjectName('PrimaryButton')
        self.apply_tirex_btn.clicked.connect(self._apply_to_tirex)
        bottom.addWidget(self.apply_tirex_btn)
        close_btn = QPushButton('Close')
        close_btn.setObjectName('GhostButton')
        close_btn.clicked.connect(self.reject)
        bottom.addWidget(close_btn)
        root.addLayout(bottom)

        self._on_mode_changed()
        self._update_applied_label()

    # ------------------------------------------------------------------ #
    def _build_increase_opts(self):
        box = QGroupBox('Increase options')
        g = QGridLayout(box)
        g.addWidget(QLabel('Target start:'), 0, 0)
        self.inc_target_cmb = QComboBox()
        g.addWidget(self.inc_target_cmb, 0, 1)

        g.addWidget(QLabel('Upstream window (nt):'), 1, 0)
        self.inc_window_spin = QSpinBox()
        self.inc_window_spin.setRange(2, 40)
        self.inc_window_spin.setValue(20)
        g.addWidget(self.inc_window_spin, 1, 1)

        self.inc_dinuc_chk = QCheckBox('Include dinucleotide substitutions')
        self.inc_dinuc_chk.setChecked(True)
        g.addWidget(self.inc_dinuc_chk, 2, 0, 1, 2)

        g.addWidget(QLabel('Search depth:'), 3, 0)
        drow = QHBoxLayout()
        self.inc_enum_rb = QRadioButton('Single+dinuc')
        self.inc_greedy_rb = QRadioButton('Greedy')
        self.inc_beam_rb = QRadioButton('Beam')
        self.inc_enum_rb.setChecked(True)
        dg = QButtonGroup(self)
        dg.addButton(self.inc_enum_rb)
        dg.addButton(self.inc_greedy_rb)
        dg.addButton(self.inc_beam_rb)
        drow.addWidget(self.inc_enum_rb)
        drow.addWidget(self.inc_greedy_rb)
        drow.addWidget(self.inc_beam_rb)
        dw = QWidget(); dw.setLayout(drow)
        g.addWidget(dw, 3, 1)

        g.addWidget(QLabel('Target fold (greedy/beam):'), 4, 0)
        self.inc_fold_spin = QDoubleSpinBox()
        self.inc_fold_spin.setRange(1.0, 1000.0)
        self.inc_fold_spin.setValue(10.0)
        g.addWidget(self.inc_fold_spin, 4, 1)

        g.addWidget(QLabel('Max rounds (greedy/beam):'), 5, 0)
        self.inc_rounds_spin = QSpinBox()
        self.inc_rounds_spin.setRange(1, 20)
        self.inc_rounds_spin.setValue(6)
        g.addWidget(self.inc_rounds_spin, 5, 1)

        g.addWidget(QLabel('Beam width:'), 6, 0)
        self.inc_beam_spin = QSpinBox()
        self.inc_beam_spin.setRange(2, 20)
        self.inc_beam_spin.setValue(5)
        g.addWidget(self.inc_beam_spin, 6, 1)
        return box

    def _build_library_opts(self):
        box = QGroupBox('RBS library / TIR ramp')
        g = QGridLayout(box)
        g.addWidget(QLabel('Target start:'), 0, 0)
        self.lib_target_cmb = QComboBox()
        g.addWidget(self.lib_target_cmb, 0, 1)

        g.addWidget(QLabel('Upstream window (nt):'), 1, 0)
        self.lib_window_spin = QSpinBox()
        self.lib_window_spin.setRange(2, 40)
        self.lib_window_spin.setValue(20)
        g.addWidget(self.lib_window_spin, 1, 1)

        self.lib_dinuc_chk = QCheckBox('Include dinucleotide substitutions')
        self.lib_dinuc_chk.setChecked(True)
        g.addWidget(self.lib_dinuc_chk, 2, 0, 1, 2)

        g.addWidget(QLabel('Number of variants:'), 3, 0)
        self.lib_n_spin = QSpinBox()
        self.lib_n_spin.setRange(2, 24)
        self.lib_n_spin.setValue(8)
        g.addWidget(self.lib_n_spin, 3, 1)

        note = QLabel('Picks variants at log-spaced TIR levels across the '
                      'achievable range (weak → strong), incl. wild-type.')
        note.setProperty('muted', True)
        note.setWordWrap(True)
        g.addWidget(note, 4, 0, 1, 2)
        return box

    def _build_decrease_opts(self):
        box = QGroupBox('Decrease options (protein preserved)')
        g = QGridLayout(box)
        g.addWidget(QLabel('Main start (keep):'), 0, 0)
        self.dec_main_cmb = QComboBox()
        self.dec_main_cmb.currentIndexChanged.connect(self._refresh_down_combo)
        g.addWidget(self.dec_main_cmb, 0, 1)

        self.dec_inframe_chk = QCheckBox('Only downstream starts in frame with main')
        self.dec_inframe_chk.setChecked(True)
        self.dec_inframe_chk.setToolTip(
            'When ticked, only list downstream start codons whose position is a '
            'multiple of 3 from the main start — i.e. internal in-frame starts '
            'of the same reading frame.')
        self.dec_inframe_chk.toggled.connect(self._refresh_down_combo)
        g.addWidget(self.dec_inframe_chk, 1, 0, 1, 2)

        g.addWidget(QLabel('Downstream start (suppress):'), 2, 0)
        self.dec_down_cmb = QComboBox()
        g.addWidget(self.dec_down_cmb, 2, 1)

        g.addWidget(QLabel('Strategy:'), 3, 0)
        self.dec_strategy_cmb = QComboBox()
        self.dec_strategy_cmb.addItems(['both', 'codon swap', 'RBS synonymous'])
        g.addWidget(self.dec_strategy_cmb, 3, 1)

        g.addWidget(QLabel('RBS window (nt):'), 4, 0)
        self.dec_rbs_spin = QSpinBox()
        self.dec_rbs_spin.setRange(3, 40)
        self.dec_rbs_spin.setValue(20)
        g.addWidget(self.dec_rbs_spin, 4, 1)
        return box

    # ------------------------------------------------------------------ #
    def _start_label(self, s):
        return f"pos {s['pos']} · {s['codon']} · TIR {s['tir']:.1f}"

    def _populate_start_combos(self):
        for cmb in (self.inc_target_cmb, self.dec_main_cmb, self.lib_target_cmb):
            cmb.blockSignals(True)
            cmb.clear()
            for i, s in enumerate(self._starts):
                cmb.addItem(self._start_label(s), i)   # store array index
            cmb.blockSignals(False)

        if not self._starts:
            self.predict_btn.setEnabled(False)
            self.dec_down_cmb.clear()
            self.status_lbl.setText('No start codons in the sequence.')
            return

        self.predict_btn.setEnabled(True)
        strongest = max(range(len(self._starts)),
                        key=lambda i: self._starts[i]['tir'])
        sel = strongest
        if self._selected_pos is not None:
            for i, s in enumerate(self._starts):
                if s['pos'] == self._selected_pos:
                    sel = i
                    break
        self.inc_target_cmb.setCurrentIndex(sel)
        self.dec_main_cmb.setCurrentIndex(sel)
        self.lib_target_cmb.setCurrentIndex(sel)
        self._refresh_down_combo()

    def _refresh_down_combo(self):
        """Rebuild the downstream-start combo from the current main start,
        optionally restricted to starts in frame with it."""
        if not hasattr(self, 'dec_down_cmb'):
            return
        self.dec_down_cmb.blockSignals(True)
        self.dec_down_cmb.clear()
        main_idx = self.dec_main_cmb.currentData()
        if main_idx is not None and 0 <= main_idx < len(self._starts):
            main_pos = self._starts[main_idx]['pos']
            only_frame = self.dec_inframe_chk.isChecked()
            for i, s in enumerate(self._starts):
                if s['pos'] <= main_pos:
                    continue
                if only_frame and (s['pos'] - main_pos) % 3 != 0:
                    continue
                self.dec_down_cmb.addItem(self._start_label(s), i)
        self.dec_down_cmb.blockSignals(False)
        if self.dec_down_cmb.count() == 0:
            self.dec_down_cmb.addItem('— none available —', None)

    def _on_mode_changed(self):
        if self.rb_increase.isChecked():
            self.opt_stack.setCurrentIndex(0)
            self._mode = 'increase'
        elif self.rb_decrease.isChecked():
            self.opt_stack.setCurrentIndex(1)
            self._mode = 'decrease'
        else:
            self.opt_stack.setCurrentIndex(2)
            self._mode = 'library'

    def _current_asd(self):
        return self.asd_cmb.currentData() or self._aSD

    def _on_asd_changed(self):
        """User switched which anti-SD the optimizer scores against → re-scan
        the current sequence so start-codon TIRs reflect the new ribosome."""
        self._aSD = self._current_asd()
        self.status_lbl.setText('Re-scanning with selected anti-SD…')
        QApplication.processEvents()
        self._rescan_starts()
        self._populate_start_combos()
        self._clear_results()
        self._update_header()
        self.status_lbl.setText(
            f'Scoring anti-SD set to {self.asd_cmb.currentText()}.')

    # ------------------------------------------------------------------ #
    def _predict(self):
        if not self._starts:
            return
        if self.rb_increase.isChecked():
            mode = 'increase'
            ti = self.inc_target_cmb.currentData()
            if ti is None:
                return
            if self.inc_greedy_rb.isChecked():
                depth = 'greedy'
            elif self.inc_beam_rb.isChecked():
                depth = 'beam'
            else:
                depth = 'enumerate'
            params = {
                'aSD': self._current_asd(),
                'target_start': self._starts[ti]['pos'],
                'upstream_len': self.inc_window_spin.value(),
                'include_dinucleotide': self.inc_dinuc_chk.isChecked(),
                'depth': depth,
                'target_fold': self.inc_fold_spin.value(),
                'max_rounds': self.inc_rounds_spin.value(),
                'beam_width': self.inc_beam_spin.value(),
            }
        elif self.rb_library.isChecked():
            mode = 'library'
            ti = self.lib_target_cmb.currentData()
            if ti is None:
                return
            params = {
                'aSD': self._current_asd(),
                'target_start': self._starts[ti]['pos'],
                'upstream_len': self.lib_window_spin.value(),
                'include_dinucleotide': self.lib_dinuc_chk.isChecked(),
                'n_variants': self.lib_n_spin.value(),
            }
        else:
            mode = 'decrease'
            mi, di = self.dec_main_cmb.currentData(), self.dec_down_cmb.currentData()
            if mi is None or di is None:
                QMessageBox.warning(
                    self, 'No downstream start',
                    'No downstream start codon is available for the current '
                    'selection. Try unticking “Only in frame with main”.')
                return
            if mi == di:
                QMessageBox.warning(self, 'Pick two starts',
                                    'Main and downstream start must differ.')
                return
            main = self._starts[mi]['pos']
            down = self._starts[di]['pos']
            if down <= main:
                QMessageBox.warning(self, 'Order',
                                    'The downstream start must be after the main start.')
                return
            params = {
                'aSD': self._current_asd(), 'main_start': main, 'down_start': down,
                'strategy': self.dec_strategy_cmb.currentText(),
                'rbs_window': self.dec_rbs_spin.value(),
            }

        self._mode = mode
        self.predict_btn.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.status_lbl.setText('Searching for edits…')

        self._worker = TunerWorker(self._seq, mode, params, parent=self)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished_ok.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _on_progress(self, done, total, msg):
        if total:
            self.progress.setRange(0, total)
            self.progress.setValue(done)
        self.status_lbl.setText(f'{msg}  ({done}/{total})' if total else msg)

    def _on_failed(self, msg):
        self.progress.setVisible(False)
        self.predict_btn.setEnabled(True)
        QMessageBox.critical(self, 'Search error', msg)
        self.status_lbl.setText('Search failed.')

    def _on_done(self, candidates, stats):
        self.progress.setVisible(False)
        self.predict_btn.setEnabled(True)
        self._candidates = candidates
        self._fill_table(candidates)
        self.plot_widget.set_candidates(candidates, self._seq, self._mode)
        self.status_lbl.setText(
            f"{stats['candidates']} candidate edit(s) · "
            f"{stats['ostir_calls']} OSTIR folds · {stats['cache_hits']} cache hits.")

    # ------------------------------------------------------------------ #
    def _fill_table(self, candidates):
        self._populating = True
        decrease = (self._mode == 'decrease')
        cols = (['Apply', '#', 'Edit', 'Pos', 'Change', 'Pred TIR', 'ΔTIR', 'Fold']
                + (['Main ΔTIR'] if decrease else [])
                + ['Warnings', 'Notes', 'Copy seq'])
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.setRowCount(len(candidates))
        self._copy_col = len(cols) - 1

        for row, c in enumerate(candidates):
            vals = [
                '',  # Apply checkbox placeholder
                str(row + 1), c.edit_type, str(c.position), c.change,
                f'{c.tir:.2f}', f'{c.delta:+.2f}',
                ('∞' if c.fold == float('inf') else f'{c.fold:.2f}×'),
            ]
            if decrease:
                md = c.main_delta
                vals.append('—' if md is None else f'{md:+.2f}')
            vals += [c.warning_text or '—', c.notes, '\U0001F4CB Copy']
            for col, v in enumerate(vals):
                if col == 0:
                    item = QTableWidgetItem()
                    item.setFlags(Qt.ItemFlag.ItemIsEnabled
                                  | Qt.ItemFlag.ItemIsUserCheckable
                                  | Qt.ItemFlag.ItemIsSelectable)
                    item.setCheckState(Qt.CheckState.Unchecked)
                    item.setData(Qt.ItemDataRole.UserRole, c.new_sequence)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    item.setToolTip('Apply this edit and re-analyze in TIRex')
                    self.table.setItem(row, col, item)
                    continue
                item = QTableWidgetItem(v)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                if col == self._copy_col:
                    item.setForeground(QBrush(QColor('#2563eb')))
                    item.setData(Qt.ItemDataRole.UserRole, c.new_sequence)
                if cols[col] == 'ΔTIR':
                    item.setForeground(QBrush(QColor('#16a34a') if c.delta >= 0
                                              else QColor('#e11d48')))
                if cols[col] == 'Warnings' and c.warnings:
                    item.setForeground(QBrush(QColor('#b45309')))
                    item.setToolTip('\n'.join(c.warnings))
                self.table.setItem(row, col, item)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents)
        baseline = candidates[0].baseline_tir if candidates else 0.0
        self.result_title.setText(
            f'Results — wild-type TIR at scored start = {baseline:.2f}')
        self._populating = False

    # ------------------------------------------------------------------ #
    def _on_item_changed(self, item):
        """A row's Apply box was ticked → apply the edit in place (the dialog
        stays open and updates with the new sequence)."""
        if self._populating or item.column() != 0:
            return
        if item.checkState() != Qt.CheckState.Checked:
            return
        row = item.row()
        if 0 <= row < len(self._candidates):
            self._apply_edit_in_place(self._candidates[row])

    def _on_cell_clicked(self, row, col):
        if col != getattr(self, '_copy_col', -1):
            return
        item = self.table.item(row, col)
        seq = item.data(Qt.ItemDataRole.UserRole) if item else None
        if seq:
            QGuiApplication.clipboard().setText(seq)
            self.status_lbl.setText(
                f'Copied mutated sequence ({len(seq)} nt) to clipboard.')

    # ================================================================== #
    #  Iterative editing
    # ================================================================== #
    def _apply_edit_in_place(self, c):
        """Adopt candidate `c`'s sequence as the new working sequence, log it,
        re-scan start codons, and reset the results for the next round."""
        self._applied.append({
            'desc': f'{c.edit_type}: {c.change} @ pos {c.position}',
            'before': self._seq,
            'after': c.new_sequence,
            'tir': c.tir,
            'baseline': c.baseline_tir,
        })
        self._seq = c.new_sequence
        self._selected_pos = None
        self.status_lbl.setText('Applied edit · re-scanning start codons…')
        QApplication.processEvents()
        self._rescan_starts()
        self._populate_start_combos()
        self._clear_results()
        self._update_header()
        self._update_applied_label()
        self.status_lbl.setText(
            f'Applied: {self._applied[-1]["desc"]}. Sequence now '
            f'{len(self._seq)} bp — predict again to stack the next change.')

    def _rescan_starts(self):
        """Re-run the OSTIR start-codon scan on the current working sequence."""
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            scorer = OstirScorer(aSD=self._aSD)
            raw = scorer.find_starts(self._seq)
        except Exception:
            raw = []
        finally:
            QApplication.restoreOverrideCursor()
        self._starts = sorted(
            [{'pos': int(r.get('start_position', 0)),
              'codon': r.get('start_codon', '?'),
              'tir': float(r.get('expression') or 0.0),
              'orf_index': None}
             for r in raw if r.get('start_position')],
            key=lambda s: s['pos'])

    def _clear_results(self):
        self._populating = True
        self.table.setRowCount(0)
        self._populating = False
        self._candidates = []
        self.plot_widget.set_candidates([], self._seq, self._mode)
        self.result_title.setText('Results')

    def _update_header(self):
        extra = (f' · {len(self._applied)} edit(s) applied'
                 if self._applied else '')
        self._sub_lbl.setText(
            f'{len(self._seq)} bp sequence · {len(self._starts)} start '
            f'codon(s) available · scored with OSTIR{extra}')

    def _update_applied_label(self):
        self.applied_list.clear()
        for i, a in enumerate(self._applied, 1):
            self.applied_list.addItem(f'{i}. {a["desc"]}')
        has = bool(self._applied)
        self.undo_btn.setEnabled(has)
        self.compare_btn.setEnabled(has)
        self.report_btn.setEnabled(has)
        self.apply_tirex_btn.setEnabled(has)

    def _undo_last(self):
        if not self._applied:
            return
        last = self._applied.pop()
        self._seq = last['before']
        self._selected_pos = None
        self.status_lbl.setText('Undoing last edit · re-scanning…')
        QApplication.processEvents()
        self._rescan_starts()
        self._populate_start_combos()
        self._clear_results()
        self._update_header()
        self._update_applied_label()
        self.status_lbl.setText(
            f'Reverted: {last["desc"]}. Sequence now {len(self._seq)} bp.')

    def _apply_to_tirex(self):
        """Send the final tuned sequence to the main window and close."""
        if not self._applied:
            QMessageBox.information(
                self, 'No edits applied',
                'Tick an “Apply” box on an edit first, then send the result '
                'to TIRex.')
            return
        self.apply_sequence.emit(self._seq)
        self.accept()

    # ================================================================== #
    #  Old-vs-new comparison
    # ================================================================== #
    def _show_diff(self):
        if not self._applied:
            return
        dlg = _SequenceDiffDialog(self._original_seq, self._seq, parent=self)
        dlg.exec()

    def _export_report(self):
        if not self._applied:
            QMessageBox.information(self, 'No edits',
                                    'Apply at least one edit first.')
            return
        from PyQt6.QtWidgets import QFileDialog
        from core.primers import design_primers
        from core.report import build_report
        import io

        path, _ = QFileDialog.getSaveFileName(
            self, 'Export tuning report', 'tirex_tuning_report.html',
            'HTML report (*.html);;All files (*)')
        if not path:
            return
        try:
            # Per-edit Q5 primers.
            primers = []
            for a in self._applied:
                p = design_primers(a['before'], a['after'])
                if p:
                    d = {'edit': a['desc']}
                    d.update(p.as_dict())
                    primers.append(d)
            # Embed current heatmap if present.
            heat = None
            if self._candidates:
                buf = io.BytesIO()
                try:
                    self.plot_widget.figure.savefig(buf, format='png', dpi=110,
                                                    bbox_inches='tight')
                    heat = buf.getvalue()
                except Exception:
                    heat = None
            baseline = self._applied[0].get('baseline')
            final_tir = self._applied[-1].get('tir')
            html = build_report(
                title='TIRex tuning report',
                original=self._original_seq, final=self._seq,
                applied=self._applied, primers=primers or None,
                baseline_tir=baseline, final_tir=final_tir,
                heatmap_png=heat,
                meta={'mode': self._mode,
                      'scoring anti-SD': self.asd_cmb.currentText(),
                      'edits applied': len(self._applied)})
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write(html)
            self.status_lbl.setText(f'Report written to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Report error', str(exc))


# ====================================================================== #
#  Old-vs-new sequence comparison dialog
# ====================================================================== #
class _SequenceDiffDialog(QDialog):
    """Side-by-side, position-aware comparison of the original vs the tuned
    sequence, with substitutions / insertions / deletions colour-highlighted."""

    _COLORS = {
        'equal':   None,
        'replace': '#fde68a',   # amber  — substitution
        'insert':  '#bbf7d0',   # green  — inserted base (in new)
        'delete':  '#fecaca',   # red    — deleted base (was in old)
    }

    def __init__(self, original, current, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Original vs. tuned sequence')
        self.resize(900, 620)
        self._a = original
        self._b = current
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 12, 12, 12)
        lay.setSpacing(8)

        legend = QLabel(
            'Legend:  '
            '<span style="background:#fde68a;">&nbsp;substitution&nbsp;</span>&nbsp;&nbsp;'
            '<span style="background:#bbf7d0;">&nbsp;insertion&nbsp;</span>&nbsp;&nbsp;'
            '<span style="background:#fecaca;">&nbsp;deletion&nbsp;</span>&nbsp;&nbsp;'
            '· gaps shown as <b>·</b>')
        legend.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(legend)

        self.summary = QLabel()
        self.summary.setProperty('muted', True)
        lay.addWidget(self.summary)

        self.view = QTextEdit()
        self.view.setReadOnly(True)
        self.view.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self.view.setFont(QFont('Consolas', 11))
        lay.addWidget(self.view, 1)

        row = QHBoxLayout()
        row.addStretch()
        copy_btn = QPushButton('Copy tuned sequence')
        copy_btn.setObjectName('GhostButton')
        copy_btn.clicked.connect(
            lambda: (QGuiApplication.clipboard().setText(self._b)))
        row.addWidget(copy_btn)
        close_btn = QPushButton('Close')
        close_btn.setObjectName('GhostButton')
        close_btn.clicked.connect(self.accept)
        row.addWidget(close_btn)
        lay.addLayout(row)

        self._render()

    def _render(self):
        import difflib
        a, b = self._a, self._b
        sm = difflib.SequenceMatcher(None, a, b, autojunk=False)

        ao, bo, tags = [], [], []
        n_sub = n_ins = n_del = 0
        for tag, i1, i2, j1, j2 in sm.get_opcodes():
            if tag == 'equal':
                for k in range(i2 - i1):
                    ao.append(a[i1 + k]); bo.append(b[j1 + k]); tags.append('equal')
            elif tag == 'replace':
                la, lb = i2 - i1, j2 - j1
                n_sub += max(la, lb)
                for k in range(max(la, lb)):
                    ao.append(a[i1 + k] if k < la else '-')
                    bo.append(b[j1 + k] if k < lb else '-')
                    tags.append('replace')
            elif tag == 'delete':
                n_del += (i2 - i1)
                for k in range(i2 - i1):
                    ao.append(a[i1 + k]); bo.append('-'); tags.append('delete')
            elif tag == 'insert':
                n_ins += (j2 - j1)
                for k in range(j2 - j1):
                    ao.append('-'); bo.append(b[j1 + k]); tags.append('insert')

        # running 1-based positions per aligned column
        oidx, midx = [], []
        oi = mi = 0
        for k in range(len(tags)):
            if ao[k] != '-': oi += 1
            if bo[k] != '-': mi += 1
            oidx.append(oi); midx.append(mi)

        def cell(ch, tag):
            color = self._COLORS[tag]
            disp = ch if ch != '-' else '·'
            if color:
                return f'<span style="background:{color};">{disp}</span>'
            return disp

        width = 60
        blocks = []
        for start in range(0, len(tags), width):
            end = min(start + width, len(tags))
            o_pos = (oidx[start - 1] if start > 0 else 0) + 1
            m_pos = (midx[start - 1] if start > 0 else 0) + 1
            orig_line = ''.join(cell(ao[k], tags[k]) for k in range(start, end))
            new_line = ''.join(cell(bo[k], tags[k]) for k in range(start, end))
            blocks.append(
                f'<div style="color:#94a3b8;">old {o_pos:>5}</div>'
                f'<div>{orig_line}</div>'
                f'<div>{new_line}</div>'
                f'<div style="color:#94a3b8;">new {m_pos:>5}</div>'
                f'<div>&nbsp;</div>')

        html = ('<div style="font-family:Consolas,monospace; '
                'white-space:pre; line-height:135%;">'
                + ''.join(blocks) + '</div>')
        self.view.setHtml(html)

        net = len(b) - len(a)
        self.summary.setText(
            f'{len(a)} bp → {len(b)} bp  (Δ {net:+d})   ·   '
            f'{n_sub} substitution(s), {n_ins} insertion(s), {n_del} deletion(s)')
