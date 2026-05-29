"""
Left-side input panel: sequence entry, FASTA loader, all OSTIR options, and Run button.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit, QPushButton,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QComboBox,
    QGroupBox, QFormLayout, QFileDialog, QSizePolicy, QFrame,
    QInputDialog, QMessageBox,
)
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtGui import QFont, QColor, QPalette

from core.asd_presets import (
    load_presets, save_preset, derive_asd_from_rbs, clean_nt,
    NATIVE_NAME, NATIVE_ASD,
)


class InputPanel(QWidget):
    """
    Emits run_requested(sequence, params_dict) when the user clicks Run.
    """

    run_requested = pyqtSignal(str, dict)

    # ------------------------------------------------------------------ #
    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    # ------------------------------------------------------------------ #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(10)

        # ── Header band ───────────────────────────────────────────── #
        header = QFrame()
        header.setObjectName('BrandHeader')
        header.setStyleSheet(
            '#BrandHeader {'
            '  background: qlineargradient(x1:0, y1:0, x2:1, y2:1,'
            '    stop:0 #0f172a, stop:1 #1e3a8a);'
            '  border-radius: 12px;'
            '}'
        )
        hlay = QVBoxLayout(header)
        hlay.setContentsMargins(14, 12, 14, 12)
        hlay.setSpacing(1)

        title = QLabel('TIRex')
        title.setFont(QFont('Segoe UI', 20, QFont.Weight.Bold))
        title.setStyleSheet('color: #ffffff; letter-spacing: 1px;')
        hlay.addWidget(title)

        subtitle = QLabel('Translation Initiation Rate Explorer')
        subtitle.setFont(QFont('Segoe UI', 8))
        subtitle.setStyleSheet('color: #93c5fd;')
        hlay.addWidget(subtitle)

        root.addWidget(header)

        # ── Sequence input ─────────────────────────────────────────── #
        seq_group = QGroupBox('Input Sequence')
        seq_layout = QVBoxLayout(seq_group)
        seq_layout.setSpacing(4)

        self.seq_edit = QTextEdit()
        self.seq_edit.setPlaceholderText(
            'Paste a DNA/RNA sequence here…\nor load a FASTA file below.'
        )
        self.seq_edit.setFont(QFont('Consolas', 9))
        self.seq_edit.setMinimumHeight(120)
        self.seq_edit.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        seq_layout.addWidget(self.seq_edit)

        btn_row = QHBoxLayout()
        self.load_btn = QPushButton('  Load FASTA…')
        self.load_btn.setObjectName('GhostButton')
        self.load_btn.clicked.connect(self._load_fasta)
        self.clear_btn = QPushButton('Clear')
        self.clear_btn.setObjectName('GhostButton')
        self.clear_btn.clicked.connect(self.seq_edit.clear)
        btn_row.addWidget(self.load_btn)
        btn_row.addWidget(self.clear_btn)
        seq_layout.addLayout(btn_row)

        self.seq_name_edit = QLineEdit()
        self.seq_name_edit.setPlaceholderText('Sequence name (optional)')
        seq_layout.addWidget(self.seq_name_edit)

        root.addWidget(seq_group)

        # ── OSTIR Options ──────────────────────────────────────────── #
        opt_group = QGroupBox('OSTIR Options')
        opt_form = QFormLayout(opt_group)
        opt_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        opt_form.setSpacing(4)

        # Start position
        self.start_spin = QSpinBox()
        self.start_spin.setRange(0, 999_999)
        self.start_spin.setSpecialValueText('Auto')
        self.start_spin.setValue(0)
        self.start_spin.setToolTip(
            "Most 5' position to consider a start codon (1-indexed). 0 = auto."
        )
        opt_form.addRow('Start pos:', self.start_spin)

        # End position
        self.end_spin = QSpinBox()
        self.end_spin.setRange(0, 999_999)
        self.end_spin.setSpecialValueText('Auto')
        self.end_spin.setValue(0)
        self.end_spin.setToolTip(
            "Most 3' position to consider a start codon (1-indexed). 0 = auto."
        )
        opt_form.addRow('End pos:', self.end_spin)

        # Circular
        self.circular_cb = QCheckBox('Treat as circular')
        self.circular_cb.setToolTip('Enable for plasmids or circular sequences.')
        opt_form.addRow('', self.circular_cb)

        # Threads
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 64)
        self.threads_spin.setValue(1)
        self.threads_spin.setToolTip('Number of parallel threads for OSTIR.')
        opt_form.addRow('Threads:', self.threads_spin)

        # Decimal places
        self.decimals_spin = QSpinBox()
        self.decimals_spin.setRange(1, 8)
        self.decimals_spin.setValue(4)
        opt_form.addRow('Decimals:', self.decimals_spin)

        # Verbosity
        self.verbosity_combo = QComboBox()
        self.verbosity_combo.addItems(['0 – Silent', '1 – Normal', '2 – Verbose'])
        self.verbosity_combo.setCurrentIndex(0)
        opt_form.addRow('Verbosity:', self.verbosity_combo)

        # Folding constraints
        self.constraints_edit = QLineEdit()
        self.constraints_edit.setPlaceholderText('Optional')
        self.constraints_edit.setFont(QFont('Consolas', 9))
        self.constraints_edit.setToolTip(
            'ViennaRNA folding constraints string (starting from the '
            '35 bp window upstream of each start codon).'
        )
        opt_form.addRow('Constraints:', self.constraints_edit)

        root.addWidget(opt_group)

        # ── Anti-SD (ribosome) ─────────────────────────────────────── #
        root.addWidget(self._build_asd_group())

        # ── Run button ────────────────────────────────────────────── #
        self.run_btn = QPushButton('▶  Run OSTIR')
        self.run_btn.setObjectName('PrimaryButton')
        self.run_btn.setMinimumHeight(42)
        self.run_btn.clicked.connect(self._on_run)
        root.addWidget(self.run_btn)

        root.addStretch()

    # ------------------------------------------------------------------ #
    #  Anti-SD group                                                       #
    # ------------------------------------------------------------------ #
    def _build_asd_group(self):
        self._presets = load_presets()
        box = QGroupBox('Anti-SD (ribosome)')
        v = QVBoxLayout(box)
        v.setSpacing(4)

        lbl1 = QLabel('Primary (always scored):')
        lbl1.setProperty('muted', True)
        v.addWidget(lbl1)
        self.primary_cmb = QComboBox()
        self.primary_cmb.currentTextChanged.connect(self._on_primary_preset)
        v.addWidget(self.primary_cmb)
        self.primary_edit = QLineEdit()
        self.primary_edit.setFont(QFont('Consolas', 9))
        self.primary_edit.setToolTip(
            "Anti-SD (3' 16S rRNA tail) used for the primary TIR column.")
        v.addWidget(self.primary_edit)

        self.secondary_chk = QCheckBox('Also score with a secondary anti-SD')
        self.secondary_chk.setChecked(True)
        self.secondary_chk.setToolTip(
            'When enabled, every ORF is also scored against this second '
            'anti-SD, giving a side-by-side TIR comparison column.')
        self.secondary_chk.toggled.connect(self._on_secondary_toggle)
        v.addWidget(self.secondary_chk)
        self.secondary_cmb = QComboBox()
        self.secondary_cmb.currentTextChanged.connect(self._on_secondary_preset)
        v.addWidget(self.secondary_cmb)
        self.secondary_edit = QLineEdit()
        self.secondary_edit.setFont(QFont('Consolas', 9))
        v.addWidget(self.secondary_edit)

        brow = QHBoxLayout()
        self.derive_btn = QPushButton('Derive from RBS…')
        self.derive_btn.setObjectName('GhostButton')
        self.derive_btn.clicked.connect(self._derive_from_rbs)
        self.save_asd_btn = QPushButton('Save…')
        self.save_asd_btn.setObjectName('GhostButton')
        self.save_asd_btn.clicked.connect(self._save_secondary_preset)
        brow.addWidget(self.derive_btn)
        brow.addWidget(self.save_asd_btn)
        v.addLayout(brow)

        self._refill_asd_combos(primary=NATIVE_NAME, secondary='pET T7 (derived)')
        self._on_secondary_toggle(self.secondary_chk.isChecked())
        return box

    def _refill_asd_combos(self, primary=None, secondary=None):
        self._presets = load_presets()
        names = list(self._presets.keys())
        if primary is None:
            primary = self.primary_cmb.currentText() or NATIVE_NAME
        if secondary is None:
            secondary = self.secondary_cmb.currentText()
        for cmb in (self.primary_cmb, self.secondary_cmb):
            cmb.blockSignals(True)
            cmb.clear()
            cmb.addItems(names)
            cmb.blockSignals(False)
        if primary in self._presets:
            self.primary_cmb.setCurrentText(primary)
            self.primary_edit.setText(self._presets[primary])
        else:
            self.primary_edit.setText(NATIVE_ASD)
        if secondary in self._presets:
            self.secondary_cmb.setCurrentText(secondary)
            self.secondary_edit.setText(self._presets[secondary])

    def _on_primary_preset(self, name):
        if name in getattr(self, '_presets', {}):
            self.primary_edit.setText(self._presets[name])

    def _on_secondary_preset(self, name):
        if name in getattr(self, '_presets', {}):
            self.secondary_edit.setText(self._presets[name])

    def _on_secondary_toggle(self, on):
        self.secondary_cmb.setEnabled(on)
        self.secondary_edit.setEnabled(on)

    def _derive_from_rbs(self):
        text, ok = QInputDialog.getMultiLineText(
            self, 'Derive anti-SD from RBS',
            "Paste an mRNA RBS / 5' leader — the SD core is detected and its "
            "reverse complement becomes the anti-SD:", '')
        if not ok or not text.strip():
            return
        asd = derive_asd_from_rbs(text)
        if not asd:
            QMessageBox.warning(self, 'Derive failed',
                                'Could not detect a Shine-Dalgarno core.')
            return
        name, ok2 = QInputDialog.getText(
            self, 'Name preset', 'Save derived anti-SD as:',
            text=f'Custom ({asd})')
        if ok2 and name.strip():
            self._presets = save_preset(name.strip(), asd)
            self._refill_asd_combos(secondary=name.strip())
        else:
            self.secondary_edit.setText(asd)
        self.secondary_chk.setChecked(True)

    def _save_secondary_preset(self):
        seq = clean_nt(self.secondary_edit.text())
        if not seq:
            QMessageBox.warning(self, 'Nothing to save',
                                'Enter a secondary anti-SD sequence first.')
            return
        name, ok = QInputDialog.getText(self, 'Save preset', 'Preset name:')
        if ok and name.strip():
            self._presets = save_preset(name.strip(), seq)
            self._refill_asd_combos(secondary=name.strip())

    def primary_asd(self):
        return clean_nt(self.primary_edit.text()) or NATIVE_ASD

    def secondary_asd(self):
        if not self.secondary_chk.isChecked():
            return None
        return clean_nt(self.secondary_edit.text()) or None

    # ------------------------------------------------------------------ #
    def _load_fasta(self):
        from PyQt6.QtWidgets import QMessageBox
        from core.seq_import import load_sequence, SUPPORTED_FILTER
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open sequence file', '', SUPPORTED_FILTER)
        if not path:
            return
        try:
            imp = load_sequence(path)
        except Exception as exc:
            QMessageBox.critical(self, 'Error loading file', str(exc))
            return
        if not imp.sequence:
            QMessageBox.warning(self, 'Empty sequence',
                                'No nucleotide sequence found in that file.')
            return
        self.seq_edit.setPlainText(imp.sequence)
        if not self.seq_name_edit.text():
            self.seq_name_edit.setText(imp.name)
        if imp.cds:
            QMessageBox.information(
                self, 'Sequence loaded',
                f'Loaded "{imp.name}" ({len(imp.sequence)} nt) from '
                f'{imp.source_format.upper()}.\n\n{len(imp.cds)} forward CDS '
                'feature(s) found — run OSTIR to score their start codons.')

    # ------------------------------------------------------------------ #
    def _on_run(self):
        import re as _re
        from PyQt6.QtWidgets import QMessageBox

        raw_text = self.seq_edit.toPlainText().strip()
        if not raw_text:
            QMessageBox.warning(self, 'No input', 'Please enter or load a sequence.')
            return

        # Strip FASTA headers; keep only letters to inspect what was provided.
        body = self._strip_headers(raw_text)
        letters = _re.sub(r'[^A-Za-z]', '', body).upper()
        sequence = _re.sub(r'[^ACGTU]', '', letters).replace('U', 'T')
        if not sequence:
            QMessageBox.warning(self, 'Invalid sequence',
                                'Could not extract a valid nucleotide sequence.')
            return

        # Warn (but proceed) if ambiguity/invalid bases had to be removed —
        # OSTIR only accepts A/C/G/T/U and would otherwise crash.
        removed = sorted(set(letters) - set('ACGTU'))
        if removed:
            QMessageBox.information(
                self, 'Non-standard bases removed',
                'These characters are not A/C/G/T/U and were removed before '
                f'analysis: {", ".join(removed)}.\n\n'
                'Note this shifts downstream positions; edit the input if that '
                'matters.')

        # Validate anti-SDs (OSTIR requires exactly 9 A/C/G/T bases).
        primary = self.primary_asd()
        if not _re.fullmatch(r'[ACGT]{9}', primary or ''):
            QMessageBox.warning(
                self, 'Invalid primary anti-SD',
                f'The primary anti-SD must be exactly 9 bases (A/C/G/T).\n\n'
                f'Got: "{primary}" ({len(primary or "")} nt). '
                'Fix it before running.')
            return
        secondary = self.secondary_asd()
        if secondary is not None and not _re.fullmatch(r'[ACGT]{9}', secondary):
            QMessageBox.warning(
                self, 'Invalid secondary anti-SD',
                f'The secondary anti-SD must be exactly 9 bases (A/C/G/T); '
                f'got "{secondary}". It will be skipped for this run.')
            secondary = None

        params = {
            'name': self.seq_name_edit.text() or 'sequence',
            'start': self.start_spin.value() or None,
            'end': self.end_spin.value() or None,
            'aSD': primary,
            'aSD_secondary': secondary,
            'circular': self.circular_cb.isChecked(),
            'threads': self.threads_spin.value(),
            'decimal_places': self.decimals_spin.value(),
            'verbosity': self.verbosity_combo.currentIndex(),
            'constraints': self.constraints_edit.text().strip() or None,
        }

        self.run_requested.emit(sequence, params)

    # ------------------------------------------------------------------ #
    @staticmethod
    def _strip_headers(text: str) -> str:
        """Drop FASTA '>' header lines, join the rest."""
        lines = text.strip().splitlines()
        return ''.join(l for l in lines if not l.startswith('>'))

    @staticmethod
    def _parse_input(text: str) -> str:
        """Accept a raw sequence string or FASTA. Returns a clean A/C/G/T
        sequence (U->T, all non-nucleotide characters removed)."""
        import re as _re
        body = InputPanel._strip_headers(text)
        return _re.sub(r'[^ACGTU]', '', body.upper().replace('U', 'T'))

    # ------------------------------------------------------------------ #
    def set_running(self, running: bool):
        """Disable/enable the run button and input while OSTIR is running."""
        self.run_btn.setEnabled(not running)
        self.run_btn.setText('⏳  Running…' if running else '▶  Run OSTIR')
