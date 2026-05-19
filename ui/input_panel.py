"""
Left-side input panel: sequence entry, FASTA loader, all OSTIR options, and Run button.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTextEdit, QPushButton,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QComboBox,
    QGroupBox, QFormLayout, QFileDialog, QSizePolicy, QFrame,
)
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtGui import QFont, QColor, QPalette


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
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # ── Title ──────────────────────────────────────────────────── #
        title = QLabel('TIRex')
        title.setFont(QFont('Segoe UI', 13, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(title)

        subtitle = QLabel('Translation Initiation Rate Predictor')
        subtitle.setFont(QFont('Segoe UI', 8))
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setStyleSheet('color: #888;')
        root.addWidget(subtitle)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet('color: #ddd;')
        root.addWidget(sep)

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
        self.load_btn = QPushButton('Load FASTA…')
        self.load_btn.clicked.connect(self._load_fasta)
        self.clear_btn = QPushButton('Clear')
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

        # Anti-Shine-Dalgarno
        self.asd_edit = QLineEdit('ACCTCCTTA')
        self.asd_edit.setMaxLength(9)
        self.asd_edit.setFont(QFont('Consolas', 9))
        self.asd_edit.setToolTip(
            '9 bp anti-Shine-Dalgarno sequence (3\' end of 16S rRNA). '
            'Default is E. coli ACCTCCTTA.'
        )
        opt_form.addRow('Anti-SD:', self.asd_edit)

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

        # ── Run button ────────────────────────────────────────────── #
        self.run_btn = QPushButton('▶  Run OSTIR')
        self.run_btn.setMinimumHeight(38)
        self.run_btn.setFont(QFont('Segoe UI', 10, QFont.Weight.Bold))
        self.run_btn.setStyleSheet(
            'QPushButton { background-color: #2e7d32; color: white; border-radius: 5px; }'
            'QPushButton:hover { background-color: #388e3c; }'
            'QPushButton:disabled { background-color: #aaa; }'
        )
        self.run_btn.clicked.connect(self._on_run)
        root.addWidget(self.run_btn)

        root.addStretch()

    # ------------------------------------------------------------------ #
    def _load_fasta(self):
        path, _ = QFileDialog.getOpenFileName(
            self, 'Open FASTA file', '',
            'FASTA files (*.fasta *.fa *.fna *.txt);;All files (*)'
        )
        if not path:
            return
        try:
            with open(path, 'r', encoding='utf-8') as fh:
                content = fh.read()
            self.seq_edit.setPlainText(content)

            # Auto-fill name from filename
            basename = os.path.splitext(os.path.basename(path))[0]
            if not self.seq_name_edit.text():
                self.seq_name_edit.setText(basename)
        except Exception as exc:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.critical(self, 'Error loading file', str(exc))

    # ------------------------------------------------------------------ #
    def _on_run(self):
        raw_text = self.seq_edit.toPlainText().strip()
        if not raw_text:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, 'No input', 'Please enter or load a sequence.')
            return

        sequence = self._parse_input(raw_text)
        if not sequence:
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, 'Invalid sequence',
                                'Could not extract a valid nucleotide sequence.')
            return

        params = {
            'name': self.seq_name_edit.text() or 'sequence',
            'start': self.start_spin.value() or None,
            'end': self.end_spin.value() or None,
            'aSD': self.asd_edit.text().strip() or None,
            'circular': self.circular_cb.isChecked(),
            'threads': self.threads_spin.value(),
            'decimal_places': self.decimals_spin.value(),
            'verbosity': self.verbosity_combo.currentIndex(),
            'constraints': self.constraints_edit.text().strip() or None,
        }

        self.run_requested.emit(sequence, params)

    # ------------------------------------------------------------------ #
    @staticmethod
    def _parse_input(text: str) -> str:
        """Accept a raw sequence string or FASTA. Returns clean sequence."""
        lines = text.strip().splitlines()
        if lines and lines[0].startswith('>'):
            # FASTA: take first entry
            seq_lines = [l.strip() for l in lines[1:] if not l.startswith('>')]
            return ''.join(seq_lines).upper().replace(' ', '').replace('\t', '')
        else:
            # Plain sequence
            return ''.join(lines).upper().replace(' ', '').replace('\t', '')

    # ------------------------------------------------------------------ #
    def set_running(self, running: bool):
        """Disable/enable the run button and input while OSTIR is running."""
        self.run_btn.setEnabled(not running)
        self.run_btn.setText('⏳  Running…' if running else '▶  Run OSTIR')
