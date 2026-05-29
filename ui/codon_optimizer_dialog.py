"""
Standalone codon-optimization dialog.

Pick an ORF (or paste a CDS), choose a host, and synonymously rewrite the CDS
to maximise CAI without changing the protein. The first N codons can be
preserved (they overlap the initiation region and affect TIR).
"""

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox,
    QSpinBox, QPushButton, QPlainTextEdit, QGroupBox, QMessageBox, QFrame,
)
from PyQt6.QtGui import QFont, QGuiApplication
from PyQt6.QtCore import Qt

from core.codon_opt import optimize_cds, CODON_USAGE
from core.tuner.genetic_code import translate

MONO = QFont('Consolas', 10)


class CodonOptimizerDialog(QDialog):
    def __init__(self, parent=None, results=None):
        super().__init__(parent)
        self.setWindowTitle('Codon optimizer (CAI)')
        self.resize(820, 640)
        self._results = list(results or [])
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        head = QFrame()
        head.setObjectName('CoHead')
        head.setStyleSheet('#CoHead{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,'
                           'stop:0 #0f172a,stop:1 #6a1b9a);border-radius:10px;}')
        hl = QVBoxLayout(head); hl.setContentsMargins(14, 9, 14, 9)
        t = QLabel('Codon optimizer'); t.setFont(QFont('Segoe UI', 14, QFont.Weight.Bold))
        t.setStyleSheet('color:white;'); hl.addWidget(t)
        sub = QLabel('Synonymous CDS rewrite to raise CAI — protein unchanged')
        sub.setStyleSheet('color:#e9d5ff;font-size:8pt;'); hl.addWidget(sub)
        root.addWidget(head)

        # Controls
        box = QGroupBox('Input')
        g = QGridLayout(box)
        g.addWidget(QLabel('ORF / CDS source:'), 0, 0)
        self.orf_cmb = QComboBox()
        self.orf_cmb.addItem('— paste a CDS below —', None)
        for r in self._results:
            cds = r.get('dna_sequence') or ''
            if not cds:
                continue
            lbl = (f"ORF #{(r.get('orf_index', 0) + 1)} · pos {r.get('start_position')}"
                   f" · {r.get('protein_length', len(cds)//3)} aa")
            self.orf_cmb.addItem(lbl, cds)
        self.orf_cmb.currentIndexChanged.connect(self._on_orf_changed)
        g.addWidget(self.orf_cmb, 0, 1, 1, 3)

        g.addWidget(QLabel('Host:'), 1, 0)
        self.host_cmb = QComboBox()
        self.host_cmb.addItems(list(CODON_USAGE.keys()))
        g.addWidget(self.host_cmb, 1, 1)
        g.addWidget(QLabel('Preserve first N codons:'), 1, 2)
        self.preserve_spin = QSpinBox()
        self.preserve_spin.setRange(0, 50)
        self.preserve_spin.setValue(0)
        self.preserve_spin.setToolTip('These codons overlap the RBS/initiation '
                                      'region and affect TIR; leave them as-is.')
        g.addWidget(self.preserve_spin, 1, 3)
        root.addWidget(box)

        g.addWidget(QLabel('CDS (frame 0, start→stop):'), 2, 0)
        self.cds_edit = QPlainTextEdit()
        self.cds_edit.setFont(MONO)
        self.cds_edit.setFixedHeight(80)
        self.cds_edit.setPlaceholderText('Paste a coding sequence, or pick an ORF above.')
        g.addWidget(self.cds_edit, 3, 0, 1, 4)

        self.opt_btn = QPushButton('Optimize  ▶')
        self.opt_btn.setObjectName('AccentButton')
        self.opt_btn.clicked.connect(self._optimize)
        root.addWidget(self.opt_btn)

        self.metrics_lbl = QLabel('')
        self.metrics_lbl.setTextFormat(Qt.TextFormat.RichText)
        root.addWidget(self.metrics_lbl)

        out_box = QGroupBox('Optimized CDS')
        ol = QVBoxLayout(out_box)
        self.out_edit = QPlainTextEdit()
        self.out_edit.setFont(MONO)
        self.out_edit.setReadOnly(True)
        ol.addWidget(self.out_edit)
        root.addWidget(out_box, 1)

        bottom = QHBoxLayout()
        bottom.addStretch()
        self.copy_btn = QPushButton('Copy optimized CDS')
        self.copy_btn.setObjectName('GhostButton')
        self.copy_btn.clicked.connect(self._copy)
        bottom.addWidget(self.copy_btn)
        close_btn = QPushButton('Close')
        close_btn.setObjectName('GhostButton')
        close_btn.clicked.connect(self.accept)
        bottom.addWidget(close_btn)
        root.addLayout(bottom)

    def _on_orf_changed(self):
        cds = self.orf_cmb.currentData()
        if cds:
            self.cds_edit.setPlainText(cds)

    def _optimize(self):
        import re
        cds = re.sub(r'[^ACGTUacgtu]', '', self.cds_edit.toPlainText().upper()).replace('U', 'T')
        if len(cds) < 6:
            QMessageBox.warning(self, 'No CDS', 'Provide a coding sequence (≥ 2 codons).')
            return
        host = self.host_cmb.currentText()
        opt, st = optimize_cds(cds, host=host,
                               preserve_first=self.preserve_spin.value())
        ok = translate(cds) == translate(opt)
        self.out_edit.setPlainText(opt)
        fold = (st['cai_after'] / st['cai_before']) if st['cai_before'] else 1.0
        self.metrics_lbl.setText(
            f"<b>CAI</b> {st['cai_before']:.3f} → "
            f"<b style='color:#16a34a'>{st['cai_after']:.3f}</b> ({fold:.2f}×) &nbsp;·&nbsp; "
            f"<b>GC</b> {st['gc_before']:.1f}% → {st['gc_after']:.1f}% &nbsp;·&nbsp; "
            f"<b>{st['codons_changed']}</b>/{st['codons_total']} codons changed &nbsp;·&nbsp; "
            f"protein preserved: <b style='color:{'#16a34a' if ok else '#e11d48'}'>"
            f"{'yes' if ok else 'NO'}</b>")

    def _copy(self):
        seq = self.out_edit.toPlainText().strip()
        if seq:
            QGuiApplication.clipboard().setText(seq)
