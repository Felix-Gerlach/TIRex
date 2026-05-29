"""
Batch mode: run OSTIR over every record in a multi-FASTA (or several files)
and collect a combined, exportable table — one row per (record, start codon).
Runs in a background thread so the UI stays responsive.
"""

import os
import csv

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QProgressBar,
    QFileDialog, QMessageBox, QLineEdit, QFrame,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ostir import run_ostir

_COLS = ['Record', 'Start pos', 'Codon', 'TIR', 'dG total', 'RBS dist']


def _read_fasta_records(path):
    """Yield (name, sequence) from a (multi-)FASTA file."""
    name, buf = None, []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(buf)
                name = line[1:].strip().split()[0] or 'seq'
                buf = []
            else:
                buf.append(line)
    if name is not None:
        yield name, ''.join(buf)


class _BatchWorker(QThread):
    progress = pyqtSignal(int, int, str)
    row_ready = pyqtSignal(str, dict)
    done = pyqtSignal(int)
    failed = pyqtSignal(str)

    def __init__(self, files, aSD, parent=None):
        super().__init__(parent)
        self.files = files
        self.aSD = aSD

    def run(self):
        import re
        try:
            records = []
            for path in self.files:
                for name, seq in _read_fasta_records(path):
                    records.append((name, re.sub(r'[^ACGTU]', '',
                                    seq.upper().replace('U', 'T'))))
            total = len(records)
            n_rows = 0
            for i, (name, seq) in enumerate(records):
                self.progress.emit(i + 1, total, f'{name} ({len(seq)} nt)')
                if len(seq) < 10:
                    continue
                try:
                    res = run_ostir(seq, aSD=self.aSD or None, threads=1,
                                    decimal_places=4, verbosity=0) or []
                except Exception:
                    res = []
                for r in res:
                    self.row_ready.emit(name, r)
                    n_rows += 1
            self.done.emit(n_rows)
        except Exception as exc:
            import traceback
            self.failed.emit(f'{exc}\n{traceback.format_exc()}')


class BatchDialog(QDialog):
    def __init__(self, parent=None, default_asd='ACCTCCTTA'):
        super().__init__(parent)
        self.setWindowTitle('Batch OSTIR — multi-FASTA')
        self.resize(900, 640)
        self._files = []
        self._worker = None
        self._default_asd = default_asd
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        head = QFrame(); head.setObjectName('BHead')
        head.setStyleSheet('#BHead{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,'
                           'stop:0 #0f172a,stop:1 #0d9488);border-radius:10px;}')
        hl = QVBoxLayout(head); hl.setContentsMargins(14, 9, 14, 9)
        t = QLabel('Batch OSTIR'); t.setFont(QFont('Segoe UI', 14, QFont.Weight.Bold))
        t.setStyleSheet('color:white;'); hl.addWidget(t)
        sub = QLabel('Score every record in one or more FASTA files')
        sub.setStyleSheet('color:#99f6e4;font-size:8pt;'); hl.addWidget(sub)
        root.addWidget(head)

        row = QHBoxLayout()
        self.add_btn = QPushButton('Add FASTA…')
        self.add_btn.clicked.connect(self._add_files)
        row.addWidget(self.add_btn)
        self.files_lbl = QLabel('No files selected.')
        self.files_lbl.setProperty('muted', True)
        row.addWidget(self.files_lbl, 1)
        row.addWidget(QLabel('Anti-SD:'))
        self.asd_edit = QLineEdit(self._default_asd)
        self.asd_edit.setFixedWidth(110)
        self.asd_edit.setFont(QFont('Consolas', 9))
        row.addWidget(self.asd_edit)
        self.run_btn = QPushButton('Run batch  ▶')
        self.run_btn.setObjectName('AccentButton')
        self.run_btn.clicked.connect(self._run)
        row.addWidget(self.run_btn)
        root.addLayout(row)

        self.table = QTableWidget()
        self.table.setColumnCount(len(_COLS))
        self.table.setHorizontalHeaderLabels(_COLS)
        self.table.setSortingEnabled(True)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents)
        root.addWidget(self.table, 1)

        bottom = QHBoxLayout()
        self.status = QLabel('Add FASTA files and run.')
        self.status.setProperty('muted', True)
        bottom.addWidget(self.status, 1)
        self.progress = QProgressBar(); self.progress.setFixedWidth(220)
        self.progress.setVisible(False)
        bottom.addWidget(self.progress)
        self.export_btn = QPushButton('Export CSV…')
        self.export_btn.setObjectName('GhostButton')
        self.export_btn.clicked.connect(self._export)
        bottom.addWidget(self.export_btn)
        close_btn = QPushButton('Close')
        close_btn.setObjectName('GhostButton')
        close_btn.clicked.connect(self.accept)
        bottom.addWidget(close_btn)
        root.addLayout(bottom)

    def _add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, 'Add FASTA files', '',
            'FASTA (*.fasta *.fa *.fna *.txt);;All files (*)')
        if paths:
            self._files.extend(p for p in paths if p not in self._files)
            self.files_lbl.setText(f'{len(self._files)} file(s) selected.')

    def _run(self):
        if not self._files:
            QMessageBox.information(self, 'No files', 'Add at least one FASTA file.')
            return
        self.table.setRowCount(0)
        self.table.setSortingEnabled(False)
        self.run_btn.setEnabled(False)
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self._worker = _BatchWorker(self._files, self.asd_edit.text().strip(),
                                    parent=self)
        self._worker.progress.connect(self._on_progress)
        self._worker.row_ready.connect(self._on_row)
        self._worker.done.connect(self._on_done)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _on_progress(self, i, total, msg):
        self.progress.setRange(0, total)
        self.progress.setValue(i)
        self.status.setText(f'Scoring {msg}  ({i}/{total})')

    def _on_row(self, name, r):
        row = self.table.rowCount()
        self.table.insertRow(row)
        vals = [name, str(r.get('start_position', '')),
                r.get('start_codon', ''),
                f"{float(r.get('expression') or 0):.2f}",
                (f"{float(r['dG_total']):.2f}" if r.get('dG_total') is not None else '—'),
                str(r.get('RBS_distance_bp', ''))]
        for c, v in enumerate(vals):
            it = QTableWidgetItem(v)
            it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.table.setItem(row, c, it)

    def _on_done(self, n_rows):
        self.progress.setVisible(False)
        self.run_btn.setEnabled(True)
        self.table.setSortingEnabled(True)
        self.status.setText(f'Done — {n_rows} start codon(s) across '
                            f'{len(self._files)} file(s).')

    def _on_failed(self, msg):
        self.progress.setVisible(False)
        self.run_btn.setEnabled(True)
        QMessageBox.critical(self, 'Batch error', msg)

    def _export(self):
        if self.table.rowCount() == 0:
            QMessageBox.information(self, 'Nothing to export', 'Run a batch first.')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save CSV', 'tirex_batch.csv', 'CSV files (*.csv)')
        if not path:
            return
        try:
            with open(path, 'w', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh)
                w.writerow(_COLS)
                for row in range(self.table.rowCount()):
                    w.writerow([self.table.item(row, c).text()
                                if self.table.item(row, c) else ''
                                for c in range(len(_COLS))])
            self.status.setText(f'Exported {self.table.rowCount()} rows to {path}')
        except Exception as exc:
            QMessageBox.critical(self, 'Export error', str(exc))
