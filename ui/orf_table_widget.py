"""
ORF results table with per-ORF visibility toggles, sorting, and a per-row
"Copy AA" button that copies the full amino-acid sequence to the clipboard.

No detail panel / composition chart — the table fills the bottom area.
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QHeaderView, QAbstractItemView, QPushButton, QMenu,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont, QColor, QBrush, QGuiApplication

# Columns: (header, result_key_or_None, format_fn_or_None, numeric?)
# The '#' column uses a synthetic display rank (populated per-table) rather
# than the raw orf_index, so it restarts from 1 when the TIR filter hides rows.
COLUMNS = [
    ('Show',          None,                        None,                                         False),
    ('Copy AA',       None,                        None,                                         False),
    ('#',             None,                        None,                                         True),
    ('Position',      'start_position',            str,                                          True),
    ('Codon',         'start_codon',               str,                                          False),
    ('TIR (1° aSD)',  'expression',                lambda v: f'{v:.4f}' if v else '—',           True),
    ('TIR (2° aSD)',  'expression_secondary',      lambda v: f'{v:.4f}' if v else '—',           True),
    ('RBS dist',      'RBS_distance_bp',           lambda v: str(v),                             True),
    ('dG total',      'dG_total',                  lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('dG rRNA:mRNA',  'dG_rRNA:mRNA',              lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('dG mRNA',       'dG_mRNA',                   lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('dG spacing',    'dG_spacing',                lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('dG standby',    'dG_standby',                lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('dG start',      'dG_start_codon',            lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('Length (aa)',   'protein_length',            lambda v: str(v) if v else '—',               True),
    ('MW (kDa)',      'protein_molecular_weight',  lambda v: f'{v:.3f}' if v else '—',           True),
    ('pI',            'protein_isoelectric_point', lambda v: f'{v:.2f}' if v else '—',           True),
    ('GRAVY',         'protein_gravy',             lambda v: f'{v:.3f}' if v is not None else '—', True),
    ('Stop?',         'has_stop',                  lambda v: '✓' if v else '✗',                  False),
]

COL_IDX = {c[0]: i for i, c in enumerate(COLUMNS)}
SHOW_COL = COL_IDX['Show']
COPY_COL = COL_IDX['Copy AA']

CODON_COLORS = {
    'ATG': QColor('#e53935'),
    'GTG': QColor('#fb8c00'),
    'TTG': QColor('#fdd835'),
}


class _NumericItem(QTableWidgetItem):
    """Table item that sorts by a stored numeric value (UserRole+1) instead of
    its display text, so 67.9 sorts before 4385.2 rather than after it.
    Items with no numeric value (e.g. '—') sort to the bottom ascending."""

    def __lt__(self, other):
        a = self.data(Qt.ItemDataRole.UserRole + 1)
        b = other.data(Qt.ItemDataRole.UserRole + 1)
        if a is None and b is None:
            return super().__lt__(other)
        if a is None:
            return False
        if b is None:
            return True
        try:
            return float(a) < float(b)
        except (TypeError, ValueError):
            return super().__lt__(other)


class ORFTableWidget(QWidget):

    visibility_changed = pyqtSignal(int, bool)
    selection_changed  = pyqtSignal(int)
    aa_copied          = pyqtSignal(int, int)   # (orf_index, aa_length)
    target_requested   = pyqtSignal(int)        # orf_index of target
    target_cleared     = pyqtSignal()
    open_gel_requested = pyqtSignal()
    tune_requested     = pyqtSignal(int)        # orf_index to open in TIR-Tuner

    # ------------------------------------------------------------------ #
    def __init__(self, parent=None):
        super().__init__(parent)
        self._results = []
        self._updating = False    # guard against itemChanged loops
        self._tir_min = 0.0
        self._tir_max = float('inf')
        self._frame_anchor = None   # target start position; hide out-of-frame
        self._build_ui()

    # ------------------------------------------------------------------ #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        title_row = QHBoxLayout()
        title_lbl = QLabel('ORF Results')
        title_lbl.setFont(QFont('Segoe UI', 10, QFont.Weight.Bold))
        title_row.addWidget(title_lbl)
        title_row.addStretch()

        self.show_all_btn = QPushButton('Show all')
        self.hide_all_btn = QPushButton('Hide all')
        self.clear_target_btn = QPushButton('Clear target')
        self.gel_btn = QPushButton('SDS-PAGE…')
        for btn in (self.show_all_btn, self.hide_all_btn,
                    self.clear_target_btn, self.gel_btn):
            btn.setFixedHeight(22)
            btn.setFont(QFont('Segoe UI', 8))
        self.show_all_btn.clicked.connect(lambda: self._set_all_visibility(True))
        self.hide_all_btn.clicked.connect(lambda: self._set_all_visibility(False))
        self.clear_target_btn.clicked.connect(self.target_cleared.emit)
        self.gel_btn.clicked.connect(self.open_gel_requested.emit)
        title_row.addWidget(self.show_all_btn)
        title_row.addWidget(self.hide_all_btn)
        title_row.addWidget(self.clear_target_btn)
        title_row.addWidget(self.gel_btn)
        root.addLayout(title_row)

        # Table (fills remaining space)
        self.table = QTableWidget()
        self.table.setColumnCount(len(COLUMNS))
        self.table.setHorizontalHeaderLabels([c[0] for c in COLUMNS])
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setDefaultSectionSize(20)
        self.table.verticalHeader().setVisible(False)
        self.table.setFont(QFont('Segoe UI', 8))
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        self.table.cellClicked.connect(self._on_cell_clicked)
        self.table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._on_context_menu)
        root.addWidget(self.table)

    # ------------------------------------------------------------------ #
    #  Public API                                                          #
    # ------------------------------------------------------------------ #
    def update_results(self, results: list):
        self._results = results
        self._frame_anchor = None     # fresh data → no target frame filter
        self._populate_table()

    def clear(self):
        self._results = []
        self.table.setRowCount(0)

    def set_tir_range(self, lo: float, hi: float):
        """Called by the visualization widget when the TIR range filter
        changes. Rows with TIR outside [lo, hi] are hidden from the table
        and the remaining rows are renumbered from 1."""
        self._tir_min = float(lo)
        # QVariant can't carry math.inf cleanly over a signal boundary, so
        # treat any value <= 0 as "no upper limit".
        self._tir_max = float(hi) if hi and hi > 0 else float('inf')
        self._populate_table()

    def set_frame_filter(self, target_start):
        """Hide every fragment NOT in the target's reading frame (frame =
        start position modulo 3). Pass None to clear."""
        self._frame_anchor = (int(target_start)
                              if target_start is not None else None)
        self._populate_table()

    def clear_frame_filter(self):
        self._frame_anchor = None
        self._populate_table()

    # ------------------------------------------------------------------ #
    #  Table population                                                    #
    # ------------------------------------------------------------------ #
    def _populate_table(self):
        self._updating = True
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)

        # ------------------------------------------------------------
        # Apply TIR + reading-frame filters; compute display ranks.
        # ------------------------------------------------------------
        def passes(r):
            t = r.get('expression') or 0
            if not (self._tir_min <= t <= self._tir_max):
                return False
            if self._frame_anchor is not None:
                sp = r.get('start_position')
                if sp is None or (sp - self._frame_anchor) % 3 != 0:
                    return False
            return True

        filtered = [r for r in self._results if passes(r)]
        # Rank by start position so # matches the visualization's numbering.
        ranked = sorted(filtered, key=lambda r: r.get('start_position', 0))
        display_rank = {
            r.get('orf_index'): i + 1 for i, r in enumerate(ranked)
        }

        # Preserve the original OSTIR order within the filtered subset for
        # the table rows — users can still click any column header to re-sort.
        rows_to_show = [r for r in self._results if passes(r)]

        self.table.setRowCount(len(rows_to_show))

        for row_idx, r in enumerate(rows_to_show):
            for col_idx, (header, key, fmt, numeric) in enumerate(COLUMNS):
                if header == 'Show':
                    item = QTableWidgetItem()
                    item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled
                        | Qt.ItemFlag.ItemIsSelectable
                        | Qt.ItemFlag.ItemIsUserCheckable
                    )
                    item.setCheckState(
                        Qt.CheckState.Checked if r.get('visible', True)
                        else Qt.CheckState.Unchecked
                    )
                    item.setData(Qt.ItemDataRole.UserRole,
                                 r.get('orf_index', row_idx))
                    self.table.setItem(row_idx, col_idx, item)
                    continue

                if header == '#':
                    rank = display_rank.get(r.get('orf_index'), row_idx + 1)
                    item = _NumericItem(str(rank))
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    item.setFlags(
                        Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                    )
                    # Store numeric sort key
                    item.setData(Qt.ItemDataRole.UserRole + 1, float(rank))
                    # Bold so the identifier stands out
                    item.setFont(QFont('Segoe UI', 8, QFont.Weight.Bold))
                    self.table.setItem(row_idx, col_idx, item)
                    continue

                if header == 'Copy AA':
                    aa = r.get('aa_sequence', '') or ''
                    item = QTableWidgetItem('📋 Copy' if aa else '—')
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                    item.setData(Qt.ItemDataRole.UserRole,
                                 r.get('orf_index', row_idx))
                    item.setData(Qt.ItemDataRole.UserRole + 2, aa)
                    if aa:
                        item.setForeground(QBrush(QColor('#1565c0')))
                        item.setToolTip(
                            f'Click to copy {len(aa)} aa to clipboard\n\n'
                            f'{aa[:80]}{"…" if len(aa) > 80 else ""}'
                        )
                    self.table.setItem(row_idx, col_idx, item)
                    continue

                val = r.get(key, None)
                if val is None or val == '':
                    text = '—'
                else:
                    try:
                        text = fmt(val) if fmt else str(val)
                    except Exception:
                        text = str(val)

                item = _NumericItem(text) if numeric else QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)

                # Numeric sort key
                if numeric and val is not None and val != '':
                    try:
                        item.setData(Qt.ItemDataRole.UserRole + 1, float(val))
                    except (TypeError, ValueError):
                        pass

                # Codon colour
                if header == 'Codon':
                    col = CODON_COLORS.get(val)
                    if col:
                        item.setForeground(QBrush(col))
                        item.setFont(QFont('Consolas', 8, QFont.Weight.Bold))

                # TIR background tint (both anti-SD columns)
                if header.startswith('TIR') and val:
                    try:
                        tv = float(val)
                        if tv > 10000:
                            item.setBackground(QBrush(QColor('#e8f5e9')))
                        elif tv > 1000:
                            item.setBackground(QBrush(QColor('#fff9c4')))
                    except Exception:
                        pass

                # Stop colour
                if header == 'Stop?':
                    item.setForeground(QBrush(
                        QColor('#2e7d32') if r.get('has_stop') else QColor('#c62828')
                    ))

                self.table.setItem(row_idx, col_idx, item)

        self.table.setSortingEnabled(True)
        self._updating = False

    # ------------------------------------------------------------------ #
    #  Visibility toggle                                                   #
    # ------------------------------------------------------------------ #
    def _on_item_changed(self, item: QTableWidgetItem):
        if self._updating or item.column() != SHOW_COL:
            return
        orf_index = item.data(Qt.ItemDataRole.UserRole)
        if orf_index is None:
            return
        visible = (item.checkState() == Qt.CheckState.Checked)
        for r in self._results:
            if r.get('orf_index') == orf_index:
                r['visible'] = visible
                break
        self.visibility_changed.emit(int(orf_index), visible)

    def _set_all_visibility(self, visible: bool):
        self._updating = True
        for row in range(self.table.rowCount()):
            item = self.table.item(row, SHOW_COL)
            if item:
                item.setCheckState(
                    Qt.CheckState.Checked if visible else Qt.CheckState.Unchecked
                )
                orf_idx = item.data(Qt.ItemDataRole.UserRole)
                if orf_idx is not None:
                    for r in self._results:
                        if r.get('orf_index') == orf_idx:
                            r['visible'] = visible
                            break
        self._updating = False
        for r in self._results:
            self.visibility_changed.emit(r.get('orf_index', 0), visible)

    # ------------------------------------------------------------------ #
    #  Copy AA                                                             #
    # ------------------------------------------------------------------ #
    def _on_cell_clicked(self, row: int, col: int):
        if col != COPY_COL:
            return
        item = self.table.item(row, col)
        if item is None:
            return
        aa = item.data(Qt.ItemDataRole.UserRole + 2) or ''
        if not aa:
            return
        QGuiApplication.clipboard().setText(aa)
        orf_idx = item.data(Qt.ItemDataRole.UserRole)
        self.aa_copied.emit(int(orf_idx) if orf_idx is not None else -1, len(aa))

    # ------------------------------------------------------------------ #
    #  Row selection                                                       #
    # ------------------------------------------------------------------ #
    def _on_selection_changed(self):
        items = self.table.selectedItems()
        if not items:
            return
        row = items[0].row()
        show_item = self.table.item(row, SHOW_COL)
        if show_item is None:
            return
        orf_index = show_item.data(Qt.ItemDataRole.UserRole)
        if orf_index is None:
            return
        self.selection_changed.emit(int(orf_index))

    # ------------------------------------------------------------------ #
    #  Context menu                                                        #
    # ------------------------------------------------------------------ #
    def _on_context_menu(self, pos):
        item = self.table.itemAt(pos)
        if item is None:
            return
        row = item.row()
        show_item = self.table.item(row, SHOW_COL)
        orf_index = (show_item.data(Qt.ItemDataRole.UserRole)
                     if show_item is not None else None)

        menu = QMenu(self.table)
        if orf_index is not None:
            act_target = menu.addAction('Set as target')
            act_target.triggered.connect(
                lambda: self.target_requested.emit(int(orf_index))
            )
            act_tune = menu.addAction('Tune translation rate…')
            act_tune.triggered.connect(
                lambda: self.tune_requested.emit(int(orf_index))
            )
        act_clear = menu.addAction('Clear target')
        act_clear.triggered.connect(self.target_cleared.emit)
        menu.addSeparator()
        act_gel = menu.addAction('Simulate SDS-PAGE…')
        act_gel.triggered.connect(self.open_gel_requested.emit)
        menu.exec(self.table.viewport().mapToGlobal(pos))

    def refresh(self):
        """Repopulate the table from the current results (e.g. after the
        main window has mutated `visible` flags)."""
        self._populate_table()
