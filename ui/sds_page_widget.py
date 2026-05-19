"""
SDS-PAGE gel simulation dialog.

A standalone QDialog accepting:
  * the current ORF results
  * the selected target ORF (optional)
  * user-configurable lanes, marker, acrylamide, intensity mode, contaminants
and rendering an SDS-PAGE gel image (black bands on white, Coomassie style).
"""

import math
from dataclasses import asdict

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt6.QtWidgets import (
    QDialog, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel,
    QComboBox, QPushButton, QToolButton, QCheckBox, QRadioButton,
    QButtonGroup, QListWidget, QListWidgetItem, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QGroupBox,
    QSlider, QSplitter, QDoubleSpinBox, QLineEdit, QFileDialog,
    QMessageBox, QSizePolicy, QFrame,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont

from core.gel_simulator import (
    MARKERS, ACRYLAMIDE_OPTIONS, compute_rf,
    GelLane, GelBand, LANE_KINDS,
    default_contaminants,
)


# ======================================================================
#  Dialog
# ======================================================================
class SDSPageDialog(QDialog):

    def __init__(self, parent=None, results=None, target_orf_index=None,
                 tir_min: float = 0.0, tir_max: float = float('inf')):
        super().__init__(parent)
        self.setWindowTitle('SDS-PAGE Gel Simulator')
        self.resize(1400, 900)

        self._results = list(results or [])
        self._target_orf_index = target_orf_index
        self._tir_min = float(tir_min)
        self._tir_max = float(tir_max) if tir_max and tir_max > 0 else float('inf')
        self._contaminants = default_contaminants()

        # Default 4-lane layout: marker | contaminants | target group | mix
        self._lanes = [
            GelLane(name='Marker',           kind='Marker'),
            GelLane(name='Contaminants',     kind='Contaminants only'),
            GelLane(name='Target group',     kind='Target group (incl. alt starts)'),
            GelLane(name='Mixture',          kind='Target group + contaminants'),
        ]

        self._build_ui()
        self._refresh_lane_list()
        self._refresh_contaminant_table()
        self._rebuild_and_render()
        # Open maximized so the gel uses the full window
        self.showMaximized()

    # ------------------------------------------------------------------ #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # -------- Top toolbar: marker / acrylamide / intensity mode ------
        top = QHBoxLayout()
        top.setSpacing(10)

        top.addWidget(QLabel('Marker:'))
        self.marker_cmb = QComboBox()
        self.marker_cmb.addItems(MARKERS.keys())
        self.marker_cmb.setCurrentText('PageRuler Prestained Plus')
        self.marker_cmb.currentTextChanged.connect(self._rebuild_and_render)
        top.addWidget(self.marker_cmb)

        top.addSpacing(10)
        top.addWidget(QLabel('Acrylamide:'))
        self.acryl_cmb = QComboBox()
        self.acryl_cmb.addItems(ACRYLAMIDE_OPTIONS)
        self.acryl_cmb.setCurrentText('12%')
        self.acryl_cmb.currentTextChanged.connect(self._rebuild_and_render)
        top.addWidget(self.acryl_cmb)

        top.addSpacing(14)
        top.addWidget(QLabel('Band intensity:'))
        self.int_uniform_rb = QRadioButton('Uniform')
        self.int_tir_rb = QRadioButton('Proportional to TIR')
        self.int_uniform_rb.setChecked(True)
        grp = QButtonGroup(self)
        grp.addButton(self.int_uniform_rb)
        grp.addButton(self.int_tir_rb)
        self.int_uniform_rb.toggled.connect(self._rebuild_and_render)
        self.int_tir_rb.toggled.connect(self._rebuild_and_render)
        top.addWidget(self.int_uniform_rb)
        top.addWidget(self.int_tir_rb)

        top.addSpacing(14)
        top.addWidget(QLabel('Contaminant master:'))
        self.master_slider = QSlider(Qt.Orientation.Horizontal)
        self.master_slider.setRange(0, 200)
        self.master_slider.setValue(100)
        self.master_slider.setFixedWidth(160)
        self.master_lbl = QLabel('100%')
        self.master_lbl.setFixedWidth(42)
        self.master_slider.valueChanged.connect(
            lambda v: (self.master_lbl.setText(f'{v}%'), self._rebuild_and_render())
        )
        top.addWidget(self.master_slider)
        top.addWidget(self.master_lbl)

        top.addStretch()
        root.addLayout(top)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet('color: #ddd;')
        root.addWidget(sep)

        # -------- Main split: controls (left) | gel canvas (right) -------
        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)

        # ----- Left side: lanes + contaminants --------------------------
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.setSpacing(8)

        # Lanes group
        lanes_box = QGroupBox('Lanes')
        lanes_lay = QVBoxLayout(lanes_box)
        self.lane_list = QListWidget()
        self.lane_list.currentRowChanged.connect(self._on_lane_selected)
        lanes_lay.addWidget(self.lane_list)

        lane_btn_row = QHBoxLayout()
        self.add_lane_btn = QPushButton('+ Add')
        self.add_lane_btn.clicked.connect(self._add_lane)
        self.rem_lane_btn = QPushButton('− Remove')
        self.rem_lane_btn.clicked.connect(self._remove_lane)
        self.up_lane_btn = QPushButton('↑')
        self.up_lane_btn.clicked.connect(lambda: self._move_lane(-1))
        self.dn_lane_btn = QPushButton('↓')
        self.dn_lane_btn.clicked.connect(lambda: self._move_lane(+1))
        for b in (self.add_lane_btn, self.rem_lane_btn,
                   self.up_lane_btn, self.dn_lane_btn):
            b.setFixedHeight(24)
            b.setFont(QFont('Segoe UI', 8))
            lane_btn_row.addWidget(b)
        lanes_lay.addLayout(lane_btn_row)

        lane_cfg_grid = QGridLayout()
        lane_cfg_grid.addWidget(QLabel('Name:'), 0, 0)
        self.lane_name_edit = QLineEdit()
        self.lane_name_edit.editingFinished.connect(self._apply_lane_edits)
        lane_cfg_grid.addWidget(self.lane_name_edit, 0, 1)

        lane_cfg_grid.addWidget(QLabel('Contents:'), 1, 0)
        self.lane_kind_cmb = QComboBox()
        self.lane_kind_cmb.addItems(LANE_KINDS)
        self.lane_kind_cmb.currentTextChanged.connect(self._apply_lane_edits)
        lane_cfg_grid.addWidget(self.lane_kind_cmb, 1, 1)

        lanes_lay.addLayout(lane_cfg_grid)
        left_lay.addWidget(lanes_box)

        # Contaminants group
        cont_box = QGroupBox('E. coli His-tag contaminants')
        cont_lay = QVBoxLayout(cont_box)
        self.cont_table = QTableWidget()
        self.cont_table.setColumnCount(4)
        self.cont_table.setHorizontalHeaderLabels(
            ['On', 'Protein', 'MW (kDa)', 'Intensity (%)']
        )
        self.cont_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self.cont_table.verticalHeader().setVisible(False)
        self.cont_table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.cont_table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.SelectedClicked
        )
        self.cont_table.itemChanged.connect(self._on_cont_item_changed)
        cont_lay.addWidget(self.cont_table)
        left_lay.addWidget(cont_box)

        left.setMaximumWidth(360)
        left.setMinimumWidth(280)
        split.addWidget(left)

        # ----- Right: gel canvas ----------------------------------------
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure(facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        right_lay.addWidget(self.canvas)
        split.addWidget(right)

        split.setSizes([320, 1500])
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        root.addWidget(split, stretch=1)

        # -------- Bottom: export + close -------------------------------
        bot = QHBoxLayout()
        bot.addStretch()
        self.export_btn = QPushButton('Export as PNG…')
        self.export_btn.clicked.connect(self._export_png)
        bot.addWidget(self.export_btn)
        close_btn = QPushButton('Close')
        close_btn.clicked.connect(self.close)
        bot.addWidget(close_btn)
        root.addLayout(bot)

    # ================================================================== #
    #  Lane management
    # ================================================================== #
    def _refresh_lane_list(self):
        prev_row = self.lane_list.currentRow()
        self.lane_list.blockSignals(True)
        self.lane_list.clear()
        for lane in self._lanes:
            item = QListWidgetItem(f'{lane.name}  ({lane.kind})')
            self.lane_list.addItem(item)
        self.lane_list.blockSignals(False)
        if self._lanes:
            self.lane_list.setCurrentRow(
                min(max(prev_row, 0), len(self._lanes) - 1)
            )

    def _on_lane_selected(self, row: int):
        if not (0 <= row < len(self._lanes)):
            self.lane_name_edit.blockSignals(True)
            self.lane_kind_cmb.blockSignals(True)
            self.lane_name_edit.clear()
            self.lane_kind_cmb.setCurrentIndex(0)
            self.lane_name_edit.blockSignals(False)
            self.lane_kind_cmb.blockSignals(False)
            return
        lane = self._lanes[row]
        self.lane_name_edit.blockSignals(True)
        self.lane_kind_cmb.blockSignals(True)
        self.lane_name_edit.setText(lane.name)
        self.lane_kind_cmb.setCurrentText(lane.kind)
        self.lane_name_edit.blockSignals(False)
        self.lane_kind_cmb.blockSignals(False)

    def _apply_lane_edits(self):
        row = self.lane_list.currentRow()
        if not (0 <= row < len(self._lanes)):
            return
        lane = self._lanes[row]
        lane.name = self.lane_name_edit.text().strip() or f'Lane {row + 1}'
        lane.kind = self.lane_kind_cmb.currentText()
        self._refresh_lane_list()
        self._rebuild_and_render()

    def _add_lane(self):
        self._lanes.append(GelLane(
            name=f'Lane {len(self._lanes) + 1}',
            kind='Target + contaminants',
        ))
        self._refresh_lane_list()
        self.lane_list.setCurrentRow(len(self._lanes) - 1)
        self._rebuild_and_render()

    def _remove_lane(self):
        row = self.lane_list.currentRow()
        if not (0 <= row < len(self._lanes)):
            return
        if len(self._lanes) <= 1:
            return
        del self._lanes[row]
        self._refresh_lane_list()
        self._rebuild_and_render()

    def _move_lane(self, delta: int):
        row = self.lane_list.currentRow()
        new = row + delta
        if not (0 <= row < len(self._lanes)) or not (0 <= new < len(self._lanes)):
            return
        self._lanes[row], self._lanes[new] = self._lanes[new], self._lanes[row]
        self._refresh_lane_list()
        self.lane_list.setCurrentRow(new)
        self._rebuild_and_render()

    # ================================================================== #
    #  Contaminant table
    # ================================================================== #
    def _refresh_contaminant_table(self):
        self.cont_table.blockSignals(True)
        self.cont_table.setRowCount(len(self._contaminants))
        for i, c in enumerate(self._contaminants):
            # On/off checkbox as native item check state
            on = QTableWidgetItem()
            on.setFlags(
                Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsUserCheckable
            )
            on.setCheckState(
                Qt.CheckState.Checked if c.enabled else Qt.CheckState.Unchecked
            )
            self.cont_table.setItem(i, 0, on)

            name = QTableWidgetItem(c.name)
            name.setFlags(Qt.ItemFlag.ItemIsEnabled)
            name.setToolTip(c.note)
            self.cont_table.setItem(i, 1, name)

            mw = QTableWidgetItem(f'{c.mw_kda:.1f}')
            mw.setFlags(Qt.ItemFlag.ItemIsEnabled)
            mw.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.cont_table.setItem(i, 2, mw)

            it = QTableWidgetItem(f'{int(c.default_intensity * 100)}')
            it.setFlags(
                Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsEditable
            )
            it.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            self.cont_table.setItem(i, 3, it)
        self.cont_table.blockSignals(False)

    def _on_cont_item_changed(self, item: QTableWidgetItem):
        row = item.row()
        if not (0 <= row < len(self._contaminants)):
            return
        c = self._contaminants[row]
        col = item.column()
        if col == 0:
            c.enabled = (item.checkState() == Qt.CheckState.Checked)
        elif col == 3:
            try:
                pct = float(item.text())
                c.default_intensity = max(0.0, min(2.0, pct / 100.0))
            except ValueError:
                self.cont_table.blockSignals(True)
                item.setText(f'{int(c.default_intensity * 100)}')
                self.cont_table.blockSignals(False)
                return
        self._rebuild_and_render()

    # ================================================================== #
    #  Protein / ORF selection
    # ================================================================== #
    def _is_active(self, r) -> bool:
        """A result is 'active' if its checkbox is on AND its TIR sits in
        the user-selected display range."""
        if not r.get('visible', True):
            return False
        tir = r.get('expression') or 0
        return self._tir_min <= tir <= self._tir_max

    def _target_result(self):
        if self._target_orf_index is None:
            return None
        for r in self._results:
            if r.get('orf_index') == self._target_orf_index:
                return r if self._is_active(r) else None
        return None

    def _target_group(self):
        """ORFs sharing the target's stop codon — active ones only."""
        if self._target_orf_index is None:
            return []
        # Find target end even if target itself was filtered out, so the
        # group can still be computed.
        t_raw = next((r for r in self._results
                      if r.get('orf_index') == self._target_orf_index), None)
        if t_raw is None:
            return []
        tend = t_raw.get('end_position')
        return [r for r in self._results
                if r.get('end_position') == tend and self._is_active(r)]

    def _all_visible_orfs(self):
        return [r for r in self._results if self._is_active(r)]

    def _build_target_bands(self, kind: str):
        """Return list of GelBand for the given lane kind (non-marker)."""
        # Decide which ORFs to include
        if kind == 'Marker':
            return []
        if kind == 'Contaminants only':
            orfs = []
        elif kind == 'Target only':
            t = self._target_result()
            orfs = [t] if t else []
        elif kind == 'Target group (incl. alt starts)':
            orfs = self._target_group()
        elif kind == 'Target + contaminants':
            t = self._target_result()
            orfs = [t] if t else []
        elif kind == 'Target group + contaminants':
            orfs = self._target_group()
        elif kind in ('All visible ORFs', 'All visible ORFs + contaminants'):
            orfs = self._all_visible_orfs()
        else:
            orfs = []

        orfs = [r for r in orfs if r]  # drop None

        bands = []
        use_tir = self.int_tir_rb.isChecked()
        # Normalisation for TIR-proportional mode
        tir_max = max((r.get('expression') or 0) for r in orfs) if orfs else 1
        tir_max = max(tir_max, 1e-6)

        for r in orfs:
            mw = r.get('protein_molecular_weight')
            if not mw:
                continue
            if use_tir:
                tir = r.get('expression') or 0
                # log scale so low TIRs are still visible
                intensity = (math.log10(max(tir, 1)) /
                             max(math.log10(max(tir_max, 10)), 1.0))
                intensity = max(0.15, min(1.0, intensity))
            else:
                intensity = 0.9
            label = f'#{r.get("orf_index", 0) + 1} · {mw:.1f} kDa'
            bands.append(GelBand(
                mw_kda=float(mw),
                intensity=float(intensity),
                label=label,
                color='#111',
            ))

        # Add contaminants if appropriate
        if kind in ('Contaminants only', 'Target + contaminants',
                    'Target group + contaminants',
                    'All visible ORFs + contaminants'):
            master = self.master_slider.value() / 100.0
            for c in self._contaminants:
                if not c.enabled:
                    continue
                bands.append(GelBand(
                    mw_kda=c.mw_kda,
                    intensity=max(0.05, min(1.0, c.default_intensity * master)),
                    label=c.name,
                    color='#37474f',
                ))
        return bands

    # ================================================================== #
    #  Gel rendering
    # ================================================================== #
    def _rebuild_and_render(self):
        marker_name = self.marker_cmb.currentText()
        acryl = self.acryl_cmb.currentText()

        # Populate each lane's bands
        for lane in self._lanes:
            if lane.kind == 'Marker':
                lane.bands = [
                    GelBand(mw_kda=mw, intensity=0.95,
                            label=f'{int(mw) if mw >= 10 else mw}', color=col)
                    for mw, col in MARKERS.get(marker_name, [])
                ]
            else:
                lane.bands = self._build_target_bands(lane.kind)

        self._render_gel(acryl)

    def _render_gel(self, acrylamide: str):
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        n = max(1, len(self._lanes))
        lane_w = 0.82
        ax.set_xlim(0.2, n + 0.8)
        ax.set_ylim(-0.04, 1.08)
        ax.set_aspect('auto')
        ax.set_facecolor('#F7F7F7')

        # Resolving gel plate (light grey)
        plate = mpatches.Rectangle(
            (0.25, 0.02), n + 0.5, 0.96,
            facecolor='#fafafa', edgecolor='#999',
            linewidth=0.8, zorder=1,
        )
        ax.add_patch(plate)

        # Stacking gel tint at top
        stack = mpatches.Rectangle(
            (0.25, 0.92), n + 0.5, 0.06,
            facecolor='#eceff1', edgecolor='none', alpha=0.8, zorder=2,
        )
        ax.add_patch(stack)

        # Wells
        for i in range(n):
            well_x = (i + 1) - lane_w / 2 + 0.05
            well = mpatches.Rectangle(
                (well_x, 0.955), lane_w - 0.1, 0.025,
                facecolor='#cfd8dc', edgecolor='#607d8b',
                linewidth=0.5, zorder=3,
            )
            ax.add_patch(well)

        # Dye front line
        ax.plot([0.3, n + 0.7], [0.035, 0.035],
                 color='#42a5f5', lw=0.5, ls='--', alpha=0.6, zorder=2)

        # --- Bands -----------------------------------------------------
        band_h = 0.012
        for lane_idx, lane in enumerate(self._lanes):
            cx = lane_idx + 1
            x_left = cx - lane_w / 2
            x_right = cx + lane_w / 2

            # Lane vertical separator (subtle)
            if lane_idx > 0:
                ax.plot([x_left - 0.03] * 2, [0.04, 0.96],
                         color='#e0e0e0', lw=0.4, zorder=2)

            # Lane title
            ax.text(cx, 1.03, lane.name,
                    ha='center', va='bottom',
                    fontsize=9, fontweight='bold', color='#263238',
                    zorder=5)
            ax.text(cx, 1.005, lane.kind,
                    ha='center', va='bottom',
                    fontsize=7, color='#78909c', zorder=5)

            # Accumulate per-position intensity so overlapping bands merge
            # additively (simple merge: we still draw each rectangle).
            for band in lane.bands:
                rf = compute_rf(band.mw_kda, acrylamide)
                y = 1 - rf          # matplotlib up = top of gel
                # Convert to gel plate range [0.04, 0.94]
                y_in_gel = 0.04 + 0.90 * y / 1.0
                alpha = max(0.08, min(1.0, band.intensity))
                # Slight blur: draw 3 thin overlapping rects
                for dy, a_mul in ((-band_h * 0.4, 0.35),
                                   (0.0,           1.00),
                                   (+band_h * 0.4, 0.35)):
                    r = mpatches.Rectangle(
                        (x_left + 0.015, y_in_gel - band_h / 2 + dy),
                        (x_right - x_left) - 0.03, band_h,
                        facecolor=band.color, edgecolor='none',
                        alpha=alpha * a_mul, zorder=4,
                    )
                    ax.add_patch(r)

            # Band labels
            if lane.kind == 'Marker':
                # MW labels to the LEFT of the marker lane
                for band in lane.bands:
                    rf = compute_rf(band.mw_kda, acrylamide)
                    y_in_gel = 0.04 + 0.90 * (1 - rf)
                    ax.text(
                        x_left - 0.05, y_in_gel,
                        band.label,
                        ha='right', va='center',
                        fontsize=7, color=band.color,
                        fontweight='bold', zorder=6,
                    )
            else:
                # Band-specific labels to the right
                # Group labels that land on nearly-identical Rfs to avoid stack
                placed = []
                for band in sorted(lane.bands, key=lambda b: -b.mw_kda):
                    rf = compute_rf(band.mw_kda, acrylamide)
                    y_in_gel = 0.04 + 0.90 * (1 - rf)
                    # Nudge label if it collides with a previous one
                    y_lbl = y_in_gel
                    while any(abs(y_lbl - py) < 0.022 for py in placed):
                        y_lbl -= 0.022
                    placed.append(y_lbl)
                    ax.text(
                        x_right + 0.02, y_lbl,
                        band.label,
                        ha='left', va='center',
                        fontsize=6.5, color='#37474f', zorder=6,
                    )

        # ----- Footer annotations --------------------------------------
        ax.text(
            0.5 * (n + 1), -0.025,
            f'{acrylamide}  SDS-PAGE   •   '
            f'{"TIR-weighted" if self.int_tir_rb.isChecked() else "uniform"} band intensities',
            ha='center', va='top', fontsize=8, color='#546e7a',
        )

        ax.set_xticks([])
        ax.set_yticks([])
        for side in ('top', 'right', 'left', 'bottom'):
            ax.spines[side].set_visible(False)

        self.figure.subplots_adjust(left=0.04, right=0.98,
                                    top=0.94, bottom=0.06)
        self.canvas.draw_idle()

    # ================================================================== #
    #  Export
    # ================================================================== #
    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save gel image', 'sds_page_simulation.png',
            'PNG files (*.png);;SVG files (*.svg);;PDF files (*.pdf);;All (*)'
        )
        if not path:
            return
        try:
            self.figure.savefig(path, dpi=250, bbox_inches='tight',
                                 facecolor='white')
        except Exception as exc:
            QMessageBox.critical(self, 'Export error', str(exc))
