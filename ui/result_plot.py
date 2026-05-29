"""
Graphical view of the substitution scan.

A saturation-mutagenesis-style heatmap:
  * x-axis  = sequence position (contiguous over the scanned window)
  * y-axis  = the four bases (A / C / G / T)
  * colour  = ΔTIR of substituting that base at that position
  * the wild-type base in each column is outlined and labelled
An annotation track above the heatmap shows user-defined, labelled regions
(SD, spacer, start codon, …) that can be added / removed at runtime.

Only single-nucleotide substitutions are plotted (the user asked for the
substitution results); deletions / insertions / multi-base edits stay in the
table.
"""

import re
from itertools import cycle

import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavToolbar
import matplotlib.patches as mpatches

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, QLineEdit,
    QPushButton, QListWidget, QListWidgetItem, QGroupBox, QColorDialog,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

BASES = ['A', 'C', 'G', 'T']
# A on top, T at bottom -> row index for pcolormesh (row 0 = bottom)
BASE_ROW = {'A': 3, 'C': 2, 'G': 1, 'T': 0}

REGION_COLORS = cycle([
    '#1976d2', '#e53935', '#2e7d32', '#f9a825', '#6a1b9a',
    '#00838f', '#d84315', '#5d4037',
])

# Accept 'G->A' or '74:G->A'
_SUB_RE = re.compile(r'^(?:\d+:)?([ACGT])->([ACGT])$')


class ResultPlotWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._sequence = ''
        self._mode = 'increase'
        self._candidates = []
        self._min_pos = 1
        self._max_pos = 1
        self._regions = []     # list of dicts: {from, to, label, color}
        self._build_ui()

    # ------------------------------------------------------------------ #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(2, 2, 2, 2)

        self.figure = Figure(facecolor='white')
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavToolbar(self.canvas, self)
        root.addWidget(self.toolbar)
        root.addWidget(self.canvas, stretch=1)

        # ----- Region-label controls -----------------------------------
        ann = QGroupBox('Region labels')
        albox = QVBoxLayout(ann)

        row = QHBoxLayout()
        row.addWidget(QLabel('From:'))
        self.from_spin = QSpinBox()
        self.from_spin.setRange(1, 1_000_000)
        row.addWidget(self.from_spin)
        row.addWidget(QLabel('To:'))
        self.to_spin = QSpinBox()
        self.to_spin.setRange(1, 1_000_000)
        row.addWidget(self.to_spin)
        row.addWidget(QLabel('Label:'))
        self.label_edit = QLineEdit()
        self.label_edit.setPlaceholderText('e.g. SD / spacer / start codon')
        row.addWidget(self.label_edit, 1)
        self._pending_color = next(REGION_COLORS)
        self.color_btn = QPushButton('Colour')
        self.color_btn.clicked.connect(self._pick_color)
        self._sync_color_btn()
        row.addWidget(self.color_btn)
        self.add_btn = QPushButton('Add')
        self.add_btn.clicked.connect(self._add_region)
        row.addWidget(self.add_btn)
        albox.addLayout(row)

        row2 = QHBoxLayout()
        self.region_list = QListWidget()
        self.region_list.setFixedHeight(70)
        row2.addWidget(self.region_list, 1)
        btns = QVBoxLayout()
        self.del_btn = QPushButton('Remove selected')
        self.del_btn.clicked.connect(self._remove_region)
        self.clear_btn = QPushButton('Clear all')
        self.clear_btn.clicked.connect(self._clear_regions)
        btns.addWidget(self.del_btn)
        btns.addWidget(self.clear_btn)
        btns.addStretch()
        row2.addLayout(btns)
        albox.addLayout(row2)

        root.addWidget(ann)

    # ------------------------------------------------------------------ #
    def _sync_color_btn(self):
        self.color_btn.setStyleSheet(
            f'background:{self._pending_color}; color:white;')

    def _pick_color(self):
        c = QColorDialog.getColor(QColor(self._pending_color), self,
                                  'Region colour')
        if c.isValid():
            self._pending_color = c.name()
            self._sync_color_btn()

    # ------------------------------------------------------------------ #
    #  Public API
    # ------------------------------------------------------------------ #
    def set_candidates(self, candidates, sequence, mode):
        self._candidates = candidates or []
        self._sequence = sequence or ''
        self._mode = mode
        subs = [c for c in self._candidates if self._parse_sub(c) is not None]
        if subs:
            positions = [c.position for c in subs]
            self._min_pos = min(positions)
            self._max_pos = max(positions)
            # Pre-fill the region spin defaults to the scanned window.
            self.from_spin.setValue(self._min_pos)
            self.to_spin.setValue(self._max_pos)
        self._render()

    # ------------------------------------------------------------------ #
    def _parse_sub(self, cand):
        """Return (position, new_base) if cand is a single substitution, else None."""
        if cand.edit_type not in ('substitution', 'codon swap'):
            return None
        m = _SUB_RE.match(cand.change.strip())
        if not m:
            return None
        return (cand.position, m.group(2))

    # ------------------------------------------------------------------ #
    #  Region management
    # ------------------------------------------------------------------ #
    def _add_region(self):
        a, b = self.from_spin.value(), self.to_spin.value()
        if b < a:
            a, b = b, a
        label = self.label_edit.text().strip() or f'{a}–{b}'
        self._regions.append({'from': a, 'to': b, 'label': label,
                              'color': self._pending_color})
        item = QListWidgetItem(f'{label}  [{a}–{b}]')
        item.setForeground(QColor(self._pending_color))
        self.region_list.addItem(item)
        self.label_edit.clear()
        self._pending_color = next(REGION_COLORS)
        self._sync_color_btn()
        self._render()

    def _remove_region(self):
        row = self.region_list.currentRow()
        if 0 <= row < len(self._regions):
            del self._regions[row]
            self.region_list.takeItem(row)
            self._render()

    def _clear_regions(self):
        self._regions.clear()
        self.region_list.clear()
        self._render()

    # ------------------------------------------------------------------ #
    #  Rendering
    # ------------------------------------------------------------------ #
    def _render(self):
        self.figure.clear()
        subs = [c for c in self._candidates if self._parse_sub(c) is not None]
        if not subs:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5,
                    'No single-substitution results to plot.\n'
                    '(Run a prediction; deletions/insertions stay in the table.)',
                    ha='center', va='center', color='#777')
            ax.axis('off')
            self.canvas.draw_idle()
            return

        lo, hi = self._min_pos, self._max_pos
        positions = list(range(lo, hi + 1))
        ncol = len(positions)
        col_of = {p: i for i, p in enumerate(positions)}

        # Build ΔTIR matrix (4 rows x ncol), NaN where untested.
        mat = np.full((4, ncol), np.nan)
        baseline = subs[0].baseline_tir
        for c in subs:
            pos, base = self._parse_sub(c)
            if pos not in col_of:
                continue
            mat[BASE_ROW[base], col_of[pos]] = c.delta

        masked = np.ma.masked_invalid(mat)
        vmax = np.nanmax(np.abs(mat)) if np.isfinite(mat).any() else 1.0
        vmax = vmax or 1.0
        cmap = matplotlib.colormaps.get_cmap('RdYlGn').copy()
        cmap.set_bad('#eeeeee')
        norm = matplotlib.colors.Normalize(vmin=-vmax, vmax=vmax)

        # Layout: heatmap + region track on top.
        ax = self.figure.add_subplot(111)
        mesh = ax.pcolormesh(np.arange(ncol + 1), np.arange(5), masked,
                             cmap=cmap, norm=norm, edgecolors='#ffffff',
                             linewidth=0.4)

        # Wild-type cell markers + base letters.
        for i, p in enumerate(positions):
            seq_idx = p - 1
            if 0 <= seq_idx < len(self._sequence):
                wt = self._sequence[seq_idx].upper()
                if wt in BASE_ROW:
                    r = BASE_ROW[wt]
                    ax.add_patch(mpatches.Rectangle(
                        (i, r), 1, 1, fill=False, edgecolor='#111',
                        linewidth=1.6, zorder=5))
                    ax.text(i + 0.5, r + 0.5, wt, ha='center', va='center',
                            fontsize=8, fontweight='bold', color='#111',
                            zorder=6)

        # Annotate ΔTIR values when the grid is small enough.
        if ncol <= 28:
            for c in subs:
                pos, base = self._parse_sub(c)
                if pos not in col_of:
                    continue
                i = col_of[pos]
                r = BASE_ROW[base]
                ax.text(i + 0.5, r + 0.5, f'{c.delta:+.0f}',
                        ha='center', va='center', fontsize=6,
                        color='#222', zorder=4)

        ax.set_yticks([r + 0.5 for r in range(4)])
        ax.set_yticklabels(['T', 'G', 'C', 'A'])   # row0=T..row3=A
        # x ticks: position numbers (thin out if dense)
        step = max(1, ncol // 30)
        ax.set_xticks([i + 0.5 for i in range(0, ncol, step)])
        ax.set_xticklabels([str(positions[i]) for i in range(0, ncol, step)],
                           rotation=90, fontsize=7)
        ax.set_xlabel('Sequence position (nt)')
        ax.set_ylabel('Substituted base')
        ax.set_xlim(0, ncol)

        # ---- Region annotation track (above the heatmap) --------------
        track_y0, track_h = 4.25, 0.55
        top = 4.25
        for reg in self._regions:
            a = max(reg['from'], lo)
            b = min(reg['to'], hi)
            if b < a:
                continue
            x0 = col_of[a]
            x1 = col_of[b] + 1
            ax.add_patch(mpatches.Rectangle(
                (x0, track_y0), x1 - x0, track_h,
                facecolor=reg['color'], edgecolor='none', alpha=0.85,
                clip_on=False, zorder=7))
            ax.text((x0 + x1) / 2, track_y0 + track_h / 2, reg['label'],
                    ha='center', va='center', fontsize=8, color='white',
                    fontweight='bold', clip_on=False, zorder=8)
            top = track_y0 + track_h
        ax.set_ylim(0, max(4.0, top + 0.15) if self._regions else 4.0)

        title = ('ΔTIR per substitution  ·  '
                 f'{"increase" if self._mode == "increase" else "decrease"} mode'
                 f'  ·  wild-type TIR = {baseline:.1f}')
        ax.set_title(title, fontsize=10)

        cbar = self.figure.colorbar(mesh, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label('ΔTIR vs wild-type')

        self.figure.tight_layout()
        self.canvas.draw_idle()
