"""
Clean, minimal sequence visualization.

Design principles
-----------------
One main track of arrow-shaped ORFs, stacked vertically only when they
actually overlap.  RBS position (from OSTIR's RBS_distance_bp) and start-codon
type are rendered as small indicators on the same ORF element, so the whole
prediction reads as one integrated glyph rather than four separate tracks.
A single subtle backbone line and a single ruler tie everything together.
"""

import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.path import Path
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QButtonGroup, QFrame, QSizePolicy, QDoubleSpinBox, QCheckBox,
    QToolButton,
)
from PyQt6.QtCore import QTimer, Qt, pyqtSignal
from PyQt6.QtGui import QFont

matplotlib.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.facecolor': '#ffffff',
    'figure.facecolor': '#ffffff',
    'path.simplify': True,
    'path.simplify_threshold': 1.0,
    'agg.path.chunksize': 10000,
})

CODON_COLORS = {
    'ATG': '#c62828',   # red
    'GTG': '#ef6c00',   # orange
    'TTG': '#f9a825',   # amber
}
RBS_COLOR = '#FFC107'
RBS_EDGE = '#F57F17'
BACKBONE_COLOR = '#90a4ae'

# Per-nucleotide box colors (common bioinformatics scheme).
NT_COLORS = {
    'A': '#4CAF50',   # green
    'T': '#E53935',   # red
    'U': '#E53935',   # red (RNA)
    'G': '#FDD835',   # yellow
    'C': '#1E88E5',   # blue
    'N': '#BDBDBD',   # grey fallback
}
# Letter color: white on dark boxes, near-black on yellow.
NT_FG = {
    'A': 'white',
    'T': 'white',
    'U': 'white',
    'G': '#222',
    'C': 'white',
    'N': '#222',
}
AA_BOX_COLOR = '#FAFAFA'   # near-white, semi-transparent
AA_BOX_EDGE  = '#263238'   # dark outline
AA_TEXT_COLOR = '#111'

STOP_CAP_COLOR = '#37474f'   # dark slate — neutral "stop" marker

# Palette for "group by shared ORF" mode — 12 distinguishable hues.
GROUP_PALETTE = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#7f7f7f',
    '#aec7e8', '#ffbb78',
]
GROUP_SINGLETON_COLOR = '#bdbdbd'   # ORFs not sharing a stop


# =====================================================================
#  Geometry helpers
# =====================================================================
def _compute_norm(results, tir_min=0.0, tir_max=float('inf')):
    vals = []
    for r in results:
        if not r.get('visible', True):
            continue
        t = r.get('expression') or 0
        if not (tir_min <= t <= tir_max):
            continue
        vals.append(math.log10(max(t or 1e-9, 1e-9)))
    if not vals:
        return Normalize(vmin=0, vmax=1)
    lo, hi = min(vals), max(vals)
    if lo == hi:
        lo, hi = lo - 0.5, hi + 0.5
    return Normalize(vmin=lo, vmax=hi)


def _tir_color(tir_value, norm, cmap):
    return cmap(norm(math.log10(max(tir_value or 1e-9, 1e-9))))


def _arrow_path(x_start, x_end, y_bot, y_top, head_bp):
    """Return a Path for a right-pointing arrow polygon."""
    head_bp = max(1.0, min(head_bp, (x_end - x_start) * 0.45))
    body_end = x_end - head_bp
    if body_end <= x_start:
        # Tiny ORF → just a triangle
        mid_y = (y_bot + y_top) / 2
        verts = [(x_start, y_bot), (x_end, mid_y),
                 (x_start, y_top), (x_start, y_bot)]
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        return Path(verts, codes)

    mid_y = (y_bot + y_top) / 2
    verts = [
        (x_start, y_bot),
        (body_end, y_bot),
        (x_end, mid_y),
        (body_end, y_top),
        (x_start, y_top),
        (x_start, y_bot),
    ]
    codes = [Path.MOVETO] + [Path.LINETO] * 4 + [Path.CLOSEPOLY]
    return Path(verts, codes)


def _assign_stack_levels(results, min_gap, predicate=None):
    """
    Pack visible ORFs onto the minimum number of y-levels so that
    overlapping ones are stacked instead of overlapping visually.

    `predicate(r) -> bool` optionally filters which ORFs participate in
    stacking (others get level 0).

    Returns a {orf_index: level} map and the total level count.
    """
    visible = [r for r in results
               if r.get('visible', True)
               and (predicate is None or predicate(r))]
    visible.sort(key=lambda r: (r.get('start_position', 0),
                                 r.get('end_position', 0)))

    # Each element in `levels` is the right-edge (in bp) currently occupied.
    levels_end = []
    assignment = {}

    for r in visible:
        s = r.get('start_position', 0)
        e = r.get('end_position', s + 3)
        placed = False
        for idx, right in enumerate(levels_end):
            if s >= right + min_gap:
                levels_end[idx] = e
                assignment[r.get('orf_index')] = idx
                placed = True
                break
        if not placed:
            levels_end.append(e)
            assignment[r.get('orf_index')] = len(levels_end) - 1

    return assignment, max(1, len(levels_end))


# =====================================================================
#  Widget
# =====================================================================
class VisualizationWidget(QWidget):

    # Notifies listeners (the results table) that the TIR range filter changed.
    # Arguments: (min_tir, max_tir). max_tir may be math.inf.
    tir_range_changed = pyqtSignal(float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._sequence = ''
        self._results = []
        self._view = 'linear'
        self._highlighted = -1
        self._cached_norm = None
        self._cmap = plt.cm.viridis

        # TIR filter (inclusive)
        self._tir_min = 0.0
        self._tir_max = float('inf')
        # Group-by-shared-stop color mode
        self._group_by_orf = False
        self._group_colors = {}   # orf_index -> hex color (only set when grouping)
        self._has_groups = False

        # Zoom-dependent sequence overlay (populated during _draw_linear)
        self._seq_artists = []     # text objects to clear on each zoom update
        self._linear_ctx = None    # dict with levels / geometry from last draw
        self._ax_linear = None     # current linear axis (for callback lookup)

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(50)
        self._debounce.timeout.connect(self._actual_redraw)

        self._build_ui()

    # ---------------------------------------------------------------- #
    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 4)
        view_lbl = QLabel('View:')
        view_lbl.setFont(QFont('Segoe UI', 9))
        bar.addWidget(view_lbl)

        self.linear_btn = QPushButton('Linear')
        self.circular_btn = QPushButton('Circular')
        for b in (self.linear_btn, self.circular_btn):
            b.setCheckable(True)
            b.setFont(QFont('Segoe UI', 9))
            b.setFixedHeight(26)
            b.setMinimumWidth(70)
        self.linear_btn.setChecked(True)
        self._style_btn(self.linear_btn, True)
        self._style_btn(self.circular_btn, False)
        self.linear_btn.toggled.connect(lambda c: self._toggle_view('linear', c))
        self.circular_btn.toggled.connect(lambda c: self._toggle_view('circular', c))

        grp = QButtonGroup(self)
        grp.setExclusive(True)
        grp.addButton(self.linear_btn)
        grp.addButton(self.circular_btn)

        bar.addWidget(self.linear_btn)
        bar.addWidget(self.circular_btn)

        # --- TIR range filter ---------------------------------------------
        bar.addSpacing(16)
        tir_lbl = QLabel('TIR range:')
        tir_lbl.setFont(QFont('Segoe UI', 9))
        tir_lbl.setToolTip('Only ORFs with TIR inside this range are drawn.')
        bar.addWidget(tir_lbl)

        self.tir_min_spin = self._make_tir_spin(0.0)
        self.tir_min_spin.setToolTip('Minimum TIR (inclusive). 0 = no lower cutoff.')
        self.tir_max_spin = self._make_tir_spin(1e9)
        self.tir_max_spin.setToolTip('Maximum TIR (inclusive). Large value = no upper cutoff.')
        dash = QLabel('–')
        dash.setFont(QFont('Segoe UI', 9))
        bar.addWidget(self.tir_min_spin)
        bar.addWidget(dash)
        bar.addWidget(self.tir_max_spin)

        self.tir_reset_btn = QToolButton()
        self.tir_reset_btn.setText('⟲')
        self.tir_reset_btn.setToolTip('Reset TIR range to full span of results')
        self.tir_reset_btn.setFixedHeight(24)
        self.tir_reset_btn.clicked.connect(self._reset_tir_range)
        bar.addWidget(self.tir_reset_btn)

        self.tir_min_spin.valueChanged.connect(self._on_tir_range_changed)
        self.tir_max_spin.valueChanged.connect(self._on_tir_range_changed)

        # --- Group-by-shared-stop ------------------------------------------
        bar.addSpacing(12)
        self.group_chk = QCheckBox('Color by ORF group')
        self.group_chk.setFont(QFont('Segoe UI', 9))
        self.group_chk.setToolTip(
            'Highlight start codons that share the same stop codon\n'
            '(i.e. alternative starts of the same protein) in distinct colors.\n'
            'Overrides the TIR color gradient.'
        )
        self.group_chk.toggled.connect(self._on_group_toggled)
        bar.addWidget(self.group_chk)

        bar.addStretch()

        self.info_label = QLabel('')
        self.info_label.setFont(QFont('Segoe UI', 8))
        self.info_label.setStyleSheet('color: #555;')
        bar.addWidget(self.info_label)
        root.addLayout(bar)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet('color: #ddd;')
        root.addWidget(sep)

        self.figure = Figure(constrained_layout=False)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding,
                                   QSizePolicy.Policy.Expanding)
        root.addWidget(self.canvas)

        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        self.mpl_toolbar.setStyleSheet('QToolBar { border: none; }')
        root.addWidget(self.mpl_toolbar)

        # Mouse-wheel zoom centered on cursor
        self.canvas.mpl_connect('scroll_event', self._on_scroll)

        self._draw_empty()

    # ---------------------------------------------------------------- #
    def _on_scroll(self, event):
        """Zoom in/out on mouse-wheel rotation, centered on the cursor.
        - Linear view: zoom the x-axis only (leaves vertical layout fixed).
        - Circular view: zoom both axes uniformly (true pinch-style zoom).
        """
        if event.inaxes is None or event.xdata is None:
            return

        base_scale = 1.25
        if event.button == 'up':
            factor = 1.0 / base_scale
        elif event.button == 'down':
            factor = base_scale
        else:
            return

        ax = event.inaxes
        if self._view == 'linear':
            x0, x1 = ax.get_xlim()
            width = x1 - x0
            new_w = max(width * factor, 6.0)         # don't zoom past ~6 bp
            new_w = min(new_w, len(self._sequence) * 1.2)  # cap zoom-out
            rel = (event.xdata - x0) / width if width else 0.5
            new_x0 = event.xdata - new_w * rel
            new_x1 = event.xdata + new_w * (1.0 - rel)
            ax.set_xlim(new_x0, new_x1)
        else:
            if event.ydata is None:
                return
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            wx, wy = x1 - x0, y1 - y0
            # Clamp zoom
            if factor < 1 and wx < 0.05:
                return
            if factor > 1 and wx > 6.0:
                return
            new_wx = wx * factor
            new_wy = wy * factor
            rel_x = (event.xdata - x0) / wx if wx else 0.5
            rel_y = (event.ydata - y0) / wy if wy else 0.5
            ax.set_xlim(event.xdata - new_wx * rel_x,
                        event.xdata + new_wx * (1.0 - rel_x))
            ax.set_ylim(event.ydata - new_wy * rel_y,
                        event.ydata + new_wy * (1.0 - rel_y))
        self.canvas.draw_idle()

    @staticmethod
    def _style_btn(btn, active):
        if active:
            btn.setStyleSheet(
                'QPushButton { background: #1565c0; color: white; '
                'border-radius: 4px; padding: 2px 8px; }'
            )
        else:
            btn.setStyleSheet(
                'QPushButton { background: #e0e0e0; color: #333; '
                'border-radius: 4px; padding: 2px 8px; }'
                'QPushButton:hover { background: #bdbdbd; }'
            )

    @staticmethod
    def _make_tir_spin(default_value):
        sb = QDoubleSpinBox()
        sb.setDecimals(2)
        sb.setRange(0.0, 1e12)
        sb.setValue(default_value)
        sb.setFixedHeight(24)
        sb.setFixedWidth(95)
        sb.setFont(QFont('Segoe UI', 8))
        sb.setAlignment(Qt.AlignmentFlag.AlignRight)
        sb.setKeyboardTracking(False)   # don't fire valueChanged on each keystroke
        return sb

    def _on_tir_range_changed(self, *_):
        lo = self.tir_min_spin.value()
        hi = self.tir_max_spin.value()
        if hi < lo:
            hi = lo
            self.tir_max_spin.blockSignals(True)
            self.tir_max_spin.setValue(hi)
            self.tir_max_spin.blockSignals(False)
        self._tir_min = lo
        self._tir_max = hi if hi > 0 else float('inf')
        self._cached_norm = _compute_norm(self._results, self._tir_min, self._tir_max)
        self._update_info()
        self.tir_range_changed.emit(float(self._tir_min), float(self._tir_max))
        self._schedule_redraw()

    def _reset_tir_range(self):
        if not self._results:
            return
        vals = [r.get('expression') or 0 for r in self._results]
        vals = [v for v in vals if v > 0]
        if not vals:
            return
        lo, hi = min(vals), max(vals)
        self.tir_min_spin.blockSignals(True)
        self.tir_max_spin.blockSignals(True)
        self.tir_min_spin.setValue(0.0)
        self.tir_max_spin.setValue(hi * 1.01)
        self.tir_min_spin.blockSignals(False)
        self.tir_max_spin.blockSignals(False)
        self._tir_min = 0.0
        self._tir_max = float('inf')
        self._cached_norm = _compute_norm(self._results)
        self._update_info()
        self.tir_range_changed.emit(float(self._tir_min), float(self._tir_max))
        self._schedule_redraw()

    def _on_group_toggled(self, checked: bool):
        self._group_by_orf = checked
        self._recompute_groups()
        self._schedule_redraw()

    def _recompute_groups(self):
        """Group ORFs by shared stop codon (same end_position)."""
        by_end = {}
        for r in self._results:
            end = r.get('end_position')
            if end is None:
                continue
            by_end.setdefault(end, []).append(r.get('orf_index'))

        self._group_colors = {}
        color_idx = 0
        self._has_groups = False
        # Deterministic order: by end_position
        for end_pos in sorted(by_end.keys()):
            members = by_end[end_pos]
            if len(members) >= 2:
                col = GROUP_PALETTE[color_idx % len(GROUP_PALETTE)]
                color_idx += 1
                self._has_groups = True
                for oi in members:
                    self._group_colors[oi] = col
            else:
                for oi in members:
                    self._group_colors[oi] = GROUP_SINGLETON_COLOR

    def _passes_tir_filter(self, r) -> bool:
        tir = r.get('expression') or 0
        return self._tir_min <= tir <= self._tir_max

    def _compute_display_ranks(self) -> dict:
        """Return {orf_index: display_rank_1based} for ORFs passing the TIR
        filter, ranked by start position."""
        passing = [r for r in self._results if self._passes_tir_filter(r)]
        passing.sort(key=lambda r: r.get('start_position', 0))
        return {r.get('orf_index'): i + 1 for i, r in enumerate(passing)}

    def _effective_face_color(self, r, tir_norm):
        """Pick fill color for one ORF, respecting group-mode override."""
        if self._group_by_orf:
            return self._group_colors.get(r.get('orf_index'), GROUP_SINGLETON_COLOR)
        tir = r.get('expression') or 0
        return _tir_color(tir, tir_norm, self._cmap)

    def _toggle_view(self, view, checked):
        if not checked:
            return
        self._view = view
        self._style_btn(self.linear_btn, view == 'linear')
        self._style_btn(self.circular_btn, view == 'circular')
        self._schedule_redraw()

    # ---------------------------------------------------------------- #
    #  Public API                                                        #
    # ---------------------------------------------------------------- #
    def set_data(self, sequence, results):
        self._sequence = sequence
        self._results = results
        self._highlighted = -1

        # Reset TIR filter to span observed values
        vals = [r.get('expression') or 0 for r in results]
        vals = [v for v in vals if v > 0]
        hi = max(vals) * 1.01 if vals else 1e6
        self.tir_min_spin.blockSignals(True)
        self.tir_max_spin.blockSignals(True)
        self.tir_min_spin.setValue(0.0)
        self.tir_max_spin.setValue(hi)
        self.tir_min_spin.blockSignals(False)
        self.tir_max_spin.blockSignals(False)
        self._tir_min = 0.0
        self._tir_max = float('inf')

        self._recompute_groups()
        self._cached_norm = _compute_norm(results, self._tir_min, self._tir_max)
        self._update_info()
        self.tir_range_changed.emit(float(self._tir_min), float(self._tir_max))
        self._schedule_redraw()

    def update_visibility(self, orf_index, visible):
        for r in self._results:
            if r.get('orf_index') == orf_index:
                r['visible'] = visible
                break
        self._cached_norm = _compute_norm(
            self._results, self._tir_min, self._tir_max
        )
        self._update_info()
        self._schedule_redraw()

    def highlight_orf(self, orf_index):
        if self._highlighted == orf_index:
            return
        self._highlighted = orf_index
        self._schedule_redraw()

    def redraw(self):
        self._schedule_redraw()

    def _schedule_redraw(self):
        self._debounce.start()

    def _actual_redraw(self):
        if self._view == 'linear':
            self._draw_linear()
        else:
            self._draw_circular()

    def _update_info(self):
        total = len(self._results)
        vis = sum(
            1 for r in self._results
            if r.get('visible', True) and self._passes_tir_filter(r)
        )
        parts = [f'{len(self._sequence)} bp', f'{vis}/{total} ORFs shown']
        if self._tir_min > 0 or self._tir_max != float('inf'):
            parts.append(
                f'TIR {self._tir_min:g}–'
                f'{"∞" if self._tir_max == float("inf") else f"{self._tir_max:g}"}'
            )
        self.info_label.setText('  ·  '.join(parts))

    # ---------------------------------------------------------------- #
    def _draw_empty(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, 'Run OSTIR to see results',
                ha='center', va='center',
                fontsize=12, color='#aaa', transform=ax.transAxes)
        ax.set_axis_off()
        self.canvas.draw_idle()

    # ================================================================ #
    #  LINEAR VIEW                                                       #
    # ================================================================ #
    def _draw_linear(self):
        self.figure.clear()
        if not self._sequence:
            self._draw_empty()
            return

        seq_len = len(self._sequence)
        results = self._results
        norm = self._cached_norm or _compute_norm(results)

        gs = self.figure.add_gridspec(
            1, 2, width_ratios=[60, 1], wspace=0.012,
            left=0.035, right=0.965, top=0.92, bottom=0.14,
        )
        ax = self.figure.add_subplot(gs[0, 0])
        cax = self.figure.add_subplot(gs[0, 1])

        # Stack overlapping ORFs — only those passing the TIR filter
        min_gap = max(seq_len * 0.003, 1.0)
        levels, n_levels = _assign_stack_levels(
            results, min_gap, predicate=self._passes_tir_filter
        )

        arrow_h = 0.7                  # height of one arrow in y-units
        level_step = 0.95              # vertical distance between stacks
        y_backbone = 0.0               # single backbone line

        # Display rank (1..N) for ORFs passing the filter, ordered by start pos
        display_rank = self._compute_display_ranks()

        # --------------------------------------------------------------
        # Backbone
        # --------------------------------------------------------------
        ax.axhline(y_backbone, color=BACKBONE_COLOR, linewidth=1.6, zorder=2)

        # --------------------------------------------------------------
        # ORF arrows (+ RBS pill + start/stop codon caps embedded)
        # --------------------------------------------------------------
        # Coordinate convention: position p (1-based) ↔ x-interval [p, p+1].
        # • An ORF occupies positions [start, end_pos_incl] → x ∈ [start, end_pos_incl + 1]
        # • Start codon = positions [start, start+2] → x ∈ [start, start+3]
        # • Stop  codon = positions [end_pos_incl-2, end_pos_incl] → x ∈ [end_pos_incl-2, end_pos_incl+1]
        arrow_head_bp = max(seq_len * 0.012, 6)

        for r in results:
            visible = r.get('visible', True)
            orf_idx = r.get('orf_index', -1)
            is_hl = (orf_idx == self._highlighted)

            # TIR filter: hide entirely
            if visible and not self._passes_tir_filter(r):
                continue

            start = r.get('start_position', 1)
            # end_pos_incl = last included position (last base of stop codon).
            # OSTIR's end_position is the 1-based position of the first base
            # AFTER the stop codon, so subtract 1.
            end_pos_incl = min(r.get('end_position', start + 3) - 1, seq_len)
            if end_pos_incl <= start:
                end_pos_incl = start + 2
            x_left  = start
            x_right = end_pos_incl + 1
            codon = r.get('start_codon', 'ATG')
            rbs_dist = r.get('RBS_distance_bp', 5) or 0

            level = levels.get(orf_idx, 0) if visible else 0
            y_bot = y_backbone + 0.25 + level * level_step
            y_top = y_bot + arrow_h

            if visible:
                fc = self._effective_face_color(r, norm)
                alpha = 0.95
                edge_col = '#111' if is_hl else '#555'
                edge_lw = 2.0 if is_hl else 0.7
            else:
                fc = '#eeeeee'
                alpha = 0.4
                edge_col = '#bbb'
                edge_lw = 0.5

            # Arrow body — covers the full ORF including stop codon
            path = _arrow_path(x_left, x_right, y_bot, y_top, arrow_head_bp)
            patch = mpatches.PathPatch(
                path, facecolor=fc, edgecolor=edge_col,
                linewidth=edge_lw, alpha=alpha, zorder=4,
                joinstyle='round',
            )
            ax.add_patch(patch)

            if not visible:
                continue

            # Start-codon cap: 3-bp colored band at positions [start, start+2]
            cap_col = CODON_COLORS.get(codon, '#c62828')
            start_cap = mpatches.Rectangle(
                (x_left, y_bot), 3.0, y_top - y_bot,
                facecolor=cap_col, edgecolor='none', alpha=0.95, zorder=5,
            )
            ax.add_patch(start_cap)

            # Stop-codon cap: 3-bp dark band at the last 3 positions
            if r.get('has_stop', True):
                stop_cap_x = x_right - 3
                stop_cap = mpatches.Rectangle(
                    (stop_cap_x, y_bot), 3.0, y_top - y_bot,
                    facecolor=STOP_CAP_COLOR, edgecolor='none',
                    alpha=0.9, zorder=5,
                )
                ax.add_patch(stop_cap)

            # RBS marker — 6-bp yellow pill just above the backbone.
            # SD spans positions [start-rbs_dist-6, start-rbs_dist-1]
            # → x ∈ [start-rbs_dist-6, start-rbs_dist]
            sd_x_right = start - rbs_dist
            sd_x_left = max(1, sd_x_right - 6)
            if sd_x_right > sd_x_left:
                rbs_y = y_backbone + 0.05
                pill = mpatches.FancyBboxPatch(
                    (sd_x_left, rbs_y), sd_x_right - sd_x_left, 0.14,
                    boxstyle='round,pad=0.3',
                    facecolor=RBS_COLOR, edgecolor=RBS_EDGE,
                    linewidth=0.5, alpha=0.95, zorder=3,
                )
                ax.add_patch(pill)
                # faint tie-line from SD right-edge up to the arrow's left edge
                ax.plot([sd_x_right, x_left],
                        [rbs_y + 0.07, (y_bot + y_top) / 2],
                        color=RBS_EDGE, linewidth=0.5,
                        linestyle=':', alpha=0.6, zorder=2)

            # Label — display rank (not raw orf_index).
            # Big arrows get the label inside; small arrows get a callout above
            # so they always show their number.
            rank = display_rank.get(orf_idx, orf_idx + 1)
            label_text = f'#{rank}'
            arrow_width = x_right - x_left
            body_mid_x = (x_left + max(x_right - arrow_head_bp, x_left + 1)) / 2
            center_x = (x_left + x_right) / 2

            if arrow_width > seq_len * 0.035:
                # Fits inside
                ax.text(
                    body_mid_x, (y_bot + y_top) / 2,
                    label_text,
                    ha='center', va='center',
                    fontsize=8, color='white', fontweight='bold', zorder=6,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground='#222')],
                )
            else:
                # Callout above the arrow tip
                callout_y = y_top + 0.18
                ax.plot(
                    [center_x, center_x], [y_top, callout_y - 0.02],
                    color='#555', linewidth=0.5, zorder=5, clip_on=True,
                )
                ax.text(
                    center_x, callout_y,
                    label_text,
                    ha='center', va='bottom',
                    fontsize=7, color='#222', fontweight='bold', zorder=6,
                    bbox=dict(boxstyle='round,pad=0.15',
                               facecolor='white', edgecolor='#555',
                               linewidth=0.4),
                    clip_on=True,
                )

        # --------------------------------------------------------------
        # Ruler at the bottom (x-axis)
        # --------------------------------------------------------------
        ax.set_ylim(-0.55, 0.25 + n_levels * level_step + arrow_h + 0.1)
        margin = seq_len * 0.015
        # Positions are 1..seq_len, and position p occupies x ∈ [p, p+1],
        # so the overall span is 1 .. seq_len+1 in x-coords.
        ax.set_xlim(1 - margin, seq_len + 1 + margin)
        ax.set_yticks([])
        ax.spines[['left', 'right', 'top']].set_visible(False)
        ax.tick_params(axis='x', labelsize=8, colors='#555', length=4)
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel('position (bp)', fontsize=8, color='#555')

        # --------------------------------------------------------------
        # Right side — TIR colorbar, OR ORF-group legend in group mode
        # --------------------------------------------------------------
        if self._group_by_orf:
            self._render_group_legend(cax)
        else:
            sm = ScalarMappable(cmap=self._cmap, norm=norm)
            sm.set_array([])
            cbar = self.figure.colorbar(sm, cax=cax)
            cbar.set_label('log₁₀(TIR)', fontsize=7, color='#555')
            cbar.ax.tick_params(labelsize=6, colors='#555')
            cbar.outline.set_linewidth(0.3)

        # --------------------------------------------------------------
        # Compact legend (codon colors + RBS)
        # --------------------------------------------------------------
        handles = [
            mpatches.Patch(facecolor=CODON_COLORS['ATG'], label='ATG'),
            mpatches.Patch(facecolor=CODON_COLORS['GTG'], label='GTG'),
            mpatches.Patch(facecolor=CODON_COLORS['TTG'], label='TTG'),
            mpatches.Patch(facecolor=STOP_CAP_COLOR, label='stop'),
            mpatches.Patch(facecolor=RBS_COLOR, edgecolor=RBS_EDGE, label='RBS'),
        ]
        ax.legend(
            handles=handles, loc='upper left',
            fontsize=7, frameon=False, ncol=5,
            bbox_to_anchor=(0.0, 1.08),
            handlelength=0.8, handleheight=0.8, columnspacing=1.0,
            labelcolor='#444',
        )

        # --------------------------------------------------------------
        # Zoom-dependent sequence overlay
        # --------------------------------------------------------------
        self._seq_artists = []
        self._linear_ctx = dict(
            levels=levels,
            n_levels=n_levels,
            arrow_h=arrow_h,
            level_step=level_step,
            y_backbone=y_backbone,
            arrow_head_bp=arrow_head_bp,
        )
        self._ax_linear = ax
        ax.callbacks.connect('xlim_changed', self._on_xlim_changed)
        # Initial pass in case view already starts zoomed (usually not)
        self._refresh_sequence_overlay(ax)

        self.canvas.draw_idle()

    # ---------------------------------------------------------------- #
    def _on_xlim_changed(self, ax):
        self._refresh_sequence_overlay(ax)
        self.canvas.draw_idle()

    def _refresh_sequence_overlay(self, ax):
        """Draw nucleotide letters along the backbone and AA letters inside
        each ORF arrow. Each letter sits inside a colored box; font sizes
        adapt to available pixel-space so nothing overlaps."""
        # Clear previous artists (text + background boxes)
        for a in self._seq_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._seq_artists.clear()

        ctx = self._linear_ctx
        seq = self._sequence
        if not ctx or not seq or self._view != 'linear':
            return

        x0, x1 = ax.get_xlim()
        span = max(x1 - x0, 1e-6)
        seq_len = len(seq)

        # Pixels available per base pair — drives font-size & visibility.
        try:
            ax_px = ax.get_window_extent().width
        except Exception:
            ax_px = 1000.0
        px_per_bp = ax_px / span

        # ---------- Nucleotide letters + colored boxes ----------
        # One box = 1 bp wide. Need ~7 px of width to draw a bold letter
        # without crowding — below that we skip the overlay entirely.
        if px_per_bp >= 7.0:
            # Bold letter width ≈ 0.55 × fontsize (pt ≈ px at 96 dpi).
            # Use 80 % of box width for the glyph.
            nt_size = min(16.0, max(8.0, px_per_bp * 0.80 / 0.55))
            y_backbone = ctx['y_backbone']
            # Box sits just below backbone, height scales mildly with font
            box_h = 0.32
            y_box_bot = y_backbone - 0.06 - box_h
            y_center = y_box_bot + box_h / 2

            # Position p (1-based) sits in x-interval [p, p+1].
            # Map x-coords back to the nearest positions to iterate.
            x_lo = max(1, int(math.floor(x0)) - 1)
            x_hi = min(seq_len, int(math.ceil(x1)) + 1)

            for p in range(x_lo, x_hi + 1):
                i = p - 1   # 0-based index into the string
                if i < 0 or i >= seq_len:
                    continue
                nt = seq[i].upper()
                fc = NT_COLORS.get(nt, NT_COLORS['N'])
                fg = NT_FG.get(nt, '#222')
                box = mpatches.Rectangle(
                    (p, y_box_bot), 1.0, box_h,
                    facecolor=fc, edgecolor='white',
                    linewidth=0.5, zorder=8, clip_on=True,
                )
                ax.add_patch(box)
                self._seq_artists.append(box)
                t = ax.text(
                    p + 0.5, y_center, nt,
                    ha='center', va='center',
                    fontsize=nt_size, fontweight='bold',
                    family='DejaVu Sans Mono',
                    color=fg, zorder=9, clip_on=True,
                )
                self._seq_artists.append(t)

        # ---------- Amino-acid letters + neutral boxes inside arrows ----------
        # Codon = 3 bp wide box. Need ~15 px per codon to draw boldly.
        if px_per_bp * 3 >= 15.0:
            aa_size = min(15.0, max(8.0, px_per_bp * 3 * 0.80 / 0.55 * 0.9))
            levels = ctx['levels']
            arrow_h = ctx['arrow_h']
            level_step = ctx['level_step']
            y_backbone = ctx['y_backbone']
            arrow_head_bp = ctx['arrow_head_bp']

            # Box dimensions (inside the arrow body)
            aa_box_h = min(arrow_h * 0.75, 0.55)

            for r in self._results:
                if not r.get('visible', True):
                    continue
                if not self._passes_tir_filter(r):
                    continue
                aa = r.get('aa_sequence') or ''
                if not aa:
                    continue

                start = r.get('start_position', 1)
                end_pos_incl = min(
                    r.get('end_position', start + 3) - 1, seq_len
                )
                x_right = end_pos_incl + 1
                orf_idx = r.get('orf_index', -1)
                level = levels.get(orf_idx, 0)
                y_bot = y_backbone + 0.25 + level * level_step
                y_top = y_bot + arrow_h
                y_mid = (y_bot + y_top) / 2
                y_box_bot = y_mid - aa_box_h / 2

                # AA boxes should stop before the arrow head begins.
                body_end = x_right - arrow_head_bp

                for i, aa_char in enumerate(aa):
                    codon_left = start + i * 3
                    codon_right = codon_left + 3
                    if codon_right < x0 - 3:
                        continue
                    if codon_left > x1 + 3:
                        break
                    # Stop when codon spills into arrow head area
                    if codon_right > body_end:
                        break
                    box = mpatches.FancyBboxPatch(
                        (codon_left + 0.2, y_box_bot),
                        2.6, aa_box_h,
                        boxstyle='round,pad=0.02,rounding_size=0.4',
                        facecolor=AA_BOX_COLOR, edgecolor=AA_BOX_EDGE,
                        linewidth=0.6, alpha=0.92,
                        zorder=7, clip_on=True,
                    )
                    ax.add_patch(box)
                    self._seq_artists.append(box)
                    t = ax.text(
                        codon_left + 1.5, y_mid, aa_char,
                        ha='center', va='center',
                        fontsize=aa_size, fontweight='bold',
                        family='DejaVu Sans Mono',
                        color=AA_TEXT_COLOR, zorder=8, clip_on=True,
                    )
                    self._seq_artists.append(t)

    # ---------------------------------------------------------------- #
    def _render_group_legend(self, cax):
        """Replace the TIR colorbar with a legend of ORF-group colors."""
        cax.clear()
        cax.set_axis_off()

        display_rank = self._compute_display_ranks()

        # Build list of unique group colors + their member ORF #s.
        # Only include members that pass the TIR filter (i.e. are currently drawn).
        groups = {}   # color -> [orf_idx, ...]
        for oi, col in self._group_colors.items():
            if col == GROUP_SINGLETON_COLOR:
                continue
            if oi not in display_rank:
                continue
            groups.setdefault(col, []).append(oi)

        if not groups:
            cax.text(0.5, 0.5, 'No\nshared\nstops',
                     ha='center', va='center', fontsize=7, color='#888',
                     transform=cax.transAxes)
            return

        # Draw swatches top-to-bottom inside cax using axes coords
        # Sort groups by the smallest display-rank they contain
        items = sorted(
            groups.items(),
            key=lambda kv: min(display_rank.get(m, 1_000_000) for m in kv[1]),
        )
        n = len(items)
        cax.text(0.5, 1.02, 'ORF groups',
                 ha='center', va='bottom', fontsize=7, color='#555',
                 transform=cax.transAxes)
        for i, (col, members) in enumerate(items):
            y = 1.0 - (i + 0.5) / n
            cax.add_patch(mpatches.Rectangle(
                (0.05, y - 0.35 / n), 0.28, 0.7 / n,
                facecolor=col, edgecolor='#333', linewidth=0.4,
                transform=cax.transAxes, clip_on=False,
            ))
            ranks = sorted(display_rank.get(m, 0) for m in members)
            label = ','.join(f'#{r}' for r in ranks)
            if len(label) > 14:
                label = label[:12] + '…'
            cax.text(0.38, y, label,
                     ha='left', va='center', fontsize=6.5, color='#333',
                     transform=cax.transAxes)

    # ================================================================ #
    #  CIRCULAR VIEW                                                     #
    # ================================================================ #
    def _draw_circular(self):
        self.figure.clear()
        if not self._sequence:
            self._draw_empty()
            return

        seq_len = len(self._sequence)
        results = self._results
        norm = self._cached_norm or _compute_norm(results)

        gs = self.figure.add_gridspec(
            1, 2, width_ratios=[40, 1], wspace=0.02,
            left=0.02, right=0.96, top=0.98, bottom=0.04,
        )
        ax = self.figure.add_subplot(gs[0, 0], aspect='equal')
        cax = self.figure.add_subplot(gs[0, 1])
        ax.set_axis_off()

        R_BACK_O = 0.995
        R_BACK_I = 0.975

        # Thin backbone ring
        backbone = mpatches.Wedge(
            (0, 0), R_BACK_O, 0, 360,
            width=R_BACK_O - R_BACK_I,
            facecolor=BACKBONE_COLOR, edgecolor='none', alpha=0.9,
        )
        ax.add_patch(backbone)

        # Position ticks
        n_ticks = 12
        tick_step = max(1, seq_len // n_ticks)
        for pos in range(0, seq_len + 1, tick_step):
            a_deg = 90 - (pos / seq_len) * 360
            a_rad = math.radians(a_deg)
            r_in, r_out = R_BACK_O + 0.005, R_BACK_O + 0.04
            ax.plot(
                [r_in * math.cos(a_rad), r_out * math.cos(a_rad)],
                [r_in * math.sin(a_rad), r_out * math.sin(a_rad)],
                '-', color='#888', lw=0.6,
            )
            if pos > 0:
                ax.text(
                    (R_BACK_O + 0.095) * math.cos(a_rad),
                    (R_BACK_O + 0.095) * math.sin(a_rad),
                    str(pos), ha='center', va='center',
                    fontsize=6.5, color='#555',
                )

        # Stack overlapping ORFs (polar) — only those passing the TIR filter
        min_gap = max(seq_len * 0.003, 1.0)
        levels, n_levels = _assign_stack_levels(
            results, min_gap, predicate=self._passes_tir_filter
        )
        n_levels = max(n_levels, 1)

        # Each ORF ring is this thick; stacks grow inward
        ring_thickness = 0.12
        ring_gap = 0.015
        first_ring_outer = R_BACK_I - 0.02

        # Display ranks for the circular labels
        display_rank = self._compute_display_ranks()

        # Helper to convert a 1-based position boundary to a polar angle.
        # Position p occupies x-interval [p, p+1]; its *start* boundary is at p,
        # so the angle is 90° - (p - 1)/seq_len * 360°.
        def pos_angle(p):
            return 90 - ((p - 1) / seq_len) * 360

        for r in results:
            visible = r.get('visible', True)
            orf_idx = r.get('orf_index', -1)
            is_hl = (orf_idx == self._highlighted)

            # TIR filter: hide entirely
            if visible and not self._passes_tir_filter(r):
                continue

            start = r.get('start_position', 1)
            end_pos_incl = min(
                r.get('end_position', start + 3) - 1, seq_len
            )
            if end_pos_incl <= start:
                end_pos_incl = start + 2
            codon = r.get('start_codon', 'ATG')
            rbs_dist = r.get('RBS_distance_bp', 5) or 0

            level = levels.get(orf_idx, 0) if visible else 0
            r_out = first_ring_outer - level * (ring_thickness + ring_gap)
            r_in = r_out - ring_thickness

            if visible:
                fc = self._effective_face_color(r, norm)
                alpha = 0.92
                edge_col = '#111' if is_hl else 'white'
                edge_lw = 1.6 if is_hl else 0.5
            else:
                fc = '#eeeeee'
                alpha = 0.3
                edge_col = '#ccc'
                edge_lw = 0.3

            # Arc covers full ORF: positions [start, end_pos_incl]
            theta_start = pos_angle(start)          # leading (5') edge
            theta_end   = pos_angle(end_pos_incl + 1)  # trailing (3') edge
            arc = mpatches.Wedge(
                (0, 0), r_out, theta_end, theta_start,
                width=r_out - r_in,
                facecolor=fc, edgecolor=edge_col,
                linewidth=edge_lw, alpha=alpha, zorder=4,
            )
            ax.add_patch(arc)

            if not visible:
                continue

            # Start-codon colored cap: 3 bp at positions [start, start+2]
            t_cap_a = pos_angle(start + 3)   # trailing angle of the cap
            t_cap_b = pos_angle(start)       # leading angle
            cap = mpatches.Wedge(
                (0, 0), r_out, t_cap_a, t_cap_b,
                width=r_out - r_in,
                facecolor=CODON_COLORS.get(codon, '#c62828'),
                edgecolor='none', alpha=0.95, zorder=5,
            )
            ax.add_patch(cap)

            # Stop-codon dark cap: last 3 bp at positions [end-2, end]
            if r.get('has_stop', True):
                t_stop_a = pos_angle(end_pos_incl + 1)
                t_stop_b = pos_angle(end_pos_incl - 2)
                stop_cap = mpatches.Wedge(
                    (0, 0), r_out, t_stop_a, t_stop_b,
                    width=r_out - r_in,
                    facecolor=STOP_CAP_COLOR, edgecolor='none',
                    alpha=0.9, zorder=5,
                )
                ax.add_patch(stop_cap)

            # RBS marker — 6-bp yellow arc just outside the backbone
            sd_pos_right = start - rbs_dist - 1  # last SD position
            sd_pos_left  = sd_pos_right - 5      # first SD position
            if sd_pos_left >= 1:
                t_sd_a = pos_angle(sd_pos_right + 1)
                t_sd_b = pos_angle(sd_pos_left)
                rbs_r_out = R_BACK_I - 0.003
                rbs_arc = mpatches.Wedge(
                    (0, 0), rbs_r_out, t_sd_a, t_sd_b,
                    width=0.02,
                    facecolor=RBS_COLOR, edgecolor=RBS_EDGE,
                    linewidth=0.3, alpha=0.95, zorder=6,
                )
                ax.add_patch(rbs_arc)

            # Label — always drawn. Big arcs get it inside; small arcs get
            # a callout just outside the ring so every ORF is identifiable.
            rank = display_rank.get(orf_idx, orf_idx + 1)
            orf_frac = (end_pos_incl - start + 1) / seq_len
            mid_p = (start + end_pos_incl + 1) / 2
            mid_deg = pos_angle(mid_p)
            mid_rad = math.radians(mid_deg)

            if orf_frac > 0.03:
                r_label = (r_out + r_in) / 2
                ax.text(
                    r_label * math.cos(mid_rad), r_label * math.sin(mid_rad),
                    f'#{rank}', ha='center', va='center',
                    fontsize=6.5, color='white', fontweight='bold', zorder=7,
                    path_effects=[pe.withStroke(linewidth=1.4, foreground='#222')],
                )
            else:
                # Callout: a thin radial line + label sitting outside backbone
                r_inner = r_out
                r_outer = R_BACK_O + 0.12
                ax.plot(
                    [r_inner * math.cos(mid_rad), r_outer * math.cos(mid_rad)],
                    [r_inner * math.sin(mid_rad), r_outer * math.sin(mid_rad)],
                    '-', color='#666', lw=0.4, zorder=6,
                )
                ax.text(
                    (r_outer + 0.04) * math.cos(mid_rad),
                    (r_outer + 0.04) * math.sin(mid_rad),
                    f'#{rank}', ha='center', va='center',
                    fontsize=6, color='#222', fontweight='bold', zorder=7,
                    bbox=dict(boxstyle='round,pad=0.12',
                               facecolor='white', edgecolor='#666',
                               linewidth=0.35),
                )

        # Centre label
        ax.text(0, 0.08, f'{seq_len} bp',
                ha='center', va='center',
                fontsize=11, color='#37474f', fontweight='bold')
        vis_n = sum(
            1 for r in results
            if r.get('visible', True) and self._passes_tir_filter(r)
        )
        ax.text(0, -0.07, f'{vis_n} ORF{"s" if vis_n != 1 else ""} shown',
                ha='center', va='center',
                fontsize=8, color='#888')

        # Right side — TIR colorbar, OR ORF-group legend in group mode
        if self._group_by_orf:
            self._render_group_legend(cax)
        else:
            sm = ScalarMappable(cmap=self._cmap, norm=norm)
            sm.set_array([])
            cbar = self.figure.colorbar(sm, cax=cax)
            cbar.set_label('log₁₀(TIR)', fontsize=7, color='#555')
            cbar.ax.tick_params(labelsize=6, colors='#555')
            cbar.outline.set_linewidth(0.3)

        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.25, 1.25)
        self.canvas.draw_idle()

    # ---------------------------------------------------------------- #
    def save_figure(self, path):
        self.figure.savefig(path, dpi=200, bbox_inches='tight')
