"""
Centralised visual theme for TIRex.

A single light, flat, modern stylesheet applied app-wide. Colours live in
PALETTE; the QSS template references them via @TOKEN@ placeholders so the whole
look can be retuned from one place.

Usage (in main.py):
    from ui.theme import apply_theme
    apply_theme(app)
"""

from PyQt6.QtGui import QPalette, QColor, QFont

# ---------------------------------------------------------------------------
#  Palette
# ---------------------------------------------------------------------------
PALETTE = {
    # surfaces
    'BG_APP':       '#eef2f7',   # window background
    'BG_SURFACE':   '#ffffff',   # cards / inputs
    'BG_SUBTLE':    '#f5f8fc',   # subtle fills, alt rows
    'BG_HEADER':    '#0f172a',   # dark header band
    'BG_ELEVATED':  '#fbfdff',

    # borders
    'BORDER':       '#dbe3ee',
    'BORDER_STRONG':'#c2cedd',

    # text
    'TEXT':         '#1f2a37',
    'TEXT_MUTED':   '#64748b',
    'TEXT_FAINT':   '#94a3b8',
    'TEXT_ON_DARK': '#e2e8f0',

    # brand / accents
    'PRIMARY':      '#3b82f6',
    'PRIMARY_HOVER':'#2563eb',
    'PRIMARY_PRESS':'#1d4ed8',
    'PRIMARY_SOFT': '#e8f1ff',

    'SUCCESS':      '#16a34a',
    'SUCCESS_HOVER':'#15803d',
    'SUCCESS_PRESS':'#166534',

    'ACCENT':       '#0d9488',
    'DANGER':       '#e11d48',
    'WARNING':      '#f59e0b',

    # selection
    'SEL_BG':       '#dbeafe',
    'SEL_FG':       '#1e3a8a',
}


# ---------------------------------------------------------------------------
#  Stylesheet template (@TOKEN@ replaced from PALETTE)
# ---------------------------------------------------------------------------
_QSS = """
* {
    font-family: "Segoe UI", "Inter", system-ui, sans-serif;
    font-size: 10pt;
    color: @TEXT@;
    outline: none;
}

QMainWindow, QDialog {
    background: @BG_APP@;
}

QToolTip {
    background: @BG_HEADER@;
    color: @TEXT_ON_DARK@;
    border: none;
    border-radius: 6px;
    padding: 6px 9px;
    font-size: 9pt;
}

/* ---------------- Group boxes (cards) ---------------- */
QGroupBox {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER@;
    border-radius: 10px;
    margin-top: 16px;
    padding: 10px 10px 10px 10px;
    font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    top: 2px;
    padding: 2px 8px;
    color: @PRIMARY_PRESS@;
    background: @PRIMARY_SOFT@;
    border-radius: 6px;
    font-size: 9pt;
    font-weight: 700;
}

/* ---------------- Labels ---------------- */
QLabel { background: transparent; }
QLabel[muted="true"] { color: @TEXT_MUTED@; }

/* ---------------- Buttons ---------------- */
QPushButton {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER_STRONG@;
    border-radius: 8px;
    padding: 6px 14px;
    color: @TEXT@;
    font-weight: 600;
}
QPushButton:hover {
    background: @BG_SUBTLE@;
    border-color: @PRIMARY@;
}
QPushButton:pressed {
    background: @PRIMARY_SOFT@;
}
QPushButton:disabled {
    background: @BG_SUBTLE@;
    color: @TEXT_FAINT@;
    border-color: @BORDER@;
}

/* Primary call-to-action (e.g. Run) */
QPushButton#PrimaryButton {
    background: @SUCCESS@;
    border: none;
    color: #ffffff;
    padding: 9px 16px;
    font-size: 10.5pt;
    font-weight: 700;
}
QPushButton#PrimaryButton:hover  { background: @SUCCESS_HOVER@; }
QPushButton#PrimaryButton:pressed{ background: @SUCCESS_PRESS@; }
QPushButton#PrimaryButton:disabled { background: #9bb8a6; color: #eef5ef; }

/* Accent button (blue) */
QPushButton#AccentButton {
    background: @PRIMARY@;
    border: none;
    color: #ffffff;
    font-weight: 700;
}
QPushButton#AccentButton:hover  { background: @PRIMARY_HOVER@; }
QPushButton#AccentButton:pressed{ background: @PRIMARY_PRESS@; }

/* Ghost / subtle button */
QPushButton#GhostButton {
    background: transparent;
    border: 1px solid @BORDER@;
    color: @TEXT_MUTED@;
    font-weight: 600;
}
QPushButton#GhostButton:hover {
    background: @BG_SUBTLE@;
    color: @TEXT@;
    border-color: @BORDER_STRONG@;
}

/* ---------------- Text inputs ---------------- */
QLineEdit, QTextEdit, QPlainTextEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER_STRONG@;
    border-radius: 8px;
    padding: 5px 8px;
    selection-background-color: @SEL_BG@;
    selection-color: @SEL_FG@;
}
QLineEdit:focus, QTextEdit:focus, QPlainTextEdit:focus,
QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border: 1px solid @PRIMARY@;
    background: @BG_ELEVATED@;
}
QLineEdit:disabled, QTextEdit:disabled, QPlainTextEdit:disabled {
    background: @BG_SUBTLE@;
    color: @TEXT_FAINT@;
}
QLineEdit::placeholder { color: @TEXT_FAINT@; }

/* ComboBox */
QComboBox::drop-down {
    border: none;
    width: 22px;
    subcontrol-origin: padding;
    subcontrol-position: center right;
}
QComboBox::down-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 6px solid @TEXT_MUTED@;
    margin-right: 8px;
}
QComboBox QAbstractItemView {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER_STRONG@;
    border-radius: 8px;
    padding: 4px;
    selection-background-color: @PRIMARY_SOFT@;
    selection-color: @PRIMARY_PRESS@;
}

/* SpinBox buttons */
QSpinBox::up-button, QDoubleSpinBox::up-button,
QSpinBox::down-button, QDoubleSpinBox::down-button {
    width: 16px;
    border: none;
    background: transparent;
}
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-bottom: 5px solid @TEXT_MUTED@;
}
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid @TEXT_MUTED@;
}

/* ---------------- Check / radio ---------------- */
QCheckBox, QRadioButton { spacing: 7px; background: transparent; }
QCheckBox::indicator, QRadioButton::indicator {
    width: 17px; height: 17px;
    border: 1px solid @BORDER_STRONG@;
    background: @BG_SURFACE@;
}
QCheckBox::indicator { border-radius: 5px; }
QRadioButton::indicator { border-radius: 9px; }
QCheckBox::indicator:hover, QRadioButton::indicator:hover {
    border-color: @PRIMARY@;
}
QCheckBox::indicator:checked {
    background: @PRIMARY@;
    border-color: @PRIMARY@;
    image: none;
}
QRadioButton::indicator:checked {
    background: qradialgradient(cx:0.5, cy:0.5, radius:0.5,
                fx:0.5, fy:0.5, stop:0 #ffffff, stop:0.45 #ffffff,
                stop:0.5 @PRIMARY@, stop:1 @PRIMARY@);
    border-color: @PRIMARY@;
}

/* ---------------- Tables ---------------- */
QTableWidget, QTableView {
    background: @BG_SURFACE@;
    alternate-background-color: @BG_SUBTLE@;
    gridline-color: @BORDER@;
    border: 1px solid @BORDER@;
    border-radius: 10px;
    selection-background-color: @SEL_BG@;
    selection-color: @SEL_FG@;
}
QTableView::item { padding: 3px 4px; }
QTableView::item:selected { background: @SEL_BG@; color: @SEL_FG@; }
QHeaderView::section {
    background: @BG_SUBTLE@;
    color: @TEXT_MUTED@;
    padding: 6px 8px;
    border: none;
    border-right: 1px solid @BORDER@;
    border-bottom: 1px solid @BORDER_STRONG@;
    font-weight: 700;
    font-size: 9pt;
}
QHeaderView::section:hover { color: @PRIMARY_PRESS@; }
QTableCornerButton::section { background: @BG_SUBTLE@; border: none; }

/* ---------------- Tabs ---------------- */
QTabWidget::pane {
    border: 1px solid @BORDER@;
    border-radius: 10px;
    top: -1px;
    background: @BG_SURFACE@;
}
QTabBar::tab {
    background: transparent;
    color: @TEXT_MUTED@;
    padding: 7px 16px;
    margin-right: 4px;
    border: none;
    border-bottom: 2px solid transparent;
    font-weight: 600;
}
QTabBar::tab:hover { color: @TEXT@; }
QTabBar::tab:selected {
    color: @PRIMARY_PRESS@;
    border-bottom: 2px solid @PRIMARY@;
}

/* ---------------- Scrollbars ---------------- */
QScrollBar:vertical {
    background: transparent; width: 11px; margin: 2px;
}
QScrollBar::handle:vertical {
    background: @BORDER_STRONG@; border-radius: 5px; min-height: 28px;
}
QScrollBar::handle:vertical:hover { background: @TEXT_FAINT@; }
QScrollBar:horizontal {
    background: transparent; height: 11px; margin: 2px;
}
QScrollBar::handle:horizontal {
    background: @BORDER_STRONG@; border-radius: 5px; min-width: 28px;
}
QScrollBar::handle:horizontal:hover { background: @TEXT_FAINT@; }
QScrollBar::add-line, QScrollBar::sub-line { height: 0; width: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }

/* ---------------- Progress ---------------- */
QProgressBar {
    background: @BG_SUBTLE@;
    border: 1px solid @BORDER@;
    border-radius: 7px;
    text-align: center;
    color: @TEXT_MUTED@;
    font-size: 8pt;
}
QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                stop:0 @PRIMARY@, stop:1 @ACCENT@);
    border-radius: 6px;
}

/* ---------------- Sliders ---------------- */
QSlider::groove:horizontal {
    height: 5px; background: @BORDER@; border-radius: 3px;
}
QSlider::sub-page:horizontal {
    background: @PRIMARY@; border-radius: 3px;
}
QSlider::handle:horizontal {
    background: @BG_SURFACE@; border: 2px solid @PRIMARY@;
    width: 15px; height: 15px; margin: -6px 0; border-radius: 8px;
}
QSlider::handle:horizontal:hover { background: @PRIMARY_SOFT@; }

/* ---------------- Menu bar / menus ---------------- */
QMenuBar {
    background: @BG_SURFACE@;
    border-bottom: 1px solid @BORDER@;
    padding: 2px 4px;
}
QMenuBar::item {
    background: transparent; padding: 6px 12px; border-radius: 6px;
    color: @TEXT@;
}
QMenuBar::item:selected { background: @PRIMARY_SOFT@; color: @PRIMARY_PRESS@; }
QMenu {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER_STRONG@;
    border-radius: 8px;
    padding: 5px;
}
QMenu::item { padding: 6px 22px 6px 14px; border-radius: 6px; }
QMenu::item:selected { background: @PRIMARY_SOFT@; color: @PRIMARY_PRESS@; }
QMenu::separator { height: 1px; background: @BORDER@; margin: 5px 8px; }

/* ---------------- Status bar ---------------- */
QStatusBar {
    background: @BG_SURFACE@;
    border-top: 1px solid @BORDER@;
    color: @TEXT_MUTED@;
}
QStatusBar::item { border: none; }

/* ---------------- Splitter ---------------- */
QSplitter::handle { background: transparent; }
QSplitter::handle:horizontal { width: 6px; }
QSplitter::handle:vertical   { height: 6px; }
QSplitter::handle:hover { background: @PRIMARY_SOFT@; }

/* ---------------- Lists ---------------- */
QListWidget {
    background: @BG_SURFACE@;
    border: 1px solid @BORDER@;
    border-radius: 8px;
    padding: 3px;
}
QListWidget::item { padding: 5px 7px; border-radius: 6px; }
QListWidget::item:selected { background: @PRIMARY_SOFT@; color: @PRIMARY_PRESS@; }
QListWidget::item:hover { background: @BG_SUBTLE@; }
"""


def build_stylesheet() -> str:
    qss = _QSS
    for token, value in PALETTE.items():
        qss = qss.replace(f'@{token}@', value)
    return qss


def apply_theme(app):
    """Apply the TIRex theme to a QApplication."""
    app.setStyle('Fusion')

    # Base palette so native-drawn bits match the QSS.
    pal = app.palette()
    pal.setColor(QPalette.ColorRole.Window, QColor(PALETTE['BG_APP']))
    pal.setColor(QPalette.ColorRole.Base, QColor(PALETTE['BG_SURFACE']))
    pal.setColor(QPalette.ColorRole.AlternateBase, QColor(PALETTE['BG_SUBTLE']))
    pal.setColor(QPalette.ColorRole.Text, QColor(PALETTE['TEXT']))
    pal.setColor(QPalette.ColorRole.WindowText, QColor(PALETTE['TEXT']))
    pal.setColor(QPalette.ColorRole.Button, QColor(PALETTE['BG_SURFACE']))
    pal.setColor(QPalette.ColorRole.ButtonText, QColor(PALETTE['TEXT']))
    pal.setColor(QPalette.ColorRole.Highlight, QColor(PALETTE['SEL_BG']))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor(PALETTE['SEL_FG']))
    pal.setColor(QPalette.ColorRole.ToolTipBase, QColor(PALETTE['BG_HEADER']))
    pal.setColor(QPalette.ColorRole.ToolTipText, QColor(PALETTE['TEXT_ON_DARK']))
    app.setPalette(pal)

    app.setFont(QFont('Segoe UI', 10))
    app.setStyleSheet(build_stylesheet())
