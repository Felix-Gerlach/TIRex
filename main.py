"""
TIRex — entry point.

Run with:
    python main.py

Package to .exe with:
    build.bat
"""

import sys
import os

# Make sure the project root is on the Python path (needed for PyInstaller builds)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication, QMessageBox
from PyQt6.QtCore import Qt


def check_viennarna() -> list:
    """Return list of missing ViennaRNA CLI tools."""
    from shutil import which
    return [t for t in ('RNAfold', 'RNAsubopt', 'RNAeval') if which(t) is None]


def main():
    # High-DPI support
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )

    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setApplicationName('TIRex')
    app.setApplicationVersion('1.0.0')

    # Dependency check
    missing = check_viennarna()
    if missing:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Icon.Critical)
        msg.setWindowTitle('Missing ViennaRNA')
        msg.setText(
            'The following ViennaRNA executables were not found in PATH:\n\n'
            + '\n'.join(f'  • {t}' for t in missing)
            + '\n\nPlease install ViennaRNA and add its folder to your '
              'system PATH, then restart the application.\n\n'
              'Download: https://www.tbi.univie.ac.at/RNA/#download'
        )
        msg.exec()
        return 1

    # Lazy import after path is set
    from ui.main_window import MainWindow
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
