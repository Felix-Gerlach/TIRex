@echo off
REM ─────────────────────────────────────────────────────────────────────────
REM  TIRex  –  PyInstaller build script
REM  Run this from inside the TIRex\ folder:
REM      build.bat
REM
REM  Output will be in:  dist\TIRex\   (one-folder bundle)
REM  The .exe is at:     dist\TIRex\TIRex.exe
REM
REM  IMPORTANT: Before building, make sure ViennaRNA is installed and that
REM  "RNAfold.exe" is reachable. The build bundles the Python code only.
REM  The end user must still have ViennaRNA in PATH when running the .exe.
REM ─────────────────────────────────────────────────────────────────────────

echo [BUILD] Installing / upgrading PyInstaller...
pip install --upgrade pyinstaller

echo [BUILD] Running PyInstaller...
pyinstaller ^
    --name "TIRex" ^
    --windowed ^
    --onedir ^
    --noconfirm ^
    --clean ^
    --add-data "ui;ui" ^
    --add-data "core;core" ^
    --hidden-import "ostir" ^
    --hidden-import "ostir.ostir" ^
    --hidden-import "ostir.ostir_factory" ^
    --hidden-import "ostir.ViennaRNA" ^
    --hidden-import "ostir.shortcuts" ^
    --hidden-import "Bio" ^
    --hidden-import "Bio.SeqUtils" ^
    --hidden-import "Bio.SeqUtils.ProtParam" ^
    --hidden-import "pandas" ^
    --hidden-import "matplotlib" ^
    --hidden-import "matplotlib.backends.backend_qtagg" ^
    --hidden-import "PyQt6.QtSvg" ^
    --hidden-import "PyQt6.QtPrintSupport" ^
    main.py

echo.
if exist "dist\TIRex\TIRex.exe" (
    echo [BUILD] SUCCESS!
    echo [BUILD] Executable: dist\TIRex\TIRex.exe
) else (
    echo [BUILD] Build may have failed. Check output above.
)
pause
