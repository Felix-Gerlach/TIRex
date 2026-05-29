# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[('ui', 'ui'), ('core', 'core')],
    hiddenimports=['ostir', 'ostir.ostir', 'ostir.ostir_factory', 'ostir.ViennaRNA', 'ostir.shortcuts', 'Bio', 'Bio.SeqUtils', 'Bio.SeqUtils.ProtParam', 'pandas', 'matplotlib', 'matplotlib.backends.backend_qtagg', 'PyQt6.QtSvg', 'PyQt6.QtPrintSupport'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='TIRex',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
