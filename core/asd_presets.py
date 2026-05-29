"""
Anti-Shine-Dalgarno (anti-SD) preset management.

Stores a small library of named anti-SD sequences (the 3' rRNA tail OSTIR
folds the mRNA SD against). Two presets ship by default:

  * 'E. coli (native)'   — ACCTCCTTA, the canonical E. coli 16S 3' end.
  * 'pET T7 (derived)'   — derived from the pET/T7 RBS leader by taking the
                           reverse complement of its Shine-Dalgarno core.

Custom presets are persisted to ``asd_presets.json`` next to the project so
they survive between sessions.

Also provides:
  * reverse_complement(seq)
  * derive_asd_from_rbs(rbs)  — turn an mRNA RBS/leader into the anti-SD that
                                pairs with its SD core.
"""

import os
import re
import json
from typing import Dict

_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A', 'N': 'N'}

# Canonical Shine-Dalgarno consensus used to locate the SD core in a leader.
_SD_CONSENSUS = 'AGGAGG'

DEFAULT_PRESETS: Dict[str, str] = {
    'E. coli (native)': 'ACCTCCTTA',
    # Reverse complement of the pET/T7 leader's SD core (AAGAAGGAG):
    'pET T7 (derived)': 'CTCCTTCTT',
}

NATIVE_NAME = 'E. coli (native)'
NATIVE_ASD = DEFAULT_PRESETS[NATIVE_NAME]

_STORE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'asd_presets.json',
)


# ---------------------------------------------------------------------- #
def reverse_complement(seq: str) -> str:
    seq = re.sub(r'[^ACGTUacgtu]', '', seq).upper()
    return ''.join(_COMP.get(b, 'N') for b in reversed(seq))


def clean_nt(seq: str) -> str:
    return re.sub(r'[^ACGTU]', '', (seq or '').upper()).replace('U', 'T')


def derive_asd_from_rbs(rbs: str, length: int = 9) -> str:
    """Derive an anti-SD from an mRNA RBS/leader.

    Finds the window best matching the SD consensus (AGGAGG), grabs a
    ``length``-nt window around it and returns its reverse complement — i.e.
    the rRNA tail that would base-pair with that SD.
    """
    s = clean_nt(rbs)
    if len(s) < 6:
        return ''
    best_i, best_score = 0, -1
    for i in range(0, len(s) - 5):
        w = s[i:i + 6]
        score = sum(1 for a, b in zip(w, _SD_CONSENSUS) if a == b)
        if score > best_score:
            best_score, best_i = score, i
    start = max(0, best_i - 1)
    window = s[start:start + length]
    if len(window) < length:                      # clamp to the 3' end
        window = s[max(0, len(s) - length):]
    return reverse_complement(window)


# ---------------------------------------------------------------------- #
def load_presets() -> Dict[str, str]:
    """Return the preset dict (defaults + any persisted custom entries)."""
    presets = dict(DEFAULT_PRESETS)
    try:
        if os.path.exists(_STORE_PATH):
            with open(_STORE_PATH, 'r', encoding='utf-8') as fh:
                user = json.load(fh)
            if isinstance(user, dict):
                for k, v in user.items():
                    if isinstance(v, str) and clean_nt(v):
                        presets[str(k)] = clean_nt(v)
    except Exception:
        pass
    return presets


def save_preset(name: str, seq: str) -> Dict[str, str]:
    """Add/overwrite a custom preset and persist. Returns the full dict."""
    name = (name or '').strip()
    seq = clean_nt(seq)
    if not name or not seq:
        return load_presets()
    # Read current custom store (defaults are implicit, not written).
    custom = {}
    try:
        if os.path.exists(_STORE_PATH):
            with open(_STORE_PATH, 'r', encoding='utf-8') as fh:
                data = json.load(fh)
            if isinstance(data, dict):
                custom = {str(k): str(v) for k, v in data.items()}
    except Exception:
        custom = {}
    custom[name] = seq
    try:
        with open(_STORE_PATH, 'w', encoding='utf-8') as fh:
            json.dump(custom, fh, indent=2)
    except Exception:
        pass
    return load_presets()


def delete_preset(name: str) -> Dict[str, str]:
    """Remove a custom preset (defaults can't be deleted)."""
    if name in DEFAULT_PRESETS:
        return load_presets()
    try:
        if os.path.exists(_STORE_PATH):
            with open(_STORE_PATH, 'r', encoding='utf-8') as fh:
                data = json.load(fh)
            if isinstance(data, dict) and name in data:
                del data[name]
                with open(_STORE_PATH, 'w', encoding='utf-8') as fh:
                    json.dump(data, fh, indent=2)
    except Exception:
        pass
    return load_presets()
