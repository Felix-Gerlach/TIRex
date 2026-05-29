"""
Session save/load for TIRex (.tirex = JSON).

Persists everything needed to restore a working session: the input sequence,
the run parameters (incl. both anti-SDs), the enriched OSTIR results, the
target ORF, and TIR-range filter. Results may contain numpy/None values, so we
sanitise to JSON-native types on save.
"""

import json
from typing import Any, Dict

FORMAT_VERSION = 1


def _jsonable(v: Any):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, (str, int, bool)) or v is None:
        return v
    try:
        f = float(v)
        return f
    except (TypeError, ValueError):
        return str(v)


def save_session(path: str, *, sequence: str, params: Dict,
                 results: list, target_orf_index=None,
                 tir_min: float = 0.0, tir_max: float = 0.0) -> None:
    data = {
        'format': 'tirex-session',
        'version': FORMAT_VERSION,
        'sequence': sequence or '',
        'params': _jsonable(params or {}),
        'results': _jsonable(results or []),
        'target_orf_index': target_orf_index,
        'tir_min': float(tir_min or 0.0),
        'tir_max': float(tir_max or 0.0),
    }
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(data, fh, indent=2)


def load_session(path: str) -> Dict:
    with open(path, 'r', encoding='utf-8') as fh:
        data = json.load(fh)
    if data.get('format') != 'tirex-session':
        raise ValueError('Not a TIRex session file.')
    return data
