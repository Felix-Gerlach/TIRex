"""
Sequence-edit constraint checks (soft-flag).

Given a proposed mutated sequence and the reference it came from, return a list
of human-readable warning strings. Nothing is rejected here — the tuner shows
the warnings in a column and lets the user decide.

Checks
------
  * newly-created restriction sites (default common cloning enzymes, editable)
  * newly-created start codons (ATG/GTG/TTG)  — competing initiation
  * newly-created in-frame stop codons         — only meaningful with a frame
  * newly-created Shine-Dalgarno cores (AGGAGG-like)
  * homopolymer runs introduced (default > 4)
  * GC content outside a window in the edited region
"""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# Common cloning enzymes (name -> recognition site, IUPAC expanded to regex).
DEFAULT_ENZYMES: Dict[str, str] = {
    'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG', 'NdeI': 'CATATG', 'NotI': 'GCGGCCGC',
    'XbaI': 'TCTAGA', 'SpeI': 'ACTAGT', 'PstI': 'CTGCAG',
    'SalI': 'GTCGAC', 'KpnI': 'GGTACC', 'SacI': 'GAGCTC',
    'NcoI': 'CCATGG', 'BglII': 'AGATCT', 'NheI': 'GCTAGC',
}

_START_CODONS = ('ATG', 'GTG', 'TTG')
_STOP_CODONS = ('TAA', 'TAG', 'TGA')
_SD_RE = re.compile(r'AGGAGG|GGAGG|AGGAG')   # SD-core-ish


@dataclass
class ConstraintConfig:
    enzymes: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_ENZYMES))
    check_restriction: bool = True
    check_new_start: bool = True
    check_new_stop: bool = True
    check_new_sd: bool = True
    check_homopolymer: bool = True
    homopolymer_max: int = 4
    check_gc: bool = True
    gc_low: float = 25.0
    gc_high: float = 75.0


def _gc(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in 'GC')
    return 100.0 * gc / len(seq)


def _longest_run(seq: str) -> int:
    best = run = 0
    prev = ''
    for b in seq:
        run = run + 1 if b == prev else 1
        prev = b
        best = max(best, run)
    return best


def _edit_span(ref: str, new: str):
    """Return (lo, hi) index span in `new` covering the changed region, with a
    little padding. Robust to indels via a common-prefix/suffix diff."""
    if ref == new:
        return None
    n = len(new)
    # common prefix
    p = 0
    while p < min(len(ref), len(new)) and ref[p] == new[p]:
        p += 1
    # common suffix
    s = 0
    while (s < min(len(ref), len(new)) - p
           and ref[len(ref) - 1 - s] == new[len(new) - 1 - s]):
        s += 1
    lo = max(0, p - 2)
    hi = min(n, n - s + 2)
    return lo, hi


def find_sites(seq: str, enzymes: Dict[str, str]) -> Dict[str, int]:
    """Count occurrences of each enzyme site in seq."""
    counts = {}
    for name, site in enzymes.items():
        c = len(re.findall(f'(?={site})', seq))
        if c:
            counts[name] = c
    return counts


def analyze(new_seq: str, ref_seq: str,
            frame_start: Optional[int] = None,
            cfg: Optional[ConstraintConfig] = None) -> List[str]:
    """Return a list of warning strings for `new_seq` relative to `ref_seq`.

    `frame_start` (1-based) enables in-frame stop-codon detection.
    """
    cfg = cfg or ConstraintConfig()
    warnings: List[str] = []
    if not new_seq or new_seq == ref_seq:
        return warnings

    span = _edit_span(ref_seq, new_seq) or (0, len(new_seq))
    lo, hi = span
    win = new_seq[lo:hi]

    # --- newly created restriction sites ---------------------------------
    if cfg.check_restriction:
        ref_sites = find_sites(ref_seq, cfg.enzymes)
        new_sites = find_sites(new_seq, cfg.enzymes)
        created = [name for name, c in new_sites.items()
                   if c > ref_sites.get(name, 0)]
        for name in created:
            warnings.append(f'creates {name} site')

    # --- newly created start codons --------------------------------------
    if cfg.check_new_start:
        # Look a few nt around the edit for a start triplet not present before.
        a = max(0, lo - 2)
        region_new = new_seq[a:hi + 2]
        # Map region back to ref by alignment is fiddly with indels; use a
        # simple heuristic: flag if a start codon now appears in the window
        # but the same-length ref window had fewer.
        ref_region = ref_seq[a:a + len(region_new)]
        def _count_starts(s):
            return sum(s[i:i + 3] in _START_CODONS
                       for i in range(max(0, len(s) - 2)))
        if _count_starts(region_new) > _count_starts(ref_region):
            warnings.append('new start codon (ATG/GTG/TTG) near edit')

    # --- newly created in-frame stop codons ------------------------------
    if cfg.check_new_stop and frame_start is not None:
        f0 = frame_start - 1
        new_stops = _inframe_stops(new_seq, f0)
        ref_stops = _inframe_stops(ref_seq, f0)
        if new_stops > ref_stops:
            warnings.append('new in-frame stop codon')

    # --- newly created SD core -------------------------------------------
    if cfg.check_new_sd:
        a = max(0, lo - 5)
        region_new = new_seq[a:hi + 5]
        ref_region = ref_seq[a:a + len(region_new)]
        if (len(_SD_RE.findall(region_new))
                > len(_SD_RE.findall(ref_region))):
            warnings.append('new Shine-Dalgarno-like motif near edit')

    # --- homopolymer run -------------------------------------------------
    if cfg.check_homopolymer:
        run_new = _longest_run(new_seq[max(0, lo - 4):hi + 4])
        run_ref = _longest_run(ref_seq[max(0, lo - 4):hi + 4])
        if run_new > cfg.homopolymer_max and run_new > run_ref:
            warnings.append(f'homopolymer run of {run_new}')

    # --- GC content of edited window -------------------------------------
    if cfg.check_gc and len(win) >= 6:
        gc = _gc(win)
        if gc < cfg.gc_low or gc > cfg.gc_high:
            warnings.append(f'local GC {gc:.0f}%')

    return warnings


def _inframe_stops(seq: str, frame0: int) -> int:
    n = 0
    for i in range(frame0, len(seq) - 2, 3):
        if seq[i:i + 3] in _STOP_CODONS:
            n += 1
    return n
