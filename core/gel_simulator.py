"""
SDS-PAGE gel simulation model.

Contains:
  * MW ladder definitions (MARKERS)
  * E. coli His-tag co-purifying contaminants (CONTAMINANTS)
  * Acrylamide-percent-dependent migration function (compute_rf)
  * Lane/band dataclasses consumed by the renderer in ui/sds_page_widget.py

Everything is *pure Python* — no Qt / matplotlib imports — so this module
is easy to unit-test.
"""

import math
from dataclasses import dataclass, field
from typing import List

# ----------------------------------------------------------------------
# MW markers (ladders)
# ----------------------------------------------------------------------
# Each marker: ordered list of (mw_kda, band_color_hex). Colors preserve
# prestained-ladder colors so the marker lane looks realistic.

MARKERS = {
    'PageRuler Prestained Plus': [
        (250, '#1f77b4'),
        (130, '#1f77b4'),
        (100, '#1f77b4'),
        (70,  '#e74c3c'),   # orange/red highlight band
        (55,  '#1f77b4'),
        (35,  '#1f77b4'),
        (25,  '#2ecc71'),   # green highlight band
        (15,  '#1f77b4'),
        (10,  '#1f77b4'),
    ],
    'PageRuler Unstained Broad Range': [
        (mw, '#333') for mw in (
            250, 200, 150, 120, 100, 85, 70, 60,
            50, 40, 30, 25, 20, 15, 10, 5,
        )
    ],
    'Precision Plus All Blue': [
        (mw, '#1976d2') for mw in (
            250, 150, 100, 75, 50, 37, 25, 20, 15, 10,
        )
    ],
    'Generic 10-band': [
        (mw, '#222') for mw in (
            250, 180, 130, 100, 70, 50, 35, 25, 15, 10,
        )
    ],
    'Tricolor Prestained': [
        (180, '#1976d2'),
        (130, '#1976d2'),
        (95,  '#1976d2'),
        (72,  '#2e7d32'),   # green highlight
        (55,  '#1976d2'),
        (43,  '#1976d2'),
        (34,  '#c62828'),   # red highlight
        (26,  '#1976d2'),
        (17,  '#1976d2'),
        (11,  '#1976d2'),
    ],
}


# ----------------------------------------------------------------------
# E. coli His-tag co-purifying contaminants
# ----------------------------------------------------------------------

@dataclass
class Contaminant:
    name: str
    mw_kda: float
    default_intensity: float     # 0..1
    note: str = ''
    enabled: bool = True


def default_contaminants() -> List[Contaminant]:
    """Return a fresh list of the canonical His-tag co-purifying contaminants."""
    return [
        Contaminant('SlyD',  21.0, 0.75, 'PPIase, intrinsic His-rich'),
        Contaminant('Crr',   18.0, 0.50, 'PTS EIIA-Glc'),
        Contaminant('Can',   25.0, 0.35, 'Carbonic anhydrase'),
        Contaminant('SodB',  21.0, 0.30, 'Superoxide dismutase'),
        Contaminant('GlmS',  67.0, 0.45, 'Glucosamine-6-P synthase'),
        Contaminant('ArnA',  74.0, 0.55, 'UDP-arabinose dehydrogenase'),
        Contaminant('DnaK',  70.0, 0.55, 'Chaperone'),
        Contaminant('GroEL', 60.0, 0.65, 'Chaperonin'),
        Contaminant('GlgB',  84.0, 0.30, 'Glycogen branching enzyme'),
    ]


# ----------------------------------------------------------------------
# Acrylamide-dependent migration
# ----------------------------------------------------------------------

ACRYLAMIDE_OPTIONS = ['8%', '10%', '12%', '15%', '4–20% gradient']

# Effective resolving window (kDa_lo, kDa_hi) for each gel.
_ACRYL_RANGE = {
    '8%':  (25, 200),
    '10%': (15, 120),
    '12%': (10, 70),
    '15%': (5,  50),
    '4–20% gradient': (3, 250),
}


def compute_rf(mw_kda: float, acrylamide: str) -> float:
    """
    Return relative mobility Rf for a protein of MW ``mw_kda`` on a gel of
    ``acrylamide`` percent. 0 = top of resolving gel (larger proteins), 1 =
    dye front (smallest). A log-linear fit inside the gel's resolving
    window, saturating gently outside it.
    """
    lo, hi = _ACRYL_RANGE.get(acrylamide, _ACRYL_RANGE['12%'])
    log_lo, log_hi = math.log10(lo), math.log10(hi)
    log_mw = math.log10(max(mw_kda, 0.5))
    rf_top, rf_bot = 0.06, 0.94

    if log_mw >= log_hi:
        excess = log_mw - log_hi
        return max(0.01, rf_top - 0.05 * excess)
    if log_mw <= log_lo:
        deficit = log_lo - log_mw
        return min(0.99, rf_bot + 0.05 * deficit)

    return rf_top + (rf_bot - rf_top) * (log_hi - log_mw) / (log_hi - log_lo)


# ----------------------------------------------------------------------
# Lane / band model
# ----------------------------------------------------------------------

@dataclass
class GelBand:
    mw_kda: float
    intensity: float                 # 0..1 (drives alpha)
    label: str = ''
    color: str = '#111'


@dataclass
class GelLane:
    name: str
    kind: str                        # see LANE_KINDS
    bands: List[GelBand] = field(default_factory=list)


LANE_KINDS = [
    'Marker',
    'Target only',
    'Target group (incl. alt starts)',
    'Contaminants only',
    'Target + contaminants',
    'Target group + contaminants',
    'All visible ORFs',
    'All visible ORFs + contaminants',
]
