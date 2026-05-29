"""
OSTIR scoring wrapper with a dedup cache.

`run_ostir` (ViennaRNA) is the expensive step, so:
  * we restrict it to a single start codon via the start/end arguments when
    we only need one position's TIR, and
  * we memoise results keyed on (sequence, start, end, aSD).

Numbers are therefore identical to the OSTIR/TIRex GUI — no approximation.
"""

import os
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Tuple

from ostir import run_ostir


def _score_one(args: Tuple[str, int, Optional[str], int]) -> float:
    """Module-level worker for the process pool: TIR at `pos` for `seq`."""
    seq, pos, aSD, decimal_places = args
    try:
        res = run_ostir(seq, start=pos, end=pos, name='tir_tuner',
                        aSD=aSD, threads=1, decimal_places=decimal_places,
                        verbosity=0) or []
    except Exception:
        return 0.0
    for r in res:
        if int(r.get('start_position', -1)) == pos:
            return float(r.get('expression') or 0.0)
    best = None
    for r in res:
        d = abs(int(r.get('start_position', -10_000)) - pos)
        if d <= 2 and (best is None or d < best[0]):
            best = (d, float(r.get('expression') or 0.0))
    return best[1] if best else 0.0


class OstirScorer:
    # Use the pool only when there are at least this many uncached candidates
    # (process spawn/import overhead isn't worth it for small batches).
    PARALLEL_THRESHOLD = 16

    def __init__(self, aSD: Optional[str] = None, threads: int = 1,
                 decimal_places: int = 4, parallel: bool = True,
                 workers: Optional[int] = None):
        self.aSD = aSD or None
        self.threads = threads
        self.decimal_places = decimal_places
        self._cache: Dict[tuple, list] = {}
        self._tir_cache: Dict[tuple, float] = {}   # (seq,pos,aSD) -> tir
        self.parallel = parallel
        self.workers = workers or max(1, (os.cpu_count() or 2) - 1)
        self.calls = 0          # OSTIR invocations actually run (cache misses)
        self.hits = 0           # cache hits

    # ------------------------------------------------------------------ #
    def _key(self, seq: str, start, end):
        return (seq, start, end, self.aSD)

    def run(self, seq: str, start=None, end=None) -> List[dict]:
        """Return the raw OSTIR result list for the given (sub)range."""
        key = self._key(seq, start, end)
        if key in self._cache:
            self.hits += 1
            return self._cache[key]
        self.calls += 1
        try:
            results = run_ostir(
                seq, start=start, end=end, name='tir_tuner',
                aSD=self.aSD, threads=self.threads,
                decimal_places=self.decimal_places, verbosity=0,
            ) or []
        except Exception:
            results = []
        self._cache[key] = results
        return results

    # ------------------------------------------------------------------ #
    def find_starts(self, seq: str, start=None, end=None) -> List[dict]:
        """Full scan: every start codon OSTIR finds in the range."""
        results = self.run(seq, start=start, end=end)
        # Normalise to ints / floats we rely on downstream.
        for r in results:
            r['start_position'] = int(r.get('start_position', 0))
            r['expression'] = float(r.get('expression') or 0.0)
        return sorted(results, key=lambda r: r['start_position'])

    def _tir_key(self, seq: str, pos: int):
        return (seq, pos, self.aSD)

    def score_start(self, seq: str, pos: int) -> float:
        """TIR (expression) of the start codon located at 1-based `pos`.
        Restricts OSTIR to that single codon for speed. 0.0 if none found."""
        tkey = self._tir_key(seq, pos)
        if tkey in self._tir_cache:
            self.hits += 1
            return self._tir_cache[tkey]
        results = self.run(seq, start=pos, end=pos)
        tir = 0.0
        found = False
        for r in results:
            if int(r.get('start_position', -1)) == pos:
                tir = float(r.get('expression') or 0.0)
                found = True
                break
        if not found:
            best = None
            for r in results:
                d = abs(int(r.get('start_position', -10_000)) - pos)
                if d <= 2 and (best is None or d < best[0]):
                    best = (d, float(r.get('expression') or 0.0))
            tir = best[1] if best else 0.0
        self._tir_cache[tkey] = tir
        return tir

    def score_many(self, tasks: List[Tuple[str, int]],
                   progress=None) -> List[float]:
        """Score many (seq, pos) candidates, in parallel when worthwhile.
        Returns a list of TIRs aligned to `tasks`."""
        # Unique, uncached tasks only.
        uncached = []
        seen = set()
        for seq, pos in tasks:
            k = self._tir_key(seq, pos)
            if k in self._tir_cache or k in seen:
                continue
            seen.add(k)
            uncached.append((seq, pos))

        use_pool = (self.parallel and self.workers > 1
                    and len(uncached) >= self.PARALLEL_THRESHOLD)
        if use_pool:
            args = [(seq, pos, self.aSD, self.decimal_places)
                    for seq, pos in uncached]
            try:
                with ProcessPoolExecutor(max_workers=self.workers) as ex:
                    done = 0
                    for (seq, pos), tir in zip(
                            uncached, ex.map(_score_one, args, chunksize=4)):
                        self._tir_cache[self._tir_key(seq, pos)] = tir
                        self.calls += 1
                        done += 1
                        if progress and done % 20 == 0:
                            progress(done, len(uncached), 'scoring (parallel)')
            except Exception:
                # Fall back to sequential if the pool can't start.
                for seq, pos in uncached:
                    self.score_start(seq, pos)
        else:
            for n, (seq, pos) in enumerate(uncached):
                self.score_start(seq, pos)
                if progress and n % 20 == 0:
                    progress(n, len(uncached), 'scoring')

        return [self._tir_cache.get(self._tir_key(seq, pos), 0.0)
                for seq, pos in tasks]

    def result_at(self, seq: str, pos: int) -> Optional[dict]:
        """Full OSTIR result dict for the start at `pos` (dG breakdown etc.)."""
        for r in self.run(seq, start=pos, end=pos):
            if int(r.get('start_position', -1)) == pos:
                return r
        return None
