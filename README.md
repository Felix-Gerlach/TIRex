# TIRex — Translation Initiation Rate Explorer

Desktop GUI for predicting and **engineering** translation initiation rates (TIR) of bacterial mRNAs. Built on the [OSTIR](https://github.com/barricklab/ostir) prediction engine, TIRex adds interactive visualization, ORF analysis, protein property calculations, a sequence‑editing optimizer (**TIR‑Tuner**), a codon optimizer, a virtual SDS‑PAGE gel, batch scoring, and session/report export.

---

## Table of contents

- [Feature overview](#feature-overview)
- [Requirements & install](#requirements--install)
- [Running / building](#running--building)
- [User manual](#user-manual)
  - [1. Loading a sequence](#1-loading-a-sequence)
  - [2. Anti‑SD configuration & dual scoring](#2-anti-sd-configuration--dual-scoring)
  - [3. OSTIR options & running](#3-ostir-options--running)
  - [4. Visualization](#4-visualization)
  - [5. ORF results table](#5-orf-results-table)
  - [6. Target protein & in‑frame filtering](#6-target-protein--in-frame-filtering)
  - [7. TIR‑Tuner — sequence‑editing optimizer](#7-tir-tuner--sequence-editing-optimizer)
  - [8. Codon optimizer](#8-codon-optimizer)
  - [9. SDS‑PAGE gel simulator](#9-sds-page-gel-simulator)
  - [10. Batch OSTIR](#10-batch-ostir)
  - [11. Sessions & exports](#11-sessions--exports)
- [Menu reference](#menu-reference)
- [Project layout](#project-layout)
- [Notes, calibration & limitations](#notes-calibration--limitations)
- [Built with](#built-with)
- [License](#license)

---

## Feature overview

**Prediction & analysis**
- Runs OSTIR to find every start codon (ATG/GTG/TTG) and predict its TIR plus the full free‑energy breakdown.
- **Dual anti‑SD scoring** — score every ORF against two ribosomes (e.g. native *E. coli* + a custom host) and compare side‑by‑side.
- Per‑ORF protein properties: length, MW, pI, GRAVY, amino‑acid composition (Biopython).
- Interactive **linear / circular** sequence map with zoom‑dependent nucleotide & amino‑acid lettering.
- Sortable results table (correct **numeric sorting**), TIR‑range filter, ORF‑group coloring.

**Engineering (new)**
- **TIR‑Tuner** — proposes nucleotide edits to **increase** initiation at a target start, **decrease** a downstream internal start without changing the protein, or build an **RBS library / TIR ramp**. Includes greedy and **beam search**, **parallel scoring**, constraint **warnings**, iterative in‑place editing, old‑vs‑new diff, and **HTML report** export with **Q5 mutagenesis primers**.
- **Codon optimizer** — synonymous CDS rewrite to raise CAI for a host without changing the protein.
- **SDS‑PAGE gel simulator** — virtual Coomassie gel of the predicted protein mixture.
- **Batch OSTIR** — score every record across one or more multi‑FASTA files into one table.

**I/O**
- Import **FASTA**, **GenBank** (`.gb/.gbk`, with CDS features), and **SnapGene** (`.dna`).
- Save/restore complete **sessions** (`.tirex`).
- Export results CSV, protein FASTA, visualization PNG/SVG, gel image, and the tuning report.

---

## Requirements & install

- **Python 3.10+**
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download)** — the executables `RNAfold`, `RNAsubopt`, `RNAeval` must be on your system **PATH**. TIRex checks for them on startup and refuses to run if they're missing.

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

| Package | Purpose |
|---------|---------|
| PyQt6 | GUI framework |
| matplotlib | Sequence visualization, ΔTIR heatmap, gel rendering |
| biopython | Protein properties (MW/pI/GRAVY) **and** GenBank/SnapGene import |
| pandas | CSV export |
| ostir | Translation‑initiation‑rate prediction engine |

---

## Running / building

```bash
python main.py
```

Build a standalone Windows executable:

```bash
build.bat        # output in dist\TIRex\TIRex.exe
```

ViennaRNA must still be on PATH on the target machine. (The launcher calls `multiprocessing.freeze_support()` so the optimizer's parallel scoring works in frozen builds.)

---

## User manual

### 1. Loading a sequence

In the left **Input** panel you can either paste a sequence or load a file with **“Load sequence…”** (`Ctrl+O`). Supported formats:

| Format | Extensions | Notes |
|--------|-----------|-------|
| FASTA | `.fasta .fa .fna .txt` | First record used; header becomes the name |
| GenBank | `.gb .gbk .genbank` | Sequence + forward‑strand **CDS features** are detected and reported |
| SnapGene | `.dna` | Read via Biopython's `snapgene` parser |
| VectorBee | `.vbee` | **Not yet supported** — shows a clear message; export to GenBank/FASTA for now |

**Input sanitization (robust against the old crash):** pasted/loaded sequence is cleaned to `A/C/G/T` (U→T); anything else (digits, `N`, dashes, whitespace) is removed and you're notified. This prevents OSTIR's internal validator from raising the cryptic `print() got an unexpected keyword 'style'` error.

### 2. Anti‑SD configuration & dual scoring

The **Anti‑SD (ribosome)** group controls which ribosome(s) the TIR is scored against. The anti‑SD is the 9‑nt 3′ tail of the 16S rRNA that base‑pairs with the mRNA Shine‑Dalgarno.

- **Primary (always scored)** — dropdown of saved presets + an editable field. Default `E. coli (native) = ACCTCCTTA`.
- **“Also score with a secondary anti‑SD”** — when ticked, every ORF is **also** scored against this second anti‑SD, producing a side‑by‑side comparison column. Default preset `pET T7 (derived) = CTCCTTCTT`.
- **Derive from RBS…** — paste an mRNA RBS/5′ leader; TIRex finds the SD core and stores its reverse complement as a new named anti‑SD preset.
- **Save…** — persist the current secondary anti‑SD as a named preset.

Presets are stored in `asd_presets.json` next to the program and survive between sessions. Anti‑SDs must be exactly **9 A/C/G/T** bases — TIRex validates this before running and warns clearly otherwise.

> Result: the table shows **`TIR (1° aSD)`** and **`TIR (2° aSD)`** columns. Both numbers come straight from OSTIR (two passes), not an approximation.

### 3. OSTIR options & running

| Option | Description | Default |
|--------|-------------|---------|
| Start pos / End pos | Restrict the start‑codon scan window (1‑indexed; 0 = auto) | Auto |
| Circular | Treat the sequence as circular (plasmids) | Off |
| Threads | OSTIR worker threads | 1 |
| Decimals | TIR decimal places | 4 |
| Constraints | ViennaRNA folding‑constraint string (advanced) | Empty |

Click **Run OSTIR**. For each start codon TIRex finds the downstream ORF (to the next in‑frame stop), translates it, and computes protein properties.

### 4. Visualization

Top panel; toggle **Linear / Circular**.

- ORFs drawn as arrows/arcs, colored by TIR (log scale); start caps colored by codon (ATG red, GTG orange, TTG amber); stop caps dark; RBS pills with a connector to the start codon; every ORF numbered.
- **Mouse‑wheel zoom**; at high zoom, per‑base letters (A green / T red / G yellow / C blue) and amino‑acid letters appear.
- **TIR range filter** (min/max + reset) — hides out‑of‑range ORFs in both the map and the table and renumbers.
- **Color by ORF group** — start codons sharing a stop codon get the same color.
- Matplotlib toolbar for pan/zoom; export via the Export menu.

### 5. ORF results table

| Column | Meaning |
|--------|---------|
| Show | Toggle ORF visibility on the map |
| Copy AA | Copy the translated amino‑acid sequence |
| # | Display rank by start position (restarts at 1 when filtered) |
| Position / Codon | Start position & codon type |
| **TIR (1° aSD)** | TIR against the primary anti‑SD |
| **TIR (2° aSD)** | TIR against the secondary anti‑SD (— if disabled) |
| RBS dist | SD↔start spacing (bp) |
| dG total / rRNA:mRNA / mRNA / spacing / standby / start | OSTIR free‑energy components |
| Length (aa) / MW (kDa) / pI / GRAVY | Protein properties |
| Stop? | In‑frame stop found |

- **Sorting** is now **numeric** on numeric columns (e.g. `67.9` sorts before `4385.2`); blanks sink to the bottom. Click a header to sort, again to reverse.
- **Show all / Hide all** buttons; **Clear target** and **SDS‑PAGE…** buttons.
- **Right‑click a row:** *Set as target*, *Tune translation rate…*, *Clear target*, *Simulate SDS‑PAGE…*.

### 6. Target protein & in‑frame filtering

Right‑click an ORF → **Set as target**. TIRex then hides every fragment **not in the target's reading frame** (`(start − target_start) mod 3 ≠ 0`) from both the table and the map, and renumbers the rest. This isolates the alternative start sites of one protein. **Clear target** (button or Tools menu) restores all ORFs. The target also seeds the SDS‑PAGE “target group” lanes.

### 7. TIR‑Tuner — sequence‑editing optimizer

Open via **Tools → Tune translation rate… (`Ctrl+T`)** or a row's right‑click menu (seeds that ORF as the target). The Tuner proposes concrete nucleotide edits, scores each with OSTIR, and lets you apply them iteratively.

**Scoring anti‑SD** — a dropdown at the top selects which ribosome the optimizer scores against (your Primary and Secondary anti‑SDs). **Switching it re‑scans** the current sequence so all numbers reflect that ribosome.

#### Modes

**(a) Increase initiation at a target start**
Edits the upstream window between the SD and the start codon. Options:
- *Upstream window (nt)* and *Include dinucleotide substitutions*.
- **Search depth:**
  - *Single+dinuc* — enumerate all single substitutions/insertions/deletions + dinucleotide substitutions (depth‑1), ranked by TIR.
  - *Greedy* — repeatedly apply the single best edit until a *target fold* or *max rounds* is reached.
  - *Beam* — keep the best *beam‑width* partial solutions each round and expand each; finds stronger multi‑edit combinations than greedy.

**(b) Decrease a downstream start (protein‑preserving)**
Suppresses an internal start codon **without changing the main protein** — every edit is synonymous in the main reading frame. Options:
- *Main start (keep)* and *Downstream start (suppress)*.
- **“Only downstream starts in frame with main”** — restricts the downstream list to internal in‑frame starts (position a multiple of 3 from the main start).
- *Strategy:* `codon swap` (mutate the codon(s) overlapping the internal start so it's no longer ATG/GTG/TTG), `RBS synonymous` (synonymously weaken the internal SD/spacer), or `both`.
- *RBS window (nt)*. Each candidate reports the downstream ΔTIR **and** the main‑start ΔTIR (which should stay ~0).

**(c) Build RBS library / TIR ramp**
Produces a panel of *N* variants whose predicted TIRs are **log‑spaced** across the achievable range (weak → strong), including wild‑type — an expression dilution series for experiments.

#### Constraint warnings (soft‑flag)

Every candidate is checked and any issues appear in a **Warnings** column (amber, with tooltips) — nothing is rejected, you decide:

- newly‑created **restriction site** (default set: EcoRI, BamHI, HindIII, XhoI, NdeI, NotI, XbaI, SpeI, PstI, SalI, KpnI, SacI, NcoI, BglII, NheI),
- newly‑created **start codon** (ATG/GTG/TTG) near the edit,
- newly‑created **in‑frame stop codon**,
- newly‑created **Shine‑Dalgarno‑like motif**,
- **homopolymer run** > 4,
- local **GC** outside 25–75 %.

#### Results & iterative editing

The results table mirrors the main table plus **Apply**, **Warnings**, **Notes**, and **Copy seq** columns, and (decrease mode) **Main ΔTIR**.

- **Tick “Apply”** on a row to adopt that edit. The Tuner **stays open**, updates to the new sequence, **re‑scans** start codons, logs the edit under **Applied edits (this session)**, and clears the results so you can predict the next change. This enables rapid stacking of edits.
- **↩ Undo last** reverts the most recent applied edit.
- **Compare…** opens a position‑aware **old‑vs‑new diff** (amber = substitution, green = insertion, red = deletion).
- **✓ Apply to TIRex & close** sends the final tuned sequence back to the main window, which loads it (auto‑versioning the name `…_v1, _v2`) and re‑runs OSTIR.

Under the hood, candidate scoring runs in **parallel** (process pool, workers = CPU−1) for large batches, with a sequential fallback, and an OSTIR result cache avoids re‑folding duplicates.

#### ΔTIR heatmap & region labels

The **Graph** tab shows a saturation‑mutagenesis heatmap: position × base, colored by **ΔTIR** (green = boosts initiation, red = lowers it), wild‑type cells outlined. You can add **labeled regions** (e.g. *SD*, *spacer*, *start codon*) as colored bars above the heatmap (From/To/label/color), remove, or clear them.

#### Export HTML report

**Export report…** writes a self‑contained HTML file containing: the wild‑type → final TIR and fold change, the run metadata, the list of applied edits, the colour‑coded sequence diff, **Q5 / NEBaseChanger‑style mutagenesis primers** for each applied edit, and the embedded ΔTIR heatmap. Opens in any browser.

> Primer Tm is an estimate — verify with the NEB Tm calculator before ordering.

### 8. Codon optimizer

**Tools → Codon optimizer (CAI)…**. Pick an ORF (uses its CDS) or paste a CDS, choose a **host** (ships *E. coli K‑12*; extensible via `CODON_USAGE`), and set **Preserve first N codons** (they overlap the initiation region and affect TIR). **Optimize** rewrites the CDS to the highest‑usage synonymous codon at each position **without changing the protein**, and reports **CAI before→after**, **GC before→after**, codons changed, and confirms the protein is preserved. Copy the optimized CDS to the clipboard.

### 9. SDS‑PAGE gel simulator

**Tools → Simulate SDS‑PAGE gel…**, the table's **SDS‑PAGE…** button, or a row's right‑click menu. Renders a virtual Coomassie gel of the predicted protein mixture:
- Selectable **marker ladder** and **acrylamide %** (affects migration).
- Configurable **lanes**: marker, target only, target group (incl. alt starts), contaminants only, combinations, custom.
- Editable list of common **E. coli His‑tag co‑purifying contaminants**, with a master‑intensity slider and uniform/TIR‑proportional band intensities.
- Respects the active **TIR filter** and only draws **selected** (visible) proteins. Export the gel image.

### 10. Batch OSTIR

**Tools → Batch OSTIR (multi‑FASTA)…**. Add one or more (multi‑)FASTA files, set the anti‑SD, and **Run batch**. Every record is scored in a background thread and collected into one sortable table (Record, Start pos, Codon, TIR, dG total, RBS dist). **Export CSV** for the combined results.

### 11. Sessions & exports

- **File → Save session… (`Ctrl+S`)** / **Open session…** — persist/restore the whole working state (sequence, anti‑SD params, full results, target, TIR filter) as a `.tirex` JSON file.
- **Export → Results as CSV (`Ctrl+Shift+S`)** — full ORF table (both TIR columns + composition).
- **Export → Protein sequences as FASTA** — visible ORFs' amino‑acid sequences.
- **Export → Visualization as PNG / SVG** — the current map.
- Tuner **Export report…** (HTML) and SDS‑PAGE image export as described above.

---

## Menu reference

| Menu | Items |
|------|-------|
| **File** | Load sequence (FASTA/GenBank/.dna)… `Ctrl+O` · Save session… `Ctrl+S` · Open session… · Quit `Ctrl+Q` |
| **Export** | Results as CSV… `Ctrl+Shift+S` · Protein sequences as FASTA… · Visualization as PNG… · Visualization as SVG… |
| **Tools** | Tune translation rate (TIR‑Tuner)… `Ctrl+T` · Codon optimizer (CAI)… · Batch OSTIR (multi‑FASTA)… · Simulate SDS‑PAGE gel… · Clear target protein |
| **Help** | About TIRex |

---

## Project layout

```
TIRex/
  main.py                       app entry (dependency check, theme, freeze_support)
  core/
    ostir_runner.py             QThread wrapper; dual‑aSD scoring; input guards
    orf_finder.py               ORF extent + translation
    protein_analysis.py         MW / pI / GRAVY (Biopython)
    gel_simulator.py            SDS‑PAGE physics, markers, contaminants
    asd_presets.py              anti‑SD preset store + derive‑from‑RBS
    codon_opt.py                CAI + host codon usage + synonymous optimizer
    primers.py                  Q5 / NEBaseChanger primer design
    seq_import.py               FASTA / GenBank / SnapGene loaders
    session.py                  .tirex save/load
    report.py                   self‑contained HTML report builder
    tuner/
      genetic_code.py           codon table, synonymous map, start/stop codons
      scoring.py                OstirScorer: caching + parallel batch scoring
      mutations.py              edit generators + synonymous enumeration
      candidate.py              Candidate dataclass (TIR, Δ, warnings, …)
      constraints.py            motif/synthesis soft‑flag checks
      increase_engine.py        enumerate / greedy / beam
      decrease_engine.py        codon swap / RBS synonymous
      library.py                RBS library / TIR ramp
  ui/
    theme.py                    central palette + QSS
    main_window.py              window, menus, wiring
    input_panel.py              sequence input, anti‑SD group, options
    visualization_widget.py     linear/circular map
    orf_table_widget.py         results table (numeric sort, frame filter)
    sds_page_widget.py          gel dialog
    tuner_dialog.py             TIR‑Tuner (+ diff dialog)
    tuner_worker.py             tuner background thread
    result_plot.py              ΔTIR heatmap + region labels
    codon_optimizer_dialog.py   codon optimizer UI
    batch_dialog.py             batch OSTIR UI
```

---

## Notes, calibration & limitations

- **OSTIR is calibrated for *E. coli*.** TIRs are in arbitrary (log‑scaled) units; treat them as relative. For other hosts, set the appropriate **anti‑SD** before scoring/tuning — RBS edits are only meaningful against the correct ribosome. Not applicable to eukaryotes or leaderless mRNAs.
- **Codon usage** ships *E. coli K‑12* only; add hosts in `core/codon_opt.py::CODON_USAGE`.
- **Met internal starts** can't be removed by codon swap (Met has a single codon) — use the RBS‑synonymous strategy instead. Decrease‑mode edits are bounded to the main CDS so they stay synonymous.
- **Primer Tm** is an estimate; verify with NEB's calculator. Primers assume whole‑plasmid PCR + KLD (NEBaseChanger style).
- **`.vbee` import is not implemented** (format undocumented here) — please share a sample file to add support; meanwhile export to GenBank/FASTA.
- Constraint checks are **advisory** (soft‑flagged), not hard filters.

---

## Built with

| Tool | Purpose |
|------|---------|
| [OSTIR](https://github.com/barricklab/ostir) | TIR prediction engine |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | RNA secondary‑structure prediction |
| [PyQt6](https://www.riverbankcomputing.com/software/pyqt/) | GUI framework |
| [Matplotlib](https://matplotlib.org/) | Visualization, heatmap, gel |
| [Biopython](https://biopython.org/) | Protein analysis + sequence file parsing |
| [pandas](https://pandas.pydata.org/) | CSV export |
| [Claude Code](https://claude.ai/claude-code) | AI‑assisted development |

## License

GPL‑3.0
