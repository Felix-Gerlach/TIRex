# TIRex — Translation Initiation Rate Explorer

Desktop GUI for predicting translation initiation rates (TIR) of bacterial mRNA sequences. Built on the [OSTIR](https://github.com/barricklab/ostir) prediction engine, TIRex adds interactive visualization, ORF analysis, protein property calculations, and a virtual SDS-PAGE gel simulator.

## Features

- **TIR prediction** — runs OSTIR on any DNA/RNA sequence to find all start codons and predict their translation initiation rates
- **Interactive sequence map** — linear and circular views showing ORFs as color-coded arrows with start codon type, stop codon, and RBS position
- **Zoom-dependent detail** — zooming into the linear view reveals individual nucleotide letters along the backbone and amino acid letters inside each ORF
- **ORF results table** — sortable table with TIR, free energy components (dG total, dG rRNA:mRNA, dG mRNA, dG spacing, dG standby, dG start codon), protein length, MW, pI, and GRAVY
- **TIR range filter** — slider to show only ORFs within a specific TIR range; affects both the visualization and the table
- **ORF group coloring** — highlights start codons that share the same stop codon (alternative starts of the same protein) in distinct colors
- **Protein analysis** — molecular weight (kDa), isoelectric point, GRAVY hydrophobicity score, and amino acid composition via Biopython
- **Copy amino acid sequences** — one-click clipboard copy of any ORF's translated amino acid sequence
- **SDS-PAGE gel simulator** — virtual Coomassie-stained gel with configurable lanes, marker ladders, acrylamide percentage, E. coli contaminant bands, and export to image
- **FASTA import** — load sequences from FASTA files or paste directly
- **Circular sequence support** — for plasmid and circular genome analysis

## Requirements

- Python 3.10+
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download) — must be installed and the executables (`RNAfold`, `RNAsubopt`, `RNAeval`) must be in your system PATH

### Python dependencies

Install with:

```bash
pip install -r requirements.txt
```

| Package | Purpose |
|---------|---------|
| PyQt6 | GUI framework |
| matplotlib | Sequence visualization and gel rendering |
| biopython | Protein property analysis (MW, pI, GRAVY) |
| pandas | Data handling and CSV export |
| ostir | Translation initiation rate prediction engine |

## Usage

```bash
python main.py
```

Or build a standalone Windows executable:

```bash
build.bat
```

The `.exe` will be in `dist\TIRex\TIRex.exe`. ViennaRNA must still be in PATH on the target machine.

---

## User Manual

### 1. Entering a sequence

The left panel contains the input area:

- **Paste** a raw DNA or RNA sequence directly into the text box
- **Load FASTA** — click the "Load FASTA..." button to import a `.fasta`, `.fa`, `.fna`, or `.txt` file. The first sequence in the file is used and the filename is auto-filled as the sequence name
- **Sequence name** — optional label shown in the results (defaults to "sequence")

### 2. Configuring OSTIR options

Below the sequence input, adjust these parameters before running:

| Option | Description | Default |
|--------|-------------|---------|
| **Start pos** | Most 5' position to scan for start codons (1-indexed). 0 = scan entire sequence | 0 (Auto) |
| **End pos** | Most 3' position to scan. 0 = scan to end | 0 (Auto) |
| **Anti-SD** | 9 bp anti-Shine-Dalgarno sequence (3' end of 16S rRNA) | `ACCTCCTTA` (E. coli) |
| **Circular** | Treat the sequence as circular (for plasmids) | Off |
| **Threads** | Number of parallel threads for OSTIR | 1 |
| **Decimals** | Decimal places in TIR output | 4 |
| **Constraints** | ViennaRNA folding constraints string (advanced) | Empty |

For most use cases, the defaults work well. Change the **Anti-SD** sequence only when working with non-E. coli organisms.

### 3. Running the analysis

Click **"Run OSTIR"**. The button changes to "Running..." and the status bar shows progress. OSTIR scans for all start codons (ATG, GTG, TTG) and predicts each one's translation initiation rate. For each start codon, TIRex then:

1. Finds the downstream ORF (to the next in-frame stop codon)
2. Translates the ORF to an amino acid sequence
3. Calculates protein properties (MW, pI, GRAVY) via Biopython

### 4. Reading the visualization

The top panel shows the sequence map. Toggle between views using the **Linear** and **Circular** buttons.

#### Linear view

- **Backbone** — grey horizontal line representing the full sequence
- **ORF arrows** — right-pointing arrows, one per start codon. Length = full ORF from start to stop codon
- **Color gradient** — arrows are colored by TIR on a log10 scale (viridis colormap). Higher TIR = brighter/yellow
- **Start codon caps** — colored band at the left edge of each arrow: red = ATG, orange = GTG, amber = TTG
- **Stop codon caps** — dark band at the right tip of each arrow
- **RBS markers** — yellow pills above the backbone showing the Shine-Dalgarno position. A dotted line connects the RBS to its start codon
- **Labels** — each ORF is numbered (#1, #2, ...) ranked by position. Large ORFs show the label inside; small ORFs show it as a callout above
- **Zoom** — scroll the mouse wheel to zoom in. At high zoom, individual nucleotide letters (color-coded: A=green, T=red, G=yellow, C=blue) appear along the backbone, and amino acid letters appear inside the ORF arrows
- **Pan** — use the matplotlib toolbar at the bottom to pan or reset the view

#### Circular view

Same information rendered on a circular map. ORFs are arcs inside the backbone ring. Position ticks and labels ring the outside. The center shows sequence length and number of ORFs displayed.

#### Toolbar options

- **TIR range** — set minimum and maximum TIR values. Only ORFs within this range are drawn and shown in the table. Click the reset button to restore the full range
- **Color by ORF group** — when checked, start codons sharing the same stop codon (i.e. alternative start sites of the same protein) are given the same color. Useful for identifying which start codons lead to the same protein

### 5. Using the ORF results table

The bottom panel lists every detected ORF with these columns:

| Column | Meaning |
|--------|---------|
| **Show** | Checkbox — toggle visibility of this ORF on the map |
| **Copy AA** | Click to copy the amino acid sequence to clipboard |
| **#** | Display rank (by start position, restarts from 1 when filtering) |
| **Position** | Start codon position (1-indexed) |
| **Codon** | Start codon type (ATG/GTG/TTG) |
| **TIR** | Predicted translation initiation rate |
| **RBS dist** | Distance in bp between the Shine-Dalgarno sequence and start codon |
| **dG total** | Total free energy of translation initiation |
| **dG rRNA:mRNA** | Free energy of rRNA-mRNA hybridization |
| **dG mRNA** | Free energy of mRNA folding near the start codon |
| **dG spacing** | Free energy penalty for non-optimal RBS-start spacing |
| **dG standby** | Free energy of the standby site |
| **dG start** | Free energy contribution of the start codon identity |
| **Length (aa)** | Protein length in amino acids |
| **MW (kDa)** | Molecular weight in kilodaltons |
| **pI** | Isoelectric point |
| **GRAVY** | Grand average of hydropathy (positive = hydrophobic, negative = hydrophilic) |
| **Stop?** | Whether a stop codon was found in-frame |

**Sorting** — click any column header to sort. Click again to reverse.

**Show/Hide all** — buttons to toggle all ORFs at once.

**Right-click context menu:**
- *Set as target* — marks an ORF as the target protein for the SDS-PAGE simulator
- *Clear target* — removes the target selection
- *Simulate SDS-PAGE* — opens the gel simulator

### 6. SDS-PAGE gel simulator

Access via the **"SDS-PAGE..."** button in the table header or the right-click menu.

The simulator renders a virtual Coomassie-stained gel based on the predicted proteins:

- **Marker lane** — choose from preset ladders (PageRuler Prestained Plus, Precision Plus Dual Color, etc.)
- **Acrylamide %** — select gel percentage (affects band migration)
- **Lane types:**
  - *Marker* — molecular weight standard
  - *Target only* — just the target protein
  - *Contaminants only* — common E. coli His-tag co-purifying contaminants
  - *Target group (incl. alt starts)* — the target protein plus all alternative start-site variants sharing the same stop codon
  - *Target group + contaminants* — combined lane
  - *Custom* — add arbitrary bands
- **Add/remove lanes** — build your gel layout
- **Contaminant table** — view and configure the E. coli contaminant band list
- **Export** — save the gel image to file

### 7. Exporting results

- **File → Export CSV** — save the full ORF results table as a CSV file
- **File → Save Figure** — save the current visualization (linear or circular) as a PNG image
- **Copy AA** — copy individual amino acid sequences from the table
- **SDS-PAGE export** — save the simulated gel image from within the gel dialog

---

## Built with

| Tool | Purpose |
|------|---------|
| [OSTIR](https://github.com/barricklab/ostir) | Translation initiation rate prediction engine |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | RNA secondary structure prediction |
| [PyQt6](https://www.riverbankcomputing.com/software/pyqt/) | GUI framework |
| [Matplotlib](https://matplotlib.org/) | Sequence visualization and gel rendering |
| [Biopython](https://biopython.org/) | Protein property analysis |
| [pandas](https://pandas.pydata.org/) | Data handling and CSV export |
| [Claude Code](https://claude.ai/claude-code) | AI-assisted development |

## License

GPL-3.0
