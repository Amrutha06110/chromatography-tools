# Chromatography Peak Analysis App

An interactive **Python / Streamlit** application for visualizing, detecting, and integrating peaks in chromatographic data.
Upload one or multiple CSV files containing `(time, intensity)` pairs and the app performs smoothing, peak detection, trapezoidal integration, and percentage area calculation.

The core library uses an **object-oriented class hierarchy** so that technique-specific defaults (axis labels, smoothing parameters, detection thresholds) are built-in for HPLC, GC, SEC/GPC, and Ion Chromatography — while all shared signal-processing logic lives in a single abstract base class.

> **📖 New to the app?** See the [Step-by-Step User Guide](docs/USER_GUIDE.md).

---

## Features

- **Multiple chromatography types** via OOP hierarchy
  - `HPLCChromatogram` — UV/Vis absorbance detection
  - `GCChromatogram` — FID / TCD detection (sharper peaks)
  - `SECChromatogram` — Size-exclusion / GPC with optional MW calibration
  - `IonChromatogram` — Conductivity detection (broader peaks)

- **Peak detection & integration**
  - Height-only threshold (% of max signal)
  - Savitzky–Golay smoothing for noise reduction
  - Trapezoidal peak integration
  - Automated relative area (%) calculation

- **Visualization** (Plotly)
  - Overlay or separate-panel modes
  - Automatic viridis-derived color palette
  - Peak labels with optional % areas

- **Export**
  - Summary CSV with peak times and relative areas
  - Plot export as **PNG**, **SVG**, or **PDF**

---

## Installation

### Clone the repository
```bash
git clone https://github.com/andrewvgrassetti/chromatography-tools.git
cd chromatography-tools
```

### Install dependencies

```bash
pip install -r requirements.txt
```

Or install as a package (editable):

```bash
pip install -e ".[dev]"
```

### Running the App

```bash
streamlit run app.py
```

A browser window will open automatically.
If not, copy the printed URL (e.g., http://localhost:8501) into your browser.

## Input Format

Each chromatogram file must be a 2-column CSV:

| time | intensity |
|------|-----------|
| 0.00 | 12.3     |
| 0.01 | 12.9     |
| ...  | ...      |

- Time in minutes (or any consistent unit)
- Intensity in arbitrary units
- Multiple files may be uploaded simultaneously.

## Output

### Summary CSV

One row per file including:

- technique type
- total number of peaks
- per-peak times
- per-peak relative area (%)

### Plot export

Available formats:

- PNG
- SVG
- PDF

with adjustable size.

## Testing

Run automated tests:

```bash
pytest
```

Expected output includes checks for:

- OOP hierarchy and technique-specific defaults
- peak count accuracy
- integration correctness
- stability with noise
- SEC molecular-weight calibration
- method chaining

## Project Structure

```
chromatography-tools/
├── app.py                        # Streamlit web application
├── chromatography/
│   ├── __init__.py               # Package exports
│   └── core.py                   # OOP class hierarchy
│       ├── BaseChromatogram      #   Abstract base (smoothing, peaks, integration)
│       ├── HPLCChromatogram      #   HPLC / UV-Vis
│       ├── GCChromatogram        #   Gas chromatography
│       ├── SECChromatogram       #   Size-exclusion / GPC (+ MW calibration)
│       └── IonChromatogram       #   Ion chromatography
├── tests/
│   ├── test_chromatogram.py      # Comprehensive pytest suite
│   └── test_edit_bounds_chart.py # Interactive chart config tests
├── docs/
│   └── USER_GUIDE.md             # Step-by-step usage guide
├── pyproject.toml                # Python project metadata
├── requirements.txt              # Dependencies
└── README.md
```

## Class Hierarchy

```
BaseChromatogram (ABC)
├── HPLCChromatogram    — Chromatogram alias
├── GCChromatogram      — narrower smoothing window, lower height threshold
├── SECChromatogram     — elution-volume axis, MW calibration support
└── IonChromatogram     — wider smoothing window, conductivity axis
```

All subclasses inherit the full signal-processing pipeline (smooth → find_peaks → integrate_peaks → summarize_peaks) and only override technique-specific metadata and defaults.
