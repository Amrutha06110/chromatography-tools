"""Streamlit application for interactive chromatogram peak analysis.

Upload one or more CSV files containing ``(time, intensity)`` columns and the
app performs smoothing, peak detection, trapezoidal integration, and
relative-area calculation.  Supports multiple chromatography techniques
(HPLC, GC, SEC/GPC, Ion Chromatography) through a proper OOP class hierarchy.
"""

from __future__ import annotations

import io
import json
import os
import pathlib
import re
from typing import Any

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    import tkinter  # noqa: F401
    _HAS_TKINTER = True
except ImportError:
    _HAS_TKINTER = False

from chromatography.core import (
    BaseChromatogram,
    GCChromatogram,
    HPLCChromatogram,
    IonChromatogram,
    SECChromatogram,
)

# ------------------------------------------------------------------
# Technique → class mapping
# ------------------------------------------------------------------

TECHNIQUE_MAP: dict[str, type[BaseChromatogram]] = {
    "HPLC": HPLCChromatogram,
    "GC": GCChromatogram,
    "SEC/GPC": SECChromatogram,
    "Ion Chromatography": IonChromatogram,
}


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

_AB_PATTERN = re.compile(r"(?<![A-Za-z])(AB[\s\-_]?\d+)", re.IGNORECASE)


def _extract_ab_number(*sources: str) -> str | None:
    """Return the first ``AB<digits>`` token found in any of *sources*.

    Each source is a free-form string — it can be a filename, a folder
    name, or a relative path.  Supports optional separators between
    ``AB`` and the digits (e.g. ``AB734``, ``AB-734``, ``AB_734``).
    Returns the normalised, uppercased identifier (e.g. ``"AB734"``)
    on the first match, or ``None`` when nothing is found.

    Sources are checked left-to-right, so callers should pass
    higher-priority strings first.
    """
    for text in sources:
        m = _AB_PATTERN.search(text)
        if m:
            # Normalise: strip any separator between "AB" and the digits.
            return re.sub(r"[\s\-_]", "", m.group(1)).upper()
    return None


def pick_folder() -> str:
    """Open a native OS folder picker dialog and return the selected path.

    Only call this when ``_HAS_TKINTER`` is ``True``.
    """
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes("-topmost", True)
    folder_path = filedialog.askdirectory(title="Select experiment folder")
    root.destroy()
    return folder_path


def extract_ab_number(folder_name: str) -> str | None:
    """Extract an AB<digits> identifier from a folder name string.

    The folder name format is typically something like
    ``"20260424 151551SYSTEM (SYSTEM)AB628"``.
    Returns the normalised, uppercased identifier (e.g. ``"AB628"``),
    or ``None`` when no match is found.
    """
    match = re.search(r"AB\d+", folder_name, re.IGNORECASE)
    return match.group(0).upper() if match else None


def read_two_col_csv(uploaded_file: Any) -> pd.DataFrame | None:
    """Read a 2-column CSV into a ``time`` / ``intensity`` DataFrame."""
    try:
        df = pd.read_csv(uploaded_file, header=None)
    except Exception:
        return None
    if df.shape[1] < 2:
        return None
    df = df.iloc[:, :2]
    df.columns = ["time", "intensity"]
    df = df.apply(pd.to_numeric, errors="coerce").dropna()
    if df.empty:
        return None
    return df


def _scan_local_directory(dir_path: str) -> list[dict[str, Any]]:
    """Scan a local directory for ``DAD1A*.csv`` files and return datasets.

    Each dataset dict has keys ``filename``, ``label``, and ``df``.  The
    label is the AB identifier extracted from the file's ancestor
    directories, falling back to the parent directory name.

    **Security note:** the *dir_path* is intentionally user-provided —
    this is a local-first analysis tool where the user pastes the path
    to their own data.  Only CSV files are parsed (via ``pd.read_csv``);
    non-CSV content is silently skipped.
    """
    root = pathlib.Path(dir_path).expanduser().resolve()  # noqa: S108
    if not root.is_dir():
        return []

    results: list[dict[str, Any]] = []
    # Find DAD1A CSV files anywhere under *root* (case-insensitive).
    dad1a_files = sorted(
        p for p in root.rglob("*")
        if p.is_file()
        and "DAD1A" in p.stem.upper()
        and p.suffix.lower() == ".csv"
    )
    for csv_path in dad1a_files:
        df = read_two_col_csv(csv_path)
        if df is None:
            continue

        # Walk up from the file's parent toward *root* looking for an AB
        # identifier in any ancestor directory name.
        ab_number: str | None = None
        for ancestor in csv_path.relative_to(root).parents:
            if ancestor.name:  # skip the empty root "."
                ab_number = _extract_ab_number(ancestor.name)
                if ab_number:
                    break

        fname = str(csv_path.relative_to(root))
        label = ab_number or csv_path.parent.name
        results.append({"filename": fname, "label": label, "df": df})

    return results


VIRIDIS_COLORS = [
    "#440154", "#482878", "#3e4989", "#31688e", "#26828e",
    "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725",
]


def _color_palette(n: int) -> list[str]:
    if n <= len(VIRIDIS_COLORS):
        step = max(1, len(VIRIDIS_COLORS) // n)
        return [VIRIDIS_COLORS[i * step % len(VIRIDIS_COLORS)] for i in range(n)]
    return [VIRIDIS_COLORS[i % len(VIRIDIS_COLORS)] for i in range(n)]


def _hex_to_rgba(hex_color: str, alpha: float = 0.15) -> str:
    """Convert a hex color string to an rgba() string with the given alpha."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _label_yshifts(peaks: list, y_max: float) -> list[float]:
    """Return per-peak y-shift values that avoid vertical label collisions.

    Labels are placed above each peak's height.  When two neighbouring
    labels would be too close together (within *threshold* of each other
    on the y-axis), the lower one is nudged upward so they don't overlap.
    """
    if not peaks:
        return []

    # Base offset above each peak (in data units)
    base_offset = 0.04 * y_max
    min_gap = 0.06 * y_max  # minimum vertical gap between labels

    # Work with (index, peak_height) pairs sorted by x position
    items = sorted(enumerate(peaks), key=lambda t: t[1].time)
    offsets = [base_offset] * len(peaks)

    for i in range(1, len(items)):
        cur_idx, cur_pk = items[i]
        prev_idx, prev_pk = items[i - 1]
        cur_y = cur_pk.height + offsets[cur_idx]
        prev_y = prev_pk.height + offsets[prev_idx]
        if abs(cur_y - prev_y) < min_gap:
            # Nudge the current label upward
            offsets[cur_idx] = (prev_y + min_gap) - cur_pk.height

    return offsets


# ------------------------------------------------------------------
# Streamlit UI
# ------------------------------------------------------------------

st.set_page_config(page_title="Chromatography Peak Analysis", layout="wide")
st.title("Chromatography Peak Analysis")

# ---- Sidebar ----
with st.sidebar:
    st.header("Data")

    data_dir = st.text_input(
        "Data directory (local path)",
        value="",
        help=(
            "Paste the full path to the folder that contains your sample "
            "subfolders (e.g. '/Users/you/Data/SequenceRun'). The app "
            "recursively finds DAD1A CSV files and automatically labels "
            "each one with the AB number from its parent folder name."
        ),
    )

    if _HAS_TKINTER:
        if st.button("📁 Select Experiment Folder"):
            folder_path = pick_folder()
            if folder_path:
                st.session_state["chrom_folder_path"] = folder_path

    st.header("Technique")
    technique = st.selectbox(
        "Chromatography type",
        list(TECHNIQUE_MAP.keys()),
        help="Selects technique-specific defaults (axis labels, smoothing window, etc.)",
    )

    st.header("Display")
    display_mode = st.radio("Mode", ["Overlay all", "Separate panels"], horizontal=True)
    show_peaks = st.checkbox("Show peaks", value=True)
    show_bounds = st.checkbox("Show peak bounds", value=False)
    label_peaks = st.checkbox("Label peaks", value=True)
    show_area_pct = st.checkbox("Include % area in labels", value=True)

    st.header("Peak Picking")
    min_height_pct = st.slider("Min height (% of max)", 0, 100, 10, step=1)

    st.header("Manual Peak Bounds")
    edit_bounds = st.checkbox(
        "Edit peak bounds",
        value=False,
        help=(
            "Enable to manually adjust integration bounds. "
            "Click on the plot to set a bound for the nearest peak. "
            "Clicks left of a peak's apex set the left bound; "
            "clicks to the right set the right bound."
        ),
    )

    st.header("Export")
    plot_fmt = st.selectbox("Save plot as", ["PNG", "SVG", "PDF"])
    plot_w = st.number_input("Width (px)", value=1000, min_value=200, step=50)
    plot_h = st.number_input("Height (px)", value=600, min_value=200, step=50)


# ---- Load data ----
datasets: list[dict[str, Any]] = []

# Priority 1: local directory path (reads filesystem directly, so all
# ancestor folder names — including AB numbers — are available).
if data_dir and data_dir.strip():
    datasets = _scan_local_directory(data_dir.strip())

# Priority 2: tkinter folder picker (native OS dialog gives full path,
# so the AB number can be extracted directly from the parent folder name).
if not datasets:
    folder_path = st.session_state.get("chrom_folder_path", None)
    if folder_path and os.path.isdir(folder_path):
        folder_name = os.path.basename(folder_path)
        ab_number = extract_ab_number(folder_name)

        if ab_number:
            st.success(f"Folder: `{folder_name}` → AB number: **{ab_number}**")
        else:
            st.warning(f"No AB number found in folder name: `{folder_name}`")
            ab_number = folder_name  # fall back to full folder name

        # Find all DAD1A CSV files inside the folder
        all_files = os.listdir(folder_path)
        dad1a_paths = sorted([
            os.path.join(folder_path, f)
            for f in all_files
            if "DAD1A" in f.upper() and f.lower().endswith(".csv")
        ])

        st.info(
            f"Found **{len(dad1a_paths)}** DAD1A file(s) out of "
            f"{len(all_files)} total file(s) in folder."
        )

        if not dad1a_paths:
            st.error("No DAD1A CSV files found in this folder.")
        else:
            # Assign nickname: AB number + counter if multiple DAD1A files
            nicknames: dict[str, str] = {}
            for i, path in enumerate(dad1a_paths):
                fname = os.path.basename(path)
                nicknames[fname] = (
                    ab_number if len(dad1a_paths) == 1
                    else f"{ab_number}_{i + 1}"
                )

            # Show preview table
            preview = pd.DataFrame([
                {
                    "File": os.path.basename(p),
                    "Nickname": nicknames[os.path.basename(p)],
                }
                for p in dad1a_paths
            ])
            st.dataframe(preview, use_container_width=True, hide_index=True)

            # Load each DAD1A file using the existing CSV parsing function
            for path in dad1a_paths:
                fname = os.path.basename(path)
                nickname = nicknames[fname]
                df = read_two_col_csv(path)
                if df is not None:
                    datasets.append({
                        "filename": nickname,
                        "label": nickname,
                        "df": df,
                    })

if not datasets:
    _hint = (
        "Paste a local data directory path in the sidebar to automatically "
        "find DAD1A files and label them with AB numbers"
    )
    if _HAS_TKINTER:
        _hint += ", or click '📁 Select Experiment Folder' to pick a folder."
    else:
        _hint += "."
    st.info(_hint)
    st.stop()

# ---- Nicknames ----
with st.sidebar:
    st.header("Nicknames")
    for ds in datasets:
        ds["label"] = st.text_input(f"Nickname for {ds['filename']}", value=ds["label"], key=f"nick_{ds['filename']}")


# ---- Session state for manual peak-bound overrides ----
# Key: (filename, peak_time) → {"left_time": float | None, "right_time": float | None}
if "manual_bounds" not in st.session_state:
    st.session_state.manual_bounds: dict[tuple[str, float], dict[str, float | None]] = {}


# ---- Build Chromatogram objects (using the selected technique class) ----
ChromClass = TECHNIQUE_MAP[technique]
chromatograms: list[dict[str, Any]] = []
for ds in datasets:
    ch = ChromClass(ds["df"]["time"].values, ds["df"]["intensity"].values)
    ch.smooth()  # uses technique-specific defaults

    ref_signal = ch.smoothed if ch.smoothed is not None else ch.intensity
    min_h_abs = (min_height_pct / 100) * float(np.nanmax(ref_signal))

    ch.find_peaks(min_height=min_h_abs)

    # Apply any saved manual bounds before integration
    for (fname, peak_t), bounds in st.session_state.manual_bounds.items():
        if fname == ds["filename"] and ch.peaks:
            left = bounds.get("left_time")
            right = bounds.get("right_time")
            if left is not None or right is not None:
                ch.set_manual_bounds(
                    peak_t, left_time=left, right_time=right,
                )

    if ch.peaks:
        ch.integrate_peaks()

    chromatograms.append({"label": ds["label"], "filename": ds["filename"], "chrom": ch})

# Axis labels from the first chromatogram instance
x_label = chromatograms[0]["chrom"].x_axis_label
y_label = chromatograms[0]["chrom"].y_axis_label

# ---- Plotting ----
colors = _color_palette(len(chromatograms))
# Mapping from Plotly trace index → chromatogram list index, used in
# separate-panels mode to identify which panel was clicked.
trace_to_chrom_idx: dict[int, int] = {}

if display_mode == "Overlay all" or len(chromatograms) == 1:
    fig = go.Figure()
    for idx, entry in enumerate(chromatograms):
        ch = entry["chrom"]
        fig.add_trace(go.Scatter(
            x=ch.time, y=ch.intensity, mode="lines",
            name=entry["label"], line=dict(color=colors[idx], width=2),
        ))
        if show_peaks and ch.peaks:
            if show_bounds:
                for pi, pk in enumerate(ch.peaks):
                    s = max(0, min(pk.start, len(ch.time) - 1))
                    e = max(0, min(pk.end, len(ch.time) - 1))
                    alpha = 0.25 if pi % 2 == 0 else 0.12
                    fig.add_trace(go.Scatter(
                        x=ch.time[s : e + 1],
                        y=ch.intensity[s : e + 1],
                        fill="tozeroy",
                        fillcolor=_hex_to_rgba(colors[idx], alpha=alpha),
                        line=dict(width=0),
                        showlegend=False,
                        hoverinfo="skip",
                    ))
            if label_peaks:
                y_max = float(np.nanmax(ch.intensity))
                yshifts = _label_yshifts(ch.peaks, y_max)
                for pi, pk in enumerate(ch.peaks):
                    txt = f"t={pk.time:.2f}"
                    if show_area_pct and pk.rel_area_pct is not None:
                        txt += f", {pk.rel_area_pct:.1f}%"
                    fig.add_annotation(x=pk.time, y=pk.height + yshifts[pi], text=txt, showarrow=False, font=dict(size=10))

    fig.update_layout(
        xaxis_title=x_label, yaxis_title=y_label,
        template="plotly_white", width=plot_w, height=plot_h,
        legend_title="Sample",
    )
else:
    n = len(chromatograms)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=[e["label"] for e in chromatograms])
    trace_count = 0
    for idx, entry in enumerate(chromatograms):
        r, c = divmod(idx, ncols)
        ch = entry["chrom"]
        trace_to_chrom_idx[trace_count] = idx
        fig.add_trace(go.Scatter(
            x=ch.time, y=ch.intensity, mode="lines",
            name=entry["label"], line=dict(color=colors[idx], width=2), showlegend=False,
        ), row=r + 1, col=c + 1)
        trace_count += 1
        if show_peaks and ch.peaks:
            if show_bounds:
                for pi, pk in enumerate(ch.peaks):
                    s = max(0, min(pk.start, len(ch.time) - 1))
                    e = max(0, min(pk.end, len(ch.time) - 1))
                    alpha = 0.25 if pi % 2 == 0 else 0.12
                    fig.add_trace(go.Scatter(
                        x=ch.time[s : e + 1],
                        y=ch.intensity[s : e + 1],
                        fill="tozeroy",
                        fillcolor=_hex_to_rgba(colors[idx], alpha=alpha),
                        line=dict(width=0),
                        showlegend=False,
                        hoverinfo="skip",
                    ), row=r + 1, col=c + 1)
                    trace_count += 1
            if label_peaks:
                y_max = float(np.nanmax(ch.intensity))
                yshifts = _label_yshifts(ch.peaks, y_max)
                for pi, pk in enumerate(ch.peaks):
                    txt = f"t={pk.time:.2f}"
                    if show_area_pct and pk.rel_area_pct is not None:
                        txt += f", {pk.rel_area_pct:.1f}%"
                    fig.add_annotation(x=pk.time, y=pk.height + yshifts[pi], text=txt, showarrow=False, font=dict(size=10), xref=f"x{idx + 1}" if idx else "x", yref=f"y{idx + 1}" if idx else "y")

    fig.update_layout(template="plotly_white", width=plot_w, height=plot_h)

# ---- Render chart (with click handling when editing bounds) ----
if edit_bounds:
    # Disable drag-based selection so that a simple single click is all
    # the user needs.  "event+select" ensures clicks on data points fire
    # selection events that Streamlit can capture.
    fig.update_layout(dragmode=False, clickmode="event+select")

    # Show a crosshair cursor instead of the default pan/zoom arrow
    # when editing bounds.  The Plotly drag-layer element (.nsewdrag)
    # controls the cursor shown over the plot area.
    st.markdown(
        "<style>.nsewdrag { cursor: crosshair !important; }</style>",
        unsafe_allow_html=True,
    )

    # Add invisible markers to the main line traces so that clicks near
    # the line snap to the nearest data-point — the same x-axis buffer
    # the hover already provides.
    for trace_data in fig.data:
        if getattr(trace_data, "mode", None) == "lines":
            trace_data.update(
                mode="lines+markers",
                marker=dict(size=8, opacity=0),
            )

    event = st.plotly_chart(
        fig,
        use_container_width=True,
        on_select="rerun",
        selection_mode=("points",),
        key="peak_bound_chart",
    )
    # Extract the clicked x-coordinate and curve number from the point selection.
    clicked_time: float | None = None
    clicked_curve: int | None = None
    if event and hasattr(event, "selection") and event.selection:
        sel = event.selection
        sel_points = sel.get("points", [])
        if sel_points:
            cx = sel_points[0].get("x")
            if cx is not None:
                clicked_time = float(cx)
            cn = sel_points[0].get("curve_number")
            if cn is not None:
                clicked_curve = int(cn)

    if clicked_time is not None:
        overlay_active = display_mode == "Overlay all" and len(chromatograms) > 1
        if overlay_active:
            # In overlay mode, simultaneously set bounds for the nearest
            # peak in every chromatogram at the same x-axis position.
            for entry in chromatograms:
                ch = entry["chrom"]
                if not ch.peaks:
                    continue
                nearest_pk = min(ch.peaks, key=lambda pk: abs(pk.time - clicked_time))
                key = (entry["filename"], nearest_pk.time)
                existing = st.session_state.manual_bounds.get(key, {})
                if clicked_time < nearest_pk.time:
                    existing["left_time"] = clicked_time
                else:
                    existing["right_time"] = clicked_time
                st.session_state.manual_bounds[key] = existing
            st.rerun()
        else:
            # Separate-panels mode: use the curve_number from the click
            # event to identify which panel was clicked, then set bounds
            # only for that chromatogram's nearest peak.
            target_entry = None
            if clicked_curve is not None and trace_to_chrom_idx and clicked_curve in trace_to_chrom_idx:
                target_entry = chromatograms[trace_to_chrom_idx[clicked_curve]]

            if target_entry is not None:
                ch = target_entry["chrom"]
                if ch.peaks:
                    nearest_pk = min(ch.peaks, key=lambda pk: abs(pk.time - clicked_time))
                    key = (target_entry["filename"], nearest_pk.time)
                    existing = st.session_state.manual_bounds.get(key, {})
                    if clicked_time < nearest_pk.time:
                        existing["left_time"] = clicked_time
                    else:
                        existing["right_time"] = clicked_time
                    st.session_state.manual_bounds[key] = existing
                    st.rerun()
            else:
                # Fallback when curve_number is unavailable (single file
                # or unexpected event shape): find the nearest peak
                # across all chromatograms.
                best_entry = None
                best_peak = None
                best_dist = float("inf")
                for entry in chromatograms:
                    ch = entry["chrom"]
                    for pk in ch.peaks:
                        d = abs(pk.time - clicked_time)
                        if d < best_dist:
                            best_dist = d
                            best_peak = pk
                            best_entry = entry

                if best_peak is not None and best_entry is not None:
                    key = (best_entry["filename"], best_peak.time)
                    existing = st.session_state.manual_bounds.get(key, {})
                    if clicked_time < best_peak.time:
                        existing["left_time"] = clicked_time
                    else:
                        existing["right_time"] = clicked_time
                    st.session_state.manual_bounds[key] = existing
                    st.rerun()
else:
    st.plotly_chart(fig, use_container_width=True)

# ---- Show active manual overrides and allow clearing ----
if st.session_state.manual_bounds:
    with st.sidebar:
        st.subheader("Active Bound Overrides")
        keys_to_remove: list[tuple[str, float]] = []
        for (fname, peak_t), bounds in st.session_state.manual_bounds.items():
            left_val = bounds.get("left_time")
            right_val = bounds.get("right_time")
            left_str = f"{left_val:.2f}" if left_val is not None else "auto"
            right_str = f"{right_val:.2f}" if right_val is not None else "auto"
            col1, col2 = st.columns([3, 1])
            col1.write(f"**{fname}** peak @ {peak_t:.2f}  \nL: {left_str} | R: {right_str}")
            if col2.button("✕", key=f"rm_{fname}_{peak_t}"):
                keys_to_remove.append((fname, peak_t))
        for k in keys_to_remove:
            del st.session_state.manual_bounds[k]
        if keys_to_remove:
            st.rerun()
        if st.button("Clear all overrides"):
            st.session_state.manual_bounds.clear()
            st.rerun()


# ---- Plot export ----
fmt_map = {"PNG": "png", "SVG": "svg", "PDF": "pdf"}


@st.cache_data(show_spinner=False)
def _render_plot_image(fig_json: str, fmt: str, width: int, height: int) -> bytes:
    """Render a Plotly figure to static image bytes (cached).

    Kaleido-based rendering is expensive (~1-2 s per call).  By caching
    on the serialised figure JSON + export settings, subsequent reruns
    with an unchanged figure return instantly.
    """
    restored_fig = go.Figure(json.loads(fig_json))
    buf = io.BytesIO()
    restored_fig.write_image(buf, format=fmt, width=width, height=height)
    return buf.getvalue()


@st.fragment
def _plot_download_button() -> None:
    """Render the download button for the plot image."""
    if edit_bounds:
        # While interactively editing bounds, skip the expensive Kaleido
        # static-image render (~1-2 s per call) so that the click → rerun
        # cycle stays responsive.  The download button reappears as soon
        # as the user un-checks "Edit peak bounds".
        return
    image_data = _render_plot_image(
        fig.to_json(), fmt_map[plot_fmt], plot_w, plot_h,
    )
    st.download_button(
        label=f"Download plot ({plot_fmt})",
        data=image_data,
        file_name=f"chromatograms.{fmt_map[plot_fmt]}",
        mime=f"image/{fmt_map[plot_fmt]}",
    )


_plot_download_button()


# ---- Summary CSV ----
@st.fragment
def _summary_csv_section() -> None:
    """Render the peak summary table and CSV download as an isolated fragment.

    Running inside ``@st.fragment`` prevents the download-button click from
    triggering an expensive full-page rerun (which would include the Kaleido
    static-image render when edit_bounds is off).
    """
    rows: list[dict[str, Any]] = []
    for entry in chromatograms:
        ch = entry["chrom"]
        base: dict[str, Any] = {
            "filename": entry["filename"],
            "nickname": entry["label"],
            "technique": ch.technique_name,
            "total_peaks": len(ch.peaks),
        }
        if ch.peaks:
            for i, pk in enumerate(ch.peaks, 1):
                base[f"peak_{i}_time"] = pk.time
                base[f"peak_{i}_rel_area"] = pk.rel_area_pct
        rows.append(base)

    summary_df = pd.DataFrame(rows)
    st.subheader("Peak Summary")
    st.dataframe(summary_df, use_container_width=True)

    csv_buf = io.StringIO()
    summary_df.to_csv(csv_buf, index=False)
    st.download_button(
        label="Download summary CSV",
        data=csv_buf.getvalue(),
        file_name="chromatogram_summary.csv",
        mime="text/csv",
    )


_summary_csv_section()
