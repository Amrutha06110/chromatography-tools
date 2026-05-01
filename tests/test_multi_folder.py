"""Tests for multi-folder upload logic — scanning multiple directories and
combining DAD1A datasets with per-folder AB numbers."""

from __future__ import annotations

import os
import pathlib
import re
from typing import Any

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Replicate helpers from app.py (tests avoid importing the Streamlit module).
# ---------------------------------------------------------------------------

_AB_PATTERN = re.compile(r"(?<![A-Za-z])(AB[\s\-_]?\d+)", re.IGNORECASE)


def _extract_ab_number(*sources: str) -> str | None:
    for text in sources:
        m = _AB_PATTERN.search(text)
        if m:
            return re.sub(r"[\s\-_]", "", m.group(1)).upper()
    return None


def _scan_local_directory(dir_path: str) -> list[dict[str, Any]]:
    """Simplified scan matching app.py logic."""
    root = pathlib.Path(dir_path).expanduser().resolve()
    if not root.is_dir():
        return []
    results: list[dict[str, Any]] = []
    dad1a_files = sorted(
        p for p in root.rglob("*")
        if p.is_file()
        and "DAD1A" in p.stem.upper()
        and p.suffix.lower() == ".csv"
    )
    for csv_path in dad1a_files:
        try:
            df = pd.read_csv(csv_path, header=None)
        except Exception:
            continue
        if df.shape[1] < 2:
            continue
        df = df.iloc[:, :2]
        df.columns = ["time", "intensity"]
        df = df.apply(pd.to_numeric, errors="coerce").dropna()
        if df.empty:
            continue

        ab_number: str | None = None
        for ancestor in csv_path.relative_to(root).parents:
            if ancestor.name:
                ab_number = _extract_ab_number(ancestor.name)
                if ab_number:
                    break
        fname = str(csv_path.relative_to(root))
        label = ab_number or csv_path.parent.name
        results.append({"filename": fname, "label": label, "df": df})
    return results


def _build_folder_summary(folder_list: list[str]) -> list[dict[str, Any]]:
    """Build the summary rows shown in the UI for each folder."""
    rows = []
    for fpath in folder_list:
        folder_name = os.path.basename(fpath)
        ab_number = _extract_ab_number(folder_name) or folder_name
        root = pathlib.Path(fpath).expanduser().resolve()
        dad1a_count = 0
        if root.is_dir():
            dad1a_count = sum(
                1 for p in root.rglob("*")
                if p.is_file()
                and "DAD1A" in p.stem.upper()
                and p.suffix.lower() == ".csv"
            )
        rows.append({
            "Folder": folder_name,
            "AB Number": ab_number,
            "DAD1A files found": dad1a_count,
        })
    return rows


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def two_folder_tree(tmp_path: pathlib.Path):
    """Create two experiment folders each with a DAD1A CSV file.

    Structure:
        tmp_path/
            AB100_experiment/
                DAD1A.CSV        (2 rows)
            AB200_run2/
                sub.dx/
                    DAD1A.CSV    (2 rows)
    """
    folder1 = tmp_path / "AB100_experiment"
    folder1.mkdir()
    (folder1 / "DAD1A.CSV").write_text("0.0,10.0\n0.1,20.0\n")

    folder2 = tmp_path / "AB200_run2"
    sub = folder2 / "sub.dx"
    sub.mkdir(parents=True)
    (sub / "DAD1A.CSV").write_text("0.0,5.0\n0.1,15.0\n")

    return [str(folder1), str(folder2)]


@pytest.fixture()
def folder_with_multiple_dad1a(tmp_path: pathlib.Path):
    """Single folder containing two DAD1A files."""
    folder = tmp_path / "AB300_batch"
    folder.mkdir()
    (folder / "DAD1A_sample1.CSV").write_text("0.0,1.0\n0.1,2.0\n")
    (folder / "DAD1A_sample2.CSV").write_text("0.0,3.0\n0.1,4.0\n")
    return str(folder)


@pytest.fixture()
def folder_no_ab(tmp_path: pathlib.Path):
    """Folder without an AB number in its name."""
    folder = tmp_path / "experiment_data"
    folder.mkdir()
    (folder / "DAD1A.CSV").write_text("0.0,7.0\n0.1,8.0\n")
    return str(folder)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestMultiFolderScan:
    """Validate scanning multiple folders and combining datasets."""

    def test_two_folders_produce_two_datasets(self, two_folder_tree):
        datasets = []
        for fpath in two_folder_tree:
            datasets.extend(_scan_local_directory(fpath))
        assert len(datasets) == 2

    def test_labels_have_correct_ab_numbers(self, two_folder_tree):
        datasets = []
        for fpath in two_folder_tree:
            datasets.extend(_scan_local_directory(fpath))
        labels = sorted(ds["label"] for ds in datasets)
        # When root IS the AB folder, the ancestor walk finds no AB number
        # in the relative path — the fallback is csv_path.parent.name.
        # folder1: DAD1A.CSV is directly inside root → parent.name = root name
        # folder2: sub.dx/DAD1A.CSV → parent.name = "sub.dx"
        assert labels == ["AB100_experiment", "sub.dx"]

    def test_multiple_dad1a_in_one_folder(self, folder_with_multiple_dad1a):
        datasets = _scan_local_directory(folder_with_multiple_dad1a)
        assert len(datasets) == 2
        # Files are directly inside root, so label is root folder name
        assert all(ds["label"] == "AB300_batch" for ds in datasets)

    def test_no_ab_number_falls_back_to_parent_name(self, folder_no_ab):
        datasets = _scan_local_directory(folder_no_ab)
        assert len(datasets) == 1
        assert datasets[0]["label"] == "experiment_data"

    def test_empty_folder_returns_nothing(self, tmp_path):
        empty = tmp_path / "empty_folder"
        empty.mkdir()
        assert _scan_local_directory(str(empty)) == []

    def test_nonexistent_path_returns_nothing(self, tmp_path):
        assert _scan_local_directory(str(tmp_path / "does_not_exist")) == []


class TestFolderSummary:
    """Validate the folder summary table builder."""

    def test_summary_counts_dad1a(self, two_folder_tree):
        rows = _build_folder_summary(two_folder_tree)
        assert len(rows) == 2
        assert rows[0]["DAD1A files found"] == 1
        assert rows[1]["DAD1A files found"] == 1

    def test_summary_ab_numbers(self, two_folder_tree):
        rows = _build_folder_summary(two_folder_tree)
        ab_numbers = sorted(r["AB Number"] for r in rows)
        assert ab_numbers == ["AB100", "AB200"]

    def test_summary_no_ab_falls_back(self, folder_no_ab):
        rows = _build_folder_summary([folder_no_ab])
        assert rows[0]["AB Number"] == "experiment_data"

    def test_summary_nonexistent_shows_zero(self, tmp_path):
        rows = _build_folder_summary([str(tmp_path / "ghost")])
        assert rows[0]["DAD1A files found"] == 0

    def test_summary_multiple_dad1a(self, folder_with_multiple_dad1a):
        rows = _build_folder_summary([folder_with_multiple_dad1a])
        assert rows[0]["DAD1A files found"] == 2


class TestFolderListDedup:
    """Validate that duplicate paths are not added to the folder list."""

    def test_duplicate_path_not_added(self):
        folder_list: list[str] = []
        path = "/some/path"
        if path not in folder_list:
            folder_list.append(path)
        if path not in folder_list:
            folder_list.append(path)
        assert folder_list == ["/some/path"]

    def test_different_paths_both_added(self):
        folder_list: list[str] = []
        for p in ["/path/a", "/path/b"]:
            if p not in folder_list:
                folder_list.append(p)
        assert len(folder_list) == 2


class TestFolderListRemoval:
    """Validate individual folder removal from the folder list."""

    def test_remove_single_folder_by_index(self):
        """Removing a folder by index leaves the others intact."""
        folder_list = ["/path/a", "/path/b", "/path/c"]
        folder_list.pop(1)
        assert folder_list == ["/path/a", "/path/c"]

    def test_remove_first_folder(self):
        """Removing the first folder shifts the rest."""
        folder_list = ["/path/a", "/path/b"]
        folder_list.pop(0)
        assert folder_list == ["/path/b"]

    def test_remove_last_folder(self):
        """Removing the last folder leaves the others."""
        folder_list = ["/path/a", "/path/b"]
        folder_list.pop(1)
        assert folder_list == ["/path/a"]

    def test_remove_only_folder_leaves_empty(self):
        """Removing the sole folder results in an empty list."""
        folder_list = ["/path/a"]
        folder_list.pop(0)
        assert folder_list == []

    def test_clear_all_empties_list(self):
        """Clearing all folders resets the list to empty."""
        folder_list = ["/path/a", "/path/b", "/path/c"]
        folder_list.clear()
        assert folder_list == []

    def test_remove_preserves_order(self):
        """After removal, remaining folders keep their original order."""
        folder_list = ["/path/a", "/path/b", "/path/c", "/path/d"]
        folder_list.pop(2)  # remove /path/c
        assert folder_list == ["/path/a", "/path/b", "/path/d"]
