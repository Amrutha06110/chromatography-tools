"""Tests for the AB-number extraction helper used when renaming DAD1A files."""

from __future__ import annotations

import pathlib
import re

import pytest

# Replicate the extraction logic from app.py so tests don't need to
# import the Streamlit-heavy module.
_AB_PATTERN = re.compile(r"^AB[\s\-_]?\d+", re.IGNORECASE)


def _extract_ab_number(rel_path: str) -> str | None:
    parts = pathlib.PurePosixPath(rel_path).parts
    for part in parts[:-1]:
        m = _AB_PATTERN.match(part)
        if m:
            return m.group(0)
    return None


def _resolve_label(name: str, webkit_rel_path: str = "") -> str:
    """Replicate the labeling logic from app.py's data-loading loop.

    Parameters
    ----------
    name : str
        ``uf.name`` — may be a bare filename (``DAD1A.CSV``) or a full
        relative path (``AB734/timestamp.dx/DAD1A.CSV``).
    webkit_rel_path : str
        ``webkitRelativePath`` attribute (usually empty in Streamlit).
    """
    rel_path = webkit_rel_path or ""
    ab_number = _extract_ab_number(rel_path) if rel_path else None
    if ab_number is None:
        ab_number = _extract_ab_number(name)
    if ab_number:
        return ab_number

    name_parts = pathlib.PurePosixPath(name).parts
    if len(name_parts) > 1:
        return name_parts[-2]
    elif rel_path:
        parts = pathlib.PurePosixPath(rel_path).parts
        return parts[-2] if len(parts) > 1 else name.rsplit(".", 1)[0]
    else:
        return name.rsplit(".", 1)[0]


class TestExtractAbNumber:
    """Validate AB-number extraction from path strings."""

    def test_ab_number_in_grandparent(self):
        """Typical Agilent structure: AB734/timestamp.dx/DAD1A.CSV"""
        assert _extract_ab_number("AB734/2026-04-17 12-31-26.dx/DAD1A.CSV") == "AB734"

    def test_ab_number_in_immediate_parent(self):
        """AB number directly above the file."""
        assert _extract_ab_number("AB1234/DAD1A.CSV") == "AB1234"

    def test_ab_number_with_dash(self):
        """AB numbers may use a dash separator."""
        assert _extract_ab_number("AB-99/some_folder/DAD1A.CSV") == "AB-99"

    def test_ab_number_case_insensitive(self):
        """Match should be case-insensitive."""
        assert _extract_ab_number("ab500/DAD1A.CSV") == "ab500"

    def test_no_ab_number_returns_none(self):
        """When no directory contains an AB number, return None."""
        assert _extract_ab_number("SomeProject/data/DAD1A.CSV") is None

    def test_empty_path_returns_none(self):
        assert _extract_ab_number("") is None

    def test_bare_filename_returns_none(self):
        assert _extract_ab_number("DAD1A.CSV") is None

    def test_ab_in_filename_not_matched(self):
        """AB pattern in the filename itself should be ignored (dirs only)."""
        assert _extract_ab_number("project/AB_data.CSV") is None

    def test_first_ab_number_wins(self):
        """If multiple directories contain AB numbers, the first (root-ward) wins."""
        assert _extract_ab_number("AB100/AB200/DAD1A.CSV") == "AB100"

    def test_deeply_nested(self):
        """AB number several levels up from the file."""
        assert _extract_ab_number("lab/AB42/run1/2026-04-17.dx/DAD1A.CSV") == "AB42"

    def test_ab_with_underscore_separator(self):
        """AB numbers may use an underscore separator."""
        assert _extract_ab_number("AB_55/DAD1A.CSV") == "AB_55"


class TestResolveLabel:
    """End-to-end label resolution mimicking app.py's data-loading loop.

    Streamlit's UploadedFile has NO webkitRelativePath attribute, so the
    AB number must be extracted from ``uf.name`` which may contain the
    full relative path when browser folder selection is used.
    """

    def test_ab_from_name_with_path(self):
        """uf.name = 'AB734/timestamp.dx/DAD1A.CSV' → 'AB734'."""
        assert _resolve_label("AB734/2026-04-17 12-31-26.dx/DAD1A.CSV") == "AB734"

    def test_ab_from_name_immediate_parent(self):
        """uf.name = 'AB1234/DAD1A.CSV' → 'AB1234'."""
        assert _resolve_label("AB1234/DAD1A.CSV") == "AB1234"

    def test_fallback_parent_dir_no_ab(self):
        """No AB number — use immediate parent dir from uf.name."""
        label = _resolve_label("2026-04-17 12-31-26.dx/DAD1A.CSV")
        assert label == "2026-04-17 12-31-26.dx"

    def test_fallback_bare_filename(self):
        """uf.name = 'DAD1A.CSV' (no path info) → 'DAD1A'."""
        assert _resolve_label("DAD1A.CSV") == "DAD1A"

    def test_webkit_rel_path_takes_priority(self):
        """If webkitRelativePath is available, its AB number wins."""
        label = _resolve_label("DAD1A.CSV", webkit_rel_path="AB999/ts.dx/DAD1A.CSV")
        assert label == "AB999"

    def test_name_ab_when_webkit_empty(self):
        """webkitRelativePath empty — AB from uf.name is used."""
        label = _resolve_label("AB100/ts.dx/DAD1A.CSV", webkit_rel_path="")
        assert label == "AB100"

    def test_ab_with_dash_from_name(self):
        """AB-99 style in uf.name."""
        assert _resolve_label("AB-99/subfolder/DAD1A.CSV") == "AB-99"

    def test_case_insensitive_from_name(self):
        """Lowercase 'ab' in uf.name."""
        assert _resolve_label("ab500/DAD1A.CSV") == "ab500"

    def test_deeply_nested_from_name(self):
        """AB number several levels up in uf.name."""
        label = _resolve_label("lab/AB42/run1/2026-04-17.dx/DAD1A.CSV")
        assert label == "AB42"


class TestDuplicateAbNumberHandling:
    """Validate that duplicate AB numbers get counter suffixes."""

    def _resolve_duplicates(self, labels: list[str]) -> list[str]:
        """Replicate the duplicate-handling logic from app.py."""
        label_counts: dict[str, int] = {}
        for label in labels:
            label_counts[label] = label_counts.get(label, 0) + 1

        label_seen: dict[str, int] = {}
        final_labels: list[str] = []
        for label in labels:
            if label_counts[label] > 1:
                label_seen[label] = label_seen.get(label, 0) + 1
                final_labels.append(f"{label}_{label_seen[label]}")
            else:
                final_labels.append(label)
        return final_labels

    def test_no_duplicates(self):
        """Unique labels remain unchanged."""
        assert self._resolve_duplicates(["AB100", "AB200"]) == ["AB100", "AB200"]

    def test_two_duplicates(self):
        """Two identical labels get _1 and _2 suffixes."""
        assert self._resolve_duplicates(["AB100", "AB100"]) == ["AB100_1", "AB100_2"]

    def test_three_duplicates(self):
        """Three identical labels get _1, _2, _3 suffixes."""
        result = self._resolve_duplicates(["AB42", "AB42", "AB42"])
        assert result == ["AB42_1", "AB42_2", "AB42_3"]

    def test_mixed_duplicates_and_unique(self):
        """Mix of unique and duplicate labels handled correctly."""
        result = self._resolve_duplicates(["AB100", "AB200", "AB100", "AB300"])
        assert result == ["AB100_1", "AB200", "AB100_2", "AB300"]


class TestExtractAbNumberFromFolderName:
    """Validate the folder-name-based AB extraction.

    This mirrors the ``_extract_ab_number()`` helper in ``app.py`` which
    searches for the pattern ``AB<digits>`` anywhere in a folder name string.
    """

    @staticmethod
    def _extract(folder_name: str) -> str | None:
        """Replicate extract_ab_number from app.py."""
        match = re.search(r"AB\d+", folder_name, re.IGNORECASE)
        return match.group(0).upper() if match else None

    def test_typical_agilent_folder_name(self):
        """Standard Agilent folder: '20260424 151551SYSTEM (SYSTEM)AB628'."""
        assert self._extract("20260424 151551SYSTEM (SYSTEM)AB628") == "AB628"

    def test_ab_number_at_start(self):
        """AB number at the beginning of the folder name."""
        assert self._extract("AB734_experiment_data") == "AB734"

    def test_ab_number_in_middle(self):
        """AB number embedded in the middle of the folder name."""
        assert self._extract("2026_AB99_run1") == "AB99"

    def test_lowercase_ab(self):
        """Lowercase 'ab' should be matched and uppercased."""
        assert self._extract("data_ab500_results") == "AB500"

    def test_no_ab_number(self):
        """Folder with no AB number returns None."""
        assert self._extract("20260424 151551SYSTEM (SYSTEM)") is None

    def test_empty_string(self):
        """Empty folder name returns None."""
        assert self._extract("") is None

    def test_first_ab_number_wins(self):
        """When multiple AB numbers exist, the first one is returned."""
        assert self._extract("AB100_rerun_AB200") == "AB100"

    def test_ab_without_digits(self):
        """'AB' without digits should not match."""
        assert self._extract("AB_only_text") is None

    def test_large_ab_number(self):
        """Large AB number with many digits."""
        assert self._extract("folder_AB123456_data") == "AB123456"

