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

