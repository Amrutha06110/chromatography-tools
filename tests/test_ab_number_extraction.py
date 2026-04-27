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


class TestExtractAbNumber:
    """Validate AB-number extraction from webkitRelativePath strings."""

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

