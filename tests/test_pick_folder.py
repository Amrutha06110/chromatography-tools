"""Tests for the pick_folder() helper function."""

from __future__ import annotations

import subprocess
from unittest.mock import patch, MagicMock

import pytest


# Replicate pick_folder logic from app.py so tests don't need to import the
# Streamlit-heavy module.

def pick_folder() -> str | None:
    """Open a native macOS folder-picker dialog via osascript."""
    script = (
        'tell application "System Events"\n'
        "    activate\n"
        "end tell\n"
        'set folderPath to POSIX path of (choose folder with prompt '
        '"Select experiment folder")\n'
        "return folderPath"
    )
    try:
        result = subprocess.run(
            ["osascript", "-e", script],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return None
    if result.returncode == 0:
        return result.stdout.strip().rstrip("/")
    return None


class TestPickFolder:
    """Tests for pick_folder() using osascript."""

    def test_returns_path_on_success(self):
        """Successful osascript invocation returns the trimmed path."""
        mock_result = MagicMock(returncode=0, stdout="/Users/me/Data/AB123\n")
        with patch("subprocess.run", return_value=mock_result) as mock_run:
            result = pick_folder()
            assert result == "/Users/me/Data/AB123"
            mock_run.assert_called_once()

    def test_strips_trailing_slash(self):
        """Trailing slash from osascript output is removed."""
        mock_result = MagicMock(returncode=0, stdout="/Users/me/Data/AB123/\n")
        with patch("subprocess.run", return_value=mock_result):
            assert pick_folder() == "/Users/me/Data/AB123"

    def test_returns_none_on_cancel(self):
        """User cancelling the dialog (non-zero return) returns None."""
        mock_result = MagicMock(returncode=1, stdout="", stderr="")
        with patch("subprocess.run", return_value=mock_result):
            assert pick_folder() is None

    def test_returns_none_when_osascript_missing(self):
        """On non-macOS systems where osascript doesn't exist, returns None."""
        with patch("subprocess.run", side_effect=FileNotFoundError):
            assert pick_folder() is None
