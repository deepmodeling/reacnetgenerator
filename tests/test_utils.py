# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test utility helpers."""

import asyncio
import hashlib
from io import BytesIO

import pytest
import requests

from reacnetgenerator.utils import download_file


class _DownloadResponse:
    """Provide the response behavior used by downloader regression tests."""

    def __init__(self, content=b"", status_error=None):
        self.raw = BytesIO(content)
        self.status_error = status_error

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def raise_for_status(self):
        """Raise the configured HTTP status error, if any."""
        if self.status_error is not None:
            raise self.status_error


def _sha256(content):
    """Return the SHA256 digest for a test payload."""
    return hashlib.sha256(content).hexdigest()


def test_download_file_retries_after_http_error(tmp_path, mocker):
    """An unsuccessful HTTP response should fall back to the next URL."""
    content = b"valid trajectory"
    session = mocker.patch("reacnetgenerator.utils.requests.Session").return_value
    session.get.side_effect = [
        _DownloadResponse(status_error=requests.HTTPError("404 Not Found")),
        _DownloadResponse(content),
    ]
    output = tmp_path / "trajectory.dump"

    result = asyncio.run(
        download_file(["https://bad", "https://good"], str(output), _sha256(content))
    )

    assert result == str(output)
    assert output.read_bytes() == content
    assert session.get.call_count == 2


def test_download_file_retries_after_hash_mismatch(tmp_path, mocker):
    """A corrupt payload should be removed before trying another URL."""
    content = b"valid trajectory"
    session = mocker.patch("reacnetgenerator.utils.requests.Session").return_value
    session.get.side_effect = [
        _DownloadResponse(b"corrupt trajectory"),
        _DownloadResponse(content),
    ]
    output = tmp_path / "trajectory.dump"

    asyncio.run(
        download_file(
            ["https://corrupt", "https://good"], str(output), _sha256(content)
        )
    )

    assert output.read_bytes() == content
    assert session.get.call_count == 2


def test_download_file_removes_rejected_payload(tmp_path, mocker):
    """No file should remain when all downloaded payloads fail validation."""
    session = mocker.patch("reacnetgenerator.utils.requests.Session").return_value
    # Include non-text bytes to ensure checksum rejection never tries to decode
    # arbitrary downloaded content for logging.
    session.get.return_value = _DownloadResponse(b"\xffcorrupt trajectory")
    output = tmp_path / "trajectory.dump"

    with pytest.raises(RuntimeError, match="Cannot download"):
        asyncio.run(download_file("https://corrupt", str(output), _sha256(b"valid")))

    assert not output.exists()
