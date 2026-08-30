# SPDX-License-Identifier: LGPL-3.0-or-later
"""Test utility helpers."""

import asyncio
import hashlib
from io import BytesIO

import pytest
import requests
from urllib3.exceptions import ProtocolError

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

    def iter_content(self, chunk_size):
        """Yield the response body using the Requests streaming interface."""
        while chunk := self.raw.read(chunk_size):
            yield chunk


class _InterruptedRaw:
    """Simulate an urllib3 stream that fails after yielding partial data."""

    def stream(self, chunk_size, decode_content=True):
        yield b"partial"
        raise ProtocolError("connection interrupted")

    def close(self):
        """Support cleanup by requests.Response's context manager."""


def _interrupted_response():
    """Return a real Requests response backed by an interrupted raw stream."""
    response = requests.Response()
    response.status_code = 200
    response.raw = _InterruptedRaw()
    return response


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


def test_download_file_retries_after_stream_interruption(tmp_path, mocker):
    """A partial streaming response should be removed before mirror fallback."""
    content = b"valid trajectory"
    session = mocker.patch("reacnetgenerator.utils.requests.Session").return_value
    session.get.side_effect = [
        _interrupted_response(),
        _DownloadResponse(content),
    ]
    output = tmp_path / "trajectory.dump"

    asyncio.run(
        download_file(
            ["https://interrupted", "https://good"],
            str(output),
            _sha256(content),
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
