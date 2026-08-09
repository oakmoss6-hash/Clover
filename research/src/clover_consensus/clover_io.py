"""Readers for original Clover input files.

Clover txt inputs are whitespace-delimited records where column 1 is a stable
read identifier/tag/index and column 2 is the DNA sequence. The original parser
uses ``line.split()`` and reads exactly these first two fields.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from .models import ReadRecord


class DuplicateReadIdError(ValueError):
    """Raised when a Clover input file contains the same read id twice."""


class MalformedInputError(ValueError):
    """Raised when a Clover input row cannot be parsed."""


def iter_clover_input(path: str | Path) -> Iterator[ReadRecord]:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 2:
                raise MalformedInputError(
                    f"{path}:{line_number}: expected at least 2 whitespace-delimited columns"
                )
            read_id, sequence = fields[0], fields[1]
            if not read_id:
                raise MalformedInputError(f"{path}:{line_number}: missing read_id")
            yield ReadRecord(read_id=read_id, sequence=sequence, tag=read_id)


def read_clover_input(path: str | Path) -> dict[str, ReadRecord]:
    records: dict[str, ReadRecord] = {}
    for record in iter_clover_input(path):
        if record.read_id in records:
            raise DuplicateReadIdError(f"duplicate read_id in Clover input: {record.read_id}")
        records[record.read_id] = record
    return records
