"""Adapter from Clover membership output to ClusterRecord objects.

This module intentionally does not infer Clover core sequences from member
reads. Core sequences must come from a reliable explicit source.
"""

from __future__ import annotations

import ast
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .clover_io import read_clover_input
from .models import ClusterRecord, ReadRecord


class MalformedMembershipError(ValueError):
    """Raised when a Clover membership output file has invalid syntax."""


class DuplicateMembershipError(ValueError):
    """Raised when a read appears more than once in membership output."""


class ConflictingMembershipError(ValueError):
    """Raised when a read is assigned to more than one cluster/core."""


class CoreSequenceUnavailable(ValueError):
    """Raised when ClusterRecord construction lacks reliable core sequences."""


@dataclass(frozen=True)
class AdapterStats:
    total_input_reads: int
    assigned_reads: int
    unassigned_reads: int
    cluster_count: int


@dataclass(frozen=True)
class ClusterBuildResult:
    clusters: list[ClusterRecord]
    stats: AdapterStats


def _parse_membership_literal(text: str, path: Path):
    try:
        parsed = ast.literal_eval(text.strip())
    except (SyntaxError, ValueError) as exc:
        raise MalformedMembershipError(f"{path}: malformed membership literal") from exc
    if not isinstance(parsed, list):
        raise MalformedMembershipError(f"{path}: membership output must be a list")
    return parsed


def read_membership(path: str | Path) -> dict[str, str]:
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        return {}
    parsed = _parse_membership_literal(text, path)
    membership: dict[str, str] = {}
    for index, item in enumerate(parsed):
        if not isinstance(item, tuple) or len(item) != 2:
            raise MalformedMembershipError(
                f"{path}: membership item {index} must be a 2-tuple"
            )
        read_id, cluster_id = str(item[0]), str(item[1])
        if not read_id:
            raise MalformedMembershipError(f"{path}: membership item {index} has empty read_id")
        if not cluster_id:
            raise MalformedMembershipError(f"{path}: membership item {index} has empty cluster_id")
        if read_id in membership:
            if membership[read_id] != cluster_id:
                raise ConflictingMembershipError(
                    f"read_id {read_id} assigned to both {membership[read_id]} and {cluster_id}"
                )
            raise DuplicateMembershipError(f"duplicate membership for read_id {read_id}")
        membership[read_id] = cluster_id
    return membership


def read_core_sequences(path: str | Path) -> dict[str, str]:
    path = Path(path)
    cores: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected cluster_id and core_sequence")
            cluster_id, core_sequence = fields
            if cluster_id in cores:
                raise ValueError(f"{path}:{line_number}: duplicate core id {cluster_id}")
            cores[cluster_id] = core_sequence
    return cores


def build_cluster_records(
    input_path: str | Path,
    membership_path: str | Path,
    core_sequences: Mapping[str, str] | None = None,
) -> ClusterBuildResult:
    reads_by_id = read_clover_input(input_path)
    membership = read_membership(membership_path)

    grouped: dict[str, list[ReadRecord]] = defaultdict(list)
    for read_id, cluster_id in membership.items():
        if read_id not in reads_by_id:
            raise ValueError(f"membership references unknown read_id: {read_id}")
        grouped[cluster_id].append(reads_by_id[read_id])

    missing_core_ids = sorted(cid for cid in grouped if core_sequences is None or cid not in core_sequences)
    if missing_core_ids:
        raise CoreSequenceUnavailable(
            "CORE_SEQUENCE_UNAVAILABLE: missing reliable core_sequence for cluster(s): "
            + ", ".join(missing_core_ids)
        )

    clusters = [
        ClusterRecord(cluster_id=cluster_id, core_sequence=core_sequences[cluster_id], reads=members)
        for cluster_id, members in sorted(grouped.items())
    ]
    assigned = len(membership)
    stats = AdapterStats(
        total_input_reads=len(reads_by_id),
        assigned_reads=assigned,
        unassigned_reads=len(reads_by_id) - assigned,
        cluster_count=len(clusters),
    )
    return ClusterBuildResult(clusters=clusters, stats=stats)
