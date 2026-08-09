"""Real-data evaluation utilities for CASPR reconstruction."""

from __future__ import annotations

import csv
import json
import statistics
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .metrics import evaluate_sequence
from .models import ClusterRecord, validate_dna_sequence
from .reconstruct import reconstruct_cluster


class TruthFormatError(ValueError):
    """Raised when a truth TSV record is malformed."""


class DuplicateTruthClusterError(ValueError):
    """Raised when a truth TSV repeats a cluster ID."""


@dataclass(frozen=True)
class ClusterBenchmarkRecord:
    """Evaluation result for one reconstructed cluster."""

    cluster_id: str
    truth_sequence: str
    core_sequence: str
    consensus: str

    raw_read_count: int
    unique_sequence_count: int
    pairwise_alignment_count: int

    core_exact: bool
    core_edit_distance: int
    core_normalized_edit_distance: float

    consensus_exact: bool
    consensus_edit_distance: int
    consensus_normalized_edit_distance: float

    consensus_substitutions: int
    consensus_insertions: int
    consensus_deletions: int

    length_error: int
    edit_distance_improvement: int

    reconstruction_backend: str
    metric_backend: str
    normalize_indels: bool

    runtime_seconds: float


@dataclass(frozen=True)
class BenchmarkSummary:
    """Aggregate statistics across benchmarked clusters."""

    cluster_count: int

    core_exact_count: int
    core_exact_rate: float

    consensus_exact_count: int
    consensus_exact_rate: float

    improved_count: int
    unchanged_count: int
    worsened_count: int

    mean_core_edit_distance: float
    mean_consensus_edit_distance: float

    median_core_edit_distance: float
    median_consensus_edit_distance: float

    mean_core_normalized_edit_distance: float
    mean_consensus_normalized_edit_distance: float

    median_core_normalized_edit_distance: float
    median_consensus_normalized_edit_distance: float

    total_runtime_seconds: float
    mean_runtime_seconds: float

    total_pairwise_alignment_count: int

    mean_raw_read_count: float
    mean_unique_sequence_count: float


BENCHMARK_RECORD_COLUMNS = (
    "cluster_id",
    "truth_sequence",
    "core_sequence",
    "consensus",
    "raw_read_count",
    "unique_sequence_count",
    "pairwise_alignment_count",
    "core_exact",
    "core_edit_distance",
    "core_normalized_edit_distance",
    "consensus_exact",
    "consensus_edit_distance",
    "consensus_normalized_edit_distance",
    "consensus_substitutions",
    "consensus_insertions",
    "consensus_deletions",
    "length_error",
    "edit_distance_improvement",
    "reconstruction_backend",
    "metric_backend",
    "normalize_indels",
    "runtime_seconds",
)


def read_truth_tsv(
    path: str | Path,
) -> dict[str, str]:
    """
    Load headerless cluster_id<TAB>truth_sequence records.

    Empty lines are ignored. Cluster IDs must be unique and non-empty.
    Truth sequences use the established A/C/G/T/N alphabet.
    """
    truth_by_cluster_id: dict[str, str] = {}

    with Path(path).open(
        "r",
        encoding="utf-8",
    ) as source:
        for line_number, raw_line in enumerate(
            source,
            start=1,
        ):
            line = raw_line.rstrip("\r\n")

            if not line.strip():
                continue

            fields = line.split("\t")

            if len(fields) != 2:
                raise TruthFormatError(
                    f"truth TSV line {line_number} must contain "
                    "exactly two tab-separated fields"
                )

            cluster_id, truth_sequence = fields

            if not cluster_id or not cluster_id.strip():
                raise TruthFormatError(
                    f"truth TSV line {line_number} has an "
                    "empty cluster ID"
                )

            if not truth_sequence:
                raise TruthFormatError(
                    f"truth TSV line {line_number} has an "
                    "empty truth sequence"
                )

            try:
                validate_dna_sequence(
                    truth_sequence,
                    "truth_sequence",
                )
            except ValueError as exc:
                raise TruthFormatError(
                    f"truth TSV line {line_number}: {exc}"
                ) from exc

            if cluster_id in truth_by_cluster_id:
                raise DuplicateTruthClusterError(
                    "duplicate truth cluster ID on line "
                    f"{line_number}: {cluster_id}"
                )

            truth_by_cluster_id[cluster_id] = (
                truth_sequence
            )

    return truth_by_cluster_id


def _validate_max_clusters(
    max_clusters: int | None,
) -> None:
    if max_clusters is None:
        return

    if (
        isinstance(max_clusters, bool)
        or not isinstance(max_clusters, int)
        or max_clusters <= 0
    ):
        raise ValueError(
            "max_clusters must be None or a positive integer"
        )


def benchmark_clusters(
    clusters: Iterable[ClusterRecord],
    truth_by_cluster_id: Mapping[str, str],
    *,
    reconstruction_backend: str = "edlib",
    metric_backend: str = "edlib",
    normalize_indels: bool = False,
    max_clusters: int | None = None,
) -> list[ClusterBenchmarkRecord]:
    """
    Reconstruct and evaluate clusters in supplied iterable order.

    runtime_seconds measures reconstruct_cluster only. Ground-truth
    metric evaluation begins after the reconstruction timer stops.
    """
    _validate_max_clusters(max_clusters)

    records: list[ClusterBenchmarkRecord] = []

    for cluster in clusters:
        if (
            max_clusters is not None
            and len(records) >= max_clusters
        ):
            break

        if cluster.cluster_id not in truth_by_cluster_id:
            raise ValueError(
                "missing ground truth for processed cluster: "
                f"{cluster.cluster_id}"
            )

        truth_sequence = truth_by_cluster_id[
            cluster.cluster_id
        ]

        started_at = time.perf_counter()

        reconstruction = reconstruct_cluster(
            cluster,
            backend=reconstruction_backend,
            normalize_indels=normalize_indels,
        )

        runtime_seconds = (
            time.perf_counter() - started_at
        )

        core_metrics = evaluate_sequence(
            truth_sequence,
            cluster.core_sequence,
            backend=metric_backend,
        )

        consensus_metrics = evaluate_sequence(
            truth_sequence,
            reconstruction.consensus,
            backend=metric_backend,
        )

        records.append(
            ClusterBenchmarkRecord(
                cluster_id=cluster.cluster_id,
                truth_sequence=truth_sequence,
                core_sequence=cluster.core_sequence,
                consensus=reconstruction.consensus,
                raw_read_count=(
                    reconstruction.raw_read_count
                ),
                unique_sequence_count=(
                    reconstruction.unique_sequence_count
                ),
                pairwise_alignment_count=(
                    reconstruction.pairwise_alignment_count
                ),
                core_exact=core_metrics.exact_match,
                core_edit_distance=(
                    core_metrics.edit_distance
                ),
                core_normalized_edit_distance=(
                    core_metrics.normalized_edit_distance
                ),
                consensus_exact=(
                    consensus_metrics.exact_match
                ),
                consensus_edit_distance=(
                    consensus_metrics.edit_distance
                ),
                consensus_normalized_edit_distance=(
                    consensus_metrics.normalized_edit_distance
                ),
                consensus_substitutions=(
                    consensus_metrics.substitutions
                ),
                consensus_insertions=(
                    consensus_metrics.insertions
                ),
                consensus_deletions=(
                    consensus_metrics.deletions
                ),
                length_error=(
                    consensus_metrics.length_error
                ),
                edit_distance_improvement=(
                    core_metrics.edit_distance
                    - consensus_metrics.edit_distance
                ),
                reconstruction_backend=(
                    reconstruction_backend
                ),
                metric_backend=metric_backend,
                normalize_indels=normalize_indels,
                runtime_seconds=runtime_seconds,
            )
        )

    return records


def summarize_benchmark(
    records: Iterable[ClusterBenchmarkRecord],
) -> BenchmarkSummary:
    """Build deterministic aggregate benchmark statistics."""
    record_list = list(records)

    if not record_list:
        raise ValueError(
            "cannot summarize an empty benchmark"
        )

    cluster_count = len(record_list)

    core_exact_count = sum(
        record.core_exact
        for record in record_list
    )

    consensus_exact_count = sum(
        record.consensus_exact
        for record in record_list
    )

    improvements = [
        record.edit_distance_improvement
        for record in record_list
    ]

    core_edit_distances = [
        record.core_edit_distance
        for record in record_list
    ]

    consensus_edit_distances = [
        record.consensus_edit_distance
        for record in record_list
    ]

    core_normalized_distances = [
        record.core_normalized_edit_distance
        for record in record_list
    ]

    consensus_normalized_distances = [
        record.consensus_normalized_edit_distance
        for record in record_list
    ]

    runtimes = [
        record.runtime_seconds
        for record in record_list
    ]

    return BenchmarkSummary(
        cluster_count=cluster_count,
        core_exact_count=core_exact_count,
        core_exact_rate=(
            core_exact_count / cluster_count
        ),
        consensus_exact_count=(
            consensus_exact_count
        ),
        consensus_exact_rate=(
            consensus_exact_count / cluster_count
        ),
        improved_count=sum(
            improvement > 0
            for improvement in improvements
        ),
        unchanged_count=sum(
            improvement == 0
            for improvement in improvements
        ),
        worsened_count=sum(
            improvement < 0
            for improvement in improvements
        ),
        mean_core_edit_distance=float(
            statistics.mean(
                core_edit_distances
            )
        ),
        mean_consensus_edit_distance=float(
            statistics.mean(
                consensus_edit_distances
            )
        ),
        median_core_edit_distance=float(
            statistics.median(
                core_edit_distances
            )
        ),
        median_consensus_edit_distance=float(
            statistics.median(
                consensus_edit_distances
            )
        ),
        mean_core_normalized_edit_distance=float(
            statistics.mean(
                core_normalized_distances
            )
        ),
        mean_consensus_normalized_edit_distance=float(
            statistics.mean(
                consensus_normalized_distances
            )
        ),
        median_core_normalized_edit_distance=float(
            statistics.median(
                core_normalized_distances
            )
        ),
        median_consensus_normalized_edit_distance=float(
            statistics.median(
                consensus_normalized_distances
            )
        ),
        total_runtime_seconds=sum(runtimes),
        mean_runtime_seconds=float(
            statistics.mean(runtimes)
        ),
        total_pairwise_alignment_count=sum(
            record.pairwise_alignment_count
            for record in record_list
        ),
        mean_raw_read_count=float(
            statistics.mean(
                record.raw_read_count
                for record in record_list
            )
        ),
        mean_unique_sequence_count=float(
            statistics.mean(
                record.unique_sequence_count
                for record in record_list
            )
        ),
    )


def _record_to_csv_row(
    record: ClusterBenchmarkRecord,
) -> dict[str, object]:
    row = asdict(record)

    row["core_exact"] = (
        "true" if record.core_exact else "false"
    )

    row["consensus_exact"] = (
        "true" if record.consensus_exact else "false"
    )

    row["normalize_indels"] = (
        "true" if record.normalize_indels else "false"
    )

    return {
        column: row[column]
        for column in BENCHMARK_RECORD_COLUMNS
    }


def write_benchmark_records_csv(
    path: str | Path,
    records: Iterable[ClusterBenchmarkRecord],
) -> None:
    """Write records using a stable explicit column order."""
    with Path(path).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=BENCHMARK_RECORD_COLUMNS,
            extrasaction="raise",
        )

        writer.writeheader()

        for record in records:
            writer.writerow(
                _record_to_csv_row(record)
            )


def write_benchmark_summary_json(
    path: str | Path,
    summary: BenchmarkSummary,
) -> None:
    """Write a benchmark summary as UTF-8 JSON."""
    with Path(path).open(
        "w",
        encoding="utf-8",
    ) as destination:
        json.dump(
            asdict(summary),
            destination,
            indent=2,
            ensure_ascii=True,
        )

        destination.write("\n")
