"""Generic real-data benchmark CLI for CASPR."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


RESEARCH_ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(
    0,
    str(RESEARCH_ROOT / "src"),
)

from clover_consensus.benchmark import (
    benchmark_clusters,
    read_truth_tsv,
    summarize_benchmark,
    write_benchmark_records_csv,
    write_benchmark_summary_json,
)
from clover_consensus.cluster_adapter import (
    build_cluster_records,
    read_core_sequences,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate existing CASPR reconstruction against "
            "generic cluster truth TSV data."
        )
    )

    parser.add_argument(
        "--input",
        dest="input_path",
        required=True,
        help="Clover read/input data path",
    )

    parser.add_argument(
        "--membership",
        required=True,
        help="Clover cluster membership output path",
    )

    parser.add_argument(
        "--core-export",
        required=True,
        help="Passive initial-core export TSV path",
    )

    parser.add_argument(
        "--truth",
        required=True,
        help="Headerless cluster_id<TAB>truth_sequence TSV",
    )

    parser.add_argument(
        "--output-records",
        required=True,
        help="Output per-cluster CSV path",
    )

    parser.add_argument(
        "--output-summary",
        required=True,
        help="Output summary JSON path",
    )

    parser.add_argument(
        "--backend",
        default="edlib",
        help="CASPR reconstruction alignment backend",
    )

    parser.add_argument(
        "--metric-backend",
        default="edlib",
        help="Truth evaluation alignment backend",
    )

    parser.add_argument(
        "--normalize-indels",
        action="store_true",
        help="Enable deterministic indel normalization",
    )

    parser.add_argument(
        "--max-clusters",
        type=int,
        default=None,
        help="Optional positive cluster processing limit",
    )

    return parser


def main(
    argv: list[str] | None = None,
) -> int:
    args = build_parser().parse_args(argv)

    core_sequences = read_core_sequences(
        args.core_export
    )

    cluster_result = build_cluster_records(
        args.input_path,
        args.membership,
        core_sequences,
    )

    truth_by_cluster_id = read_truth_tsv(
        args.truth
    )

    records = benchmark_clusters(
        cluster_result.clusters,
        truth_by_cluster_id,
        reconstruction_backend=args.backend,
        metric_backend=args.metric_backend,
        normalize_indels=args.normalize_indels,
        max_clusters=args.max_clusters,
    )

    summary = summarize_benchmark(records)

    Path(args.output_records).parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    Path(args.output_summary).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    write_benchmark_records_csv(
        args.output_records,
        records,
    )

    write_benchmark_summary_json(
        args.output_summary,
        summary,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
