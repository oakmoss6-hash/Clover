#!/usr/bin/env python3
"""Phase 14.3B: truth-blind Clover clustering/reconstruction accuracy audit."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path


BENCHMARK_DIR = Path(__file__).resolve().parent
RESEARCH_ROOT = BENCHMARK_DIR.parent
ROOT = RESEARCH_ROOT.parent

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(RESEARCH_ROOT / "src"))
sys.path.insert(0, str(BENCHMARK_DIR))

import run_phase14_clover_e2e as e2e

from clover_consensus.clover_worker_reconstruction import (
    CloverWorkerReconstructor,
)
from clover_consensus.metrics import evaluate_sequence


@dataclass
class EvaluationWorkerReconstructor(CloverWorkerReconstructor):
    """
    Record read IDs only after Clover has already decided membership.

    No truth labels or references are available inside the worker.
    """

    read_memberships: dict[int, list[str]] = field(
        default_factory=dict
    )

    def record_membership(
        self,
        *,
        core_index: int,
        sequence: str,
        read_id: str | None = None,
        is_new_core: bool = False,
    ) -> None:
        if read_id is None:
            raise RuntimeError(
                "Phase 14.3 requires read_id metadata"
            )

        super().record_membership(
            core_index=core_index,
            sequence=sequence,
            read_id=read_id,
            is_new_core=is_new_core,
        )

        self.read_memberships.setdefault(
            core_index,
            [],
        ).append(read_id)

    def finalize(self) -> dict:
        for core_index, state in self.adapter.states.items():
            ids = self.read_memberships.get(core_index, [])

            if len(ids) != state.raw_read_count:
                raise RuntimeError(
                    "membership count mismatch: "
                    f"worker={self.worker_name!r} "
                    f"cluster={core_index} "
                    f"ids={len(ids)} "
                    f"M={state.raw_read_count}"
                )

        output = super().finalize()
        prefix = self.worker_name

        output[prefix + "phase14_read_memberships"] = [
            (
                core_index,
                tuple(self.read_memberships[core_index]),
            )
            for core_index in sorted(self.adapter.states)
        ]

        output[prefix + "phase14_routing_cores"] = [
            (
                core_index,
                self.adapter.states[core_index].routing_core,
            )
            for core_index in sorted(self.adapter.states)
        ]

        return output


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--references", required=True)

    parser.add_argument(
        "--max-reads",
        type=int,
        default=100000,
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        choices=[1, 4, 16, 64],
    )

    parser.add_argument(
        "--backend",
        default="wfa",
        choices=["wfa", "edlib", "nw"],
    )

    parser.add_argument(
        "--backbone-policy",
        default="core",
        choices=[
            "core",
            "support_length",
            "max_span",
        ],
    )

    parser.add_argument(
        "--metric-backend",
        default="edlib",
    )

    parser.add_argument(
        "--expected-reference-count",
        type=int,
        default=72000,
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
    )

    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    parts: list[str] = []

    with path.open(
        "r",
        encoding="utf-8",
    ) as source:
        for raw in source:
            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(parts)

                header = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.upper())

    if header is not None:
        records[header] = "".join(parts)

    return records


def load_truth_for_read_ids(
    audit_path: Path,
    wanted_ids: set[str],
    input_sequences: dict[str, str],
):
    """
    Scan the complete unordered audit file.

    Truth is loaded only after Clover/reconstruction have finished.
    """
    truth: dict[str, str] = {}

    scanned = 0
    duplicates = 0
    conflicts = 0
    sequence_mismatches = 0

    with audit_path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as source:
        header = source.readline().rstrip("\n").split("\t")

        if header[:3] != [
            "read_id",
            "label",
            "sequence",
        ]:
            raise RuntimeError(
                f"unexpected audit header: {header!r}"
            )

        for raw in source:
            scanned += 1
            fields = raw.rstrip("\n").split("\t")

            if len(fields) != 9:
                raise RuntimeError(
                    f"audit row {scanned + 1}: "
                    f"expected 9 fields, got {len(fields)}"
                )

            read_id = fields[0]

            if read_id not in wanted_ids:
                continue

            label = fields[1]
            audit_sequence = fields[2].upper()

            if input_sequences[read_id] != audit_sequence:
                sequence_mismatches += 1

            previous = truth.get(read_id)

            if previous is None:
                truth[read_id] = label
            else:
                duplicates += 1

                if previous != label:
                    conflicts += 1

    if sequence_mismatches:
        raise RuntimeError(
            "FASTQ/audit sequence mismatch count: "
            f"{sequence_mismatches}"
        )

    if conflicts:
        raise RuntimeError(
            f"conflicting truth labels: {conflicts}"
        )

    return truth, {
        "audit_rows_scanned": scanned,
        "matched_assigned_reads": len(truth),
        "duplicate_rows": duplicates,
        "conflicting_labels": conflicts,
        "sequence_mismatches": sequence_mismatches,
    }


def sequence_metric(
    truth: str,
    sequence: str,
    backend: str,
):
    metric = evaluate_sequence(
        truth,
        sequence,
        backend=backend,
    )

    return {
        "exact": int(metric.exact_match),
        "edit_distance": metric.edit_distance,
        "substitutions": metric.substitutions,
        "insertions": metric.insertions,
        "deletions": metric.deletions,
    }


def summarize(rows: list[dict]) -> dict:
    if not rows:
        return {
            "cluster_count": 0,
        }

    n = len(rows)

    core_exact = sum(
        row["core_exact"]
        for row in rows
    )

    consensus_exact = sum(
        row["consensus_exact"]
        for row in rows
    )

    core_mean_ed = statistics.fmean(
        row["core_edit_distance"]
        for row in rows
    )

    consensus_mean_ed = statistics.fmean(
        row["consensus_edit_distance"]
        for row in rows
    )

    better = sum(
        row["comparison"] == "better"
        for row in rows
    )

    same = sum(
        row["comparison"] == "same"
        for row in rows
    )

    worse = sum(
        row["comparison"] == "worse"
        for row in rows
    )

    labeled = sum(
        row["labeled_read_count"]
        for row in rows
    )

    majority = sum(
        row["majority_count"]
        for row in rows
    )

    raw = sum(
        row["raw_read_count"]
        for row in rows
    )

    return {
        "cluster_count": n,
        "raw_read_count": raw,
        "labeled_read_count": labeled,
        "label_coverage": (
            labeled / raw
            if raw else 0.0
        ),
        "mean_purity": statistics.fmean(
            row["purity"]
            for row in rows
        ),
        "weighted_purity": (
            majority / labeled
            if labeled else 0.0
        ),
        "core_exact_count": core_exact,
        "core_exact_rate": core_exact / n,
        "consensus_exact_count": consensus_exact,
        "consensus_exact_rate": (
            consensus_exact / n
        ),
        "exact_rate_gain": (
            consensus_exact - core_exact
        ) / n,
        "core_mean_edit_distance": core_mean_ed,
        "consensus_mean_edit_distance": (
            consensus_mean_ed
        ),
        "mean_edit_distance_reduction": (
            core_mean_ed - consensus_mean_ed
        ),
        "better": better,
        "same": same,
        "worse": worse,
    }


def main():
    args = parse_args()

    # ==========================================================
    # Stage 1: original reads only. No truth loaded.
    # ==========================================================
    records = []

    for read_id, sequence in e2e.iter_fastq(
        args.input
    ):
        if (
            args.max_reads is not None
            and len(records) >= args.max_reads
        ):
            break

        records.append(
            (
                read_id,
                sequence,
            )
        )

    input_sequences = dict(records)

    names, shards, invalid_prefix_reads = (
        e2e.partition_records(
            records,
            args.workers,
        )
    )

    print(
        f"loaded {len(records):,} original FASTQ reads"
    )
    print(
        "truth labels/references are NOT loaded yet"
    )

    # ==========================================================
    # Stage 2: Clover + reconstruction, truth-blind.
    # ==========================================================
    MyProcess, _ = e2e.import_components()

    outputs, algorithm_wall = e2e.run_pass(
        MyProcess=MyProcess,
        CloverWorkerReconstructor=(
            EvaluationWorkerReconstructor
        ),
        names=names,
        shards=shards,
        reconstruct=True,
        backend=args.backend,
        backbone_policy=args.backbone_policy,
    )

    reconstruction = e2e.collect_reconstruction(
        outputs,
        names,
    )

    cluster_rows = {}

    for row in reconstruction["rows"]:
        key = (
            row["worker"],
            int(row["cluster_id"]),
        )

        cluster_rows[key] = dict(row)

    memberships = {}
    routing_cores = {}

    for name in names:
        output = outputs[name]

        for core_index, read_ids in output[
            name + "phase14_read_memberships"
        ]:
            memberships[
                (
                    name,
                    int(core_index),
                )
            ] = list(read_ids)

        for core_index, core in output[
            name + "phase14_routing_cores"
        ]:
            routing_cores[
                (
                    name,
                    int(core_index),
                )
            ] = core

    if set(cluster_rows) != set(memberships):
        raise RuntimeError(
            "cluster/result membership keys differ"
        )

    if set(cluster_rows) != set(routing_cores):
        raise RuntimeError(
            "cluster/routing-core keys differ"
        )

    assigned_ids = set()
    membership_total = 0

    for key, ids in memberships.items():
        row = cluster_rows[key]

        if len(ids) != row["raw_read_count"]:
            raise RuntimeError(
                f"{key}: read-id count "
                f"{len(ids)} != "
                f"M={row['raw_read_count']}"
            )

        membership_total += len(ids)
        assigned_ids.update(ids)

    if len(assigned_ids) != membership_total:
        raise RuntimeError(
            "same read_id assigned to multiple clusters"
        )

    prefix = Path(args.output_prefix)

    prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    membership_path = Path(
        str(prefix) + ".memberships.tsv"
    )

    # Materialize truth-blind memberships before loading labels.
    with membership_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "worker",
                "cluster_id",
                "read_id",
            ]
        )

        for key in sorted(memberships):
            worker, cluster_id = key

            for read_id in memberships[key]:
                writer.writerow(
                    [
                        worker,
                        cluster_id,
                        read_id,
                    ]
                )

    print()
    print("Clover + reconstruction finished")
    print("truth still not loaded")
    print(
        "algorithm wall:",
        f"{algorithm_wall:.4f}s",
    )
    print(
        "clusters C:",
        f"{reconstruction['cluster_count']:,}",
    )
    print(
        "assigned M:",
        f"{reconstruction['raw_clustered_reads']:,}",
    )
    print(
        "unique U:",
        f"{reconstruction['unique_sequences']:,}",
    )
    print(
        "pairwise A:",
        f"{reconstruction['pairwise_alignments']:,}",
    )

    # ==========================================================
    # Stage 3: post-hoc truth loading.
    # ==========================================================
    print()
    print(
        "Now loading truth labels and references..."
    )

    truth_started = time.perf_counter()

    read_truth, truth_join_stats = (
        load_truth_for_read_ids(
            Path(args.audit),
            assigned_ids,
            input_sequences,
        )
    )

    references = read_fasta(
        Path(args.references)
    )

    truth_join_seconds = (
        time.perf_counter()
        - truth_started
    )

    if len(references) != args.expected_reference_count:
        raise RuntimeError(
            f"reference count "
            f"{len(references)} != "
            f"{args.expected_reference_count}"
        )

    # ==========================================================
    # Stage 4: cluster purity + core/consensus accuracy.
    # ==========================================================
    evaluation_started = time.perf_counter()

    rows = []

    for key in sorted(cluster_rows):
        worker, cluster_id = key

        reconstruction_row = cluster_rows[key]
        ids = memberships[key]

        labels = [
            read_truth[read_id]
            for read_id in ids
            if read_id in read_truth
        ]

        counts = Counter(labels)

        raw_count = reconstruction_row[
            "raw_read_count"
        ]

        labeled_count = len(labels)

        base_row = {
            "worker": worker,
            "cluster_id": cluster_id,
            "raw_read_count": raw_count,
            "unique_sequence_count": (
                reconstruction_row[
                    "unique_sequence_count"
                ]
            ),
            "pairwise_alignment_count": (
                reconstruction_row[
                    "pairwise_alignment_count"
                ]
            ),
            "labeled_read_count": (
                labeled_count
            ),
            "label_coverage": (
                labeled_count / raw_count
                if raw_count else 0.0
            ),
            "fully_labeled": int(
                labeled_count == raw_count
            ),
            "distinct_truth_labels": (
                len(counts)
            ),
            "majority_label": "",
            "majority_count": 0,
            "purity": "",
            "evaluation_status": "",
            "routing_core": (
                routing_cores[key]
            ),
            "backbone": (
                reconstruction_row[
                    "backbone"
                ]
            ),
            "consensus": (
                reconstruction_row[
                    "consensus"
                ]
            ),
            "truth": "",
            "core_exact": "",
            "core_edit_distance": "",
            "consensus_exact": "",
            "consensus_edit_distance": "",
            "comparison": "",
        }

        if not counts:
            base_row[
                "evaluation_status"
            ] = "no_labeled_reads"
            rows.append(base_row)
            continue

        majority_count = max(
            counts.values()
        )

        majority_labels = sorted(
            label
            for label, count
            in counts.items()
            if count == majority_count
        )

        purity = (
            majority_count / labeled_count
        )

        base_row[
            "majority_count"
        ] = majority_count

        base_row[
            "purity"
        ] = purity

        if len(majority_labels) != 1:
            base_row[
                "evaluation_status"
            ] = "majority_tie"

            base_row[
                "majority_label"
            ] = "|".join(
                majority_labels
            )

            rows.append(base_row)
            continue

        truth_label = majority_labels[0]

        truth_sequence = references.get(
            truth_label
        )

        base_row[
            "majority_label"
        ] = truth_label

        if truth_sequence is None:
            base_row[
                "evaluation_status"
            ] = "missing_reference"

            rows.append(base_row)
            continue

        core_metric = sequence_metric(
            truth_sequence,
            routing_cores[key],
            args.metric_backend,
        )

        consensus_metric = sequence_metric(
            truth_sequence,
            reconstruction_row[
                "consensus"
            ],
            args.metric_backend,
        )

        if (
            consensus_metric["edit_distance"]
            < core_metric["edit_distance"]
        ):
            comparison = "better"
        elif (
            consensus_metric["edit_distance"]
            > core_metric["edit_distance"]
        ):
            comparison = "worse"
        else:
            comparison = "same"

        base_row.update(
            {
                "evaluation_status": "ok",
                "truth": truth_sequence,
                "core_exact": (
                    core_metric["exact"]
                ),
                "core_edit_distance": (
                    core_metric[
                        "edit_distance"
                    ]
                ),
                "consensus_exact": (
                    consensus_metric[
                        "exact"
                    ]
                ),
                "consensus_edit_distance": (
                    consensus_metric[
                        "edit_distance"
                    ]
                ),
                "comparison": comparison,
            }
        )

        rows.append(base_row)

    evaluation_seconds = (
        time.perf_counter()
        - evaluation_started
    )

    ok_rows = [
        row
        for row in rows
        if row["evaluation_status"] == "ok"
    ]

    groups = {
        "all_evaluable": ok_rows,

        "multi_unique": [
            row
            for row in ok_rows
            if row[
                "unique_sequence_count"
            ] > 1
        ],

        "pure": [
            row
            for row in ok_rows
            if row["purity"] == 1.0
        ],

        "pure_multi_unique": [
            row
            for row in ok_rows
            if (
                row["purity"] == 1.0
                and row[
                    "unique_sequence_count"
                ] > 1
            )
        ],

        "mixed": [
            row
            for row in ok_rows
            if row["purity"] < 1.0
        ],

        "fully_labeled_pure_multi_unique": [
            row
            for row in ok_rows
            if (
                row["purity"] == 1.0
                and row["fully_labeled"] == 1
                and row[
                    "unique_sequence_count"
                ] > 1
            )
        ],
    }

    summaries = {
        name: summarize(group)
        for name, group in groups.items()
    }

    labeled_rows = [
        row
        for row in rows
        if row["labeled_read_count"] > 0
    ]

    total_labeled = sum(
        row["labeled_read_count"]
        for row in labeled_rows
    )

    total_majority = sum(
        row["majority_count"]
        for row in labeled_rows
    )

    clustering = {
        "clusters_with_labeled_reads": (
            len(labeled_rows)
        ),
        "assigned_labeled_reads": (
            total_labeled
        ),
        "assigned_label_coverage": (
            total_labeled
            / reconstruction[
                "raw_clustered_reads"
            ]
            if reconstruction[
                "raw_clustered_reads"
            ]
            else 0.0
        ),
        "weighted_purity": (
            total_majority
            / total_labeled
            if total_labeled
            else 0.0
        ),
        "majority_tie_clusters": sum(
            row["evaluation_status"]
            == "majority_tie"
            for row in rows
        ),
    }

    summary = {
        "phase": "14.3B",
        "processed_reads": len(records),
        "workers": args.workers,
        "backend": args.backend,
        "backbone_policy": (
            args.backbone_policy
        ),
        "metric_backend": (
            args.metric_backend
        ),
        "truth_usage": (
            "truth loaded only after Clover "
            "clustering and reconstruction completed"
        ),
        "algorithm_wall_seconds": (
            algorithm_wall
        ),
        "truth_join_seconds": (
            truth_join_seconds
        ),
        "evaluation_seconds": (
            evaluation_seconds
        ),
        "invalid_prefix_reads": (
            invalid_prefix_reads
        ),
        "cluster_count_C": (
            reconstruction[
                "cluster_count"
            ]
        ),
        "clustered_raw_reads_M": (
            reconstruction[
                "raw_clustered_reads"
            ]
        ),
        "unique_sequences_U": (
            reconstruction[
                "unique_sequences"
            ]
        ),
        "pairwise_alignments_A": (
            reconstruction[
                "pairwise_alignments"
            ]
        ),
        "expected_U_minus_C": (
            reconstruction[
                "unique_sequences"
            ]
            - reconstruction[
                "cluster_count"
            ]
        ),
        "pairwise_invariant_ok": (
            reconstruction[
                "pairwise_alignments"
            ]
            == reconstruction[
                "unique_sequences"
            ]
            - reconstruction[
                "cluster_count"
            ]
        ),
        "truth_join": truth_join_stats,
        "clustering": clustering,
        "groups": summaries,
    }

    cluster_path = Path(
        str(prefix) + ".clusters.tsv"
    )

    summary_path = Path(
        str(prefix) + ".summary.json"
    )

    fields = [
        "worker",
        "cluster_id",
        "raw_read_count",
        "unique_sequence_count",
        "pairwise_alignment_count",
        "labeled_read_count",
        "label_coverage",
        "fully_labeled",
        "distinct_truth_labels",
        "majority_label",
        "majority_count",
        "purity",
        "evaluation_status",
        "routing_core",
        "backbone",
        "consensus",
        "truth",
        "core_exact",
        "core_edit_distance",
        "consensus_exact",
        "consensus_edit_distance",
        "comparison",
    ]

    with cluster_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    with summary_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            summary,
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    print()
    print(
        "=== Phase 14.3B Accuracy ==="
    )

    print(
        "assigned label coverage:",
        f"{clustering['assigned_label_coverage']:.6%}",
    )

    print(
        "labeled weighted purity:",
        f"{clustering['weighted_purity']:.6%}",
    )

    print(
        "majority tie clusters:",
        clustering[
            "majority_tie_clusters"
        ],
    )

    for name in (
        "all_evaluable",
        "multi_unique",
        "pure",
        "pure_multi_unique",
        "mixed",
        "fully_labeled_pure_multi_unique",
    ):
        result = summaries[name]

        print()
        print(f"[{name}]")
        print(
            "clusters:",
            f"{result['cluster_count']:,}",
        )

        if result["cluster_count"] == 0:
            continue

        print(
            "core exact:",
            f"{result['core_exact_rate']:.6%}",
        )
        print(
            "consensus exact:",
            f"{result['consensus_exact_rate']:.6%}",
        )
        print(
            "exact gain:",
            f"{result['exact_rate_gain']:+.6%}",
        )
        print(
            "core mean ED:",
            f"{result['core_mean_edit_distance']:.6f}",
        )
        print(
            "consensus mean ED:",
            f"{result['consensus_mean_edit_distance']:.6f}",
        )
        print(
            "mean ED reduction:",
            f"{result['mean_edit_distance_reduction']:+.6f}",
        )
        print(
            "better / same / worse:",
            result["better"],
            "/",
            result["same"],
            "/",
            result["worse"],
        )

    print()
    print(
        "algorithm wall:",
        f"{algorithm_wall:.4f}s",
    )

    print(
        "truth join:",
        f"{truth_join_seconds:.4f}s",
    )

    print(
        "evaluation:",
        f"{evaluation_seconds:.4f}s",
    )

    print(
        "summary:",
        summary_path,
    )

    print(
        "clusters:",
        cluster_path,
    )

    print(
        "truth-blind memberships:",
        membership_path,
    )


if __name__ == "__main__":
    main()
