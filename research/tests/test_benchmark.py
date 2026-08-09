import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.benchmark import (
    BENCHMARK_RECORD_COLUMNS,
    ClusterBenchmarkRecord,
    DuplicateTruthClusterError,
    TruthFormatError,
    benchmark_clusters,
    read_truth_tsv,
    summarize_benchmark,
    write_benchmark_records_csv,
    write_benchmark_summary_json,
)
from clover_consensus.metrics import evaluate_sequence
from clover_consensus.models import (
    ClusterRecord,
    ReadRecord,
)
from clover_consensus.reconstruct import reconstruct_cluster


def make_cluster(
    cluster_id: str,
    core: str,
    variants: list[tuple[str, str, int]] | None = None,
) -> ClusterRecord:
    reads = [
        ReadRecord(
            read_id=f"{cluster_id}-core",
            sequence=core,
        )
    ]

    for read_id, sequence, count in variants or []:
        reads.append(
            ReadRecord(
                read_id=read_id,
                sequence=sequence,
                count=count,
            )
        )

    return ClusterRecord(
        cluster_id=cluster_id,
        core_sequence=core,
        reads=reads,
        core_read_id=f"{cluster_id}-core",
    )


def make_record(
    cluster_id: str,
    *,
    core_distance: int,
    consensus_distance: int,
    raw_read_count: int = 3,
    unique_sequence_count: int = 2,
    pairwise_alignment_count: int = 1,
    runtime_seconds: float = 0.25,
) -> ClusterBenchmarkRecord:
    return ClusterBenchmarkRecord(
        cluster_id=cluster_id,
        truth_sequence="ACGT",
        core_sequence="ACGT",
        consensus="ACGT",
        raw_read_count=raw_read_count,
        unique_sequence_count=unique_sequence_count,
        pairwise_alignment_count=(
            pairwise_alignment_count
        ),
        core_exact=core_distance == 0,
        core_edit_distance=core_distance,
        core_normalized_edit_distance=(
            core_distance / 4
        ),
        consensus_exact=consensus_distance == 0,
        consensus_edit_distance=(
            consensus_distance
        ),
        consensus_normalized_edit_distance=(
            consensus_distance / 4
        ),
        consensus_substitutions=consensus_distance,
        consensus_insertions=0,
        consensus_deletions=0,
        length_error=0,
        edit_distance_improvement=(
            core_distance - consensus_distance
        ),
        reconstruction_backend="nw",
        metric_backend="nw",
        normalize_indels=False,
        runtime_seconds=runtime_seconds,
    )


def write_temp(text: str) -> Path:
    temporary = tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        encoding="utf-8",
    )

    temporary.write(text)
    temporary.close()

    return Path(temporary.name)


class TestBenchmarkRunner(unittest.TestCase):
    def test_one_exact_cluster(self):
        records = benchmark_clusters(
            [make_cluster("c1", "ACGT")],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )

        self.assertEqual(len(records), 1)
        self.assertTrue(records[0].core_exact)
        self.assertTrue(records[0].consensus_exact)
        self.assertEqual(
            records[0].consensus_edit_distance,
            0,
        )

    def test_consensus_improves_over_core(self):
        cluster = make_cluster(
            "c1",
            "AGGT",
            [("r1", "ACGT", 2)],
        )

        record = benchmark_clusters(
            [cluster],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertEqual(record.core_edit_distance, 1)
        self.assertEqual(
            record.consensus_edit_distance,
            0,
        )
        self.assertEqual(
            record.edit_distance_improvement,
            1,
        )

    def test_consensus_unchanged(self):
        record = benchmark_clusters(
            [make_cluster("c1", "AGGT")],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertEqual(
            record.core_edit_distance,
            record.consensus_edit_distance,
        )
        self.assertEqual(
            record.edit_distance_improvement,
            0,
        )

    def test_consensus_worsens(self):
        cluster = make_cluster(
            "c1",
            "ACGT",
            [("r1", "AGGT", 2)],
        )

        record = benchmark_clusters(
            [cluster],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertEqual(record.core_edit_distance, 0)
        self.assertEqual(
            record.consensus_edit_distance,
            1,
        )
        self.assertEqual(
            record.edit_distance_improvement,
            -1,
        )

    def test_core_and_consensus_metrics_are_independent(self):
        cluster = make_cluster(
            "c1",
            "AGGT",
            [("r1", "ACGT", 2)],
        )

        record = benchmark_clusters(
            [cluster],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertFalse(record.core_exact)
        self.assertTrue(record.consensus_exact)
        self.assertEqual(record.core_edit_distance, 1)
        self.assertEqual(
            record.consensus_edit_distance,
            0,
        )

    def test_pairwise_alignment_count_is_copied(self):
        cluster = make_cluster(
            "c1",
            "ACGT",
            [
                ("r1", "AGGT", 1),
                ("r2", "AGGT", 3),
            ],
        )

        record = benchmark_clusters(
            [cluster],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertEqual(
            record.pairwise_alignment_count,
            1,
        )

    def test_runtime_is_non_negative(self):
        record = benchmark_clusters(
            [make_cluster("c1", "ACGT")],
            {"c1": "ACGT"},
            reconstruction_backend="nw",
            metric_backend="nw",
        )[0]

        self.assertGreaterEqual(
            record.runtime_seconds,
            0.0,
        )

    def test_summary_counts_and_rates(self):
        records = [
            make_record(
                "unchanged",
                core_distance=0,
                consensus_distance=0,
            ),
            make_record(
                "improved",
                core_distance=2,
                consensus_distance=0,
            ),
            make_record(
                "worsened",
                core_distance=1,
                consensus_distance=3,
            ),
        ]

        summary = summarize_benchmark(records)

        self.assertEqual(summary.cluster_count, 3)
        self.assertEqual(summary.improved_count, 1)
        self.assertEqual(summary.unchanged_count, 1)
        self.assertEqual(summary.worsened_count, 1)
        self.assertEqual(summary.core_exact_count, 1)
        self.assertEqual(
            summary.consensus_exact_count,
            2,
        )
        self.assertAlmostEqual(
            summary.core_exact_rate,
            1 / 3,
        )
        self.assertAlmostEqual(
            summary.consensus_exact_rate,
            2 / 3,
        )

    def test_summary_mean_metrics(self):
        records = [
            make_record(
                "c1",
                core_distance=0,
                consensus_distance=0,
            ),
            make_record(
                "c2",
                core_distance=1,
                consensus_distance=1,
            ),
            make_record(
                "c3",
                core_distance=2,
                consensus_distance=3,
            ),
        ]

        summary = summarize_benchmark(records)

        self.assertAlmostEqual(
            summary.mean_core_edit_distance,
            1.0,
        )
        self.assertAlmostEqual(
            summary.mean_consensus_edit_distance,
            4 / 3,
        )
        self.assertAlmostEqual(
            summary.mean_core_normalized_edit_distance,
            0.25,
        )
        self.assertAlmostEqual(
            summary.mean_consensus_normalized_edit_distance,
            1 / 3,
        )

    def test_summary_median_metrics(self):
        records = [
            make_record(
                "c1",
                core_distance=0,
                consensus_distance=0,
            ),
            make_record(
                "c2",
                core_distance=1,
                consensus_distance=1,
            ),
            make_record(
                "c3",
                core_distance=4,
                consensus_distance=3,
            ),
        ]

        summary = summarize_benchmark(records)

        self.assertEqual(
            summary.median_core_edit_distance,
            1.0,
        )
        self.assertEqual(
            summary.median_consensus_edit_distance,
            1.0,
        )
        self.assertEqual(
            summary.median_core_normalized_edit_distance,
            0.25,
        )
        self.assertEqual(
            summary.median_consensus_normalized_edit_distance,
            0.25,
        )

    def test_missing_truth_raises(self):
        with self.assertRaisesRegex(
            ValueError,
            "missing ground truth",
        ):
            benchmark_clusters(
                [make_cluster("missing", "ACGT")],
                {},
                reconstruction_backend="nw",
                metric_backend="nw",
            )

    def test_unused_truth_entry_does_not_raise(self):
        records = benchmark_clusters(
            [make_cluster("c1", "ACGT")],
            {
                "c1": "ACGT",
                "unused": "TGCA",
            },
            reconstruction_backend="nw",
            metric_backend="nw",
        )

        self.assertEqual(len(records), 1)

    def test_max_clusters_limits_and_preserves_order(self):
        clusters = [
            make_cluster("c3", "ACGT"),
            make_cluster("c1", "ACGT"),
            make_cluster("c2", "ACGT"),
        ]

        records = benchmark_clusters(
            clusters,
            {
                "c1": "ACGT",
                "c2": "ACGT",
                "c3": "ACGT",
            },
            reconstruction_backend="nw",
            metric_backend="nw",
            max_clusters=2,
        )

        self.assertEqual(
            [record.cluster_id for record in records],
            ["c3", "c1"],
        )

    def test_invalid_max_clusters_is_rejected(self):
        for value in (0, -1):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    benchmark_clusters(
                        [],
                        {},
                        max_clusters=value,
                    )

    def test_empty_summary_is_rejected(self):
        with self.assertRaises(ValueError):
            summarize_benchmark([])

    def test_options_reach_algorithm_layers(self):
        cluster = make_cluster("c1", "ACGT")

        with patch(
            "clover_consensus.benchmark.reconstruct_cluster",
            wraps=reconstruct_cluster,
        ) as mocked_reconstruct:
            with patch(
                "clover_consensus.benchmark.evaluate_sequence",
                wraps=evaluate_sequence,
            ) as mocked_evaluate:
                benchmark_clusters(
                    [cluster],
                    {"c1": "ACGT"},
                    reconstruction_backend="nw",
                    metric_backend="nw",
                    normalize_indels=True,
                )

        self.assertEqual(
            mocked_reconstruct.call_args.kwargs["backend"],
            "nw",
        )
        self.assertTrue(
            mocked_reconstruct.call_args.kwargs[
                "normalize_indels"
            ]
        )
        self.assertEqual(mocked_evaluate.call_count, 2)
        self.assertTrue(
            all(
                call.kwargs["backend"] == "nw"
                for call in mocked_evaluate.call_args_list
            )
        )

    def test_csv_writer_emits_stable_columns(self):
        record = make_record(
            "c1",
            core_distance=0,
            consensus_distance=0,
        )

        path = Path(
            tempfile.NamedTemporaryFile(
                delete=False,
            ).name
        )

        write_benchmark_records_csv(
            path,
            [record],
        )

        with path.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as source:
            rows = list(csv.reader(source))

        self.assertEqual(
            tuple(rows[0]),
            BENCHMARK_RECORD_COLUMNS,
        )
        self.assertEqual(
            rows[1][
                BENCHMARK_RECORD_COLUMNS.index(
                    "core_exact"
                )
            ],
            "true",
        )

    def test_summary_json_writer_produces_valid_json(self):
        summary = summarize_benchmark(
            [
                make_record(
                    "c1",
                    core_distance=0,
                    consensus_distance=0,
                )
            ]
        )

        path = Path(
            tempfile.NamedTemporaryFile(
                delete=False,
            ).name
        )

        write_benchmark_summary_json(
            path,
            summary,
        )

        with path.open(
            "r",
            encoding="utf-8",
        ) as source:
            decoded = json.load(source)

        self.assertEqual(
            decoded["cluster_count"],
            1,
        )

    def test_truth_tsv_loader(self):
        path = write_temp(
            "c1\tACGT\n"
            "c2\tTGCA\n"
        )

        self.assertEqual(
            read_truth_tsv(path),
            {
                "c1": "ACGT",
                "c2": "TGCA",
            },
        )

    def test_duplicate_truth_cluster_is_rejected(self):
        path = write_temp(
            "c1\tACGT\n"
            "c1\tTGCA\n"
        )

        with self.assertRaises(
            DuplicateTruthClusterError
        ):
            read_truth_tsv(path)

    def test_invalid_truth_dna_is_rejected(self):
        path = write_temp("c1\tACGX\n")

        with self.assertRaises(
            TruthFormatError
        ):
            read_truth_tsv(path)

    def test_blank_truth_line_is_ignored(self):
        path = write_temp(
            "\n"
            "c1\tACGT\n"
            "   \n"
            "c2\tTGCA\n"
        )

        self.assertEqual(
            list(read_truth_tsv(path)),
            ["c1", "c2"],
        )


if __name__ == "__main__":
    unittest.main()
