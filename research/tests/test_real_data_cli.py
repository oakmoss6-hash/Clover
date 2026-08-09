import argparse
import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import ModuleType
from unittest.mock import patch


RESEARCH_ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = (
    RESEARCH_ROOT
    / "benchmarks"
    / "run_real_data.py"
)


def load_runner_module() -> ModuleType:
    specification = (
        importlib.util.spec_from_file_location(
            "caspr_run_real_data",
            RUNNER_PATH,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            "could not load real-data CLI module"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    sys.modules[specification.name] = module
    specification.loader.exec_module(module)

    return module


run_real_data = load_runner_module()


class TestRealDataCli(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = (
            tempfile.TemporaryDirectory()
        )

        self.root = Path(
            self.temporary_directory.name
        )

        self.input_path = (
            self.root / "reads.txt"
        )

        self.membership_path = (
            self.root / "membership.txt"
        )

        self.core_export_path = (
            self.root / "cores.tsv"
        )

        self.truth_path = (
            self.root / "truth.tsv"
        )

        self.records_path = (
            self.root / "records.csv"
        )

        self.summary_path = (
            self.root / "summary.json"
        )

        self.input_path.write_text(
            "r1 ACGT\n",
            encoding="utf-8",
        )

        self.membership_path.write_text(
            "[('r1', 'c1')]",
            encoding="utf-8",
        )

        self.core_export_path.write_text(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n",
            encoding="utf-8",
        )

        self.truth_path.write_text(
            "c1\tACGT\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def arguments(
        self,
        *extra: str,
    ) -> list[str]:
        return [
            "--input",
            str(self.input_path),
            "--membership",
            str(self.membership_path),
            "--core-export",
            str(self.core_export_path),
            "--truth",
            str(self.truth_path),
            "--output-records",
            str(self.records_path),
            "--output-summary",
            str(self.summary_path),
            "--backend",
            "nw",
            "--metric-backend",
            "nw",
            *extra,
        ]

    def test_parser_accepts_all_required_arguments(self):
        parser = run_real_data.build_parser()

        args = parser.parse_args(
            self.arguments(
                "--normalize-indels",
                "--max-clusters",
                "1",
            )
        )

        self.assertIsInstance(
            parser,
            argparse.ArgumentParser,
        )
        self.assertEqual(
            args.input_path,
            str(self.input_path),
        )
        self.assertEqual(
            args.membership,
            str(self.membership_path),
        )
        self.assertEqual(
            args.core_export,
            str(self.core_export_path),
        )
        self.assertEqual(
            args.truth,
            str(self.truth_path),
        )
        self.assertEqual(args.backend, "nw")
        self.assertEqual(
            args.metric_backend,
            "nw",
        )
        self.assertTrue(args.normalize_indels)
        self.assertEqual(args.max_clusters, 1)

    def test_exact_cluster_completes_and_writes_outputs(self):
        return_code = run_real_data.main(
            self.arguments()
        )

        self.assertEqual(return_code, 0)
        self.assertTrue(self.records_path.is_file())
        self.assertTrue(self.summary_path.is_file())

        with self.records_path.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as source:
            records = list(
                csv.DictReader(source)
            )

        self.assertEqual(len(records), 1)
        self.assertEqual(
            records[0]["cluster_id"],
            "c1",
        )
        self.assertEqual(
            records[0]["core_exact"],
            "true",
        )
        self.assertEqual(
            records[0]["consensus_exact"],
            "true",
        )

        with self.summary_path.open(
            "r",
            encoding="utf-8",
        ) as source:
            summary = json.load(source)

        self.assertEqual(
            summary["cluster_count"],
            1,
        )

    def test_cli_reaches_actual_adapter_workflow(self):
        with patch.object(
            run_real_data,
            "read_core_sequences",
            wraps=run_real_data.read_core_sequences,
        ) as mocked_read_cores:
            with patch.object(
                run_real_data,
                "build_cluster_records",
                wraps=run_real_data.build_cluster_records,
            ) as mocked_build_clusters:
                run_real_data.main(
                    self.arguments()
                )

        mocked_read_cores.assert_called_once_with(
            str(self.core_export_path)
        )

        mocked_build_clusters.assert_called_once()

        build_args = mocked_build_clusters.call_args.args

        self.assertEqual(
            build_args[0],
            str(self.input_path),
        )
        self.assertEqual(
            build_args[1],
            str(self.membership_path),
        )

        core_records = build_args[2]

        self.assertIsInstance(core_records, dict)
        self.assertEqual(
            set(core_records),
            {"c1"},
        )
        self.assertEqual(
            core_records["c1"].cluster_id,
            "c1",
        )
        self.assertEqual(
            core_records["c1"].core_read_id,
            "r1",
        )
        self.assertEqual(
            core_records["c1"].sequence,
            "ACGT",
        )

    def test_max_clusters_reaches_benchmark_layer(self):
        with patch.object(
            run_real_data,
            "benchmark_clusters",
            wraps=run_real_data.benchmark_clusters,
        ) as mocked_benchmark:
            run_real_data.main(
                self.arguments(
                    "--max-clusters",
                    "1",
                )
            )

        self.assertEqual(
            mocked_benchmark.call_args.kwargs[
                "max_clusters"
            ],
            1,
        )

    def test_normalize_indels_reaches_benchmark_layer(self):
        with patch.object(
            run_real_data,
            "benchmark_clusters",
            wraps=run_real_data.benchmark_clusters,
        ) as mocked_benchmark:
            run_real_data.main(
                self.arguments(
                    "--normalize-indels",
                )
            )

        self.assertTrue(
            mocked_benchmark.call_args.kwargs[
                "normalize_indels"
            ]
        )

    def test_missing_truth_raises_clearly(self):
        self.truth_path.write_text(
            "other\tACGT\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(
            ValueError,
            "missing ground truth",
        ):
            run_real_data.main(
                self.arguments()
            )


if __name__ == "__main__":
    unittest.main()
