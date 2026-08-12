from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
RESEARCH_SRC = ROOT / "research" / "src"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(RESEARCH_SRC))

from clover_consensus.clover_worker_reconstruction import (
    CloverWorkerReconstructor,
)


HAS_PYWFA = importlib.util.find_spec("pywfa") is not None


class RecordingQueue:
    def __init__(self):
        self.values = []

    def put(self, value):
        self.values.append(value)


def make_process(q_output):
    with patch.object(sys, "argv", [sys.argv[0]]):
        from clover.main import MyProcess

        process = MyProcess(
            "phase14-worker",
            [],
            q_output,
        )

    return process


class TestWorkerFinalizer(unittest.TestCase):
    def test_finalizer_disabled_preserves_output(self):
        q = RecordingQueue()
        process = make_process(q)

        process.num_dict = {
            "existing": 7,
        }

        process._finalize_worker_outputs()

        self.assertEqual(
            process.num_dict,
            {"existing": 7},
        )

    def test_finalizer_merges_without_overwrite(self):
        q = RecordingQueue()
        process = make_process(q)

        process.num_dict = {
            "existing": 7,
        }

        process.worker_finalize_observer = (
            lambda: {"extra": 11}
        )

        process._finalize_worker_outputs()

        self.assertEqual(
            process.num_dict,
            {
                "existing": 7,
                "extra": 11,
            },
        )

    def test_finalizer_rejects_key_collision(self):
        q = RecordingQueue()
        process = make_process(q)

        process.num_dict = {
            "existing": 7,
        }

        process.worker_finalize_observer = (
            lambda: {"existing": 99}
        )

        with self.assertRaises(RuntimeError):
            process._finalize_worker_outputs()

    def test_reconstructor_attaches_both_worker_hooks(self):
        q = RecordingQueue()
        process = make_process(q)

        reconstructor = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
        )

        reconstructor.attach_to_process(process)

        self.assertIsNotNone(
            process.cluster_membership_observer
        )
        self.assertIsNotNone(
            process.worker_finalize_observer
        )

        self.assertEqual(
            process.cluster_membership_observer.__self__,
            reconstructor,
        )
        self.assertEqual(
            process.worker_finalize_observer.__self__,
            reconstructor,
        )

    def test_reconstructor_rejects_double_attachment(self):
        q = RecordingQueue()
        process = make_process(q)

        first = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
        )
        second = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
        )

        first.attach_to_process(process)

        with self.assertRaises(RuntimeError):
            second.attach_to_process(process)

    def test_singleton_cluster_skips_reconstruction_engine(self):
        reconstructor = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
            backbone_policy="core",
        )

        # M = 3 but U = 1.
        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAA",
            is_new_core=True,
        )
        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAA",
        )
        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAA",
        )

        with patch(
            "clover_consensus.clover_worker_reconstruction."
            "reconstruct_state"
        ) as mocked_reconstruct:
            output = reconstructor.finalize()

        mocked_reconstruct.assert_not_called()

        self.assertEqual(
            output["Areconstruction_cluster_count"],
            1,
        )
        self.assertEqual(
            output["Areconstruction_raw_read_count"],
            3,
        )
        self.assertEqual(
            output["Areconstruction_unique_sequence_count"],
            1,
        )
        self.assertEqual(
            output["Areconstruction_pairwise_alignment_count"],
            0,
        )
        self.assertEqual(
            output["Areconstruction_singleton_cluster_count"],
            1,
        )
        self.assertEqual(
            output["Areconstruction_multi_unique_cluster_count"],
            0,
        )

        rows = output["Areconstruction_results"]

        self.assertEqual(
            rows,
            [
                (
                    1,
                    "AAAAAA",
                    "AAAAAA",
                    3,
                    1,
                    0,
                )
            ],
        )

    @unittest.skipUnless(
        HAS_PYWFA,
        "optional dependency pywfa is not installed",
    )
    def test_worker_reconstructor_returns_compact_u_minus_one_result(self):
        reconstructor = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
            backbone_policy="core",
        )

        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAA",
            is_new_core=True,
        )

        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAT",
        )

        reconstructor.record_membership(
            core_index=1,
            sequence="AAAAAT",
        )

        output = reconstructor.finalize()

        self.assertEqual(
            output["Areconstruction_cluster_count"],
            1,
        )
        self.assertEqual(
            output["Areconstruction_raw_read_count"],
            3,
        )
        self.assertEqual(
            output["Areconstruction_unique_sequence_count"],
            2,
        )
        self.assertEqual(
            output["Areconstruction_pairwise_alignment_count"],
            1,
        )

        rows = output["Areconstruction_results"]

        self.assertEqual(len(rows), 1)

        (
            core_index,
            consensus,
            backbone,
            raw_count,
            unique_count,
            pairwise_count,
        ) = rows[0]

        self.assertEqual(core_index, 1)
        self.assertEqual(consensus, "AAAAAT")
        self.assertEqual(backbone, "AAAAAA")
        self.assertEqual(raw_count, 3)
        self.assertEqual(unique_count, 2)
        self.assertEqual(pairwise_count, 1)


if __name__ == "__main__":
    unittest.main()
