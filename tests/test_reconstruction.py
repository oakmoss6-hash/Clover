from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from clover.reconstruction import (
    CloverWorkerReconstructor,
    ClusterState,
    reconstruct_state,
    write_reconstruction_output,
)
from clover.reconstruction.backbone import select_backbone


class TestProductionReconstruction(unittest.TestCase):
    def test_duplicate_compression_and_u_minus_one(self):
        state = ClusterState("1", "AAAAAA")
        state.add_sequence("AAAAAA")
        state.add_sequence("AAAAAT")
        state.add_sequence("AAAAAT")

        result = reconstruct_state(state, backend="nw")

        self.assertEqual(state.raw_read_count, 3)
        self.assertEqual(state.unique_sequence_count, 2)
        self.assertEqual(result.pairwise_alignment_count, 1)
        self.assertEqual(result.consensus, "AAAAAT")
        self.assertEqual(result.routing_core, "AAAAAA")
        self.assertEqual(result.backbone, "AAAAAA")

    def test_backbone_policy_is_configurable(self):
        weights = {"AAA": 1, "AAAA": 3, "AAAT": 3}
        self.assertEqual(
            select_backbone(
                weights,
                core_sequence="AAA",
                policy="support_length",
            ),
            "AAAA",
        )

    def test_singleton_worker_fast_path(self):
        worker = CloverWorkerReconstructor(
            worker_name="A",
            backend="nw",
        )
        worker.record_membership(
            core_index=1,
            sequence="ACGT",
            is_new_core=True,
        )
        worker.record_membership(
            core_index=1,
            sequence="ACGT",
        )

        output = worker.finalize()

        self.assertEqual(output["Areconstruction_pairwise_alignment_count"], 0)
        self.assertEqual(
            output["Areconstruction_results"],
            [(1, "ACGT", "ACGT", "ACGT", 2, 1, 0)],
        )

    def test_writer_checks_global_invariant(self):
        count_dict = {
            "Areconstruction_results": [
                (1, "AAAA", "AAAA", "AAAT", 3, 2, 1),
            ],
            "Areconstruction_cluster_count": 1,
            "Areconstruction_raw_read_count": 3,
            "Areconstruction_unique_sequence_count": 2,
            "Areconstruction_pairwise_alignment_count": 1,
        }

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "consensus.tsv"
            summary = write_reconstruction_output(
                count_dict,
                ["A"],
                path,
            )
            self.assertEqual(summary.pairwise_alignment_count, 1)
            header = path.read_text().splitlines()[0]
            self.assertIn("routing_core", header)
            self.assertIn("backbone", header)
            self.assertIn("consensus", header)


if __name__ == "__main__":
    unittest.main()
