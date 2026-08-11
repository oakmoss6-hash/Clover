from __future__ import annotations

import unittest

from clover_consensus.backbone import (
    select_backbone,
    weighted_median_length,
)
from clover_consensus.cluster_state import ClusterState
from clover_consensus.models import ClusterRecord, ReadRecord
from clover_consensus.reconstruct import (
    reconstruct_cluster,
    reconstruct_state,
)


class TestClusterState(unittest.TestCase):
    def make_cluster(self) -> ClusterRecord:
        return ClusterRecord(
            cluster_id="cluster",
            core_sequence="ACGT",
            core_read_id="core",
            reads=[
                ReadRecord("core", "ACGT", count=2),
                ReadRecord("a", "ACGTT", count=3),
                ReadRecord("b", "ACGTT", count=4),
                ReadRecord("c", "ACGG", count=5),
            ],
        )

    def test_state_compresses_duplicates_while_preserving_weight(self):
        state = ClusterState.from_cluster(self.make_cluster())
        self.assertEqual(state.raw_read_count, 14)
        self.assertEqual(state.total_weight, 14)
        self.assertEqual(state.unique_sequence_count, 3)
        self.assertEqual(state.sequence_weights["ACGTT"], 7)

    def test_streaming_add_sequence_matches_cluster_adapter(self):
        state = ClusterState("cluster", "ACGT")
        state.add_sequence("ACGT", 2)
        state.add_sequence("ACGTT", 3)
        state.add_sequence("ACGTT", 4)
        state.add_sequence("ACGG", 5)
        state.validate()
        expected = ClusterState.from_cluster(self.make_cluster())
        self.assertEqual(state.sequence_weights, expected.sequence_weights)
        self.assertEqual(state.raw_read_count, expected.raw_read_count)

    def test_reconstruct_cluster_preserves_historical_default(self):
        cluster = self.make_cluster()
        adapted = reconstruct_cluster(cluster, backend="nw")
        direct = reconstruct_state(
            ClusterState.from_cluster(cluster),
            backend="nw",
            backbone_policy="core",
        )
        self.assertEqual(adapted.backbone, "ACGT")
        self.assertEqual(adapted.consensus, direct.consensus)
        self.assertEqual(
            adapted.pairwise_alignment_count,
            adapted.unique_sequence_count - 1,
        )


class TestBackboneSelection(unittest.TestCase):
    def test_weighted_median_uses_raw_multiplicity(self):
        self.assertEqual(
            weighted_median_length({"AA": 1, "CCCC": 5, "GGGGGG": 1}),
            4,
        )

    def test_support_length_is_truth_blind_and_deterministic(self):
        # Highest support ties at 2. The weighted median length is 4, so
        # the length-4 modal candidate is selected without a target length.
        weights = {"AA": 2, "CCCC": 2, "TTTT": 1}
        self.assertEqual(
            select_backbone(
                weights,
                core_sequence="AA",
                policy="support_length",
            ),
            "CCCC",
        )

    def test_max_span_prefers_longer_then_more_supported(self):
        weights = {"AAAA": 10, "CCCCCC": 1, "GGGGGG": 3}
        self.assertEqual(
            select_backbone(
                weights,
                core_sequence="AAAA",
                policy="max_span",
            ),
            "GGGGGG",
        )

    def test_reconstruction_keeps_u_minus_one_under_new_policy(self):
        cluster = ClusterRecord(
            cluster_id="cluster",
            core_sequence="ACG",
            core_read_id="core",
            reads=[
                ReadRecord("core", "ACG", count=1),
                ReadRecord("a", "ACGT", count=3),
                ReadRecord("b", "ACGT", count=2),
                ReadRecord("c", "ACGG", count=2),
            ],
        )
        result = reconstruct_cluster(
            cluster,
            backend="nw",
            backbone_policy="support_length",
        )
        self.assertEqual(result.backbone, "ACGT")
        self.assertEqual(
            result.pairwise_alignment_count,
            result.unique_sequence_count - 1,
        )


if __name__ == "__main__":
    unittest.main()
