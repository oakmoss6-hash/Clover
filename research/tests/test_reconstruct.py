import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.alignment_models import AlignmentResult
from clover_consensus.aligners import align_global
from clover_consensus.models import (
    ClusterRecord,
    ReadRecord,
)
from clover_consensus.reconstruct import reconstruct_cluster


def make_cluster(
    core: str,
    reads: list[ReadRecord],
    cluster_id: str = "c1",
) -> ClusterRecord:
    return ClusterRecord(
        cluster_id=cluster_id,
        core_sequence=core,
        reads=reads,
        core_read_id="core",
    )


def make_right_aligned_deletion() -> AlignmentResult:
    return AlignmentResult.from_gapped_alignment(
        reference="AAAA",
        query="AAA",
        aligned_reference="AAAA",
        aligned_query="AAA-",
        backend="manual",
        backend_score=1,
        backend_score_name="unit_edit_distance",
    )


class TestCasprReconstruction(unittest.TestCase):
    def test_core_only_cluster(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [ReadRecord("core", "ACGT")],
            )
        )

        self.assertEqual(result.consensus, "ACGT")
        self.assertEqual(
            result.pairwise_alignment_count,
            0,
        )

    def test_identical_reads_leave_consensus_unchanged(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "ACGT"),
                    ReadRecord(
                        "r2",
                        "ACGT",
                        count=3,
                    ),
                ],
            )
        )

        self.assertEqual(result.consensus, "ACGT")
        self.assertEqual(result.total_weight, 5)
        self.assertEqual(
            result.unique_sequence_count,
            1,
        )
        self.assertEqual(
            result.pairwise_alignment_count,
            0,
        )

    def test_substitution_majority_repairs_backbone(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "AGGT"),
                    ReadRecord("r2", "AGGT"),
                ],
            )
        )

        self.assertEqual(result.consensus, "AGGT")
        self.assertEqual(
            result.changed_base_count,
            1,
        )

    def test_deletion_majority_removes_backbone_position(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "AGT"),
                    ReadRecord("r2", "AGT"),
                ],
            )
        )

        self.assertEqual(result.consensus, "AGT")
        self.assertEqual(
            result.deleted_base_count,
            1,
        )

    def test_insertion_majority_adds_inserted_base(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "ACTGT"),
                    ReadRecord("r2", "ACTGT"),
                ],
            )
        )

        self.assertEqual(result.consensus, "ACTGT")
        self.assertEqual(
            result.inserted_base_count,
            1,
        )

    def test_leading_insertion_uses_slot_zero(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "TACGT"),
                    ReadRecord("r2", "TACGT"),
                ],
            )
        )

        self.assertEqual(result.consensus, "TACGT")
        self.assertEqual(
            result.insertion_profiles[
                0
            ].non_empty_counts,
            {"T": 2},
        )

    def test_trailing_insertion_uses_slot_l(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "ACGTA"),
                    ReadRecord("r2", "ACGTA"),
                ],
            )
        )

        self.assertEqual(result.consensus, "ACGTA")
        self.assertEqual(
            result.insertion_profiles[
                4
            ].non_empty_counts,
            {"A": 2},
        )

    def test_duplicate_count_weighting_changes_consensus(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord(
                        "core",
                        "ACGT",
                        count=2,
                    ),
                    ReadRecord(
                        "variant",
                        "AGGT",
                        count=3,
                    ),
                ],
            )
        )

        self.assertEqual(result.consensus, "AGGT")
        self.assertEqual(
            result.raw_read_count,
            5,
        )
        self.assertEqual(result.total_weight, 5)

    def test_duplicate_sequences_are_aligned_once(self):
        cluster = make_cluster(
            "ACGT",
            [
                ReadRecord("core", "ACGT"),
                ReadRecord("r1", "AGGT"),
                ReadRecord(
                    "r2",
                    "AGGT",
                    count=4,
                ),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            wraps=align_global,
        ) as mocked_align:
            result = reconstruct_cluster(cluster)

        self.assertEqual(
            mocked_align.call_count,
            1,
        )
        self.assertEqual(
            result.pairwise_alignment_count,
            1,
        )
        self.assertEqual(
            result.unique_sequence_count,
            2,
        )

    def test_pairwise_call_bound_with_backbone_represented(self):
        cluster = make_cluster(
            "ACGT",
            [
                ReadRecord("core", "ACGT"),
                ReadRecord("r1", "AGGT"),
                ReadRecord("r2", "ACGTT"),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            wraps=align_global,
        ) as mocked_align:
            result = reconstruct_cluster(cluster)

        self.assertLessEqual(
            mocked_align.call_count,
            result.unique_sequence_count - 1,
        )
        self.assertEqual(
            mocked_align.call_count,
            result.pairwise_alignment_count,
        )

    def test_base_tie_preserves_backbone_base(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord(
                        "variant",
                        "AGGT",
                    ),
                ],
            )
        )

        self.assertEqual(result.consensus, "ACGT")
        self.assertEqual(
            result.base_profiles[
                1
            ].consensus_symbol(),
            "C",
        )

    def test_insertion_tie_prefers_empty(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord(
                        "variant",
                        "ACTGT",
                    ),
                ],
            )
        )

        self.assertEqual(result.consensus, "ACGT")
        self.assertEqual(
            result.insertion_profiles[
                2
            ].consensus_insertion(
                result.total_weight
            ),
            "",
        )

    def test_base_weight_is_conserved_at_every_position(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord(
                        "core",
                        "ACGT",
                        count=2,
                    ),
                    ReadRecord(
                        "sub",
                        "AGGT",
                        count=3,
                    ),
                    ReadRecord(
                        "del",
                        "AGT",
                        count=4,
                    ),
                    ReadRecord(
                        "ins",
                        "ACTGT",
                        count=5,
                    ),
                ],
            )
        )

        for profile in result.base_profiles:
            self.assertEqual(
                sum(profile.counts.values()),
                result.total_weight,
            )

    def test_insertion_implicit_empty_weight(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord(
                        "core",
                        "ACGT",
                        count=2,
                    ),
                    ReadRecord(
                        "inserted",
                        "ACTGT",
                        count=3,
                    ),
                ],
            )
        )

        profile = result.insertion_profiles[2]

        self.assertEqual(
            profile.non_empty_counts,
            {"T": 3},
        )
        self.assertEqual(
            profile.empty_weight(
                result.total_weight
            ),
            2,
        )

    def test_nw_backend_reconstructs_cluster(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "AGGT"),
                    ReadRecord("r2", "AGGT"),
                ],
            ),
            backend="nw",
        )

        self.assertEqual(result.backend, "nw")
        self.assertEqual(result.consensus, "AGGT")

    @unittest.skipUnless(
        importlib.util.find_spec("edlib") is not None,
        "optional dependency edlib is not installed",
    )
    def test_edlib_backend_reconstructs_cluster(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACGT",
                [
                    ReadRecord("core", "ACGT"),
                    ReadRecord("r1", "AGGT"),
                    ReadRecord("r2", "AGGT"),
                ],
            ),
            backend="edlib",
        )

        self.assertEqual(result.backend, "edlib")
        self.assertEqual(result.consensus, "AGGT")

    def test_multi_base_insertion_is_one_string_vote(self):
        result = reconstruct_cluster(
            make_cluster(
                "ACAC",
                [
                    ReadRecord("core", "ACAC"),
                    ReadRecord("r1", "ACGTAC"),
                    ReadRecord("r2", "ACGTAC"),
                ],
            )
        )

        profile = result.insertion_profiles[2]

        self.assertEqual(
            result.consensus,
            "ACGTAC",
        )
        self.assertEqual(
            profile.non_empty_counts,
            {"GT": 2},
        )
        self.assertNotIn(
            "G",
            profile.non_empty_counts,
        )
        self.assertNotIn(
            "T",
            profile.non_empty_counts,
        )

    def test_no_read_to_read_alignment_path(self):
        backbone = "ACGT"

        cluster = make_cluster(
            backbone,
            [
                ReadRecord("core", backbone),
                ReadRecord("sub", "AGGT"),
                ReadRecord("del", "AGT"),
                ReadRecord("ins", "ACTGT"),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            wraps=align_global,
        ) as mocked_align:
            result = reconstruct_cluster(cluster)

        self.assertEqual(
            mocked_align.call_count,
            3,
        )
        self.assertEqual(
            result.pairwise_alignment_count,
            3,
        )
        self.assertTrue(
            all(
                call.args[0] == backbone
                for call in mocked_align.call_args_list
            )
        )
        self.assertEqual(
            {
                call.args[1]
                for call in mocked_align.call_args_list
            },
            {
                "AGGT",
                "AGT",
                "ACTGT",
            },
        )

    def test_missing_backbone_sequence_is_rejected(self):
        cluster = make_cluster(
            "ACGT",
            [
                ReadRecord("r1", "AGGT"),
                ReadRecord("r2", "ACGTT"),
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            "CASPR/Phase-4.5 contract",
        ):
            reconstruct_cluster(cluster)

    def test_normalization_false_preserves_v0_projection(self):
        cluster = make_cluster(
            "AAAA",
            [
                ReadRecord("core", "AAAA"),
                ReadRecord("short", "AAA"),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            return_value=make_right_aligned_deletion(),
        ):
            default_result = reconstruct_cluster(cluster)

        with patch(
            "clover_consensus.reconstruct.align_global",
            return_value=make_right_aligned_deletion(),
        ):
            explicit_result = reconstruct_cluster(
                cluster,
                normalize_indels=False,
            )

        self.assertEqual(
            default_result,
            explicit_result,
        )
        self.assertEqual(
            explicit_result.base_profiles[
                0
            ].counts["D"],
            0,
        )
        self.assertEqual(
            explicit_result.base_profiles[
                3
            ].counts["D"],
            1,
        )

    def test_normalization_precedes_profile_projection(self):
        cluster = make_cluster(
            "AAAA",
            [
                ReadRecord("core", "AAAA"),
                ReadRecord("short", "AAA"),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            return_value=make_right_aligned_deletion(),
        ):
            result = reconstruct_cluster(
                cluster,
                normalize_indels=True,
            )

        self.assertEqual(
            result.base_profiles[
                0
            ].counts["D"],
            1,
        )
        self.assertEqual(
            result.base_profiles[
                3
            ].counts["D"],
            0,
        )

    def test_normalization_does_not_add_alignment_calls(self):
        cluster = make_cluster(
            "AAAA",
            [
                ReadRecord("core", "AAAA"),
                ReadRecord("short", "AAA"),
            ],
        )

        with patch(
            "clover_consensus.reconstruct.align_global",
            return_value=make_right_aligned_deletion(),
        ) as mocked_align:
            result = reconstruct_cluster(
                cluster,
                normalize_indels=True,
            )

        self.assertEqual(
            mocked_align.call_count,
            1,
        )
        self.assertEqual(
            result.pairwise_alignment_count,
            1,
        )


if __name__ == "__main__":
    unittest.main()
