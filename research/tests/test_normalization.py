import importlib.util
import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.alignment_models import AlignmentResult
from clover_consensus.aligners import align_global
from clover_consensus.normalization import (
    normalize_indels_left,
)


def make_alignment(
    *,
    reference: str,
    query: str,
    aligned_reference: str,
    aligned_query: str,
    backend: str = "manual",
    backend_score: int | float | None = None,
    backend_score_name: str | None = None,
) -> AlignmentResult:
    return AlignmentResult.from_gapped_alignment(
        reference=reference,
        query=query,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        backend=backend,
        backend_score=backend_score,
        backend_score_name=backend_score_name,
    )


class TestIndelLeftNormalization(unittest.TestCase):
    def test_alignment_without_indels_is_unchanged(self):
        alignment = make_alignment(
            reference="ACGT",
            query="ACGT",
            aligned_reference="ACGT",
            aligned_query="ACGT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertIs(normalized, alignment)

    def test_substitution_only_alignment_is_unchanged(self):
        alignment = make_alignment(
            reference="ACGT",
            query="AGGT",
            aligned_reference="ACGT",
            aligned_query="AGGT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertIs(normalized, alignment)
        self.assertEqual(
            normalized.aligned_reference,
            "ACGT",
        )
        self.assertEqual(
            normalized.aligned_query,
            "AGGT",
        )

    def test_homopolymer_deletion_is_left_normalized(self):
        alignment = make_alignment(
            reference="AAAA",
            query="AAA",
            aligned_reference="AAAA",
            aligned_query="AAA-",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_reference,
            "AAAA",
        )
        self.assertEqual(
            normalized.aligned_query,
            "-AAA",
        )
        self.assertEqual(
            normalized.cigar,
            "1D3=",
        )

    def test_homopolymer_insertion_is_left_normalized(self):
        alignment = make_alignment(
            reference="AAA",
            query="AAAA",
            aligned_reference="AAA-",
            aligned_query="AAAA",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_reference,
            "-AAA",
        )
        self.assertEqual(
            normalized.aligned_query,
            "AAAA",
        )
        self.assertEqual(
            normalized.cigar,
            "1I3=",
        )

    def test_multi_base_deletion_is_left_normalized(self):
        alignment = make_alignment(
            reference="ATATAT",
            query="ATAT",
            aligned_reference="ATATAT",
            aligned_query="ATAT--",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_reference,
            "ATATAT",
        )
        self.assertEqual(
            normalized.aligned_query,
            "--ATAT",
        )
        self.assertEqual(
            normalized.cigar,
            "2D4=",
        )

    def test_multi_base_insertion_is_left_normalized(self):
        alignment = make_alignment(
            reference="ATAT",
            query="ATATAT",
            aligned_reference="ATAT--",
            aligned_query="ATATAT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_reference,
            "--ATAT",
        )
        self.assertEqual(
            normalized.aligned_query,
            "ATATAT",
        )
        self.assertEqual(
            normalized.cigar,
            "2I4=",
        )

    def test_non_equivalent_gap_does_not_move(self):
        alignment = make_alignment(
            reference="ACGT",
            query="AGT",
            aligned_reference="ACGT",
            aligned_query="A-GT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertIs(normalized, alignment)
        self.assertEqual(
            normalized.aligned_query,
            "A-GT",
        )

    def test_normalized_alignment_reconstructs_reference(self):
        alignment = make_alignment(
            reference="ATATAT",
            query="ATAT",
            aligned_reference="ATATAT",
            aligned_query="ATAT--",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_reference.replace("-", ""),
            alignment.reference,
        )

    def test_normalized_alignment_reconstructs_query(self):
        alignment = make_alignment(
            reference="ATAT",
            query="ATATAT",
            aligned_reference="ATAT--",
            aligned_query="ATATAT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.aligned_query.replace("-", ""),
            alignment.query,
        )

    def test_edit_distance_is_preserved(self):
        alignment = make_alignment(
            reference="AACGT",
            query="AGGT",
            aligned_reference="AACGT",
            aligned_query="A-GGT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.edit_distance,
            alignment.edit_distance,
        )

    def test_substitution_count_is_preserved(self):
        alignment = make_alignment(
            reference="AACGT",
            query="AGGT",
            aligned_reference="AACGT",
            aligned_query="A-GGT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            alignment.substitutions,
            1,
        )
        self.assertEqual(
            normalized.substitutions,
            alignment.substitutions,
        )

    def test_insertion_count_is_preserved(self):
        alignment = make_alignment(
            reference="ATAT",
            query="ATATAT",
            aligned_reference="ATAT--",
            aligned_query="ATATAT",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.insertions,
            alignment.insertions,
        )

    def test_deletion_count_is_preserved(self):
        alignment = make_alignment(
            reference="ATATAT",
            query="ATAT",
            aligned_reference="ATATAT",
            aligned_query="ATAT--",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.deletions,
            alignment.deletions,
        )

    def test_normalization_is_idempotent(self):
        alignment = make_alignment(
            reference="ATATAT",
            query="ATAT",
            aligned_reference="ATATAT",
            aligned_query="ATAT--",
        )

        once = normalize_indels_left(alignment)
        twice = normalize_indels_left(once)

        self.assertEqual(twice, once)
        self.assertIs(twice, once)

    def test_backend_metadata_is_preserved(self):
        alignment = make_alignment(
            reference="AAAA",
            query="AAA",
            aligned_reference="AAAA",
            aligned_query="AAA-",
            backend="custom-backend",
            backend_score=17,
            backend_score_name="custom-score",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.backend,
            "custom-backend",
        )
        self.assertEqual(
            normalized.backend_score,
            17,
        )
        self.assertEqual(
            normalized.backend_score_name,
            "custom-score",
        )

    def test_nw_normalized_result_is_valid(self):
        alignment = align_global(
            "AAAA",
            "AAA",
            backend="nw",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.reference,
            "AAAA",
        )
        self.assertEqual(
            normalized.query,
            "AAA",
        )
        self.assertEqual(
            normalized.aligned_query,
            "-AAA",
        )
        self.assertEqual(
            normalized.edit_distance,
            1,
        )

    @unittest.skipUnless(
        importlib.util.find_spec("edlib") is not None,
        "optional dependency edlib is not installed",
    )
    def test_edlib_normalized_result_is_valid(self):
        alignment = align_global(
            "AAAA",
            "AAA",
            backend="edlib",
        )

        normalized = normalize_indels_left(alignment)

        self.assertEqual(
            normalized.reference,
            "AAAA",
        )
        self.assertEqual(
            normalized.query,
            "AAA",
        )
        self.assertEqual(
            normalized.aligned_query,
            "-AAA",
        )
        self.assertEqual(
            normalized.edit_distance,
            1,
        )

    def test_equivalent_alignments_normalize_identically(self):
        right_aligned = make_alignment(
            reference="AAAA",
            query="AAA",
            aligned_reference="AAAA",
            aligned_query="AAA-",
        )

        left_aligned = make_alignment(
            reference="AAAA",
            query="AAA",
            aligned_reference="AAAA",
            aligned_query="-AAA",
        )

        normalized_right = normalize_indels_left(
            right_aligned
        )
        normalized_left = normalize_indels_left(
            left_aligned
        )

        self.assertEqual(
            normalized_right.aligned_reference,
            normalized_left.aligned_reference,
        )
        self.assertEqual(
            normalized_right.aligned_query,
            normalized_left.aligned_query,
        )
        self.assertEqual(
            normalized_right.cigar,
            normalized_left.cigar,
        )


if __name__ == "__main__":
    unittest.main()
