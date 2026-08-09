import random
import importlib.util
import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.alignment_models import (
    operations_to_cigar,
)
from clover_consensus.aligners import (
    align_global,
    align_global_nw,
)


class TestAlignmentModels(unittest.TestCase):

    def test_operations_to_cigar(self):
        operations = (
            "=",
            "=",
            "X",
            "I",
            "I",
            "D",
            "=",
        )

        self.assertEqual(
            operations_to_cigar(operations),
            "2=1X2I1D1=",
        )

    def test_empty_operations_cigar(self):
        self.assertEqual(
            operations_to_cigar(tuple()),
            "",
        )


class TestNeedlemanWunsch(unittest.TestCase):

    def test_identical_sequences(self):
        result = align_global_nw(
            "ACGT",
            "ACGT",
        )

        self.assertEqual(
            result.aligned_reference,
            "ACGT",
        )

        self.assertEqual(
            result.aligned_query,
            "ACGT",
        )

        self.assertEqual(
            result.cigar,
            "4=",
        )

        self.assertEqual(
            result.matches,
            4,
        )

        self.assertEqual(
            result.substitutions,
            0,
        )

        self.assertEqual(
            result.insertions,
            0,
        )

        self.assertEqual(
            result.deletions,
            0,
        )

        self.assertEqual(
            result.edit_distance,
            0,
        )

        self.assertEqual(
            result.normalized_edit_distance,
            0.0,
        )

    def test_single_substitution(self):
        result = align_global_nw(
            "ACGT",
            "AGGT",
        )

        self.assertEqual(
            result.cigar,
            "1=1X2=",
        )

        self.assertEqual(
            result.substitutions,
            1,
        )

        self.assertEqual(
            result.insertions,
            0,
        )

        self.assertEqual(
            result.deletions,
            0,
        )

        self.assertEqual(
            result.edit_distance,
            1,
        )

        self.assertAlmostEqual(
            result.normalized_edit_distance,
            0.25,
        )

    def test_single_insertion(self):
        result = align_global_nw(
            "ACGT",
            "ACGTT",
        )

        self.assertEqual(
            result.insertions,
            1,
        )

        self.assertEqual(
            result.deletions,
            0,
        )

        self.assertEqual(
            result.substitutions,
            0,
        )

        self.assertEqual(
            result.edit_distance,
            1,
        )

        self.assertEqual(
            result.aligned_reference.replace("-", ""),
            "ACGT",
        )

        self.assertEqual(
            result.aligned_query.replace("-", ""),
            "ACGTT",
        )

    def test_single_deletion(self):
        result = align_global_nw(
            "ACGT",
            "AGT",
        )

        self.assertEqual(
            result.cigar,
            "1=1D2=",
        )

        self.assertEqual(
            result.deletions,
            1,
        )

        self.assertEqual(
            result.insertions,
            0,
        )

        self.assertEqual(
            result.substitutions,
            0,
        )

        self.assertEqual(
            result.edit_distance,
            1,
        )

    def test_empty_reference(self):
        result = align_global_nw(
            "",
            "ACG",
        )

        self.assertEqual(
            result.cigar,
            "3I",
        )

        self.assertEqual(
            result.insertions,
            3,
        )

        self.assertEqual(
            result.edit_distance,
            3,
        )

        self.assertEqual(
            result.normalized_edit_distance,
            1.0,
        )

    def test_empty_query(self):
        result = align_global_nw(
            "ACG",
            "",
        )

        self.assertEqual(
            result.cigar,
            "3D",
        )

        self.assertEqual(
            result.deletions,
            3,
        )

        self.assertEqual(
            result.edit_distance,
            3,
        )

        self.assertEqual(
            result.normalized_edit_distance,
            1.0,
        )

    def test_both_empty(self):
        result = align_global_nw(
            "",
            "",
        )

        self.assertEqual(
            result.cigar,
            "",
        )

        self.assertEqual(
            result.edit_distance,
            0,
        )

        self.assertEqual(
            result.normalized_edit_distance,
            0.0,
        )

        self.assertEqual(
            result.identity,
            1.0,
        )

    def test_edit_distance_is_symmetric(self):
        forward = align_global_nw(
            "ACGT",
            "AGT",
        )

        reverse = align_global_nw(
            "AGT",
            "ACGT",
        )

        self.assertEqual(
            forward.edit_distance,
            reverse.edit_distance,
        )

    def test_backend_entry_point(self):
        result = align_global(
            "ACGT",
            "AGGT",
            backend="nw",
        )

        self.assertEqual(
            result.backend,
            "nw",
        )

        self.assertEqual(
            result.edit_distance,
            1,
        )

    def test_unknown_backend_rejected(self):
        with self.assertRaises(ValueError):
            align_global(
                "ACGT",
                "ACGT",
                backend="unknown",
            )

    def test_invalid_reference_symbol_rejected(self):
        with self.assertRaises(ValueError):
            align_global_nw(
                "ACGU",
                "ACGT",
            )

    def test_invalid_query_symbol_rejected(self):
        with self.assertRaises(ValueError):
            align_global_nw(
                "ACGT",
                "ACGX",
            )

    def test_iter_columns_for_insertion(self):
        result = align_global_nw(
            "ACG",
            "ATCG",
        )

        columns = list(
            result.iter_columns()
        )

        insertion_columns = [
            column
            for column in columns
            if column.operation == "I"
        ]

        self.assertEqual(
            len(insertion_columns),
            1,
        )

        insertion = insertion_columns[0]

        self.assertIsNone(
            insertion.ref_index,
        )

        self.assertIsNotNone(
            insertion.query_index,
        )

        self.assertEqual(
            insertion.ref_base,
            "-",
        )

        self.assertEqual(
            insertion.query_base,
            "T",
        )

    def test_backend_score_separate_from_canonical_metrics(self):
        result = align_global_nw(
            "ACGT",
            "AGGT",
        )

        self.assertEqual(
            result.backend_score_name,
            "unit_edit_distance",
        )

        self.assertEqual(
            result.backend_score,
            1,
        )

        self.assertEqual(
            result.edit_distance,
            1,
        )



@unittest.skipUnless(
    importlib.util.find_spec("edlib") is not None,
    "optional dependency edlib is not installed",
)
class TestEdlibBackend(unittest.TestCase):

    def test_edlib_identical(self):
        result = align_global(
            "ACGT",
            "ACGT",
            backend="edlib",
        )

        self.assertEqual(result.backend, "edlib")
        self.assertEqual(result.edit_distance, 0)
        self.assertEqual(result.cigar, "4=")

    def test_edlib_substitution(self):
        result = align_global(
            "ACGT",
            "AGGT",
            backend="edlib",
        )

        self.assertEqual(result.substitutions, 1)
        self.assertEqual(result.edit_distance, 1)

    def test_edlib_insertion_direction(self):
        result = align_global(
            "ACGT",
            "ACGTT",
            backend="edlib",
        )

        self.assertEqual(result.insertions, 1)
        self.assertEqual(result.deletions, 0)
        self.assertEqual(result.edit_distance, 1)

    def test_edlib_deletion_direction(self):
        result = align_global(
            "ACGT",
            "AGT",
            backend="edlib",
        )

        self.assertEqual(result.deletions, 1)
        self.assertEqual(result.insertions, 0)
        self.assertEqual(result.edit_distance, 1)

    def test_edlib_backend_score_matches_canonical_distance(self):
        result = align_global(
            "ACGT",
            "AGT",
            backend="edlib",
        )

        self.assertEqual(
            result.backend_score,
            result.edit_distance,
        )

        self.assertEqual(
            result.backend_score_name,
            "unit_edit_distance",
        )

    def test_nw_edlib_random_differential(self):
        rng = random.Random(20260809)
        alphabet = "ACGT"

        for _ in range(500):
            reference_length = rng.randint(0, 25)
            query_length = rng.randint(0, 25)

            reference = "".join(
                rng.choice(alphabet)
                for _ in range(reference_length)
            )

            query = "".join(
                rng.choice(alphabet)
                for _ in range(query_length)
            )

            nw = align_global(
                reference,
                query,
                backend="nw",
            )

            edlib_result = align_global(
                reference,
                query,
                backend="edlib",
            )

            self.assertEqual(
                nw.edit_distance,
                edlib_result.edit_distance,
                msg=(
                    f"reference={reference!r}, "
                    f"query={query!r}"
                ),
            )

            self.assertEqual(
                edlib_result.aligned_reference.replace("-", ""),
                reference,
            )

            self.assertEqual(
                edlib_result.aligned_query.replace("-", ""),
                query,
            )

if __name__ == "__main__":
    unittest.main()
