import importlib.util
import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.metrics import evaluate_sequence


class TestSequenceMetrics(unittest.TestCase):
    def test_exact_sequence(self):
        metrics = evaluate_sequence(
            "ACGT",
            "ACGT",
            backend="nw",
        )

        self.assertTrue(metrics.exact_match)
        self.assertEqual(metrics.edit_distance, 0)
        self.assertEqual(
            metrics.normalized_edit_distance,
            0.0,
        )
        self.assertEqual(metrics.substitutions, 0)
        self.assertEqual(metrics.insertions, 0)
        self.assertEqual(metrics.deletions, 0)

    def test_substitution(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGGT",
            backend="nw",
        )

        self.assertFalse(metrics.exact_match)
        self.assertEqual(metrics.edit_distance, 1)
        self.assertEqual(metrics.substitutions, 1)
        self.assertEqual(metrics.insertions, 0)
        self.assertEqual(metrics.deletions, 0)

    def test_insertion(self):
        metrics = evaluate_sequence(
            "ACGT",
            "ACGTA",
            backend="nw",
        )

        self.assertEqual(metrics.edit_distance, 1)
        self.assertEqual(metrics.insertions, 1)
        self.assertEqual(metrics.deletions, 0)

    def test_deletion(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGT",
            backend="nw",
        )

        self.assertEqual(metrics.edit_distance, 1)
        self.assertEqual(metrics.deletions, 1)
        self.assertEqual(metrics.insertions, 0)

    def test_normalized_edit_distance(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGGT",
            backend="nw",
        )

        self.assertAlmostEqual(
            metrics.normalized_edit_distance,
            0.25,
        )

    def test_positive_length_error(self):
        metrics = evaluate_sequence(
            "ACGT",
            "ACGTA",
            backend="nw",
        )

        self.assertEqual(metrics.reference_length, 4)
        self.assertEqual(metrics.prediction_length, 5)
        self.assertEqual(metrics.length_error, 1)

    def test_negative_length_error(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGT",
            backend="nw",
        )

        self.assertEqual(metrics.reference_length, 4)
        self.assertEqual(metrics.prediction_length, 3)
        self.assertEqual(metrics.length_error, -1)

    def test_nw_metric_backend(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGGT",
            backend="nw",
        )

        self.assertEqual(metrics.edit_distance, 1)

    @unittest.skipUnless(
        importlib.util.find_spec("edlib")
        is not None,
        "optional dependency edlib is not installed",
    )
    def test_edlib_metric_backend(self):
        metrics = evaluate_sequence(
            "ACGT",
            "AGGT",
            backend="edlib",
        )

        self.assertEqual(metrics.edit_distance, 1)
        self.assertEqual(metrics.substitutions, 1)

    def test_invalid_backend_follows_aligner_behavior(self):
        with self.assertRaises(ValueError):
            evaluate_sequence(
                "ACGT",
                "ACGT",
                backend="unknown",
            )


if __name__ == "__main__":
    unittest.main()
