import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.mutation import (
    mutate_ids_channel,
)


class TestIDSChannel(unittest.TestCase):

    def test_zero_error_channel(self):
        result = mutate_ids_channel(
            "ACGTACGT",
            seed=1,
        )

        self.assertEqual(
            result.mutated,
            "ACGTACGT",
        )

        self.assertEqual(
            result.injected_errors,
            0,
        )

        self.assertEqual(
            result.events,
            tuple(),
        )

    def test_reproducible_with_seed(self):
        first = mutate_ids_channel(
            "ACGTACGTACGT",
            substitution_rate=0.1,
            insertion_rate=0.1,
            deletion_rate=0.1,
            seed=12345,
        )

        second = mutate_ids_channel(
            "ACGTACGTACGT",
            substitution_rate=0.1,
            insertion_rate=0.1,
            deletion_rate=0.1,
            seed=12345,
        )

        self.assertEqual(
            first,
            second,
        )

    def test_all_transmitted_bases_substituted(self):
        source = "ACGTACGT"

        result = mutate_ids_channel(
            source,
            substitution_rate=1.0,
            insertion_rate=0.0,
            deletion_rate=0.0,
            seed=10,
        )

        self.assertEqual(
            result.substitutions,
            len(source),
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
            len(result.mutated),
            len(source),
        )

        for original, mutated in zip(
            source,
            result.mutated,
        ):
            self.assertNotEqual(
                original,
                mutated,
            )

    def test_all_bases_deleted(self):
        source = "ACGTACGT"

        result = mutate_ids_channel(
            source,
            substitution_rate=0.0,
            insertion_rate=0.0,
            deletion_rate=1.0,
            seed=10,
        )

        self.assertEqual(
            result.mutated,
            "",
        )

        self.assertEqual(
            result.deletions,
            len(source),
        )

        self.assertEqual(
            result.injected_errors,
            len(source),
        )

    def test_insertions_preserve_source_as_subsequence(self):
        source = "ACGT"

        result = mutate_ids_channel(
            source,
            substitution_rate=0.0,
            insertion_rate=0.8,
            deletion_rate=0.0,
            seed=42,
        )

        self.assertGreater(
            result.insertions,
            0,
        )

        iterator = iter(
            result.mutated
        )

        self.assertTrue(
            all(
                base in iterator
                for base in source
            )
        )

    def test_length_accounting(self):
        result = mutate_ids_channel(
            "ACGTACGTACGTACGT",
            substitution_rate=0.2,
            insertion_rate=0.2,
            deletion_rate=0.2,
            seed=2026,
        )

        expected_length = (
            len(result.source)
            + result.insertions
            - result.deletions
        )

        self.assertEqual(
            len(result.mutated),
            expected_length,
        )

    def test_event_count_accounting(self):
        result = mutate_ids_channel(
            "ACGT" * 20,
            substitution_rate=0.2,
            insertion_rate=0.1,
            deletion_rate=0.1,
            seed=99,
        )

        self.assertEqual(
            len(result.events),
            result.injected_errors,
        )

        self.assertEqual(
            sum(
                event.operation == "X"
                for event in result.events
            ),
            result.substitutions,
        )

        self.assertEqual(
            sum(
                event.operation == "I"
                for event in result.events
            ),
            result.insertions,
        )

        self.assertEqual(
            sum(
                event.operation == "D"
                for event in result.events
            ),
            result.deletions,
        )

    def test_invalid_rate_rejected(self):
        with self.assertRaises(ValueError):
            mutate_ids_channel(
                "ACGT",
                substitution_rate=1.1,
            )

    def test_nonterminating_insertion_rate_rejected(self):
        with self.assertRaises(ValueError):
            mutate_ids_channel(
                "ACGT",
                insertion_rate=1.0,
            )

    def test_insertion_deletion_sum_rejected(self):
        with self.assertRaises(ValueError):
            mutate_ids_channel(
                "ACGT",
                insertion_rate=0.7,
                deletion_rate=0.4,
            )

    def test_invalid_dna_rejected(self):
        with self.assertRaises(ValueError):
            mutate_ids_channel(
                "ACGU",
            )


if __name__ == "__main__":
    unittest.main()
