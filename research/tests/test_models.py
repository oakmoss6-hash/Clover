import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from clover_consensus.models import ClusterRecord, ReadRecord


class TestModels(unittest.TestCase):
    def test_read_record_validates_sequence_and_count(self):
        with self.assertRaisesRegex(ValueError, "non-empty"):
            ReadRecord("r1", "")
        with self.assertRaisesRegex(ValueError, "invalid DNA"):
            ReadRecord("r1", "ACGX")
        with self.assertRaisesRegex(ValueError, ">= 1"):
            ReadRecord("r1", "ACGT", count=0)

    def test_cluster_properties(self):
        cluster = ClusterRecord(
            cluster_id="c1",
            core_sequence="ACGT",
            reads=[ReadRecord("r1", "ACGT"), ReadRecord("r2", "ACGA", count=2)],
        )
        self.assertEqual(cluster.raw_read_count, 3)
        self.assertEqual(cluster.unique_sequence_count, 2)
        self.assertEqual(cluster.sequence_lengths, (4, 4))

    def test_empty_cluster_prevention(self):
        with self.assertRaisesRegex(ValueError, "non-empty"):
            ClusterRecord("c1", "ACGT", [])


if __name__ == "__main__":
    unittest.main()
