import tempfile
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from clover_consensus.cluster_adapter import (
    ConflictingMembershipError,
    CoreSequenceUnavailable,
    DuplicateMembershipError,
    MalformedMembershipError,
    build_cluster_records,
    read_membership,
)
from clover_consensus.models import CoreRecord


class TestClusterAdapter(unittest.TestCase):

    def write_temp(self, text: str) -> Path:
        tmp = tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            encoding="utf-8",
        )
        tmp.write(text)
        tmp.close()
        return Path(tmp.name)

    def test_single_cluster_many_reads(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
            "r2 ACGT\n"
            "r3 ACGA\n"
        )

        membership_path = self.write_temp(
            "[('r1', 'c1'), ('r2', 'c1'), ('r3', 'c1')]"
        )

        cores = {
            "c1": CoreRecord(
                cluster_id="c1",
                core_read_id="r1",
                sequence="ACGT",
            )
        }

        result = build_cluster_records(
            input_path,
            membership_path,
            cores,
        )

        self.assertEqual(result.stats.total_input_reads, 3)
        self.assertEqual(result.stats.assigned_reads, 3)
        self.assertEqual(result.stats.unassigned_reads, 0)
        self.assertEqual(result.stats.cluster_count, 1)

        cluster = result.clusters[0]

        self.assertEqual(cluster.raw_read_count, 3)
        self.assertEqual(cluster.cluster_id, "c1")
        self.assertEqual(cluster.core_read_id, "r1")
        self.assertEqual(cluster.core_sequence, "ACGT")

    def test_two_clusters_and_unassigned_read(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
            "r2 TGCA\n"
            "r3 GGGG\n"
        )

        membership_path = self.write_temp(
            "[('r1', 'c1'), ('r2', 'c2')]"
        )

        cores = {
            "c1": CoreRecord(
                cluster_id="c1",
                core_read_id="r1",
                sequence="ACGT",
            ),
            "c2": CoreRecord(
                cluster_id="c2",
                core_read_id="r2",
                sequence="TGCA",
            ),
        }

        result = build_cluster_records(
            input_path,
            membership_path,
            cores,
        )

        self.assertEqual(
            [cluster.cluster_id for cluster in result.clusters],
            ["c1", "c2"],
        )

        self.assertEqual(result.stats.unassigned_reads, 1)

        self.assertEqual(
            [cluster.core_read_id for cluster in result.clusters],
            ["r1", "r2"],
        )

    def test_malformed_membership_is_rejected(self):
        path = self.write_temp("not a list")

        with self.assertRaises(MalformedMembershipError):
            read_membership(path)

    def test_conflicting_membership_is_rejected(self):
        path = self.write_temp(
            "[('r1', 'c1'), ('r1', 'c2')]"
        )

        with self.assertRaises(ConflictingMembershipError):
            read_membership(path)

    def test_duplicate_membership_is_rejected(self):
        path = self.write_temp(
            "[('r1', 'c1'), ('r1', 'c1')]"
        )

        with self.assertRaises(DuplicateMembershipError):
            read_membership(path)

    def test_missing_read_id_in_membership_is_rejected(self):
        path = self.write_temp(
            "[('', 'c1')]"
        )

        with self.assertRaises(MalformedMembershipError):
            read_membership(path)

    def test_core_sequence_unavailable(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
        )

        membership_path = self.write_temp(
            "[('r1', 'c1')]"
        )

        with self.assertRaisesRegex(
            CoreSequenceUnavailable,
            "CORE_SEQUENCE_UNAVAILABLE",
        ):
            build_cluster_records(
                input_path,
                membership_path,
            )


if __name__ == "__main__":
    unittest.main()
