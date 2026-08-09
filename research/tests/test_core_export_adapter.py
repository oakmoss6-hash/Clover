import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src"),
)

from clover_consensus.cluster_adapter import (
    CoreSequenceMismatchError,
    DuplicateCoreRecordError,
    MalformedCoreExportError,
    build_cluster_records,
    read_core_sequences,
)


class TestCoreExportAdapter(unittest.TestCase):

    def write_temp(self, text: str) -> Path:
        tmp = tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            encoding="utf-8",
        )
        tmp.write(text)
        tmp.close()
        return Path(tmp.name)

    def test_one_cluster_core_export(self):
        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n"
        )

        cores = read_core_sequences(core_path)

        self.assertEqual(cores["c1"].cluster_id, "c1")
        self.assertEqual(cores["c1"].core_read_id, "r1")
        self.assertEqual(cores["c1"].sequence, "ACGT")

    def test_two_cluster_core_export(self):
        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n"
            "c2\tr3\tTGCA\n"
        )

        cores = read_core_sequences(core_path)

        self.assertEqual(sorted(cores), ["c1", "c2"])
        self.assertEqual(cores["c1"].core_read_id, "r1")
        self.assertEqual(cores["c2"].core_read_id, "r3")

    def test_duplicate_core_rejected(self):
        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n"
            "c1\tr2\tTGCA\n"
        )

        with self.assertRaises(DuplicateCoreRecordError):
            read_core_sequences(core_path)

    def test_malformed_core_tsv_rejected(self):
        core_path = self.write_temp(
            "cluster\tread\tseq\n"
            "c1\tr1\tACGT\n"
        )

        with self.assertRaises(MalformedCoreExportError):
            read_core_sequences(core_path)

    def test_sequence_mismatch_rejected(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
            "r2 ACGA\n"
        )

        membership_path = self.write_temp(
            "[('r1', 'c1'), ('r2', 'c1')]"
        )

        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tTGCA\n"
        )

        cores = read_core_sequences(core_path)

        with self.assertRaises(CoreSequenceMismatchError):
            build_cluster_records(
                input_path,
                membership_path,
                cores,
            )

    def test_core_already_in_membership_does_not_duplicate(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
            "r2 ACGA\n"
        )

        membership_path = self.write_temp(
            "[('r1', 'c1'), ('r2', 'c1')]"
        )

        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n"
        )

        result = build_cluster_records(
            input_path,
            membership_path,
            read_core_sequences(core_path),
        )

        cluster = result.clusters[0]

        self.assertEqual(
            [read.read_id for read in cluster.reads],
            ["r1", "r2"],
        )
        self.assertEqual(cluster.core_read_id, "r1")

    def test_core_absent_from_membership_is_added(self):
        input_path = self.write_temp(
            "r1 ACGT\n"
            "r2 ACGA\n"
        )

        membership_path = self.write_temp(
            "[('r2', 'c1')]"
        )

        core_path = self.write_temp(
            "cluster_id\tcore_read_id\tcore_sequence\n"
            "c1\tr1\tACGT\n"
        )

        result = build_cluster_records(
            input_path,
            membership_path,
            read_core_sequences(core_path),
        )

        cluster = result.clusters[0]

        self.assertEqual(
            [read.read_id for read in cluster.reads],
            ["r1", "r2"],
        )
        self.assertEqual(cluster.core_read_id, "r1")


if __name__ == "__main__":
    unittest.main()
