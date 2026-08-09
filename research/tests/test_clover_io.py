import tempfile
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from clover_consensus.clover_io import DuplicateReadIdError, MalformedInputError, read_clover_input


class TestCloverInput(unittest.TestCase):
    def write_temp(self, text):
        tmp = tempfile.NamedTemporaryFile("w", delete=False)
        tmp.write(text)
        tmp.close()
        return Path(tmp.name)

    def test_reads_two_column_input(self):
        path = self.write_temp("r1 ACGT\nr2 TGCA\n")
        reads = read_clover_input(path)
        self.assertEqual(reads["r1"].sequence, "ACGT")
        self.assertEqual(reads["r1"].tag, "r1")

    def test_duplicate_sequences_with_different_ids_are_allowed(self):
        path = self.write_temp("r1 ACGT\nr2 ACGT\n")
        reads = read_clover_input(path)
        self.assertEqual(len(reads), 2)

    def test_duplicate_read_id_is_rejected(self):
        path = self.write_temp("r1 ACGT\nr1 TGCA\n")
        with self.assertRaises(DuplicateReadIdError):
            read_clover_input(path)

    def test_missing_read_id_or_sequence_is_rejected(self):
        path = self.write_temp("ACGT\n")
        with self.assertRaises(MalformedInputError):
            read_clover_input(path)


if __name__ == "__main__":
    unittest.main()
