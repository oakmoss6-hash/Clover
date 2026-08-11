from __future__ import annotations

import importlib.util
import random
import sys
import unittest
from pathlib import Path

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESEARCH_ROOT / 'src'))

from clover_consensus.aligners import align_global
from clover_consensus.models import ClusterRecord, ReadRecord
from clover_consensus.reconstruct import reconstruct_cluster

HAS_PYWFA = importlib.util.find_spec('pywfa') is not None
HAS_EDLIB = importlib.util.find_spec('edlib') is not None


@unittest.skipUnless(HAS_PYWFA, 'optional dependency pywfa is not installed')
class TestWfaBackend(unittest.TestCase):
    def test_exact(self):
        result = align_global('ACGT', 'ACGT', backend='wfa')
        self.assertEqual(result.edit_distance, 0)
        self.assertEqual(result.backend, 'wfa')
        self.assertEqual(result.backend_score, 0)

    def test_substitution(self):
        result = align_global('ACGT', 'ACGA', backend='wfa')
        self.assertEqual(result.edit_distance, 1)
        self.assertEqual(result.substitutions, 1)

    def test_insertion_direction(self):
        result = align_global('ACGT', 'ACGTT', backend='wfa')
        self.assertEqual(result.edit_distance, 1)
        self.assertEqual(result.insertions, 1)
        self.assertEqual(result.deletions, 0)

    def test_deletion_direction(self):
        result = align_global('ACGTT', 'ACGT', backend='wfa')
        self.assertEqual(result.edit_distance, 1)
        self.assertEqual(result.insertions, 0)
        self.assertEqual(result.deletions, 1)

    @unittest.skipUnless(HAS_EDLIB, 'optional dependency edlib is not installed')
    def test_random_distance_matches_edlib(self):
        rng = random.Random(20260810)
        alphabet = 'ACGT'

        def mutate(sequence: str) -> str:
            out = []
            for base in sequence:
                if rng.random() < 0.03:
                    continue
                if rng.random() < 0.03:
                    out.append(rng.choice([x for x in alphabet if x != base]))
                else:
                    out.append(base)
                if rng.random() < 0.03:
                    out.append(rng.choice(alphabet))
            return ''.join(out)

        for index in range(250):
            reference = ''.join(rng.choice(alphabet) for _ in range(150))
            query = mutate(reference)
            wfa = align_global(reference, query, backend='wfa')
            edlib = align_global(reference, query, backend='edlib')
            with self.subTest(index=index):
                self.assertEqual(wfa.edit_distance, edlib.edit_distance)
                self.assertEqual(wfa.backend_score, edlib.edit_distance)

    def test_reconstruction_accepts_wfa_backend(self):
        cluster = ClusterRecord(
            cluster_id='phase11',
            core_sequence='ACGT',
            core_read_id='core',
            reads=[
                ReadRecord('core', 'ACGT'),
                ReadRecord('a', 'AGGT'),
                ReadRecord('b', 'AGGT'),
            ],
        )
        result = reconstruct_cluster(cluster, backend='wfa')
        self.assertEqual(result.backend, 'wfa')
        self.assertEqual(result.consensus, 'AGGT')
        self.assertEqual(result.pairwise_alignment_count, 1)


if __name__ == '__main__':
    unittest.main()
