from __future__ import annotations

import importlib.util
import multiprocessing as mp
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
RESEARCH_SRC = ROOT / "research" / "src"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(RESEARCH_SRC))

from clover_consensus.clover_worker_reconstruction import (
    CloverWorkerReconstructor,
)


HAS_PYWFA = importlib.util.find_spec("pywfa") is not None


@unittest.skipUnless(
    HAS_PYWFA,
    "optional dependency pywfa is not installed",
)
class TestPhase14MultiprocessSmoke(unittest.TestCase):
    def test_real_clover_process_reconstructs_locally(self):
        # 152-nt DNA-storage-like reads matching Clover defaults.
        core = "A" * 152
        variant = ("A" * 151) + "T"

        data = [
            f"r1 {core}",
            f"r2 {variant}",
            f"r3 {variant}",
        ]

        q_output = mp.Queue(2)

        # Clover's config loader reads sys.argv.
        with patch.object(sys, "argv", [sys.argv[0]]):
            from clover.main import MyProcess

            process = MyProcess(
                "A",
                data,
                q_output,
            )

        reconstructor = CloverWorkerReconstructor(
            worker_name="A",
            backend="wfa",
            backbone_policy="core",
        )

        reconstructor.attach_to_process(process)

        process.start()

        output = q_output.get(timeout=15)

        process.join(timeout=15)

        if process.is_alive():
            process.terminate()
            process.join()
            self.fail("Clover worker did not terminate")

        self.assertEqual(process.exitcode, 0)

        # Real child process must have reconstructed one locally-owned cluster.
        self.assertEqual(
            output["Areconstruction_cluster_count"],
            1,
        )

        # M = 3 raw reads.
        self.assertEqual(
            output["Areconstruction_raw_read_count"],
            3,
        )

        # U = {core, variant} = 2.
        self.assertEqual(
            output["Areconstruction_unique_sequence_count"],
            2,
        )

        # Exact star invariant: U - C = 2 - 1 = 1.
        self.assertEqual(
            output["Areconstruction_pairwise_alignment_count"],
            1,
        )

        rows = output["Areconstruction_results"]

        self.assertEqual(len(rows), 1)

        (
            cluster_id,
            consensus,
            backbone,
            raw_count,
            unique_count,
            pairwise_count,
        ) = rows[0]

        self.assertEqual(raw_count, 3)
        self.assertEqual(unique_count, 2)
        self.assertEqual(pairwise_count, 1)

        # Routing core remains the original first read.
        self.assertEqual(backbone, core)

        # Two copies of the variant outvote the core.
        self.assertEqual(consensus, variant)


if __name__ == "__main__":
    unittest.main()
