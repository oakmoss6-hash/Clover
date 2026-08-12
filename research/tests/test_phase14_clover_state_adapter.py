from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
RESEARCH_SRC = ROOT / "research" / "src"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(RESEARCH_SRC))

from clover_consensus.clover_state_adapter import (
    CloverClusterStateAdapter,
)
from clover_consensus.reconstruct import reconstruct_state

HAS_PYWFA = importlib.util.find_spec("pywfa") is not None


class StubTree:
    def __init__(self, results):
        self.results = list(results)
        self.calls = []
        self.inserted = []

    def fuzz_fin(self, word, max_value):
        self.calls.append((word, max_value))

        if not self.results:
            raise AssertionError("unexpected fuzz_fin call")

        return self.results.pop(0)

    def insert(self, word, value):
        self.inserted.append((word, value))


def make_process():
    # Clover's config loader inspects sys.argv.
    with patch.object(sys, "argv", [sys.argv[0]]):
        from clover.main import MyProcess
        process = MyProcess("phase14-test", [], [])

    process.config_dict["Virtual_mode"] = True
    process.config_dict["read_len_min"] = 4
    process.config_dict["other_tree_nums"] = 2
    process.config_dict["tree_threshold"] = 10

    process.align_swicth = False
    process.read_len = 6
    process.dna_tree_nums = 2
    process.fuzz_list = [1, 1, 2]
    process.loc_nums = [-1, 0, 1]
    process.fuzz_tree_nums = 4
    process.now_clust_threshold = 8

    return process


def configure_first_read_as_new_core(process):
    process.a_tree = StubTree([
        ["", 1000],
    ])
    process.b_tree = StubTree([
        ["", 1000],
    ])

    process.c_tree = StubTree([
        ["", 1000],
        ["", 1000],
        ["", 1000],
    ])
    process.d_tree = StubTree([
        ["", 1000],
        ["", 1000],
        ["", 1000],
    ])


class TestCloverClusterStateAdapter(unittest.TestCase):
    def test_adapter_compresses_duplicates(self):
        adapter = CloverClusterStateAdapter()

        adapter.record_membership(
            core_index=7,
            sequence="AAAAAA",
            is_new_core=True,
        )
        adapter.record_membership(
            core_index=7,
            sequence="AAAAAT",
        )
        adapter.record_membership(
            core_index=7,
            sequence="AAAAAT",
        )

        state = adapter.get_state(7)

        self.assertEqual(state.routing_core, "AAAAAA")
        self.assertEqual(state.raw_read_count, 3)
        self.assertEqual(state.unique_sequence_count, 2)
        self.assertEqual(
            state.sequence_weights,
            {
                "AAAAAA": 1,
                "AAAAAT": 2,
            },
        )

        state.validate()

    def test_real_clover_cluster_calls_build_state(self):
        process = make_process()
        adapter = CloverClusterStateAdapter()

        process.cluster_membership_observer = (
            adapter.record_membership
        )

        # Read 1 fails every lookup and becomes Clover core 1.
        configure_first_read_as_new_core(process)
        process.cluster("r1 AAAAAA")

        # Reads 2 and 3 route to that core through the front tree.
        process.a_tree = StubTree([
            [1, 0],
            [1, 0],
        ])

        process.cluster("r2 AAAAAT")
        process.cluster("r3 AAAAAT")

        state = adapter.get_state(1)

        self.assertEqual(adapter.cluster_count, 1)
        self.assertEqual(state.routing_core, "AAAAAA")
        self.assertEqual(state.raw_read_count, 3)
        self.assertEqual(state.unique_sequence_count, 2)
        self.assertEqual(state.sequence_weights["AAAAAA"], 1)
        self.assertEqual(state.sequence_weights["AAAAAT"], 2)

        # Observer does not replace Clover's own cluster bookkeeping.
        self.assertEqual(
            process.ref_dict[1],
            ["r1", "r2", "r3"],
        )

    def test_observer_disabled_preserves_clover(self):
        process = make_process()

        self.assertIsNone(
            process.cluster_membership_observer
        )

        configure_first_read_as_new_core(process)
        process.cluster("r1 AAAAAA")

        self.assertEqual(
            process.ref_dict[1],
            ["r1"],
        )

    @unittest.skipUnless(
        HAS_PYWFA,
        "optional dependency pywfa is not installed",
    )
    def test_clover_state_reconstructs_with_u_minus_one_wfa(self):
        adapter = CloverClusterStateAdapter()

        adapter.record_membership(
            core_index=1,
            sequence="AAAAAA",
            is_new_core=True,
        )
        adapter.record_membership(
            core_index=1,
            sequence="AAAAAT",
        )
        adapter.record_membership(
            core_index=1,
            sequence="AAAAAT",
        )

        state = adapter.get_state(1)

        result = reconstruct_state(
            state,
            backend="wfa",
            backbone_policy="core",
        )

        self.assertEqual(result.raw_read_count, 3)
        self.assertEqual(result.unique_sequence_count, 2)

        # U = 2, therefore exactly U - 1 = 1 pairwise call.
        self.assertEqual(
            result.pairwise_alignment_count,
            1,
        )

        self.assertEqual(
            result.consensus,
            "AAAAAT",
        )


if __name__ == "__main__":
    unittest.main()
