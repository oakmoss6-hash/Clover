import unittest

from clover import main
from clover.routing import RoutingHint


class StubTree:
    """Deterministic stand-in for Clover Trie.fuzz_fin()."""

    def __init__(self, results):
        self.results = list(results)
        self.calls = []

    def fuzz_fin(self, word, max_value):
        self.calls.append((word, max_value))
        if not self.results:
            raise AssertionError("unexpected fuzz_fin call")
        return self.results.pop(0)


class TestPassiveRoutingHint(unittest.TestCase):

    def make_process(self, capture):
        process = main.MyProcess("test", [], [])
        process.now_clust_threshold = 0
        process.read_len = 6
        process.dna_tree_nums = 4
        process.fuzz_list = [1, 1, 4]
        process.loc_nums = [0]
        process.fuzz_tree_nums = 1
        process.config_dict["read_len_min"] = 4
        process.capture_routing_hints = capture
        return process

    def test_capture_disabled_preserves_historical_behavior(self):
        process = self.make_process(False)

        process.cluster("1 AAAAAA")
        process.cluster("2 AAAAAA")

        self.assertEqual(process.ref_dict[1], ["1", "2"])
        self.assertEqual(process.routing_hints, {})

    def test_front_route_is_captured(self):
        process = self.make_process(True)

        process.cluster("1 AAAAAA")
        process.cluster("2 AAAAAA")

        self.assertEqual(process.ref_dict[1], ["1", "2"])

        hint = process.routing_hints[1]["AAAAAA"]

        self.assertIsInstance(hint, RoutingHint)
        self.assertEqual(hint.tree_kind, "front")
        self.assertEqual(hint.horizontal_drifts, 0)
        self.assertEqual(hint.query_shift, 0)

    def test_duplicate_sequence_has_one_hint(self):
        process = self.make_process(True)

        process.cluster("1 AAAAAA")
        process.cluster("2 AAAAAA")
        process.cluster("3 AAAAAA")

        self.assertEqual(process.ref_dict[1], ["1", "2", "3"])
        self.assertEqual(len(process.routing_hints[1]), 1)

    def test_capture_on_off_produces_identical_clusters(self):
        without = self.make_process(False)
        with_capture = self.make_process(True)

        reads = [
            "1 AAAAAA",
            "2 AAAAAA",
            "3 AAAAAT",
            "4 AAAAAA",
        ]

        for read in reads:
            without.cluster(read)
            with_capture.cluster(read)

        self.assertEqual(without.ref_dict, with_capture.ref_dict)
        self.assertEqual(without.index_list, with_capture.index_list)


    def make_middle_process(self):
        process = self.make_process(True)

        # Pretend core 1 already exists.  The stub trees control only the
        # routing decision; no additional alignment is involved.
        process.ref_dict = {1: ["core"]}
        process.ref_list = {1: "AAAAAA"}

        process.loc_nums = [-1, 0, 1]
        process.fuzz_tree_nums = 4
        process.now_clust_threshold = 10

        # Force front/back routes to fail so Clover reaches middle search.
        process.a_tree = StubTree([
            ["", 1000],
        ])
        process.b_tree = StubTree([
            ["", 1000],
        ])

        return process

    def test_middle_c_records_winning_query_shift(self):
        process = self.make_middle_process()

        # i=-1, 0, +1 respectively.
        process.c_tree = StubTree([
            [1, 3],
            [1, 0],
            [1, 2],
        ])
        process.d_tree = StubTree([
            [1, 3],
            [1, 2],
            [1, 3],
        ])

        process.cluster("2 AAAAAT")

        hint = process.routing_hints[1]["AAAAAT"]

        self.assertEqual(hint.tree_kind, "middle_c")
        self.assertEqual(hint.horizontal_drifts, 0)
        self.assertEqual(hint.query_shift, 0)

    def test_middle_d_records_winning_query_shift(self):
        process = self.make_middle_process()

        process.c_tree = StubTree([
            [1, 3],
            [1, 2],
            [1, 3],
        ])
        process.d_tree = StubTree([
            [1, 3],
            [1, 0],
            [1, 2],
        ])

        process.cluster("2 AAAAAT")

        hint = process.routing_hints[1]["AAAAAT"]

        self.assertEqual(hint.tree_kind, "middle_d")
        self.assertEqual(hint.horizontal_drifts, 0)
        self.assertEqual(hint.query_shift, 0)

    def test_middle_tie_keeps_first_winner(self):
        process = self.make_middle_process()

        # The first candidate reaches drift=1 at i=-1.
        # Every later result ties at 1.  Clover uses "<", not "<=",
        # so the first middle_c candidate must remain the winner.
        process.c_tree = StubTree([
            [1, 1],
            [1, 1],
            [1, 1],
        ])
        process.d_tree = StubTree([
            [1, 1],
            [1, 1],
            [1, 1],
        ])

        process.cluster("2 AAAAAT")

        hint = process.routing_hints[1]["AAAAAT"]

        self.assertEqual(hint.tree_kind, "middle_c")
        self.assertEqual(hint.horizontal_drifts, 1)
        self.assertEqual(hint.query_shift, -1)


if __name__ == "__main__":
    unittest.main()
