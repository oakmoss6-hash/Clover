import unittest

from clover import main
from clover.routing import RoutingHint


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


if __name__ == "__main__":
    unittest.main()
