"""Passive routing metadata produced by Clover clustering.

Routing hints describe how Clover found a cluster candidate. They are
acceleration hints only; they are not alignment constraints or ground truth.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class RoutingHint:
    """Minimal metadata for one successful Clover routing decision.

    Attributes:
        tree_kind:
            Clover tree that produced the accepted match:
            ``front``, ``back``, ``middle_c`` or ``middle_d``.
        horizontal_drifts:
            Fuzzy-search drift count returned by ``Trie.fuzz_fin``.
        query_shift:
            Query-window shift used by Clover for a middle-tree search.
            This is not assumed to be a global-alignment offset.
    """

    tree_kind: str
    horizontal_drifts: int
    query_shift: int = 0
