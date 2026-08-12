"""Incremental adapter from Clover cluster membership to ClusterState."""

from __future__ import annotations

from dataclasses import dataclass, field

from .cluster_state import ClusterState


@dataclass
class CloverClusterStateAdapter:
    """
    Build compact reconstruction states directly from Clover decisions.

    Clover remains responsible for routing.  This adapter only receives
    already-decided cluster memberships and never modifies Clover's tries,
    thresholds, cores, or routing decisions.
    """

    states: dict[int, ClusterState] = field(default_factory=dict)

    def record_membership(
        self,
        *,
        core_index: int,
        sequence: str,
        read_id: str | None = None,
        is_new_core: bool = False,
    ) -> None:
        if is_new_core:
            if core_index in self.states:
                raise RuntimeError(
                    f"duplicate creation of Clover cluster {core_index}"
                )

            state = ClusterState(
                cluster_id=str(core_index),
                routing_core=sequence,
            )
            state.add_sequence(sequence)
            self.states[core_index] = state
            return

        state = self.states.get(core_index)

        if state is None:
            raise RuntimeError(
                "Clover membership arrived before its routing core: "
                f"cluster={core_index}"
            )

        state.add_sequence(sequence)

    def get_state(self, core_index: int) -> ClusterState:
        return self.states[core_index]

    def validate_all(self) -> None:
        for state in self.states.values():
            state.validate()

    @property
    def cluster_count(self) -> int:
        return len(self.states)

    @property
    def raw_read_count(self) -> int:
        return sum(
            state.raw_read_count
            for state in self.states.values()
        )

    @property
    def unique_sequence_count(self) -> int:
        return sum(
            state.unique_sequence_count
            for state in self.states.values()
        )
