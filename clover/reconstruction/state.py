"""Compact streaming state for Clover cluster reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ClusterState:
    """Unique sequence evidence and multiplicities for one Clover cluster."""

    cluster_id: str
    routing_core: str
    sequence_weights: dict[str, int] = field(default_factory=dict)
    raw_read_count: int = 0

    def add_sequence(self, sequence: str, count: int = 1) -> None:
        if not sequence:
            raise ValueError("sequence must be non-empty")
        if count < 1:
            raise ValueError("count must be positive")
        self.sequence_weights[sequence] = (
            self.sequence_weights.get(sequence, 0) + count
        )
        self.raw_read_count += count

    @property
    def unique_sequence_count(self) -> int:
        return len(self.sequence_weights)

    @property
    def total_weight(self) -> int:
        return sum(self.sequence_weights.values())

    def validate(self) -> None:
        if not self.cluster_id:
            raise ValueError("cluster_id must be non-empty")
        if not self.routing_core:
            raise ValueError("routing_core must be non-empty")
        if self.routing_core not in self.sequence_weights:
            raise ValueError(
                "routing_core must be represented in sequence evidence"
            )
        if self.raw_read_count != self.total_weight:
            raise ValueError(
                "raw_read_count must equal summed sequence multiplicity"
            )


@dataclass
class CloverClusterStateAdapter:
    """Observe already-decided Clover memberships without affecting routing."""

    states: dict[int, ClusterState] = field(default_factory=dict)

    def record_membership(
        self,
        *,
        core_index: int,
        sequence: str,
        read_id: str | None = None,
        is_new_core: bool = False,
    ) -> None:
        # read_id is deliberately accepted for evaluation hooks but is not
        # retained by the reconstruction algorithm.
        del read_id

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

    def validate_all(self) -> None:
        for state in self.states.values():
            state.validate()

    @property
    def cluster_count(self) -> int:
        return len(self.states)
