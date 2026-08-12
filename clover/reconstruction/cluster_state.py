"""Compact, streaming-friendly evidence state for one Clover cluster."""

from __future__ import annotations

from dataclasses import dataclass, field

from .models import ClusterRecord


@dataclass
class ClusterState:
    """
    Minimal reconstruction state owned by one Clover cluster.

    Identical reads are represented once with a multiplicity.  This state can
    be built from the existing ``ClusterRecord`` API or incrementally by a
    future Clover integration without retaining duplicate read objects.
    """

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
                "CASPR/Phase-4.5 contract: CASPR/Phase-4.5 contract: CASPR/Phase-4.5 contract: routing_core must be represented in sequence evidence"
            )
        if self.raw_read_count != self.total_weight:
            raise ValueError(
                "raw_read_count must equal summed sequence multiplicity"
            )

    @classmethod
    def from_cluster(cls, cluster: ClusterRecord) -> "ClusterState":
        state = cls(
            cluster_id=cluster.cluster_id,
            routing_core=cluster.core_sequence,
        )
        for read in cluster.reads:
            state.add_sequence(read.sequence, read.count)
        state.validate()
        return state
