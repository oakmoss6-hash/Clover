"""Core data models for cluster-level reconstruction inputs."""

from __future__ import annotations

from dataclasses import dataclass

ALLOWED_DNA_ALPHABET = frozenset("ACGTN")


def validate_dna_sequence(sequence: str, field_name: str = "sequence") -> None:
    if not isinstance(sequence, str) or not sequence:
        raise ValueError(f"{field_name} must be a non-empty string")
    invalid = sorted(set(sequence) - ALLOWED_DNA_ALPHABET)
    if invalid:
        raise ValueError(
            f"{field_name} contains invalid DNA symbols: {''.join(invalid)}; "
            f"allowed alphabet is {''.join(sorted(ALLOWED_DNA_ALPHABET))}"
        )


@dataclass(frozen=True)
class ReadRecord:
    read_id: str
    sequence: str
    tag: str | None = None
    quality: str | None = None
    count: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.read_id, str) or not self.read_id:
            raise ValueError("read_id must be a non-empty string")
        validate_dna_sequence(self.sequence)
        if self.count < 1:
            raise ValueError("count must be >= 1")


@dataclass
class ClusterRecord:
    cluster_id: str
    core_sequence: str
    reads: list[ReadRecord]

    def __post_init__(self) -> None:
        if not isinstance(self.cluster_id, str) or not self.cluster_id:
            raise ValueError("cluster_id must be a non-empty string")
        validate_dna_sequence(self.core_sequence, "core_sequence")
        if not self.reads:
            raise ValueError("ClusterRecord reads must be non-empty")
        for read in self.reads:
            if not isinstance(read, ReadRecord):
                raise ValueError("ClusterRecord reads must contain ReadRecord objects")

    @property
    def raw_read_count(self) -> int:
        return sum(read.count for read in self.reads)

    @property
    def unique_sequence_count(self) -> int:
        return len({read.sequence for read in self.reads})

    @property
    def sequence_lengths(self) -> tuple[int, ...]:
        return tuple(len(read.sequence) for read in self.reads)
