"""Clover-aware multi-read global reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field

from .alignment import GlobalAlignment, align_global
from .backbone import BackbonePolicy, select_backbone
from .state import ClusterState


BASE_SYMBOLS = ("A", "C", "G", "T", "N", "D")
BASE_SYMBOL_SET = frozenset(BASE_SYMBOLS)
DNA_SYMBOL_SET = frozenset("ACGTN")


@dataclass
class _BaseProfile:
    backbone_base: str
    counts: dict[str, int] = field(
        default_factory=lambda: {symbol: 0 for symbol in BASE_SYMBOLS}
    )

    @property
    def total_weight(self) -> int:
        return sum(self.counts.values())

    def add(self, symbol: str, weight: int) -> None:
        if symbol not in BASE_SYMBOL_SET:
            raise ValueError(f"unsupported base profile symbol: {symbol}")
        self.counts[symbol] += weight

    def consensus(self) -> str:
        maximum = max(self.counts.values())
        tied = {
            symbol for symbol, weight in self.counts.items()
            if weight == maximum
        }
        if self.backbone_base in tied:
            return self.backbone_base
        return next(symbol for symbol in BASE_SYMBOLS if symbol in tied)


@dataclass
class _InsertionProfile:
    counts: dict[str, int] = field(default_factory=dict)

    def add(self, insertion: str, weight: int) -> None:
        if not insertion:
            raise ValueError("insertion vote must be non-empty")
        invalid = set(insertion) - DNA_SYMBOL_SET
        if invalid:
            raise ValueError(
                "insertion contains invalid DNA symbols: "
                + "".join(sorted(invalid))
            )
        self.counts[insertion] = self.counts.get(insertion, 0) + weight

    def consensus(self, total_weight: int) -> str:
        non_empty_weight = sum(self.counts.values())
        if non_empty_weight > total_weight:
            raise RuntimeError("insertion profile exceeds total voting weight")

        candidates = {
            "": total_weight - non_empty_weight,
            **self.counts,
        }
        maximum = max(candidates.values())
        tied = [
            insertion for insertion, weight in candidates.items()
            if weight == maximum
        ]
        return min(
            tied,
            key=lambda insertion: (
                insertion != "",
                len(insertion),
                insertion,
            ),
        )


@dataclass(frozen=True)
class ReconstructionResult:
    cluster_id: str
    routing_core: str
    backbone: str
    consensus: str
    backend: str
    raw_read_count: int
    unique_sequence_count: int
    pairwise_alignment_count: int


def _project_alignment(
    alignment: GlobalAlignment,
    weight: int,
    base_profiles: list[_BaseProfile],
    insertion_profiles: list[_InsertionProfile],
) -> None:
    backbone_position = 0
    insertion_parts: list[str] = []

    def flush_insertion() -> None:
        if insertion_parts:
            insertion_profiles[backbone_position].add(
                "".join(insertion_parts),
                weight,
            )
            insertion_parts.clear()

    for ref_base, query_base in zip(
        alignment.aligned_reference,
        alignment.aligned_query,
    ):
        if ref_base == "-":
            insertion_parts.append(query_base)
            continue

        flush_insertion()
        if backbone_position >= len(base_profiles):
            raise RuntimeError("alignment consumed too many backbone bases")

        if query_base == "-":
            base_profiles[backbone_position].add("D", weight)
        else:
            base_profiles[backbone_position].add(query_base, weight)
        backbone_position += 1

    flush_insertion()
    if backbone_position != len(base_profiles):
        raise RuntimeError("alignment did not cover the complete backbone")


def _render_consensus(
    base_profiles: list[_BaseProfile],
    insertion_profiles: list[_InsertionProfile],
    total_weight: int,
) -> str:
    parts: list[str] = []
    for position, profile in enumerate(base_profiles):
        parts.append(insertion_profiles[position].consensus(total_weight))
        symbol = profile.consensus()
        if symbol != "D":
            parts.append(symbol)
    parts.append(insertion_profiles[len(base_profiles)].consensus(total_weight))
    return "".join(parts)


def reconstruct_state(
    state: ClusterState,
    *,
    backend: str = "wfa",
    backbone_policy: BackbonePolicy | str = "core",
) -> ReconstructionResult:
    """
    Reconstruct one Clover cluster from compact unique-sequence evidence.

    Every unique non-backbone sequence is globally aligned exactly once, so
    a cluster with U unique sequences performs exactly U-1 pairwise global
    alignments.
    """
    state.validate()
    backbone = select_backbone(
        state.sequence_weights,
        core_sequence=state.routing_core,
        policy=backbone_policy,
    )
    total_weight = state.total_weight

    base_profiles = [_BaseProfile(base) for base in backbone]
    insertion_profiles = [
        _InsertionProfile() for _ in range(len(backbone) + 1)
    ]

    pairwise_alignment_count = 0
    for sequence, weight in state.sequence_weights.items():
        if sequence == backbone:
            for index, base in enumerate(backbone):
                base_profiles[index].add(base, weight)
            continue

        alignment = align_global(backbone, sequence, backend=backend)
        pairwise_alignment_count += 1
        _project_alignment(
            alignment,
            weight,
            base_profiles,
            insertion_profiles,
        )

    expected = state.unique_sequence_count - 1
    if pairwise_alignment_count != expected:
        raise RuntimeError(
            "star-alignment invariant failed: "
            f"{pairwise_alignment_count} != U-1 {expected}"
        )

    for profile in base_profiles:
        if profile.total_weight != total_weight:
            raise RuntimeError("base profile weight conservation failed")

    return ReconstructionResult(
        cluster_id=state.cluster_id,
        routing_core=state.routing_core,
        backbone=backbone,
        consensus=_render_consensus(
            base_profiles,
            insertion_profiles,
            total_weight,
        ),
        backend=backend,
        raw_read_count=state.raw_read_count,
        unique_sequence_count=state.unique_sequence_count,
        pairwise_alignment_count=pairwise_alignment_count,
    )
