"""Clover-Aware Star Profile Reconstruction."""

from __future__ import annotations

from collections import defaultdict

from .aligners import align_global
from .alignment_models import AlignmentResult
from .models import ClusterRecord
from .profile_models import (
    BaseProfile,
    InsertionProfile,
    ReconstructionResult,
)


def _compress_sequences(
    cluster: ClusterRecord,
) -> dict[str, int]:
    """Aggregate read multiplicities by identical sequence."""
    sequence_weights: dict[str, int] = defaultdict(int)

    for read in cluster.reads:
        sequence_weights[read.sequence] += read.count

    return dict(sequence_weights)


def _project_alignment(
    alignment: AlignmentResult,
    weight: int,
    base_profiles: list[BaseProfile],
    insertion_profiles: list[InsertionProfile],
) -> None:
    """Project one read-to-backbone alignment into shared profiles."""
    backbone_position = 0
    insertion_parts: list[str] = []

    def flush_insertion() -> None:
        if not insertion_parts:
            return

        insertion_profiles[
            backbone_position
        ].add_vote(
            "".join(insertion_parts),
            weight,
        )

        insertion_parts.clear()

    for column in alignment.iter_columns():
        if column.operation == "I":
            insertion_parts.append(column.query_base)
            continue

        flush_insertion()

        if backbone_position >= len(base_profiles):
            raise ValueError(
                "alignment consumes more backbone bases "
                "than expected"
            )

        if column.ref_index != backbone_position:
            raise ValueError(
                "alignment reference coordinates are not contiguous"
            )

        if column.operation in {"=", "X"}:
            base_profiles[
                backbone_position
            ].add_vote(
                column.query_base,
                weight,
            )

        elif column.operation == "D":
            base_profiles[
                backbone_position
            ].add_vote(
                "D",
                weight,
            )

        else:
            raise ValueError(
                "unsupported canonical operation: "
                f"{column.operation}"
            )

        backbone_position += 1

    flush_insertion()

    if backbone_position != len(base_profiles):
        raise ValueError(
            "alignment does not cover every backbone position"
        )


def _build_consensus(
    base_profiles: list[BaseProfile],
    insertion_profiles: list[InsertionProfile],
    total_weight: int,
) -> str:
    consensus_parts: list[str] = [
        insertion_profiles[0].consensus_insertion(
            total_weight
        )
    ]

    for index, profile in enumerate(base_profiles):
        symbol = profile.consensus_symbol()

        if symbol != "D":
            consensus_parts.append(symbol)

        consensus_parts.append(
            insertion_profiles[
                index + 1
            ].consensus_insertion(
                total_weight
            )
        )

    return "".join(consensus_parts)


def reconstruct_cluster(
    cluster: ClusterRecord,
    *,
    backend: str = "nw",
) -> ReconstructionResult:
    """
    Reconstruct one cluster using a fixed-backbone star profile.

    The fixed backbone is cluster.core_sequence. Duplicate reads are
    compressed by identical sequence and their ReadRecord.count values
    are summed. Every unique non-backbone sequence is aligned against
    the backbone exactly once.

    No read-to-read comparison, all-pairs alignment, progressive
    merging, or U-by-U distance structure is created.

    Duplicate compression is linear in the total number of input
    sequence bases. Alignment cost is the sum of unique read-to-
    backbone alignment costs. Profile projection is linear in total
    alignment output size.

    With the reference Needleman-Wunsch backend and sequences of length
    approximately L, alignment costs approximately O(U * L**2). For
    fixed short DNA length L, scaling in unique sequence count U is
    linear.
    """
    backbone = cluster.core_sequence
    sequence_weights = _compress_sequences(cluster)

    if backbone not in sequence_weights:
        raise ValueError(
            "ClusterRecord violates the CASPR/Phase-4.5 contract: "
            "cluster.core_sequence must be represented in cluster.reads "
            "so the initial Clover core participates as evidence"
        )

    total_weight = sum(sequence_weights.values())

    base_profiles = [
        BaseProfile(
            backbone_index=index,
            backbone_base=base,
        )
        for index, base in enumerate(backbone)
    ]

    insertion_profiles = [
        InsertionProfile(slot=slot)
        for slot in range(len(backbone) + 1)
    ]

    pairwise_alignment_count = 0

    for sequence, weight in sequence_weights.items():
        if sequence == backbone:
            for index, base in enumerate(backbone):
                base_profiles[index].add_vote(
                    base,
                    weight,
                )

            continue

        alignment = align_global(
            backbone,
            sequence,
            backend=backend,
        )

        pairwise_alignment_count += 1

        _project_alignment(
            alignment,
            weight,
            base_profiles,
            insertion_profiles,
        )

    for profile in base_profiles:
        if profile.total_weight != total_weight:
            raise RuntimeError(
                "base profile weight conservation failed"
            )

    consensus = _build_consensus(
        base_profiles,
        insertion_profiles,
        total_weight,
    )

    return ReconstructionResult(
        cluster_id=cluster.cluster_id,
        backbone=backbone,
        consensus=consensus,
        backend=backend,
        raw_read_count=cluster.raw_read_count,
        total_weight=total_weight,
        unique_sequence_count=len(sequence_weights),
        pairwise_alignment_count=pairwise_alignment_count,
        base_profiles=tuple(base_profiles),
        insertion_profiles=tuple(insertion_profiles),
    )
