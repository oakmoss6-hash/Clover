"""Clover-aware sparse star-profile reconstruction."""

from __future__ import annotations

from .aligners import align_global
from .alignment_models import AlignmentResult
from .backbone import BackbonePolicy, select_backbone
from .cluster_state import ClusterState
from .local_repair import (
    LocalRepairConfig,
    choose_local_repairs,
    core_projection,
    detect_repeat_windows,
    project_read_alignment,
    render_consensus,
)
from .models import ClusterRecord
from .normalization import normalize_indels_left
from .profile_models import (
    BaseProfile,
    InsertionProfile,
    ReconstructionResult,
)


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
        insertion_profiles[backbone_position].add_vote(
            "".join(insertion_parts), weight
        )
        insertion_parts.clear()

    for column in alignment.iter_columns():
        if column.operation == "I":
            insertion_parts.append(column.query_base)
            continue

        flush_insertion()
        if backbone_position >= len(base_profiles):
            raise ValueError(
                "alignment consumes more backbone bases than expected"
            )
        if column.ref_index != backbone_position:
            raise ValueError(
                "alignment reference coordinates are not contiguous"
            )

        if column.operation in {"=", "X"}:
            base_profiles[backbone_position].add_vote(
                column.query_base, weight
            )
        elif column.operation == "D":
            base_profiles[backbone_position].add_vote("D", weight)
        else:
            raise ValueError(
                f"unsupported canonical operation: {column.operation}"
            )
        backbone_position += 1

    flush_insertion()
    if backbone_position != len(base_profiles):
        raise ValueError(
            "alignment does not cover every backbone position"
        )


def _profile_decisions(
    base_profiles: list[BaseProfile],
    insertion_profiles: list[InsertionProfile],
    total_weight: int,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    return (
        tuple(profile.consensus_symbol() for profile in base_profiles),
        tuple(
            profile.consensus_insertion(total_weight)
            for profile in insertion_profiles
        ),
    )


def _target_length_guard(
    baseline: str,
    repaired: str,
    target_length: int | None,
) -> str:
    """Reject a repair only when it moves farther from a known target length."""
    if target_length is None:
        return repaired
    if target_length < 1:
        raise ValueError("target_length must be positive")
    baseline_error = abs(len(baseline) - target_length)
    repaired_error = abs(len(repaired) - target_length)
    return baseline if repaired_error > baseline_error else repaired


def reconstruct_state(
    state: ClusterState,
    *,
    backend: str = "nw",
    normalize_indels: bool = False,
    local_repeat_repair: bool = False,
    target_length: int | None = None,
    local_repair_config: LocalRepairConfig | None = None,
    backbone_policy: BackbonePolicy | str = "core",
) -> ReconstructionResult:
    """
    Reconstruct one compact Clover cluster state.

    The coordinate backbone is selected without reference truth. Every unique
    non-backbone sequence is globally aligned exactly once, so the pairwise
    call count remains U-1 independently of raw cluster coverage.
    """
    if normalize_indels and local_repeat_repair:
        raise ValueError(
            "local repeat repair cannot be combined with indel normalization"
        )

    state.validate()
    sequence_weights = state.sequence_weights
    backbone = select_backbone(
        sequence_weights,
        core_sequence=state.routing_core,
        policy=backbone_policy,
    )
    total_weight = state.total_weight

    base_profiles = [
        BaseProfile(backbone_index=index, backbone_base=base)
        for index, base in enumerate(backbone)
    ]
    insertion_profiles = [
        InsertionProfile(slot=slot)
        for slot in range(len(backbone) + 1)
    ]

    repair_config = local_repair_config or LocalRepairConfig()
    windows = (
        detect_repeat_windows(
            backbone,
            flank=repair_config.flank,
            max_window_length=repair_config.max_window_length,
            homopolymer_min_length=(
                repair_config.homopolymer_min_length
            ),
            tandem_periods=repair_config.tandem_periods,
            tandem_min_copies=repair_config.tandem_min_copies,
        )
        if local_repeat_repair else []
    )
    local_projections = (
        {backbone: core_projection(backbone)}
        if windows else {}
    )
    pairwise_alignment_count = 0

    for sequence, weight in sequence_weights.items():
        if sequence == backbone:
            for index, base in enumerate(backbone):
                base_profiles[index].add_vote(base, weight)
            continue

        alignment = align_global(backbone, sequence, backend=backend)
        pairwise_alignment_count += 1

        if windows:
            local_projections[sequence] = project_read_alignment(alignment)
        if normalize_indels:
            alignment = normalize_indels_left(alignment)

        _project_alignment(
            alignment,
            weight,
            base_profiles,
            insertion_profiles,
        )

    expected_pairwise = state.unique_sequence_count - 1
    if pairwise_alignment_count != expected_pairwise:
        raise RuntimeError(
            "star-alignment invariant failed: "
            f"{pairwise_alignment_count} != U-1 {expected_pairwise}"
        )

    for profile in base_profiles:
        if profile.total_weight != total_weight:
            raise RuntimeError("base profile weight conservation failed")

    base_decisions, insertion_decisions = _profile_decisions(
        base_profiles, insertion_profiles, total_weight
    )
    repairs = (
        choose_local_repairs(
            backbone,
            windows,
            sequence_weights,
            local_projections,
            base_decisions,
            insertion_decisions,
            config=repair_config,
        )
        if windows else []
    )
    baseline_consensus = render_consensus(
        base_decisions, insertion_decisions
    )
    repaired_consensus = render_consensus(
        base_decisions, insertion_decisions, repairs
    )
    consensus = _target_length_guard(
        baseline_consensus,
        repaired_consensus,
        target_length,
    )

    return ReconstructionResult(
        cluster_id=state.cluster_id,
        backbone=backbone,
        consensus=consensus,
        backend=backend,
        raw_read_count=state.raw_read_count,
        total_weight=total_weight,
        unique_sequence_count=state.unique_sequence_count,
        pairwise_alignment_count=pairwise_alignment_count,
        base_profiles=tuple(base_profiles),
        insertion_profiles=tuple(insertion_profiles),
    )


def reconstruct_cluster(
    cluster: ClusterRecord,
    *,
    backend: str = "nw",
    normalize_indels: bool = False,
    local_repeat_repair: bool = False,
    target_length: int | None = None,
    local_repair_config: LocalRepairConfig | None = None,
    backbone_policy: BackbonePolicy | str = "core",
) -> ReconstructionResult:
    """
    Backward-compatible adapter from the existing Clover ``ClusterRecord``.

    Future Clover integration can aggregate reads directly into
    ``ClusterState`` and call ``reconstruct_state`` to avoid retaining and
    recompressing duplicate read objects.
    """
    return reconstruct_state(
        ClusterState.from_cluster(cluster),
        backend=backend,
        normalize_indels=normalize_indels,
        local_repeat_repair=local_repeat_repair,
        target_length=target_length,
        local_repair_config=local_repair_config,
        backbone_policy=backbone_policy,
    )
