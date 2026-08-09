"""Sequence-level evaluation metrics for CASPR benchmarks."""

from __future__ import annotations

from dataclasses import dataclass

from .aligners import align_global


@dataclass(frozen=True)
class SequenceMetrics:
    """Backend-independent metrics for one predicted sequence."""

    exact_match: bool
    edit_distance: int
    normalized_edit_distance: float
    substitutions: int
    insertions: int
    deletions: int
    reference_length: int
    prediction_length: int
    length_error: int


def evaluate_sequence(
    reference: str,
    prediction: str,
    *,
    backend: str = "edlib",
) -> SequenceMetrics:
    """Evaluate a prediction against a ground-truth reference."""
    alignment = align_global(
        reference,
        prediction,
        backend=backend,
    )

    return SequenceMetrics(
        exact_match=reference == prediction,
        edit_distance=alignment.edit_distance,
        normalized_edit_distance=(
            alignment.normalized_edit_distance
        ),
        substitutions=alignment.substitutions,
        insertions=alignment.insertions,
        deletions=alignment.deletions,
        reference_length=len(reference),
        prediction_length=len(prediction),
        length_error=len(prediction) - len(reference),
    )
