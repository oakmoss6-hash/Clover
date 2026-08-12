"""Truth-blind reconstruction-backbone selection."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Literal


BackbonePolicy = Literal["core", "support_length", "max_span"]
SUPPORTED_BACKBONES = ("core", "support_length", "max_span")


def _validate(sequence_weights: Mapping[str, int]) -> None:
    if not sequence_weights:
        raise ValueError("cannot select a backbone from an empty cluster")
    for sequence, weight in sequence_weights.items():
        if not sequence:
            raise ValueError("backbone candidates must be non-empty")
        if weight < 1:
            raise ValueError("sequence multiplicities must be positive")


def weighted_median_length(sequence_weights: Mapping[str, int]) -> float:
    """Read-count-weighted median observed sequence length."""
    _validate(sequence_weights)
    length_weights: dict[int, int] = {}
    total = 0
    for sequence, weight in sequence_weights.items():
        length_weights[len(sequence)] = (
            length_weights.get(len(sequence), 0) + weight
        )
        total += weight

    def observation_at(rank: int) -> int:
        cumulative = 0
        for length in sorted(length_weights):
            cumulative += length_weights[length]
            if cumulative >= rank:
                return length
        raise RuntimeError("weighted median calculation failed")

    if total % 2:
        return float(observation_at(total // 2 + 1))
    return (
        observation_at(total // 2)
        + observation_at(total // 2 + 1)
    ) / 2.0


def select_backbone(
    sequence_weights: Mapping[str, int],
    *,
    core_sequence: str,
    policy: BackbonePolicy | str = "core",
) -> str:
    """Choose a backbone using only observed cluster evidence."""
    _validate(sequence_weights)
    if core_sequence not in sequence_weights:
        raise ValueError("core_sequence must be represented in cluster reads")

    if policy == "core":
        return core_sequence

    items = list(sequence_weights.items())

    if policy == "support_length":
        maximum_weight = max(weight for _, weight in items)
        candidates = [
            (index, sequence)
            for index, (sequence, weight) in enumerate(items)
            if weight == maximum_weight
        ]
        if len(candidates) == 1:
            return candidates[0][1]
        center = weighted_median_length(sequence_weights)
        return min(
            candidates,
            key=lambda item: (
                abs(len(item[1]) - center),
                item[0],
            ),
        )[1]

    if policy == "max_span":
        return min(
            enumerate(items),
            key=lambda item: (
                -len(item[1][0]),
                -item[1][1],
                item[0],
            ),
        )[1][0]

    raise ValueError(
        f"unknown backbone policy: {policy!r}; "
        f"expected one of {', '.join(SUPPORTED_BACKBONES)}"
    )
