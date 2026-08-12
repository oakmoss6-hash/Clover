"""Truth-blind reconstruction-backbone selection for Clover clusters."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Literal


BackbonePolicy = Literal["core", "support_length", "max_span"]


def _validate_sequence_weights(sequence_weights: Mapping[str, int]) -> None:
    if not sequence_weights:
        raise ValueError("cannot select a backbone from an empty cluster")
    for sequence, weight in sequence_weights.items():
        if not sequence:
            raise ValueError("backbone candidates must be non-empty")
        if weight < 1:
            raise ValueError("sequence multiplicities must be positive")


def weighted_median_length(sequence_weights: Mapping[str, int]) -> float:
    """Return the read-count-weighted median observed sequence length."""
    _validate_sequence_weights(sequence_weights)
    length_weights: dict[int, int] = {}
    total_weight = 0
    for sequence, weight in sequence_weights.items():
        length = len(sequence)
        length_weights[length] = length_weights.get(length, 0) + weight
        total_weight += weight

    def observation_at(rank: int) -> int:
        cumulative = 0
        for length in sorted(length_weights):
            cumulative += length_weights[length]
            if cumulative >= rank:
                return length
        raise RuntimeError("weighted median length calculation failed")

    if total_weight % 2:
        return float(observation_at(total_weight // 2 + 1))
    lower = observation_at(total_weight // 2)
    upper = observation_at(total_weight // 2 + 1)
    return (lower + upper) / 2.0


def select_backbone(
    sequence_weights: Mapping[str, int],
    *,
    core_sequence: str,
    policy: BackbonePolicy | str = "core",
) -> str:
    """
    Select a reconstruction coordinate backbone without reference truth.

    ``core`` preserves the historical Clover/CASPR behavior.

    ``support_length`` first selects the most abundant observed sequence.
    If several sequences share the maximum multiplicity, it chooses the one
    whose length is closest to the read-count-weighted median cluster length.
    Remaining ties preserve first-observed mapping order.

    ``max_span`` is an experimental policy that selects the longest observed
    sequence, then higher multiplicity, then first-observed order.  It is
    retained for generality experiments and is not assumed to be the final
    production default.
    """
    _validate_sequence_weights(sequence_weights)
    if core_sequence not in sequence_weights:
        raise ValueError("core_sequence must be represented in cluster reads")

    if policy == "core":
        return core_sequence

    items = list(sequence_weights.items())

    if policy == "support_length":
        maximum_weight = max(weight for _, weight in items)
        candidates = [
            (index, sequence, weight)
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
        "unknown backbone policy: "
        f"{policy!r}; expected core, support_length, or max_span"
    )
