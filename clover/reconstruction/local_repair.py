"""Conservative cluster-level repair for short repeat-rich regions."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Mapping

from .alignment_models import AlignmentResult


Window = tuple[int, int]
Repair = tuple[int, int, str]


@dataclass(frozen=True)
class LocalRepairConfig:
    """Data-size-independent controls for local repeat repair."""

    flank: int = 2
    max_window_length: int | None = None
    homopolymer_min_length: int = 4
    tandem_periods: tuple[int, ...] = (2, 3, 4)
    tandem_min_copies: int = 3
    minimum_top_fraction: float = 0.35
    minimum_margin_fraction: float = 0.10
    minimum_position_fraction: float = 0.60
    minimum_indel_support_weight: int | None = None

    def __post_init__(self) -> None:
        if self.flank < 0:
            raise ValueError("flank must be non-negative")
        if self.max_window_length is not None and self.max_window_length < 1:
            raise ValueError("max_window_length must be positive or None")
        if self.homopolymer_min_length < 2:
            raise ValueError("homopolymer_min_length must be at least 2")
        if not self.tandem_periods or any(
            period < 1 for period in self.tandem_periods
        ):
            raise ValueError(
                "tandem_periods must contain positive integers"
            )
        if self.tandem_min_copies < 2:
            raise ValueError("tandem_min_copies must be at least 2")
        for name, value in (
            ("minimum_top_fraction", self.minimum_top_fraction),
            ("minimum_margin_fraction", self.minimum_margin_fraction),
            ("minimum_position_fraction", self.minimum_position_fraction),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        if (
            self.minimum_indel_support_weight is not None
            and self.minimum_indel_support_weight < 1
        ):
            raise ValueError(
                "minimum_indel_support_weight must be positive or None"
            )


def _merge_intervals(intervals: Iterable[Window]) -> list[Window]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start >= merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def detect_repeat_windows(
    core: str,
    *,
    flank: int = 2,
    max_window_length: int | None = None,
    homopolymer_min_length: int = 4,
    tandem_periods: tuple[int, ...] = (2, 3, 4),
    tandem_min_copies: int = 3,
) -> list[Window]:
    """Find repeat windows without assuming a fixed strand length."""
    config = LocalRepairConfig(
        flank=flank,
        max_window_length=max_window_length,
        homopolymer_min_length=homopolymer_min_length,
        tandem_periods=tandem_periods,
        tandem_min_copies=tandem_min_copies,
    )

    intervals: list[Window] = []
    index = 0
    while index < len(core):
        end = index + 1
        while end < len(core) and core[end] == core[index]:
            end += 1
        if end - index >= config.homopolymer_min_length:
            intervals.append((index, end))
        index = end

    for period in config.tandem_periods:
        minimum = config.tandem_min_copies * period
        for start in range(len(core) - minimum + 1):
            motif = core[start:start + period]
            if (
                core[start:start + minimum]
                != motif * config.tandem_min_copies
            ):
                continue
            end = start + minimum
            while (
                end < len(core)
                and core[end] == motif[(end - start) % period]
            ):
                end += 1
            intervals.append((start, end))

    expanded = [
        (
            max(0, start - config.flank),
            min(len(core), end + config.flank),
        )
        for start, end in _merge_intervals(intervals)
    ]
    windows = _merge_intervals(expanded)
    if config.max_window_length is None:
        return windows
    return [
        window
        for window in windows
        if window[1] - window[0] <= config.max_window_length
    ]


def core_projection(core: str) -> dict:
    return {
        "bases": tuple(core),
        "insertions": tuple("" for _ in range(len(core) + 1)),
        "indels": tuple(False for _ in core),
    }


def project_read_alignment(alignment: AlignmentResult) -> dict:
    """Project one raw core-to-read alignment for reuse by all windows."""
    core_length = len(alignment.reference)
    bases: list[str | None] = [None] * core_length
    insertion_parts: list[list[str]] = [[] for _ in range(core_length + 1)]
    indels = [False] * core_length
    position = 0

    for column in alignment.iter_columns():
        if column.operation == "I":
            insertion_parts[position].append(column.query_base)
            if position < core_length:
                indels[position] = True
            continue

        if position >= core_length:
            raise RuntimeError("alignment consumes too many core bases")
        if column.operation == "D":
            indels[position] = True
        else:
            bases[position] = column.query_base
        position += 1

    if position != core_length:
        raise RuntimeError("alignment does not consume the complete core")

    return {
        "bases": tuple(bases),
        "insertions": tuple("".join(parts) for parts in insertion_parts),
        "indels": tuple(indels),
    }


def extract_local_fragment(
    projection: Mapping[str, tuple],
    start: int,
    end: int,
) -> tuple[str, bool]:
    """Extract slots start..end-1 and bases start..end-1."""
    bases = projection["bases"]
    insertions = projection["insertions"]
    indels = projection["indels"]
    if not 0 <= start < end <= len(bases):
        raise ValueError("local window is outside the core")

    parts: list[str] = []
    for position in range(start, end):
        parts.append(insertions[position])
        if bases[position] is not None:
            parts.append(bases[position])
    return "".join(parts), any(indels[start:end])


def collect_window_votes(
    sequence_weights: Mapping[str, int],
    projections: Mapping[str, Mapping[str, tuple]],
    window: Window,
) -> tuple[Counter[str], int]:
    """Collect exact local fragments plus total weighted indel evidence."""
    votes: Counter[str] = Counter()
    indel_support_weight = 0
    start, end = window

    for sequence, weight in sequence_weights.items():
        fragment, has_indel = extract_local_fragment(
            projections[sequence], start, end
        )
        votes[fragment] += weight
        if has_indel:
            indel_support_weight += weight
    return votes, indel_support_weight


def _length_votes(votes: Mapping[str, int]) -> Counter[int]:
    grouped: Counter[int] = Counter()
    for fragment, weight in votes.items():
        grouped[len(fragment)] += weight
    return grouped


def consensus_for_length_group(
    votes: Mapping[str, int],
    length: int,
    *,
    minimum_position_fraction: float = 0.60,
) -> str | None:
    """Build a weighted positional consensus within one fragment length."""
    members = [(fragment, weight) for fragment, weight in votes.items()
               if len(fragment) == length]
    group_weight = sum(weight for _, weight in members)
    if group_weight < 1:
        return None

    consensus: list[str] = []
    for position in range(length):
        symbols: Counter[str] = Counter()
        for fragment, weight in members:
            symbols[fragment[position]] += weight
        symbol, support = min(
            symbols.items(),
            key=lambda item: (-item[1], item[0]),
        )
        if support / group_weight < minimum_position_fraction:
            return None
        consensus.append(symbol)
    return "".join(consensus)


def select_length_consensus(
    votes: Mapping[str, int],
    *,
    core_window_length: int,
    total_weight: int,
    indel_support_weight: int,
    minimum_indel_support_weight: int | None = None,
    minimum_top_fraction: float = 0.35,
    minimum_margin_fraction: float = 0.10,
    minimum_position_fraction: float = 0.60,
) -> str | None:
    """
    Select a dominant local length state, then consensus within that state.

    This merges fragments that support the same local indel length even when
    substitutions make their complete strings different.
    """
    if not votes or total_weight < 1:
        return None
    if indel_support_weight <= 0:
        return None
    if (
        minimum_indel_support_weight is not None
        and indel_support_weight < minimum_indel_support_weight
    ):
        return None

    length_votes = _length_votes(votes)
    ranked = sorted(
        length_votes,
        key=lambda length: (
            -length_votes[length],
            abs(length - core_window_length),
            length,
        ),
    )
    winner_length = ranked[0]
    top_weight = length_votes[winner_length]
    second_weight = max(
        (length_votes[length] for length in ranked[1:]),
        default=0,
    )
    required_margin = max(
        1,
        math.ceil(minimum_margin_fraction * total_weight),
    )

    if top_weight / total_weight < minimum_top_fraction:
        return None
    if top_weight - second_weight < required_margin:
        return None

    return consensus_for_length_group(
        votes,
        winner_length,
        minimum_position_fraction=minimum_position_fraction,
    )


def baseline_local_fragment(
    base_decisions: tuple[str, ...],
    insertion_decisions: tuple[str, ...],
    start: int,
    end: int,
) -> str:
    parts: list[str] = []
    for position in range(start, end):
        parts.append(insertion_decisions[position])
        if base_decisions[position] != "D":
            parts.append(base_decisions[position])
    return "".join(parts)


def choose_local_repairs(
    core: str,
    windows: Iterable[Window],
    sequence_weights: Mapping[str, int],
    projections: Mapping[str, Mapping[str, tuple]],
    base_decisions: tuple[str, ...],
    insertion_decisions: tuple[str, ...],
    *,
    config: LocalRepairConfig | None = None,
) -> list[Repair]:
    """Choose conservative length-stratified repairs for repeat windows."""
    config = config or LocalRepairConfig()
    total_weight = sum(sequence_weights.values())
    repairs: list[Repair] = []

    for start, end in windows:
        votes, indel_weight = collect_window_votes(
            sequence_weights, projections, (start, end)
        )
        winner = select_length_consensus(
            votes,
            core_window_length=end - start,
            total_weight=total_weight,
            indel_support_weight=indel_weight,
            minimum_indel_support_weight=(
                config.minimum_indel_support_weight
            ),
            minimum_top_fraction=config.minimum_top_fraction,
            minimum_margin_fraction=config.minimum_margin_fraction,
            minimum_position_fraction=config.minimum_position_fraction,
        )
        baseline = baseline_local_fragment(
            base_decisions, insertion_decisions, start, end
        )
        if winner is not None and winner != baseline:
            repairs.append((start, end, winner))
    return repairs


def render_consensus(
    base_decisions: tuple[str, ...],
    insertion_decisions: tuple[str, ...],
    repairs: Iterable[Repair] = (),
) -> str:
    repairs_by_start = {
        start: (end, fragment) for start, end, fragment in repairs
    }
    parts: list[str] = []
    position = 0

    while position < len(base_decisions):
        repair = repairs_by_start.get(position)
        if repair is not None:
            end, fragment = repair
            parts.append(fragment)
            position = end
            continue

        parts.append(insertion_decisions[position])
        if base_decisions[position] != "D":
            parts.append(base_decisions[position])
        position += 1

    parts.append(insertion_decisions[len(base_decisions)])
    return "".join(parts)
