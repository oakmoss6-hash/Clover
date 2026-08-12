"""Weighted profile models for CASPR reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field


BASE_SYMBOLS = ("A", "C", "G", "T", "N", "D")
BASE_SYMBOL_SET = frozenset(BASE_SYMBOLS)
DNA_SYMBOL_SET = frozenset("ACGTN")


@dataclass
class BaseProfile:
    """Weighted base and deletion votes at one backbone position."""

    backbone_index: int
    backbone_base: str
    counts: dict[str, int] = field(
        default_factory=lambda: {
            symbol: 0
            for symbol in BASE_SYMBOLS
        }
    )

    def __post_init__(self) -> None:
        if self.backbone_index < 0:
            raise ValueError("backbone_index must be >= 0")

        if self.backbone_base not in DNA_SYMBOL_SET:
            raise ValueError(
                "backbone_base must be one of A, C, G, T, N"
            )

        if set(self.counts) != BASE_SYMBOL_SET:
            raise ValueError(
                "counts must contain exactly A, C, G, T, N, D"
            )

        if any(
            not isinstance(weight, int) or weight < 0
            for weight in self.counts.values()
        ):
            raise ValueError(
                "base vote weights must be non-negative integers"
            )

    @property
    def total_weight(self) -> int:
        return sum(self.counts.values())

    def add_vote(
        self,
        symbol: str,
        weight: int,
    ) -> None:
        if symbol not in BASE_SYMBOL_SET:
            raise ValueError(
                f"unsupported base profile symbol: {symbol}"
            )

        if not isinstance(weight, int) or weight < 1:
            raise ValueError(
                "vote weight must be a positive integer"
            )

        self.counts[symbol] += weight

    def consensus_symbol(self) -> str:
        """Select a symbol using the deterministic CASPR tie rule."""
        maximum_weight = max(self.counts.values())

        tied_symbols = {
            symbol
            for symbol, weight in self.counts.items()
            if weight == maximum_weight
        }

        if self.backbone_base in tied_symbols:
            return self.backbone_base

        for symbol in BASE_SYMBOLS:
            if symbol in tied_symbols:
                return symbol

        raise RuntimeError("base profile has no consensus symbol")


@dataclass
class InsertionProfile:
    """Sparse complete-string insertion votes for one backbone slot."""

    slot: int
    non_empty_counts: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.slot < 0:
            raise ValueError("slot must be >= 0")

        for insertion, weight in self.non_empty_counts.items():
            self._validate_insertion(insertion)

            if not isinstance(weight, int) or weight < 0:
                raise ValueError(
                    "insertion vote weights must be "
                    "non-negative integers"
                )

    @staticmethod
    def _validate_insertion(insertion: str) -> None:
        if not insertion:
            raise ValueError(
                "non-empty insertion profiles cannot store "
                "an empty string"
            )

        invalid = set(insertion) - DNA_SYMBOL_SET

        if invalid:
            raise ValueError(
                "insertion contains invalid DNA symbols: "
                + "".join(sorted(invalid))
            )

    def add_vote(
        self,
        insertion: str,
        weight: int,
    ) -> None:
        self._validate_insertion(insertion)

        if not isinstance(weight, int) or weight < 1:
            raise ValueError(
                "vote weight must be a positive integer"
            )

        self.non_empty_counts[insertion] = (
            self.non_empty_counts.get(insertion, 0)
            + weight
        )

    def empty_weight(self, total_weight: int) -> int:
        if (
            not isinstance(total_weight, int)
            or total_weight < 0
        ):
            raise ValueError(
                "total_weight must be a non-negative integer"
            )

        non_empty_weight = sum(self.non_empty_counts.values())

        if non_empty_weight > total_weight:
            raise ValueError(
                "non-empty insertion weight exceeds "
                "total voting weight"
            )

        return total_weight - non_empty_weight

    def consensus_insertion(self, total_weight: int) -> str:
        """Select an insertion using deterministic CASPR tie rules."""
        candidates = {
            "": self.empty_weight(total_weight),
            **self.non_empty_counts,
        }

        maximum_weight = max(candidates.values())

        tied_insertions = [
            insertion
            for insertion, weight in candidates.items()
            if weight == maximum_weight
        ]

        return min(
            tied_insertions,
            key=lambda insertion: (
                insertion != "",
                len(insertion),
                insertion,
            ),
        )


@dataclass(frozen=True)
class ReconstructionResult:
    """Consensus and inspectable weighted profiles for one cluster."""

    cluster_id: str
    backbone: str
    consensus: str
    backend: str
    raw_read_count: int
    total_weight: int
    unique_sequence_count: int
    pairwise_alignment_count: int
    base_profiles: tuple[BaseProfile, ...]
    insertion_profiles: tuple[InsertionProfile, ...]

    @property
    def changed_base_count(self) -> int:
        return sum(
            profile.consensus_symbol()
            not in {profile.backbone_base, "D"}
            for profile in self.base_profiles
        )

    @property
    def deleted_base_count(self) -> int:
        return sum(
            profile.consensus_symbol() == "D"
            for profile in self.base_profiles
        )

    @property
    def inserted_base_count(self) -> int:
        return sum(
            len(
                profile.consensus_insertion(
                    self.total_weight
                )
            )
            for profile in self.insertion_profiles
        )
