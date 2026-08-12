from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator


VALID_OPERATIONS = frozenset({"=", "X", "I", "D"})


def operations_to_cigar(operations: tuple[str, ...]) -> str:
    """
    Convert canonical per-column alignment operations to extended CIGAR.

    Canonical operation semantics:

        = : match
        X : substitution / mismatch
        I : insertion in query relative to reference
        D : deletion in query relative to reference
    """
    if not operations:
        return ""

    for op in operations:
        if op not in VALID_OPERATIONS:
            raise ValueError(f"unsupported alignment operation: {op}")

    parts: list[str] = []
    current = operations[0]
    count = 1

    for op in operations[1:]:
        if op == current:
            count += 1
        else:
            parts.append(f"{count}{current}")
            current = op
            count = 1

    parts.append(f"{count}{current}")
    return "".join(parts)


@dataclass(frozen=True)
class AlignmentColumn:
    """
    One column of a global alignment.

    ref_index/query_index use zero-based coordinates in the original
    ungapped sequences. A gap is represented by None.
    """

    ref_index: int | None
    query_index: int | None
    ref_base: str
    query_base: str
    operation: str


@dataclass(frozen=True)
class AlignmentResult:
    """
    Backend-independent representation of a pairwise global alignment.

    The canonical internal CIGAR uses extended operations:

        = : identical bases
        X : substitution
        I : query contains an inserted base relative to reference
        D : query lacks a reference base

    `edit_distance` is always the canonical unit-cost edit distance:

        substitutions + insertions + deletions

    `backend_score` is deliberately separate because different backends
    may use different scoring models.
    """

    reference: str
    query: str

    aligned_reference: str
    aligned_query: str

    operations: tuple[str, ...]
    cigar: str

    matches: int
    substitutions: int
    insertions: int
    deletions: int

    edit_distance: int
    normalized_edit_distance: float

    backend: str

    backend_score: float | int | None = None
    backend_score_name: str | None = None

    def __post_init__(self) -> None:
        if len(self.aligned_reference) != len(self.aligned_query):
            raise ValueError(
                "aligned_reference and aligned_query must have equal length"
            )

        if len(self.operations) != len(self.aligned_reference):
            raise ValueError(
                "operations length must equal alignment length"
            )

        if any(op not in VALID_OPERATIONS for op in self.operations):
            raise ValueError("operations contain unsupported symbols")

        if operations_to_cigar(self.operations) != self.cigar:
            raise ValueError("cigar does not match operations")

        expected_edit_distance = (
            self.substitutions
            + self.insertions
            + self.deletions
        )

        if self.edit_distance != expected_edit_distance:
            raise ValueError(
                "edit_distance must equal substitutions + "
                "insertions + deletions"
            )

        expected_alignment_length = (
            self.matches
            + self.substitutions
            + self.insertions
            + self.deletions
        )

        if expected_alignment_length != len(self.operations):
            raise ValueError(
                "operation counts do not equal alignment length"
            )

        denominator = max(
            len(self.reference),
            len(self.query),
        )

        expected_ned = (
            0.0
            if denominator == 0
            else self.edit_distance / denominator
        )

        if abs(self.normalized_edit_distance - expected_ned) > 1e-12:
            raise ValueError(
                "normalized_edit_distance is inconsistent with "
                "canonical edit distance"
            )

        reconstructed_reference = self.aligned_reference.replace("-", "")
        reconstructed_query = self.aligned_query.replace("-", "")

        if reconstructed_reference != self.reference:
            raise ValueError(
                "aligned_reference does not reconstruct reference"
            )

        if reconstructed_query != self.query:
            raise ValueError(
                "aligned_query does not reconstruct query"
            )

    @property
    def alignment_length(self) -> int:
        return len(self.operations)

    @property
    def identity(self) -> float:
        """
        Fraction of alignment columns that are exact matches.
        """
        if self.alignment_length == 0:
            return 1.0
        return self.matches / self.alignment_length

    def iter_columns(self) -> Iterator[AlignmentColumn]:
        """
        Iterate through aligned columns while exposing original
        zero-based sequence coordinates.

        This interface will later be used by the star-profile
        reconstruction algorithm.
        """
        ref_index = 0
        query_index = 0

        for ref_base, query_base, op in zip(
            self.aligned_reference,
            self.aligned_query,
            self.operations,
        ):
            if op in {"=", "X"}:
                yield AlignmentColumn(
                    ref_index=ref_index,
                    query_index=query_index,
                    ref_base=ref_base,
                    query_base=query_base,
                    operation=op,
                )
                ref_index += 1
                query_index += 1

            elif op == "I":
                yield AlignmentColumn(
                    ref_index=None,
                    query_index=query_index,
                    ref_base="-",
                    query_base=query_base,
                    operation=op,
                )
                query_index += 1

            elif op == "D":
                yield AlignmentColumn(
                    ref_index=ref_index,
                    query_index=None,
                    ref_base=ref_base,
                    query_base="-",
                    operation=op,
                )
                ref_index += 1

            else:
                raise RuntimeError(f"unexpected operation: {op}")

    @classmethod
    def from_gapped_alignment(
        cls,
        *,
        reference: str,
        query: str,
        aligned_reference: str,
        aligned_query: str,
        backend: str,
        backend_score: float | int | None = None,
        backend_score_name: str | None = None,
    ) -> "AlignmentResult":
        """
        Construct the canonical result from two gapped alignment strings.

        This constructor is intentionally backend-independent and will
        later be reused to normalize Edlib/WFA outputs.
        """
        if len(aligned_reference) != len(aligned_query):
            raise ValueError(
                "aligned sequences must have equal length"
            )

        operations: list[str] = []

        matches = 0
        substitutions = 0
        insertions = 0
        deletions = 0

        for ref_base, query_base in zip(
            aligned_reference,
            aligned_query,
        ):
            if ref_base == "-" and query_base == "-":
                raise ValueError(
                    "alignment column cannot contain two gaps"
                )

            if ref_base == "-":
                operations.append("I")
                insertions += 1

            elif query_base == "-":
                operations.append("D")
                deletions += 1

            elif ref_base == query_base:
                operations.append("=")
                matches += 1

            else:
                operations.append("X")
                substitutions += 1

        operation_tuple = tuple(operations)
        cigar = operations_to_cigar(operation_tuple)

        edit_distance = (
            substitutions
            + insertions
            + deletions
        )

        denominator = max(
            len(reference),
            len(query),
        )

        normalized_edit_distance = (
            0.0
            if denominator == 0
            else edit_distance / denominator
        )

        return cls(
            reference=reference,
            query=query,
            aligned_reference=aligned_reference,
            aligned_query=aligned_query,
            operations=operation_tuple,
            cigar=cigar,
            matches=matches,
            substitutions=substitutions,
            insertions=insertions,
            deletions=deletions,
            edit_distance=edit_distance,
            normalized_edit_distance=normalized_edit_distance,
            backend=backend,
            backend_score=backend_score,
            backend_score_name=backend_score_name,
        )
