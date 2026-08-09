"""Backend-independent deterministic indel normalization."""

from __future__ import annotations

from .alignment_models import AlignmentResult


def _gap_kind_at(
    aligned_reference: list[str],
    aligned_query: list[str],
    index: int,
) -> str | None:
    if aligned_reference[index] == "-":
        return "I"

    if aligned_query[index] == "-":
        return "D"

    return None


def _is_exact_match(
    aligned_reference: list[str],
    aligned_query: list[str],
    index: int,
) -> bool:
    return (
        aligned_reference[index] != "-"
        and aligned_query[index] != "-"
        and aligned_reference[index]
        == aligned_query[index]
    )


def _maximum_left_shift(
    aligned_reference: list[str],
    aligned_query: list[str],
    start: int,
    end: int,
    gap_kind: str,
) -> int:
    shift = 0

    while start - shift > 0:
        left_index = start - shift - 1
        right_index = end - shift

        if not _is_exact_match(
            aligned_reference,
            aligned_query,
            left_index,
        ):
            break

        if gap_kind == "D":
            if (
                aligned_query[left_index]
                != aligned_reference[right_index]
            ):
                break

        elif gap_kind == "I":
            if (
                aligned_reference[left_index]
                != aligned_query[right_index]
            ):
                break

        else:
            raise ValueError(
                f"unsupported gap kind: {gap_kind}"
            )

        shift += 1

    return shift


def _apply_left_shift(
    aligned_reference: list[str],
    aligned_query: list[str],
    start: int,
    end: int,
    shift: int,
    gap_kind: str,
) -> int:
    if shift == 0:
        return start

    gap_length = end - start + 1
    segment_start = start - shift

    if gap_kind == "D":
        moved_query_bases = aligned_query[
            segment_start:start
        ]

        aligned_query[
            segment_start:end + 1
        ] = (
            ["-"] * gap_length
            + moved_query_bases
        )

    elif gap_kind == "I":
        moved_reference_bases = aligned_reference[
            segment_start:start
        ]

        aligned_reference[
            segment_start:end + 1
        ] = (
            ["-"] * gap_length
            + moved_reference_bases
        )

    else:
        raise ValueError(
            f"unsupported gap kind: {gap_kind}"
        )

    return segment_start


def _left_normalize_gap_blocks(
    aligned_reference: list[str],
    aligned_query: list[str],
) -> None:
    alignment_length = len(aligned_reference)
    index = 0

    while index < alignment_length:
        gap_kind = _gap_kind_at(
            aligned_reference,
            aligned_query,
            index,
        )

        if gap_kind is None:
            index += 1
            continue

        start = index
        end = start

        while (
            end + 1 < alignment_length
            and _gap_kind_at(
                aligned_reference,
                aligned_query,
                end + 1,
            )
            == gap_kind
        ):
            end += 1

        shift = _maximum_left_shift(
            aligned_reference,
            aligned_query,
            start,
            end,
            gap_kind,
        )

        new_start = _apply_left_shift(
            aligned_reference,
            aligned_query,
            start,
            end,
            shift,
            gap_kind,
        )

        if (
            new_start > 0
            and _gap_kind_at(
                aligned_reference,
                aligned_query,
                new_start - 1,
            )
            == gap_kind
        ):
            index = new_start - 1

            while (
                index > 0
                and _gap_kind_at(
                    aligned_reference,
                    aligned_query,
                    index - 1,
                )
                == gap_kind
            ):
                index -= 1

            continue

        index = end + 1


def normalize_indels_left(
    alignment: AlignmentResult,
) -> AlignmentResult:
    """
    Return a deterministic left-normalized alignment.

    Normalization operates only on the existing gapped alignment.
    It does not run dynamic programming or call an alignment backend.

    A gap block is shifted left only across exact-match columns and
    only when the shifted bases remain exact matches. Substitution
    columns are therefore never moved or altered.

    The common-case cost is linear in alignment length. Re-scanning
    occurs only when previously separate same-kind gap blocks merge.
    """
    if alignment.insertions == 0 and alignment.deletions == 0:
        return alignment

    aligned_reference = list(
        alignment.aligned_reference
    )
    aligned_query = list(
        alignment.aligned_query
    )

    _left_normalize_gap_blocks(
        aligned_reference,
        aligned_query,
    )

    normalized_reference = "".join(
        aligned_reference
    )
    normalized_query = "".join(
        aligned_query
    )

    if (
        normalized_reference
        == alignment.aligned_reference
        and normalized_query
        == alignment.aligned_query
    ):
        return alignment

    normalized = AlignmentResult.from_gapped_alignment(
        reference=alignment.reference,
        query=alignment.query,
        aligned_reference=normalized_reference,
        aligned_query=normalized_query,
        backend=alignment.backend,
        backend_score=alignment.backend_score,
        backend_score_name=(
            alignment.backend_score_name
        ),
    )

    if normalized.reference != alignment.reference:
        raise RuntimeError(
            "indel normalization changed the reference"
        )

    if normalized.query != alignment.query:
        raise RuntimeError(
            "indel normalization changed the query"
        )

    if normalized.edit_distance != alignment.edit_distance:
        raise RuntimeError(
            "indel normalization changed the edit distance"
        )

    if (
        normalized.substitutions
        != alignment.substitutions
    ):
        raise RuntimeError(
            "indel normalization changed substitution count"
        )

    if normalized.insertions != alignment.insertions:
        raise RuntimeError(
            "indel normalization changed insertion count"
        )

    if normalized.deletions != alignment.deletions:
        raise RuntimeError(
            "indel normalization changed deletion count"
        )

    return normalized
