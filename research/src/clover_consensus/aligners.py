from __future__ import annotations

from .alignment_models import AlignmentResult


DNA_ALPHABET = frozenset("ACGTN")


def _validate_dna_sequence(sequence: str, name: str) -> None:
    if not isinstance(sequence, str):
        raise TypeError(f"{name} must be a string")

    invalid = set(sequence) - DNA_ALPHABET

    if invalid:
        invalid_text = "".join(sorted(invalid))
        raise ValueError(
            f"{name} contains invalid DNA symbols: {invalid_text}"
        )


def align_global_nw(
    reference: str,
    query: str,
) -> AlignmentResult:
    """
    Exact pairwise global alignment using a Needleman-Wunsch-style
    dynamic programming recurrence with unit edit costs.

    Objective
    ---------
    Minimize:

        substitutions + insertions + deletions

    Costs
    -----
    match        : 0
    substitution : 1
    insertion    : 1
    deletion     : 1

    Complexity
    ----------
    Time:  O(m * n)
    Space: O(m * n)

    where m=len(reference), n=len(query).

    This implementation is intentionally simple and deterministic.
    It is a correctness/reference backend, not the final high-speed
    production aligner.

    Tie breaking
    ------------
    If multiple optimal traceback paths exist:

        diagonal > deletion > insertion

    This guarantees deterministic output. More sophisticated indel
    normalization will be handled later as a separate layer rather
    than being hidden inside this reference aligner.
    """
    _validate_dna_sequence(reference, "reference")
    _validate_dna_sequence(query, "query")

    m = len(reference)
    n = len(query)

    # Traceback codes.
    DIAGONAL = 0
    DELETE = 1
    INSERT = 2

    # dp[i][j] =
    # minimum edit distance between reference[:i] and query[:j]
    dp = [
        [0] * (n + 1)
        for _ in range(m + 1)
    ]

    trace = [
        [DIAGONAL] * (n + 1)
        for _ in range(m + 1)
    ]

    for i in range(1, m + 1):
        dp[i][0] = i
        trace[i][0] = DELETE

    for j in range(1, n + 1):
        dp[0][j] = j
        trace[0][j] = INSERT

    for i in range(1, m + 1):
        ref_base = reference[i - 1]

        for j in range(1, n + 1):
            query_base = query[j - 1]

            substitution_cost = (
                0 if ref_base == query_base else 1
            )

            diagonal_cost = (
                dp[i - 1][j - 1]
                + substitution_cost
            )

            deletion_cost = (
                dp[i - 1][j]
                + 1
            )

            insertion_cost = (
                dp[i][j - 1]
                + 1
            )

            best_cost = min(
                diagonal_cost,
                deletion_cost,
                insertion_cost,
            )

            dp[i][j] = best_cost

            # Deterministic tie-breaking:
            # diagonal > deletion > insertion
            if diagonal_cost == best_cost:
                trace[i][j] = DIAGONAL
            elif deletion_cost == best_cost:
                trace[i][j] = DELETE
            else:
                trace[i][j] = INSERT

    aligned_reference: list[str] = []
    aligned_query: list[str] = []

    i = m
    j = n

    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == DIAGONAL:
            aligned_reference.append(
                reference[i - 1]
            )
            aligned_query.append(
                query[j - 1]
            )
            i -= 1
            j -= 1

        elif i > 0 and (
            j == 0
            or trace[i][j] == DELETE
        ):
            aligned_reference.append(
                reference[i - 1]
            )
            aligned_query.append("-")
            i -= 1

        else:
            aligned_reference.append("-")
            aligned_query.append(
                query[j - 1]
            )
            j -= 1

    aligned_reference.reverse()
    aligned_query.reverse()

    aligned_reference_text = "".join(
        aligned_reference
    )

    aligned_query_text = "".join(
        aligned_query
    )

    result = AlignmentResult.from_gapped_alignment(
        reference=reference,
        query=query,
        aligned_reference=aligned_reference_text,
        aligned_query=aligned_query_text,
        backend="nw",
        backend_score=dp[m][n],
        backend_score_name="unit_edit_distance",
    )

    # The DP optimum and canonical result must agree.
    if result.edit_distance != dp[m][n]:
        raise RuntimeError(
            "NW traceback is inconsistent with DP optimum"
        )

    return result


def align_global(
    reference: str,
    query: str,
    *,
    backend: str = "nw",
) -> AlignmentResult:
    """
    Unified global-alignment entry point.

    More backends will be registered here in Phase 5.3.
    """
    if backend == "nw":
        return align_global_nw(
            reference,
            query,
        )

    raise ValueError(
        f"unsupported global alignment backend: {backend}"
    )
