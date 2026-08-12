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



def align_global_edlib(
    reference: str,
    query: str,
) -> AlignmentResult:
    """
    Exact unit-cost global alignment using Edlib.

    Canonical operation semantics:

        = : match
        X : substitution
        I : insertion in query relative to reference
        D : deletion in query relative to reference
    """
    _validate_dna_sequence(reference, "reference")
    _validate_dna_sequence(query, "query")

    try:
        import edlib
    except ImportError as exc:
        raise RuntimeError(
            "Edlib backend requested but Python package 'edlib' "
            "is not installed"
        ) from exc

    raw_result = edlib.align(
        query,
        reference,
        mode="NW",
        task="path",
    )

    edit_distance = raw_result["editDistance"]

    if edit_distance < 0:
        raise RuntimeError(
            "Edlib failed to produce a global alignment"
        )

    # Edlib may return an empty CIGAR when either sequence is empty.
    # Handle these global-alignment boundary cases explicitly so that
    # canonical I/D semantics remain backend-independent.
    if reference == "" and query == "":
        return AlignmentResult.from_gapped_alignment(
            reference="",
            query="",
            aligned_reference="",
            aligned_query="",
            backend="edlib",
            backend_score=0,
            backend_score_name="unit_edit_distance",
        )

    if reference == "":
        return AlignmentResult.from_gapped_alignment(
            reference=reference,
            query=query,
            aligned_reference="-" * len(query),
            aligned_query=query,
            backend="edlib",
            backend_score=len(query),
            backend_score_name="unit_edit_distance",
        )

    if query == "":
        return AlignmentResult.from_gapped_alignment(
            reference=reference,
            query=query,
            aligned_reference=reference,
            aligned_query="-" * len(reference),
            backend="edlib",
            backend_score=len(reference),
            backend_score_name="unit_edit_distance",
        )

    nice = edlib.getNiceAlignment(
        raw_result,
        query,
        reference,
    )

    result = AlignmentResult.from_gapped_alignment(
        reference=reference,
        query=query,
        aligned_reference=nice["target_aligned"],
        aligned_query=nice["query_aligned"],
        backend="edlib",
        backend_score=edit_distance,
        backend_score_name="unit_edit_distance",
    )

    if result.edit_distance != edit_distance:
        raise RuntimeError(
            "Edlib alignment path is inconsistent with "
            "Edlib edit distance"
        )

    return result

def _align_global_existing(
    reference: str,
    query: str,
    *,
    backend: str = "nw",
) -> AlignmentResult:
    """
    Unified global-alignment entry point.
    """
    if backend == "nw":
        return align_global_nw(
            reference,
            query,
        )

    if backend == "edlib":
        return align_global_edlib(
            reference,
            query,
        )

    raise ValueError(
        f"unsupported global alignment backend: {backend}"
    )
# ---- Phase 11 WFA backend -------------------------------------------------
_WFA_ALIGNER = None


def _get_wfa_aligner():
    global _WFA_ALIGNER
    if _WFA_ALIGNER is None:
        try:
            from pywfa import WavefrontAligner
        except ImportError as exc:
            raise ImportError(
                "backend='wfa' requires pywfa; install pywfa==0.5.1"
            ) from exc
        _WFA_ALIGNER = WavefrontAligner(
            distance="affine",
            span="end-to-end",
            scope="full",
            heuristic=None,
            mismatch=1,
            gap_opening=0,
            gap_extension=1,
        )
    return _WFA_ALIGNER


def _validate_wfa_sequence(sequence: str, name: str) -> None:
    if not isinstance(sequence, str):
        raise TypeError(f"{name} must be a string")
    invalid = set(sequence) - set("ACGT")
    if invalid:
        raise ValueError(
            f"{name} contains invalid DNA symbols: {sorted(invalid)}"
        )


def _wfa_gapped_alignment(reference: str, query: str, cigar_tuples):
    ref_parts = []
    query_parts = []
    ref_index = 0
    query_index = 0

    for operation, length in cigar_tuples:
        if operation in {0, 7, 8}:  # M, =, X
            ref_end = ref_index + length
            query_end = query_index + length
            ref_parts.append(reference[ref_index:ref_end])
            query_parts.append(query[query_index:query_end])
            ref_index = ref_end
            query_index = query_end
        elif operation == 1:  # I: query only
            query_end = query_index + length
            ref_parts.append("-" * length)
            query_parts.append(query[query_index:query_end])
            query_index = query_end
        elif operation == 2:  # D: reference only
            ref_end = ref_index + length
            ref_parts.append(reference[ref_index:ref_end])
            query_parts.append("-" * length)
            ref_index = ref_end
        else:
            raise RuntimeError(
                f"unsupported WFA CIGAR operation code: {operation}"
            )

    if ref_index != len(reference):
        raise RuntimeError("WFA CIGAR did not consume the complete reference")
    if query_index != len(query):
        raise RuntimeError("WFA CIGAR did not consume the complete query")

    aligned_reference = "".join(ref_parts)
    aligned_query = "".join(query_parts)
    if len(aligned_reference) != len(aligned_query):
        raise RuntimeError("WFA produced unequal gapped alignment lengths")
    return aligned_reference, aligned_query


def _align_wfa(reference: str, query: str) -> AlignmentResult:
    """Exact end-to-end WFA under the same unit-edit objective as Edlib."""
    _validate_wfa_sequence(reference, "reference")
    _validate_wfa_sequence(query, "query")
    result = _get_wfa_aligner()(query, reference)
    if result.status != 0:
        raise RuntimeError(f"WFA alignment failed with status {result.status}")

    aligned_reference, aligned_query = _wfa_gapped_alignment(
        reference, query, result.cigartuples
    )
    distance = -int(result.score)
    if distance < 0:
        raise RuntimeError(f"unexpected positive WFA score: {result.score}")

    alignment = AlignmentResult.from_gapped_alignment(
        reference=reference,
        query=query,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        backend="wfa",
        backend_score=distance,
        backend_score_name="unit_edit_distance",
    )
    if alignment.edit_distance != distance:
        raise RuntimeError(
            "WFA score and canonical edit distance disagree: "
            f"{distance} != {alignment.edit_distance}"
        )
    return alignment


def align_global(reference: str, query: str, *, backend: str = "nw") -> AlignmentResult:
    if backend == "wfa":
        return _align_wfa(reference, query)
    return _align_global_existing(reference, query, backend=backend)
# ---- End Phase 11 WFA backend --------------------------------------------
