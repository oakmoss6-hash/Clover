"""Pairwise global-alignment backends used by Clover reconstruction."""

from __future__ import annotations

from dataclasses import dataclass


DNA_ALPHABET = frozenset("ACGTN")
WFA_DNA_ALPHABET = frozenset("ACGT")
SUPPORTED_BACKENDS = ("wfa", "edlib", "nw")


def _validate_sequence(sequence: str, name: str, alphabet=DNA_ALPHABET) -> None:
    if not isinstance(sequence, str):
        raise TypeError(f"{name} must be a string")
    invalid = set(sequence) - set(alphabet)
    if invalid:
        raise ValueError(
            f"{name} contains invalid DNA symbols: {''.join(sorted(invalid))}"
        )


@dataclass(frozen=True)
class GlobalAlignment:
    """Backend-independent gapped global alignment."""

    reference: str
    query: str
    aligned_reference: str
    aligned_query: str
    edit_distance: int
    backend: str

    def __post_init__(self) -> None:
        if len(self.aligned_reference) != len(self.aligned_query):
            raise ValueError("aligned sequences must have equal length")
        if self.aligned_reference.replace("-", "") != self.reference:
            raise ValueError("aligned_reference does not reconstruct reference")
        if self.aligned_query.replace("-", "") != self.query:
            raise ValueError("aligned_query does not reconstruct query")

        distance = 0
        for ref_base, query_base in zip(
            self.aligned_reference,
            self.aligned_query,
        ):
            if ref_base == "-" and query_base == "-":
                raise ValueError("alignment column cannot contain two gaps")
            if ref_base != query_base:
                distance += 1

        if distance != self.edit_distance:
            raise ValueError(
                "edit_distance is inconsistent with the gapped alignment"
            )


def _make_alignment(
    reference: str,
    query: str,
    aligned_reference: str,
    aligned_query: str,
    *,
    backend: str,
) -> GlobalAlignment:
    distance = sum(
        ref_base != query_base
        for ref_base, query_base in zip(aligned_reference, aligned_query)
    )
    return GlobalAlignment(
        reference=reference,
        query=query,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        edit_distance=distance,
        backend=backend,
    )


def _align_nw(reference: str, query: str) -> GlobalAlignment:
    """Reference Needleman-Wunsch backend with unit edit costs."""
    _validate_sequence(reference, "reference")
    _validate_sequence(query, "query")

    m = len(reference)
    n = len(query)
    diagonal, delete, insert = 0, 1, 2

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[diagonal] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = i
        trace[i][0] = delete
    for j in range(1, n + 1):
        dp[0][j] = j
        trace[0][j] = insert

    for i in range(1, m + 1):
        ref_base = reference[i - 1]
        for j in range(1, n + 1):
            query_base = query[j - 1]
            diag_cost = dp[i - 1][j - 1] + (ref_base != query_base)
            del_cost = dp[i - 1][j] + 1
            ins_cost = dp[i][j - 1] + 1
            best = min(diag_cost, del_cost, ins_cost)
            dp[i][j] = best
            if diag_cost == best:
                trace[i][j] = diagonal
            elif del_cost == best:
                trace[i][j] = delete
            else:
                trace[i][j] = insert

    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    i, j = m, n

    while i > 0 or j > 0:
        if i > 0 and j > 0 and trace[i][j] == diagonal:
            aligned_reference.append(reference[i - 1])
            aligned_query.append(query[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or trace[i][j] == delete):
            aligned_reference.append(reference[i - 1])
            aligned_query.append("-")
            i -= 1
        else:
            aligned_reference.append("-")
            aligned_query.append(query[j - 1])
            j -= 1

    aligned_reference.reverse()
    aligned_query.reverse()

    result = _make_alignment(
        reference,
        query,
        "".join(aligned_reference),
        "".join(aligned_query),
        backend="nw",
    )
    if result.edit_distance != dp[m][n]:
        raise RuntimeError("NW traceback is inconsistent with DP optimum")
    return result


def _align_edlib(reference: str, query: str) -> GlobalAlignment:
    """Exact unit-cost global alignment through Edlib."""
    _validate_sequence(reference, "reference")
    _validate_sequence(query, "query")

    try:
        import edlib
    except ImportError as exc:
        raise RuntimeError(
            "backend='edlib' requires the optional Python package 'edlib'"
        ) from exc

    if not reference and not query:
        return _make_alignment("", "", "", "", backend="edlib")
    if not reference:
        return _make_alignment(
            reference,
            query,
            "-" * len(query),
            query,
            backend="edlib",
        )
    if not query:
        return _make_alignment(
            reference,
            query,
            reference,
            "-" * len(reference),
            backend="edlib",
        )

    raw = edlib.align(query, reference, mode="NW", task="path")
    if raw["editDistance"] < 0:
        raise RuntimeError("Edlib failed to produce a global alignment")

    nice = edlib.getNiceAlignment(raw, query, reference)
    result = _make_alignment(
        reference,
        query,
        nice["target_aligned"],
        nice["query_aligned"],
        backend="edlib",
    )
    if result.edit_distance != raw["editDistance"]:
        raise RuntimeError("Edlib path and edit distance disagree")
    return result


_WFA_ALIGNER = None


def _get_wfa_aligner():
    global _WFA_ALIGNER
    if _WFA_ALIGNER is None:
        try:
            from pywfa import WavefrontAligner
        except ImportError as exc:
            raise RuntimeError(
                "backend='wfa' requires pywfa==0.5.1"
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


def _wfa_to_gapped(reference: str, query: str, cigar_tuples):
    ref_parts: list[str] = []
    query_parts: list[str] = []
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
        elif operation == 1:  # insertion in query
            query_end = query_index + length
            ref_parts.append("-" * length)
            query_parts.append(query[query_index:query_end])
            query_index = query_end
        elif operation == 2:  # deletion in query
            ref_end = ref_index + length
            ref_parts.append(reference[ref_index:ref_end])
            query_parts.append("-" * length)
            ref_index = ref_end
        else:
            raise RuntimeError(
                f"unsupported WFA CIGAR operation code: {operation}"
            )

    if ref_index != len(reference) or query_index != len(query):
        raise RuntimeError("WFA CIGAR did not consume both complete sequences")

    return "".join(ref_parts), "".join(query_parts)


def _align_wfa(reference: str, query: str) -> GlobalAlignment:
    """Exact end-to-end WFA under unit mismatch/indel cost."""
    _validate_sequence(reference, "reference", WFA_DNA_ALPHABET)
    _validate_sequence(query, "query", WFA_DNA_ALPHABET)

    raw = _get_wfa_aligner()(query, reference)
    if raw.status != 0:
        raise RuntimeError(f"WFA alignment failed with status {raw.status}")

    aligned_reference, aligned_query = _wfa_to_gapped(
        reference,
        query,
        raw.cigartuples,
    )
    result = _make_alignment(
        reference,
        query,
        aligned_reference,
        aligned_query,
        backend="wfa",
    )

    distance = -int(raw.score)
    if result.edit_distance != distance:
        raise RuntimeError(
            "WFA score and canonical edit distance disagree: "
            f"{distance} != {result.edit_distance}"
        )
    return result


def align_global(
    reference: str,
    query: str,
    *,
    backend: str = "wfa",
) -> GlobalAlignment:
    """Run one exact end-to-end pairwise global alignment."""
    backend = backend.lower()
    if backend == "wfa":
        return _align_wfa(reference, query)
    if backend == "edlib":
        return _align_edlib(reference, query)
    if backend == "nw":
        return _align_nw(reference, query)
    raise ValueError(
        f"unsupported reconstruction backend: {backend!r}; "
        f"expected one of {', '.join(SUPPORTED_BACKENDS)}"
    )
