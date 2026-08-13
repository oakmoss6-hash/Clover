"""Reed multi-read global reconstruction for Clover clusters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


DNA_ALPHABET = frozenset("ACGTN")
WFA_DNA_ALPHABET = frozenset("ACGT")
SUPPORTED_BACKENDS = ("wfa", "edlib", "nw")
SUPPORTED_BACKBONES = ("core", "support_length", "max_span")
BackbonePolicy = Literal["core", "support_length", "max_span"]

BASE_SYMBOLS = ("A", "C", "G", "T", "N", "D")
BASE_SYMBOL_SET = frozenset(BASE_SYMBOLS)
DNA_SYMBOL_SET = frozenset("ACGTN")


@dataclass
class ClusterState:
    """Unique sequence evidence and multiplicities for one Clover cluster."""

    cluster_id: str
    routing_core: str
    sequence_weights: dict[str, int] = field(default_factory=dict)
    raw_read_count: int = 0

    def add_sequence(self, sequence: str, count: int = 1) -> None:
        if not sequence:
            raise ValueError("sequence must be non-empty")
        if count < 1:
            raise ValueError("count must be positive")
        self.sequence_weights[sequence] = (
            self.sequence_weights.get(sequence, 0) + count
        )
        self.raw_read_count += count

    @property
    def unique_sequence_count(self) -> int:
        return len(self.sequence_weights)

    @property
    def total_weight(self) -> int:
        return sum(self.sequence_weights.values())

    def validate(self) -> None:
        if not self.cluster_id:
            raise ValueError("cluster_id must be non-empty")
        if not self.routing_core:
            raise ValueError("routing_core must be non-empty")
        if self.routing_core not in self.sequence_weights:
            raise ValueError(
                "routing_core must be represented in sequence evidence"
            )
        if self.raw_read_count != self.total_weight:
            raise ValueError(
                "raw_read_count must equal summed sequence multiplicity"
            )


@dataclass
class CloverClusterStateAdapter:
    """Observe already-decided Clover memberships without affecting routing."""

    states: dict[int, ClusterState] = field(default_factory=dict)

    def record_membership(
        self,
        *,
        core_index: int,
        sequence: str,
        read_id: str | None = None,
        is_new_core: bool = False,
    ) -> None:
        # read_id is deliberately accepted for evaluation hooks but is not
        # retained by the reconstruction algorithm.
        del read_id

        if is_new_core:
            if core_index in self.states:
                raise RuntimeError(
                    f"duplicate creation of Clover cluster {core_index}"
                )
            state = ClusterState(
                cluster_id=str(core_index),
                routing_core=sequence,
            )
            state.add_sequence(sequence)
            self.states[core_index] = state
            return

        state = self.states.get(core_index)
        if state is None:
            raise RuntimeError(
                "Clover membership arrived before its routing core: "
                f"cluster={core_index}"
            )
        state.add_sequence(sequence)

    def validate_all(self) -> None:
        for state in self.states.values():
            state.validate()

    @property
    def cluster_count(self) -> int:
        return len(self.states)


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


@dataclass(frozen=True)
class ReconstructionResult:
    cluster_id: str
    routing_core: str
    backbone: str
    consensus: str
    backend: str
    raw_read_count: int
    unique_sequence_count: int
    pairwise_alignment_count: int


@dataclass(frozen=True)
class ReconstructionSummary:
    cluster_count: int
    raw_read_count: int
    unique_sequence_count: int
    pairwise_alignment_count: int

    @property
    def expected_pairwise_alignment_count(self) -> int:
        return self.unique_sequence_count - self.cluster_count


def _validate_sequence(sequence: str, name: str, alphabet=DNA_ALPHABET) -> None:
    if not isinstance(sequence, str):
        raise TypeError(f"{name} must be a string")
    invalid = set(sequence) - set(alphabet)
    if invalid:
        raise ValueError(
            f"{name} contains invalid DNA symbols: {''.join(sorted(invalid))}"
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


def _validate_backbone_candidates(
    sequence_weights: Mapping[str, int],
) -> None:
    if not sequence_weights:
        raise ValueError("cannot select a backbone from an empty cluster")
    for sequence, weight in sequence_weights.items():
        if not sequence:
            raise ValueError("backbone candidates must be non-empty")
        if weight < 1:
            raise ValueError("sequence multiplicities must be positive")


def weighted_median_length(sequence_weights: Mapping[str, int]) -> float:
    """Read-count-weighted median observed sequence length."""
    _validate_backbone_candidates(sequence_weights)
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
    _validate_backbone_candidates(sequence_weights)
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


@dataclass
class _BaseProfile:
    backbone_base: str
    counts: dict[str, int] = field(
        default_factory=lambda: {symbol: 0 for symbol in BASE_SYMBOLS}
    )

    @property
    def total_weight(self) -> int:
        return sum(self.counts.values())

    def add(self, symbol: str, weight: int) -> None:
        if symbol not in BASE_SYMBOL_SET:
            raise ValueError(f"unsupported base profile symbol: {symbol}")
        self.counts[symbol] += weight

    def consensus(self) -> str:
        maximum = max(self.counts.values())
        tied = {
            symbol for symbol, weight in self.counts.items()
            if weight == maximum
        }
        if self.backbone_base in tied:
            return self.backbone_base
        return next(symbol for symbol in BASE_SYMBOLS if symbol in tied)


@dataclass
class _InsertionProfile:
    counts: dict[str, int] = field(default_factory=dict)

    def add(self, insertion: str, weight: int) -> None:
        if not insertion:
            raise ValueError("insertion vote must be non-empty")
        invalid = set(insertion) - DNA_SYMBOL_SET
        if invalid:
            raise ValueError(
                "insertion contains invalid DNA symbols: "
                + "".join(sorted(invalid))
            )
        self.counts[insertion] = self.counts.get(insertion, 0) + weight

    def consensus(self, total_weight: int) -> str:
        non_empty_weight = sum(self.counts.values())
        if non_empty_weight > total_weight:
            raise RuntimeError("insertion profile exceeds total voting weight")

        candidates = {
            "": total_weight - non_empty_weight,
            **self.counts,
        }
        maximum = max(candidates.values())
        tied = [
            insertion for insertion, weight in candidates.items()
            if weight == maximum
        ]
        return min(
            tied,
            key=lambda insertion: (
                insertion != "",
                len(insertion),
                insertion,
            ),
        )


def _project_alignment(
    alignment: GlobalAlignment,
    weight: int,
    base_profiles: list[_BaseProfile],
    insertion_profiles: list[_InsertionProfile],
) -> None:
    backbone_position = 0
    insertion_parts: list[str] = []

    def flush_insertion() -> None:
        if insertion_parts:
            insertion_profiles[backbone_position].add(
                "".join(insertion_parts),
                weight,
            )
            insertion_parts.clear()

    for ref_base, query_base in zip(
        alignment.aligned_reference,
        alignment.aligned_query,
    ):
        if ref_base == "-":
            insertion_parts.append(query_base)
            continue

        flush_insertion()
        if backbone_position >= len(base_profiles):
            raise RuntimeError("alignment consumed too many backbone bases")

        if query_base == "-":
            base_profiles[backbone_position].add("D", weight)
        else:
            base_profiles[backbone_position].add(query_base, weight)
        backbone_position += 1

    flush_insertion()
    if backbone_position != len(base_profiles):
        raise RuntimeError("alignment did not cover the complete backbone")


def _render_consensus(
    base_profiles: list[_BaseProfile],
    insertion_profiles: list[_InsertionProfile],
    total_weight: int,
) -> str:
    parts: list[str] = []
    for position, profile in enumerate(base_profiles):
        parts.append(insertion_profiles[position].consensus(total_weight))
        symbol = profile.consensus()
        if symbol != "D":
            parts.append(symbol)
    parts.append(insertion_profiles[len(base_profiles)].consensus(total_weight))
    return "".join(parts)


def reconstruct_cluster(
    state: ClusterState,
    *,
    backend: str = "wfa",
    backbone_policy: BackbonePolicy | str = "core",
) -> ReconstructionResult:
    """
    Reconstruct one Clover cluster from compact unique-sequence evidence.

    Every unique non-backbone sequence is globally aligned exactly once, so
    a cluster with U unique sequences performs exactly U-1 pairwise global
    alignments.
    """
    state.validate()
    backbone = select_backbone(
        state.sequence_weights,
        core_sequence=state.routing_core,
        policy=backbone_policy,
    )
    total_weight = state.total_weight

    base_profiles = [_BaseProfile(base) for base in backbone]
    insertion_profiles = [
        _InsertionProfile() for _ in range(len(backbone) + 1)
    ]

    pairwise_alignment_count = 0
    for sequence, weight in state.sequence_weights.items():
        if sequence == backbone:
            for index, base in enumerate(backbone):
                base_profiles[index].add(base, weight)
            continue

        alignment = align_global(backbone, sequence, backend=backend)
        pairwise_alignment_count += 1
        _project_alignment(
            alignment,
            weight,
            base_profiles,
            insertion_profiles,
        )

    expected = state.unique_sequence_count - 1
    if pairwise_alignment_count != expected:
        raise RuntimeError(
            "star-alignment invariant failed: "
            f"{pairwise_alignment_count} != U-1 {expected}"
        )

    for profile in base_profiles:
        if profile.total_weight != total_weight:
            raise RuntimeError("base profile weight conservation failed")

    return ReconstructionResult(
        cluster_id=state.cluster_id,
        routing_core=state.routing_core,
        backbone=backbone,
        consensus=_render_consensus(
            base_profiles,
            insertion_profiles,
            total_weight,
        ),
        backend=backend,
        raw_read_count=state.raw_read_count,
        unique_sequence_count=state.unique_sequence_count,
        pairwise_alignment_count=pairwise_alignment_count,
    )


reconstruct_state = reconstruct_cluster


@dataclass
class CloverWorkerReconstructor:
    """Own compact Reed evidence inside one Clover worker."""

    worker_name: str
    backend: str = "wfa"
    backbone_policy: str = "core"
    adapter: CloverClusterStateAdapter = field(
        default_factory=CloverClusterStateAdapter
    )

    @classmethod
    def from_config(
        cls,
        worker_name: str,
        config: Mapping[str, object],
    ) -> "CloverWorkerReconstructor":
        return cls(
            worker_name=worker_name,
            backend=str(config.get("reconstruct_backend", "wfa")),
            backbone_policy=str(
                config.get("reconstruct_backbone", "core")
            ),
        )

    def attach_to_process(self, process) -> None:
        if process.cluster_membership_observer is not None:
            raise RuntimeError(
                "Clover worker already has a cluster membership observer"
            )
        if process.worker_finalize_observer is not None:
            raise RuntimeError("Clover worker already has a finalizer observer")

        process.cluster_membership_observer = self.record_membership
        process.worker_finalize_observer = self.finalize

    def record_membership(
        self,
        *,
        core_index: int,
        sequence: str,
        read_id: str | None = None,
        is_new_core: bool = False,
    ) -> None:
        self.adapter.record_membership(
            core_index=core_index,
            sequence=sequence,
            read_id=read_id,
            is_new_core=is_new_core,
        )

    def finalize(self) -> dict:
        """Reconstruct all clusters owned by this worker."""
        self.adapter.validate_all()
        rows = []
        raw_reads = 0
        unique_sequences = 0
        pairwise_alignments = 0

        for core_index in sorted(self.adapter.states):
            state = self.adapter.states[core_index]
            raw_reads += state.raw_read_count
            unique_sequences += state.unique_sequence_count

            if state.unique_sequence_count == 1:
                rows.append(
                    (
                        core_index,
                        state.routing_core,
                        state.routing_core,
                        state.routing_core,
                        state.raw_read_count,
                        1,
                        0,
                    )
                )
                continue

            result = reconstruct_cluster(
                state,
                backend=self.backend,
                backbone_policy=self.backbone_policy,
            )
            pairwise_alignments += result.pairwise_alignment_count
            rows.append(
                (
                    core_index,
                    result.routing_core,
                    result.backbone,
                    result.consensus,
                    result.raw_read_count,
                    result.unique_sequence_count,
                    result.pairwise_alignment_count,
                )
            )

        expected = unique_sequences - self.adapter.cluster_count
        if pairwise_alignments != expected:
            raise RuntimeError(
                "worker star-alignment invariant failed: "
                f"{pairwise_alignments} != {expected}"
            )

        prefix = self.worker_name
        return {
            prefix + "reconstruction_results": rows,
            prefix + "reconstruction_cluster_count": self.adapter.cluster_count,
            prefix + "reconstruction_raw_read_count": raw_reads,
            prefix + "reconstruction_unique_sequence_count": unique_sequences,
            prefix + "reconstruction_pairwise_alignment_count": pairwise_alignments,
        }


def write_reconstruction_output(
    count_dict: Mapping[str, object],
    worker_names,
    output_path: str | Path,
) -> ReconstructionSummary:
    """Write compact cluster consensus TSV and verify A = U - C."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_clusters = 0
    total_raw = 0
    total_unique = 0
    total_pairwise = 0

    with output_path.open("w", encoding="utf-8", newline="") as output:
        output.write(
            "worker\tcluster_id\trouting_core\tbackbone\tconsensus\t"
            "raw_read_count\tunique_sequence_count\t"
            "pairwise_alignment_count\n"
        )

        for worker_name in worker_names:
            prefix = str(worker_name)
            rows = count_dict[prefix + "reconstruction_results"]
            total_clusters += int(
                count_dict[prefix + "reconstruction_cluster_count"]
            )
            total_raw += int(
                count_dict[prefix + "reconstruction_raw_read_count"]
            )
            total_unique += int(
                count_dict[prefix + "reconstruction_unique_sequence_count"]
            )
            total_pairwise += int(
                count_dict[prefix + "reconstruction_pairwise_alignment_count"]
            )

            for row in rows:
                output.write("\t".join(map(str, (prefix, *row))) + "\n")

    summary = ReconstructionSummary(
        cluster_count=total_clusters,
        raw_read_count=total_raw,
        unique_sequence_count=total_unique,
        pairwise_alignment_count=total_pairwise,
    )
    if (
        summary.pairwise_alignment_count
        != summary.expected_pairwise_alignment_count
    ):
        raise RuntimeError(
            "global reconstruction invariant failed: "
            f"A={summary.pairwise_alignment_count} != "
            f"U-C={summary.expected_pairwise_alignment_count}"
        )
    return summary


__all__ = [
    "ClusterState",
    "CloverWorkerReconstructor",
    "GlobalAlignment",
    "ReconstructionResult",
    "ReconstructionSummary",
    "align_global",
    "reconstruct_cluster",
    "reconstruct_state",
    "select_backbone",
    "weighted_median_length",
    "write_reconstruction_output",
]
