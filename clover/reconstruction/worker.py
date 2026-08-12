"""Worker-local integration between Clover routing and reconstruction."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

from .reconstruct import reconstruct_state
from .state import CloverClusterStateAdapter


@dataclass(frozen=True)
class ReconstructionSummary:
    cluster_count: int
    raw_read_count: int
    unique_sequence_count: int
    pairwise_alignment_count: int

    @property
    def expected_pairwise_alignment_count(self) -> int:
        return self.unique_sequence_count - self.cluster_count


@dataclass
class CloverWorkerReconstructor:
    """Own compact reconstruction evidence inside one Clover worker."""

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

            result = reconstruct_state(
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
