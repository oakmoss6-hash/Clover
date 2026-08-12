"""Worker-local Clover reconstruction integration."""

from __future__ import annotations

from dataclasses import dataclass, field

from .clover_state_adapter import CloverClusterStateAdapter
from .reconstruct import reconstruct_state


@dataclass
class CloverWorkerReconstructor:
    """
    Own reconstruction evidence inside one Clover worker.

    Raw reads never leave the worker.  Only compact reconstruction results
    are returned when ``finalize()`` is called.
    """

    worker_name: str
    backend: str = "wfa"
    backbone_policy: str = "core"
    local_repeat_repair: bool = False

    adapter: CloverClusterStateAdapter = field(
        default_factory=CloverClusterStateAdapter
    )

    def attach_to_process(self, process) -> None:
        """Attach this reconstructor to one already-created Clover worker."""
        if process.cluster_membership_observer is not None:
            raise RuntimeError(
                "Clover worker already has a cluster membership observer"
            )

        if process.worker_finalize_observer is not None:
            raise RuntimeError(
                "Clover worker already has a finalizer observer"
            )

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
        # read_id is intentionally ignored by reconstruction itself.
        # It exists for passive downstream evaluation/auditing.
        self.adapter.record_membership(
            core_index=core_index,
            sequence=sequence,
            read_id=read_id,
            is_new_core=is_new_core,
        )

    def finalize(self) -> dict:
        """Reconstruct all locally owned Clover clusters."""

        self.adapter.validate_all()

        rows = []

        raw_reads = 0
        unique_sequences = 0
        pairwise_alignments = 0

        singleton_clusters = 0
        multi_unique_clusters = 0

        for core_index in sorted(self.adapter.states):
            state = self.adapter.states[core_index]

            # Exact singleton fast path.
            #
            # If U_c == 1, every read in the cluster is the same sequence
            # after exact duplicate aggregation.  Consensus and backbone are
            # therefore already known and no profile or pairwise alignment
            # needs to be constructed.
            if state.unique_sequence_count == 1:
                sequence = state.routing_core

                singleton_clusters += 1
                raw_reads += state.raw_read_count
                unique_sequences += 1

                rows.append(
                    (
                        core_index,
                        sequence,
                        sequence,
                        state.raw_read_count,
                        1,
                        0,
                    )
                )
                continue

            multi_unique_clusters += 1

            result = reconstruct_state(
                state,
                backend=self.backend,
                backbone_policy=self.backbone_policy,
                local_repeat_repair=self.local_repeat_repair,
            )

            raw_reads += result.raw_read_count
            unique_sequences += result.unique_sequence_count
            pairwise_alignments += result.pairwise_alignment_count

            # Tuple instead of full ReconstructionResult:
            # keep multiprocessing output compact.
            rows.append(
                (
                    core_index,
                    result.consensus,
                    result.backbone,
                    result.raw_read_count,
                    result.unique_sequence_count,
                    result.pairwise_alignment_count,
                )
            )

        # Global invariant over independently reconstructed clusters:
        # sum(U_c - 1) = total_U - number_of_clusters
        expected_pairwise = (
            unique_sequences - self.adapter.cluster_count
        )

        if pairwise_alignments != expected_pairwise:
            raise RuntimeError(
                "worker star-alignment invariant failed: "
                f"{pairwise_alignments} != "
                f"{unique_sequences} - {self.adapter.cluster_count}"
            )

        prefix = self.worker_name

        return {
            prefix + "reconstruction_results": rows,
            prefix + "reconstruction_cluster_count": (
                self.adapter.cluster_count
            ),
            prefix + "reconstruction_raw_read_count": raw_reads,
            prefix + "reconstruction_unique_sequence_count": (
                unique_sequences
            ),
            prefix + "reconstruction_pairwise_alignment_count": (
                pairwise_alignments
            ),
            prefix + "reconstruction_singleton_cluster_count": (
                singleton_clusters
            ),
            prefix + "reconstruction_multi_unique_cluster_count": (
                multi_unique_clusters
            ),
        }
