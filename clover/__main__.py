"""Command-line entry point for Clover."""

import multiprocessing as mp
import sys
import time

from clover.load_config import HELP_TEXT
from clover.main import run_cli
from clover.reed import CloverWorkerReconstructor, reconstruct_cluster


_PARALLEL_STATES = None
_PARALLEL_BACKEND = "wfa"
_PARALLEL_BACKBONE = "core"


def _parallel_reconstruct_one(core_index):
    state = _PARALLEL_STATES[core_index]
    result = reconstruct_cluster(
        state,
        backend=_PARALLEL_BACKEND,
        backbone_policy=_PARALLEL_BACKBONE,
    )
    return (
        core_index,
        result.routing_core,
        result.backbone,
        result.consensus,
        result.raw_read_count,
        result.unique_sequence_count,
        result.pairwise_alignment_count,
    )


def _make_parallel_finalize(worker_count):
    def finalize(self):
        self.adapter.validate_all()
        rows = []
        raw_reads = 0
        unique_sequences = 0
        pairwise_alignments = 0
        multi_ids = []

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
            else:
                multi_ids.append(core_index)

        # Larger clusters first gives better dynamic load balance.
        multi_ids.sort(
            key=lambda core_index: self.adapter.states[
                core_index
            ].unique_sequence_count,
            reverse=True,
        )

        global _PARALLEL_STATES
        global _PARALLEL_BACKEND
        global _PARALLEL_BACKBONE
        _PARALLEL_STATES = self.adapter.states
        _PARALLEL_BACKEND = self.backend
        _PARALLEL_BACKBONE = self.backbone_policy

        started = time.perf_counter()
        if worker_count == 1 or len(multi_ids) <= 1:
            for core_index in multi_ids:
                rows.append(_parallel_reconstruct_one(core_index))
        else:
            # Linux/Ubuntu: fork inherits the frozen cluster-state dictionary
            # copy-on-write, so tasks only send integer cluster ids.
            context = mp.get_context("fork")
            with context.Pool(processes=worker_count) as pool:
                rows.extend(
                    pool.imap_unordered(
                        _parallel_reconstruct_one,
                        multi_ids,
                        chunksize=16,
                    )
                )

        reconstruction_seconds = time.perf_counter() - started
        rows.sort(key=lambda row: row[0])
        pairwise_alignments = sum(int(row[6]) for row in rows)

        expected = unique_sequences - self.adapter.cluster_count
        if pairwise_alignments != expected:
            raise RuntimeError(
                "worker star-alignment invariant failed: "
                f"{pairwise_alignments} != {expected}"
            )

        print("Reconstruction workers:", worker_count)
        print("Reconstruction-only time:", reconstruction_seconds)

        prefix = self.worker_name
        return {
            prefix + "reconstruction_results": rows,
            prefix + "reconstruction_cluster_count": self.adapter.cluster_count,
            prefix + "reconstruction_raw_read_count": raw_reads,
            prefix + "reconstruction_unique_sequence_count": unique_sequences,
            prefix + "reconstruction_pairwise_alignment_count": pairwise_alignments,
        }

    return finalize


def _extract_reconstruction_workers(argv):
    workers = 1
    cleaned = []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg == "--reconstruct-workers":
            if index + 1 >= len(argv):
                raise ValueError("--reconstruct-workers requires an integer")
            workers = int(argv[index + 1])
            index += 2
            continue
        if arg.startswith("--reconstruct-workers="):
            workers = int(arg.split("=", 1)[1])
            index += 1
            continue
        cleaned.append(arg)
        index += 1

    if workers < 1:
        raise ValueError("--reconstruct-workers must be at least 1")
    return workers, cleaned


def _force_serial_clover(argv):
    cleaned = []
    index = 0
    while index < len(argv):
        arg = argv[index]
        if arg == "-P":
            if index + 1 >= len(argv):
                raise ValueError("-P requires a value")
            index += 2
            continue
        if arg.startswith("-P") and len(arg) > 2:
            index += 1
            continue
        cleaned.append(arg)
        index += 1
    cleaned.extend(["-P", "0"])
    return cleaned


def main():
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(HELP_TEXT)
        print(
            "  --reconstruct-workers N       Parallel Reed workers after serial Clover clustering."
        )
        return

    workers, argv = _extract_reconstruction_workers(sys.argv[1:])

    if "--reconstruct" in argv:
        # Reconstruction experiments must not change Clover's clustering by
        # changing the number of Reed workers.  Freeze Clover at one worker;
        # parallelism starts only after clustering has finished.
        argv = _force_serial_clover(argv)
        CloverWorkerReconstructor.finalize = _make_parallel_finalize(workers)
        print("Clover clustering workers: 1 (forced for reconstruction)")

    sys.argv = [sys.argv[0], *argv]
    run_cli()


if __name__ == "__main__":
    main()
