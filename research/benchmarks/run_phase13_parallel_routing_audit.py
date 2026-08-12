#!/usr/bin/env python3

import argparse
import json
import math
import multiprocessing as mp
import queue
import sys
import time
import traceback
from collections import Counter
from pathlib import Path


BASE_TO_INT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Phase 13.2P parallel Clover routing audit"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-reads", type=int, default=None)
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        choices=[1, 4, 16, 64],
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5000,
    )
    return parser.parse_args()


def iter_fastq(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                return

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise RuntimeError("truncated FASTQ record")

            if not header.startswith("@"):
                raise RuntimeError(
                    f"invalid FASTQ header: {header[:80]!r}"
                )

            read_id = header[1:].strip().split()[0]
            yield read_id, sequence.strip().upper()


def prefix_length_for_workers(workers):
    if workers == 1:
        return 0

    k = round(math.log(workers, 4))

    if 4 ** k != workers:
        raise ValueError(
            "workers must be 1, 4, 16, or 64"
        )

    return k


def owner_index(sequence, prefix_length):
    """
    Clover-aware deterministic base-4 prefix partition.

    A=0, C=1, G=2, T=3.

    Invalid/short prefixes are sent to worker 0 and counted separately.
    """
    if prefix_length == 0:
        return 0, True

    if len(sequence) < prefix_length:
        return 0, False

    owner = 0

    for base in sequence[:prefix_length]:
        value = BASE_TO_INT.get(base)

        if value is None:
            return 0, False

        owner = owner * 4 + value

    return owner, True


def joint_key(tree_kind, horizontal_drifts, query_shift):
    return f"{tree_kind}|h={horizontal_drifts}|q={query_shift}"


def counter_to_json(counter):
    return {
        str(key): int(value)
        for key, value in sorted(
            counter.items(),
            key=lambda x: str(x[0]),
        )
    }


def joint_counter_to_json(counter):
    return {
        joint_key(tree, h, q): int(value)
        for (tree, h, q), value in sorted(
            counter.items(),
            key=lambda x: (
                str(x[0][0]),
                int(x[0][1]),
                int(x[0][2]),
            ),
        )
    }


def safe_put(q, value, error_event):
    while True:
        if error_event.is_set():
            raise RuntimeError(
                "a worker failed while producer was dispatching reads"
            )

        try:
            q.put(value, timeout=1.0)
            return
        except queue.Full:
            continue


def worker_main(
    worker_id,
    input_queue,
    result_queue,
    error_event,
):
    try:
        # Clover parses sys.argv while loading configuration.
        saved_argv = sys.argv[:]
        sys.argv = [saved_argv[0]]

        try:
            from clover import main as clover_main

            process = clover_main.MyProcess(
                f"phase13-audit-worker-{worker_id}",
                [],
                [],
            )
        finally:
            sys.argv = saved_argv

        process.capture_routing_hints = True

        raw_tree = Counter()
        raw_horizontal = Counter()
        raw_middle_shift = Counter()
        raw_joint = Counter()

        processed_reads = 0
        eligible_reads = 0
        routed_reads = 0
        new_core_reads = 0

        # Only O(number_of_clusters), not O(number_of_reads).
        core_sequence_by_id = {}

        last_route = {
            "core": None,
        }

        original_record = process._record_routing_hint

        def audited_record(
            core_index,
            sequence,
            tree_kind,
            horizontal_drifts,
            query_shift=0,
        ):
            nonlocal routed_reads

            routed_reads += 1
            last_route["core"] = core_index

            raw_tree[tree_kind] += 1
            raw_horizontal[horizontal_drifts] += 1
            raw_joint[
                (
                    tree_kind,
                    horizontal_drifts,
                    query_shift,
                )
            ] += 1

            if tree_kind.startswith("middle_"):
                raw_middle_shift[query_shift] += 1

            return original_record(
                core_index,
                sequence,
                tree_kind,
                horizontal_drifts,
                query_shift=query_shift,
            )

        process._record_routing_hint = audited_record

        started = time.perf_counter()

        while True:
            batch = input_queue.get()

            if batch is None:
                break

            for read_id, sequence in batch:
                processed_reads += 1

                if len(sequence) >= process.config_dict["read_len_min"]:
                    eligible_reads += 1

                last_route["core"] = None

                before_core_count = len(process.ref_dict)

                process.cluster(
                    f"{read_id} {sequence}"
                )

                after_core_count = len(process.ref_dict)

                if after_core_count > before_core_count:
                    if after_core_count != before_core_count + 1:
                        raise RuntimeError(
                            "one read created more than one Clover core"
                        )

                    core_id = next(
                        reversed(process.ref_dict)
                    )

                    core_sequence_by_id[core_id] = sequence
                    new_core_reads += 1

                    # Audit-only memory compaction:
                    # Clover Virtual_mode stores all tags in ref_dict,
                    # but old tags do not participate in later tree lookup.
                    if process.ref_dict.get(core_id):
                        process.ref_dict[core_id] = [
                            process.ref_dict[core_id][0]
                        ]

                route_core = last_route["core"]

                if route_core is not None:
                    if process.ref_dict.get(route_core):
                        process.ref_dict[route_core] = [
                            process.ref_dict[route_core][0]
                        ]

        elapsed = time.perf_counter() - started

        unique_tree = Counter()
        unique_horizontal = Counter()
        unique_middle_shift = Counter()
        unique_joint = Counter()

        routed_unique_all = 0
        noncore_unique_with_hint = 0
        clustered_unique_sequences = 0

        for core_id in process.ref_dict:
            cluster_hints = process.routing_hints.get(
                core_id,
                {},
            )

            core_sequence = core_sequence_by_id.get(
                core_id
            )

            # Distinct sequences in this Clover cluster:
            #
            # hints contains all unique routed sequences.
            # The newly-created core itself has no hint unless the same
            # sequence is observed again.
            cluster_unique = len(cluster_hints)

            if core_sequence not in cluster_hints:
                cluster_unique += 1

            clustered_unique_sequences += cluster_unique

            for sequence, hint in cluster_hints.items():
                routed_unique_all += 1

                # Exact copies of the routing core do not require a
                # non-core pairwise alignment after deduplication.
                if sequence == core_sequence:
                    continue

                noncore_unique_with_hint += 1

                unique_tree[hint.tree_kind] += 1
                unique_horizontal[
                    hint.horizontal_drifts
                ] += 1

                unique_joint[
                    (
                        hint.tree_kind,
                        hint.horizontal_drifts,
                        hint.query_shift,
                    )
                ] += 1

                if hint.tree_kind.startswith("middle_"):
                    unique_middle_shift[
                        hint.query_shift
                    ] += 1

        cluster_count = len(process.ref_dict)

        noncore_unique_sequences = max(
            0,
            clustered_unique_sequences - cluster_count,
        )

        hint_unique_coverage = (
            noncore_unique_with_hint
            / noncore_unique_sequences
            if noncore_unique_sequences
            else 1.0
        )

        front_back_h0 = sum(
            count
            for (tree, h, q), count in unique_joint.items()
            if tree in {"front", "back"} and h == 0
        )

        middle_h0_q0 = sum(
            count
            for (tree, h, q), count in unique_joint.items()
            if (
                tree.startswith("middle_")
                and h == 0
                and q == 0
            )
        )

        result = {
            "ok": True,
            "worker_id": worker_id,
            "processed_reads": processed_reads,
            "eligible_reads": eligible_reads,
            "routed_reads": routed_reads,
            "new_core_reads": new_core_reads,
            "unassigned_reads": (
                processed_reads
                - routed_reads
                - new_core_reads
            ),
            "cluster_count": cluster_count,
            "clustered_raw_reads": (
                routed_reads + new_core_reads
            ),
            "clustered_unique_sequences": (
                clustered_unique_sequences
            ),
            "routed_unique_all": routed_unique_all,
            "noncore_unique_sequences": (
                noncore_unique_sequences
            ),
            "noncore_unique_with_hint": (
                noncore_unique_with_hint
            ),
            "hint_unique_coverage": (
                hint_unique_coverage
            ),
            "elapsed_seconds": elapsed,
            "raw_tree": raw_tree,
            "raw_horizontal": raw_horizontal,
            "raw_middle_shift": raw_middle_shift,
            "raw_joint": raw_joint,
            "unique_tree": unique_tree,
            "unique_horizontal": unique_horizontal,
            "unique_middle_shift": unique_middle_shift,
            "unique_joint": unique_joint,
            "front_back_h0": front_back_h0,
            "middle_h0_q0": middle_h0_q0,
            "effective_config": {
                "read_len": process.read_len,
                "read_len_min": (
                    process.config_dict["read_len_min"]
                ),
                "end_tree_len": process.dna_tree_nums,
                "other_tree_len": process.fuzz_list[2],
                "other_tree_nums": (
                    process.config_dict["other_tree_nums"]
                ),
                "horizontal_drift_threshold": (
                    process.fuzz_tree_nums
                ),
                "vertical_drift_values": list(
                    process.loc_nums
                ),
                "tree_threshold": (
                    process.config_dict["tree_threshold"]
                ),
                "now_clust_threshold": (
                    process.now_clust_threshold
                ),
            },
        }

        result_queue.put(result)

    except Exception:
        error_event.set()

        result_queue.put(
            {
                "ok": False,
                "worker_id": worker_id,
                "traceback": traceback.format_exc(),
            }
        )


def merge_results(
    args,
    worker_results,
    wall_seconds,
    invalid_prefix_reads,
):
    raw_tree = Counter()
    raw_horizontal = Counter()
    raw_middle_shift = Counter()
    raw_joint = Counter()

    unique_tree = Counter()
    unique_horizontal = Counter()
    unique_middle_shift = Counter()
    unique_joint = Counter()

    scalar_fields = [
        "processed_reads",
        "eligible_reads",
        "routed_reads",
        "new_core_reads",
        "unassigned_reads",
        "cluster_count",
        "clustered_raw_reads",
        "clustered_unique_sequences",
        "routed_unique_all",
        "noncore_unique_sequences",
        "noncore_unique_with_hint",
        "front_back_h0",
        "middle_h0_q0",
    ]

    totals = {
        field: 0
        for field in scalar_fields
    }

    for result in worker_results:
        for field in scalar_fields:
            totals[field] += result[field]

        raw_tree.update(result["raw_tree"])
        raw_horizontal.update(
            result["raw_horizontal"]
        )
        raw_middle_shift.update(
            result["raw_middle_shift"]
        )
        raw_joint.update(result["raw_joint"])

        unique_tree.update(result["unique_tree"])
        unique_horizontal.update(
            result["unique_horizontal"]
        )
        unique_middle_shift.update(
            result["unique_middle_shift"]
        )
        unique_joint.update(result["unique_joint"])

    noncore_u = totals["noncore_unique_sequences"]

    hint_coverage = (
        totals["noncore_unique_with_hint"]
        / noncore_u
        if noncore_u
        else 1.0
    )

    M = totals["clustered_raw_reads"]
    U = totals["clustered_unique_sequences"]

    effective_configs = [
        result["effective_config"]
        for result in worker_results
    ]

    if effective_configs:
        first_config = effective_configs[0]

        for config in effective_configs[1:]:
            if config != first_config:
                raise RuntimeError(
                    "workers used different Clover configurations"
                )
    else:
        first_config = {}

    summary = {
        "phase": "13.2P/13.2C",
        "input": str(Path(args.input)),
        "max_reads": args.max_reads,
        "workers": args.workers,
        "prefix_length": (
            prefix_length_for_workers(args.workers)
        ),
        "batch_size": args.batch_size,
        "wall_seconds": wall_seconds,
        "reads_per_second": (
            totals["processed_reads"] / wall_seconds
            if wall_seconds
            else 0.0
        ),
        "invalid_prefix_reads": invalid_prefix_reads,
        **totals,
        "duplicate_compression_ratio_m_over_u": (
            M / U if U else 0.0
        ),
        "hint_unique_coverage": hint_coverage,
        "effective_clover_config": first_config,
        "raw": {
            "tree_kind": counter_to_json(raw_tree),
            "horizontal_drifts": counter_to_json(
                raw_horizontal
            ),
            "middle_query_shift": counter_to_json(
                raw_middle_shift
            ),
            "joint_route": joint_counter_to_json(
                raw_joint
            ),
        },
        "unique_noncore": {
            "tree_kind": counter_to_json(unique_tree),
            "horizontal_drifts": counter_to_json(
                unique_horizontal
            ),
            "middle_query_shift": counter_to_json(
                unique_middle_shift
            ),
            "joint_route": joint_counter_to_json(
                unique_joint
            ),
            "front_back_h0": totals["front_back_h0"],
            "middle_h0_q0": totals["middle_h0_q0"],
        },
        "workers_detail": [
            {
                "worker_id": r["worker_id"],
                "processed_reads": r["processed_reads"],
                "routed_reads": r["routed_reads"],
                "new_core_reads": r["new_core_reads"],
                "cluster_count": r["cluster_count"],
                "clustered_unique_sequences": (
                    r["clustered_unique_sequences"]
                ),
                "elapsed_seconds": r["elapsed_seconds"],
            }
            for r in sorted(
                worker_results,
                key=lambda x: x["worker_id"],
            )
        ],
    }

    return summary


def main():
    args = parse_args()

    prefix_length = prefix_length_for_workers(
        args.workers
    )

    ctx = mp.get_context("spawn")

    error_event = ctx.Event()
    result_queue = ctx.Queue()

    input_queues = [
        ctx.Queue(maxsize=4)
        for _ in range(args.workers)
    ]

    workers = []

    for worker_id in range(args.workers):
        process = ctx.Process(
            target=worker_main,
            args=(
                worker_id,
                input_queues[worker_id],
                result_queue,
                error_event,
            ),
            name=f"phase13-worker-{worker_id}",
        )

        process.start()
        workers.append(process)

    batches = [
        []
        for _ in range(args.workers)
    ]

    processed = 0
    invalid_prefix_reads = 0

    wall_started = time.perf_counter()

    try:
        for read_id, sequence in iter_fastq(args.input):
            if (
                args.max_reads is not None
                and processed >= args.max_reads
            ):
                break

            owner, valid_prefix = owner_index(
                sequence,
                prefix_length,
            )

            if not valid_prefix:
                invalid_prefix_reads += 1

            batches[owner].append(
                (read_id, sequence)
            )

            processed += 1

            if len(batches[owner]) >= args.batch_size:
                safe_put(
                    input_queues[owner],
                    batches[owner],
                    error_event,
                )
                batches[owner] = []

            if processed % 100000 == 0:
                elapsed = time.perf_counter() - wall_started

                print(
                    f"dispatched={processed:,} "
                    f"rate={processed / elapsed:,.0f} reads/s"
                )

        for worker_id, batch in enumerate(batches):
            if batch:
                safe_put(
                    input_queues[worker_id],
                    batch,
                    error_event,
                )

        for q in input_queues:
            safe_put(
                q,
                None,
                error_event,
            )

    except Exception:
        error_event.set()

        for process in workers:
            if process.is_alive():
                process.terminate()

        for process in workers:
            process.join()

        raise

    worker_results = []

    while len(worker_results) < args.workers:
        result = result_queue.get()

        if not result["ok"]:
            error_event.set()

            for process in workers:
                if process.is_alive():
                    process.terminate()

            for process in workers:
                process.join()

            raise RuntimeError(
                "worker failed:\n"
                + result["traceback"]
            )

        worker_results.append(result)

    for process in workers:
        process.join()

        if process.exitcode != 0:
            raise RuntimeError(
                f"{process.name} exited with "
                f"code {process.exitcode}"
            )

    wall_seconds = time.perf_counter() - wall_started

    summary = merge_results(
        args,
        worker_results,
        wall_seconds,
        invalid_prefix_reads,
    )

    output = Path(args.output)
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            summary,
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")

    print()
    print("=== Phase 13.2P parallel audit ===")
    print(
        "processed:",
        f"{summary['processed_reads']:,}",
    )
    print(
        "workers:",
        summary["workers"],
    )
    print(
        "wall:",
        f"{summary['wall_seconds']:.3f}s",
    )
    print(
        "throughput:",
        f"{summary['reads_per_second']:,.0f} reads/s",
    )
    print(
        "clusters:",
        f"{summary['cluster_count']:,}",
    )
    print(
        "M:",
        f"{summary['clustered_raw_reads']:,}",
    )
    print(
        "U:",
        f"{summary['clustered_unique_sequences']:,}",
    )
    print(
        "M/U:",
        f"{summary['duplicate_compression_ratio_m_over_u']:.4f}",
    )
    print(
        "hint coverage:",
        f"{summary['hint_unique_coverage']:.6f}",
    )
    print(
        "unique tree:",
        summary["unique_noncore"]["tree_kind"],
    )
    print(
        "unique horizontal:",
        summary["unique_noncore"]["horizontal_drifts"],
    )
    print(
        "unique middle shift:",
        summary["unique_noncore"]["middle_query_shift"],
    )
    print(
        "front/back h=0:",
        summary["unique_noncore"]["front_back_h0"],
    )
    print(
        "middle h=0 q=0:",
        summary["unique_noncore"]["middle_h0_q0"],
    )
    print(
        "invalid-prefix reads:",
        summary["invalid_prefix_reads"],
    )
    print(
        "output:",
        output,
    )


if __name__ == "__main__":
    main()
