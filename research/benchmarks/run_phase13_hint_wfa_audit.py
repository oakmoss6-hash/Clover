#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
import queue
import random
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
    p = argparse.ArgumentParser(
        description="Phase 13.3 Clover RoutingHint vs exact WFA audit"
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--max-reads", type=int, default=200000)
    p.add_argument(
        "--workers",
        type=int,
        default=4,
        choices=[1, 4, 16, 64],
    )
    p.add_argument("--batch-size", type=int, default=5000)

    # Per-worker reservoir sizes.
    p.add_argument("--candidate-per-worker", type=int, default=2500)
    p.add_argument("--control-per-worker", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260812)
    return p.parse_args()


def iter_fastq(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        while True:
            header = f.readline()
            if not header:
                return

            sequence = f.readline()
            plus = f.readline()
            quality = f.readline()

            if not sequence or not plus or not quality:
                raise RuntimeError("truncated FASTQ")

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
        raise ValueError("workers must be 1, 4, 16, or 64")

    return k


def owner_index(sequence, prefix_length):
    if prefix_length == 0:
        return 0

    if len(sequence) < prefix_length:
        return 0

    owner = 0

    for base in sequence[:prefix_length]:
        value = BASE_TO_INT.get(base)

        if value is None:
            return 0

        owner = owner * 4 + value

    return owner


def safe_put(q, value, error_event):
    while True:
        if error_event.is_set():
            raise RuntimeError("worker failure during dispatch")

        try:
            q.put(value, timeout=1.0)
            return
        except queue.Full:
            pass


def reservoir_add(reservoir, item, seen, limit, rng):
    seen += 1

    if limit <= 0:
        return seen

    if len(reservoir) < limit:
        reservoir.append(item)
        return seen

    index = rng.randrange(seen)

    if index < limit:
        reservoir[index] = item

    return seen


def max_abs_diagonal(aligned_reference, aligned_query):
    """
    Maximum |query_position - reference_position| along the WFA path.
    """
    if len(aligned_reference) != len(aligned_query):
        raise RuntimeError("unequal gapped alignment lengths")

    ref_pos = 0
    query_pos = 0
    maximum = 0

    for ref_base, query_base in zip(
        aligned_reference,
        aligned_query,
    ):
        if ref_base != "-":
            ref_pos += 1

        if query_base != "-":
            query_pos += 1

        maximum = max(
            maximum,
            abs(query_pos - ref_pos),
        )

    return maximum


def analyse_pair(
    worker_id,
    class_name,
    core_id,
    reference,
    query,
    hint,
    align_global,
):
    alignment = align_global(
        reference,
        query,
        backend="wfa",
    )

    same_length = len(reference) == len(query)

    hamming = None

    if same_length:
        hamming = sum(
            a != b
            for a, b in zip(reference, query)
        )

    max_diagonal = max_abs_diagonal(
        alignment.aligned_reference,
        alignment.aligned_query,
    )

    no_indel = (
        alignment.insertions == 0
        and alignment.deletions == 0
    )

    return {
        "worker_id": worker_id,
        "class": class_name,
        "core_id": core_id,
        "tree_kind": hint.tree_kind,
        "horizontal_drifts": hint.horizontal_drifts,
        "query_shift": hint.query_shift,
        "reference_length": len(reference),
        "query_length": len(query),
        "length_difference": len(query) - len(reference),
        "edit_distance": alignment.edit_distance,
        "substitutions": alignment.substitutions,
        "insertions": alignment.insertions,
        "deletions": alignment.deletions,
        "has_indel": int(not no_indel),
        "same_length": int(same_length),
        "same_length_no_indel": int(
            same_length and no_indel
        ),
        "hamming_distance": (
            "" if hamming is None else hamming
        ),
        "hamming_equals_edit": int(
            hamming is not None
            and hamming == alignment.edit_distance
        ),
        "max_abs_diagonal": max_diagonal,
        "reference": reference,
        "query": query,
    }


def worker_main(
    worker_id,
    input_queue,
    result_queue,
    error_event,
    candidate_limit,
    control_limit,
    seed,
):
    try:
        research_root = (
            Path(__file__).resolve().parents[1]
        )
        research_src = research_root / "src"

        if str(research_src) not in sys.path:
            sys.path.insert(0, str(research_src))

        from clover_consensus.aligners import align_global

        saved_argv = sys.argv[:]
        sys.argv = [saved_argv[0]]

        try:
            from clover import main as clover_main

            process = clover_main.MyProcess(
                f"phase13-wfa-worker-{worker_id}",
                [],
                [],
            )
        finally:
            sys.argv = saved_argv

        process.capture_routing_hints = True

        core_sequences = {}
        last_route_core = {"value": None}

        original_record = process._record_routing_hint

        def record_hint(
            core_index,
            sequence,
            tree_kind,
            horizontal_drifts,
            query_shift=0,
        ):
            last_route_core["value"] = core_index

            return original_record(
                core_index,
                sequence,
                tree_kind,
                horizontal_drifts,
                query_shift=query_shift,
            )

        process._record_routing_hint = record_hint

        processed_reads = 0

        while True:
            batch = input_queue.get()

            if batch is None:
                break

            for read_id, sequence in batch:
                processed_reads += 1

                before_core_count = len(process.ref_dict)
                last_route_core["value"] = None

                process.cluster(
                    f"{read_id} {sequence}"
                )

                after_core_count = len(process.ref_dict)

                if after_core_count > before_core_count:
                    if after_core_count != before_core_count + 1:
                        raise RuntimeError(
                            "one read created multiple cores"
                        )

                    core_id = next(
                        reversed(process.ref_dict)
                    )

                    core_sequences[core_id] = sequence

                    # Audit-only memory compaction.
                    if process.ref_dict.get(core_id):
                        process.ref_dict[core_id] = [
                            process.ref_dict[core_id][0]
                        ]

                routed_core = last_route_core["value"]

                if routed_core is not None:
                    if process.ref_dict.get(routed_core):
                        process.ref_dict[routed_core] = [
                            process.ref_dict[routed_core][0]
                        ]

        rng_candidate = random.Random(
            seed + worker_id * 1009 + 1
        )
        rng_control = random.Random(
            seed + worker_id * 1009 + 2
        )

        candidate_sample = []
        control_sample = []

        candidate_seen = 0
        control_seen = 0
        exact_core_duplicates = 0

        for core_id, cluster_hints in process.routing_hints.items():
            reference = core_sequences.get(core_id)

            if reference is None:
                raise RuntimeError(
                    f"missing routing-core sequence for {core_id}"
                )

            for query, hint in cluster_hints.items():
                # Deduplication removes copies identical to the core,
                # therefore they are not U-1 alignment work.
                if query == reference:
                    exact_core_duplicates += 1
                    continue

                item = (
                    core_id,
                    reference,
                    query,
                    hint,
                )

                is_candidate = (
                    hint.tree_kind in {"front", "back"}
                    and hint.horizontal_drifts == 0
                )

                if is_candidate:
                    candidate_seen = reservoir_add(
                        candidate_sample,
                        item,
                        candidate_seen,
                        candidate_limit,
                        rng_candidate,
                    )
                else:
                    control_seen = reservoir_add(
                        control_sample,
                        item,
                        control_seen,
                        control_limit,
                        rng_control,
                    )

        rows = []

        for item in candidate_sample:
            core_id, reference, query, hint = item

            rows.append(
                analyse_pair(
                    worker_id,
                    "candidate",
                    core_id,
                    reference,
                    query,
                    hint,
                    align_global,
                )
            )

        for item in control_sample:
            core_id, reference, query, hint = item

            rows.append(
                analyse_pair(
                    worker_id,
                    "control",
                    core_id,
                    reference,
                    query,
                    hint,
                    align_global,
                )
            )

        result_queue.put(
            {
                "ok": True,
                "worker_id": worker_id,
                "processed_reads": processed_reads,
                "cluster_count": len(process.ref_dict),
                "candidate_seen": candidate_seen,
                "control_seen": control_seen,
                "candidate_sampled": len(candidate_sample),
                "control_sampled": len(control_sample),
                "exact_core_duplicates": exact_core_duplicates,
                "rows": rows,
            }
        )

    except Exception:
        error_event.set()

        result_queue.put(
            {
                "ok": False,
                "worker_id": worker_id,
                "traceback": traceback.format_exc(),
            }
        )


def summarize_rows(rows):
    result = {}

    for class_name in ("candidate", "control"):
        subset = [
            row for row in rows
            if row["class"] == class_name
        ]

        n = len(subset)

        edit_counter = Counter(
            row["edit_distance"]
            for row in subset
        )
        length_counter = Counter(
            row["length_difference"]
            for row in subset
        )
        diagonal_counter = Counter(
            row["max_abs_diagonal"]
            for row in subset
        )

        if n == 0:
            result[class_name] = {
                "n": 0,
            }
            continue

        same_length = sum(
            row["same_length"]
            for row in subset
        )
        no_indel = sum(
            not row["has_indel"]
            for row in subset
        )
        same_length_no_indel = sum(
            row["same_length_no_indel"]
            for row in subset
        )
        hamming_equals_edit = sum(
            row["hamming_equals_edit"]
            for row in subset
        )

        diag0 = sum(
            row["max_abs_diagonal"] == 0
            for row in subset
        )
        diag1 = sum(
            row["max_abs_diagonal"] <= 1
            for row in subset
        )
        diag2 = sum(
            row["max_abs_diagonal"] <= 2
            for row in subset
        )

        mean_edit = sum(
            row["edit_distance"]
            for row in subset
        ) / n

        result[class_name] = {
            "n": n,
            "same_length_count": same_length,
            "same_length_rate": same_length / n,
            "no_indel_count": no_indel,
            "no_indel_rate": no_indel / n,
            "same_length_no_indel_count": (
                same_length_no_indel
            ),
            "same_length_no_indel_rate": (
                same_length_no_indel / n
            ),
            "hamming_equals_edit_count": (
                hamming_equals_edit
            ),
            "hamming_equals_edit_rate": (
                hamming_equals_edit / n
            ),
            "diag_0_rate": diag0 / n,
            "diag_le_1_rate": diag1 / n,
            "diag_le_2_rate": diag2 / n,
            "mean_edit_distance": mean_edit,
            "total_substitutions": sum(
                row["substitutions"]
                for row in subset
            ),
            "total_insertions": sum(
                row["insertions"]
                for row in subset
            ),
            "total_deletions": sum(
                row["deletions"]
                for row in subset
            ),
            "edit_distance_distribution": {
                str(k): v
                for k, v in sorted(edit_counter.items())
            },
            "length_difference_distribution": {
                str(k): v
                for k, v in sorted(length_counter.items())
            },
            "max_abs_diagonal_distribution": {
                str(k): v
                for k, v in sorted(diagonal_counter.items())
            },
        }

    return result


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
                args.candidate_per_worker,
                args.control_per_worker,
                args.seed,
            ),
        )
        process.start()
        workers.append(process)

    batches = [
        []
        for _ in range(args.workers)
    ]

    dispatched = 0
    started = time.perf_counter()

    try:
        for read_id, sequence in iter_fastq(args.input):
            if (
                args.max_reads is not None
                and dispatched >= args.max_reads
            ):
                break

            owner = owner_index(
                sequence,
                prefix_length,
            )

            batches[owner].append(
                (read_id, sequence)
            )
            dispatched += 1

            if len(batches[owner]) >= args.batch_size:
                safe_put(
                    input_queues[owner],
                    batches[owner],
                    error_event,
                )
                batches[owner] = []

            if dispatched % 100000 == 0:
                elapsed = time.perf_counter() - started

                print(
                    f"dispatched={dispatched:,} "
                    f"rate={dispatched / elapsed:,.0f} reads/s"
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
                f"worker exited with {process.exitcode}"
            )

    rows = []

    for result in worker_results:
        rows.extend(result["rows"])

    summary = {
        "phase": "13.3",
        "input": args.input,
        "max_reads": args.max_reads,
        "workers": args.workers,
        "prefix_length": prefix_length,
        "candidate_definition": (
            "tree_kind in {front,back} "
            "and horizontal_drifts == 0"
        ),
        "routing_reference": "Clover routing core",
        "wfa_backend": "clover_consensus.align_global(..., backend='wfa')",
        "candidate_seen": sum(
            r["candidate_seen"]
            for r in worker_results
        ),
        "control_seen": sum(
            r["control_seen"]
            for r in worker_results
        ),
        "candidate_sampled": sum(
            r["candidate_sampled"]
            for r in worker_results
        ),
        "control_sampled": sum(
            r["control_sampled"]
            for r in worker_results
        ),
        "exact_core_duplicates_excluded": sum(
            r["exact_core_duplicates"]
            for r in worker_results
        ),
        "wall_seconds": (
            time.perf_counter() - started
        ),
        "groups": summarize_rows(rows),
    }

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = Path(str(prefix) + ".json")
    tsv_path = Path(str(prefix) + ".tsv")

    with json_path.open(
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            summary,
            f,
            indent=2,
            sort_keys=True,
        )
        f.write("\n")

    fieldnames = [
        "worker_id",
        "class",
        "core_id",
        "tree_kind",
        "horizontal_drifts",
        "query_shift",
        "reference_length",
        "query_length",
        "length_difference",
        "edit_distance",
        "substitutions",
        "insertions",
        "deletions",
        "has_indel",
        "same_length",
        "same_length_no_indel",
        "hamming_distance",
        "hamming_equals_edit",
        "max_abs_diagonal",
        "reference",
        "query",
    ]

    with tsv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    print()
    print("=== Phase 13.3 Hint vs WFA ===")
    print(
        "candidate seen:",
        f"{summary['candidate_seen']:,}",
    )
    print(
        "control seen:",
        f"{summary['control_seen']:,}",
    )
    print(
        "candidate sampled:",
        f"{summary['candidate_sampled']:,}",
    )
    print(
        "control sampled:",
        f"{summary['control_sampled']:,}",
    )

    for name in ("candidate", "control"):
        group = summary["groups"][name]

        print()
        print(name.upper())

        if group["n"] == 0:
            print("n=0")
            continue

        print("n:", group["n"])
        print(
            "same length:",
            f"{group['same_length_rate']:.6f}",
        )
        print(
            "no indel:",
            f"{group['no_indel_rate']:.6f}",
        )
        print(
            "same length + no indel:",
            f"{group['same_length_no_indel_rate']:.6f}",
        )
        print(
            "hamming == edit:",
            f"{group['hamming_equals_edit_rate']:.6f}",
        )
        print(
            "diag=0:",
            f"{group['diag_0_rate']:.6f}",
        )
        print(
            "diag<=1:",
            f"{group['diag_le_1_rate']:.6f}",
        )
        print(
            "diag<=2:",
            f"{group['diag_le_2_rate']:.6f}",
        )
        print(
            "mean edit:",
            f"{group['mean_edit_distance']:.4f}",
        )

    print()
    print("JSON:", json_path)
    print("TSV :", tsv_path)


if __name__ == "__main__":
    main()
