#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing as mp
import queue
import sys
import time
from itertools import product
from pathlib import Path


ALPHABET = "ACGT"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Phase 14.2B real Clover -> WFA reconstruction smoke"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--max-reads", type=int, default=20000)
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        choices=[1, 4, 16, 64],
    )
    parser.add_argument(
        "--backend",
        default="wfa",
        choices=["wfa", "edlib", "nw"],
    )
    parser.add_argument(
        "--backbone-policy",
        default="core",
        choices=["core", "support_length", "max_span"],
    )
    parser.add_argument("--output-prefix", required=True)
    return parser.parse_args()


def iter_fastq(path):
    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
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
            "workers must be an exact power of four"
        )

    return k


def worker_names(workers):
    k = prefix_length_for_workers(workers)

    if k == 0:
        return [""]

    return [
        "".join(chars)
        for chars in product(ALPHABET, repeat=k)
    ]


def partition_records(records, workers):
    names = worker_names(workers)
    k = prefix_length_for_workers(workers)

    name_to_index = {
        name: index
        for index, name in enumerate(names)
    }

    shards = [
        []
        for _ in names
    ]

    invalid_prefix_reads = 0

    for read_id, sequence in records:
        if k == 0:
            owner = 0
        else:
            prefix = sequence[:k]

            if (
                len(prefix) == k
                and all(base in ALPHABET for base in prefix)
            ):
                owner = name_to_index[prefix]
            else:
                owner = 0
                invalid_prefix_reads += 1

        shards[owner].append(
            f"{read_id} {sequence}"
        )

    return names, shards, invalid_prefix_reads


def import_components():
    # Clover's configuration parser inspects sys.argv on import.
    saved_argv = sys.argv[:]
    sys.argv = [saved_argv[0]]

    try:
        from clover.main import MyProcess
        from clover_consensus.clover_worker_reconstruction import (
            CloverWorkerReconstructor,
        )
    finally:
        sys.argv = saved_argv

    return MyProcess, CloverWorkerReconstructor


def make_process(
    MyProcess,
    name,
    data,
    q_output,
):
    saved_argv = sys.argv[:]
    sys.argv = [saved_argv[0]]

    try:
        process = MyProcess(
            name,
            data,
            q_output,
        )
    finally:
        sys.argv = saved_argv

    return process


def detect_worker_name(output, names):
    matches = [
        name
        for name in names
        if name + "sum_read_num" in output
    ]

    if len(matches) != 1:
        raise RuntimeError(
            "could not uniquely identify worker output: "
            f"{matches}"
        )

    return matches[0]


def run_pass(
    *,
    MyProcess,
    CloverWorkerReconstructor,
    names,
    shards,
    reconstruct,
    backend,
    backbone_policy,
):
    q_output = mp.Queue(maxsize=max(2, len(names) * 2))
    processes = []

    for name, data in zip(names, shards):
        process = make_process(
            MyProcess,
            name,
            data,
            q_output,
        )

        if reconstruct:
            reconstructor = CloverWorkerReconstructor(
                worker_name=name,
                backend=backend,
                backbone_policy=backbone_policy,
                local_repeat_repair=False,
            )

            reconstructor.attach_to_process(process)

            # Retain an owner-side reference until process.start().
            # The child receives it through the process hooks.
            process._phase14_reconstructor = reconstructor

        processes.append(process)

    started = time.perf_counter()

    for process in processes:
        process.start()

    outputs = {}

    try:
        while len(outputs) < len(processes):
            try:
                output = q_output.get(timeout=120)
            except queue.Empty:
                states = [
                    (
                        process.name,
                        process.is_alive(),
                        process.exitcode,
                    )
                    for process in processes
                ]
                raise RuntimeError(
                    "timed out waiting for Clover worker output: "
                    f"{states}"
                )

            name = detect_worker_name(
                output,
                names,
            )

            if name in outputs:
                raise RuntimeError(
                    f"duplicate output for worker {name!r}"
                )

            outputs[name] = output

    finally:
        for process in processes:
            process.join(timeout=120)

        for process in processes:
            if process.is_alive():
                process.terminate()
                process.join()
                raise RuntimeError(
                    f"Clover worker {process.name} did not terminate"
                )

            if process.exitcode != 0:
                raise RuntimeError(
                    f"Clover worker {process.name} exited "
                    f"with code {process.exitcode}"
                )

        q_output.close()

    wall_seconds = time.perf_counter() - started

    return outputs, wall_seconds


def strip_reconstruction_fields(output):
    return {
        key: value
        for key, value in output.items()
        if "reconstruction_" not in key
    }


def compare_clover_outputs(
    baseline_outputs,
    reconstruction_outputs,
    names,
):
    differences = []

    for name in names:
        baseline = strip_reconstruction_fields(
            baseline_outputs[name]
        )
        reconstructed = strip_reconstruction_fields(
            reconstruction_outputs[name]
        )

        if baseline != reconstructed:
            differing_keys = sorted(
                key
                for key in set(baseline) | set(reconstructed)
                if baseline.get(key) != reconstructed.get(key)
            )

            differences.append(
                {
                    "worker": name,
                    "keys": differing_keys,
                }
            )

    return differences


def collect_reconstruction(
    outputs,
    names,
):
    rows = []
    worker_details = []

    total_clusters = 0
    total_raw = 0
    total_unique = 0
    total_pairwise = 0
    total_singletons = 0
    total_multi_unique = 0

    for name in names:
        output = outputs[name]

        cluster_count = output[
            name + "reconstruction_cluster_count"
        ]
        raw_count = output[
            name + "reconstruction_raw_read_count"
        ]
        unique_count = output[
            name + "reconstruction_unique_sequence_count"
        ]
        pairwise_count = output[
            name + "reconstruction_pairwise_alignment_count"
        ]
        singleton_count = output[
            name + "reconstruction_singleton_cluster_count"
        ]
        multi_unique_count = output[
            name + "reconstruction_multi_unique_cluster_count"
        ]
        worker_rows = output[
            name + "reconstruction_results"
        ]

        total_clusters += cluster_count
        total_raw += raw_count
        total_unique += unique_count
        total_pairwise += pairwise_count
        total_singletons += singleton_count
        total_multi_unique += multi_unique_count

        worker_details.append(
            {
                "worker": name,
                "input_reads": None,
                "cluster_count": cluster_count,
                "raw_clustered_reads": raw_count,
                "unique_sequences": unique_count,
                "pairwise_alignments": pairwise_count,
                "singleton_clusters": singleton_count,
                "multi_unique_clusters": multi_unique_count,
            }
        )

        for (
            cluster_id,
            consensus,
            backbone,
            raw_cluster_count,
            unique_cluster_count,
            cluster_pairwise_count,
        ) in worker_rows:
            rows.append(
                {
                    "worker": name,
                    "cluster_id": cluster_id,
                    "consensus": consensus,
                    "backbone": backbone,
                    "raw_read_count": raw_cluster_count,
                    "unique_sequence_count": unique_cluster_count,
                    "pairwise_alignment_count": (
                        cluster_pairwise_count
                    ),
                }
            )

    expected_pairwise = total_unique - total_clusters

    if total_pairwise != expected_pairwise:
        raise RuntimeError(
            "global U-C invariant failed: "
            f"{total_pairwise} != "
            f"{total_unique} - {total_clusters}"
        )

    return {
        "rows": rows,
        "worker_details": worker_details,
        "cluster_count": total_clusters,
        "raw_clustered_reads": total_raw,
        "unique_sequences": total_unique,
        "pairwise_alignments": total_pairwise,
        "singleton_clusters": total_singletons,
        "multi_unique_clusters": total_multi_unique,
    }


def main():
    args = parse_args()

    MyProcess, CloverWorkerReconstructor = (
        import_components()
    )

    records = []

    for read_id, sequence in iter_fastq(args.input):
        if (
            args.max_reads is not None
            and len(records) >= args.max_reads
        ):
            break

        records.append(
            (read_id, sequence)
        )

    # Match Clover's current input eligibility guard:
    # len(sequence) >= read_len_min and no ambiguous N.
    saved_argv = sys.argv[:]
    sys.argv = [saved_argv[0]]
    try:
        probe = MyProcess("phase14-config-probe", [], None)
        read_len_min = probe.config_dict["read_len_min"]
    finally:
        sys.argv = saved_argv

    eligible_reads = sum(
        1
        for _, sequence in records
        if len(sequence) >= read_len_min
        and "N" not in sequence
    )

    names, shards, invalid_prefix_reads = (
        partition_records(
            records,
            args.workers,
        )
    )

    print(
        f"loaded {len(records):,} reads "
        f"into {args.workers} Clover workers"
    )

    print(
        "shards:",
        {
            name or "<single>": len(shard)
            for name, shard in zip(names, shards)
        },
    )

    # ------------------------------------------------------------
    # Pass A: historical Clover only.
    # ------------------------------------------------------------
    baseline_outputs, baseline_seconds = run_pass(
        MyProcess=MyProcess,
        CloverWorkerReconstructor=CloverWorkerReconstructor,
        names=names,
        shards=shards,
        reconstruct=False,
        backend=args.backend,
        backbone_policy=args.backbone_policy,
    )

    print(
        f"baseline Clover: {baseline_seconds:.4f}s"
    )

    # ------------------------------------------------------------
    # Pass B: Clover + worker-local reconstruction.
    # ------------------------------------------------------------
    reconstruction_outputs, integrated_seconds = run_pass(
        MyProcess=MyProcess,
        CloverWorkerReconstructor=CloverWorkerReconstructor,
        names=names,
        shards=shards,
        reconstruct=True,
        backend=args.backend,
        backbone_policy=args.backbone_policy,
    )

    print(
        "Clover + reconstruction:",
        f"{integrated_seconds:.4f}s",
    )

    differences = compare_clover_outputs(
        baseline_outputs,
        reconstruction_outputs,
        names,
    )

    if differences:
        raise RuntimeError(
            "reconstruction changed historical Clover output: "
            + json.dumps(differences)
        )

    reconstruction = collect_reconstruction(
        reconstruction_outputs,
        names,
    )

    for detail, shard in zip(
        reconstruction["worker_details"],
        shards,
    ):
        detail["input_reads"] = len(shard)

    M = reconstruction["raw_clustered_reads"]
    U = reconstruction["unique_sequences"]
    C = reconstruction["cluster_count"]
    A = reconstruction["pairwise_alignments"]
    singleton_clusters = reconstruction["singleton_clusters"]
    multi_unique_clusters = reconstruction["multi_unique_clusters"]

    summary = {
        "phase": "14.2B",
        "input": str(Path(args.input)),
        "max_reads": args.max_reads,
        "processed_reads": len(records),
        "eligible_reads": eligible_reads,
        "below_eligibility_filter": (
            len(records) - eligible_reads
        ),
        "assigned_reads": M,
        "eligible_unassigned_reads": (
            eligible_reads - M
        ),
        "assignment_rate_of_processed": (
            M / len(records) if records else 0.0
        ),
        "assignment_rate_of_eligible": (
            M / eligible_reads if eligible_reads else 0.0
        ),
        "workers": args.workers,
        "backend": args.backend,
        "backbone_policy": args.backbone_policy,
        "local_repeat_repair": False,
        "target_length": None,
        "invalid_prefix_reads": invalid_prefix_reads,
        "clover_outputs_identical": True,
        "baseline_clover_wall_seconds": baseline_seconds,
        "integrated_wall_seconds": integrated_seconds,
        "reconstruction_overhead_seconds": (
            integrated_seconds - baseline_seconds
        ),
        "integrated_reads_per_second": (
            len(records) / integrated_seconds
            if integrated_seconds
            else 0.0
        ),
        "cluster_count_C": C,
        "clustered_raw_reads_M": M,
        "unique_sequences_U": U,
        "pairwise_alignments_A": A,
        "singleton_cluster_count": singleton_clusters,
        "multi_unique_cluster_count": multi_unique_clusters,
        "singleton_cluster_rate": (
            singleton_clusters / C if C else 0.0
        ),
        "expected_pairwise_U_minus_C": U - C,
        "pairwise_invariant_ok": A == U - C,
        "duplicate_compression_M_over_U": (
            M / U if U else 0.0
        ),
        "worker_details": reconstruction[
            "worker_details"
        ],
    }

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    json_path = Path(str(prefix) + ".json")
    tsv_path = Path(str(prefix) + ".consensus.tsv")

    with json_path.open(
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

    fieldnames = [
        "worker",
        "cluster_id",
        "consensus",
        "backbone",
        "raw_read_count",
        "unique_sequence_count",
        "pairwise_alignment_count",
    ]

    with tsv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(
            reconstruction["rows"]
        )

    print()
    print("=== Phase 14.2B real Clover E2E ===")
    print(
        "processed reads:",
        f"{len(records):,}",
    )
    print(
        "Clover outputs identical:",
        summary["clover_outputs_identical"],
    )
    print(
        "baseline Clover:",
        f"{baseline_seconds:.4f}s",
    )
    print(
        "integrated:",
        f"{integrated_seconds:.4f}s",
    )
    print(
        "reconstruction overhead:",
        f"{summary['reconstruction_overhead_seconds']:.4f}s",
    )
    print(
        "eligible reads:",
        f"{eligible_reads:,}",
    )
    print(
        "below eligibility filter:",
        f"{len(records) - eligible_reads:,}",
    )
    print(
        "eligible but unassigned:",
        f"{eligible_reads - M:,}",
    )
    print(
        "assignment / eligible:",
        f"{M / eligible_reads:.2%}" if eligible_reads else "n/a",
    )
    print(
        "clusters C:",
        f"{C:,}",
    )
    print(
        "clustered raw M:",
        f"{M:,}",
    )
    print(
        "unique U:",
        f"{U:,}",
    )
    print(
        "M/U:",
        f"{summary['duplicate_compression_M_over_U']:.4f}",
    )
    print(
        "singleton clusters:",
        f"{singleton_clusters:,}",
        f"({singleton_clusters / C:.2%})" if C else "",
    )
    print(
        "multi-unique clusters:",
        f"{multi_unique_clusters:,}",
    )
    print(
        "pairwise A:",
        f"{A:,}",
    )
    print(
        "U-C:",
        f"{U - C:,}",
    )
    print(
        "U-C invariant:",
        summary["pairwise_invariant_ok"],
    )
    print(
        "throughput:",
        f"{summary['integrated_reads_per_second']:,.0f} reads/s",
    )
    print(
        "summary:",
        json_path,
    )
    print(
        "consensus:",
        tsv_path,
    )


if __name__ == "__main__":
    main()
