#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def load_truth(path, expected_length):
    sequences = []
    bad_length = 0
    bad_chars = 0

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            seq = line.strip().upper()
            if not seq:
                continue

            sequences.append(seq)

            if expected_length is not None and len(seq) != expected_length:
                bad_length += 1

            if any(base not in "ACGT" for base in seq):
                bad_chars += 1

    truth = set(sequences)

    print(f"truth_total={len(sequences)}")
    print(f"truth_unique={len(truth)}")
    print(f"truth_bad_length={bad_length}")
    print(f"truth_bad_chars={bad_chars}")

    if bad_length:
        raise RuntimeError("truth contains sequences with unexpected length")
    if bad_chars:
        raise RuntimeError("truth contains non-ACGT symbols")
    if len(sequences) != len(truth):
        raise RuntimeError(
            "truth contains duplicate sequences; sequence-set recovery "
            "cannot distinguish duplicated originals"
        )

    return truth


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--expected-length", type=int, default=None)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    truth = load_truth(args.truth, args.expected_length)

    core_hit_truths = set()
    consensus_hit_truths = set()

    core_outputs = set()
    consensus_outputs = set()

    core_exact_rows = 0
    consensus_exact_rows = 0
    cluster_count = 0

    with open(args.consensus, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        required = {"routing_core", "consensus"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                "consensus TSV missing columns: " + ", ".join(sorted(missing))
            )

        for row in reader:
            cluster_count += 1

            core = row["routing_core"].strip().upper()
            consensus = row["consensus"].strip().upper()

            core_outputs.add(core)
            consensus_outputs.add(consensus)

            if core in truth:
                core_exact_rows += 1
                core_hit_truths.add(core)

            if consensus in truth:
                consensus_exact_rows += 1
                consensus_hit_truths.add(consensus)

    both = core_hit_truths & consensus_hit_truths
    consensus_only = consensus_hit_truths - core_hit_truths
    core_only = core_hit_truths - consensus_hit_truths
    neither = truth - (core_hit_truths | consensus_hit_truths)

    truth_n = len(truth)

    core_recovery = len(core_hit_truths) / truth_n
    consensus_recovery = len(consensus_hit_truths) / truth_n

    gain = len(consensus_hit_truths) - len(core_hit_truths)
    gain_pp = (consensus_recovery - core_recovery) * 100.0

    metrics = [
        ("clusters", cluster_count),
        ("routing_core_exact_rows", core_exact_rows),
        ("routing_core_exact_truths", len(core_hit_truths)),
        ("routing_core_exact_recovery_pct", f"{core_recovery * 100:.6f}"),
        ("consensus_exact_rows", consensus_exact_rows),
        ("consensus_exact_truths", len(consensus_hit_truths)),
        ("consensus_exact_recovery_pct", f"{consensus_recovery * 100:.6f}"),
        ("consensus_gain_truths", gain),
        ("consensus_gain_percentage_points", f"{gain_pp:.6f}"),
        ("both_truths", len(both)),
        ("consensus_only_truths", len(consensus_only)),
        ("core_only_truths", len(core_only)),
        ("neither_truths", len(neither)),
        ("unique_routing_core_outputs", len(core_outputs)),
        ("unique_consensus_outputs", len(consensus_outputs)),
        ("consensus_duplicate_output_rows", cluster_count - len(consensus_outputs)),
        ("consensus_exact_duplicate_rows", consensus_exact_rows - len(consensus_hit_truths)),
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as f:
        f.write(f"truth_total={truth_n}\n")
        for key, value in metrics:
            f.write(f"{key}={value}\n")

    print()
    print("===== EXACT RECOVERY =====")
    print(f"truth_total={truth_n}")
    for key, value in metrics:
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
