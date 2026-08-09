#!/usr/bin/env python3
from __future__ import annotations

import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_FILES = [ROOT / 'example' / 'example_tag_data.txt', ROOT / 'example' / 'example_index_data.txt']


def parse_rows(path: Path):
    rows = []
    for raw in path.read_text(encoding='utf-8', errors='replace').splitlines():
        parts = raw.strip().split()
        if len(parts) >= 2:
            rows.append((parts[0], parts[1].upper()))
    return rows


def median(values):
    return statistics.median(values) if values else 0


def positional_differences(seqs):
    if not seqs:
        return []
    ref = seqs[0]
    out = []
    for seq in seqs[1:]:
        if len(seq) == len(ref):
            diffs = [(i, ref[i], seq[i]) for i in range(len(ref)) if ref[i] != seq[i]]
            if diffs:
                out.append(diffs)
    return out


def summarize_file(path: Path):
    rows = parse_rows(path)
    seqs = [seq for _, seq in rows]
    lengths = [len(seq) for seq in seqs]
    counts = Counter(seqs)
    duplicates = sum(v - 1 for v in counts.values() if v > 1)
    return {
        'file': str(path.relative_to(ROOT)),
        'file_size_bytes': path.stat().st_size,
        'read_count': len(seqs),
        'unique_sequence_count': len(counts),
        'exact_duplicate_count': duplicates,
        'duplicate_ratio': duplicates / len(seqs) if seqs else 0,
        'min_length': min(lengths) if lengths else 0,
        'median_length': median(lengths),
        'max_length': max(lengths) if lengths else 0,
        'all_same_length': len(set(lengths)) <= 1,
        'alphabet': ''.join(sorted(set(''.join(seqs)))),
    }


def tag_report(path: Path):
    rows = parse_rows(path)
    grouped = defaultdict(list)
    for tag, seq in rows:
        grouped[tag].append(seq)
    lines = []
    observed_duplicate = False
    observed_same_len_var = False
    observed_diff_len = False
    potential_indel = False
    for tag in sorted(grouped):
        seqs = grouped[tag]
        lengths = [len(s) for s in seqs]
        counts = Counter(seqs)
        observed_duplicate = observed_duplicate or any(v > 1 for v in counts.values())
        observed_diff_len = observed_diff_len or len(set(lengths)) > 1
        if len(set(seqs)) > 1 and len(set(lengths)) == 1:
            observed_same_len_var = True
        if len(set(lengths)) > 1:
            potential_indel = True
        lines.append(f'tag={tag} read_count={len(seqs)} unique_read_count={len(counts)} min_length={min(lengths)} max_length={max(lengths)} all_identical={len(counts)==1} has_variant={len(counts)>1} has_different_length={len(set(lengths))>1} has_exact_duplicate={any(v > 1 for v in counts.values())}')
        diffs = positional_differences(seqs)
        if diffs:
            lines.append(f'  positional_differences_vs_first_same_length={diffs}')
    feature_lines = [
        'exact duplicate: ' + ('OBSERVED' if observed_duplicate else 'NOT OBSERVED'),
        'substitution-like same-length variation: ' + ('OBSERVED' if observed_same_len_var else 'NOT OBSERVED'),
        'different-length variation: ' + ('OBSERVED' if observed_diff_len else 'NOT OBSERVED'),
        'truncated reads: UNKNOWN',
        'insertion/deletion: ' + ('POTENTIAL' if potential_indel else 'NOT OBSERVED'),
    ]
    return len(grouped), lines, feature_lines


def main():
    out_dir = ROOT / 'research' / 'results'
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [summarize_file(path) for path in EXAMPLE_FILES]
    csv_path = out_dir / 'example_dataset_audit.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    lines = ['# Clover Example Dataset Audit', '']
    for row in rows:
        lines.append(f"## {row['file']}")
        for key, value in row.items():
            if key != 'file':
                lines.append(f'- {key}: {value}')
        lines.append('')
    tag_count, tag_lines, feature_lines = tag_report(ROOT / 'example' / 'example_tag_data.txt')
    lines.append('## example/example_tag_data.txt tags')
    lines.append(f'- number_of_tags: {tag_count}')
    lines.extend(f'- {line}' for line in tag_lines)
    lines.append('')
    lines.append('## Feature observations')
    lines.extend(f'- {line}' for line in feature_lines)
    (out_dir / 'example_dataset_audit.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
