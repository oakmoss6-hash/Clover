# Cluster Adapter Design

## Input Format

Clover txt input is parsed in `clover/main.py` with `line.split()`.
The first two whitespace-delimited columns are used:

- column 1: read identifier in non-virtual/low-memory modes; tag in virtual mode
- column 2: DNA sequence

For the adapter, column 1 is treated as the stable `read_id` / join key and is also stored as `ReadRecord.tag`. Column 2 is the DNA sequence. Extra columns are ignored by Clover's current parser and are not used by this adapter.

The input file order matters to Clover internally because `dna_num` is assigned from the processing counter, and core creation depends on stream order. The adapter does not change that order; it only joins existing records by explicit identifiers.

Duplicate read IDs are rejected. Duplicate sequences with different read IDs are allowed.

## Membership Format

Clover writes `index_list` as the Python string representation of a list of 2-tuples:

```text
[(read_id, cluster_or_core_label), ...]
```

In non-virtual mode, `read_id` comes from column 1 of the input. The cluster/core label comes from the stored first value in `ref_dict[core_id]`. The output does not include unmatched/discarded reads and does not include full read sequences. It also does not include core sequences.

The adapter parses this format with `ast.literal_eval`, requires a list of 2-tuples, and rejects malformed rows, duplicate memberships, and conflicting memberships.

## Core Sequence Availability

Existing Clover output is insufficient to reliably recover `core_sequence` without modifying Clover or providing an explicit core export. The adapter therefore refuses to guess core sequences from first/longest/most-common reads and raises `CORE_SEQUENCE_UNAVAILABLE` when no reliable core sequence source is supplied.

## Memory Modes

Mode A: in-memory. The first implementation reads input records and membership into dictionaries, groups members by cluster id, and builds `ClusterRecord` objects. This is suitable for small development tests and Phase 4 validation.

Mode B: stream / grouped export. Future large-scale experiments should consume membership sorted/grouped by cluster and stream matching input records or pre-built shards. The public API keeps parsing and construction separate so a streaming implementation can replace the in-memory grouping later.
