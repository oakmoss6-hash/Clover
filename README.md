# Clover with Multi-Read Global Reconstruction

This repository is a fork of [Guanjinqu/Clover](https://github.com/Guanjinqu/Clover), extended with an optional **cluster-level multi-read global reconstruction** stage for DNA-storage sequencing data.

The original Clover algorithm is responsible for fast clustering. The reconstruction extension does **not** replace Clover routing. Instead, it observes the reads that Clover has already assigned to each cluster, compresses exact duplicates, globally aligns the unique sequences to one shared backbone, and produces one weighted consensus sequence per cluster.

## What this fork adds

For a Clover cluster containing `M` raw reads and `U` unique sequences:

1. Aggregate exact duplicates as `sequence -> multiplicity`.
2. If `U = 1`, return the sequence directly with no alignment.
3. Select one truth-blind reconstruction backbone.
4. Globally align every other unique sequence to the backbone exactly once.
5. Project all weighted evidence onto one shared coordinate profile.
6. Call a deterministic weighted consensus.

Therefore a cluster performs exactly:

```text
U - 1
```

pairwise global alignments. Across `C` reconstructed clusters:

```text
A = sum(U_c - 1) = U - C
```

This avoids all-pairs sequence comparison inside a cluster.

## Architecture

```text
FASTQ / FASTA / Clover TXT
            |
            v
     streaming input
            |
            v
      Clover routing
            |
            v
  sequence -> multiplicity
        (M -> U)
            |
       +----+----+
       |         |
     U = 1     U > 1
       |         |
  direct output  v
          reconstruction backbone
                  |
                  v
          U - 1 global alignments
                  |
                  v
          shared weighted profile
                  |
                  v
               consensus
```

Reconstruction mode uses streaming input by default so the complete raw read set is not preloaded into the parent process.

## Installation

Clone this fork and install the dependencies from the repository root:

```bash
git clone <YOUR-FORK-URL>
cd Clover
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The production default uses `pywfa==0.5.1` for exact end-to-end Wavefront Alignment.

Optional backends:

```bash
python -m pip install edlib
```

`nw` is also available as a pure-Python correctness/reference backend and requires no additional package.

## Quick start: Clover + reconstruction

```bash
python -m clover \
  --input reads.fastq \
  -P 2 \
  --reconstruct \
  --reconstruct-backend wfa \
  --consensus-output consensus.tsv
```

`-P` follows the original Clover convention:

```text
-P 0 -> 1 worker
-P 1 -> 4 workers
-P 2 -> 16 workers
```

### Reconstruction options

| Option | Default | Description |
|---|---|---|
| `--reconstruct` | off | Enable cluster-level multi-read reconstruction. |
| `--reconstruct-backend` | `wfa` | Pairwise backend: `wfa`, `edlib`, or `nw`. |
| `--reconstruct-backbone` | `core` | Backbone policy: `core`, `support_length`, or `max_span`. |
| `--consensus-output` | `clover_consensus.tsv` | Output TSV path. |

The validated benchmark default is:

```text
backend  = wfa
backbone = core
```

`--align` is the original Clover global-matching feature and is separate from the new reconstruction stage. `--align` and `--reconstruct` are mutually exclusive.

## Input formats

The reconstruction CLI accepts:

### FASTQ

Standard four-line FASTQ records.

### FASTA

Single- or multi-line FASTA records.

### Clover text

Two whitespace-separated columns:

```text
read_id sequence
```

Reads containing `N` or shorter than Clover's configured minimum length are filtered by the normal Clover processing path.

## Output

`--consensus-output` writes one row per Clover output cluster:

```text
worker
cluster_id
routing_core
backbone
consensus
raw_read_count
unique_sequence_count
pairwise_alignment_count
```

The implementation checks the global invariant:

```text
pairwise_alignment_count = total_unique_sequences - cluster_count
```

and raises an error if it is violated.

## Full ERR1816980 result

The following result was obtained with the full assembled ERR1816980 dataset used during development of this fork.

| Metric | Value |
|---|---:|
| Input reads | 15,787,115 |
| Original source strands | 72,000 |
| Workers | 16 (`-P 2`) |
| Clover output clusters | 153,813 |
| Assigned reconstruction reads `M` | 13,402,459 |
| Unique sequence evidence `U` | 3,963,397 |
| Pairwise global alignments `A` | 3,809,584 |
| Verified invariant | `A = U - C` |
| End-to-end wall time | 140.26 s |

### Exact distinct strand recovery

A source strand counts as recovered when at least one output sequence is exactly identical to that source reference. Duplicate output clusters matching the same source strand are counted only once.

| Output used for recovery | Exact strands recovered | Recovery rate |
|---|---:|---:|
| Clover routing cores | 48,703 / 72,000 | 67.6431% |
| Clover + multi-read reconstruction | **71,823 / 72,000** | **99.7542%** |
| Net improvement | **+23,120 strands** | **+32.1111 percentage points** |

More specifically, reconstruction recovered 23,175 source strands that were absent from the exact routing-core set, while 55 strands present in the core set were absent from the final consensus set.

> **Metric note:** `99.7542%` is an **exact distinct source-strand recovery rate**. It is not the same metric as Clover's clustering accuracy reported in the original paper.

The wall-time number is a measurement from the development machine (WSL2, 24 logical CPUs) and should not be interpreted as a hardware-independent guarantee.

## Reconstruction code layout

```text
clover/reconstruction/
├── __init__.py     public reconstruction API
├── alignment.py    WFA / Edlib / Needleman-Wunsch backends
├── backbone.py     truth-blind backbone selection
├── state.py        compact streaming cluster evidence
├── reconstruct.py  shared profile and consensus algorithm
└── worker.py       Clover multiprocessing integration and TSV output
```

Experimental simulation, diagnostics, baseline comparisons, and previous local-repair experiments remain research code and are not part of the production reconstruction path.

## Original Clover usage

The original clustering options remain available. Common options include:

| Option | Description |
|---|---|
| `-I`, `--input` | Input file. |
| `-L` | Expected read length. |
| `-P` | Clover process exponent (`4^P` workers for `P > 0`). |
| `-D` | End-tree depth. |
| `-V` | Vertical drift setting. |
| `-H` | Horizontal drift setting. |
| `--no-tag` | Untagged input mode. |
| `--no-fast` | Original lower-memory Clover input mode. |
| `--low` | Original minimum-memory Clover mode. |
| `--align` | Original Clover global-matching feature. |

For the original algorithm and parameter definitions, see the upstream Clover repository and paper.

## Testing

From the repository root:

```bash
python -m unittest discover -s tests -v
```

The development branch additionally contains research regression tests used to validate backend equivalence, the `U-1` alignment bound, Clover worker integration, coverage behavior, and real-data reconstruction.

## Relationship to upstream Clover

This fork preserves Clover's clustering/routing logic and adds reconstruction as an optional downstream stage. Reconstruction receives only already-decided cluster memberships and does not use reference truth, Bowtie labels, or known source sequences during normal execution.

## License

Clover is distributed under the GNU General Public License. See `LICENSE` for the repository license terms.

## Citation

For the original Clover clustering algorithm, please cite:

> Qu G, Yan Z, Wu H. **Clover: tree structure-based efficient DNA clustering for DNA-based data storage.** Briefings in Bioinformatics. 2022;23(5):bbac336.

If you use the reconstruction extension in this fork, please also cite the corresponding reconstruction work once its citation is available.
