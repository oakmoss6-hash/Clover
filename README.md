# Clover

Clover is a tree-structure-based DNA clustering tool for DNA data storage.
This repository contains the original Clover clustering algorithm plus Reed,
an optional cluster-level multi-read global reconstruction module.

The reconstruction module runs after Clover has determined cluster
membership. It observes the reads already assigned to each cluster, compresses
exact duplicate sequences into multiplicities, aligns unique sequences to one
cluster backbone, and emits one weighted consensus sequence per cluster.

Bowtie, BLAST, source labels, and known reference strands are not part of the
normal reconstruction algorithm.

## Clover Clustering

Clover routes reads through prefix, suffix, and middle trie searches with
configurable drift thresholds. The clustering path is responsible for deciding
which reads belong to each cluster and for producing Clover's original cluster
statistics.

Common clustering options:

| Option | Description |
|---|---|
| `-I`, `--input` | Input FASTQ, FASTA, or Clover text file. |
| `-L` | Expected read length. |
| `-D` | End-tree depth. |
| `-V` | Vertical drift setting. |
| `-H` | Horizontal drift threshold. |
| `-T` | Expected tag count. |
| `-P` | Process exponent: `0` uses one worker, `N > 0` uses `4^N` workers. |
| `-O` | Write original Clover cluster index output. |
| `--no-tag` | Use untagged input mode. |
| `--no-fast` | Use Clover's lower-memory input mode. |
| `--low` | Use Clover's minimum-memory mode. |

## Multi-Read Reconstruction

Reconstruction is an optional downstream stage for cluster-level consensus
calling. It does not replace Clover clustering and does not participate in
routing decisions.

For each cluster:

1. Store observed read evidence as `sequence -> multiplicity`.
2. Return singleton clusters directly without alignment.
3. Select one truth-blind backbone.
4. Globally align each non-backbone unique sequence once.
5. Project weighted evidence onto a shared profile.
6. Call a deterministic weighted consensus.

For a cluster with `U` unique sequences, reconstruction performs exactly
`U - 1` pairwise global alignments.

## Installation

Install from the repository root:

```bash
python -m pip install -r requirements.txt
python -m pip install .
```

The package installs `clover`; Reed is implemented in `clover.reed`.

## Usage

Run Clover clustering:

```bash
python -m clover --input reads.fastq -P 0
```

Run Clover clustering followed by reconstruction:

```bash
python -m clover \
  --input reads.fastq \
  -P 0 \
  --reconstruct \
  --reconstruct-backend wfa \
  --consensus-output clover_consensus.tsv
```

FASTQ, FASTA, and Clover text inputs are accepted. Clover text input uses two
whitespace-separated columns:

```text
read_id sequence
```

Reads containing `N` or shorter than Clover's configured minimum length are
filtered by Clover's normal processing path.

## Reconstruction Options

| Option | Default | Description |
|---|---|---|
| `--reconstruct` | off | Enable cluster-level multi-read reconstruction. |
| `--reconstruct-backend` | `wfa` | Pairwise global-alignment backend: `wfa`, `edlib`, or `nw`. |
| `--reconstruct-backbone` | `core` | Backbone policy: `core`, `support_length`, or `max_span`. |
| `--consensus-output` | `clover_consensus.tsv` | Consensus TSV output path. |

## Alignment Backend

The default backend is WFA through `pywfa==0.5.1`.

`edlib` is supported as an optional backend when the Python package is
installed. The `nw` backend is a pure-Python Needleman-Wunsch reference
implementation and requires no additional package.

## Output

Reconstruction writes a TSV file with one row per reconstructed Clover
cluster:

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

The implementation verifies the global invariant:

```text
pairwise_alignment_count = total_unique_sequences - cluster_count
```

## Citation

For the original Clover clustering algorithm, please cite:

Qu G, Yan Z, Wu H. Clover: tree structure-based efficient DNA clustering for
DNA-based data storage. Briefings in Bioinformatics. 2022;23(5):bbac336.

## License

Clover is distributed under the GNU General Public License. See `LICENSE` for
the repository license terms.
