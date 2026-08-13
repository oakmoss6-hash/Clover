# Clover-Reed

**Clover-Reed** is a research implementation that extends the original **Clover** DNA clustering framework with **Reed**, a cluster-level multi-sequence global alignment and reconstruction method designed for highly homologous DNA reads.

The key assumption is different from general-purpose read mapping tools such as BLAST or Bowtie: **Clover has already grouped reads that are expected to originate from the same DNA strand**. Therefore, Reed does not search for which reference a read belongs to. Instead, it focuses on efficiently identifying substitutions, insertions, and deletions within each cluster and reconstructing a consensus sequence.

> Current status: research prototype for DNA storage sequence reconstruction.

---

## 1. Motivation

Clover is a tree-structure-based DNA sequence clustering algorithm. Its clustering stage is designed to be efficient and does not require sequence-to-sequence global alignment. The original Clover framework also leaves room for an external global matching/alignment method.

Reed is added **after Clover clustering**:

```text
DNA reads
   |
   v
Clover clustering
   |
   v
High-homology clusters
   |
   v
Reed reconstruction
   |
   +--> duplicate compression
   +--> backbone selection
   +--> global alignment
   +--> shared profile
   +--> weighted consensus
   |
   v
Reconstructed DNA sequences
```

The main research objective is to improve reconstruction accuracy while keeping the computational cost with respect to the number of sequences as low as possible.

---

## 2. What Reed Does

For each Clover cluster, Reed performs the following steps.

### 2.1 Duplicate compression

Identical reads inside a cluster are compressed into one unique sequence with a multiplicity count.

Instead of aligning the same sequence many times:

```text
ACGT... x 40
ACGG... x 12
ACGA... x  5
```

Reed stores:

```text
ACGT... -> weight 40
ACGG... -> weight 12
ACGA... -> weight 5
```

The multiplicity is later used as evidence in consensus generation.

### 2.2 Shared backbone selection

One representative sequence is selected as the backbone of the cluster.

Current backbone policies include:

- `core`
- `support_length`
- `max_span`

The default experimental configuration uses the Clover routing core as the backbone.

### 2.3 Star-shaped global alignment

Let a cluster contain \(U_i\) unique sequences.

Reed does **not** perform all-pairs alignment.

Instead, every non-backbone unique sequence is globally aligned to the shared backbone exactly once:

```text
             read 1
               |
             read 2
               |
read 3 ---- backbone ---- read 4
               |
             read 5
               |
              ...
```

Therefore, the number of pairwise global alignments for one cluster is:

\[
A_i = U_i - 1
\]

For all \(C\) clusters, if the total number of unique sequences is \(U\):

\[
A = \sum_{i=1}^{C}(U_i-1)=U-C
\]

This avoids the \(O(U_i^2)\) number of pairwise comparisons required by an all-pairs strategy.

### 2.4 Global alignment backends

Reed currently supports:

- **WFA** — default backend for large-scale reconstruction
- **Edlib** — optional backend
- **Needleman-Wunsch (NW)** — reference/fallback implementation

WFA is used as the primary backend because the reads within a Clover cluster are expected to be highly homologous.

### 2.5 Shared profile projection

Each pairwise alignment is projected onto the same backbone coordinate system.

Reed accumulates evidence for:

- matches/substitutions,
- deletions,
- insertions between backbone positions.

This is important for DNA storage reads because insertion/deletion errors shift sequence coordinates and cannot be handled reliably by simple position-wise majority voting.

### 2.6 Weighted consensus

The shared profile combines evidence from all unique sequences.

Duplicate multiplicity is retained as weight, so a sequence observed many times contributes more evidence than a sequence observed once.

The final profile is converted into one reconstructed consensus sequence per Clover cluster.

---

## 3. Why Reed Is Different from General Read Mapping

BLAST, Bowtie and similar tools are mainly designed to answer a mapping/search problem:

> Which reference or genomic position best matches this read?

Reed addresses a different problem:

> These reads have already been assigned to the same high-homology cluster. What are their differences, and what sequence best represents their common source?

Therefore, Reed is designed as a **cluster-level reconstruction method**, not a database search or read-mapping algorithm.

---

## 4. Computational Design

Assume:

- \(M\): number of raw reads entering reconstruction,
- \(U\): number of unique sequences after duplicate compression,
- \(C\): number of Clover clusters,
- \(L\): typical sequence length.

Reed performs:

```text
duplicate compression        -> one pass over cluster reads
backbone selection           -> no all-pairs distance matrix
global alignment             -> U - C pairwise alignments
profile accumulation         -> one projection per alignment
consensus generation         -> one profile scan per cluster
```

The **number of pairwise global alignments grows linearly with the number of unique sequences**:

\[
A = U-C
\]

This is the main complexity property of the current implementation.

The total runtime is not strictly \(O(U)\) in all possible cases because the cost of each global alignment also depends on sequence length, edit distance, and the selected alignment backend. Reed's design specifically targets DNA-storage clusters where sequence length is limited and within-cluster edit distance is expected to be small.

---

## 5. Current Experimental Results

### 5.1 Exact reconstruction on labeled ERR1816980 data

A truth-labeled reconstruction experiment was performed using the 72,000 known original DNA sequences of ERR1816980.

The evaluation metric is **exact sequence recovery**: a reconstructed sequence is counted as recovered only when it is exactly identical to the corresponding original DNA sequence.

| Method | Exactly recovered original sequences | Exact recovery rate |
|---|---:|---:|
| Clover routing core | 48,703 / 72,000 | 67.64% |
| Clover + Reed | **71,823 / 72,000** | **99.75%** |

Reed therefore recovered:

- **23,120 additional original DNA sequences**
- an absolute improvement of **32.11 percentage points**

over using the Clover routing core directly.

These numbers describe **reconstruction accuracy**, not Clover clustering accuracy.

### 5.2 Full large-scale run

A full local PEAR-assembled ERR1816980 input containing:

```text
15,787,115 reads
```

was processed with 16 workers and WFA reconstruction.

Observed reconstruction statistics:

| Metric | Result |
|---|---:|
| Clover clusters \(C\) | 153,813 |
| Raw reads entering Reed \(M\) | 13,402,459 |
| Unique sequences after compression \(U\) | 3,963,397 |
| Pairwise global alignments \(A\) | 3,809,584 |
| Expected \(U-C\) | 3,809,584 |
| \(A = U-C\) | **Yes** |
| End-to-end wall-clock time | about **79.96 s** |
| Workers | 16 |

The large-scale run therefore confirms that the implementation follows the intended alignment-count relation:

\[
A=U-C
\]

and does not fall back to an all-pairs alignment path.

`C = 153,813` is the number of clusters produced by Clover. It is not the number of original DNA sequences and is not generated by Reed.

---

## 6. What Has Been Improved Compared with the Original Clover Workflow

The current implementation improves the reconstruction stage in several ways.

### 6.1 From one core sequence to multi-read evidence

A Clover routing core is only one observed read and can contain sequencing errors.

Reed instead combines evidence from multiple reads belonging to the same cluster and reconstructs a consensus.

The labeled ERR1816980 experiment shows the practical effect:

```text
routing core exact recovery: 67.64%
Reed exact recovery:         99.75%
```

### 6.2 Duplicate-aware computation

Repeated reads are compressed before global alignment.

On the full large-scale run:

```text
raw clustered reads M = 13,402,459
unique sequences U    =  3,963,397
```

Thus many repeated observations contribute through multiplicity without requiring repeated pairwise global alignments.

### 6.3 Avoidance of all-pairs multiple alignment

For a cluster with \(U_i\) unique sequences:

```text
all-pairs strategy:
U_i(U_i - 1) / 2 alignments

Reed:
U_i - 1 alignments
```

This is important for scalability when cluster coverage increases.

### 6.4 Explicit handling of insertion and deletion evidence

Reed does not rely on direct position-wise voting over unaligned reads.

Global alignments are projected to a shared profile containing base, deletion and insertion evidence before consensus calling.

### 6.5 Independent reconstruction switch

Clover clustering can still be executed without Reed.

Reconstruction is enabled with:

```bash
--reconstruct
```

The global alignment backend can be selected independently:

```bash
--reconstruct-backend wfa
--reconstruct-backend edlib
--reconstruct-backend nw
```

This keeps clustering and reconstruction as separate algorithmic stages.

---

## 7. Current Limitations

Reed is currently a research prototype and still has several important limitations.

### 7.1 Dependence on Clover cluster quality

Reed assumes that a cluster contains reads from the same source sequence.

If Clover mixes reads from different original strands, the shared profile may produce an incorrect consensus.

Reed currently reconstructs **within** each Clover cluster and does not correct the upstream clustering assignment.

### 7.2 Cluster fragmentation is not solved

One original DNA strand may be split into multiple Clover clusters.

Reed currently reconstructs each cluster independently and does not merge fragmented clusters.

Therefore, the number of Reed consensus sequences can be larger than the true number of original DNA sequences.

### 7.3 Single-backbone bias

All pairwise alignments are projected onto one selected backbone.

If that backbone contains an unusual error pattern, the common coordinate system may be suboptimal.

The current backbone policies reduce this risk but do not eliminate it.

### 7.4 Difficult indel regions

Repeated bases, homopolymers and competing insertion patterns may produce ambiguous alignment representations.

A simple backbone-centered profile may not capture every equivalent indel representation optimally.

### 7.5 Consensus model is still relatively simple

The current method primarily combines alignment evidence and sequence multiplicity.

It does not yet use a complete probabilistic DNA-storage error model covering sequence context, platform-specific insertion/deletion bias, or all available base-quality information.

### 7.6 Validation is still incomplete

The strongest current exact-recovery result is from ERR1816980.

A publishable evaluation should include additional real DNA-storage datasets, simulated error-rate sweeps, different cluster sizes, and comparisons against alternative reconstruction/MSA approaches.

### 7.7 Alignment count is linear, but alignment cost is data-dependent

The number of alignments is \(U-C\), but the runtime of each WFA/NW/Edlib call depends on sequence length and divergence.

Very long sequences or high-error clusters may therefore have different scaling behavior from the current low-error DNA-storage setting.

---

## 8. Planned Algorithm Improvements

The next stage focuses on improving Reed without losing its low alignment-count property.

### Priority 1 — Better backbone selection without \(O(U^2)\)

Develop a more robust backbone score using inexpensive cluster statistics such as:

- multiplicity,
- sequence length,
- Clover routing-core information,
- k-mer support,
- approximate similarity sketches.

The goal is to improve the coordinate reference without constructing an all-pairs distance matrix.

### Priority 2 — Linear-complexity consensus refinement

Use the first Reed consensus as a refined backbone and perform at most one additional alignment pass.

For example:

```text
initial backbone
      |
      v
U-1 alignments
      |
      v
first consensus
      |
      v
selective second pass
      |
      v
final consensus
```

A fixed two-pass design would still require only a constant multiple of \(U\) alignments rather than \(O(U^2)\).

### Priority 3 — Confidence-aware reconstruction

Add confidence scores to profile positions and identify ambiguous regions.

Only low-confidence positions or suspicious reads would require additional processing.

This can concentrate computation where it is useful instead of increasing work for the entire cluster.

### Priority 4 — Better insertion/deletion modeling

Improve the insertion-slot representation and normalization of equivalent indel alignments.

For difficult local regions, a small local graph/POA-like refinement could be applied selectively while keeping the global algorithm backbone-centered.

### Priority 5 — Outlier detection

Detect reads that are inconsistent with the dominant cluster profile.

Possible signals include:

- unusually large edit distance to the backbone,
- low profile agreement,
- abnormal length,
- inconsistent insertion/deletion patterns.

This may reduce the influence of occasional misclustered or highly corrupted reads.

### Priority 6 — Quality- and error-model-aware consensus

Incorporate:

- base quality,
- position-dependent error rates,
- substitution/insertion/deletion priors,
- DNA synthesis/sequencing platform characteristics.

The objective is to replace simple evidence accumulation with better calibrated reconstruction confidence while preserving efficiency.

### Priority 7 — Broader complexity and accuracy experiments

Future experiments should explicitly vary:

```text
cluster size
coverage
sequence length
substitution rate
insertion rate
deletion rate
number of workers
```

and measure:

```text
exact reconstruction rate
edit distance
runtime
alignment count
memory
```

The most important scalability experiment is runtime versus cluster size, comparing the observed Reed trend with an all-pairs alignment baseline.

---

## 9. Installation

```bash
git clone https://github.com/oakmoss6-hash/Clover.git
cd Clover
git checkout reconstruction-v1

python3 -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
pip install -e .
```

The default large-scale alignment backend is WFA through `pywfa`.

---

## 10. Usage

### Clover clustering only

Run Clover without reconstruction:

```bash
python -m clover \
  -I input.fastq \
  -L 152 \
  -D 15 \
  -H 3 \
  -V 3 \
  -P 2
```

### Clover + Reed reconstruction

```bash
python -m clover \
  -I input.fastq \
  -L 152 \
  -D 15 \
  -H 3 \
  -V 3 \
  -P 2 \
  --reconstruct \
  --reconstruct-backend wfa \
  --reconstruct-backbone core \
  --consensus-output consensus.tsv
```

For `-P N`, Clover uses \(4^N\) worker partitions. For example:

```text
-P 0 -> 1 worker
-P 1 -> 4 workers
-P 2 -> 16 workers
```

---

## 11. Reconstruction Output

`consensus.tsv` contains one reconstruction result per Clover cluster.

Current output fields include information such as:

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

The global summary printed by the program includes:

```text
Reconstruction clusters: C
Reconstruction raw reads M: M
Reconstruction unique U: U
Pairwise global alignments A: A
U-C: U-C
```

For the current Reed implementation, a correct run should satisfy:

```text
A == U - C
```

---

## 12. Repository Structure

```text
clover/
├── __init__.py
├── __main__.py
├── input_io.py
├── load_config.py
├── main.py
├── reed.py
└── tree.py
```

Main responsibilities:

```text
tree.py
    Original Clover tree-based clustering structures.

main.py
    Clover clustering process and Reed integration.

reed.py
    Reed cluster-level multi-sequence global alignment
    and reconstruction algorithm.

input_io.py
    FASTA / FASTQ / Clover text input.

load_config.py
    Command-line and algorithm configuration.

__main__.py
    Command-line entry point.
```

---

## 13. Research Scope

The current work focuses on the following question:

> Given a set of reads that Clover has already assigned to the same high-homology cluster, how can their differences be identified and their common source sequence reconstructed accurately without using an expensive all-pairs multiple-alignment strategy?

Reed should therefore be viewed as a **DNA-storage cluster reconstruction / multi-sequence global alignment method**, rather than a replacement for BLAST, Bowtie, or general read mapping.

---

## 14. Original Clover

This project is based on the Clover DNA clustering algorithm:

> **Clover: tree structure-based efficient DNA clustering for DNA-based data storage**  
> *Briefings in Bioinformatics*, 2022.  
> DOI: `10.1093/bib/bbac336`

The original Clover clustering code and algorithm remain the foundation of the clustering stage. Reed is an experimental reconstruction extension built on top of Clover-generated clusters.

Please also see the repository `LICENSE` for licensing information.

---

## 15. Current Conclusion

The current implementation demonstrates three main results:

1. **Cluster-level multi-read reconstruction works on large real DNA-storage data.**
2. **Exact recovery on the current labeled ERR1816980 evaluation improves from 67.64% using Clover routing cores to 99.75% using Reed consensus.**
3. **The number of pairwise global alignments is exactly \(U-C\), avoiding an all-pairs \(O(U^2)\) alignment-count design.**

The next research stage is not to increase code complexity, but to improve backbone robustness, indel handling, confidence modeling and experimental validation while preserving the low alignment-count property.
