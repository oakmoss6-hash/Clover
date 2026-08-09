# Global Alignment Flow

## Files Read

- `clover/main.py`
- `clover/align.py`
- `clover/load_config.py`

## Call Sites

`clover/main.py` imports the module as:

```python
from clover import align as ag
```

`ag.global_align(self.ref_list[align_list[0]], dna_str)` is called inside `MyProcess.cluster` after a read matches an existing core sequence through one of the trie retrieval paths:

- front tree `a_tree`
- back tree `b_tree`
- middle tree `c_tree` / optional `d_tree`

## Flow

```text
input read
  -> MyProcess.cluster(read)
  -> split into dna_tag / dna_str
  -> tree retrieval against core-sequence tries
  -> existing core match selected
  -> optional global_align(core_sequence, dna_str)
  -> positional error counts in ref_error_dict
  -> possible core base correction by replacing one position
```

## Answers

1. Global alignment occurs after a cluster/core match, not before clustering and not as the primary clustering decision.

2. It is optional. The switch is `align_fuc` / `self.align_swicth`, default `False` in `clover/load_config.py`.

3. If global matching is off, Clover clustering still proceeds through trie retrieval and core-set updates.

4. For one successfully clustered read, the code calls `global_align` at most once along the chosen successful tree path.

5. It compares the read only with the selected core sequence in `self.ref_list[align_list[0]]`.

6. There is no read-vs-every-read-in-cluster alignment in the audited code.

7. Alignment output participates in:

- error marking: counts mismatch positions in `self.ref_error_dict[core_id]`.
- core correction: if a position exceeds support thresholds, `self.ref_list[core_id]` is modified at that single base.
- candidate sequence: no full candidate alignment path is stored.

8. Relevant variables:

- file: `clover/main.py`
- function: `MyProcess.cluster`
- `dna_str`: current read sequence
- `dna_tag`: read tag/index depending on mode
- `a_align`, `b_align`, `fin_align`, `align_list`: selected core id plus tree drift distance
- `self.ref_list`: core sequence storage used only when alignment mode stores core sequences
- `self.ref_dict`: core id to tag/index/member labels depending on mode
- `self.ref_error_dict`: per-core positional mismatch counts
- `self.index_list`: output mapping in non-virtual modes
