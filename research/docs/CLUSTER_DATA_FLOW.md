# Cluster Data Flow Audit

## 1. Core Sequence Set

Clover represents core sequences through trie indexes plus dictionaries in `MyProcess`.

When `align_fuc` is true and a read becomes a new core, the code stores:

```python
self.ref_list[dna_num] = dna_str
```

The tries store slices of that core sequence and map them to `dna_num`.

## 2. `ref_list`

- type: initialized as `{}` despite the name "list"
- key/index: `dna_num`, an integer assigned from the sequential processed read counter
- value: core sequence string, but only populated when `config_dict['align_fuc'] == True`

## 3. `ref_dict`

- type: dict
- key: core id / `dna_num`
- value: list of tags or indexes assigned to that core

In virtual/tag mode it appends `dna_tag`; in non-virtual mode it is used to map the matched core back to the first stored value.

## 4. Existing `cluster_id -> all reads`

NO. The current in-memory state stores tags/indexes for matched members, not the full original read strings for every cluster.

## 5. Whether the matched read remains in memory

In fast mode, input lines are preloaded in `data_dict` before processing, but `MyProcess.cluster` does not store all read sequences by cluster. In streaming/low-memory modes, reads are processed line by line and are not retained as cluster member sequences.

## 6. Low-memory mode

`--low` sets `mmr_mode=True`, `fast_mode=False`, `align_fuc=False`, `Statistical_model=False`, and `Virtual_mode=False`. It streams from `input_path`, writes output chunks when requested, and does not keep full cluster member reads.

## 7. Final Cluster Output

If `output_file` is configured:

- non-low-memory path writes `new_count_dict['index_list']` to `config_dict['output_file']`.
- low-memory path writes per-process files named `self.name + '_' + output_file`.

The output is a string representation of tuples from `self.index_list`, such as `(dna_index, core_label)`. It does not include full read sequence strings or core sequence strings.

## 8. Can current output recover `cluster_id -> list of reads`?

Partially, only if the original input file is available and output indexes can be joined back to input records. The output alone is insufficient because it does not include read sequences or core sequences.

## 9. Minimal Invasive Options

| Option | Changes clustering decision | Memory overhead | Complexity | Reproducibility | Assessment |
|---|---|---|---|---|---|
| A. modify clustering internal state | Risky | Higher | Medium | Harder to isolate | Not recommended now |
| B. collect at output/export layer | No if passive | Moderate | Medium | Good | Recommended |
| C. read cluster output after Clover completes | No | Low | Low to medium | Good, but requires joining original input | Viable fallback |
| D. wrapper / adapter | No | Low to moderate | Medium | Good | Recommended as outer interface |

Current recommendation: prioritize cluster export/interface work before pairwise aligner work, because Clover does not directly provide `cluster_id`, `core_sequence`, and `list_of_reads` as a complete reconstruction input.
