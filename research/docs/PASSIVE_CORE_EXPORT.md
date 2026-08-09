# Passive Initial Clover Core Export

Phase 4.5 exports the initial Clover core: the exact input read sequence at the moment Clover creates a new core. This is distinct from any corrected core sequence that may later appear in `ref_list` after original Clover global correction.

Schema:

```text
cluster_id<TAB>core_read_id<TAB>core_sequence
```

The export is opt-in through `--export-cluster-cores PATH`. Default behavior remains disabled.

The implementation writes at the new-core creation event in `MyProcess.cluster`, immediately after `self.ref_dict[dna_num]=[dna_tag]` and before tree insertion. It exports `dna_tag`, `dna_tag`, and `dna_str`.

Multiprocessing policy: first version is single-process only. If core export is enabled with `N_PROCESS != 1`, Clover fails fast with a clear `ValueError`.
