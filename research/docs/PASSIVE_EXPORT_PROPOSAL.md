# Passive Core Export Proposal

## Problem

Current Clover membership output does not include reliable core sequences. The research adapter must not infer core sequence from cluster members because Clover core and reconstruction backbone are different concepts.

## Minimal Passive Export

Add an explicit opt-in flag, conceptually:

```text
--export-cluster-cores path.tsv
```

The file should contain:

```text
core_id<TAB>core_sequence
```

or, if Clover exposes a stable cluster label separately:

```text
cluster_id<TAB>core_sequence
```

## Constraints

The export must not change:

- tree retrieval
- match threshold
- drift behavior
- cluster membership
- original `global_align`
- `ref_list` contents

It should only read existing core state and write an additional file when the explicit flag is provided.

## Why This Is Minimal

The adapter can already join original input records with membership output. The missing piece is a reliable mapping from cluster/core id to core sequence. A passive export provides that without changing clustering decisions.
