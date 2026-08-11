Phase 13.0 — architecture refactor
==================================

Purpose
-------
Establish the final low-complexity reconstruction interfaces before adding
Clover tree-hit anchors, alignment elimination, and parallel scheduling.

This patch intentionally preserves historical results by default:
    backbone_policy="core"

New pieces
----------
1. ClusterState
   - streaming-friendly exact duplicate aggregation
   - stores sequence -> multiplicity once
   - separates Clover routing core from reconstruction state
   - reconstruct_state() can later be called directly from Clover

2. Backbone selector
   Policies:
   - core            : historical behavior / reproducibility baseline
   - support_length  : truth-blind abundance-first, weighted-median-length tie break
   - max_span        : experimental diagnostic policy

3. Reconstruction invariant
   Every unique non-backbone sequence is aligned once:
       pairwise_alignment_count == unique_sequence_count - 1

Not included yet
----------------
- Clover tree-hit/drift anchor fast path
- WFA skipping / certified direct path
- multiprocessing/MPI/GPU scheduler

Those should be added only after the Clover source-code interface audit and
the workload-reduction experiment, so the core method stays simple.

Install
-------
From the Clover repository root:

    unzip -o phase13_0_architecture_refactor.zip -d .
    python -m unittest discover -s research/tests -v

The existing reconstruct_cluster(...) API remains compatible.
