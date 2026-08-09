# Adapter Complexity

Let:

- `N` be the number of raw input reads.
- `C` be the number of clusters in the membership output.
- `B` be the total number of input bases.

## Time

- input parsing: `O(B)` because each input line and sequence is read once.
- membership parsing: `O(N)` for one membership tuple per assigned read.
- ClusterRecord construction: `O(N)` to group assigned reads and instantiate records.

No edit distance, NW, WFA, Edlib, profile, consensus, or all-pairs comparison is called in this phase.

## Memory

The current Phase 4 implementation is in-memory:

- read table: `O(N)` records plus sequence storage from the input file
- membership map: `O(N)` assigned reads
- grouped cluster members: `O(N)` references to `ReadRecord`
- cluster list: `O(C)`

Future streaming mode should avoid one huge nested list by consuming grouped membership/export shards and yielding `ClusterRecord` one cluster at a time.

## No Quadratic Operation

The adapter only reorganizes data. It does not compare cluster members pairwise and does not introduce `O(N^2)` behavior.
