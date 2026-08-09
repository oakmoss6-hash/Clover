# Passive Core Export Complexity

For `C` created cores and strand length `L`, passive export writes one TSV row per core, so additional work is proportional to emitted bytes: `O(C * L)`.

Memory overhead is `O(1)` beyond the file handle. It performs no alignment, edit distance, consensus, benchmark, or all-pairs comparison.
