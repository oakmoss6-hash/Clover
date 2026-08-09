# Alignment Baseline Audit

## Files Read

- `clover/align.py`
- `tests/test_align.py`

## Answers

1. Global alignment function name: `global_align`.

2. Signature:

```python
def global_align(read_1, read_2):
```

3. Roles:

- `read_1`: the first compared sequence; docstring says it is in the core set.
- `read_2`: the post-match sequence/read being compared to the core sequence.

4. Return data structure: a Python list of tuples.

5. Return value semantics: each tuple is `(position, base)` where `position` is a zero-based mismatch position and `base` is the base from `read_2` at that position. It does not return an alignment path, edit distance, CIGAR, or operation labels.

6. Algorithmic machinery:

- dynamic programming matrix: NONE
- Needleman-Wunsch: NONE
- Smith-Waterman: NONE
- Levenshtein: NONE
- CIGAR: NONE

7. Supported error behavior:

- match: supported implicitly by emitting no tuple when bases match.
- substitution: supported only for same-position mismatch while iterating over `range(len(read_1))`.
- insertion: unsupported. If `read_2` has extra trailing bases, they are ignored; internal insertions shift all downstream positions and are reported as positional substitutions until the loop ends.
- deletion: unsupported. If `read_2` is shorter than `read_1`, the function can raise `IndexError`.

8. Unequal-length behavior from minimal experiment:

```text
AAA ATA -> [(1, 'T')]
AAA AAAA -> []
AAAA AAA -> IndexError: string index out of range
```

9. `tests/test_align.py` coverage:

- match: not explicitly covered
- substitution: covered by `AAA` versus `ATA`
- insertion: not covered
- deletion: not covered

## Conclusion

The original `global_align` is a positional mismatch marker, not a global alignment algorithm for IDS errors.
