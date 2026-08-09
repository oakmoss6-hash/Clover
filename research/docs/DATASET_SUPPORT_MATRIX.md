# Dataset Support Matrix

| Experiment | Support | Reason |
|---|---|---|
| upstream regression | PARTIALLY_SUPPORTED | Original tests exist, but one test module cannot run because `tqdm` is missing. |
| original Clover behavior | PARTIALLY_SUPPORTED | Core files and examples exist; full run blocked by missing dependency. |
| pairwise substitution correctness | SUPPORTED | `tests/test_align.py` covers one same-length substitution case. |
| pairwise insertion correctness | NOT_SUPPORTED | No test coverage and implementation does not handle insertion. |
| pairwise deletion correctness | NOT_SUPPORTED | No test coverage and shorter `read_2` can raise `IndexError`. |
| multi-read consensus correctness | NOT_SUPPORTED | Example tags contain exact duplicate pairs, not within-tag variants for consensus recovery. |
| IDS robustness | NOT_SUPPORTED | No observed insertion/deletion/substitution mixture in same-tag examples. |
| PCR duplicate robustness | PARTIALLY_SUPPORTED | Exact duplicates are observed, but only tiny examples. |
| contamination robustness | NOT_SUPPORTED | No observed same-tag contamination in examples. |
| runtime versus U | NOT_SUPPORTED | Example has only 3 unique sequences and tiny clusters. |
| memory versus U | NOT_SUPPORTED | Dataset too small. |
| complexity exponent alpha | NOT_SUPPORTED | Cannot fit `log(T) = alpha * log(U) + c` from this example. |
| real DNA reconstruction accuracy | UNKNOWN | Excel workbook was not parsed; examples are toy data. |

## Controlled U Scaling

The GitHub example cannot support controlled experiments for:

```text
U = 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000
```

It has only 3 unique sequences in the audited examples. Simple duplication of example reads cannot be treated as a real scaling dataset because it changes multiplicity without increasing unique read diversity or cluster difficulty.
