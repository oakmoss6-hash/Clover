# Original Datasets Audit

## Example Data

Detailed outputs:

- `research/results/example_dataset_audit.csv`
- `research/results/example_dataset_audit.txt`

Summary:

- `example/example_tag_data.txt`: 6 reads, 3 unique sequences, duplicate ratio 0.5, all length 150, alphabet ACGT.
- `example/example_index_data.txt`: 6 reads, 3 unique sequences, duplicate ratio 0.5, all length 150, alphabet ACGT.

The tag data contains 3 tags, each with 2 reads. Within each tag, reads are identical exact duplicates. No same-tag substitution-like variation and no different-length variation were observed.

## Excel Dataset

File:

```text
experiment result/raw_data.xlsx
```

- present: yes
- file size: 139758 bytes
- parsed: no
- reason: `openpyxl` is unavailable in the current Python environment, and this phase forbids installing dependencies.

Because the workbook was not parsed, this audit cannot determine whether it contains raw sequencing reads or paper experiment summaries/results.
