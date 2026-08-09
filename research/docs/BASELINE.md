# Phase 3.5 Baseline

## Source

- source: Guanjinqu/Clover GitHub ZIP
- upstream official commit: unknown
- local baseline commit: 127af066ae147ebd9e480bc81cdb91ed69f0f180
- branch: master
- OS: Linux DESKTOP-VMTF55G 6.18.33.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 18 21:54:43 UTC 2026 x86_64 GNU/Linux
- Python: Python 3.12.3
- pwd: /home/oakmoss/clover-research/Clover

## Git Status

`git status --short` reported:

```text
?? clover/__pycache__/
?? tests/__pycache__/
```

The requested baseline was expected to be clean, but the current repository contains untracked Python bytecode cache directories under upstream paths. No source file under `clover/`, `tests/`, or `example/` was modified by this audit.

## Required Files

All required files are present:

- `clover/align.py`
- `clover/main.py`
- `clover/tree.py`
- `tests/test_align.py`
- `tests/test_clust.py`
- `tests/test_tree.py`
- `example/example_tag_data.txt`
- `example/example_index_data.txt`

## Upstream SHA256

Saved to:

```text
research/results/upstream_files.sha256
```

## Original Unittest

Command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p 'test_*.py' -v
```

Full output:

```text
research/results/original_unittest.txt
```

Result:

- discovered: 7
- passed: 6
- failed: 0
- errors: 1
- skipped: 0

Failure category: C. missing dependency

Reason: `tests/test_clust.py` imports `clover.main`, and `clover/main.py` imports `tqdm`. The current environment does not have `tqdm`; per phase rules, no package was installed.
