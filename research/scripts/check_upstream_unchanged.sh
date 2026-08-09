#!/usr/bin/env bash
set -euo pipefail

if ! git diff --exit-code -- clover tests example >/tmp/clover_upstream_diff.txt; then
    cat /tmp/clover_upstream_diff.txt
    git diff --name-only -- clover tests example
    exit 1
fi

untracked=$(git status --porcelain -- clover tests example)
if [ -n "$untracked" ]; then
    echo "$untracked"
    exit 1
fi

echo "UPSTREAM CLOVER FILES UNCHANGED"
