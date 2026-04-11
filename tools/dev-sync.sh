#!/usr/bin/env bash
set -euo pipefail

echo "[1/5] Sync submodule metadata"
git submodule sync --recursive

echo "[2/5] Init/update to recorded commits"
git submodule update --init --recursive

echo "[3/5] Fetch all remotes"
git submodule foreach 'git fetch origin'

echo "[4/5] Attach to declared branches"

git submodule foreach '
  set -e

  branch=$(git config -f $toplevel/.gitmodules submodule.$name.branch)

  echo "-> $name -> $branch"

  # Create or switch branch
  if git show-ref --verify --quiet refs/heads/$branch; then
    git checkout $branch
  else
    git checkout -b $branch origin/$branch
  fi

  # Hard-align to remote (important for reproducibility of PGO builds)
  git reset --hard origin/$branch
'

echo "[5/5] Done: all submodules on correct branches"
