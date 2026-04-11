#!/usr/bin/env bash
set -euo pipefail

if ! command -v pre-commit &> /dev/null; then
    echo "pre-commit not found."
    echo "Make sure the correct environment is active."
    exit 1
fi

echo "[1/2] Install git hooks:"
pre-commit install

echo "[2/2] Install hooks in submodules:"

git submodule foreach '
  if [ -f .pre-commit-config.yaml ]; then
    echo "-> installing hooks in $name"
    pre-commit install
  fi
'

echo "Done"
