#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_DIR="tmp/smoke_align_keep"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

python main.py align -f data/arp23.fasta -o "$TMP_DIR/trimmed.aln" -c 1 --keep

test -s "$TMP_DIR/trimmed.aln"
test -s "$TMP_DIR/trimmed.aln.untrimmed"

echo "align keep smoke test passed"
