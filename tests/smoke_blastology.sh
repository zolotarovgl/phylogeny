#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_ROOT="tmp/smoke_blastology"
TEMP_DIR="$OUT_ROOT/tmp"
RESULT_DIR="$OUT_ROOT/results"
OUTFILE="$OUT_ROOT/Owenia_bcl2.tab"

python main.py blastology \
  --query data/BCL2.fasta \
  --refnames data/BCL2.names \
  --target data/sample.fasta \
  -c 2 \
  --force \
  --soi Owefus \
  --outputfile "$OUTFILE" \
  --phymethod fasttree \
  --output_dir "$RESULT_DIR" \
  --temp_dir "$TEMP_DIR"

test -s "$OUTFILE"
test "$(wc -l < "$OUTFILE")" -eq 3

for gene in \
  Owefus_OFUSG13935.2 \
  Owefus_OFUSG16636.1 \
  Owefus_OFUSG14207.1
do
  grep -q "^${gene}[[:space:]]" "$OUTFILE"
done

awk -F'\t' 'NF != 5 { exit 1 }' "$OUTFILE"

test -s "$RESULT_DIR/clusters/search.HG0.align.log"
test -s "$RESULT_DIR/clusters/search.HG0.phylogeny.log"
test -s "$RESULT_DIR/clusters/search.HG0.possvm.log"

echo "blastology smoke test passed: $OUTFILE"
