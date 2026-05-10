#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TMP_DIR="tmp/smoke_pfamscan"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

cat > "$TMP_DIR/toy_alignment.fasta" <<'EOF'
>toy_a
ACDEFGHIK
>toy_b
ACDEYGHIK
EOF

cat > "$TMP_DIR/toy_query.fasta" <<'EOF'
>toy_query_match
ACDEFGHIK
>toy_query_noise
MMMMMMMMM
EOF

hmmbuild --amino -n ToyDomain "$TMP_DIR/toy.hmm" "$TMP_DIR/toy_alignment.fasta" > "$TMP_DIR/hmmbuild.log"
hmmpress -f "$TMP_DIR/toy.hmm" > "$TMP_DIR/hmmpress.log"

python main.py pfamscan \
  --fasta "$TMP_DIR/toy_query.fasta" \
  --pfam_db "$TMP_DIR/toy.hmm" \
  --output_dir "$TMP_DIR/results" \
  --prefix toy \
  --no-cut_ga \
  --domE 1e-2

test -s "$TMP_DIR/results/toy.pfamscan.domtblout"
test -s "$TMP_DIR/results/toy.pfamscan.tsv"
test -s "$TMP_DIR/results/toy.pfamscan.domains.csv"
test -s "$TMP_DIR/results/toy.pfamscan_archs.csv"

grep -q '^sequence_id' "$TMP_DIR/results/toy.pfamscan.tsv"
grep -q '^toy_query_match[[:space:]]' "$TMP_DIR/results/toy.pfamscan.tsv"
grep -q $'^toy_query_match\t[0-9][0-9]*\t[0-9][0-9]*\tToyDomain$' "$TMP_DIR/results/toy.pfamscan.domains.csv"
grep -q $'^toy_query_match\t[0-9][0-9]*\t[0-9][0-9]*\tToyDomain$' "$TMP_DIR/results/toy.pfamscan_archs.csv"
! grep -q '/' "$TMP_DIR/results/toy.pfamscan_archs.csv"

echo "pfamscan smoke test passed"
