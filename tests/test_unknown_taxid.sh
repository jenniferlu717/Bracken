#!/bin/bash
# kmer2read_distr must not crash on a taxid missing from nodes.dmp.
# Usage: tests/test_unknown_taxid.sh [path/to/kmer2read_distr]

BIN=${1:-"$(dirname "$0")/../src/kmer2read_distr"}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

mkdir "$WORK/tax"
printf "1\t|\t1\t|\tno rank\t|\t\n10\t|\t1\t|\tgenus\t|\t\n11\t|\t10\t|\tspecies\t|\t\n12\t|\t10\t|\tspecies\t|\t\n" > "$WORK/tax/nodes.dmp"
printf "seqA\t11\n" > "$WORK/seqid2taxid.map"
# Taxid 99 is not in nodes.dmp
printf "C\tseqA\t11\t100\t99:1 11:1 12:1\n" > "$WORK/database.kraken"

"$BIN" --seqid2taxid "$WORK/seqid2taxid.map" --taxonomy "$WORK/tax" \
    --kraken "$WORK/database.kraken" --output "$WORK/out.tsv" \
    -k 31 -l 33 -t 1 > "$WORK/log.txt" 2>&1
status=$?
if [ $status -ne 0 ]; then
    echo "FAIL: kmer2read_distr exited with status $status"
    exit 1
fi

expected=$(printf "seqA\t11\t\t10:1 ")
actual=$(head -n 1 "$WORK/out.tsv")
if [ "$actual" != "$expected" ]; then
    echo "FAIL: unexpected output line: [$actual]"
    exit 1
fi
echo "PASS: taxid missing from nodes.dmp is handled"
