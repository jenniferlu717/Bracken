#!/bin/bash
# Regression test for bracken-build with kraken2 (issue #289, "Argument list too long").
# Case 1: a library with more FASTA files than fit on one command line.
# Case 2: a library with no FASTA file must stop the build with an error.
# A stub kraken2 stands in for the real one.
# Usage: tests/test_many_library_files.sh
# Exit status 0 on success, 1 on failure.

ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Enough files that their paths cannot fit in ARG_MAX bytes.
PREFIX="$WORK/db/library/added/a_deliberately_long_file_name_for_padding_"
ARG_MAX=$(getconf ARG_MAX)
N=$(( ARG_MAX / (${#PREFIX} + 10) + 2000 ))

mkdir -p "$WORK/db/library/added" "$WORK/db/taxonomy" "$WORK/bin"
touch "$WORK/db/hash.k2d" "$WORK/db/taxonomy/nodes.dmp"
for i in $(seq 1 $N); do
    printf ">seq%d\nACGT" "$i" > "$PREFIX$i.fna"
done

# Stub kraken2: reads FASTA from the files named on its command line (also
# "/dev/stdin") and prints one classification line per record. Like kraken2 it
# stops with status 64 when it gets no input file.
cat > "$WORK/bin/kraken2" <<"STUB"
#!/bin/bash
files=()
while [ $# -gt 0 ]; do
    case "$1" in
        --db|--threads) shift 2 ;;
        *) files+=("$1"); shift ;;
    esac
done
if [ ${#files[@]} -eq 0 ]; then
    echo "Need to specify input filenames!" >&2
    exit 64
fi
for f in "${files[@]}"; do
    while read -r line || [ -n "$line" ]; do
        case "$line" in
            ">"*) printf "C\t%s\t1\t4\t1:1\n" "${line#>}" ;;
        esac
    done < "$f"
done
STUB
chmod +x "$WORK/bin/kraken2"

# bracken-build goes on to run kmer2read_distr, which this test does not need.
PATH="$WORK/bin:$PATH" "$ROOT/bracken-build" -d "$WORK/db" -k 31 -l 100 -t 1 > "$WORK/log.txt" 2>&1

if [ ! -s "$WORK/db/database.kraken" ]; then
    echo "FAIL: database.kraken was not created"
    grep -i "too long" "$WORK/log.txt"
    exit 1
fi
lines=$(wc -l < "$WORK/db/database.kraken")
if [ "$lines" -ne "$N" ]; then
    echo "FAIL: database.kraken has $lines lines, expected $N"
    exit 1
fi
echo "PASS: bracken-build passed $N library files to kraken2"

# Case 2: no FASTA file in the library.
rm -rf "$WORK/db/library" "$WORK/db/database.kraken" "$WORK/db/database.kraken.tmp"
mkdir -p "$WORK/db/library/added"
PATH="$WORK/bin:$PATH" "$ROOT/bracken-build" -d "$WORK/db" -k 31 -l 100 -t 1 > "$WORK/log2.txt" 2>&1
status=$?
if [ $status -eq 0 ] || grep -q "Finished creating database.kraken" "$WORK/log2.txt"; then
    echo "FAIL: bracken-build reported success (status $status) for a library without FASTA files"
    exit 1
fi
echo "PASS: bracken-build stops when the library has no FASTA files"
