#!/bin/bash
# bracken-build with kraken2: many library files, and a library with no FASTA file.

ROOT=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

PREFIX="$WORK/db/library/added/a_deliberately_long_file_name_for_padding_"
ARG_MAX=$(getconf ARG_MAX)
N=$(( ARG_MAX / (${#PREFIX} + 10) + 2000 ))

mkdir -p "$WORK/db/library/added" "$WORK/db/taxonomy" "$WORK/bin"
touch "$WORK/db/hash.k2d" "$WORK/db/taxonomy/nodes.dmp"
for i in $(seq 1 $N); do
    printf ">seq%d\nACGT" "$i" > "$PREFIX$i.fna"
done

# Stub kraken2: one output line per FASTA record; status 64 without input files.
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

rm -rf "$WORK/db/library" "$WORK/db/database.kraken" "$WORK/db/database.kraken.tmp"
mkdir -p "$WORK/db/library/added"
PATH="$WORK/bin:$PATH" "$ROOT/bracken-build" -d "$WORK/db" -k 31 -l 100 -t 1 > "$WORK/log2.txt" 2>&1
status=$?
if [ $status -eq 0 ] || grep -q "Finished creating database.kraken" "$WORK/log2.txt"; then
    echo "FAIL: bracken-build reported success (status $status) for a library without FASTA files"
    exit 1
fi
echo "PASS: bracken-build stops when the library has no FASTA files"
