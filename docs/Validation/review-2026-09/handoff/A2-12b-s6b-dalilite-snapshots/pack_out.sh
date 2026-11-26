#!/usr/bin/env bash
# pack_out.sh OUTDIR ZIPFILE : zip OUTDIR without the third-party archives it downloaded (OUTDIR/found: academic-licence / third-party code, public repository, GitHub 100 MB limit; their name, size and sha256 are in
# found_index.tsv and report.txt, their text inspection in inspect/), without .tmp / .page partial files and without any file of 3 MB or more.
set -uo pipefail
out="$1"; zipf="$2"
rm -f "$zipf" "${zipf%.zip}.tar.gz"
list=$(mktemp)
find "$out" -type f -size -3000k ! -path "$out/found/*" ! -name '.tmp*' ! -name '.page' | LC_ALL=C sort > "$list"
if command -v zip > /dev/null; then zip -q "$zipf" -@ < "$list"; else tar czf "${zipf%.zip}.tar.gz" -T "$list"; fi
rm -f "$list"
