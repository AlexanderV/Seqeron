#!/usr/bin/env bash
# Runs inside the image: every non-comment line of rows.txt is run (60 s limit) in its own directory; standard output, error stream, exit status and written files go to /work/out/rows/NNN/.
set -u
cd /work
rm -rf out; mkdir -p out/rows
{ echo "uname -m: $(uname -m)"; seqret -version 2>&1 | head -2; sha256sum /usr/share/EMBOSS/data/DRCAT.dat; date -u; getent hosts www.ebi.ac.uk | head -1; } > out/environment.txt
cp /usr/share/EMBOSS/data/DRCAT.dat out/DRCAT.dat 2>/dev/null || true
n=0
while IFS= read -r line; do
  case "$line" in ''|\#*) continue;; esac
  n=$((n+1)); d=out/rows/$(printf %03d $n); mkdir -p $d
  echo "$line" > $d/command.txt
  ( cd $d && timeout 60 bash -c "$line" > stdout.txt 2> stderr.txt; echo $? > status.txt )
done < rows.txt
echo "$n rows" | tee out/summary.txt
if [ -n "${HOST_UID:-}" ]; then chown -R "$HOST_UID:$HOST_GID" out 2>/dev/null || true; fi
