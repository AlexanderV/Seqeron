#!/usr/bin/env bash
# Runs inside the image: every row is `patmatmotifs -auto ARGS` on the inputs of /work/inputs with the real PROSITE data; standard output, error stream, exit status and the written files go to /work/out/rows/<n>.*
set -u
cd /work
rm -rf out; mkdir -p out/rows inputs
python3 -I make_inputs.py inputs; python3 -I make_rows.py inputs/rows.txt
{ echo "uname -m: $(uname -m)"; patmatmotifs -version 2>&1 | head -2; sha256sum /opt/prosite/prosite.dat /opt/prosite/prosite.doc /opt/emboss-data/PROSITE/* 2>/dev/null | head -20; date -u; } > out/environment.txt
cp /opt/prosite/SHA256SUMS out/PROSITE.SHA256SUMS; cp /opt/emboss-data/prosextract.log out/ 2>/dev/null || true
export EMBOSS_DATA=/opt/emboss-data:/usr/share/EMBOSS/data
n=0
while IFS= read -r args; do
  [ -z "$args" ] && continue
  n=$((n+1)); d=out/rows/$(printf %03d $n); mkdir -p $d; cp -r inputs/* $d/
  ( cd $d && eval "patmatmotifs -auto $args" > stdout.txt 2> stderr.txt; echo $? > status.txt )
  echo "$args" > $d/args.txt
  # valgrind verdict (a run that reads invalid memory has no defined result)
  ( cd $d && eval "valgrind -q --error-exitcode=99 patmatmotifs -auto $args" > /dev/null 2> valgrind.txt; echo $? > valgrind_status.txt )
done < inputs/rows.txt
echo "$n rows" | tee out/summary.txt
if [ -n "${HOST_UID:-}" ]; then chown -R "$HOST_UID:$HOST_GID" out inputs 2>/dev/null || true; fi
