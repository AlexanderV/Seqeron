#!/usr/bin/env bash
# Runs inside the image (cwd /work, the package directory mounted).  Everything goes to ./out; nothing here uses account data.
set -u
cd /work
rm -rf out; mkdir -p out/sets out/direct out/http
{
  echo "uname -m: $(uname -m)"; echo "date (UTC): $(date -u +%FT%TZ)"
  echo "mafft-src sha256: $(cat /opt/mafft-src.sha256)"
  go version
  ls -l /usr/local/libexec/mafft/dash_client 2>&1; sha256sum /usr/local/libexec/mafft/dash_client 2>&1
  mafft --version 2>&1 | head -1
  [ -e /opt/DASH_CLIENT_BUILD_FAILED ] && echo "DASH_CLIENT_BUILD_FAILED (see mafft.build.log)"
} > out/versions.txt 2>&1
cp /opt/mafft.build.log out/mafft.build.log 2>/dev/null
[ -e /opt/DASH_CLIENT_BUILD_FAILED ] && { cat out/versions.txt; echo "The image has no dash_client: nothing to run." >&2; exit 1; }
python3 make_inputs.py inputs > /dev/null
PORT=8765
LOG=/work/out/http/exchanges.jsonl
: > "$LOG"
python3 dashproxy.py $PORT "$LOG" &
PROXY=$!
sleep 2
SERVER="http://127.0.0.1:$PORT/dash/REST1.0/"
SERVERNOSLASH="http://127.0.0.1:$PORT/dash/REST1.0"
export TMPDIR=/tmp
# (1) the real dash_client, on every set, exactly as the script calls it:  dash_client -url URL -i dashin -sequences dashsequences -hat3 hat3.seed
while IFS=$'\t' read -r sid fn; do
  d=/tmp/direct-$sid; rm -rf "$d"; mkdir -p "$d"; cd "$d"
  awk '{if(/^>/)print; else {gsub( /-/,"" ); print;}}' /work/inputs/$fn > dashin
  mark=$(wc -l < "$LOG")
  timeout 1800 /usr/local/libexec/mafft/dash_client -url "$SERVER" -i dashin -sequences dashsequences -hat3 hat3.seed > stdout.txt 2> stderr.txt; echo $? > rc.txt
  tail -n +$((mark+1)) "$LOG" > exchanges.jsonl
  cp /work/inputs/$fn in.fa
  tar czf /work/out/direct/$sid.tar.gz in.fa dashin stdout.txt stderr.txt rc.txt dashsequences hat3.seed exchanges.jsonl 2> /dev/null
  cd /work
done < inputs/sets.tsv
# (2) the script mafft with --dash: the first option line on every set, every option line on 5 sets
SUB="single-ubiquitin family-00 mix-00 mix-03 globin-pair"
printf '1\t2\n3\t4\n' > /work/out/merge_table.txt
while IFS=$'\t' read -r oid opts; do
  while IFS=$'\t' read -r sid fn; do
    if [ "$oid" != "dash" ]; then case " $SUB " in *" $sid "*) ;; *) continue;; esac; fi
    o=${opts//@SERVER@/$SERVER}; o=${o//@SERVERNOSLASH@/$SERVERNOSLASH}; o=${o//@MERGE@/merge_table.txt}
    d=/tmp/run-$oid-$sid; rm -rf "$d"; mkdir -p "$d"; cd "$d"
    cp /work/inputs/$fn in.fa; cp /work/out/merge_table.txt .
    mark=$(wc -l < "$LOG")
    # shellcheck disable=SC2086
    timeout 3600 mafft $o in.fa > stdout.txt 2> stderr.txt; echo $? > rc.txt
    tail -n +$((mark+1)) "$LOG" > exchanges.jsonl
    echo "$o" > options.txt
    tar czf /work/out/sets/$oid.$sid.tar.gz in.fa options.txt stdout.txt stderr.txt rc.txt exchanges.jsonl $(ls in.fa.* GuideTree 2>/dev/null) 2> /dev/null
    cd /work
  done < inputs/sets.tsv
done < inputs/options.tsv
kill $PROXY 2> /dev/null
echo "done: $(ls out/direct | wc -l) direct runs, $(ls out/sets | wc -l) mafft runs, $(wc -l < "$LOG") HTTP exchanges"
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" /work/out /work/inputs 2> /dev/null
exit 0
