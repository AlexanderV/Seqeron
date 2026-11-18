#!/usr/bin/env bash
# Runs INSIDE the container (cwd /work = this package directory mounted from the host).
# Runs the whole ORFfinder matrix of make_jobs.py and writes everything into ./out/ :
#   out/ENV.txt  out/SHA256SUMS  out/runs.tsv  out/jobs.tsv  out/runs/<job id>/{cmd.txt,stdout.txt,stderr.txt,rc.txt,meta.tsv,orfs.txt}
# Environment: JOBS (parallel processes, default = CPU count), QUICK=1 (reduced matrix), SINGLE=0 (skip one-sequence runs),
#              FORCE=1 (re-run jobs that already have out/runs/<id>/rc.txt).
set -u
cd /work || exit 1
OUT=/work/out
mkdir -p "$OUT/runs" work/single
JOBS=${JOBS:-$(nproc)}
case "$JOBS" in ''|*[!0-9]*) JOBS=2 ;; esac

{
  echo "== date (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "== uname -a: $(uname -a)"
  echo "== uname -m: $(uname -m)"
  echo "== /etc/os-release:"; cat /etc/os-release
  echo "== nproc: $(nproc)   JOBS=$JOBS QUICK=${QUICK:-0} SINGLE=${SINGLE:-1}"
  echo "== checksums (as downloaded / unpacked in the image):"; cat /opt/orffinder/SHA256SUMS
  echo "== ldd:"; cat /opt/orffinder/ldd.txt
  echo "== ORFfinder -version:";      ORFfinder -version 2>&1; echo "rc=$?"
  echo "== ORFfinder -version-full:"; ORFfinder -version-full 2>&1; echo "rc=$?"
  echo "== python3: $(python3 --version 2>&1)"
} > "$OUT/ENV.txt" 2>&1
cp /opt/orffinder/USAGE.txt /opt/orffinder/CHANGELOG.txt "$OUT/"

# Is the committed input set what the generator produces? (determinism check; a mismatch is reported, not fatal)
python3 -I make_inputs.py --check > "$OUT/inputs_check.txt" 2>&1 || echo "WARNING: inputs differ from make_inputs.py output, see out/inputs_check.txt"
(cd inputs && sha256sum -c inputs.sha256 > "$OUT/inputs_sha256_check.txt" 2>&1) || echo "WARNING: inputs/ differs from inputs.sha256"

# One-sequence FASTA files for the per-process runs.
python3 -I - <<'PY'
import os
rec = None
with open("inputs/designed.fa") as fh:
    for line in fh:
        if line.startswith(">"):
            if rec: open("work/single/%s.fa" % rec[0], "w").write(">%s\n%s" % (rec[0], "".join(rec[1])))
            rec = (line[1:].strip(), [])
        else:
            rec[1].append(line)
if rec: open("work/single/%s.fa" % rec[0], "w").write(">%s\n%s" % (rec[0], "".join(rec[1])))
PY

python3 -I make_jobs.py > "$OUT/jobs.tsv" || { echo "make_jobs.py failed"; exit 1; }
TOTAL=$(wc -l < "$OUT/jobs.tsv")
echo "$TOTAL ORFfinder runs, $JOBS in parallel"

# Smoke test: if the binary cannot even start, stop with a clear message (and keep ENV.txt).
if ! ORFfinder -version > /dev/null 2>&1; then
  echo "ORFfinder does not start (see out/ENV.txt). On Apple Silicon use: docker run --platform linux/amd64 ..." >&2
  exit 2
fi

run_job() {
  local id="$1" group="$2" infile="$3" tmo="$4" useout="$5" args="$6"
  local dir="$OUT/runs/$id"
  if [ -z "${FORCE:-}" ] && [ -f "$dir/rc.txt" ]; then return 0; fi
  mkdir -p "$dir"; cd "$dir" || return 1
  local -a a=() cmd=()
  read -r -a a <<< "$args"
  cmd=(timeout -k 5 "$tmo" ORFfinder)
  [ "$infile" != "-" ] && cmd+=(-in "/work/$infile")
  cmd+=("${a[@]}")
  [ "$useout" = "1" ] && cmd+=(-out "$dir/orfs.txt")
  printf '%q ' "${cmd[@]}" > cmd.txt; echo >> cmd.txt
  local t0 t1 rc
  t0=$(date +%s%N)
  "${cmd[@]}" > stdout.txt 2> stderr.txt < /dev/null
  rc=$?
  t1=$(date +%s%N)
  echo "$rc" > rc.txt
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$group" "$rc" "$(( (t1 - t0) / 1000000 ))" \
      "$(stat -c %s orfs.txt 2>/dev/null || echo -)" "$(stat -c %s stdout.txt)" "$(stat -c %s stderr.txt)" > meta.tsv
}

n=0; done_n=0
while IFS=$'\t' read -r id group infile tmo useout args; do
  run_job "$id" "$group" "$infile" "$tmo" "$useout" "$args" &
  n=$((n + 1))
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  if [ $((n % 100)) -eq 0 ]; then echo "  started $n / $TOTAL  ($(date -u +%H:%M:%S))"; fi
done < "$OUT/jobs.tsv"
wait
echo "all $TOTAL runs finished"

printf 'id\tgroup\trc\tmillis\torfs_bytes\tstdout_bytes\tstderr_bytes\n' > "$OUT/runs.tsv"
cat "$OUT"/runs/*/meta.tsv >> "$OUT/runs.tsv"
echo "== rc histogram:"; tail -n +2 "$OUT/runs.tsv" | cut -f3 | sort | uniq -c | tee "$OUT/rc_histogram.txt"
( cd "$OUT" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS )
if [ -n "${HOST_UID:-}" ]; then chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" "$OUT" work 2>/dev/null || true; fi
echo done
