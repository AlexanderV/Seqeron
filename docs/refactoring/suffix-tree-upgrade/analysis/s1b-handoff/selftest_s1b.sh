#!/bin/bash
# Self-test of the hand-off package scripts (work package 1b.9f): runs the entry point s1b.sh at 10^6 and 2x10^7 symbols natively (no Docker) and checks its outputs and its failure modes: a missing tool, too little
# disk, a bind-mount-like file system, a killed build and its resume (Resume on and off), a second run of the same text (identical container), a doctored form. Needs the published program:
#   S1B_BENCH=<SuffixTree.Esa.Bench.dll> S1B_EMPROTO=<SuffixTree.EmProto.dll> [S1B_TEST_DIR=<work directory>] selftest_s1b.sh [--with-tests]
# --with-tests also runs the Persistent test set through the test driver (about 40 s, needs the SDK and the repository checkout).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
S1B=$HERE/s1b.sh
W=${S1B_TEST_DIR:-$(mktemp -d)}
D=$W/data
: "${S1B_BENCH:?set S1B_BENCH}" "${S1B_EMPROTO:?set S1B_EMPROTO}"
export S1B_BENCH S1B_EMPROTO S1B_HEARTBEAT=0
PASS=0; FAIL=0
ok()   { PASS=$((PASS + 1)); echo "  PASS  $*"; }
bad()  { FAIL=$((FAIL + 1)); echo "  FAIL  $*"; }
expect() { # name condition-exit-code
  if [ "$2" = 0 ]; then ok "$1"; else bad "$1"; fi
}
run() { "$S1B" "$@" --data "$D" > "$W/last.out" 2> "$W/last.err"; echo $? > "$W/last.rc"; LAST_RC=$(cat "$W/last.rc"); }
has() { grep -q -- "$1" "$W/last.out" "$W/last.err" 2> /dev/null; }
sha_of_run() { cut -d' ' -f1 "$D/runs/$1/container.sha256"; }
rm -rf "$D"; mkdir -p "$D"

echo "1. syntax and the unit tests of the Python parts"
for f in "$HERE"/*.sh; do bash -n "$f"; expect "bash -n $(basename "$f")" $?; done
for f in "$HERE"/*.py; do python3 -I -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f"; expect "python syntax $(basename "$f")" $?; done
for t in test_s1b_check.py test_s1b_sampler.py test_s1b_gen.py; do
  if [ "$t" = test_s1b_gen.py ] && [ ! -f "$HERE/../em/make_1b9e_inputs.py" ]; then echo "  SKIP  $t needs the in-memory generators of the repository (analysis/em, analysis/bench): not in the image"; continue; fi
  python3 -I "$HERE/$t" > "$W/$t.log" 2>&1; expect "unit tests $t ($(grep -E '^Ran ' "$W/$t.log"))" $?
done

echo "2. help, plan, dry run"
run help; expect "help prints the commands" "$([ "$LAST_RC" = 0 ] && grep -q 'chain' "$W/last.out" && echo 0 || echo 1)"
run plan --n 3g; expect "plan prints the three sizes" "$(grep -q '3000000000' "$W/last.out" && grep -q '1000000000' "$W/last.out" && grep -q 'ESTIMATE' "$W/last.out" && echo 0 || echo 1)"
run chain --n 1m --dry-run; expect "chain --dry-run exits 0" "$LAST_RC"
expect "chain --dry-run creates no text and no container" "$([ ! -d "$D/text" ] && [ ! -d "$D/work" ] && echo 0 || echo 1)"
has "e8-gen" && has "sesa resume build" && has "handoff queries"; expect "chain --dry-run lists the steps" $?

echo "3. the chain at 10^6 on four kinds of text"
for kind in genome dna prot rep; do
  run chain --n 1m --kind "$kind" --structure yes --psi true
  expect "chain $kind 1m exits 0" "$LAST_RC"
  has "RESULT: PASS"; expect "chain $kind 1m: the checker says RESULT: PASS" $?
  if grep -q "^$kind	1000000	psi" "$HERE/reference_hashes.tsv"; then has "equals the reference table"; expect "chain $kind 1m: the container equals the reference container of the table" $?; fi
  [ -f "$D/runs/$kind-n1000000/build.1.sampler.json" ]; expect "chain $kind 1m: the sampler summary exists" $?
done
expect "the form names four runs" "$(python3 -I -c "import json;print(0 if len(json.load(open('$D/form.json'))['runs'])==4 else 1)")"

echo "4. a second run of the same text gives the same container; the form carries the hash"
run chain --n 1m --kind genome --tag b --structure no
expect "the second genome run exits 0" "$LAST_RC"
expect "the two genome containers are identical" "$([ "$(sha_of_run genome-n1000000)" = "$(sha_of_run genome-n1000000-b)" ] && echo 0 || echo 1)"
has "the other run genome-n1000000"; expect "the checker compared the two runs" $?
run check
has "Structure level not run" || has "NOT RUN S1" || has "NOT RUN"; expect "a skipped Structure level is reported as NOT RUN, not as a pass" $?

echo "5. a doctored form fails"
cp "$D/form.json" "$W/form.good"
python3 -I - "$D/form.json" <<'PY'
import json, sys
f = json.load(open(sys.argv[1])); r = f['runs']['dna-n1000000']
r['rss']['in_process_peak_mb'] = 400.0
json.dump(f, open(sys.argv[1], 'w'))
PY
run check; expect "a doctored RssAnon: exit 1" "$([ "$LAST_RC" = 1 ] && echo 0 || echo 1)"; has "FAIL.*R1"; expect "  and R1 says FAIL" $?
cp "$W/form.good" "$D/form.json"
python3 -I - "$D/form.json" <<'PY'
import json, sys
f = json.load(open(sys.argv[1])); f['runs']['genome-n1000000-b']['container']['sha256'] = '0' * 64
json.dump(f, open(sys.argv[1], 'w'))
PY
run check; expect "a doctored container hash: exit 1" "$([ "$LAST_RC" = 1 ] && echo 0 || echo 1)"; has "FAIL.*I1"; expect "  and I1 says FAIL" $?
cp "$W/form.good" "$D/form.json"
python3 -I - "$D/form.json" <<'PY'
import json, sys
f = json.load(open(sys.argv[1])); del f['runs']['prot-n1000000']['queries']
json.dump(f, open(sys.argv[1], 'w'))
PY
run check; expect "a missing section: exit 1" "$([ "$LAST_RC" = 1 ] && echo 0 || echo 1)"; has "missing section 'queries'"; expect "  and it names the section" $?
cp "$W/form.good" "$D/form.json"; run check; expect "the original form passes again" "$LAST_RC"

echo "6. failure modes of the pre-flight"
S1B_BENCH=/nonexistent/SuffixTree.Esa.Bench.dll run preflight --n 1m; expect "a missing program: pre-flight exits 3" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"; has "not found"; expect "  and says what is missing" $?
S1B_BENCH=/nonexistent/SuffixTree.Esa.Bench.dll run chain --n 1m --kind dna --tag nobench; expect "a missing program: the chain is refused (3)" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"
expect "  and nothing was generated" "$([ ! -f "$D/text/dna_n1000000.txt.nobench" ] && [ ! -d "$D/work/dna-n1000000-nobench" ] && echo 0 || echo 1)"
S1B_FREE_OVERRIDE=1000000000 run preflight --n 3g; expect "too little disk: pre-flight exits 3" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"; has "needs about"; expect "  and says how much" $?
S1B_FREE_OVERRIDE=1000000000 run chain --n 3g --kind genome --tag nodisk; expect "too little disk: the chain is refused (3)" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"
expect "  and no 3 GB text was written" "$([ ! -f "$D/text/genome_n3000000000.txt" ] && echo 0 || echo 1)"
mkdir -p "$W/fakeproc"; cp /proc/meminfo /proc/cpuinfo "$W/fakeproc/" 2> /dev/null
{ real=$(readlink -f "$D"); echo "virtiofs $real virtiofs rw 0 0"; } > "$W/fakeproc/mounts"
E8_PROC="$W/fakeproc" run preflight --n 1m; expect "a bind-mount-like file system: pre-flight exits 3" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"; has "bind mount"; expect "  and names the bind mount" $?
S1B_ALLOW_BADFS=1 E8_PROC="$W/fakeproc" run preflight --n 1m; expect "  S1B_ALLOW_BADFS=1 lets it through" "$LAST_RC"

echo "7. a killed build and its resume (2x10^7 symbols, genome-like, about 6 s)"
run gen --n 20m --kind genome
run build --n 20m --kind genome --tag ref --resume off
expect "the uninterrupted build exits 0" "$LAST_RC"; run hash --n 20m --kind genome --tag ref
REF=$(sha_of_run genome-n20000000-ref)
run build --n 20m --kind genome --tag kon --resume on --kill-after 2.5
if [ "$LAST_RC" = 4 ]; then
  ok "Resume on: the killed build exits 4"; has "Resume is ON"; expect "  and tells to run the same command again" $?
  [ "$(cat "$D/runs/genome-n20000000-kon/build.1.rc")" = 137 ]; expect "  exit code 137 recorded" $?
  run build --n 20m --kind genome --tag kon --resume on
  expect "  the second call resumes and exits 0" "$LAST_RC"
  grep -aP "^R\tresume\tstate_found\tyes" "$D/runs/genome-n20000000-kon/build.2.out" > /dev/null; expect "  the state of the killed attempt was found" $?
  run hash --n 20m --kind genome --tag kon
  expect "  the resumed container equals the uninterrupted one" "$([ "$(sha_of_run genome-n20000000-kon)" = "$REF" ] && echo 0 || echo 1)"
else
  echo "  SKIP  the kill landed after the end of the build (this machine is faster than the test assumes): rc $LAST_RC"
fi
# the pre-flight of a resumed chain counts what the interrupted attempt already holds: with just less free space than a whole run needs (2.54 GB at 2x10^7) the chain must still go on
run build --n 20m --kind genome --tag kfree --resume on --kill-after 2.5
if [ "$LAST_RC" = 4 ]; then
  S1B_FREE_OVERRIDE=2400000000 run chain --n 20m --kind genome --tag kfree --resume on --structure no
  expect "a resumed chain counts the files of the killed attempt: it is not refused although a whole run would not fit" "$LAST_RC"
  S1B_FREE_OVERRIDE=2400000000 run preflight --n 20m --kind genome
  expect "  while a pre-flight that does not know them refuses (3)" "$([ "$LAST_RC" = 3 ] && echo 0 || echo 1)"
else
  echo "  SKIP  the kill landed after the end of the build: rc $LAST_RC"
fi
run build --n 20m --kind genome --tag koff --resume off --kill-after 2.5
if [ "$LAST_RC" = 4 ]; then
  ok "Resume off: the killed build exits 4"; has "Resume is OFF"; expect "  and says that the next call starts again" $?
  run build --n 20m --kind genome --tag koff --resume off
  expect "  the second call starts again and exits 0" "$LAST_RC"; run hash --n 20m --kind genome --tag koff
  expect "  and gives the same container" "$([ "$(sha_of_run genome-n20000000-koff)" = "$REF" ] && echo 0 || echo 1)"
else
  echo "  SKIP  the kill landed after the end of the build: rc $LAST_RC"
fi

if [ "${1:-}" = --with-tests ]; then
  echo "8. the test driver on the Persistent set"
  S1B_SRC=${S1B_SRC:-$(git -C "$HERE" rev-parse --show-toplevel)} run tests persistent
  expect "tests persistent exits 0" "$LAST_RC"
  python3 -I -c "import json;j=json.load(open('$D/tests/persistent.json'));print(j);raise SystemExit(0 if j['status']=='pass' and j['passed']>600 and j['failed']==0 else 1)"; expect "the result file has the totals" $?
  run form; run check; has "T1 persistent"; expect "the checker reads the test set" $?
fi

echo
echo "SELFTEST: $PASS passed, $FAIL failed (work directory $W)"
[ "$FAIL" = 0 ]
