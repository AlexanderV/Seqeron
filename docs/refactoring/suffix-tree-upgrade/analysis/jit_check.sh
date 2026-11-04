#!/usr/bin/env bash
# T3 of the S1a review (section 6): the JIT'd, optimized code of the ESA kernel methods must contain no managed range-check helper
# (CORINFO_HELP_RNGCHKFAIL) and no call to the Buf<T> accessors (they must inline to raw loads and stores). Tier-0 listings are ignored (unoptimized by design);
# every optimized listing (Tier1, OSR, FullOpts) of every listed method is checked, and the run FAILS when there is none (a pattern that matches nothing proves nothing).
#   usage: jit_check.sh <input file> "<method names>" [bench args...]
#   e.g.   jit_check.sh data/dna_1m.txt "ComputePlcp WriteU8 WriteU16 CommonPrefix" --lcp-kernel pointer
set -eu
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$B/../../../.." && pwd)"
export DOTNET_ROOT=${DOTNET_ROOT:-$HOME/.dotnet} PATH=$HOME/.dotnet:$PATH DOTNET_CLI_TELEMETRY_OPTOUT=1
DLL=$ROOT/apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll
INPUT=$1; METHODS=$2; shift 2
OUT=$(mktemp)
DOTNET_JitDisasm="$METHODS" DOTNET_JitStdOutFile="$OUT" dotnet "$DLL" build --input "$INPUT" "$@" >/dev/null
python3 -I - "$OUT" <<'PY'
import re, sys
text = open(sys.argv[1]).read()
listings = re.split(r'(?m)^; Assembly listing for method ', text)[1:]
bad = checked = 0
for l in listings:
    head = l.split('\n', 1)[0]
    tier = re.search(r'\(([^)]*)\)\s*$', head)
    tier = tier.group(1) if tier else '?'
    if 'Tier0' in tier and 'OSR' not in tier and 'Tier1' not in tier:
        continue
    checked += 1
    rng = l.count('CORINFO_HELP_RNGCHKFAIL')
    calls = len(re.findall(r'call\s+\[(?:SuffixTree\.Esa\.)?Buf`1\[', l))  # a surviving call INTO the accessor struct (its indexer, loads, slices); a call that merely takes or returns a Buf (arena allocation, a helper method) is not an accessor call
    print(f'{head}: range-check helpers {rng}, accessor calls {calls}')
    if rng or calls: bad += 1
if checked == 0:
    print('FAIL: no optimized listing found (method names do not match, or the methods were never promoted to Tier1/OSR)'); sys.exit(2)
print(f'{checked} optimized listing(s) checked, {bad} with findings')
sys.exit(1 if bad else 0)
PY
status=$?
rm -f "$OUT"
exit $status
