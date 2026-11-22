#!/bin/bash
# Shared helpers of the hand-off package s1b-handoff (work package 1b.9f); sourced by s1b.sh and preflight.sh. The helpers of the E8 package (cgroup, memory, disk, file system, human) are reused, not copied.
# Nothing here sends anything anywhere.

S1B_HERE=${S1B_HERE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
E8_DIR=${E8_DIR:-$S1B_HERE/../e8}
# shellcheck source=../e8/e8lib.sh
. "$E8_DIR/e8lib.sh"

# the largest text that the reader (the Structure level and the engine) takes: SesaVerifier and EsaSuffixTree, ADR-004 1b.7f
S1B_READER_LIMIT=2147483645
S1B_GATE_MB=320

die() { echo "s1b: $*" >&2; exit 2; }
warnx() { echo "s1b: warning: $*" >&2; }

# sha256_of FILE -> hex (sha256sum on Linux, shasum on macOS)
sha256_of() {
  if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

# fsize FILE -> bytes
fsize() { stat -c %s "$1" 2> /dev/null || stat -f %z "$1"; }

# now_s -> seconds since the epoch (integer)
now_s() { date +%s; }

# bytes per symbol of the container by text kind (MEASURED at 10^7 and 10^8, see the ADR: the container with PSI; the sigma sets the bits of the packed text, the repeats set the overflow list)
container_bpc() {
  case $1 in
    dna) echo 13.25 ;;      # 2 bits per symbol: TEXT 0.25 + SA 4 + CHILD 4 + PSI 4 + LCP 1 (no overflow)
    genome) echo 13.50 ;;   # 3 bits (A C G N T) + 1.39 M overflow pairs at 10^8
    rep) echo 13.40 ;;
    prot) echo 13.70 ;;     # 5 bits for 20 letters
    *) echo 13.70 ;;
  esac
}
# scratch bytes per symbol at the peak (before the suffix array is complete); the build plans for 12 (SesaBuilder.ScratchBytesPerChar)
S1B_SCRATCH_BPC=13

# need_bytes N KIND -> the disk that one run needs at its peak: the text, the container (reserved from the start), the scratch at its peak, and 2 GB + 2 % of margin
need_bytes() { awk -v n="$1" -v c="$(container_bpc "$2")" -v s="$S1B_SCRATCH_BPC" 'BEGIN { printf "%.0f", n * (1 + c + s) * 1.02 + 2e9 }'; }

# label of a run
run_label() { echo "$1-n$2${3:+-$3}"; }

# is_reader_n N -> true when the reader can open a container of N symbols
is_reader_n() { [ "$1" -le "$S1B_READER_LIMIT" ]; }

# json_escape STRING
json_escape() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }

# phase of a Bench marker line "M<TAB>name<TAB>seconds" -> the last one in a file
last_marker() { grep -a '^M	' "$1" 2> /dev/null | tail -1 | cut -f2,3 | tr '\t' ' '; }

# portable free and total bytes of the file system of DIR (POSIX df -P, 1024-byte blocks; works on macOS and Linux)
free_bytes_of() { if [ -n "${S1B_FREE_OVERRIDE:-}" ]; then echo "$S1B_FREE_OVERRIDE"; return; fi; df -Pk "$1" 2> /dev/null | awk 'NR == 2 { printf "%.0f", $4 * 1024 }'; }   # S1B_FREE_OVERRIDE: a test hook for the disk checks
size_bytes_of() { df -Pk "$1" 2> /dev/null | awk 'NR == 2 { printf "%.0f", $2 * 1024 }'; }
