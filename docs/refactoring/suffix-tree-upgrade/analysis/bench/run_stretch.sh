#!/usr/bin/env bash
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $B
TMO=900 ./run_main.sh seqpers dna_20m 3; rm -f trees/dna_20m.st
TMO=900 ./run_main.sh seqpers dna_50m 3
echo STRETCH_DONE
