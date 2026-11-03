#!/usr/bin/env bash
B="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $B
for x in dna_1m prot_1m dna_10m rep_10m dna_20m dna_50m; do ./run_main.sh go $x 3; done
for x in dna_1m prot_1m dna_10m rep_10m; do ./run_main.sh seqmem $x 3; done
for x in dna_1m prot_1m dna_10m rep_10m; do ./run_main.sh seqpers $x 3; rm -f trees/$x.st; done
echo MAIN_DONE
