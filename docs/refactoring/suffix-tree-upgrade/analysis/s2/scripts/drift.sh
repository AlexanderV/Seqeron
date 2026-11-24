#!/usr/bin/env bash
# The machine drift control (analysis/bench/HANDOFF.md section 2, PROGRESS-S1.md "Measurement protocol"): the libsais single-thread SA build of dna_10m, 5 fresh processes.
# The S1 reference is 0.328 s; a session whose median differs by more than 10 % is annotated "drift". usage: drift.sh <label>   (appends one row to results/drift_s2.tsv)
set -u
. "$(dirname "${BASH_SOURCE[0]}")/common.sh"
label=${1:?label}
quiet
vals=()
for i in 1 2 3 4 5; do
  v=$("$LIBSAIS_BIN" "$D/data/dna_10m.txt" "$D/queries/dna_10m" | awk -F'\t' '$1=="R" && $2=="build" && $3=="time" {print $4; exit}')
  vals+=("$v")
done
printf '%s\n' "${vals[@]}" | sort -g | awk -v l="$label" -v ref=0.328 '{a[NR]=$1} END {m=a[int((NR+1)/2)]; r=m/ref; d=(r<0.9||r>1.1)?"drift":"ok"; printf "%s\t%s\tlibsais_dna_10m_build_s\t%.4f\t%.4f\t%.4f\t%.3f\t%s\n", strftime("%Y-%m-%dT%H:%M:%SZ", systime(), 1), l, m, a[1], a[NR], r, d}' | tee -a "$S2/results/drift_s2.tsv"
