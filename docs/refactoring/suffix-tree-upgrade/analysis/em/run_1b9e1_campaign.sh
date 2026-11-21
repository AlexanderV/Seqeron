#!/bin/bash
# Work package 1b.9e.1: the kill campaign of 1b.9e (run_1b9e_kill.py: `kill -9` of a separate process at random moments inside every phase of a build with Resume, optionally a second kill inside the resume, then a resume in a NEW
# process with 1, 2 or 4 threads, SHA-256 against the uninterrupted build, nothing but the target left) run again with the code of this package: the resumes now normalise the container before they go on.
# usage: run_1b9e1_campaign.sh <bench dll> <work directory> <rows directory> <data directory with dna.txt rep.txt prot.txt (10^7) and dna100.txt rep100.txt (10^8)> [per-phase at 10^7] [per-phase at 10^8]
# One heavy job at a time. The rows are rows_1b9e1_kill_<kind>.tsv; summarize_1b9e_kill.py turns them into the tables.
set -u
BENCH=$1; WORK=$2; ROWS=$3; DATA=$4; P7=${5:-2}; P8=${6:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK" "$ROWS"
for kind in dna rep prot; do
  python3 "$HERE/run_1b9e_kill.py" --bench "$BENCH" --input "$DATA/$kind.txt" --kind "${kind}_10m" --work "$WORK/$kind" --rows "$ROWS/rows_1b9e1_kill_${kind}_10m.tsv" --per-phase "$P7" --boundaries no --seed 21 --chain 0.25
  rm -rf "$WORK/$kind"
done
for kind in dna100 rep100; do
  python3 "$HERE/run_1b9e_kill.py" --bench "$BENCH" --input "$DATA/$kind.txt" --kind "${kind%100}_100m" --work "$WORK/$kind" --rows "$ROWS/rows_1b9e1_kill_${kind%100}_100m.tsv" --per-phase "$P8" --boundaries no --seed 22 --chain 0.25
  rm -rf "$WORK/$kind"
done
