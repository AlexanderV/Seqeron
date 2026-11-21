#!/bin/bash
# Work package 1b.9e.1: the build with Resume = false must not slow down. The build of the commit before the package (ff71de38: the Esa dll in <old directory>) against this commit, both with Resume off, alternating in one
# process (`sesa resume abhead`), with and without fsync; pair 0 is the warm-up. The rows are summarised by summarize_1b9e_cost.py (the median and the spread of the pair ratios).
# usage: run_1b9e1_cost.sh <bench dll> <dll directory of the build before the package> <input> <label> <work directory> <rows file> [pairs]
# One heavy job at a time.
set -eu
BENCH=$1; OLD=$2; INPUT=$3; LABEL=$4; WORK=$5; ROWS=$6; PAIRS=${7:-5}
mkdir -p "$WORK"
echo "# $(date -u +%FT%TZ) $LABEL pairs $PAIRS uptime $(uptime)" >> "$ROWS"
dotnet "$BENCH" sesa resume abhead --old "$OLD" --input "$INPUT" --work "$WORK" --label "$LABEL.head.fsync" --fsync true --pairs "$PAIRS" >> "$ROWS"
dotnet "$BENCH" sesa resume abhead --old "$OLD" --input "$INPUT" --work "$WORK" --label "$LABEL.head.nofsync" --fsync false --pairs "$PAIRS" >> "$ROWS"
echo "# $(date -u +%FT%TZ) done uptime $(uptime)" >> "$ROWS"
