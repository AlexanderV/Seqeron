#!/bin/bash
# Cost of work package 1b.9e: (1) Resume off against on, alternating in one process (the cost of recording the state: the fsyncs, the CRC pass, the state file), with and without fsync;
# (2) the build of the commit before the package against this commit with Resume off, alternating in one process (the default must not slow down).
# usage: run_1b9e_cost.sh <bench dll> <dll directory of the build before the package> <input> <label> <work directory> <rows file> [pairs]
# Every row is `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes`; summarize_1b9e_cost.py turns them into tables. One heavy job at a time.
set -eu
BENCH=$1; OLD=$2; INPUT=$3; LABEL=$4; WORK=$5; ROWS=$6; PAIRS=${7:-5}
mkdir -p "$WORK"
echo "# $(date -u +%FT%TZ) $LABEL pairs $PAIRS uptime $(uptime)" >> "$ROWS"
dotnet "$BENCH" sesa resume ab --input "$INPUT" --work "$WORK" --label "$LABEL.fsync" --fsync true --pairs "$PAIRS" >> "$ROWS"
dotnet "$BENCH" sesa resume ab --input "$INPUT" --work "$WORK" --label "$LABEL.nofsync" --fsync false --pairs "$PAIRS" >> "$ROWS"
dotnet "$BENCH" sesa resume abhead --old "$OLD" --input "$INPUT" --work "$WORK" --label "$LABEL.head.fsync" --fsync true --pairs "$PAIRS" >> "$ROWS"
dotnet "$BENCH" sesa resume abhead --old "$OLD" --input "$INPUT" --work "$WORK" --label "$LABEL.head.nofsync" --fsync false --pairs "$PAIRS" >> "$ROWS"
echo "# $(date -u +%FT%TZ) done uptime $(uptime)" >> "$ROWS"
