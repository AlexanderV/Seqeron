#!/bin/bash
# usage: run_race.sh <mem|file> <processes> <rounds> <racebench dir> <results file> [<container> <text file>]
# one process per run of the client; the exit code of each is recorded (0 clean, 1 wrong answer or another exception, 134 / 139 a crash of the process, 124 a hang that the 300 s limit ended).
set -u
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
mode=$1; procs=$2; rounds=$3; dir=$4; out=$5; container=${6:-}; text=${7:-}
: > "$out"
clean=0; wrong=0; crash=0; hang=0
for i in $(seq 1 "$procs"); do
  if [ "$mode" = mem ]; then timeout 300 dotnet "$dir/RaceBench.dll" mem "$rounds" "$i" > "$out.$i.txt" 2>&1; else timeout 300 dotnet "$dir/RaceBench.dll" file "$rounds" "$i" "$container" "$text" > "$out.$i.txt" 2>&1; fi
  rc=$?
  line=$(grep -a '^RESULT' "$out.$i.txt" | head -1)
  echo -e "$i\t$rc\t$line" >> "$out"
  case $rc in 0) clean=$((clean+1));; 1) wrong=$((wrong+1));; 124) hang=$((hang+1));; *) crash=$((crash+1));; esac
  [ "$rc" = 0 ] && rm -f "$out.$i.txt"
done
echo -e "SUMMARY\t$mode\tprocesses=$procs\trounds=$rounds\tclean=$clean\twrong_or_other=$wrong\tcrash=$crash\thang=$hang" | tee -a "$out"
