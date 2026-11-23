#!/bin/bash
export PATH="$PATH:/root/.dotnet" DOTNET_ROOT=/root/.dotnet DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
P=/tmp/stbench/b9h/allocprobe_out/allocprobe.dll; D=/home/user/Seqeron/src/SuffixTree/Algorithms/SuffixTree.Esa/bin/Release/net10.0/SuffixTree.Esa.dll
OUT=/tmp/stbench/b9h/results_1b9h_allocprobe.txt
{
echo "# work package 1b.9h, item 3: allocprobe (analysis/bench/allocprobe), Release, .NET $(dotnet --version), 4 vCPU; quiet machine unless a row says load 4 (four busy processes); the columns are those of the program"
echo "# part 1: GC.GetAllocatedBytesForCurrentThread over one window, another thread that makes collections (helper 0 none, 1 small garbage flat out, 2 GC.Collect(0)+1 ms, 4 large garbage flat out, 5 GC.Collect(2, background)+1 ms); window n nothing, s a long[11], l a long[200001]; region 0/1"
for w in n s l; do for h in 0 1 2 4 5; do for r in 0 1; do dotnet $P $h 6 $r $w | cut -c1-420; done; done; done
echo "# control: the concurrent collector off (DOTNET_gcConcurrent=0), helper 4"
for w in n l; do DOTNET_gcConcurrent=0 dotnet $P 4 6 0 $w | cut -c1-420; done
echo "# part 2: the test of 1b.3 over the real SuffixTree.Esa.dll, 200,000 symbols (the two bodies one at a time: old = the plain counter, new = three no-GC regions); helper 0 none, 5 large arrays of 100 KB per millisecond, 4 of 200 KB flat out"
for h in 0 5 4; do for w in old new; do dotnet $P real $D $h 30 $w | cut -c1-520; done; done
echo "# part 2 under load: four busy processes beside the process"
pids=""; for i in 1 2 3 4; do python3 -c "
while True: pass" & pids="$pids $!"; done
for h in 0 5 4; do for w in old new; do echo -n "load 4: "; dotnet $P real $D $h 30 $w | cut -c1-520; done; done
kill $pids
echo done
} > $OUT 2>&1
echo finished >> /tmp/stbench/b9h/allocprobe_all.done
