using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Running;
using BenchmarkDotNet.Toolchains.InProcess.NoEmit;
using SuffixTree.Benchmarks;

// ============================================================
// Two-phase JIT vs NativeAOT benchmark strategy:
//
// Phase 1 — JIT baseline (fast, ~2 min):
//   dotnet run --project apps\SuffixTree.Benchmarks -c Release -f net9.0 -- \
//     --filter "*Build_Short*" "*Build_DNA*" "*Contains*" "*LRS*"
//
// Phase 2 — Publish as NativeAOT ONCE (~5 min):
//   dotnet publish apps\SuffixTree.Benchmarks -c Release -r win-x64 -f net9.0 \
//     /p:PublishAot=true /p:OptimizationPreference=Speed \
//     /p:IlcInstructionSet=native /p:IlcFoldIdenticalMethodBodies=true \
//     /p:StripSymbols=true /p:InvariantGlobalization=true
//
// Phase 3 — Run AOT binary with InProcess toolchain (fast, ~2 min):
//   .\apps\SuffixTree.Benchmarks\bin\Release\net9.0\win-x64\publish\SuffixTree.Benchmarks.exe \
//     --inprocess --filter "*Build_Short*" "*Build_DNA*" "*Contains*" "*LRS*"
//
// Compare JIT vs AOT results from BenchmarkDotNet.Artifacts/results/
// ============================================================

var useInProcess = args.Contains("--inprocess");
var filteredArgs = args.Where(a => a != "--inprocess").ToArray();

if (useInProcess)
{
    // InProcess mode: benchmark runs inside the current process.
    // Use this when running from a pre-published NativeAOT binary —
    // no child process spawning, no re-compilation, instant start.
    var config = ManualConfig.Create(DefaultConfig.Instance)
        .AddJob(Job.ShortRun
            .WithToolchain(InProcessNoEmitToolchain.Instance)
            .WithId("NativeAOT-InProcess"));
    BenchmarkRunner.Run<SuffixTreeBenchmarks>(config, filteredArgs);
}
else
{
    BenchmarkRunner.Run<SuffixTreeBenchmarks>(args: filteredArgs);
}
