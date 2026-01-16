using BenchmarkDotNet.Running;
using SuffixTree.Benchmarks;

// Run benchmarks in Release mode:
// dotnet run -c Release
//
// To run specific category:
// dotnet run -c Release -- --filter "*Build*"
// dotnet run -c Release -- --filter "*Contains*"

BenchmarkRunner.Run<SuffixTreeBenchmarks>();
