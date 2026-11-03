using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Bench;

// Measurement harness of the ESA backend (ST-S1 work package 1a.8). See README.md.
//
//   dotnet SuffixTree.Esa.Bench.dll build   --input <file> [--name dna_10m] [--verify <sa_vectors.tsv>]
//   dotnet SuffixTree.Esa.Bench.dll steady  --input <file> [--builds 5] [--name dna_10m]
//   dotnet SuffixTree.Esa.Bench.dll profile --input <file> [--builds 3] [--name dna_10m]
//
// Every command prints `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes` lines (the protocol of lib/benchlib.py). One process = one measurement:
// the input is loaded first (excluded from the timers), the build runs from a cold JIT.
return Cli.Run(args);
