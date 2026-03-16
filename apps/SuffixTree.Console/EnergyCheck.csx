using System;
using System.Linq;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

var validHairpin57 = "GCAUAGCUAGCUAGCUAGCUAGCUA" + "GAAAUUU" + "UAGCUAGCUAGCUAGCUAGCUAUGC";
Console.WriteLine($"ValidHairpin57 length: {validHairpin57.Length}");

var results = FindPreMiRnaHairpins(validHairpin57, minHairpinLength: 55).ToList();
Console.WriteLine($"Candidates found: {results.Count}");
foreach (var r in results)
{
    Console.WriteLine($"  Start={r.Start} End={r.End} SeqLen={r.Sequence.Length} Energy={r.FreeEnergy:F4}");
    Console.WriteLine($"  Mature='{r.MatureSequence}' (len={r.MatureSequence.Length})");
    Console.WriteLine($"  Star  ='{r.StarSequence}' (len={r.StarSequence.Length})");
    Console.WriteLine($"  Struct='{r.Structure}'");
    Console.WriteLine();
}

// Also check stem20 hairpin
var stem20 = "GCAUAGCUAGCUAGCUAGCU" + "AUGCAUGCAUGCAUG" + "AGCUAGCUAGCUAGCUAUGC";
Console.WriteLine($"stem20Hairpin55 length: {stem20.Length}");
var results2 = FindPreMiRnaHairpins(stem20, minHairpinLength: 55).ToList();
Console.WriteLine($"Candidates: {results2.Count}");
foreach (var r in results2)
{
    Console.WriteLine($"  Start={r.Start} End={r.End} SeqLen={r.Sequence.Length} Energy={r.FreeEnergy:F4}");
    Console.WriteLine($"  Struct='{r.Structure}'");
}

// Also check wobble hairpin
var wobble57 = "GCAUAGCUAGCUAGCUAGCUAGCUA" + "GAAAUUU" + "UAGCUAGUUAGCUAGUUAGUUAUGU";
Console.WriteLine($"\nWobbleHairpin57 length: {wobble57.Length}");
var results3 = FindPreMiRnaHairpins(wobble57, minHairpinLength: 55).ToList();
Console.WriteLine($"Candidates: {results3.Count}");
foreach (var r in results3)
{
    Console.WriteLine($"  Start={r.Start} End={r.End} SeqLen={r.Sequence.Length} Energy={r.FreeEnergy:F4}");
    Console.WriteLine($"  Struct='{r.Structure}'");
}
