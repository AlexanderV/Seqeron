// Multiple sequence alignment benchmark (audit A6-8d, B09): the exact center-star alignment SequenceAligner.MultipleAlign against the anchor-based
// SequenceAligner.MultipleAlignAnchored on related and random DNA sets, printed as a table. It used to be an [Explicit] NUnit fixture that no test tier ran
// (wall-clock numbers are unreliable under the parallel suite); the correctness properties it checked are the deterministic tests of
// tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignAnchored_Scenarios_Tests.cs.
//
// Run on an idle machine:  dotnet run -c Release --project apps/Seqeron.Alignment.Benchmarks

using System.Diagnostics;

var scoring = SequenceAligner.SimpleDna;
var scenarios = new (string Label, List<DnaSequence> Sequences)[]
{
    ("5x100,  95%", Related(5, 100, 0.95, 10)),
    ("5x200,  90%", Related(5, 200, 0.90, 20)),
    ("10x500, 90%", Related(10, 500, 0.90, 30)),
    ("10x500, 80%", Related(10, 500, 0.80, 40)),
    ("10x1K,  90%", Related(10, 1000, 0.90, 50)),
    ("20x1K,  90%", Related(20, 1000, 0.90, 60)),
    ("10x500, rnd", Random(10, 500, 70)),
};

Console.WriteLine("| Scenario     | Classic ms | Anchor ms | Speedup |");
Console.WriteLine("|--------------|-----------:|----------:|--------:|");
foreach (var (label, seqs) in scenarios)
{
    SequenceAligner.MultipleAlign(seqs, scoring);                 // warm-up
    SequenceAligner.MultipleAlignAnchored(seqs, scoring);
    GC.Collect();
    GC.WaitForPendingFinalizers();
    GC.Collect();
    var sw = Stopwatch.StartNew();
    SequenceAligner.MultipleAlign(seqs, scoring);
    double classicMs = sw.Elapsed.TotalMilliseconds;
    GC.Collect();
    sw.Restart();
    SequenceAligner.MultipleAlignAnchored(seqs, scoring);
    double anchorMs = sw.Elapsed.TotalMilliseconds;
    Console.WriteLine($"| {label,-12} | {classicMs,10:F2} | {anchorMs,9:F2} | {classicMs / Math.Max(anchorMs, 0.001),6:F2}x |");
}

Console.WriteLine();
Console.WriteLine($"Runtime: {System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription}; OS: {System.Runtime.InteropServices.RuntimeInformation.OSDescription}");

static string Dna(int length, int seed)
{
    var rng = new Random(seed);
    return new string(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
}

// a reference and (count - 1) copies with point mutations at the rate 1 - similarity and about 1% indels
static List<DnaSequence> Related(int count, int length, double similarity, int seed)
{
    var rng = new Random(seed);
    string reference = Dna(length, seed);
    var sequences = new List<DnaSequence> { new(reference) };
    for (int s = 1; s < count; s++)
    {
        var chars = reference.ToCharArray();
        for (int m = 0; m < (int)(length * (1.0 - similarity)); m++)
        {
            int pos = rng.Next(length);
            char replacement;
            do { replacement = "ACGT"[rng.Next(4)]; } while (replacement == chars[pos]);
            chars[pos] = replacement;
        }

        var sb = new System.Text.StringBuilder(new string(chars));
        for (int d = 0; d < Math.Max(1, length / 100) && sb.Length > 10; d++)
        {
            int pos = rng.Next(sb.Length);
            if (rng.NextDouble() < 0.5)
                sb.Insert(pos, "ACGT"[rng.Next(4)]);
            else
                sb.Remove(pos, 1);
        }

        sequences.Add(new DnaSequence(sb.ToString()));
    }

    return sequences;
}

static List<DnaSequence> Random(int count, int length, int seed)
    => [.. Enumerable.Range(0, count).Select(i => new DnaSequence(Dna(length, seed + (i * 1000))))];
