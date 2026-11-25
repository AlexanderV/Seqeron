using System.Text;
using Seqeron.Genomics.Chromosome;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for ASSEMBLY-COMPLETE-001. (1) EstimateCompletenessFromKmers against a naive in-test
/// reference of the Merqury / GenomeScope definitions (meryl-style canonical counting with sets, the
/// kmerHistToPloidyDepth first-valley rule, Merqury completeness / QV, GenomeScope error rate and size).
/// (2) AssessCompleteness on planted exact marker copies, whose BUSCO class is known by construction.
/// The Python oracle (tools/review/oracles/B10/completeness_oracle.py; Biopython local alignment,
/// scipy root of GenomeScope's f1) agreed 400/400 on both (Evidence).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyCompletenessFuzzTests
{
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) sb.Append("ACGT"[rng.Next(4)]);
        return sb.ToString();
    }

    private static string Canonical(string kmer)
    {
        string rc = DnaSequence.GetReverseComplementString(kmer);
        return string.CompareOrdinal(rc, kmer) < 0 ? rc : kmer;
    }

    private static Dictionary<string, long> NaiveCanonicalCounts(IEnumerable<string> seqs, int k)
    {
        var counts = new Dictionary<string, long>();
        foreach (string s in seqs)
            for (int i = 0; i + k <= s.Length; i++)
            {
                string w = s.Substring(i, k);
                if (w.Any(c => "ACGT".IndexOf(c) < 0)) continue;
                w = Canonical(w);
                counts[w] = counts.GetValueOrDefault(w) + 1;
            }
        return counts;
    }

    [Test]
    public void EstimateCompletenessFromKmers_RandomReadSets_MatchNaiveMerquryGenomeScopeReference()
    {
        var rng = new Random(4101);
        for (int trial = 0; trial < 400; trial++)
        {
            int k = new[] { 3, 4, 5, 7, 11 }[rng.Next(5)];
            string genome = RandomDna(rng, rng.Next(30, 300));
            int coverage = rng.Next(4, 25);
            var reads = new List<string>();
            for (int r = 0; r < coverage * genome.Length / 50 + 1; r++)
            {
                int a = rng.Next(0, Math.Max(1, genome.Length - 50));
                var read = genome.Substring(a, Math.Min(50, genome.Length - a)).ToCharArray();
                if (rng.Next(2) == 0) read[rng.Next(read.Length)] = "ACGT"[rng.Next(4)];
                string s = new(read);
                reads.Add(rng.Next(2) == 0 ? s : DnaSequence.GetReverseComplementString(s));
            }
            var spectrum = NaiveCanonicalCounts(reads, k);
            var assembly = new[] { ("a", genome[..(genome.Length / 2)]), ("b", "N" + genome[(genome.Length / 2)..] + RandomDna(rng, 5)) };
            int expected = rng.Next(4) == 0 ? coverage : 0;

            // Reference.
            var hist = spectrum.Values.GroupBy(c => c).OrderBy(g => g.Key).Select(g => (Depth: (int)g.Key, Count: (double)g.Count())).ToList();
            int filt = 0;
            for (int i = 1; i < hist.Count; i++)
            {
                double prevD = i == 1 ? hist[0].Count - double.MaxValue : hist[i - 1].Count - hist[i - 2].Count;
                double d = hist[i].Count - hist[i - 1].Count;
                if (d * prevD < 0 && d > prevD) { filt = hist[i].Depth - 1; break; }
            }
            long total = spectrum.Values.Sum(), error = spectrum.Values.Where(c => c <= filt).Sum();
            var above = hist.Where(h => h.Depth > filt).ToList();
            int peak = expected > 0 ? expected : above.First(h => h.Count == above.Max(x => x.Count)).Depth;
            var solid = spectrum.Where(p => p.Value > filt).Select(p => p.Key).ToHashSet();
            var asm = NaiveCanonicalCounts(assembly.Select(x => x.Item2), k);
            long asmTotal = asm.Values.Sum(), asmOnly = asm.Where(p => !spectrum.ContainsKey(p.Key)).Sum(p => p.Value);
            double asmErr = 1 - Math.Pow(1 - asmOnly / (double)asmTotal, 1.0 / k);

            var e = GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers(
                spectrum.Select(p => (p.Key, (int)p.Value)), assembly, expected);
            string ctx = $"trial {trial} k={k}";
            Assert.Multiple(() =>
            {
                Assert.That((e.K, e.SolidThreshold, e.PeakCoverage, e.TotalKmers, e.ErrorKmers), Is.EqualTo((k, filt, peak, total, error)), ctx);
                Assert.That(e.ErrorRate, Is.EqualTo(1 - Math.Pow(1 - error / (double)total, 1.0 / k)).Within(1e-12), ctx);
                Assert.That(e.EstimatedGenomeSize, Is.EqualTo((long)Math.Round((total - error) / (double)peak, MidpointRounding.ToEven)), ctx);
                Assert.That(e.SolidKmers, Is.EqualTo(solid.Count), ctx);
                Assert.That(e.AssemblySolidKmers, Is.EqualTo(asm.Keys.Count(solid.Contains)), ctx);
                Assert.That(e.Completeness, Is.EqualTo(asm.Keys.Count(solid.Contains) / (double)solid.Count).Within(1e-12), ctx);
                Assert.That((e.AssemblyKmers, e.AssemblyOnlyKmers), Is.EqualTo((asmTotal, asmOnly)), ctx);
                Assert.That(e.AssemblyErrorRate, Is.EqualTo(asmErr).Within(1e-12), ctx);
                if (asmOnly == 0) Assert.That(double.IsPositiveInfinity(e.QualityValue), Is.True, ctx);
                else Assert.That(e.QualityValue, Is.EqualTo(-10 * Math.Log10(asmErr)).Within(1e-9), ctx);
            });
        }
    }

    [Test]
    public void AssessCompleteness_PlantedExactCopies_HaveTheirConstructedBuscoClass()
    {
        var rng = new Random(4102);
        for (int trial = 0; trial < 150; trial++)
        {
            string marker = RandomDna(rng, rng.Next(40, 200));
            int m = marker.Length;
            // 0 = absent, 1 = one complete copy, 2 = two copies on different contigs,
            // 3 = 70 % fragment only, 4 = complete + 40 % fragment.
            int kind = rng.Next(5);
            string Copy() => rng.Next(2) == 0 ? marker : DnaSequence.GetReverseComplementString(marker);
            string Flank() => RandomDna(rng, rng.Next(0, 60));
            var contigs = new List<string> { Flank() + Flank() };
            switch (kind)
            {
                case 1: contigs.Add(Flank() + Copy() + Flank()); break;
                case 2: contigs.Add(Flank() + Copy() + Flank()); contigs.Add(Flank() + Copy() + Flank()); break;
                case 3: contigs.Add(Flank() + marker[..(int)Math.Ceiling(m * 0.7)] + Flank()); break;
                case 4: contigs.Add(Flank() + Copy() + Flank()); contigs.Add(Flank() + marker[(m - (int)(m * 0.4))..] + Flank()); break;
            }
            var asm = contigs.Select((s, i) => ($"c{i}", s)).ToList();
            var r = GenomeAssemblyAnalyzer.AssessCompleteness(asm, new[] { ("g", marker) });
            var expected = kind switch
            {
                0 => (0, 0, 0, 1),
                1 or 4 => (1, 0, 0, 0),
                2 => (0, 1, 0, 0),
                _ => (0, 0, 1, 0),
            };
            string ctx = $"trial {trial} kind {kind} marker {marker}";
            Assert.That((r.CompleteSingleCopy, r.CompleteDuplicated, r.Fragmented, r.Missing), Is.EqualTo(expected), ctx);

            foreach (var h in GenomeAssemblyAnalyzer.FindMarkerHits(marker, asm))
            {
                int length = asm.Single(a => a.Item1 == h.SequenceId).Item2.Length;
                Assert.That(h.Start >= 0 && h.Start <= h.End && h.End < length, Is.True, ctx);
                Assert.That(h.MarkerStart >= 0 && h.MarkerStart <= h.MarkerEnd && h.MarkerEnd < m, Is.True, ctx);
                Assert.That(h.Coverage, Is.EqualTo((h.MarkerEnd - h.MarkerStart + 1) / (double)m).Within(1e-12), ctx);
                Assert.That(h.Identity, Is.InRange(0.0, 1.0), ctx);
            }
        }
    }
}
