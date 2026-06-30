// 08_DIFFERENTIAL_TESTING rows 140, 147, 217. Independent oracles: manual codon-usage triplet scan,
// a manual N50/L50 cumulative computation, and a manual majority-vote consensus.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class AnnotAssembly_Differential_Tests
{
    // ---- Row 217: ANNOT-CODONUSAGE-001 — GetCodonUsage vs manual triplet scan ----

    [Test]
    [Category("ANNOT-CODONUSAGE-001")]
    [TestCase("ATGAAATTTATG")]
    [TestCase("ATGNNNAAATTT")]   // NNN excluded
    [TestCase("ATGAA")]          // trailing partial ignored
    public void GetCodonUsage_MatchesManualScan(string seq)
    {
        var s = seq.ToUpperInvariant();
        var expected = new Dictionary<string, int>();
        for (int i = 0; i + 3 <= s.Length; i += 3)
        {
            string c = s.Substring(i, 3);
            if (c.All(ch => "ACGT".Contains(ch))) expected[c] = expected.GetValueOrDefault(c) + 1;
        }
        Assert.That(GenomeAnnotator.GetCodonUsage(seq), Is.EquivalentTo(expected));
    }

    // ---- Row 147: ASSEMBLY-STATS-001 — N50 / L50 vs manual cumulative ----

    private static (int n50, int l50) N50Oracle(int[] lengths)
    {
        var sorted = lengths.OrderByDescending(x => x).ToList();
        long total = sorted.Sum(x => (long)x);
        long cum = 0; int count = 0;
        foreach (int len in sorted)
        {
            cum += len; count++;
            if (cum * 100 >= total * 50) return (len, count);
        }
        return (sorted.Count > 0 ? sorted[^1] : 0, count);
    }

    [Test]
    [Category("ASSEMBLY-STATS-001")]
    [TestCase(new[] { 100, 200, 300, 400, 500 })]
    [TestCase(new[] { 50, 50, 50, 50 })]
    [TestCase(new[] { 1000, 1, 1, 1 })]
    public void N50_MatchesManualCumulative(int[] lengths)
    {
        var (n50, l50) = N50Oracle(lengths);
        Assert.That(GenomeAssemblyAnalyzer.CalculateN50(lengths), Is.EqualTo(n50));
        var nx = GenomeAssemblyAnalyzer.CalculateNx(lengths, 50);
        Assert.That(nx.Nx, Is.EqualTo(n50));
        Assert.That(nx.Lx, Is.EqualTo(l50));
    }

    // ---- Row 140: ASSEMBLY-CONSENSUS-001 — ComputeConsensus vs manual majority vote ----

    private static string ConsensusOracle(string[] reads, double threshold = 0.5)
    {
        int length = reads.Max(r => r.Length);
        var sb = new System.Text.StringBuilder();
        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            int atoms = 0;
            foreach (var r in reads)
            {
                if (pos >= r.Length) continue;
                char c = char.ToUpperInvariant(r[pos]);
                if (c == '-' || c == '.') continue;
                counts[c] = counts.GetValueOrDefault(c) + 1; atoms++;
            }
            int maxSize = 0, ties = 0; char best = 'N';
            foreach (var kv in counts)
            {
                if (kv.Value > maxSize) { maxSize = kv.Value; ties = 1; best = kv.Key; }
                else if (kv.Value == maxSize) ties++;
            }
            bool commit = ties == 1 && atoms > 0 && (double)maxSize / atoms >= threshold;
            sb.Append(commit ? best : 'N');
        }
        return sb.ToString();
    }

    [Test]
    [Category("ASSEMBLY-CONSENSUS-001")]
    public void ComputeConsensus_MatchesManualMajority()
    {
        foreach (var reads in new[]
        {
            new[] { "ACGT", "ACGT", "AGGT" },
            new[] { "AC", "GC" },          // col0 tie -> N
            new[] { "A-GT", "ACGT" },      // gap skipped
            new[] { "AAA", "AAT", "ATA" },
        })
        {
            Assert.That(SequenceAssembler.ComputeConsensus(reads), Is.EqualTo(ConsensusOracle(reads)),
                string.Join(",", reads));
        }
    }
}
