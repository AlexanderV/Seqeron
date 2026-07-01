// 08_DIFFERENTIAL_TESTING rows 210, 212, 227, 230. Independent oracles: (A−T)/(A+T) AT-skew, an RNA
// complement base table, a non-overlapping triplet codon-frequency scan, and a manual Shannon k-mer
// entropy.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class CompositionExtraDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 210: SEQ-ATSKEW-001 — (A−T)/(A+T) ----

    [Test]
    [Category("SEQ-ATSKEW-001")]
    [TestCase("AAAATTT")]
    [TestCase("AAAA")]      // no T -> +1
    [TestCase("TTTT")]      // no A -> -1
    [TestCase("GGGCCC")]    // no A/T -> 0
    [TestCase("AATTGGCC")]
    public void AtSkew_MatchesFormula(string seq)
    {
        var s = seq.ToUpperInvariant();
        int a = s.Count(c => c == 'A'), t = s.Count(c => c == 'T');
        double expected = (a + t) > 0 ? (double)(a - t) / (a + t) : 0;
        Assert.That(GcSkewCalculator.CalculateAtSkew(seq), Is.EqualTo(expected).Within(Tol));
    }

    // ---- Row 212: SEQ-RNACOMP-001 — RNA complement base table ----

    private static readonly Dictionary<char, char> RnaComp = new() { ['A'] = 'U', ['U'] = 'A', ['C'] = 'G', ['G'] = 'C' };

    [Test]
    [Category("SEQ-RNACOMP-001")]
    [TestCase("ACGU")]
    [TestCase("AAAA")]
    [TestCase("GCGCGC")]
    [TestCase("AUGCUAGC")]
    public void RnaComplement_MatchesBaseTable(string rna)
    {
        string expected = new string(rna.Select(c => RnaComp[c]).ToArray());
        Assert.That(new RnaSequence(rna).Complement().Sequence, Is.EqualTo(expected));
    }

    // ---- Row 227: SEQ-CODON-FREQ-001 — non-overlapping triplet frequencies ----

    [Test]
    [Category("SEQ-CODON-FREQ-001")]
    [TestCase("ATGAAATTT", 0)]
    [TestCase("ATGAAANTTGGG", 0)]   // NTT excluded
    [TestCase("AATGAAATTT", 1)]     // frame 1
    public void CodonFrequencies_MatchesTripletGroupBy(string seq, int frame)
    {
        var upper = seq.ToUpperInvariant();
        var counts = new Dictionary<string, int>();
        int total = 0;
        for (int i = frame; i + 3 <= upper.Length; i += 3)
        {
            string c = upper.Substring(i, 3);
            if (c.All(ch => "ATGC".Contains(ch))) { counts[c] = counts.GetValueOrDefault(c) + 1; total++; }
        }
        var expected = counts.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);

        var actual = SequenceStatistics.CalculateCodonFrequencies(seq, frame);
        Assert.That(actual.Count, Is.EqualTo(expected.Count));
        foreach (var kv in expected) Assert.That(actual[kv.Key], Is.EqualTo(kv.Value).Within(Tol), kv.Key);
    }

    // ---- Row 230: SEQ-COMPLEX-KMER-001 — k-mer entropy vs manual Shannon ----

    [Test]
    [Category("SEQ-COMPLEX-KMER-001")]
    [TestCase("ACGTACGT", 2)]
    [TestCase("AAAAA", 2)]
    [TestCase("ACGTACGTACGT", 3)]
    public void KmerEntropy_MatchesManualShannon(string seq, int k)
    {
        var s = seq.ToUpperInvariant();
        var counts = new Dictionary<string, int>();
        int total = 0;
        for (int i = 0; i + k <= s.Length; i++) { var w = s.Substring(i, k); counts[w] = counts.GetValueOrDefault(w) + 1; total++; }
        double expected = total == 0 ? 0 : -counts.Values.Sum(c => { double p = (double)c / total; return p * Math.Log2(p); });

        Assert.That(KmerAnalyzer.CalculateKmerEntropy(seq, k), Is.EqualTo(expected).Within(1e-12));
    }
}
