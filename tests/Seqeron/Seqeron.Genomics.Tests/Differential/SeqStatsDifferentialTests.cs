// 08_DIFFERENTIAL_TESTING rows 74, 121, 122, 123, 130. Independent oracles: miRNA seed substring,
// manual nucleotide composition, manual sliding dinucleotide frequencies, an independent Kyte-Doolittle
// hydropathy mean, and the closed-form Wallace/Marmur-Doty Tm.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class SeqStatsDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 74: MIRNA-SEED-001 — GetSeedSequence vs substring extraction ----

    [Test]
    [Category("MIRNA-SEED-001")]
    [TestCase("UGAGGUAGUAGGUUGUAUAGUU")]   // let-7
    [TestCase("ACGUACGU")]
    [TestCase("ACGUAC")]                     // length < 8 -> ""
    [TestCase("")]
    public void GetSeedSequence_MatchesSubstring(string mir)
    {
        string expected = mir.Length < 8 ? "" : mir.Substring(1, 7).ToUpperInvariant();
        Assert.That(MiRnaAnalyzer.GetSeedSequence(mir), Is.EqualTo(expected));
    }

    // ---- Row 121: SEQ-COMPOSITION-001 — nucleotide composition vs manual count ----

    [Test]
    [Category("SEQ-COMPOSITION-001")]
    [TestCase("ACGTACGTNNXX")]
    [TestCase("GGGGCCCC")]
    [TestCase("AUAUAU")]
    public void NucleotideComposition_MatchesManualCount(string seq)
    {
        var s = seq.ToUpperInvariant();
        int a = s.Count(c => c == 'A'), t = s.Count(c => c == 'T'), g = s.Count(c => c == 'G');
        int c2 = s.Count(c => c == 'C'), u = s.Count(c => c == 'U'), n = s.Count(c => c == 'N');
        int total = a + t + g + c2 + u;

        var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
        Assert.That((comp.CountA, comp.CountT, comp.CountG, comp.CountC, comp.CountU, comp.CountN),
            Is.EqualTo((a, t, g, c2, u, n)));
        Assert.That(comp.GcContent, Is.EqualTo(total > 0 ? (double)(g + c2) / total : 0).Within(Tol));
        Assert.That(comp.AtContent, Is.EqualTo(total > 0 ? (double)(a + t + u) / total : 0).Within(Tol));
        Assert.That(comp.GcSkew, Is.EqualTo((g + c2) > 0 ? (double)(g - c2) / (g + c2) : 0).Within(Tol));
        Assert.That(comp.AtSkew, Is.EqualTo((a + t) > 0 ? (double)(a - t) / (a + t) : 0).Within(Tol));
    }

    // ---- Row 122: SEQ-DINUC-001 — dinucleotide frequencies vs manual sliding count ----

    [Test]
    [Category("SEQ-DINUC-001")]
    [TestCase("ACGTACGT")]
    [TestCase("AAAA")]
    [TestCase("ACGTNNACGT")]
    public void DinucleotideFrequencies_MatchesManualSlidingCount(string seq)
    {
        var s = seq.ToUpperInvariant();
        var counts = new Dictionary<string, int>();
        int total = 0;
        for (int i = 0; i < s.Length - 1; i++)
        {
            string di = s.Substring(i, 2);
            if (di.All(c => "ATGCU".Contains(c))) { counts[di] = counts.GetValueOrDefault(di) + 1; total++; }
        }
        var expected = counts.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);

        var actual = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        Assert.That(actual.Count, Is.EqualTo(expected.Count));
        foreach (var kv in expected)
            Assert.That(actual[kv.Key], Is.EqualTo(kv.Value).Within(Tol), kv.Key);
    }

    // ---- Row 123: SEQ-HYDRO-001 — hydrophobicity vs independent Kyte-Doolittle mean ----

    private static readonly Dictionary<char, double> Kd = new()
    {
        ['A'] = 1.8, ['R'] = -4.5, ['N'] = -3.5, ['D'] = -3.5, ['C'] = 2.5, ['E'] = -3.5,
        ['Q'] = -3.5, ['G'] = -0.4, ['H'] = -3.2, ['I'] = 4.5, ['L'] = 3.8, ['K'] = -3.9,
        ['M'] = 1.9, ['F'] = 2.8, ['P'] = -1.6, ['S'] = -0.8, ['T'] = -0.7, ['W'] = -0.9,
        ['Y'] = -1.3, ['V'] = 4.2,
    };

    [Test]
    [Category("SEQ-HYDRO-001")]
    [TestCase("ACDEFGHIKLMNPQRSTVWY")]
    [TestCase("IIIVVVLLL")]
    [TestCase("MKLV")]
    public void Hydrophobicity_MatchesKyteDoolittleMean(string protein)
    {
        var residues = protein.ToUpperInvariant().Where(c => Kd.ContainsKey(c)).ToList();
        double expected = residues.Count > 0 ? residues.Sum(c => Kd[c]) / residues.Count : 0;
        Assert.That(SequenceStatistics.CalculateHydrophobicity(protein), Is.EqualTo(expected).Within(1e-12));
    }

    // ---- Row 130: SEQ-TM-001 — Tm vs closed-form Wallace / Marmur-Doty ----

    private static double TmOracle(string seq)
    {
        var s = seq.ToUpperInvariant();
        int a = s.Count(c => c == 'A'), t = s.Count(c => c == 'T'), g = s.Count(c => c == 'G'), c2 = s.Count(c => c == 'C');
        if (s.Length < 14) return 2.0 * (a + t) + 4.0 * (g + c2);   // Wallace
        int total = a + t + g + c2;
        return total == 0 ? 0 : 64.9 + 41.0 * (g + c2 - 16.4) / total; // Marmur-Doty
    }

    [Test]
    [Category("SEQ-TM-001")]
    [TestCase("ACGT")]
    [TestCase("GGGGCCCC")]
    [TestCase("ATGCATGCATGCATGCATGC")]
    [TestCase("GGGGCCCCGGGGCCCCGGGG")]
    public void MeltingTemperature_MatchesClosedForm(string seq)
    {
        Assert.That(SequenceStatistics.CalculateMeltingTemperature(seq), Is.EqualTo(TmOracle(seq)).Within(1e-9));
    }
}
