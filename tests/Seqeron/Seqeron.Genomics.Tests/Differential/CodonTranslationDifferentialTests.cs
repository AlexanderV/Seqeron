// 08_DIFFERENTIAL_TESTING rows 58-63 (Codon + Translation). An INDEPENDENT hard-coded standard genetic
// code (built here from the published NCBI table 1) drives the oracles: codon-table equality, manual
// triplet translation, synonymous-optimisation protein preservation, hand-derived CAI geometric mean,
// manual rare-codon filtering, and manual codon-usage counting.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class CodonTranslationDifferentialTests
{
    private const double Tol = 1e-12;

    // Independent standard genetic code (NCBI translation table 1), RNA codons -> 1-letter AA / '*'.
    private static readonly Dictionary<string, char> StdCode = BuildStandardCode();

    private static Dictionary<string, char> BuildStandardCode()
    {
        const string bases = "UCAG";
        // AA order indexed by [first][second][third] per the canonical table layout.
        const string aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        var d = new Dictionary<string, char>();
        int i = 0;
        foreach (char b1 in bases)
            foreach (char b2 in bases)
                foreach (char b3 in bases)
                    d[$"{b1}{b2}{b3}"] = aa[i++];
        return d;
    }

    private static string TranslateOracle(string seq)
    {
        var rna = seq.ToUpperInvariant().Replace('T', 'U');
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i + 3 <= rna.Length; i += 3)
            sb.Append(StdCode[rna.Substring(i, 3)]);
        return sb.ToString();
    }

    // ---- Row 62: TRANS-CODON-001 — GeneticCode.Standard vs independent table (all 64 codons) ----

    [Test]
    [Category("TRANS-CODON-001")]
    public void GeneticCode_Standard_MatchesIndependentTable()
    {
        foreach (var (codon, expected) in StdCode)
        {
            Assert.That(GeneticCode.Standard.Translate(codon), Is.EqualTo(expected), $"RNA {codon}");
            Assert.That(GeneticCode.Standard.Translate(codon.Replace('U', 'T')), Is.EqualTo(expected), $"DNA {codon}");
        }
    }

    // ---- Row 63: TRANS-PROT-001 — Translator.Translate vs manual triplet loop ----

    [Test]
    [Category("TRANS-PROT-001")]
    [TestCase("ATGGCCTTTTAA")]
    [TestCase("ATGAAATTTGGCGCATGA")]
    [TestCase("TTTCCCGGGAAACCCTAG")]
    public void Translate_MatchesManualTripletLoop(string seq)
    {
        Assert.That(Translator.Translate(seq).Sequence, Is.EqualTo(TranslateOracle(seq)));
    }

    // ---- Row 58: CODON-OPT-001 — optimized synonymous sequence translates to the same protein ----

    [Test]
    [Category("CODON-OPT-001")]
    [TestCase("ATGAAATTTGGCGCATGCCTGTAA")]
    [TestCase("ATGCTGCTGAAAGATTTCGGCTAA")]
    public void OptimizeSequence_PreservesProtein(string seq)
    {
        var result = CodonOptimizer.OptimizeSequence(seq, CodonOptimizer.EColiK12);
        // Independent: both the original and the optimized sequence must translate to the same protein.
        Assert.That(TranslateOracle(result.OptimizedSequence), Is.EqualTo(TranslateOracle(seq)));
    }

    // ---- Custom partial table for CAI / rare-codon tests ----
    // Phe {UUU,UUC} and Lys {AAA,AAG}; relative adaptiveness w = freq / max(synonymous freq).
    private static readonly CodonOptimizer.CodonUsageTable CustomTable = new(
        "Custom",
        new Dictionary<string, double> { ["UUU"] = 0.5, ["UUC"] = 1.0, ["AAA"] = 1.0, ["AAG"] = 0.25 },
        new Dictionary<string, string> { ["UUU"] = "F", ["UUC"] = "F", ["AAA"] = "K", ["AAG"] = "K" });

    // ---- Row 59: CODON-CAI-001 — CAI vs hand-derived geometric mean ----

    [Test]
    [Category("CODON-CAI-001")]
    public void Cai_MatchesGeometricMeanOfAdaptiveness()
    {
        // TTT AAA TTC AAG -> w = [0.5, 1.0, 1.0, 0.25]; CAI = (0.5*1*1*0.25)^(1/4) = 0.125^0.25.
        double expected = Math.Pow(0.5 * 1.0 * 1.0 * 0.25, 1.0 / 4.0);
        Assert.That(CodonOptimizer.CalculateCAI("TTTAAATTCAAG", CustomTable), Is.EqualTo(expected).Within(1e-12));
    }

    // ---- Row 60: CODON-RARE-001 — FindRareCodons vs manual frequency lookup ----

    [Test]
    [Category("CODON-RARE-001")]
    public void FindRareCodons_MatchesManualFrequencyFilter()
    {
        const string seq = "TTTAAATTCAAG"; // UUU(0.5) AAA(1.0) UUC(1.0) AAG(0.25)
        const double threshold = 0.6;
        var actual = CodonOptimizer.FindRareCodons(seq, CustomTable, threshold)
            .Select(r => (r.Position, r.Codon, r.Frequency)).ToList();

        var rna = seq.Replace("T", "U");
        var expected = new List<(int, string, double)>();
        for (int i = 0; i + 3 <= rna.Length; i += 3)
        {
            string c = rna.Substring(i, 3);
            double f = CustomTable.CodonFrequencies.GetValueOrDefault(c, 0);
            if (f < threshold) expected.Add((i, c, f));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 61: CODON-USAGE-001 — CalculateCodonUsage vs manual triplet scan ----

    [Test]
    [Category("CODON-USAGE-001")]
    [TestCase("TTTAAATTCAAGTTT")]
    [TestCase("ATGATGATGAAA")]
    public void CodonUsage_MatchesManualScan(string seq)
    {
        var actual = CodonOptimizer.CalculateCodonUsage(seq);

        var rna = seq.ToUpperInvariant().Replace('T', 'U');
        var expected = new Dictionary<string, int>();
        for (int i = 0; i + 3 <= rna.Length; i += 3)
        {
            string c = rna.Substring(i, 3);
            expected[c] = expected.GetValueOrDefault(c, 0) + 1;
        }
        Assert.That(actual, Is.EquivalentTo(expected));
    }
}
