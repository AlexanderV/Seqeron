using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// FIN-B03 WP14/WP15: options of the Analysis MCP wrappers over <see cref="SequenceStatistics"/>.
/// Expected numbers are the executed reference values recorded in
/// docs/Validation/review-2026-09/B03.md (F8, F13, F14, F19, F24, F34, F37, F38) — seqinr rho,
/// Biopython 1.88 gc_fraction / protein_scale, scipy entropy and the hand-traced Chou-Fasman cases —
/// not the wrapper output. Default parameters are checked bit-identical to the pre-existing library calls.
/// </summary>
[TestFixture]
public class B03SequenceStatisticsMcpOptionsTests
{
    private static void AssertMap(IReadOnlyDictionary<string, double> actual, IReadOnlyDictionary<string, double> expected, string label)
    {
        Assert.That(actual.Keys, Is.EquivalentTo(expected.Keys), label + " keys");
        foreach (var (k, v) in expected)
            Assert.That(actual[k], Is.EqualTo(v).Within(1e-12), $"{label} {k}");
    }

    private static void AssertValues(double[] actual, double[] expected, string label)
    {
        Assert.That(actual, Has.Length.EqualTo(expected.Length), label + " length");
        for (int i = 0; i < expected.Length; i++)
            Assert.That(actual[i], Is.EqualTo(expected[i]).Within(1e-10), $"{label} [{i}]");
    }

    private static string RandomSequence(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    #region dinucleotide_ratios (A25, F13/F37)

    // Karlin rho* = seqinr rho on seq + separator + reverse complement (F13, exact fractions).
    [Test]
    public void DinucleotideRatios_StrandSymmetric_MatchesKarlinRhoStar()
    {
        AssertMap(AnalysisTools.DinucleotideRatios("ATGCGCGT", strandSymmetric: true).Ratios,
            new Dictionary<string, double>
            {
                ["AC"] = 128.0 / 105, ["AT"] = 256.0 / 63, ["CA"] = 128.0 / 105, ["CG"] = 512.0 / 175,
                ["GC"] = 512.0 / 175, ["GT"] = 128.0 / 105, ["TG"] = 128.0 / 105,
            }, "rho* ATGCGCGT");

        AssertMap(AnalysisTools.DinucleotideRatios("GATCCTTAAAGGCGCATTTAGGCCCATG", strandSymmetric: true).Ratios,
            new Dictionary<string, double>
            {
                ["AA"] = 40.0 / 27, ["AG"] = 8.0 / 9, ["AT"] = 16.0 / 9, ["CA"] = 8.0 / 9,
                ["CC"] = 40.0 / 27, ["CG"] = 16.0 / 27, ["CT"] = 8.0 / 9, ["GA"] = 16.0 / 27,
                ["GC"] = 16.0 / 9, ["GG"] = 40.0 / 27, ["TA"] = 32.0 / 27, ["TC"] = 16.0 / 27,
                ["TG"] = 8.0 / 9, ["TT"] = 40.0 / 27,
            }, "rho* 28-mer");

        Assert.That(AnalysisTools.DinucleotideRatios("ATGCGCGTNNacgtRAT", strandSymmetric: true).Ratios["CG"],
            Is.EqualTo(48.0 / 11).Within(1e-12), "rho* CG with N/R/lowercase");
    }

    // F37: Karlin classes (rho* <= 0.78 under, >= 1.23 over) — independent Python rho* values.
    [Test]
    public void DinucleotideRatios_IncludeRepresentation_StrandSymmetric_MatchesKarlinClasses()
    {
        var r = AnalysisTools.DinucleotideRatios("ATGCATGCAATTGGCCATGCAAACGTTTAAGG", strandSymmetric: true, includeRepresentation: true);
        Assert.That(r.Representation, Is.Not.Null);
        var over = new[] { "AA", "AT", "CA", "GC", "TG", "TT" };
        var under = new[] { "AC", "AG", "CG", "CT", "GT", "TA" };
        var normal = new[] { "CC", "GG" };
        Assert.Multiple(() =>
        {
            Assert.That(r.Representation!.Keys, Is.EquivalentTo(over.Concat(under).Concat(normal)));
            foreach (var d in over) Assert.That(r.Representation[d], Is.EqualTo("OverRepresented"), d);
            foreach (var d in under) Assert.That(r.Representation[d], Is.EqualTo("UnderRepresented"), d);
            foreach (var d in normal) Assert.That(r.Representation[d], Is.EqualTo("Normal"), d);
            Assert.That(r.Ratios["GC"], Is.EqualTo(2.69651086240948).Within(1e-12));
            Assert.That(r.Ratios["CG"], Is.EqualTo(0.67412771560237).Within(1e-12));
            Assert.That(r.Ratios["CC"], Is.EqualTo(1.011191573403555).Within(1e-12));
        });
    }

    [Test]
    public void DinucleotideRatios_IncludeRepresentation_SingleStrand_MatchesKarlinClasses()
    {
        var r = AnalysisTools.DinucleotideRatios("AAAATTTTCCCCGGGGACGT", includeRepresentation: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.Ratios["AA"], Is.EqualTo(2.526315789473684).Within(1e-12));
            Assert.That(r.Ratios["CG"], Is.EqualTo(1.6842105263157894).Within(1e-12));
            Assert.That(r.Ratios["AC"], Is.EqualTo(0.8421052631578947).Within(1e-12));
            Assert.That(r.Representation!["CG"], Is.EqualTo("OverRepresented"));
            foreach (var d in new[] { "AC", "AT", "GA", "GT", "TC" })
                Assert.That(r.Representation[d], Is.EqualTo("Normal"), d);
            Assert.That(r.Representation.Keys, Is.EquivalentTo(r.Ratios.Keys));
        });
    }

    [Test]
    public void DinucleotideRatios_Defaults_EqualLibrarySingleStrand_AndOmitRepresentation()
    {
        var rng = new Random(20261014);
        for (int n = 0; n < 300; n++)
        {
            string seq = RandomSequence(rng, "ACGTUNacgtR-", rng.Next(1, 60));
            var r = AnalysisTools.DinucleotideRatios(seq);
            var lib = SequenceStatistics.CalculateDinucleotideRatios(seq);
            Assert.That(r.Representation, Is.Null);
            Assert.That(r.Ratios.Keys, Is.EquivalentTo(lib.Keys), seq);
            foreach (var (k, v) in lib)
                Assert.That(r.Ratios[k], Is.EqualTo(v), seq + " " + k); // bit-identical

            var sym = AnalysisTools.DinucleotideRatios(seq, strandSymmetric: true).Ratios;
            var libSym = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: true);
            Assert.That(sym.Keys, Is.EquivalentTo(libSym.Keys), seq);
            foreach (var (k, v) in libSym)
                Assert.That(sym[k], Is.EqualTo(v), seq + " rho* " + k);
        }
    }

    [Test]
    public void DinucleotideRatios_RepresentationOmittedFromJsonByDefault()
    {
        string plain = System.Text.Json.JsonSerializer.Serialize(AnalysisTools.DinucleotideRatios("ATAT"));
        string withClasses = System.Text.Json.JsonSerializer.Serialize(
            AnalysisTools.DinucleotideRatios("ATAT", includeRepresentation: true));
        Assert.Multiple(() =>
        {
            Assert.That(plain, Does.Not.Contain("Representation"));
            Assert.That(withClasses, Does.Contain("\"Representation\""));
        });
    }

    #endregion

    #region chou_fasman_assign (A26, F14/F34)

    // F14 hand-traced Standard cases.
    [TestCase("EEEEEEVVVVVV", "HHHHEEEEEEEE")]
    [TestCase("PPPP", "CCCC")]
    [TestCase("NPDG", "TTTT")]
    [TestCase("SGIAKQ", "CCCCCC")]
    [TestCase("PTIQGQ", "CCCCCC")]
    public void ChouFasmanAssign_Standard_MatchesHandTrace(string sequence, string expected)
    {
        Assert.That(AnalysisTools.ChouFasmanAssign(sequence).Assignment, Is.EqualTo(expected));
        Assert.That(AnalysisTools.ChouFasmanAssign(sequence.ToLowerInvariant(), "standard").Assignment, Is.EqualTo(expected));
    }

    // F34 hand-traced Graded1978 / Standard pairs.
    [TestCase("AAAHHS", "HHHHHH", "CCCCCC")]
    [TestCase("EEEENY", "CCCCCC", "HHHHHH")]
    [TestCase("EEESRR", "HHHHHH", "CCCCCC")]
    [TestCase("VVVSW", "CCCCC", "EEEEE")]
    [TestCase("EEEEEEP", "HHHHHHC", "HHHHHHH")]
    [TestCase("AEEPEEEEE", "CHHHHHHHH", "HHHHHHHHH")]
    [TestCase("EEEEGG", "CCTTTT", "HHTTTT")]
    public void ChouFasmanAssign_Graded1978_MatchesHandTrace(string sequence, string graded, string standard)
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.ChouFasmanAssign(sequence, "Graded1978").Assignment, Is.EqualTo(graded));
            Assert.That(AnalysisTools.ChouFasmanAssign(sequence, "graded1978").Assignment, Is.EqualTo(graded));
            Assert.That(AnalysisTools.ChouFasmanAssign(sequence, "Standard").Assignment, Is.EqualTo(standard));
        });
    }

    [Test]
    public void ChouFasmanAssign_EqualsLibrary_Random()
    {
        var rng = new Random(1978);
        for (int n = 0; n < 300; n++)
        {
            string seq = RandomSequence(rng, "ACDEFGHIKLMNPQRSTVWYXacdpg", rng.Next(1, 80));
            Assert.That(AnalysisTools.ChouFasmanAssign(seq).Assignment,
                Is.EqualTo(SequenceStatistics.PredictSecondaryStructureChouFasman(seq)), seq);
            Assert.That(AnalysisTools.ChouFasmanAssign(seq, "Graded1978").Assignment,
                Is.EqualTo(SequenceStatistics.PredictSecondaryStructureChouFasman(seq, SequenceStatistics.ChouFasmanRuleSet.Graded1978)), seq);
        }
    }

    [Test]
    public void ChouFasmanAssign_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AnalysisTools.ChouFasmanAssign(""));
            Assert.Throws<ArgumentException>(() => AnalysisTools.ChouFasmanAssign(null!));
            Assert.Throws<ArgumentException>(() => AnalysisTools.ChouFasmanAssign("EEEE", "1"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.ChouFasmanAssign("EEEE", "Graded"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.ChouFasmanAssign("EEEE", null!));
        });
    }

    #endregion
    #region gc_content_profile (A27, F19/F24)

    // Biopython 1.88 [gc_fraction(s[i:i+6], mode) for i in range(0, len(s) - 5, 3)] (re-executed 2026-10-14).
    [TestCase("remove", new[] { 0.6666666666666666, 0.4, 0.0, 0.6666666666666666, 0.5 })]
    [TestCase("ignore", new[] { 0.6666666666666666, 0.3333333333333333, 0.0, 0.3333333333333333, 0.5 })]
    [TestCase("Weighted", new[] { 0.6666666666666666, 0.4166666666666667, 0.3333333333333333, 0.5833333333333334, 0.5 })]
    public void GcContentProfile_AmbiguityMode_MatchesBiopythonGcFraction(string mode, double[] expected)
    {
        const string seq = "ACGTSSWWNNRYacgusw";
        AssertValues(AnalysisTools.GcContentProfile(seq, 6, 3, fraction: true, ambiguityMode: mode).Values, expected, mode);
        AssertValues(AnalysisTools.GcContentProfile(seq, 6, 3, ambiguityMode: mode).Values,
            expected.Select(v => v * 100).ToArray(), mode + " percent");
    }

    // Biopython [gc_fraction(s[i:i+w]) for i in range(0, len(s), step)] — trailing partial windows (F24).
    [Test]
    public void GcContentProfile_IncludePartialWindow_MatchesBiopythonRangeSlicing()
    {
        AssertValues(AnalysisTools.GcContentProfile("GGGCCCAAATTTGCGCA", 6, 4, fraction: true, includePartialWindow: true).Values,
            new[] { 1.0, 1.0 / 3, 1.0 / 3, 0.8, 0.0 }, "w6 s4");
        AssertValues(AnalysisTools.GcContentProfile("GGGCCCAAATTTGCGCA", 6, 4, fraction: true).Values,
            new[] { 1.0, 1.0 / 3, 1.0 / 3 }, "w6 s4 complete only");
        AssertValues(AnalysisTools.GcContentProfile("ACGTNSWACGN", 4, 4, fraction: true, ambiguityMode: "remove", includePartialWindow: true).Values,
            new[] { 0.5, 0.3333333333333333, 1.0 }, "remove");
        AssertValues(AnalysisTools.GcContentProfile("ACGTNSWACGN", 4, 4, fraction: true, ambiguityMode: "ignore", includePartialWindow: true).Values,
            new[] { 0.5, 0.25, 0.6666666666666666 }, "ignore");
        AssertValues(AnalysisTools.GcContentProfile("ACGTNSWACGN", 4, 4, fraction: true, ambiguityMode: "weighted", includePartialWindow: true).Values,
            new[] { 0.5, 0.375, 0.8333333333333334 }, "weighted");
        AssertValues(AnalysisTools.GcContentProfile("GCA", 10, 1, fraction: true, includePartialWindow: true).Values,
            new[] { 2.0 / 3, 0.5, 0.0 }, "window longer than sequence");
    }

    [Test]
    public void GcContentProfile_Defaults_EqualLibrary_Random()
    {
        var rng = new Random(19);
        for (int n = 0; n < 300; n++)
        {
            string seq = RandomSequence(rng, "ACGTUNSWRYacgt-", rng.Next(1, 80));
            int w = rng.Next(1, 20), st = rng.Next(1, 6);
            Assert.That(AnalysisTools.GcContentProfile(seq, w, st).Values,
                Is.EqualTo(SequenceStatistics.CalculateGcContentProfile(seq, w, st).ToArray()), seq);
            Assert.That(AnalysisTools.GcContentProfile(seq, w, st, fraction: true, includePartialWindow: true).Values,
                Is.EqualTo(SequenceStatistics.CalculateGcContentProfile(seq, w, st, true, includePartialWindow: true).ToArray()), seq);
            Assert.That(AnalysisTools.GcContentProfile(seq, w, st, ambiguityMode: "weighted", includePartialWindow: true).Values,
                Is.EqualTo(SequenceStatistics.CalculateGcContentProfile(seq, w, st, false,
                    Seqeron.Genomics.Core.SequenceExtensions.GcAmbiguityMode.Weighted, true).ToArray()), seq);
        }
    }

    [Test]
    public void GcContentProfile_InvalidAmbiguityMode_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() => AnalysisTools.GcContentProfile("ACGT", 2, 1, ambiguityMode: "1"));
        Assert.That(ex!.ParamName, Is.EqualTo("ambiguityMode"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GcContentProfile("ACGT", 2, 1, ambiguityMode: "strict"));
    }

    #endregion

    #region hydrophobicity_profile (A28, F8/F29/F30)

    private const string Ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    private const string S22 = "MKTAYIAKQRQISFVKSHFSRQ";

    // Biopython 1.88 ProteinAnalysis(seq).protein_scale(ProtParamData.kd, window, edge) (re-executed 2026-10-14).
    [Test]
    public void HydrophobicityProfile_EdgeWeight_MatchesBiopythonProteinScale()
    {
        var ub = AnalysisTools.HydrophobicityProfile(Ubiquitin, 5, edgeWeight: 0.4).Values;
        Assert.That(ub, Has.Length.EqualTo(72));
        AssertValues(ub.Take(6).ToArray(),
            new[] { 2.015625, 1.853125, 1.546875, 0.37187499999999996, 0.19687499999999997, 0.3437499999999999 }, "ubiquitin w5 e0.4");
        AssertValues(AnalysisTools.HydrophobicityProfile("FLIVAG", 5, edgeWeight: 0.4).Values,
            new[] { 3.7312499999999997, 3.1156249999999996 }, "FLIVAG w5 e0.4");
        AssertValues(AnalysisTools.HydrophobicityProfile("FLIVAG", 3, edgeWeight: 0.0).Values,
            new[] { 3.8, 4.5, 4.2, 1.8 }, "FLIVAG w3 e0");
    }

    [Test]
    public void HydrophobicityProfile_BiopythonConvention_MatchesProteinScale()
    {
        AssertValues(AnalysisTools.HydrophobicityProfile("ACDXEFGHIK", 5, convention: "Biopython").Values,
            new[] { -1.04, -0.34, -1.48, -0.22000000000000003, 0.039999999999999925, -0.040000000000000126 }, "unknown X");
        AssertValues(AnalysisTools.HydrophobicityProfile(S22, 4, convention: "biopython").Values.Take(4).ToArray(),
            new[] { -0.3199999999999999, -0.45999999999999996, 0.6, 2.2600000000000002 }, "even window");
        // Standard convention differs off-centre: AXCDE (X at offset 1) Biopython 0.16 vs standard -0.54 (F29).
        AssertValues(AnalysisTools.HydrophobicityProfile("AXCDE", 5, convention: "Biopython").Values, new[] { 0.16 }, "AXCDE bio");
        AssertValues(AnalysisTools.HydrophobicityProfile("AXCDE", 5).Values, new[] { -0.54 }, "AXCDE standard");
    }

    // Biopython protein_scale on the min-max normalised ProtParamData.kd table (ProtScale 0..1, F30).
    [Test]
    public void HydrophobicityProfile_NormalizeScale_MatchesBiopythonOnNormalisedKd()
    {
        AssertValues(AnalysisTools.HydrophobicityProfile(S22, 9, normalizeScale: true).Values,
            new[]
            {
                0.4592592592592592, 0.3802469135802469, 0.3851851851851852, 0.4493827160493827, 0.41728395061728396,
                0.4679012345679012, 0.46419753086419746, 0.39382716049382716, 0.43209876543209874, 0.4358024691358025,
                0.5259259259259259, 0.5592592592592592, 0.44814814814814824, 0.4148148148148148,
            }, "kd W9 normalised");
    }

    [Test]
    public void HydrophobicityProfile_Strict_RejectsUnknownSymbols()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("ACDX", 1, strict: true));
            Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("ACDXEFGHIK", 5, convention: "Biopython", strict: true));
            Assert.That(AnalysisTools.HydrophobicityProfile("FLIVAG", 5, strict: true).Values,
                Is.EqualTo(AnalysisTools.HydrophobicityProfile("FLIVAG", 5).Values));
        });
    }

    [Test]
    public void HydrophobicityProfile_Defaults_EqualLibrary_Random()
    {
        var rng = new Random(1982);
        for (int n = 0; n < 300; n++)
        {
            string seq = RandomSequence(rng, "ACDEFGHIKLMNPQRSTVWYXBZ*-acdk", rng.Next(1, 60));
            int w = rng.Next(1, 15);
            Assert.That(AnalysisTools.HydrophobicityProfile(seq, w).Values,
                Is.EqualTo(SequenceStatistics.CalculateHydrophobicityProfile(seq, w).ToArray()), seq);
            if (w % 2 == 1)
            {
                Assert.That(AnalysisTools.HydrophobicityProfile(seq, w, edgeWeight: 0.3).Values,
                    Is.EqualTo(SequenceStatistics.CalculateHydrophobicityProfile(seq, w, 0.3).ToArray()), seq);
            }
        }
    }

    [Test]
    public void HydrophobicityProfile_InvalidOptions_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 5, edgeWeight: 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 5, edgeWeight: -0.1));
            Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 4, edgeWeight: 0.4));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 1, convention: "Biopython"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 5, convention: "1"));
            Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("FLIVAG", 5, convention: "ProtScale"));
        });
    }

    #endregion

    #region entropy_profile (A28 tail, F38)

    // scipy.stats.entropy(counts, base=2) over s[i:i+W], i in range(0, len(s), step) (re-executed 2026-10-14).
    [Test]
    public void EntropyProfile_IncludePartialWindow_MatchesScipy()
    {
        AssertValues(AnalysisTools.EntropyProfile("ACGTACGTAAAACCCGGT", 8, 3, includePartialWindow: true).Values,
            new[] { 2.0, 1.75, 1.75, 1.561278124459133, 1.4591479170272446, 0.9182958340544894 }, "letters partial");
        AssertValues(AnalysisTools.EntropyProfile("AAAAAAAAACGT", 5, 4, includePartialWindow: true).Values,
            new[] { 0.0, 0.0, 2.0 }, "AAAAAAAAACGT");
    }

    [Test]
    public void EntropyProfile_Kmer_MatchesScipy()
    {
        AssertValues(AnalysisTools.EntropyProfile("ACGTACGTAAAACCCGGT", 8, 3, k: 2).Values,
            new[] { 1.950212064914747, 2.2359263506290326, 2.128085278891394, 2.2359263506290326 }, "k2");
        AssertValues(AnalysisTools.EntropyProfile("ACGTACGTAAAACCCGGT", 8, 3, k: 2, includePartialWindow: true).Values,
            new[] { 1.950212064914747, 2.2359263506290326, 2.128085278891394, 2.2359263506290326, 1.9219280948873625, 1.0 }, "k2 partial");
        AssertValues(AnalysisTools.EntropyProfile("AAAAAAAAACGT", 6, 2, k: 1, includePartialWindow: true).Values,
            new[] { 0.0, 0.0, 0.6500224216483541, 1.7924812503605778, 2.0, 1.0 }, "k1 partial");
        AssertValues(AnalysisTools.EntropyProfile("GATTACA", 5, 2, k: 2, includePartialWindow: true).Values,
            new[] { 2.0, 2.0, 1.0, 0.0 }, "GATTACA k2 partial");
    }

    [Test]
    public void EntropyProfile_Defaults_EqualLibrary_AndKGuard()
    {
        var rng = new Random(38);
        for (int n = 0; n < 300; n++)
        {
            string seq = RandomSequence(rng, "ACGTNacgt-", rng.Next(1, 60));
            int w = rng.Next(1, 15), st = rng.Next(1, 5);
            Assert.That(AnalysisTools.EntropyProfile(seq, w, st).Values,
                Is.EqualTo(SequenceStatistics.CalculateEntropyProfile(seq, w, st).ToArray()), seq);
            Assert.That(AnalysisTools.EntropyProfile(seq, w, st, k: 3, includePartialWindow: true).Values,
                Is.EqualTo(SequenceStatistics.CalculateEntropyProfile(seq, w, st, 3, true).ToArray()), seq);
        }

        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.EntropyProfile("ACGT", 2, 1, k: 0));
    }

    #endregion

    #region codon_frequencies normalization / genetic_code (round 2 B5, F36 → F49)

    // EMBOSS 6.6.0 cusp executed 2026-10-14 (F36; re-run for F49) on this 26-codon sequence, standard code.
    private const string CuspSequence =
        "ATGGCTGCCGCAGCGGCTAAAAAGAAACTGCTTCTCTTATTGCTGCGTCGCAGAAGGTGGTGTTGCTAATAGTGAATG";

    // codon, cusp Fraction, cusp Frequency (/1000)
    private static readonly (string Codon, double Fraction, double PerThousand)[] CuspRows =
    {
        ("GCA", 0.200, 38.462), ("GCC", 0.200, 38.462), ("GCG", 0.200, 38.462), ("GCT", 0.400, 76.923),
        ("TGC", 0.500, 38.462), ("TGT", 0.500, 38.462), ("AAA", 0.667, 76.923), ("AAG", 0.333, 38.462),
        ("CTC", 0.167, 38.462), ("CTG", 0.333, 76.923), ("CTT", 0.167, 38.462), ("TTA", 0.167, 38.462),
        ("TTG", 0.167, 38.462), ("ATG", 1.000, 76.923), ("AGA", 0.250, 38.462), ("AGG", 0.250, 38.462),
        ("CGC", 0.250, 38.462), ("CGT", 0.250, 38.462), ("TGG", 1.000, 38.462), ("TAA", 0.333, 38.462),
        ("TAG", 0.333, 38.462), ("TGA", 0.333, 38.462),
    };

    [Test]
    public void CodonFrequencies_PerAminoAcid_MatchesEmbossCuspFraction()
    {
        var f = AnalysisTools.CodonFrequencies(CuspSequence, normalization: "PerAminoAcid").Frequencies;
        var lowerRna = AnalysisTools.CodonFrequencies(CuspSequence.ToLowerInvariant().Replace('t', 'u'), normalization: "peraminoacid").Frequencies;
        Assert.That(f.Keys, Is.EquivalentTo(CuspRows.Select(r => r.Codon)));
        Assert.Multiple(() =>
        {
            foreach (var row in CuspRows)
            {
                Assert.That(Math.Round(f[row.Codon], 3), Is.EqualTo(row.Fraction), row.Codon);
                Assert.That(lowerRna[row.Codon], Is.EqualTo(f[row.Codon]), row.Codon + " lower-case RNA");
            }

            Assert.That(f["AAA"], Is.EqualTo(2.0 / 3).Within(1e-12));
            Assert.That(f["TGA"], Is.EqualTo(1.0 / 3).Within(1e-12), "stops form one '*' family");
        });
    }

    [Test]
    public void CodonFrequencies_PerThousand_MatchesEmbossCuspFrequency()
    {
        var f = AnalysisTools.CodonFrequencies(CuspSequence, 0, "PerThousand").Frequencies;
        Assert.Multiple(() =>
        {
            foreach (var row in CuspRows)
                Assert.That(Math.Round(f[row.Codon], 3), Is.EqualTo(row.PerThousand), row.Codon);
            Assert.That(f["GCT"], Is.EqualTo(2 * 1000.0 / 26).Within(1e-9));
        });
    }

    // Hand-derived (F36): NCBI table 2 moves AGA/AGG into the stop family and TGA to Trp.
    [Test]
    public void CodonFrequencies_PerAminoAcid_GeneticCode2_RegroupsFamilies()
    {
        var f = AnalysisTools.CodonFrequencies(CuspSequence, normalization: "PerAminoAcid", genetic_code: 2).Frequencies;
        Assert.Multiple(() =>
        {
            foreach (var stop in new[] { "TAA", "TAG", "AGA", "AGG" })
                Assert.That(f[stop], Is.EqualTo(0.25).Within(1e-12), stop);
            Assert.That(f["TGA"], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(f["TGG"], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(f["CGC"], Is.EqualTo(0.5).Within(1e-12));
        });
    }

    [Test]
    public void CodonFrequencies_Defaults_EqualLibrary_AndOptionsEqualLibraryOverload()
    {
        var rng = new Random(49);
        foreach (int table in new[] { 1, 2, 11 })
        {
            for (int n = 0; n < 100; n++)
            {
                string seq = RandomSequence(rng, "ACGTUacgtN-", rng.Next(1, 60));
                int frame = rng.Next(0, 3);
                var code = Seqeron.Genomics.Core.GeneticCode.GetByTableNumber(table);
                Assert.That(AnalysisTools.CodonFrequencies(seq, frame).Frequencies,
                    Is.EquivalentTo(SequenceStatistics.CalculateCodonFrequencies(seq, frame)), seq);
                Assert.That(AnalysisTools.CodonFrequencies(seq, frame, "Total", table).Frequencies,
                    Is.EquivalentTo(SequenceStatistics.CalculateCodonFrequencies(seq, frame)), seq);
                foreach (var mode in new[] { SequenceStatistics.CodonFrequencyNormalization.PerThousand, SequenceStatistics.CodonFrequencyNormalization.PerAminoAcid })
                    Assert.That(AnalysisTools.CodonFrequencies(seq, frame, mode.ToString(), table).Frequencies,
                        Is.EquivalentTo(SequenceStatistics.CalculateCodonFrequencies(seq, frame, mode, code)), $"{seq} {mode} {table}");
            }
        }
    }

    [TestCase("Fraction")]
    [TestCase("2")]
    [TestCase("")]
    public void CodonFrequencies_UnknownNormalization_Throws(string normalization) =>
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies(CuspSequence, 0, normalization));

    [TestCase(0)]
    [TestCase(7)]
    [TestCase(34)]
    public void CodonFrequencies_UnknownGeneticCode_Throws(int table)
    {
        var ex = Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies(CuspSequence, 0, "PerAminoAcid", table));
        Assert.That(ex!.Message, Does.Contain("Unknown genetic code table"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies(CuspSequence, 0, "Total", table));
    }

    #endregion
}
