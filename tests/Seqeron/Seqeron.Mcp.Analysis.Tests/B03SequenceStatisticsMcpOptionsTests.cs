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
}
