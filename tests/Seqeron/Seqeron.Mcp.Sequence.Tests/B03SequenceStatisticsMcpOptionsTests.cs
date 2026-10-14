using System.Text.Json;
using ModelContextProtocol;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// FIN-B03 WP13 (B03.md F39, audit rows A22–A24): optional SequenceStatistics parameters exposed by the
/// Sequence MCP wrappers. Reference values are Biopython 1.88 / EMBOSS 6.6.0 numbers recorded in B03.md
/// F2, F3, F6, F10, F26, F27 and F28 (re-executed 2026-10-14). Defaults must stay bit-identical.
/// </summary>
[TestFixture]
public class B03SequenceStatisticsMcpOptionsTests
{
    // ---------------- A22: nucleotide_composition N ----------------

    [Test]
    public void NucleotideComposition_CountsN_SeparatelyFromOther()
    {
        var r = SequenceTools.NucleotideComposition("ATGCNN");
        Assert.Multiple(() =>
        {
            Assert.That(r.N, Is.EqualTo(2));
            Assert.That(r.Other, Is.EqualTo(0));
            Assert.That(r.A + r.T + r.G + r.C + r.U + r.N + r.Other, Is.EqualTo(6));
        });
    }

    [Test]
    public void NucleotideComposition_LowercaseNAndOtherSymbols()
    {
        var r = SequenceTools.NucleotideComposition("acgunRx");
        Assert.Multiple(() =>
        {
            Assert.That(r.N, Is.EqualTo(1));
            Assert.That(r.Other, Is.EqualTo(2));
            Assert.That(r.A + r.T + r.G + r.C + r.U + r.N + r.Other, Is.EqualTo(r.Length));
        });
    }

    // ---------------- A23: isoelectric_point scale ----------------

    // EMBOSS 6.6.0 iep (F3): FLPVLAGLTPSIVPKLVCLLTKKC 9.5678, "A" 5.5500.
    // Biopython 1.88 IsoelectricPoint (F6): INGAR 9.750021, PETER 4.532080, AKD 6.131529.
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", "Emboss", 9.57)]
    [TestCase("A", "emboss", 5.55)]
    [TestCase("INGAR", "Bjellqvist", 9.75)]
    [TestCase("PETER", "BJELLQVIST", 4.53)]
    [TestCase("AKD", "bjellqvist", 6.13)]
    public void IsoelectricPoint_Scale_MatchesReference(string sequence, string scale, double expected)
    {
        Assert.That(SequenceTools.IsoelectricPoint(sequence, scale).PI, Is.EqualTo(expected).Within(1e-9));
    }

    [TestCase("INGAR")]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC")]
    [TestCase("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS")]
    public void IsoelectricPoint_DefaultScale_EqualsLibraryOneArgument(string sequence)
    {
        Assert.That(SequenceTools.IsoelectricPoint(sequence).PI,
            Is.EqualTo(SequenceStatistics.CalculateIsoelectricPoint(sequence)));
    }

    [TestCase("1")]
    [TestCase("")]
    [TestCase("ExPASy")]
    public void IsoelectricPoint_UnknownScale_Throws(string scale)
    {
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("INGAR", scale));
    }

    // ---------------- A23: thermodynamics selfComplementary + salt options ----------------

    private const string Owczarzy28 = "CGTTCCAAAGATGTGGGCATGAGCTTAC";

    // Biopython 1.88 Tm_NN(DNA_NN3, Na=50, dnac1=dnac2=125 nM) — F27.
    [TestCase(0, 0, 0, 0, 0, 77.61921589084807)]
    [TestCase(1, 0, 0, 0, 0, 56.02211796282597)]
    [TestCase(5, 0, 0, 0, 0, 61.9240199979082)]
    [TestCase(6, 0, 0, 0, 0, 61.375386349394205)]
    [TestCase(1, 0.020, 0.010, 0.0015, 0.0008, 65.07007431795537)]
    [TestCase(5, 0.020, 0.010, 0.0015, 0.0008, 68.32519243488775)]
    [TestCase(6, 0.020, 0.010, 0.0015, 0.0008, 70.40066503149728)]
    [TestCase(7, 0.020, 0.010, 0.0015, 0.0008, 67.73409261682195)]
    public void Thermodynamics_SaltOptions_MatchBiopythonTmNn(int method, double k, double tris, double mg, double dntps, double expectedTm)
    {
        var unrounded = SequenceTools.Thermodynamics(Owczarzy28, 0.05, 2.5e-7, false, method, k, tris, mg, dntps, roundResults: false);
        var rounded = SequenceTools.Thermodynamics(Owczarzy28, 0.05, 2.5e-7, false, method, k, tris, mg, dntps);
        Assert.Multiple(() =>
        {
            Assert.That(unrounded.MeltingTemperature, Is.EqualTo(expectedTm).Within(1e-9));
            Assert.That(rounded.MeltingTemperature, Is.EqualTo(Math.Round(expectedTm, 1, MidpointRounding.AwayFromZero)).Within(1e-9));
        });
    }

    [Test]
    public void Thermodynamics_Method7MagnesiumOnly_NaZero_MatchesBiopython()
    {
        // F27: method 7, Mg 3 mM only (Na 0): 69.36483747600732.
        var r = SequenceTools.Thermodynamics(Owczarzy28, 0, 2.5e-7, saltMethod: 7, mgConcentration: 0.003, roundResults: false);
        Assert.That(r.MeltingTemperature, Is.EqualTo(69.36483747600732).Within(1e-9));
    }

    [Test]
    public void Thermodynamics_SelfComplementary_MatchesBiopythonSelfcomp()
    {
        // F10: Biopython Tm_NN("CGCGAATTCGCG", DNA_NN3, Na=50, dnac1=250, dnac2=0, selfcomp=True) = 46.96872419752236.
        var rounded = SequenceTools.Thermodynamics("CGCGAATTCGCG", selfComplementary: true);
        var unrounded = SequenceTools.Thermodynamics("CGCGAATTCGCG", selfComplementary: true, roundResults: false);
        Assert.Multiple(() =>
        {
            Assert.That(rounded.MeltingTemperature, Is.EqualTo(47.0).Within(1e-9));
            Assert.That(unrounded.MeltingTemperature, Is.EqualTo(46.96872419752236).Within(1e-9));
        });
    }

    [Test]
    public void Thermodynamics_SelfComplementaryOwczarzy2008Mg_MatchesBiopython()
    {
        // F27: selfcomp GCGCATGCGC method 7 Na 50 Mg 2 (dnac1 250 nM) = 51.99165540136477.
        var r = SequenceTools.Thermodynamics("GCGCATGCGC", 0.05, 2.5e-7, true, 7, mgConcentration: 0.002, roundResults: false);
        Assert.That(r.MeltingTemperature, Is.EqualTo(51.99165540136477).Within(1e-9));
    }

    [Test]
    public void Thermodynamics_Defaults_BitIdenticalToLibraryThreeArgument()
    {
        var rng = new Random(13);
        const string alphabet = "ACGTNacgtu";
        for (int i = 0; i < 300; i++)
        {
            int len = rng.Next(2, 40);
            var chars = new char[len];
            for (int j = 0; j < len; j++)
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            string seq = new(chars);
            var lib = SequenceStatistics.CalculateThermodynamics(seq, 0.05, 2.5e-7);
            var mcp = SequenceTools.Thermodynamics(seq);
            Assert.That((mcp.DeltaH, mcp.DeltaS, mcp.DeltaG, mcp.MeltingTemperature),
                Is.EqualTo((lib.DeltaH, lib.DeltaS, lib.DeltaG, lib.MeltingTemperature)), seq);
        }
    }

    [Test]
    public void Thermodynamics_InosinePairWithoutParameters_Throws()
    {
        // F26: Biopython Tm_NN("ACIIGT") -> ValueError "no thermodynamic data for neighbors 'II/II' available".
        Assert.Throws<ArgumentException>(() => SequenceTools.Thermodynamics("ACIIGT"));
    }

    [Test]
    public void Thermodynamics_InvalidSaltOptions_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceTools.Thermodynamics(Owczarzy28, saltMethod: 8));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceTools.Thermodynamics(Owczarzy28, saltMethod: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceTools.Thermodynamics(Owczarzy28, mgConcentration: -0.001));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceTools.Thermodynamics(Owczarzy28, 0));
        });
    }

    // ---------------- A24: molecular_weight_nucleotide / molecular_weight_protein options ----------------

    // Biopython 1.88 molecular_weight — F2 (AGC average) and F28 (monoisotopic).
    [TestCase("AGC", true, false, false, false, 949.6095)]
    [TestCase("AGC", true, true, false, false, 1890.2057)]
    [TestCase("AGC", true, false, true, false, 931.5942)]
    [TestCase("AGC", true, true, true, false, 1854.1751)]
    [TestCase("AGC", false, true, false, false, 1972.1755)]
    [TestCase("ATGCGTAC", true, true, true, true, 4940.810128)]
    [TestCase("AUGCGUAC", false, false, false, true, 2588.343647)]
    public void MolecularWeightNucleotide_Options_MatchBiopython(
        string sequence, bool isDna, bool doubleStranded, bool circular, bool monoisotopic, double expected)
    {
        var r = SequenceTools.MolecularWeightNucleotide(sequence, isDna, doubleStranded, circular, monoisotopic);
        Assert.That(r.MolecularWeight, Is.EqualTo(expected).Within(1e-6));
    }

    [TestCase("ACGN", true)]
    [TestCase("ACGU", true)]
    [TestCase("ACGT", false)]
    public void MolecularWeightNucleotide_Strict_RejectsLikeBiopython(string sequence, bool isDna)
    {
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightNucleotide(sequence, isDna, strict: true));
        Assert.DoesNotThrow(() => SequenceTools.MolecularWeightNucleotide(sequence, isDna));
    }

    // Biopython 1.88 molecular_weight("GAGA", "protein", circular, monoisotopic) — F28.
    [TestCase(false, false, 274.2737)]
    [TestCase(false, true, 256.2584)]
    [TestCase(true, false, 274.127717)]
    [TestCase(true, true, 256.117152)]
    public void MolecularWeightProtein_Options_MatchBiopython(bool monoisotopic, bool circular, double expected)
    {
        Assert.That(SequenceTools.MolecularWeightProtein("GAGA", monoisotopic, circular).MolecularWeight,
            Is.EqualTo(expected).Within(1e-6));
    }

    [TestCase("MAXK")]
    [TestCase("MABK")]
    [TestCase("MAZK")]
    [TestCase("M*K")]
    public void MolecularWeightProtein_Strict_RejectsAmbiguousLikeBiopython(string sequence)
    {
        Assert.Throws<ArgumentException>(() => SequenceTools.MolecularWeightProtein(sequence, strict: true));
        Assert.DoesNotThrow(() => SequenceTools.MolecularWeightProtein(sequence));
    }

    [Test]
    public void MolecularWeight_Defaults_BitIdenticalToLibraryOverloads()
    {
        var rng = new Random(24);
        const string nt = "ACGTUNacgtu-";
        const string aa = "ACDEFGHIKLMNPQRSTVWYUOXBZ*acdk";
        for (int i = 0; i < 300; i++)
        {
            string n = RandomString(rng, nt, rng.Next(1, 30));
            string p = RandomString(rng, aa, rng.Next(1, 30));
            bool isDna = (i & 1) == 0;
            Assert.That(SequenceTools.MolecularWeightNucleotide(n, isDna).MolecularWeight,
                Is.EqualTo(SequenceStatistics.CalculateNucleotideMolecularWeight(n, isDna)), n);
            Assert.That(SequenceTools.MolecularWeightProtein(p).MolecularWeight,
                Is.EqualTo(SequenceStatistics.CalculateMolecularWeight(p)), p);
        }
    }

    // ---------------- Round 2 B3 (F47): isoelectric_point window / termini / modifications ----------------
    // EMBOSS 6.6.0 iep and Biopython 1.88 numbers re-executed 2026-10-14 (F43/F44).

    private const string Ifna2Human =
        "MALTFALLVALLVLSCKSSCSVGCDLPQTHSLGSRRTLMLLAQMRKISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTK" +
        "DSSAAWDETLLDKFYTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE";

    [Test]
    public void IsoelectricPoint_EmbossWindow_Ifna2Modified_MatchesIepDocumentation()
    {
        // iep -disulphides 2 -lysinemodified 2 → 5.7240 (iep documentation example).
        var r = SequenceTools.IsoelectricPoint(Ifna2Human, window: "Emboss", lysinesModified: 2, disulphides: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.PI, Is.EqualTo(5.7240).Within(5e-5));
            Assert.That(r.Found, Is.True);
        });
    }

    [Test]
    public void IsoelectricPoint_EmbossWindow_Poly40R_IsIepNone()
    {
        var r = SequenceTools.IsoelectricPoint(new string('R', 40), window: "emboss");
        Assert.Multiple(() =>
        {
            Assert.That(r.PI, Is.Null);
            Assert.That(r.Found, Is.False);
        });
        string json = JsonSerializer.Serialize(r, McpJsonUtilities.DefaultOptions);
        Assert.That(json, Is.EqualTo("{\"found\":false}"));
    }

    [TestCase("DDDD", "Bjellqvist", 4.0500284194946286)]
    [TestCase(Ifna2Human, "Bjellqvist", 6.32307643890381)]
    public void IsoelectricPoint_BiopythonWindow_MatchesBiopythonPi(string sequence, string scale, double expected)
    {
        var r = SequenceTools.IsoelectricPoint(sequence, scale, window: "BIOPYTHON");
        Assert.Multiple(() =>
        {
            Assert.That(r.PI, Is.EqualTo(expected).Within(1e-12));
            Assert.That(r.Found, Is.Null);
        });
    }

    [Test]
    public void IsoelectricPoint_NoTermini_MatchesIep()
    {
        // iep -notermini PETER → 4.0999.
        var r = SequenceTools.IsoelectricPoint("PETER", window: "Emboss", termini: false);
        Assert.That(r.PI, Is.EqualTo(4.0999).Within(5e-5));
    }

    [TestCase("PETER", "Default", true, 0, 0)]
    [TestCase("ACCKCKKA", "Default", false, 1, 1)]
    [TestCase("FLPVLAGLTPSIVPKLVCLLTKKC", "Emboss", true, 2, 1)]
    [TestCase("KDDK", "Biopython", false, 1, 0)]
    public void IsoelectricPoint_Options_EqualLibraryTryCalculate(string sequence, string window, bool termini, int lysines, int disulphides)
    {
        foreach (var scale in new[] { SequenceStatistics.PkaScale.Emboss, SequenceStatistics.PkaScale.Bjellqvist })
        {
            bool found = SequenceStatistics.TryCalculateIsoelectricPoint(sequence,
                new SequenceStatistics.IsoelectricPointOptions
                {
                    Scale = scale,
                    Window = Enum.Parse<SequenceStatistics.IsoelectricPointWindow>(window),
                    Termini = termini,
                    LysinesModified = lysines,
                    Disulphides = disulphides,
                },
                out double expected);
            var r = SequenceTools.IsoelectricPoint(sequence, scale.ToString(), window, termini, lysines, disulphides);
            Assert.That(r.PI, Is.EqualTo(found ? expected : (double?)null), $"{sequence} {scale}");
        }
    }

    [Test]
    public void IsoelectricPoint_DefaultOptions_JsonUnchanged()
    {
        string json = JsonSerializer.Serialize(SequenceTools.IsoelectricPoint("FLPVLAGLTPSIVPKLVCLLTKKC"), McpJsonUtilities.DefaultOptions);
        using var doc = JsonDocument.Parse(json);
        var props = doc.RootElement.EnumerateObject().ToArray();
        Assert.That(props, Has.Length.EqualTo(1));
        Assert.That(props[0].Value.GetDouble(), Is.EqualTo(9.57));
    }

    [Test]
    public void IsoelectricPoint_InvalidOptions_Throw()
    {
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("PETER", window: "None"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("PETER", window: "1"));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("PETER", lysinesModified: -1));
        Assert.Throws<ArgumentException>(() => SequenceTools.IsoelectricPoint("PETER", disulphides: -1));
    }

    // ---------------- Round 2 B4 (F48): melting_temperature weightAmbiguousBases ----------------
    // Biopython 1.88 Tm_Wallace(strict=False) executed 2026-10-14: 42.0, 39.0; the ≥ 14-base GC formula
    // with Biopython's non-strict GC weights (F33): 57.38333333333334.

    [TestCase("acgtNNSWacgt-ry", 42.0)]
    [TestCase("ACGTRYKMBDHVN", 39.0)]
    [TestCase("ACGTACGTACGTACGTRYKMBDHV", 57.38333333333334)]
    public void MeltingTemperature_WeightAmbiguousBases_MatchesBiopythonNonStrict(string sequence, double expected)
    {
        Assert.That(SequenceTools.MeltingTemperature(sequence, weightAmbiguousBases: true).Tm, Is.EqualTo(expected).Within(1e-9));
        Assert.That(SequenceTools.MeltingTemperature(sequence, weightAmbiguousBases: true).Tm,
            Is.EqualTo(SequenceTools.SummarizeSequence(sequence, weightAmbiguousBasesInTm: true).MeltingTemperature).Within(1e-12));
    }

    [TestCase("acgtNNSWacgt-ry")]
    [TestCase("ATGCGATCGATCGTAGCTAGCTAG")]
    public void MeltingTemperature_Default_BitIdenticalToLibrary(string sequence)
    {
        Assert.That(SequenceTools.MeltingTemperature(sequence).Tm, Is.EqualTo(SequenceStatistics.CalculateMeltingTemperature(sequence)));
        Assert.That(SequenceTools.MeltingTemperature(sequence, false).Tm, Is.EqualTo(SequenceStatistics.CalculateMeltingTemperature(sequence, false)));
    }

    private static string RandomString(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }
}
