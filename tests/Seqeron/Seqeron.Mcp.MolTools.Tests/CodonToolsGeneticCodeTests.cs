using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Optional <c>genetic_code</c> parameter of the codon-usage MCP tools (Validation 2026-09, B02 F34):
/// rscu, codon_adaptation_index, effective_number_of_codons, codon_usage_statistics.
/// Each tool must (a) keep the table-1 result when the parameter is omitted, (b) equal the
/// table-aware C# overload, (c) differ from table 1 where the code matters (NCBI table 2:
/// ATA = Met, TGA = Trp, AGA/AGG = stop), and (d) reject an unknown table number.
/// </summary>
[TestFixture]
public class CodonToolsGeneticCodeTests
{
    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    // Gene M3 of the ENC tests (CodonW 1.4.4 -enc: 41.29 under NCBI 1, 45.00 under NCBI 2).
    private static readonly string M3Gene =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // ATA×1, ATG×2, TGA×1, TGG×2, AGA×1, TAA×1.
    private const string MitoGene = "ATAATGATGTGATGGTGGAGATAA";

    private static IEnumerable<int> InvalidTables => new[] { 0, 7, 8, 17, 20, 34, -1 };

    #region rscu

    // Families from Biopython 1.88 CodonTable.unambiguous_dna_by_id[n] (forward_table + stop_codons),
    // RSCU = n·x / Σx (Sharp, Tuohy & Mosurski 1986). Table 2: Met {ATA, ATG}, Trp {TGA, TGG},
    // stops {TAA, TAG, AGA, AGG}; table 1: Ile {ATT, ATC, ATA}, Met {ATG}, stops {TAA, TAG, TGA},
    // Arg {CGN, AGA, AGG}.
    [TestCase(2, "ATA", 2.0 / 3)]
    [TestCase(2, "ATG", 4.0 / 3)]
    [TestCase(2, "TGA", 2.0 / 3)]
    [TestCase(2, "TGG", 4.0 / 3)]
    [TestCase(2, "AGA", 2.0)]
    [TestCase(2, "TAA", 2.0)]
    [TestCase(1, "ATA", 3.0)]
    [TestCase(1, "ATG", 1.0)]
    [TestCase(1, "TGA", 1.5)]
    [TestCase(1, "AGA", 6.0)]
    public void Rscu_GeneticCode_MatchesBiopythonFamilies(int table, string codon, double expected)
    {
        Assert.That(MolToolsTools.rscu(MitoGene, table).Rscu[codon], Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void Rscu_GeneticCode_EqualsCSharpOverload_DefaultIsTable1()
    {
        var mito = MolToolsTools.rscu(MitoGene, 2).Rscu;
        Assert.Multiple(() =>
        {
            Assert.That(mito, Is.EquivalentTo(CodonUsageAnalyzer.CalculateRscu(MitoGene, GeneticCode.GetByTableNumber(2))));
            Assert.That(MolToolsTools.rscu(MitoGene).Rscu, Is.EquivalentTo(CodonUsageAnalyzer.CalculateRscu(MitoGene)));
            Assert.That(MolToolsTools.rscu(MitoGene).Rscu, Is.EquivalentTo(MolToolsTools.rscu(MitoGene, 1).Rscu));
            Assert.That(mito, Is.Not.EquivalentTo(MolToolsTools.rscu(MitoGene).Rscu));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void Rscu_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.rscu(MitoGene, table), table);

    #endregion

    #region codon_adaptation_index

    // CodonW 1.4.4 cai_out with the Sharp & Li E. coli w values (-code 1 = NCBI 2 / -code 0 = NCBI 1);
    // same values locked for CodonUsageAnalyzer.CalculateCai(string, ..., GeneticCode) (F11/F18 tests).
    [TestCase("ATAATG", 2, 0.0547722557505166)]
    [TestCase("TGATGG", 2, 0.1)]
    [TestCase("AGAAGGCTG", 2, 1.0)]
    [TestCase("ATAATG", 1, 0.003)]
    [TestCase("TGATGG", 1, 0.0)]
    [TestCase("AGAAGGCTG", 1, 0.02)]
    public void CodonAdaptationIndex_GeneticCode_MatchesCodonW(string gene, int table, double expected)
    {
        double cai = MolToolsTools.codon_adaptation_index(gene, CodonUsageAnalyzer.EColiOptimalCodons, table).Cai;

        Assert.Multiple(() =>
        {
            Assert.That(cai, Is.EqualTo(expected).Within(1e-10));
            Assert.That(cai, Is.EqualTo(CodonUsageAnalyzer.CalculateCai(gene, CodonUsageAnalyzer.EColiOptimalCodons,
                GeneticCode.GetByTableNumber(table))));
        });
    }

    [Test]
    public void CodonAdaptationIndex_GeneticCode_ReferenceFamiliesFollowTheCode()
    {
        // Sharp & Li w = RSCU / max synonymous RSCU. Table 2: Met = {ATA, ATG} → w = {0.5, 1}, CAI = √0.5.
        // Table 1 (default): ATG is a single-codon family (excluded); ATA is Ile with family max 0.5 → w = 1.
        var reference = new Dictionary<string, double> { ["ATG"] = 1.0, ["ATA"] = 0.5 };

        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.codon_adaptation_index("ATAATG", reference, 2).Cai,
                Is.EqualTo(Math.Sqrt(0.5)).Within(1e-12));
            Assert.That(MolToolsTools.codon_adaptation_index("ATAATG", reference).Cai, Is.EqualTo(1.0).Within(1e-12));
            Assert.That(MolToolsTools.codon_adaptation_index("ATAATG", reference, 1).Cai, Is.EqualTo(1.0).Within(1e-12));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void CodonAdaptationIndex_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.codon_adaptation_index("ATAATG", CodonUsageAnalyzer.EColiOptimalCodons, table), table);

    #endregion

    #region effective_number_of_codons

    // CodonW 1.4.4 -enc: table 1 → 41.29, table 2 → 45.00, table 3 → 43.19, table 9 → 43.83
    // (locked for CalculateEnc(string, GeneticCode) in CodonUsageAnalyzer_CalculateEnc_Tests, F17).
    [TestCase(1, 41.288461538461526)]
    [TestCase(2, 45.0)]
    [TestCase(3, 43.18531468531467)]
    [TestCase(9, 43.83333333333333)]
    public void EffectiveNumberOfCodons_GeneticCode_MatchesCodonW(int table, double expected)
    {
        double enc = MolToolsTools.effective_number_of_codons(M3Gene, table).Enc;

        Assert.Multiple(() =>
        {
            Assert.That(enc, Is.EqualTo(expected).Within(1e-9));
            Assert.That(enc, Is.EqualTo(CodonUsageAnalyzer.CalculateEnc(M3Gene, GeneticCode.GetByTableNumber(table))));
        });
    }

    [Test]
    public void EffectiveNumberOfCodons_DefaultIsTable1()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.effective_number_of_codons(M3Gene).Enc,
                Is.EqualTo(MolToolsTools.effective_number_of_codons(M3Gene, 1).Enc));
            Assert.That(MolToolsTools.effective_number_of_codons(M3Gene).Enc,
                Is.EqualTo(CodonUsageAnalyzer.CalculateEnc(M3Gene)));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void EffectiveNumberOfCodons_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.effective_number_of_codons(M3Gene, table), table);

    #endregion

    #region codon_usage_statistics

    // CodonW 1.4.4 -gc3s (-code 0 = NCBI 1 / -code 1 = NCBI 2): ATAATG 0.000 / 0.500;
    // AGAAGGCTG 0.667 / 1.000 (locked for GetStatistics(string, GeneticCode), F18).
    [TestCase("ATAATG", 1, 0.0)]
    [TestCase("ATAATG", 2, 50.0)]
    [TestCase("AGAAGGCTG", 1, 200.0 / 3)]
    [TestCase("AGAAGGCTG", 2, 100.0)]
    public void CodonUsageStatistics_GeneticCode_Gc3sMatchesCodonW(string gene, int table, double gc3s)
    {
        var stats = MolToolsTools.codon_usage_statistics(gene, table);
        var expected = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.GetByTableNumber(table));

        Assert.Multiple(() =>
        {
            Assert.That(stats.Gc3s, Is.EqualTo(gc3s).Within(1e-9));
            Assert.That(stats.Gc3s, Is.EqualTo(expected.Gc3s));
            Assert.That(stats.Enc, Is.EqualTo(expected.Enc));
            Assert.That(stats.Rscu, Is.EquivalentTo(expected.Rscu));
            Assert.That(stats.CodonCounts, Is.EquivalentTo(expected.CodonCounts));
        });
    }

    [Test]
    public void CodonUsageStatistics_GeneticCode_RscuAndEncFollowTheCode_DefaultIsTable1()
    {
        var mito = MolToolsTools.codon_usage_statistics(M3Gene, 2);
        var standard = MolToolsTools.codon_usage_statistics(M3Gene);

        Assert.Multiple(() =>
        {
            Assert.That(mito.Enc, Is.EqualTo(45.0).Within(1e-9));
            Assert.That(standard.Enc, Is.EqualTo(41.288461538461526).Within(1e-9));
            Assert.That(standard.Rscu, Is.EquivalentTo(CodonUsageAnalyzer.GetStatistics(M3Gene).Rscu));
            Assert.That(MolToolsTools.codon_usage_statistics(M3Gene, 1).Gc3s, Is.EqualTo(standard.Gc3s));
            Assert.That(mito.TotalCodons, Is.EqualTo(standard.TotalCodons));
            Assert.That(mito.Gc3, Is.EqualTo(standard.Gc3));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void CodonUsageStatistics_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.codon_usage_statistics(M3Gene, table), table);

    #endregion

    #region count_codons

    // count_codons is code-independent and reads U as T (F10); a triplet with N is skipped frame-preservingly.
    [Test]
    public void CountCodons_RnaInput_CountedAsDna()
    {
        var counts = MolToolsTools.count_codons("augaugUUUNNNgcc").Counts;

        Assert.That(counts, Is.EquivalentTo(new Dictionary<string, int> { ["ATG"] = 2, ["TTT"] = 1, ["GCC"] = 1 }));
    }

    #endregion

    private static void AssertUnknownTable(TestDelegate call, int table)
    {
        var ex = Assert.Throws<ArgumentException>(call)!;
        Assert.Multiple(() =>
        {
            Assert.That(ex.ParamName, Is.EqualTo("genetic_code"));
            Assert.That(ex.Message, Does.Contain($"Unknown genetic code table: {table}."));
            Assert.That(ex.Message, Does.Contain(string.Join(", ", GeneticCode.SupportedTableNumbers)));
        });
    }
}
