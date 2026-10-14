// CODON-OPT/RARE/CAI — MCP: optional genetic_code on the CodonOptimizer tools (FIN-B02 WP15, audit A2-2;
// B02 F49). The tools delegate to the genetic-code-aware CodonOptimizer API of F48 (commit 353f809);
// the table-2 locks below are those of
// tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_GeneticCode_Tests.cs (DNA Chisel
// EnforceTranslation(genetic_table=…) / Biopython forward_table families). Each tool must (a) keep the
// old result when genetic_code is omitted (== explicit 1 == the C# overload without a code), (b) equal
// the C# overload under table 2, (c) differ from table 1 where the code matters, (d) reject an unknown
// table with the F34 error.
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CodonOptimizerToolsGeneticCodeTests
{
    private static readonly GeneticCode Table2 = GeneticCode.GetByTableNumber(2);
    private static readonly CodonUsageTableInput EColi = new(Preset: "EColiK12");
    private static readonly CodonUsageTableInput Human = new(Preset: "Human");

    // AUG CGU AGA CGC AGG UGA AUA UAA — table 1: M R R R R * I *, table 2: M R * R * W M *.
    private const string MixedCodeGene = "AUGCGUAGACGCAGGUGAAUAUAA";

    private static readonly Dictionary<string, double> AgaRichFrequencies = new()
    {
        { "AGA", 0.60 }, { "AGG", 0.05 }, { "CGU", 0.05 }, { "CGC", 0.20 }, { "CGA", 0.05 }, { "CGG", 0.05 },
        { "AUG", 1.00 }, { "AUA", 0.30 }, { "AUU", 0.40 }, { "AUC", 0.30 },
        { "UGG", 0.80 }, { "UGA", 0.20 },
        { "UAA", 0.60 }, { "UAG", 0.20 },
        { "CAA", 0.20 }, { "CAG", 0.30 },
        { "GAA", 0.60 }, { "GAG", 0.40 },
        { "AAA", 0.70 }, { "AAG", 0.30 },
    };

    private static CodonUsageTableInput AgaRichInput => new(OrganismName: "AGA-rich (synthetic)", CodonFrequencies: AgaRichFrequencies);
    private static CodonOptimizer.CodonUsageTable AgaRichTable => CodonOptimizer.CreateCodonUsageTable("AGA-rich (synthetic)", AgaRichFrequencies);

    private static IEnumerable<int> InvalidTables => new[] { 0, 7, 8, 17, 20, 34, -1 };

    #region optimize_codons

    [Test]
    public void OptimizeCodons_DefaultEqualsTable1_Table2IsTheLockedCSharpResult()
    {
        var dflt = MolToolsTools.optimize_codons(MixedCodeGene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var t1 = MolToolsTools.optimize_codons(MixedCodeGene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI, genetic_code: 1);
        var t2 = MolToolsTools.optimize_codons(MixedCodeGene, EColi, CodonOptimizer.OptimizationStrategy.MaximizeCAI, genetic_code: 2);
        var api1 = CodonOptimizer.OptimizeSequence(MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var api2 = CodonOptimizer.OptimizeSequence(MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI, geneticCode: Table2);

        Assert.Multiple(() =>
        {
            // F48 locks: table 1 -> AUGCGCCGCCGCCGCUGAAUUUAA; table 2 keeps AGA/AGG (stops), UGA -> UGG, AUA -> AUG.
            Assert.That(dflt.OptimizedSequence, Is.EqualTo("AUGCGCCGCCGCCGCUGAAUUUAA"));
            Assert.That(t1, Is.EqualTo(dflt).Using<OptimizationResultDto>(SameOptimization));
            Assert.That(dflt.OptimizedSequence, Is.EqualTo(api1.OptimizedSequence));
            Assert.That(t2.OptimizedSequence, Is.EqualTo("AUGCGCAGACGCAGGUGGAUGUAA"));
            Assert.That(t2.OptimizedSequence, Is.EqualTo(api2.OptimizedSequence));
            Assert.That(t2.ProteinSequence, Is.EqualTo(api2.ProteinSequence));
            Assert.That(t2.ProteinSequence, Is.EqualTo("MR*R*WM*"));
            Assert.That(t2.OptimizedCAI, Is.EqualTo(api2.OptimizedCAI));
            Assert.That(t2.ChangedCodons, Is.EqualTo(api2.ChangedCodons));
        });
    }

    [Test]
    public void OptimizeCodons_InlineTable_Table2NeverUsesAgaAggForArg()
    {
        const string argGene = "AUGCGUCGUCGUCGUUAA";
        var t1 = MolToolsTools.optimize_codons(argGene, AgaRichInput, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var t2 = MolToolsTools.optimize_codons(argGene, AgaRichInput, CodonOptimizer.OptimizationStrategy.MaximizeCAI, genetic_code: 2);
        Assert.Multiple(() =>
        {
            Assert.That(t1.OptimizedSequence, Is.EqualTo("AUGAGAAGAAGAAGAUAA"));
            Assert.That(t2.OptimizedSequence, Is.EqualTo("AUGCGCCGCCGCCGCUAA"));
            Assert.That(t2.OptimizedSequence, Is.EqualTo(CodonOptimizer.OptimizeSequence(
                argGene, AgaRichTable, CodonOptimizer.OptimizationStrategy.MaximizeCAI, geneticCode: Table2).OptimizedSequence));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void OptimizeCodons_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.optimize_codons(MixedCodeGene, EColi, genetic_code: table), table);

    #endregion

    #region cai_from_organism_table

    [Test]
    public void CaiFromOrganismTable_GeneticCode_EqualsCSharpOverload()
    {
        double dflt = MolToolsTools.cai_from_organism_table(MixedCodeGene, EColi).Cai;
        double t1 = MolToolsTools.cai_from_organism_table(MixedCodeGene, EColi, 1).Cai;
        double t2 = MolToolsTools.cai_from_organism_table(MixedCodeGene, EColi, 2).Cai;
        Assert.Multiple(() =>
        {
            Assert.That(t1, Is.EqualTo(dflt));
            Assert.That(dflt, Is.EqualTo(CodonOptimizer.CalculateCAI(MixedCodeGene, CodonOptimizer.EColiK12)));
            Assert.That(t2, Is.EqualTo(CodonOptimizer.CalculateCAI(MixedCodeGene, CodonOptimizer.EColiK12, geneticCode: Table2)));
            Assert.That(t2, Is.Not.EqualTo(dflt), "AGA/AGG drop out as stops, UGA/AUA join Trp/Met under table 2");
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void CaiFromOrganismTable_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.cai_from_organism_table(MixedCodeGene, EColi, table), table);

    #endregion

    #region remove_restriction_sites

    [Test]
    public void RemoveRestrictionSites_GeneticCode_UsesThatCodesFamilies()
    {
        const string gene = "AUGCGUCGUCGUUAA";
        string[] site = { "CGTCGT" };
        string dflt = MolToolsTools.remove_restriction_sites(gene, site, AgaRichInput).OptimizedSequence;
        string t1 = MolToolsTools.remove_restriction_sites(gene, site, AgaRichInput, 1).OptimizedSequence;
        string t2 = MolToolsTools.remove_restriction_sites(gene, site, AgaRichInput, 2).OptimizedSequence;
        Assert.Multiple(() =>
        {
            Assert.That(t1, Is.EqualTo(dflt));
            Assert.That(dflt, Is.EqualTo(CodonOptimizer.RemoveRestrictionSites(gene, site, AgaRichTable)));
            Assert.That(dflt, Does.Contain("AGA"));
            Assert.That(t2, Is.EqualTo(CodonOptimizer.RemoveRestrictionSites(gene, site, AgaRichTable, Table2)));
            Assert.That(t2, Does.Not.Contain("AGA").And.Not.Contain("AGG"));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void RemoveRestrictionSites_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.remove_restriction_sites("AUGCGUCGUUAA", new[] { "CGTCGT" }, EColi, table), table);

    #endregion

    #region reduce_secondary_structure

    [Test]
    public void ReduceSecondaryStructure_GeneticCode_EqualsCSharpOverload()
    {
        const string gene = "AUGCGUCGUCGUCAACAGGAAGAGAAAAAACGUCGUUGAUAA";
        string dflt = MolToolsTools.reduce_secondary_structure(gene, AgaRichInput, 21).OptimizedSequence;
        string t1 = MolToolsTools.reduce_secondary_structure(gene, AgaRichInput, 21, genetic_code: 1).OptimizedSequence;
        string t2 = MolToolsTools.reduce_secondary_structure(gene, AgaRichInput, 21, genetic_code: 2).OptimizedSequence;
        Assert.Multiple(() =>
        {
            Assert.That(t1, Is.EqualTo(dflt));
            Assert.That(dflt, Is.EqualTo(CodonOptimizer.ReduceSecondaryStructure(gene, AgaRichTable, 21)));
            Assert.That(t2, Is.EqualTo(CodonOptimizer.ReduceSecondaryStructure(gene, AgaRichTable, 21, geneticCode: Table2)));
            Assert.That(Translator.Translate(t2, Table2).Sequence, Is.EqualTo(Translator.Translate(gene, Table2).Sequence));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void ReduceSecondaryStructure_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.reduce_secondary_structure(MixedCodeGene, EColi, genetic_code: table), table);

    #endregion

    #region find_rare_codons

    [Test]
    public void FindRareCodons_GeneticCode_ChangesOnlyTheAminoAcid()
    {
        var dflt = MolToolsTools.find_rare_codons(MixedCodeGene, EColi, 0.30).RareCodons;
        var t1 = MolToolsTools.find_rare_codons(MixedCodeGene, EColi, 0.30, genetic_code: 1).RareCodons;
        var t2 = MolToolsTools.find_rare_codons(MixedCodeGene, EColi, 0.30, genetic_code: 2).RareCodons;
        var api2 = CodonOptimizer.FindRareCodons(MixedCodeGene, CodonOptimizer.EColiK12, 0.30, geneticCode: Table2)
            .Select(t => new RareCodon(t.Position, t.Codon, t.AminoAcid, t.Frequency));
        var api1 = CodonOptimizer.FindRareCodons(MixedCodeGene, CodonOptimizer.EColiK12, 0.30)
            .Select(t => new RareCodon(t.Position, t.Codon, t.AminoAcid, t.Frequency));
        Assert.Multiple(() =>
        {
            Assert.That(t1, Is.EqualTo(dflt));
            Assert.That(dflt, Is.EqualTo(api1));
            Assert.That(t2, Is.EqualTo(api2));
            Assert.That(t2.Select(r => r.Position), Is.EqualTo(dflt.Select(r => r.Position)));
            Assert.That(dflt.Single(r => r.Codon == "UGA").AminoAcid, Is.EqualTo("*"));
            Assert.That(t2.Single(r => r.Codon == "UGA").AminoAcid, Is.EqualTo("W"));
            Assert.That(t2.Single(r => r.Codon == "AGA").AminoAcid, Is.EqualTo("*"));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void FindRareCodons_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.find_rare_codons(MixedCodeGene, EColi, genetic_code: table), table);

    #endregion

    #region codon_minmax_profile

    [Test]
    public void CodonMinMaxProfile_GeneticCode_EqualsCSharpOverload()
    {
        const string gene = "AUGCGUAGACGCAGGUGAAUAAAACAACGUCGCAGAAGGUGGAUGAAAUAA";
        var dflt = MolToolsTools.codon_minmax_profile(gene, Human, 5).Windows.Select(w => w.PercentMinMax).ToList();
        var t1 = MolToolsTools.codon_minmax_profile(gene, Human, 5, genetic_code: 1).Windows.Select(w => w.PercentMinMax).ToList();
        var t2 = MolToolsTools.codon_minmax_profile(gene, Human, 5, genetic_code: 2).Windows.Select(w => w.PercentMinMax).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(t1, Is.EqualTo(dflt));
            Assert.That(dflt, Is.EqualTo(CodonOptimizer.CalculateMinMaxProfile(gene, CodonOptimizer.Human, 5).Select(w => w.PercentMinMax)));
            Assert.That(t2, Is.EqualTo(CodonOptimizer.CalculateMinMaxProfile(gene, CodonOptimizer.Human, 5, geneticCode: Table2).Select(w => w.PercentMinMax)));
            Assert.That(t2, Is.Not.EqualTo(dflt));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void CodonMinMaxProfile_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.codon_minmax_profile(MixedCodeGene, Human, 3, genetic_code: table), table);

    #endregion

    #region build_codon_table

    [Test]
    public void BuildCodonTable_GeneticCode_EqualsCSharpOverload()
    {
        const string reference = "AUGCGUCGUCGCAGAAGGAUAAUGUGAUGGUAAUAG";
        var dflt = MolToolsTools.build_codon_table(reference, "ref");
        var t1 = MolToolsTools.build_codon_table(reference, "ref", 1);
        var t2 = MolToolsTools.build_codon_table(reference, "ref", 2);
        var api1 = CodonOptimizer.CreateCodonTableFromSequence(reference, "ref");
        var api2 = CodonOptimizer.CreateCodonTableFromSequence(reference, "ref", Table2);
        Assert.Multiple(() =>
        {
            Assert.That(t1.CodonFrequencies, Is.EquivalentTo(dflt.CodonFrequencies));
            Assert.That(dflt.CodonFrequencies, Is.EquivalentTo(api1.CodonFrequencies));
            Assert.That(t2.CodonFrequencies, Is.EquivalentTo(api2.CodonFrequencies));
            Assert.That(t2.CodonToAminoAcid, Is.EquivalentTo(api2.CodonToAminoAcid));
            Assert.That(t2.CodonToAminoAcid["UGA"], Is.EqualTo("W"));
            Assert.That(dflt.CodonToAminoAcid["UGA"], Is.EqualTo("*"));
            Assert.That(t2.CodonFrequencies, Is.Not.EquivalentTo(dflt.CodonFrequencies));
        });
    }

    [TestCaseSource(nameof(InvalidTables))]
    public void BuildCodonTable_UnknownGeneticCode_Throws(int table) =>
        AssertUnknownTable(() => MolToolsTools.build_codon_table(MixedCodeGene, "ref", table), table);

    #endregion

    private static bool SameOptimization(OptimizationResultDto a, OptimizationResultDto b) =>
        a.OriginalSequence == b.OriginalSequence && a.OptimizedSequence == b.OptimizedSequence &&
        a.ProteinSequence == b.ProteinSequence && a.OriginalCAI.Equals(b.OriginalCAI) &&
        a.OptimizedCAI.Equals(b.OptimizedCAI) && a.GcContentOriginal.Equals(b.GcContentOriginal) &&
        a.GcContentOptimized.Equals(b.GcContentOptimized) && a.ChangedCodons == b.ChangedCodons &&
        a.Changes.SequenceEqual(b.Changes);

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
