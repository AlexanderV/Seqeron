namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for CODON-OPT-001 / CODON-RARE-001 / CODON-CAI-001: the additive, optional
/// <c>GeneticCode? geneticCode = null</c> parameter of every <c>CodonOptimizer</c> API whose
/// behaviour depends on the genetic code (review 2026-09, B02 F48, audit item A2-1).
///
/// Reference: DNA Chisel 3.2.16 <c>EnforceTranslation(genetic_table=…)</c> /
/// <c>CodonOptimize(…)</c>, whose synonymous families come from
/// <c>dnachisel.biotools.get_backtranslation_table(genetic_table)</c> — Biopython
/// <c>Bio.Data.CodonTable.unambiguous_dna_by_name[…].forward_table</c> grouped by amino acid,
/// with the table's <c>stop_codons</c> as the '*' family. Verified family sets (dnachisel 3.2.16,
/// Biopython 1.88):
/// <list type="bullet">
/// <item>table 2 "Vertebrate Mitochondrial": R = CGT CGC CGA CGG, M = ATA ATG, W = TGA TGG,
/// * = TAA TAG AGA AGG;</item>
/// <item>table 11 "Bacterial": identical <c>forward_table</c> and <c>stop_codons</c> to table 1
/// (only the start codons differ), so every optimizer output under table 11 equals the Standard one;</item>
/// <item>table 27 "Karyorelict Nuclear" (dual-coding): Q = TAA TAG CAA CAG, W = TGA TGG,
/// * = TGA — a context-dependent stop belongs to the family of the amino acid it encodes, exactly
/// as in Biopython's <c>forward_table</c> and in DNA Chisel's back-translation table.</item>
/// </list>
/// The <c>HarmonizeExpression</c> cross-check reuses the F23 approach: DNA Chisel
/// <c>DnaOptimizationProblem(constraints=[EnforceTranslation(genetic_table="Vertebrate
/// Mitochondrial")], objectives=[CodonOptimize(codon_usage_table=…, method="match_codon_usage")])</c>
/// scores our deterministic allocation −0.6, exactly the score its own randomised optimizer
/// reaches (input CGT×10: −11.8).
/// </summary>
[TestFixture]
[Category("Codon")]
[Category("CODON-OPT-001")]
public class CodonOptimizer_GeneticCode_Tests
{
    private static readonly GeneticCode Table2 = GeneticCode.GetByTableNumber(2);
    private static readonly GeneticCode Table11 = GeneticCode.GetByTableNumber(11);
    private static readonly GeneticCode Table27 = GeneticCode.GetByTableNumber(27);

    // AUG CGU AGA CGC AGG UGA AUA UAA — table 1: M R R R R * I *, table 2: M R * R * W M *.
    private const string MixedCodeGene = "AUGCGUAGACGCAGGUGAAUAUAA";

    // A usage table in which AGA is by far the most frequent Arg codon, so the Standard code picks
    // AGA while table 2 (AGA = stop) cannot: the family choice itself is observable.
    private static readonly CodonOptimizer.CodonUsageTable AgaRichArg = CodonOptimizer.CreateCodonUsageTable(
        "AGA-rich (synthetic)",
        new Dictionary<string, double>
        {
            { "AGA", 0.60 }, { "AGG", 0.05 }, { "CGU", 0.05 }, { "CGC", 0.20 }, { "CGA", 0.05 }, { "CGG", 0.05 },
            { "AUG", 1.00 }, { "AUA", 0.30 }, { "AUU", 0.40 }, { "AUC", 0.30 },
            { "UGG", 0.80 }, { "UGA", 0.20 },
            { "UAA", 0.60 }, { "UAG", 0.20 },
            { "CAA", 0.20 }, { "CAG", 0.30 },
            { "GAA", 0.60 }, { "GAG", 0.40 },
            { "AAA", 0.70 }, { "AAG", 0.30 },
        });

    #region OptimizeSequence

    [Test]
    [Description("Default (null) genetic code is bit-identical to an explicit Standard code")]
    public void OptimizeSequence_DefaultGeneticCode_EqualsExplicitStandard()
    {
        foreach (var strategy in new[]
                 {
                     CodonOptimizer.OptimizationStrategy.MaximizeCAI,
                     CodonOptimizer.OptimizationStrategy.BalancedOptimization,
                     CodonOptimizer.OptimizationStrategy.AvoidRareCodeons,
                     CodonOptimizer.OptimizationStrategy.HarmonizeExpression,
                 })
        {
            var standard = CodonOptimizer.OptimizeSequence(
                MixedCodeGene, CodonOptimizer.EColiK12, strategy, geneticCode: GeneticCode.Standard);
            var dflt = CodonOptimizer.OptimizeSequence(MixedCodeGene, CodonOptimizer.EColiK12, strategy);

            Assert.Multiple(() =>
            {
                Assert.That(standard.OptimizedSequence, Is.EqualTo(dflt.OptimizedSequence), $"{strategy} sequence");
                Assert.That(standard.ProteinSequence, Is.EqualTo(dflt.ProteinSequence), $"{strategy} protein");
                Assert.That(standard.OriginalCAI, Is.EqualTo(dflt.OriginalCAI), $"{strategy} original CAI");
                Assert.That(standard.OptimizedCAI, Is.EqualTo(dflt.OptimizedCAI), $"{strategy} optimized CAI");
                Assert.That(standard.ChangedCodons, Is.EqualTo(dflt.ChangedCodons), $"{strategy} change count");
            });
        }
    }

    [Test]
    [Description("Table 2: AGA/AGG are stops (never Arg), UGA is Trp and AUA is Met")]
    public void OptimizeSequence_VertebrateMitochondrial_UsesTable2Families()
    {
        var standard = CodonOptimizer.OptimizeSequence(
            MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var mito = CodonOptimizer.OptimizeSequence(
            MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI,
            geneticCode: Table2);

        Assert.Multiple(() =>
        {
            // Standard: every Arg codon (CGU, AGA, CGC, AGG) becomes the most frequent Arg codon
            // CGC (0.40); UGA is a stop and is kept; AUA (Ile) becomes AUU (0.51).
            Assert.That(standard.OptimizedSequence, Is.EqualTo("AUGCGCCGCCGCCGCUGAAUUUAA"));
            Assert.That(standard.ProteinSequence, Is.EqualTo("MRRRR*I*"));

            // Table 2: AGA/AGG are stops and are preserved; the Arg codons CGU/CGC go to CGC;
            // UGA is Trp and becomes UGG (1.00 > 0.29); AUA is Met and becomes AUG.
            Assert.That(mito.OptimizedSequence, Is.EqualTo("AUGCGCAGACGCAGGUGGAUGUAA"));
            Assert.That(mito.ProteinSequence, Is.EqualTo("MR*R*WM*"));

            // The protein is preserved under the code that was asked for (Translator).
            Assert.That(Translator.Translate(mito.OptimizedSequence, Table2).Sequence,
                Is.EqualTo(Translator.Translate(MixedCodeGene, Table2).Sequence));
            Assert.That(Translator.Translate(standard.OptimizedSequence, GeneticCode.Standard).Sequence,
                Is.EqualTo(Translator.Translate(MixedCodeGene, GeneticCode.Standard).Sequence));
        });
    }

    [Test]
    [Description("Table 2 never emits AGA/AGG for Arg even when the usage table makes AGA the best codon")]
    public void OptimizeSequence_VertebrateMitochondrial_NeverEmitsAgaAggForArg()
    {
        const string argGene = "AUGCGUCGUCGUCGUUAA";

        var standard = CodonOptimizer.OptimizeSequence(
            argGene, AgaRichArg, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var mito = CodonOptimizer.OptimizeSequence(
            argGene, AgaRichArg, CodonOptimizer.OptimizationStrategy.MaximizeCAI, geneticCode: Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard.OptimizedSequence, Is.EqualTo("AUGAGAAGAAGAAGAUAA"),
                "the Standard code may use AGA, the most frequent Arg codon of this table");
            Assert.That(mito.OptimizedSequence, Is.EqualTo("AUGCGCCGCCGCCGCUAA"),
                "under table 2 AGA/AGG are stop codons, so the best Arg codon is CGC");
            Assert.That(Translator.Translate(mito.OptimizedSequence, Table2).Sequence,
                Is.EqualTo("MRRRR*"));
        });
    }

    [Test]
    [Description("Table 11 has table 1's codon table, so the optimizer output is identical")]
    public void OptimizeSequence_BacterialPlastid_EqualsStandardOutput()
    {
        var standard = CodonOptimizer.OptimizeSequence(
            MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization);
        var bacterial = CodonOptimizer.OptimizeSequence(
            MixedCodeGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.BalancedOptimization,
            geneticCode: Table11);

        Assert.That(bacterial.OptimizedSequence, Is.EqualTo(standard.OptimizedSequence));
        Assert.That(bacterial.ProteinSequence, Is.EqualTo(standard.ProteinSequence));
    }

    [Test]
    [Description("HarmonizeExpression under table 2 matches the DNA Chisel match_codon_usage optimum")]
    public void OptimizeSequence_HarmonizeExpression_Table2_MatchesDnaChiselScore()
    {
        // 10 Arg codons. Table-2 Arg family CGU .38 CGC .40 CGA .06 CGG .10 (total .94):
        // quotas 4.043 / 4.255 / 0.638 / 1.064 → floors 4,4,0,1 (sum 9), the remaining seat goes
        // to the largest fraction (CGA) → CGU4 CGC4 CGA1 CGG1. The four positions that already
        // hold CGU keep it. DNA Chisel scores this allocation −0.6 under
        // EnforceTranslation(genetic_table="Vertebrate Mitochondrial") — the same score its own
        // randomised match_codon_usage optimizer reaches (input CGU×10 scores −11.8).
        const string argGene = "CGUCGUCGUCGUCGUCGUCGUCGUCGUCGU";

        var result = CodonOptimizer.OptimizeSequence(
            argGene, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.HarmonizeExpression,
            geneticCode: Table2);

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Is.EqualTo("CGUCGUCGUCGUCGCCGCCGCCGCCGACGG"));
            Assert.That(Translator.Translate(result.OptimizedSequence, Table2).Sequence, Is.EqualTo("RRRRRRRRRR"));
        });
    }

    [Test]
    [Description("Dual-coding tables 27/28/31: the codon belongs to the amino-acid family (Biopython / DNA Chisel)")]
    public void OptimizeSequence_DualCodingTable27_TreatsStopAsItsAminoAcid()
    {
        // Table 27: TAA/TAG encode Gln (and TGA is the only stop). DNA Chisel's
        // get_backtranslation_table("Karyorelict Nuclear")["Q"] == ["TAA", "TAG", "CAA", "CAG"],
        // so a Gln codon may be rewritten to UAA/UAG and vice versa; with this table the most
        // frequent Gln codon of the family is UAA (0.60).
        const string glnGene = "AUGCAACAGUGA";

        var result = CodonOptimizer.OptimizeSequence(
            glnGene, AgaRichArg, CodonOptimizer.OptimizationStrategy.MaximizeCAI, geneticCode: Table27);

        Assert.Multiple(() =>
        {
            Assert.That(result.ProteinSequence, Is.EqualTo("MQQW"),
                "under table 27 UGA translates as Trp and UAA/UAG as Gln (Biopython forward_table)");
            Assert.That(result.OptimizedSequence, Is.EqualTo("AUGUAAUAAUGG"));
            Assert.That(Translator.Translate(result.OptimizedSequence, Table27).Sequence, Is.EqualTo("MQQW"));
        });
    }

    #endregion

    #region CalculateCAI

    [Test]
    [Description("CalculateCAI: default == explicit Standard, and table N == the canonical CAI core under table N")]
    public void CalculateCAI_GeneticCode_MatchesCanonicalCoreUnderThatCode()
    {
        var reference = new Dictionary<string, double>();
        foreach (var (codon, frequency) in CodonOptimizer.EColiK12.CodonFrequencies)
            reference[codon.Replace('U', 'T')] = frequency;

        double dflt = CodonOptimizer.CalculateCAI(MixedCodeGene, CodonOptimizer.EColiK12);
        double standard = CodonOptimizer.CalculateCAI(
            MixedCodeGene, CodonOptimizer.EColiK12, geneticCode: GeneticCode.Standard);
        double mito = CodonOptimizer.CalculateCAI(
            MixedCodeGene, CodonOptimizer.EColiK12, geneticCode: Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard, Is.EqualTo(dflt), "the default must be the Standard code");
            Assert.That(dflt, Is.EqualTo(
                CodonUsageAnalyzer.CalculateCai(MixedCodeGene, reference, GeneticCode.Standard)).Within(1e-12));
            Assert.That(mito, Is.EqualTo(
                CodonUsageAnalyzer.CalculateCai(MixedCodeGene, reference, Table2)).Within(1e-12),
                "CAI under table N must equal the CodonW-locked canonical core under table N");
            Assert.That(mito, Is.Not.EqualTo(dflt).Within(1e-9),
                "table 2 changes the stops and the Met/Trp families, so the scored codon set differs");
        });
    }

    #endregion

    #region RemoveRestrictionSites

    [Test]
    [Description("RemoveRestrictionSites: default == Standard; table 2 substitutions stay inside table-2 families")]
    public void RemoveRestrictionSites_GeneticCode_UsesThatCodesFamilies()
    {
        const string gene = "AUGCGUCGUCGUUAA";
        string[] site = { "CGTCGT" };

        string dflt = CodonOptimizer.RemoveRestrictionSites(gene, site, AgaRichArg);
        string standard = CodonOptimizer.RemoveRestrictionSites(gene, site, AgaRichArg, GeneticCode.Standard);
        string mito = CodonOptimizer.RemoveRestrictionSites(gene, site, AgaRichArg, Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard, Is.EqualTo(dflt));
            Assert.That(dflt, Does.Contain("AGA"), "the Standard code removes the site with the AGA codon (0.60)");
            Assert.That(mito, Does.Not.Contain("AGA").And.Not.Contain("AGG"),
                "under table 2 AGA/AGG are stops and must never replace an Arg codon");
            Assert.That(Translator.Translate(mito, Table2).Sequence,
                Is.EqualTo(Translator.Translate(gene, Table2).Sequence),
                "the protein under the requested code is preserved");
            Assert.That(mito.Replace('U', 'T'), Does.Not.Contain("CGTCGT"));
        });
    }

    #endregion

    #region ReduceSecondaryStructure

    [Test]
    [Description("ReduceSecondaryStructure: default == Standard; table 2 keeps the table-2 protein")]
    public void ReduceSecondaryStructure_GeneticCode_PreservesProteinUnderThatCode()
    {
        // 48 nt, Arg/Gln/Glu/Lys codons plus AUG and a UGA (Trp under table 2).
        const string gene = "AUGCGUCGUCGUCAACAGGAAGAGAAAAAACGUCGUUGAUAA";

        string dflt = CodonOptimizer.ReduceSecondaryStructure(gene, AgaRichArg, windowSize: 21);
        string standard = CodonOptimizer.ReduceSecondaryStructure(
            gene, AgaRichArg, windowSize: 21, geneticCode: GeneticCode.Standard);
        string mito = CodonOptimizer.ReduceSecondaryStructure(
            gene, AgaRichArg, windowSize: 21, geneticCode: Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard, Is.EqualTo(dflt));
            Assert.That(Translator.Translate(dflt, GeneticCode.Standard).Sequence,
                Is.EqualTo(Translator.Translate(gene, GeneticCode.Standard).Sequence));
            Assert.That(Translator.Translate(mito, Table2).Sequence,
                Is.EqualTo(Translator.Translate(gene, Table2).Sequence),
                "only synonymous codons of the requested code may be used");
            for (int i = 0; i + 3 <= mito.Length; i += 3)
            {
                string codon = mito.Substring(i, 3);
                if (Table2.Translate(codon) == 'R')
                    Assert.That(codon, Is.AnyOf("CGU", "CGC", "CGA", "CGG"), $"codon {i / 3}");
            }
        });
    }

    #endregion

    #region FindRareCodons / %MinMax

    [Test]
    [Description("FindRareCodons: the rare set is code-independent, the reported amino acid is not")]
    public void FindRareCodons_GeneticCode_ChangesOnlyTheAminoAcid()
    {
        var dflt = CodonOptimizer.FindRareCodons(MixedCodeGene, CodonOptimizer.EColiK12, 0.30).ToList();
        var standard = CodonOptimizer.FindRareCodons(
            MixedCodeGene, CodonOptimizer.EColiK12, 0.30, geneticCode: GeneticCode.Standard).ToList();
        var mito = CodonOptimizer.FindRareCodons(
            MixedCodeGene, CodonOptimizer.EColiK12, 0.30, geneticCode: Table2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(standard, Is.EqualTo(dflt));
            Assert.That(mito.Select(r => r.Position), Is.EqualTo(dflt.Select(r => r.Position)),
                "rarity is a per-codon table lookup and does not depend on the code");
            // UGA: frequency 0.29 < 0.30 → reported; '*' under table 1, 'W' under table 2.
            Assert.That(dflt.Single(r => r.Codon == "UGA").AminoAcid, Is.EqualTo("*"));
            Assert.That(mito.Single(r => r.Codon == "UGA").AminoAcid, Is.EqualTo("W"));
            Assert.That(mito.Single(r => r.Codon == "AGA").AminoAcid, Is.EqualTo("*"));
        });
    }

    [Test]
    [Description("%MinMax: default == Standard; table 2 changes the families and hence the profile")]
    public void CalculateMinMaxProfile_GeneticCode_UsesThatCodesFamilies()
    {
        const string gene = "AUGCGUAGACGCAGGUGAAUAAAACAACGUCGCAGAAGGUGGAUGAAAUAA";

        var dflt = CodonOptimizer.CalculateMinMaxProfile(gene, CodonOptimizer.Human, 5);
        var standard = CodonOptimizer.CalculateMinMaxProfile(
            gene, CodonOptimizer.Human, 5, geneticCode: GeneticCode.Standard);
        var mito = CodonOptimizer.CalculateMinMaxProfile(gene, CodonOptimizer.Human, 5, geneticCode: Table2);
        var basisAndCode = CodonOptimizer.CalculateMinMaxProfile(
            gene, CodonOptimizer.Human, 5, CodonOptimizer.MinMaxUsageBasis.PerThousandWhenAvailable, Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard.Select(w => w.PercentMinMax), Is.EqualTo(dflt.Select(w => w.PercentMinMax)));
            Assert.That(basisAndCode.Select(w => w.PercentMinMax), Is.EqualTo(mito.Select(w => w.PercentMinMax)),
                "the 5-argument overload must agree with the 4-argument one at the default basis");
            Assert.That(mito, Has.Count.EqualTo(dflt.Count));
            Assert.That(mito.Select(w => w.PercentMinMax), Is.Not.EqualTo(dflt.Select(w => w.PercentMinMax)),
                "the Arg family loses AGA/AGG and Met/Trp gain a second codon under table 2");
        });
    }

    #endregion

    #region Table construction

    [Test]
    [Description("CreateCodonTableFromSequence: families and codon→AA map follow the requested code")]
    public void CreateCodonTableFromSequence_GeneticCode_NormalisesWithinThatCodesFamilies()
    {
        const string reference = "AUGCGUCGUCGCAGAAGGAUAAUGUGAUGGUAAUAG";

        var dflt = CodonOptimizer.CreateCodonTableFromSequence(reference, "default");
        var standard = CodonOptimizer.CreateCodonTableFromSequence(reference, "standard", GeneticCode.Standard);
        var mito = CodonOptimizer.CreateCodonTableFromSequence(reference, "mito", Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard.CodonFrequencies, Is.EqualTo(dflt.CodonFrequencies));
            Assert.That(standard.CodonToAminoAcid, Is.EqualTo(dflt.CodonToAminoAcid));

            // Table 2 codon → amino acid.
            Assert.That(mito.CodonToAminoAcid["AGA"], Is.EqualTo("*"));
            Assert.That(mito.CodonToAminoAcid["AUA"], Is.EqualTo("M"));
            Assert.That(mito.CodonToAminoAcid["UGA"], Is.EqualTo("W"));

            // Arg family of table 2 has four codons and its frequencies sum to 1.
            string[] mitoArg = { "CGU", "CGC", "CGA", "CGG" };
            Assert.That(mitoArg.Sum(c => mito.CodonFrequencies[c]), Is.EqualTo(1.0).Within(1e-12));
            // Standard Arg family has six and includes AGA/AGG.
            string[] standardArg = { "CGU", "CGC", "CGA", "CGG", "AGA", "AGG" };
            Assert.That(standardArg.Sum(c => standard.CodonFrequencies[c]), Is.EqualTo(1.0).Within(1e-12));
            // Table-2 stop family {UAA, UAG, AGA, AGG} sums to 1.
            string[] mitoStops = { "UAA", "UAG", "AGA", "AGG" };
            Assert.That(mitoStops.Sum(c => mito.CodonFrequencies[c]), Is.EqualTo(1.0).Within(1e-12));
            // Table-2 Trp family {UGA, UGG}: both observed once → 0.5 each.
            Assert.That(mito.CodonFrequencies["UGA"], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(mito.CodonFrequencies["UGG"], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(mito.CodonFrequencies, Has.Count.EqualTo(64));
        });
    }

    [Test]
    [Description("CreateCodonUsageTable: the codon→AA map comes from the requested code")]
    public void CreateCodonUsageTable_GeneticCode_StoresThatCodesMapping()
    {
        var frequencies = new Dictionary<string, double> { { "AGA", 1.0 } };

        var dflt = CodonOptimizer.CreateCodonUsageTable("default", frequencies);
        var standard = CodonOptimizer.CreateCodonUsageTable("standard", frequencies, GeneticCode.Standard);
        var mito = CodonOptimizer.CreateCodonUsageTable("mito", frequencies, Table2);

        Assert.Multiple(() =>
        {
            Assert.That(standard.CodonToAminoAcid, Is.EqualTo(dflt.CodonToAminoAcid));
            Assert.That(dflt.CodonToAminoAcid["AGA"], Is.EqualTo("R"));
            Assert.That(mito.CodonToAminoAcid["AGA"], Is.EqualTo("*"));
            Assert.That(mito.CodonToAminoAcid["UGA"], Is.EqualTo("W"));
        });
    }

    #endregion
}
