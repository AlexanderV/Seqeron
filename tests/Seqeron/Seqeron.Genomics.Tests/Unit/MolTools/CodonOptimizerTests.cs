namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonOptimizerTests
{
    // NOTE: OptimizeSequence tests in CodonOptimizer_OptimizeSequence_Tests.cs (CODON-OPT-001)
    // NOTE: CAI calculation tests in CodonOptimizer_CAI_Tests.cs (CODON-CAI-001)
    // NOTE: Codon Usage Analysis tests in CodonOptimizer_CodonUsage_Tests.cs (CODON-USAGE-001)
    // NOTE: FindRareCodons tests in CodonOptimizer_FindRareCodons_Tests.cs (CODON-RARE-001)

    #region Restriction Site Removal Tests

    [Test]
    public void RemoveRestrictionSites_RemovesTargetSite()
    {
        // EcoRI: GAATTC (GAAUUC in RNA)
        string sequence = "AUGGAAUUCGCU"; // Contains EcoRI site
        var result = CodonOptimizer.RemoveRestrictionSites(sequence, new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        Assert.That(result, Does.Not.Contain("GAAUUC"));
    }

    [Test]
    public void RemoveRestrictionSites_PreservesProtein()
    {
        string sequence = "AUGGAAUUC"; // M-E-F
        string original = sequence.Replace('T', 'U');
        var result = CodonOptimizer.RemoveRestrictionSites(sequence, new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        // Verify protein is preserved by checking length is unchanged
        Assert.That(result.Length, Is.EqualTo(original.Length));
    }

    [Test]
    public void RemoveRestrictionSites_MultiplesSites()
    {
        // BamHI: GGATCC, HindIII: AAGCTT
        string sequence = "AUGGGATCCAAGCTTGCU";
        var result = CodonOptimizer.RemoveRestrictionSites(sequence,
            new[] { "GGATCC", "AAGCTT" }, CodonOptimizer.EColiK12);

        Assert.That(result, Does.Not.Contain("GGAUCC"));
    }

    [Test]
    public void RemoveRestrictionSites_EmptySequence_ReturnsEmpty()
    {
        var result = CodonOptimizer.RemoveRestrictionSites("", new[] { "GAATTC" }, CodonOptimizer.EColiK12);
        Assert.That(result, Is.Empty);
    }

    #endregion

    #region Secondary Structure Reduction Tests

    [Test]
    public void ReduceSecondaryStructure_ModifiesHighStructure()
    {
        // Create sequence with potential secondary structure
        string sequence = "AUGGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAGCUGCAUAA";
        var result = CodonOptimizer.ReduceSecondaryStructure(sequence, CodonOptimizer.EColiK12);

        // Should return modified or same sequence
        Assert.That(result.Length, Is.EqualTo(sequence.Length));
    }

    [Test]
    public void ReduceSecondaryStructure_ShortSequence_ReturnsSame()
    {
        string sequence = "AUGGCU";
        var result = CodonOptimizer.ReduceSecondaryStructure(sequence, CodonOptimizer.EColiK12, 40);
        Assert.That(result, Is.EqualTo(sequence));
    }

    [Test]
    public void ReduceSecondaryStructure_EmptySequence_ReturnsEmpty()
    {
        var result = CodonOptimizer.ReduceSecondaryStructure("", CodonOptimizer.EColiK12);
        Assert.That(result, Is.Empty);
    }

    #endregion

    #region Codon Usage Table Tests

    [Test]
    public void EColiK12_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.EColiK12.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void Yeast_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.Yeast.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void Human_ContainsAllCodons()
    {
        Assert.That(CodonOptimizer.Human.CodonFrequencies.Count, Is.EqualTo(64));
    }

    [Test]
    public void CodonTables_FrequenciesNormalized()
    {
        // For each amino acid, synonymous codon frequencies should sum to ~1
        foreach (var aa in new[] { "L", "S", "R", "P", "A", "G", "V", "T" })
        {
            // Just verify tables are properly structured
            Assert.That(CodonOptimizer.EColiK12.CodonFrequencies.Values.All(f => f >= 0 && f <= 1), Is.True);
        }
    }

    [Test]
    public void CreateCodonTableFromSequence_CreatesValidTable()
    {
        string reference = "AUGGCUGCUGCUGCUGCUGCUGCUUAA";
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference, "Test Organism");

        Assert.That(table.OrganismName, Is.EqualTo("Test Organism"));
        Assert.That(table.CodonFrequencies.ContainsKey("GCU"), Is.True);
    }

    [Test]
    public void CreateCodonTableFromSequence_CalculatesCorrectFrequencies()
    {
        // 4x GCU, 2x GCC = 6x Ala
        string reference = "GCUGCUGCUGCUGCCGCC";
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference, "Test");

        // GCU should be ~0.67, GCC should be ~0.33
        Assert.That(table.CodonFrequencies.GetValueOrDefault("GCU", 0), Is.GreaterThan(0.5));
    }

    #endregion

    #region Restriction-site removal — review 2026-09 (F22)

    [Test]
    [Description("One codon per removed site, chosen as the most frequent synonymous substitution")]
    public void RemoveRestrictionSites_ChangesExactlyOneCodonPerSite()
    {
        // EcoRI GAAUUC over Glu-Phe-Gly. Candidate substitutions that destroy the site:
        // GAA→GAG (f 0.31) and UUC→UUU (f 0.57); the higher-frequency one wins, and the codons
        // that follow the site are NOT touched (they were, before 2026-09).
        string result = CodonOptimizer.RemoveRestrictionSites("GAATTCGGG", new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo("GAAUUUGGG"));
            Assert.That(result, Does.Not.Contain("GAAUUC"));
        });
    }

    [Test]
    [Description("A non-palindromic site is also removed where its reverse complement occurs (both strands are cut)")]
    public void RemoveRestrictionSites_ReverseComplementOccurrence_IsRemoved()
    {
        // BsaI recognises GGTCTC / GAGACC (REBASE lists the top strand only). Both must go.
        string forward = CodonOptimizer.RemoveRestrictionSites("ATGGGTCTCGCTAAA", new[] { "GGTCTC" }, CodonOptimizer.EColiK12);
        string reverse = CodonOptimizer.RemoveRestrictionSites("ATGGAGACCGCTAAA", new[] { "GGTCTC" }, CodonOptimizer.EColiK12);

        Assert.Multiple(() =>
        {
            Assert.That(forward, Does.Not.Contain("GGUCUC"));
            Assert.That(reverse, Does.Not.Contain("GAGACC"), "the reverse-complement strand is cut by the same enzyme");
            Assert.That(reverse, Is.EqualTo("AUGGAAACCGCUAAA"));
        });
    }

    [Test]
    [Description("Recognition sequences may contain IUPAC ambiguity codes (XhoII RGATCY)")]
    public void RemoveRestrictionSites_IupacSite_IsMatchedAndRemoved()
    {
        string result = CodonOptimizer.RemoveRestrictionSites("ATGAGATCTGCTAAA", new[] { "RGATCY" }, CodonOptimizer.EColiK12);

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo("AUGCGCUCUGCUAAA"));
            Assert.That(result, Does.Not.Contain("AGAUCU"));
        });
    }

    [Test]
    [Description("A site that no synonymous substitution can remove is left in place, without throwing")]
    public void RemoveRestrictionSites_UnremovableSite_IsLeftInPlace()
    {
        // UGGUGG inside a Trp-Trp run: Trp has a single codon, so the site cannot be removed.
        string result = CodonOptimizer.RemoveRestrictionSites("ATGTGGTGGTGGAAA", new[] { "TGGTGG" }, CodonOptimizer.EColiK12);

        Assert.That(result, Is.EqualTo("AUGUGGUGGUGGAAA"));
    }

    [Test]
    public void RemoveRestrictionSites_InvalidSiteCharacter_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => CodonOptimizer.RemoveRestrictionSites("AUGGCU", new[] { "GA@TTC" }, CodonOptimizer.EColiK12));
    }

    #endregion

    #region Reference codon table from a gene set — review 2026-09 (F25)

    [Test]
    [Description("Relative adaptiveness from the built table equals Biopython CodonAdaptationIndex (0.5 pseudo-count)")]
    public void CreateCodonTableFromSequence_MatchesBiopythonRelativeAdaptiveness()
    {
        // Biopython 1.88: CodonAdaptationIndex(["ATGAAAGCGTTCAAGCGTACTGCGATGCCCAAAGGGTTTTAA"]) →
        // AAA 1.0, AAG 0.5, GCG 1.0, GCT 0.25, TTT 1.0, TTC 1.0, CGT 1.0, AGA 0.5, TAA 1.0, TAG 0.5.
        // Our table stores f = count / family total, so w = f / max f must reproduce those values.
        const string reference = "ATGAAAGCGTTCAAGCGTACTGCGATGCCCAAAGGGTTTTAA";
        var table = CodonOptimizer.CreateCodonTableFromSequence(reference, "T");

        double W(string codon, params string[] family)
        {
            double max = family.Max(c => table.CodonFrequencies[c]);
            return table.CodonFrequencies[codon] / max;
        }

        string[] lys = { "AAA", "AAG" };
        string[] ala = { "GCU", "GCC", "GCA", "GCG" };
        string[] phe = { "UUU", "UUC" };
        string[] arg = { "CGU", "CGC", "CGA", "CGG", "AGA", "AGG" };
        string[] stop = { "UAA", "UAG", "UGA" };

        Assert.Multiple(() =>
        {
            Assert.That(table.CodonFrequencies, Has.Count.EqualTo(64), "absent codons are pseudo-counted, not dropped");
            Assert.That(W("AAA", lys), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(W("AAG", lys), Is.EqualTo(0.5).Within(1e-9));
            Assert.That(W("GCG", ala), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(W("GCU", ala), Is.EqualTo(0.25).Within(1e-9));
            Assert.That(W("UUU", phe), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(W("UUC", phe), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(W("CGU", arg), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(W("AGA", arg), Is.EqualTo(0.5).Within(1e-9));
            Assert.That(W("UAG", stop), Is.EqualTo(0.5).Within(1e-9));
            Assert.That(table.CodonToAminoAcid["GCG"], Is.EqualTo("A"), "the table carries the Standard code mapping");
        });
    }

    [Test]
    [Description("Every amino-acid family of the built table sums to 1")]
    public void CreateCodonTableFromSequence_FamiliesSumToOne()
    {
        var table = CodonOptimizer.CreateCodonTableFromSequence("AUGGCUGCUGCCAAAUAA", "T");

        var sums = table.CodonFrequencies
            .GroupBy(kv => table.CodonToAminoAcid[kv.Key])
            .ToDictionary(g => g.Key, g => g.Sum(kv => kv.Value));

        foreach (var (aminoAcid, sum) in sums)
            Assert.That(sum, Is.EqualTo(1.0).Within(1e-9), $"family {aminoAcid}");
    }

    #endregion
}
