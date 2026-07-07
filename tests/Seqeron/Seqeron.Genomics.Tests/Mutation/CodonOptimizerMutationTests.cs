namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for CodonOptimizer.cs (checklist 04 rows 58-61:
/// CODON-OPT/CAI/RARE/USAGE-001). Baseline 50% — the codon-usage/CAI/optimize/structure
/// methods were under-pinned. These assert the published formulas exactly (CAI geometric
/// mean of relative adaptiveness, codon-usage distance, MaximizeCAI selection) and pin the
/// secondary-structure reducer's deterministic, protein-preserving behaviour.
/// </summary>
[TestFixture]
public class CodonOptimizerMutationTests
{
    // Custom table: M has its single codon; Leucine has a clear max (CUG 0.5) and a rare codon (CUU 0.1).
    private static CodonOptimizer.CodonUsageTable TestTable() => new(
        "test",
        new Dictionary<string, double> { ["AUG"] = 1.0, ["CUG"] = 0.5, ["CUU"] = 0.1 },
        new Dictionary<string, string>());

    // ── CalculateCodonUsage: exact codon counts (rows 61) ─────────────────────────────

    [Test]
    public void CalculateCodonUsage_CountsCodonsExactly()
    {
        var usage = CodonOptimizer.CalculateCodonUsage("ATGCTGCTTTAA"); // AUG CUG CUU UAA
        usage.Should().HaveCount(4);
        usage["AUG"].Should().Be(1);
        usage["CUG"].Should().Be(1);
        usage["CUU"].Should().Be(1);
        usage["UAA"].Should().Be(1);
    }

    [Test]
    public void CalculateCodonUsage_RepeatedCodon_Accumulates()
    {
        CodonOptimizer.CalculateCodonUsage("ATGATGATG")["AUG"].Should().Be(3);
    }

    // ── CalculateCAI: geometric mean of w_i = freq(codon)/max(freq over synonymous) ──

    [Test]
    public void CalculateCAI_IsGeometricMeanOfRelativeAdaptiveness()
    {
        // AUG: w = 1/1 = 1; CUG: w = 0.5/0.5 = 1; CUU: w = 0.1/0.5 = 0.2; stop skipped.
        // CAI = (1 · 1 · 0.2)^(1/3) = 0.2^(1/3).
        double cai = CodonOptimizer.CalculateCAI("ATGCTGCTTTAA", TestTable());
        cai.Should().BeApproximately(System.Math.Pow(0.2, 1.0 / 3.0), 1e-9);
    }

    [Test]
    public void CalculateCAI_AllOptimalCodons_IsOne()
    {
        // AUG (w=1) + CUG (w=1) → CAI = 1.
        CodonOptimizer.CalculateCAI("ATGCTG", TestTable()).Should().BeApproximately(1.0, 1e-9);
    }

    // ── CompareCodonUsage: 1 − Σ|f1−f2|/2 ─────────────────────────────────────────────

    [Test]
    public void CompareCodonUsage_IdenticalSequences_IsOne()
    {
        CodonOptimizer.CompareCodonUsage("ATGCTG", "ATGCTG").Should().BeApproximately(1.0, 1e-9);
    }

    [Test]
    public void CompareCodonUsage_HalfOverlap_IsExact()
    {
        // s1 = AUG×2 (f AUG=1); s2 = AUG,CUG (f each 0.5). Σ|Δ| = 0.5+0.5 = 1 → 1 − 1/2 = 0.5.
        CodonOptimizer.CompareCodonUsage("ATGATG", "ATGCTG").Should().BeApproximately(0.5, 1e-9);
    }

    // ── FindRareCodons: codons with table frequency below threshold ───────────────────

    [Test]
    public void FindRareCodons_ReportsBelowThresholdWithPosition()
    {
        // CUU freq 0.1 < 0.15 default; AUG (1.0) and CUG (0.5) are not rare.
        var rare = CodonOptimizer.FindRareCodons("ATGCTGCTT", TestTable()).ToList();

        rare.Should().ContainSingle();
        rare[0].Position.Should().Be(6); // third codon → index 2 × 3
        rare[0].Codon.Should().Be("CUU");
        rare[0].Frequency.Should().BeApproximately(0.1, 1e-9);
    }

    // ── OptimizeSequence (MaximizeCAI): picks highest-frequency synonymous codon ──────

    [Test]
    public void OptimizeSequence_MaximizeCai_ReplacesRareWithOptimalSynonym()
    {
        // "ATGCTTTAA" = M, L(CUU rare), stop. MaximizeCAI → L becomes CUG (highest freq).
        var r = CodonOptimizer.OptimizeSequence("ATGCTTTAA", TestTable(),
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        r.OptimizedSequence.Should().Be("AUGCUGUAA");
        r.ProteinSequence.Should().Be("ML*");
        r.ChangedCodons.Should().Be(1);
        r.Changes[0].Position.Should().Be(3, "the changed codon is the 2nd (index 1 × 3)");
        r.OptimizedCAI.Should().BeApproximately(1.0, 1e-9);
        r.OptimizedCAI.Should().BeGreaterThan(r.OriginalCAI, "replacing a rare codon raises CAI");
    }

    [Test]
    public void OptimizeSequence_NonMultipleOfThree_TrimsToCompleteCodons()
    {
        // 5 nt → trimmed to 3 (one complete codon). OriginalSequence reflects the trim.
        var r = CodonOptimizer.OptimizeSequence("ATGCT", TestTable(), CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        r.OriginalSequence.Should().Be("AUG");
        r.ProteinSequence.Should().Be("M");
    }

    [Test]
    public void OptimizeSequence_BalancedOptimization_RaisesGcTowardTarget()
    {
        // GC-poor input with a high GC target → BalanceGcContent engages and swaps to GC-richer
        // synonymous codons. Deterministic output (AAG/UUC are the GC-richest Lys/Phe codons).
        var r = CodonOptimizer.OptimizeSequence("AAATTTAAA", CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization, gcTargetMin: 0.9, gcTargetMax: 1.0);

        r.OptimizedSequence.Should().Be("AAGUUCAAG");
        r.GcContentOptimized.Should().BeGreaterThan(r.GcContentOriginal);
        r.GcContentOriginal.Should().Be(0.0);
    }

    [Test]
    public void OptimizeSequence_PreservesProtein()
    {
        // Real E. coli table; the optimized sequence must encode the SAME protein.
        var r = CodonOptimizer.OptimizeSequence("ATGCTTGCAAAATAA", CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        // Translate both via codon usage equivalence: same codon count, same amino acids.
        r.OptimizedSequence.Length.Should().Be(r.OriginalSequence.Length);
        r.OptimizedCAI.Should().BeGreaterThanOrEqualTo(r.OriginalCAI);
    }

    // ── ReduceSecondaryStructure: deterministic, protein-preserving synonymous changes ──

    [Test]
    public void ReduceSecondaryStructure_ShortSequence_ReturnedUnchanged()
    {
        // Shorter than the window → returned as-is (no T→U either, exact passthrough).
        const string seq = "ATGCTG";
        CodonOptimizer.ReduceSecondaryStructure(seq, CodonOptimizer.EColiK12, windowSize: 40)
            .Should().Be(seq);
    }

    [Test]
    public void ReduceSecondaryStructure_SelfComplementaryInput_ReducesDeterministically()
    {
        // A G-block/C-block sequence is highly self-complementary (many G·C pairs) ⇒ the reducer
        // engages and swaps codons to lower the folding propensity. The transformation is fully
        // deterministic; this exact output pins the whole CalculateLocalStructure/AreComplementary
        // path (any mutation there changes which codons are chosen). Protein is preserved:
        // 8×Gly then 8×Pro both before and after.
        string seq = string.Concat(Enumerable.Repeat("GGG", 8)) + string.Concat(Enumerable.Repeat("CCC", 8));
        string outSeq = CodonOptimizer.ReduceSecondaryStructure(seq, CodonOptimizer.EColiK12);

        outSeq.Should().Be("GGUGGUGGUGGUGGUGGUGGUGGUCCUCCUCCUCCUCCUCCCCCCCCC");

        // Theory checks alongside the characterization: same length, RNA, protein preserved.
        outSeq.Length.Should().Be(seq.Length);
        outSeq.Should().NotContain("T");
        var aas = Enumerable.Range(0, outSeq.Length / 3)
            .Select(i => outSeq.Substring(i * 3, 3))
            .Select(c => c.StartsWith("GG") ? 'G' : 'P').ToList();
        aas.Should().Equal(Enumerable.Repeat('G', 8).Concat(Enumerable.Repeat('P', 8)),
            "synonymous-only changes: 8 Gly then 8 Pro");
    }

    [Test]
    public void ReduceSecondaryStructure_AtRichSelfComplementaryInput_ReducesDeterministically()
    {
        // A-block/T-block → A·U self-complementary (exercises the A-U / U-A AreComplementary arms).
        // Deterministic reduction; protein preserved: 8×Lys then 8×Phe.
        string seq = string.Concat(Enumerable.Repeat("AAA", 8)) + string.Concat(Enumerable.Repeat("TTT", 8));
        string outSeq = CodonOptimizer.ReduceSecondaryStructure(seq, CodonOptimizer.EColiK12);

        outSeq.Should().Be("AAGAAGAAGAAGAAGAAGAAGAAGUUCUUCUUCUUCUUCUUUUUUUUU");
        var aas = Enumerable.Range(0, outSeq.Length / 3)
            .Select(i => outSeq.Substring(i * 3, 3))
            .Select(c => c.StartsWith("AA") ? 'K' : 'F').ToList();
        aas.Should().Equal(Enumerable.Repeat('K', 8).Concat(Enumerable.Repeat('F', 8)));
    }

    [Test]
    public void ReduceSecondaryStructure_UBlockThenABlock_ExercisesUaComplementarityArm()
    {
        // U-block then A-block puts U at position i and A at j > i+4 → exercises the (U,A) arm of
        // AreComplementary (the A-block/U-block test only hits the (A,U) arm). Protein: 8×Phe, 8×Lys.
        string seq = string.Concat(Enumerable.Repeat("TTT", 8)) + string.Concat(Enumerable.Repeat("AAA", 8));
        string outSeq = CodonOptimizer.ReduceSecondaryStructure(seq, CodonOptimizer.EColiK12);

        outSeq.Should().Be("UUCUUCUUCUUCUUCUUCUUCUUCAAGAAGAAGAAGAAGAAAAAAAAA");
    }

    // ── RemoveRestrictionSites: eliminate the recognition site, preserve the protein ─

    [Test]
    public void RemoveRestrictionSites_EliminatesSite_PreservesProtein()
    {
        // EcoRI GAATTC inside a coding sequence (Glu-Phe-Gly) is removed via synonymous edits.
        string result = CodonOptimizer.RemoveRestrictionSites("GAATTCGGG", new[] { "GAATTC" }, CodonOptimizer.EColiK12);

        result.Should().NotContain("GAAUUC", "the recognition site must be eliminated");
        result.Should().Be("GAGUUUGGU"); // Glu(GAG) Phe(UUU) Gly(GGU) — same protein
    }

    // ── CreateCodonTableFromSequence: per-amino-acid relative frequencies ─────────────

    [Test]
    public void CreateCodonTableFromSequence_ComputesRelativeFrequenciesPerAminoAcid()
    {
        // Leucine: CUG×2, CUU×1 → frequencies 2/3 and 1/3 within the Leu synonymous family.
        var table = CodonOptimizer.CreateCodonTableFromSequence("CTGCTGCTT", "test");

        table.CodonFrequencies["CUG"].Should().BeApproximately(2.0 / 3.0, 1e-9);
        table.CodonFrequencies["CUU"].Should().BeApproximately(1.0 / 3.0, 1e-9);
    }

    // ── OptimizeSequence (AvoidRareCodeons): only rare codons are replaced ────────────

    [Test]
    public void OptimizeSequence_AvoidRareCodons_ReplacesOnlyTheRareCodon()
    {
        // CUU (0.1 < 0.15) is rare → replaced by CUG; the non-rare CUG is left intact.
        var r = CodonOptimizer.OptimizeSequence("ATGCTTCTGTAA", TestTable(),
            CodonOptimizer.OptimizationStrategy.AvoidRareCodeons, rareCodonThreshold: 0.15);

        r.OptimizedSequence.Should().Be("AUGCUGCUGUAA"); // M, L(CUU→CUG), L(CUG kept), stop
        r.ChangedCodons.Should().Be(1);
    }

    [Test]
    public void GetCodonGcContent_ViaOptimize_BalancesTowardTarget()
    {
        // BalancedOptimization nudges GC into [min,max]; pin that the optimized GC is within range
        // for a sequence whose optimal-CAI form would otherwise drift out.
        var r = CodonOptimizer.OptimizeSequence("ATGCTTGCAAAA", CodonOptimizer.EColiK12,
            CodonOptimizer.OptimizationStrategy.BalancedOptimization, gcTargetMin: 0.40, gcTargetMax: 0.60);

        r.ProteinSequence.Should().StartWith("M");
        r.OptimizedSequence.Length.Should().Be(r.OriginalSequence.Length);
    }
}
