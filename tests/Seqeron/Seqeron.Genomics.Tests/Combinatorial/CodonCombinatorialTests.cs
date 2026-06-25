namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Codon area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Codon")]
public class CodonCombinatorialTests
{
    // Curated sense codons (each has synonyms ⇒ real optimisation opportunities; no stop codons).
    private static readonly string[] SenseCodons = { "CTT", "GTT", "TCT", "CCT", "ACT", "GCT", "CGT", "GGT", "AAA", "GAA" };

    private static string CodingSequence(int nCodons) =>
        string.Concat(Enumerable.Range(0, nCodons).Select(i => SenseCodons[(i * 7 + 3) % SenseCodons.Length]));

    private static string Protein(string seq) =>
        Translator.Translate(seq.Replace('U', 'T'), GeneticCode.Standard).Sequence;

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CODON-OPT-001 — Codon optimization (Codon)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 58.
    // Spec: tests/TestSpecs/CODON-OPT-001.md (canonical CodonOptimizer.OptimizeSequence).
    // Dimensions: organism(3) × strategy(2: freq/balanced) × seqLen(3). Grid 3×2×3 = 18.
    //
    // Model (codon optimization): replace each codon with a synonymous codon favoured by the
    // target organism's usage table to raise the Codon Adaptation Index (CAI; Sharp & Li 1987)
    // WITHOUT changing the encoded protein. The "freq/random" strategy axis maps to the two
    // implemented strategies most relevant here: MaximizeCAI (frequency-greedy) and
    // BalancedOptimization (frequency + GC balancing).
    //
    // The combinatorial point: organism table, strategy and length interact, yet the cardinal
    // invariant — the optimized sequence encodes the identical protein — holds in every cell,
    // CAI stays in [0,1], and the frequency-greedy strategy never lowers CAI.
    // ═══════════════════════════════════════════════════════════════════════

    public enum OptStrategy { MaximizeCai, Balanced }

    private static CodonOptimizer.CodonUsageTable Table(int organism) => organism switch
    {
        0 => CodonOptimizer.EColiK12,
        1 => CodonOptimizer.Yeast,
        _ => CodonOptimizer.Human,
    };

    [Test, Combinatorial]
    public void CodonOpt_PreservesProtein_AndImprovesAdaptation(
        [Values(0, 1, 2)] int organism,
        [Values(OptStrategy.MaximizeCai, OptStrategy.Balanced)] OptStrategy strategy,
        [Values(10, 30, 60)] int nCodons)
    {
        string coding = CodingSequence(nCodons);
        var s = strategy == OptStrategy.MaximizeCai
            ? CodonOptimizer.OptimizationStrategy.MaximizeCAI
            : CodonOptimizer.OptimizationStrategy.BalancedOptimization;

        var r = CodonOptimizer.OptimizeSequence(coding, Table(organism), s);

        // Cardinal invariant: synonymous substitution preserves the protein.
        Protein(r.OptimizedSequence).Should().Be(Protein(coding), "optimization is synonymous");
        r.ProteinSequence.Should().Be(Protein(coding), "the reported protein matches the input");
        r.OptimizedSequence.Length.Should().Be(coding.Length, "codon count is preserved");

        r.OriginalCAI.Should().BeInRange(0.0, 1.0);
        r.OptimizedCAI.Should().BeInRange(0.0, 1.0);

        if (strategy == OptStrategy.MaximizeCai)
            r.OptimizedCAI.Should().BeGreaterThanOrEqualTo(r.OriginalCAI - 1e-9,
                "frequency-greedy optimization never lowers CAI");
    }

    /// <summary>
    /// Interaction witness: MaximizeCAI strictly improves a deliberately rare-codon sequence —
    /// the optimized CAI exceeds the original for the E. coli table.
    /// </summary>
    [Test]
    public void CodonOpt_MaximizeCai_RaisesRareCodonSequence()
    {
        // Repeated Leu codon CTA is rare in E. coli; optimization should pick a frequent synonym.
        string rare = string.Concat(Enumerable.Repeat("CTA", 20));
        var r = CodonOptimizer.OptimizeSequence(rare, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        r.OptimizedCAI.Should().BeGreaterThan(r.OriginalCAI, "a rare-codon gene gains CAI");
        Protein(r.OptimizedSequence).Should().Be(Protein(rare), "still all-leucine");
    }

    /// <summary>
    /// Interaction witness: re-optimising an already CAI-maximised sequence does not lower its
    /// CAI (optimisation is stable at the optimum).
    /// </summary>
    [Test]
    public void CodonOpt_MaximizeCai_IsStableAtOptimum()
    {
        string coding = CodingSequence(30);
        var first = CodonOptimizer.OptimizeSequence(coding, CodonOptimizer.Human, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
        var second = CodonOptimizer.OptimizeSequence(first.OptimizedSequence.Replace('U', 'T'), CodonOptimizer.Human,
            CodonOptimizer.OptimizationStrategy.MaximizeCAI);

        second.OptimizedCAI.Should().BeApproximately(first.OptimizedCAI, 1e-9, "the optimum is a fixed point");
        Protein(second.OptimizedSequence).Should().Be(Protein(coding));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CODON-RARE-001 — Rare codon detection (Codon)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 60.
    // Spec: tests/TestSpecs/CODON-RARE-001.md (canonical FindRareCodons).
    // Dimensions: threshold(3) × referenceTable(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model: a codon is "rare" for an organism when its usage frequency in that organism's table
    // is below a threshold (rare codons slow/stall translation). FindRareCodons reports every
    // such codon with its position and frequency.
    //
    // The combinatorial point: threshold, reference organism and length interact — the flagged
    // set is exactly the codons below threshold in the chosen table (verified against the table),
    // it grows monotonically with the threshold, and the same codon can be rare in one organism
    // but common in another.
    // ═══════════════════════════════════════════════════════════════════════

    private static List<string> RnaCodons(string coding)
    {
        string rna = coding.ToUpperInvariant().Replace('T', 'U');
        return Enumerable.Range(0, rna.Length / 3).Select(i => rna.Substring(i * 3, 3)).ToList();
    }

    [Test, Combinatorial]
    public void CodonRare_FlagsExactlyBelowThreshold(
        [Values(0.10, 0.20, 0.40)] double threshold,
        [Values(0, 1, 2)] int organism,
        [Values(10, 30, 60)] int nCodons)
    {
        string coding = CodingSequence(nCodons);
        var table = Table(organism);

        var rare = CodonOptimizer.FindRareCodons(coding, table, threshold).ToList();

        var expected = RnaCodons(coding)
            .Select((c, i) => (Pos: i * 3, Codon: c, Freq: table.CodonFrequencies.GetValueOrDefault(c, 0)))
            .Where(x => x.Freq < threshold)
            .ToList();

        rare.Should().HaveCount(expected.Count, "exactly the sub-threshold codons are flagged");
        for (int j = 0; j < expected.Count; j++)
        {
            rare[j].Position.Should().Be(expected[j].Pos);
            rare[j].Codon.Should().Be(expected[j].Codon);
            rare[j].Frequency.Should().BeApproximately(expected[j].Freq, 1e-12);
        }
        rare.Should().OnlyContain(r => r.Frequency < threshold);
    }

    /// <summary>
    /// Interaction witness: raising the threshold can only add rare codons — the flagged
    /// position set at a lower threshold is a subset of that at a higher one.
    /// </summary>
    [Test]
    public void CodonRare_Threshold_IsMonotone()
    {
        string coding = CodingSequence(40);
        var low = CodonOptimizer.FindRareCodons(coding, CodonOptimizer.EColiK12, 0.10).Select(r => r.Position).ToHashSet();
        var high = CodonOptimizer.FindRareCodons(coding, CodonOptimizer.EColiK12, 0.30).Select(r => r.Position).ToHashSet();
        low.Should().BeSubsetOf(high, "a higher threshold flags at least as many codons");
    }

    /// <summary>
    /// Interaction witness: rarity is organism-specific — a codon below threshold in one table
    /// can be at or above threshold in another, so the flagged sets differ by reference table.
    /// </summary>
    [Test]
    public void CodonRare_IsOrganismSpecific()
    {
        string coding = CodingSequence(60);
        var eColi = CodonOptimizer.FindRareCodons(coding, CodonOptimizer.EColiK12, 0.20).Select(r => r.Position).ToHashSet();
        var human = CodonOptimizer.FindRareCodons(coding, CodonOptimizer.Human, 0.20).Select(r => r.Position).ToHashSet();
        eColi.Should().NotBeEquivalentTo(human, "codon rarity depends on the organism's usage table");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CODON-ENC-001 — Effective Number of Codons (Wright 1990) (Codon)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 213.
    // Spec: tests/TestSpecs/CODON-ENC-001.md (canonical CodonUsageAnalyzer.CalculateEnc). ADVANCED §10.
    // Dimensions: geneticCode(3) × seqLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Wright 1990): the Effective Number of Codons Nc ∈ [20, 61] measures synonymous-codon
    // bias — 20 when each amino acid uses a single codon (maximal bias), 61 when all synonymous
    // codons are used equally (no bias).
    //
    // Axis mapping (documented — CalculateEnc fixes the standard code): geneticCode → the codon-bias
    // regime {Biased (one codon used), Mixed, Uniform (all sense codons cycled)}; seqLen → number of
    // codons. The combinatorial point: Nc stays within the Wright bound [20,61] at every cell, and a
    // biased sequence has a strictly lower Nc than an unbiased one (witness).
    // ═══════════════════════════════════════════════════════════════════════

    public enum CodonBias { Biased, Mixed, Uniform }

    private static readonly string[] EncSenseCodons = BuildEncSenseCodons();

    private static string[] BuildEncSenseCodons()
    {
        const string b = "ACGT";
        var stops = new HashSet<string> { "TAA", "TAG", "TGA" };
        var list = new List<string>();
        foreach (char x in b) foreach (char y in b) foreach (char z in b)
        {
            string c = $"{x}{y}{z}";
            if (!stops.Contains(c)) list.Add(c);
        }
        return list.ToArray();
    }

    private static string EncSequence(CodonBias bias, int nCodons)
    {
        string[] palette = bias switch
        {
            CodonBias.Biased => new[] { "GCT" },                 // one codon (Ala) ⇒ maximal bias
            CodonBias.Mixed => EncSenseCodons.Take(10).ToArray(),
            _ => EncSenseCodons,                                    // all 61 sense codons ⇒ minimal bias
        };
        var sb = new System.Text.StringBuilder(nCodons * 3);
        for (int i = 0; i < nCodons; i++) sb.Append(palette[i % palette.Length]);
        return sb.ToString();
    }

    [Test, Combinatorial]
    public void CodonEnc_WithinWrightBound_AcrossBiasAndLength(
        [Values(CodonBias.Biased, CodonBias.Mixed, CodonBias.Uniform)] CodonBias bias,
        [Values(60, 150, 300)] int nCodons)
    {
        double enc = CodonUsageAnalyzer.CalculateEnc(EncSequence(bias, nCodons));
        enc.Should().BeInRange(20.0, 61.0, "Nc lies in the Wright [20,61] bound");
    }

    /// <summary>
    /// Interaction witness — codon bias lowers Nc: a single-codon sequence is near the 20 floor,
    /// while an all-codons sequence approaches the 61 ceiling.
    /// </summary>
    [Test]
    public void CodonEnc_BiasLowersEffectiveNumber()
    {
        double biased = CodonUsageAnalyzer.CalculateEnc(EncSequence(CodonBias.Biased, 200));
        double uniform = CodonUsageAnalyzer.CalculateEnc(EncSequence(CodonBias.Uniform, 200));
        biased.Should().BeLessThan(uniform, "stronger codon bias means fewer effective codons");
        uniform.Should().BeGreaterThan(biased + 1.0, "an unbiased sequence uses many more effective codons");
    }
}
