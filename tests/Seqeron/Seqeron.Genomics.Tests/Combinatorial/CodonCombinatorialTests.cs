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
}
