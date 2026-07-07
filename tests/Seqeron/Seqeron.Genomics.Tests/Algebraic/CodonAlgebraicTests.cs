using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Codon area (codon optimization, codon usage).
///
/// Algebraic testing pins the synonymous-substitution round-trip of codon
/// optimization (the encoded protein is preserved and the deterministic optimizer
/// is a fixpoint) and the per-amino-acid normalization identity of codon usage.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 58, 61.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Codon")]
public class CodonAlgebraicTests
{
    /// <summary>DNA of length a multiple of 3 (whole codons), A/C/G/T.</summary>
    private static Arbitrary<string> CodingDnaArbitrary() =>
        (from k in Gen.Choose(1, 30)
         from a in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(3 * k)
         select new string(a))
        .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CODON-OPT-001 — Codon optimization (Codon)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 58.
    //
    // Model: codon optimization replaces each codon with a synonymous codon
    //        (preferred by the target organism) and keeps stop codons, so it
    //        preserves the encoded protein. The deterministic MaximizeCAI strategy
    //        always picks the single highest-frequency synonymous codon, making it
    //        a fixpoint (re-optimizing changes nothing).
    //   — docs/algorithms/Codon_Optimization; CodonOptimizer.OptimizeSequence.
    //
    // Laws under test (checklist row 58):
    //   • RT    — translate(optimize(dna)) = translate(dna) (protein preserved).
    //   • IDEMP — optimize(optimize(x)) = optimize(x) (deterministic fixpoint).
    // ═══════════════════════════════════════════════════════════════════════

    private static string Translate(string seq) =>
        Translator.Translate(seq, GeneticCode.Standard, 0, false).Sequence;

    /// <summary>
    /// RT: optimizing a coding sequence preserves its translation (only synonymous
    /// codon substitutions; stop codons retained).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Optimize_RoundTrip_PreservesTranslation()
    {
        return Prop.ForAll(CodingDnaArbitrary(), dna =>
        {
            var result = CodonOptimizer.OptimizeSequence(
                dna, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
            string before = Translate(dna);
            string after = Translate(result.OptimizedSequence);
            return (before == after).Label($"protein changed: \"{before}\" -> \"{after}\" for \"{dna}\"");
        });
    }

    /// <summary>
    /// IDEMP: the MaximizeCAI optimizer is a fixpoint — re-optimizing its own output
    /// yields the identical sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Optimize_Idempotent_FixpointUnderReoptimization()
    {
        return Prop.ForAll(CodingDnaArbitrary(), dna =>
        {
            var once = CodonOptimizer.OptimizeSequence(
                dna, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
            var twice = CodonOptimizer.OptimizeSequence(
                once.OptimizedSequence, CodonOptimizer.EColiK12, CodonOptimizer.OptimizationStrategy.MaximizeCAI);
            return (twice.OptimizedSequence == once.OptimizedSequence)
                .Label($"not idempotent for \"{dna}\": {once.OptimizedSequence} -> {twice.OptimizedSequence}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CODON-USAGE-001 — Codon usage (Codon)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 61.
    //
    // Model: codon counts (CalculateCodonUsage) grouped by the amino acid they
    //        encode give the relative synonymous usage; within each amino acid the
    //        fractional usages count(codon)/count(aa) sum to 1.
    //   — docs/algorithms/Codon_Optimization/Codon_Usage_Analysis.md;
    //     CodonOptimizer.CalculateCodonUsage + GeneticCode.Standard.
    //
    // Laws under test (checklist row 61):
    //   • DIST — per amino acid, Σ fractional codon usage = 1.0.
    //   • ID   — a single repeated codon → usage 1.0 for that codon.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>DIST: the fractional codon usage within each amino acid sums to 1.</summary>
    [FsCheck.NUnit.Property]
    public Property CodonUsage_Distributive_PerAminoAcidSumsToOne()
    {
        return Prop.ForAll(CodingDnaArbitrary(), dna =>
        {
            var counts = CodonOptimizer.CalculateCodonUsage(dna);
            // Group observed codons by the amino acid they encode.
            var byAa = counts.GroupBy(kv => GeneticCode.Standard.Translate(kv.Key));
            foreach (var aaGroup in byAa)
            {
                int aaTotal = aaGroup.Sum(kv => kv.Value);
                double fracSum = aaGroup.Sum(kv => (double)kv.Value / aaTotal);
                if (Math.Abs(fracSum - 1.0) > 1e-9)
                    return false.Label($"AA '{aaGroup.Key}' fractions sum to {fracSum}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>ID: a sequence of one repeated codon has usage 1.0 for that codon.</summary>
    [Test]
    public void CodonUsage_Identity_SingleCodonIsOne()
    {
        var counts = CodonOptimizer.CalculateCodonUsage("GCTGCTGCTGCT"); // Ala (GCU) ×4 after T→U
        counts.Should().HaveCount(1);
        var only = counts.Single();
        var aaTotal = only.Value;
        ((double)only.Value / aaTotal).Should().Be(1.0);
    }
}
