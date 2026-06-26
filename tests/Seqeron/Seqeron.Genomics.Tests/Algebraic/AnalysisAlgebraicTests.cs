using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Analysis area (genomic similarity).
///
/// Algebraic testing pins the reflexive maximum and symmetry of the k-mer
/// Jaccard similarity.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 179.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Analysis")]
public class AnalysisAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-SIMILARITY-001 — k-mer Jaccard similarity (Analysis), row 179.
    //
    // Model: similarity = |A∩B| / |A∪B| × 100 over the two sequences' k-mer sets —
    //        a symmetric Jaccard index whose self-comparison is the maximum (100).
    //   — docs/algorithms/Genomic_Analysis; GenomicAnalyzer.CalculateSimilarity.
    //
    // Laws (row 179): ID — sim(x, x) = 100 (max; Jaccard 1 on a sequence's own
    //                 k-mer set).  COMM — sim(a, b) = sim(b, a).
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property Similarity_Identity_SelfIsMax()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 6), seq =>
        {
            double sim = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(seq), new DnaSequence(seq));
            return (System.Math.Abs(sim - 100.0) < 1e-9).Label($"sim(x,x)={sim} for \"{seq}\"");
        });
    }

    [FsCheck.NUnit.Property]
    public Property Similarity_Commutative_Symmetric()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 6), DnaArbitrary(minLen: 6), (a, b) =>
        {
            double ab = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b));
            double ba = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(b), new DnaSequence(a));
            return (System.Math.Abs(ab - ba) < 1e-9).Label($"sim(a,b)={ab} != sim(b,a)={ba}");
        });
    }

    // An RNA alphabet sampler (the pseudoknot predictor reads T as U).
    private static Arbitrary<string> RnaArbitrary(int minLen, int maxLen) =>
        (from len in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf(len)
         select new string(chars)).ToArbitrary();

    // The designed canonical H-type pseudoknot from the RNA-PKPREDICT-001 evidence (M1):
    // two crossing 4-bp helices that beat the plain-MFE hairpin.
    private const string DesignedHType = "GGGGAACCCCAACCCCAAGGGG";

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-PKPREDICT-001 — Canonical H-type pseudoknot prediction (Analysis), row 236.
    //
    // Model: pknotsRG canonical simple-recursive pseudoknot (Reeder & Giegerich 2004) —
    //        a single H-type fold accepted only if it strictly beats the plain MFE.
    //   — RnaSecondaryStructure.PredictStructurePseudoknot; tests/TestSpecs/RNA-PKPREDICT-001.md.
    //
    // Laws (row 236): ID — no pairable bases → empty structure (the neutral, all-unpaired
    //                 input maps to no base pairs, ΔG = 0, no pseudoknot).
    //                 IDEMP — prediction is a pure, deterministic function: f(x) = f(x).
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Pseudoknot_Identity_NoPairableBasesIsEmptyStructure()
    {
        // A homopolymer has no Watson–Crick/GU pair anywhere, so the neutral input folds to nothing.
        var s = RnaSecondaryStructure.PredictStructurePseudoknot(new string('A', 24));
        s.BasePairs.Should().BeEmpty();
        s.HasPseudoknot.Should().BeFalse();
        s.FreeEnergy.Should().Be(0.0);
        s.DotBracket.Should().Be(new string('.', 24));
    }

    [Test]
    public void Pseudoknot_Identity_EmptyInputIsEmptyStructure()
    {
        var s = RnaSecondaryStructure.PredictStructurePseudoknot(string.Empty);
        s.BasePairs.Should().BeEmpty();
        s.HasPseudoknot.Should().BeFalse();
        s.FreeEnergy.Should().Be(0.0);
    }

    [FsCheck.NUnit.Property]
    public Property Pseudoknot_Idempotent_Deterministic()
    {
        return Prop.ForAll(RnaArbitrary(minLen: 11, maxLen: 30), seq =>
        {
            var a = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            var b = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            bool same = a.DotBracket == b.DotBracket
                        && a.FreeEnergy == b.FreeEnergy
                        && a.HasPseudoknot == b.HasPseudoknot
                        && a.BasePairs.SequenceEqual(b.BasePairs);
            return same.Label($"non-deterministic prediction for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-PKRECURSIVE-001 — Recursive pseudoknot folding (Analysis), row 237.
    //
    // Model: the recursive pknotsRG decomposition F(i,j) chains pseudoknot-free blocks
    //        and H-type knots whose loops fold by F again. For an interval that holds a
    //        SINGLE H-type knot it must collapse to the non-recursive single-knot result.
    //   — RnaSecondaryStructure.PredictStructurePseudoknotRecursive; RNA-PKRECURSIVE-001.md.
    //
    // Laws (row 237): ID — a single H-type input → recursive result equals the
    //                 non-recursive PredictStructurePseudoknot.
    //                 IDEMP — recursive prediction is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void RecursivePseudoknot_Identity_SingleHTypeEqualsNonRecursive()
    {
        var nonRecursive = RnaSecondaryStructure.PredictStructurePseudoknot(DesignedHType);
        var recursive = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(DesignedHType);

        // The single canonical knot must be recovered by both routes.
        nonRecursive.HasPseudoknot.Should().BeTrue();
        recursive.HasPseudoknot.Should().BeTrue();
        recursive.DotBracket.Should().Be(nonRecursive.DotBracket);
        recursive.FreeEnergy.Should().Be(nonRecursive.FreeEnergy);
        recursive.BasePairs.Should().Equal(nonRecursive.BasePairs);
    }

    [FsCheck.NUnit.Property]
    public Property RecursivePseudoknot_Idempotent_Deterministic()
    {
        return Prop.ForAll(RnaArbitrary(minLen: 11, maxLen: 30), seq =>
        {
            var a = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq);
            var b = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq);
            bool same = a.DotBracket == b.DotBracket
                        && a.FreeEnergy == b.FreeEnergy
                        && a.HasPseudoknot == b.HasPseudoknot
                        && a.BasePairs.SequenceEqual(b.BasePairs);
            return same.Label($"non-deterministic recursive prediction for \"{seq}\"");
        });
    }
}
