using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for RNA secondary structure prediction.
/// Verifies structural, stem-loop, and energy invariants using FsCheck.
///
/// Test Units: RNA-STRUCT-001, RNA-STEMLOOP-001, RNA-ENERGY-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class RnaStructureProperties
{
    #region Generators

    /// <summary>
    /// Generates random RNA sequences from valid bases (A, C, G, U).
    /// </summary>
    private static Arbitrary<string> RnaArbitrary(int minLen = 10) =>
        Gen.Elements('A', 'C', 'G', 'U')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Generates RNA sequences with known stem-loop potential (palindromic arms).
    /// Pattern: complementary_arm + loop(≥3nt) + reverse_arm.
    /// </summary>
    private static Arbitrary<string> StemLoopRnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'U')
            .ArrayOf()
            .Where(a => a.Length >= 4 && a.Length <= 8)
            .SelectMany(arm =>
                Gen.Elements('A', 'C', 'G', 'U').ArrayOf()
                    .Where(loop => loop.Length >= 3 && loop.Length <= 8)
                    .Select(loop =>
                    {
                        var comp = arm.Reverse().Select(RnaComplement).ToArray();
                        return new string(arm) + new string(loop) + new string(comp);
                    }))
            .ToArbitrary();

    /// <summary>
    /// GC-rich RNA: forms more stable structures (lower ΔG).
    /// Evidence: GC pairs have 3 hydrogen bonds vs 2 for AU.
    /// </summary>
    private static string GcRichRna =>
        "GGGCCCGGGAAACCCGGGCCC";

    /// <summary>
    /// AU-rich RNA: forms less stable structures (higher ΔG).
    /// </summary>
    private static string AuRichRna =>
        "AAAUUUAAAGGGAAAUUUAAA";

    private static char RnaComplement(char b) => b switch
    {
        'A' => 'U',
        'U' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => 'A'
    };

    #endregion

    #region RNA-STRUCT-001: R: pairs ≤ len/2; P: no crossing pairs; P: complementary; D: deterministic

    /// <summary>
    /// INV-1: Number of base pairs ≤ len/2.
    /// Evidence: Each nucleotide can participate in at most one base pair.
    /// Source: Nussinov &amp; Jacobson (1980) — maximum matching constraint.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictStructure_BasePairCount_AtMostHalfLength()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            var structure = RnaSecondaryStructure.PredictStructure(seq);
            return (structure.BasePairs.Count <= seq.Length / 2)
                .Label($"BasePairs={structure.BasePairs.Count} must be ≤ {seq.Length / 2}");
        });
    }

    /// <summary>
    /// INV-2: No crossing pairs in predicted structure (planar graph / Nussinov constraint).
    /// Evidence: Standard RNA secondary structure excludes pseudoknots in the
    /// base prediction. For any two pairs (i,j) and (k,l) with i &lt; k,
    /// either j &lt; k (nested or disjoint) or l &lt; j (nested inside).
    /// Source: Nussinov et al. (1978) — dynamic programming forbids crossing pairs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictStructure_NoCrossingPairs()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            var structure = RnaSecondaryStructure.PredictStructure(seq);
            var pairs = structure.BasePairs;
            for (int a = 0; a < pairs.Count; a++)
            {
                int i = Math.Min(pairs[a].Position1, pairs[a].Position2);
                int j = Math.Max(pairs[a].Position1, pairs[a].Position2);
                for (int b = a + 1; b < pairs.Count; b++)
                {
                    int k = Math.Min(pairs[b].Position1, pairs[b].Position2);
                    int l = Math.Max(pairs[b].Position1, pairs[b].Position2);
                    // Crossing: i < k < j < l
                    if (i < k && k < j && j < l)
                        return false.Label($"Crossing pairs: ({i},{j}) and ({k},{l})");
                }
            }
            return true.Label("No crossing pairs");
        });
    }

    /// <summary>
    /// INV-3: All predicted base pairs are complementary (A-U, G-C, or G-U wobble).
    /// Evidence: RNA folding only allows Watson-Crick and wobble pairs.
    /// Source: Turner &amp; Mathews (2010) — NNDB nearest-neighbor model.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictStructure_AllBasePairs_AreComplementary()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            var structure = RnaSecondaryStructure.PredictStructure(seq);
            return structure.BasePairs.All(bp =>
                RnaSecondaryStructure.CanPair(bp.Base1, bp.Base2))
                .Label("All base pairs must be complementary (WC or wobble)");
        });
    }

    /// <summary>
    /// INV-4: Dot-bracket length equals sequence length.
    /// Evidence: Each position maps to exactly one character in dot-bracket notation.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictStructure_DotBracketLength_EqualsSequenceLength()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            var structure = RnaSecondaryStructure.PredictStructure(seq);
            return (structure.DotBracket.Length == seq.Length)
                .Label($"DotBracket length={structure.DotBracket.Length}, seqLen={seq.Length}");
        });
    }

    /// <summary>
    /// INV-5: PredictStructure is deterministic — same input yields identical output.
    /// Evidence: Zuker-style DP is a pure function with no randomness.
    /// Source: Zuker &amp; Stiegler (1981) — deterministic energy minimization.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictStructure_IsDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            var s1 = RnaSecondaryStructure.PredictStructure(seq);
            var s2 = RnaSecondaryStructure.PredictStructure(seq);
            bool same = s1.DotBracket == s2.DotBracket &&
                        Math.Abs(s1.MinimumFreeEnergy - s2.MinimumFreeEnergy) < 1e-10 &&
                        s1.BasePairs.Count == s2.BasePairs.Count;
            return same.Label("PredictStructure must be deterministic");
        });
    }

    #endregion

    #region RNA-STEMLOOP-001: R: stem len > 0; P: loop len ≥ minLoop; P: stem arms complementary; D: deterministic

    /// <summary>
    /// INV-6: Every found stem-loop has stem length > 0.
    /// Evidence: A stem-loop requires at least one base pair in the stem.
    /// Source: Wikipedia Stem-loop — hairpin requires a double-stranded stem region.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindStemLoops_StemLength_GreaterThanZero()
    {
        return Prop.ForAll(StemLoopRnaArbitrary(), seq =>
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq, minStemLength: 3).ToList();
            return stemLoops.All(sl => sl.Stem.Length > 0)
                .Label("All stem-loops must have stem length > 0");
        });
    }

    /// <summary>
    /// INV-7: Loop size ≥ minLoop (default 3).
    /// Evidence: Loops fewer than 3 bases are sterically impossible.
    /// Source: NNDB Turner 2004 — "The nearest neighbor rules prohibit hairpin
    /// loops with fewer than 3 nucleotides."
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindStemLoops_LoopSize_AtLeastMinLoop()
    {
        return Prop.ForAll(StemLoopRnaArbitrary(), seq =>
        {
            const int minLoop = 3;
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq, minLoopSize: minLoop).ToList();
            return stemLoops.All(sl => sl.Loop.Size >= minLoop)
                .Label($"All loops must have size ≥ {minLoop}");
        });
    }

    /// <summary>
    /// INV-8: Stem arm base pairs are all complementary (A-U, G-C, or G-U wobble).
    /// Evidence: Stem formation requires hydrogen bonding between paired bases.
    /// Source: Turner &amp; Mathews (2010) — allowed pairs in nearest-neighbor model.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindStemLoops_StemArms_AreComplementary()
    {
        return Prop.ForAll(StemLoopRnaArbitrary(), seq =>
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq, minStemLength: 3).ToList();
            return stemLoops.All(sl =>
                sl.Stem.BasePairs.All(bp => RnaSecondaryStructure.CanPair(bp.Base1, bp.Base2)))
                .Label("All stem base pairs must be complementary");
        });
    }

    /// <summary>
    /// INV-9: FindStemLoops is deterministic.
    /// Evidence: Deterministic scanning algorithm with no randomness.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindStemLoops_IsDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(12), seq =>
        {
            var sl1 = RnaSecondaryStructure.FindStemLoops(seq).ToList();
            var sl2 = RnaSecondaryStructure.FindStemLoops(seq).ToList();
            bool same = sl1.Count == sl2.Count &&
                        sl1.Zip(sl2).All(p => p.First.Start == p.Second.Start &&
                                               p.First.End == p.Second.End &&
                                               Math.Abs(p.First.TotalFreeEnergy - p.Second.TotalFreeEnergy) < 1e-10);
            return same.Label("FindStemLoops must be deterministic");
        });
    }

    /// <summary>
    /// INV-10: Stem-loop positions are within sequence bounds.
    /// Evidence: Start and End indices reference positions within the input string.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindStemLoops_Positions_WithinBounds()
    {
        return Prop.ForAll(StemLoopRnaArbitrary(), seq =>
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq).ToList();
            return stemLoops.All(sl => sl.Start >= 0 && sl.End < seq.Length)
                .Label("Stem-loop positions must be within [0, seqLen)");
        });
    }

    /// <summary>
    /// Poly-A has no stem-loops — single base type cannot form pairs.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PolyA_HasNoStemLoops()
    {
        string rna = "AAAAAAAAAAAAAAAAAAA";
        var stemLoops = RnaSecondaryStructure.FindStemLoops(rna).ToList();
        Assert.That(stemLoops, Is.Empty, "Poly-A should have no stem-loops");
    }

    #endregion

    #region RNA-ENERGY-001: R: ΔG ≤ 0 for stable structures; M: more GC → lower energy; D: deterministic

    /// <summary>
    /// INV-11: MFE ≤ 0 for sequences capable of forming structure.
    /// Evidence: Folded structures are energetically favorable (negative ΔG) or neutral.
    /// Source: Turner 2004 — all stacking energies are negative; hairpin initiation is positive
    /// but is overcome when sufficient stem stacking exists.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MinimumFreeEnergy_ForStructuredRna_IsNonPositive()
    {
        // A sequence with a clear stem-loop: GGGAAACCC folds into (((...)))
        string rna = "GGGAAACCCUUUAAAGGGCCC";
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);
        Assert.That(mfe, Is.LessThanOrEqualTo(0.001),
            $"MFE={mfe} should be ≤ 0 for a structured RNA");
    }

    /// <summary>
    /// INV-12: GC-rich RNA has lower (more negative) MFE than AU-rich RNA.
    /// Evidence: GC base pairs contribute 3 hydrogen bonds and stronger stacking
    /// (e.g., GC/CG = -3.42 kcal/mol vs AU/UA = -1.10 kcal/mol in Turner 2004).
    /// Source: NNDB Turner 2004 stacking parameters.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MinimumFreeEnergy_GcRich_LowerThanAuRich()
    {
        double mfeGc = RnaSecondaryStructure.CalculateMinimumFreeEnergy(GcRichRna);
        double mfeAu = RnaSecondaryStructure.CalculateMinimumFreeEnergy(AuRichRna);

        Assert.That(mfeGc, Is.LessThan(mfeAu),
            $"GC-rich MFE ({mfeGc}) should be < AU-rich MFE ({mfeAu})");
    }

    /// <summary>
    /// INV-13: CalculateMinimumFreeEnergy is deterministic.
    /// Evidence: Zuker-style DP is a pure deterministic algorithm.
    /// Source: Zuker &amp; Stiegler (1981).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MinimumFreeEnergy_IsDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            double mfe1 = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            double mfe2 = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            return (Math.Abs(mfe1 - mfe2) < 1e-10)
                .Label($"MFE must be deterministic: {mfe1} vs {mfe2}");
        });
    }

    /// <summary>
    /// INV-14: MFE is finite for any valid RNA input.
    /// Evidence: DP algorithm produces bounded numerical results for valid inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MinimumFreeEnergy_IsFinite()
    {
        return Prop.ForAll(RnaArbitrary(10), seq =>
        {
            double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            return double.IsFinite(mfe).Label($"MFE={mfe} must be finite");
        });
    }

    /// <summary>
    /// INV-15: Stem-loop total free energy is finite.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property StemLoop_Energy_IsFinite()
    {
        return Prop.ForAll(StemLoopRnaArbitrary(), seq =>
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq).ToList();
            return stemLoops.All(sl => double.IsFinite(sl.TotalFreeEnergy))
                .Label("All stem-loop energies must be finite");
        });
    }

    #endregion
}
