using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for RNA secondary structure prediction.
/// Verifies structural, stem-loop, and energy invariants using FsCheck.
///
/// Test Units: RNA-STRUCT-001, RNA-STEMLOOP-001, RNA-ENERGY-001, RNA-DOTBRACKET-001, RNA-HAIRPIN-001, RNA-INVERT-001, RNA-MFE-001, RNA-PAIR-001, RNA-PARTITION-001, RNA-PSEUDOKNOT-001
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

    #region RNA-DOTBRACKET-001: RT: parse∘format = identity; P: balanced brackets → valid pairs; R: pair count ≤ len/2; D: deterministic

    // ParseDotBracket extracts base pairs from dot-bracket notation (ViennaRNA/WUSS); ValidateDotBracket
    // tests well-formedness. For nested round-bracket structures, parsing then re-rendering the pairs
    // reproduces the original string.

    /// <summary>Generates a well-formed, nested (no pseudoknot) dot-bracket string over '(', ')', '.'.</summary>
    private static Arbitrary<string> BalancedDotBracketArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int steps = 5 + rng.Next(20);
            var chars = new List<char>(steps);
            int open = 0;
            for (int i = 0; i < steps; i++)
            {
                int r = rng.Next(3);
                if (r == 0) { chars.Add('('); open++; }
                else if (r == 1 && open > 0) { chars.Add(')'); open--; }
                else chars.Add('.');
            }
            while (open-- > 0) chars.Add(')');
            return new string(chars.ToArray());
        }).ToArbitrary();

    private static string RenderPairs(int length, IEnumerable<(int Position1, int Position2)> pairs)
    {
        var arr = new char[length];
        Array.Fill(arr, '.');
        foreach (var (p1, p2) in pairs) { arr[p1] = '('; arr[p2] = ')'; }
        return new string(arr);
    }

    /// <summary>
    /// INV-1 (RT): For a nested round-bracket structure, parsing to base pairs and re-rendering those
    /// pairs reproduces the original dot-bracket string exactly (parse ∘ format = identity).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotBracket_ParseFormat_RoundTrips()
    {
        return Prop.ForAll(BalancedDotBracketArbitrary(), s =>
        {
            var pairs = RnaSecondaryStructure.ParseDotBracket(s).ToList();
            return (RenderPairs(s.Length, pairs) == s)
                .Label($"round-trip failed: '{RenderPairs(s.Length, pairs)}' ≠ '{s}'");
        });
    }

    /// <summary>
    /// INV-2 (P): A balanced structure validates, and every parsed pair opens before it closes with
    /// each position used at most once (a proper matching).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotBracket_Balanced_YieldsValidPairs()
    {
        return Prop.ForAll(BalancedDotBracketArbitrary(), s =>
        {
            var pairs = RnaSecondaryStructure.ParseDotBracket(s).ToList();
            bool ordered = pairs.All(p => p.Position1 < p.Position2);
            var positions = pairs.SelectMany(p => new[] { p.Position1, p.Position2 }).ToList();
            bool disjoint = positions.Distinct().Count() == positions.Count;
            return (RnaSecondaryStructure.ValidateDotBracket(s) && ordered && disjoint)
                .Label("balanced string failed validation or produced overlapping/inverted pairs");
        });
    }

    /// <summary>
    /// INV-3 (R): The number of base pairs is at most ⌊len/2⌋ (each pair occupies two positions).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotBracket_PairCount_AtMostHalfLength()
    {
        return Prop.ForAll(BalancedDotBracketArbitrary(), s =>
        {
            int pairs = RnaSecondaryStructure.ParseDotBracket(s).Count();
            return (2 * pairs <= s.Length).Label($"{pairs} pairs > len/2 ({s.Length})");
        });
    }

    /// <summary>
    /// INV-4 (D): Parsing is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DotBracket_Parse_IsDeterministic()
    {
        return Prop.ForAll(BalancedDotBracketArbitrary(), s =>
            RnaSecondaryStructure.ParseDotBracket(s).SequenceEqual(RnaSecondaryStructure.ParseDotBracket(s))
                .Label("ParseDotBracket must be deterministic"));
    }

    /// <summary>
    /// INV-5 (golden/negative): nested pairs are recovered; unbalanced or family-mismatched strings
    /// are rejected by validation.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DotBracket_GoldenAndInvalidCases()
    {
        Assert.Multiple(() =>
        {
            var pairs = RnaSecondaryStructure.ParseDotBracket("(())").OrderBy(p => p.Position1).ToList();
            Assert.That(pairs, Is.EqualTo(new[] { (0, 3), (1, 2) }), "nested pairs recovered");
            Assert.That(RnaSecondaryStructure.ValidateDotBracket("((..))"), Is.True);
            Assert.That(RnaSecondaryStructure.ValidateDotBracket("(()"), Is.False, "unbalanced rejected");
            Assert.That(RnaSecondaryStructure.ValidateDotBracket("(]"), Is.False, "family mismatch rejected");
            Assert.That(RnaSecondaryStructure.ValidateDotBracket("())"), Is.False, "extra close rejected");
        });
    }

    #endregion

    #region RNA-HAIRPIN-001: R: sub-minimal loops are prohibitive; M: larger destabilising loop → higher energy; D: deterministic

    // CalculateHairpinLoopEnergy returns the Turner 2004 hairpin loop ΔG. Loops shorter than the 3-nt
    // steric minimum return a prohibitive +100 kcal/mol; for very large loops the Jacobson-Stockmayer
    // extrapolation ΔG(n) = ΔG(9) + 1.75·RT·ln(n/9) grows with loop size.

    /// <summary>
    /// INV-1 (R): A loop below the 3-nt steric minimum is assigned a prohibitive (large positive)
    /// energy, while a feasible loop (≥ 3 nt) has a finite, far smaller energy.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Hairpin_SubMinimalLoop_IsProhibitive()
    {
        double tooShort = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AA", 'G', 'C');   // size 2
        double feasible = RnaSecondaryStructure.CalculateHairpinLoopEnergy("AAA", 'G', 'C');  // size 3

        Assert.Multiple(() =>
        {
            Assert.That(tooShort, Is.EqualTo(100.0), "loops < 3 nt must be prohibitive");
            Assert.That(feasible, Is.LessThan(100.0), "a 3-nt loop must be feasible");
        });
    }

    /// <summary>
    /// INV-2 (M): In the Jacobson-Stockmayer regime (loops &gt; 30 nt) a larger loop is more
    /// destabilising — energy is non-decreasing in loop size when composition and closing pair are
    /// held fixed (all-A loop, G-C closure keep the sequence-dependent terms constant).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_LargerLoop_HasHigherEnergy()
    {
        var sizes = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int a = 31 + rng.Next(120);
            int b = a + 1 + rng.Next(50);
            return (a, b);
        }).ToArbitrary();

        return Prop.ForAll(sizes, ab =>
        {
            double ea = RnaSecondaryStructure.CalculateHairpinLoopEnergy(new string('A', ab.a), 'G', 'C');
            double eb = RnaSecondaryStructure.CalculateHairpinLoopEnergy(new string('A', ab.b), 'G', 'C');
            return (eb >= ea).Label($"loop {ab.b} (ΔG={eb}) not ≥ loop {ab.a} (ΔG={ea})");
        });
    }

    /// <summary>
    /// INV-3 (D): Hairpin loop energy is deterministic and finite for feasible loops.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_Energy_IsDeterministicAndFinite()
    {
        return Prop.ForAll(RnaArbitrary(4), loop =>
        {
            double e1 = RnaSecondaryStructure.CalculateHairpinLoopEnergy(loop, 'G', 'C');
            double e2 = RnaSecondaryStructure.CalculateHairpinLoopEnergy(loop, 'G', 'C');
            return (e1 == e2 && double.IsFinite(e1))
                .Label($"hairpin energy non-deterministic or non-finite: {e1}");
        });
    }

    #endregion

    #region RNA-INVERT-001: P: arms reverse-complementary; R: positions valid; D: deterministic

    // FindInvertedRepeats reports W G W̄ᴿ patterns: a left arm, a loop, and a right arm equal to the
    // reverse complement of the left (antiparallel stem). The stem extends outward while
    // GetComplement(seq[q+k]) == seq[p−k] (Alamro et al. 2021, IUPACpal).

    private static string RandRna(Random rng, int len)
    {
        const string bases = "ACGU";
        var chars = new char[len];
        for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
        return new string(chars);
    }

    /// <summary>Builds left-arm + loop + reverseComplement(left-arm) so a perfect inverted repeat is present.</summary>
    private static Arbitrary<string> EmbeddedInvertedRepeatArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            string left = RandRna(rng, 4 + rng.Next(5));      // arm length 4..8
            string right = new string(left.Reverse().Select(RnaSecondaryStructure.GetComplement).ToArray());
            string loop = RandRna(rng, 3 + rng.Next(4));       // loop length 3..6
            return left + loop + right;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P + R): every reported inverted repeat has valid, non-overlapping arm coordinates of
    /// equal length with a loop in [minSpacing, maxSpacing], and its two arms are antiparallel
    /// reverse complements (GetComplement(seq[Start2+m]) == seq[End1−m]). A sequence built to contain
    /// a perfect stem yields at least one such repeat (non-vacuous).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeats_AreValidReverseComplementaryArms()
    {
        return Prop.ForAll(EmbeddedInvertedRepeatArbitrary(), seq =>
        {
            var reps = RnaSecondaryStructure.FindInvertedRepeats(seq).ToList();
            if (reps.Count == 0)
                return false.Label("expected at least one inverted repeat in a constructed stem");

            bool allValid = reps.All(r =>
                0 <= r.Start1 && r.Start1 <= r.End1 && r.End1 < r.Start2 && r.Start2 <= r.End2 && r.End2 < seq.Length
                && r.End1 - r.Start1 + 1 == r.Length && r.End2 - r.Start2 + 1 == r.Length
                && r.Length >= 4
                && (r.Start2 - r.End1 - 1) >= 3 && (r.Start2 - r.End1 - 1) <= 100
                && Enumerable.Range(0, r.Length)
                    .All(m => RnaSecondaryStructure.GetComplement(seq[r.Start2 + m]) == seq[r.End1 - m]));
            return allValid.Label("a reported repeat had invalid coordinates or non-complementary arms");
        });
    }

    /// <summary>
    /// INV-2 (D): Inverted-repeat detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeats_AreDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(20), seq =>
            RnaSecondaryStructure.FindInvertedRepeats(seq).SequenceEqual(RnaSecondaryStructure.FindInvertedRepeats(seq))
                .Label("FindInvertedRepeats must be deterministic"));
    }

    /// <summary>
    /// INV-3 (golden): the constructed stem AAGG-UUU-CCUU is detected as a single inverted repeat
    /// with reverse-complementary arms.
    /// </summary>
    [Test]
    [Category("Property")]
    public void InvertedRepeats_GoldenStem_IsDetected()
    {
        const string seq = "AAGGUUUCCUU"; // arm AAGG, loop UUU, arm CCUU = revcomp(AAGG)
        var reps = RnaSecondaryStructure.FindInvertedRepeats(seq).ToList();

        Assert.That(reps, Is.Not.Empty);
        var r = reps[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.Length, Is.GreaterThanOrEqualTo(4));
            for (int m = 0; m < r.Length; m++)
                Assert.That(RnaSecondaryStructure.GetComplement(seq[r.Start2 + m]), Is.EqualTo(seq[r.End1 - m]),
                    $"arm position {m} not reverse-complementary");
        });
    }

    #endregion

    #region RNA-MFE-001: R: MFE ≤ 0; M: more GC pairs → lower energy; D: deterministic

    // CalculateMinimumFreeEnergy is a Zuker-style DP with Turner 2004 parameters. The unfolded
    // structure always scores 0, so the MFE can never be positive; a GC stem (≈ −3 kcal/mol/stack)
    // is more stable than an equivalent AU stem (≈ −1 to −2).

    /// <summary>
    /// INV-1 (R): The MFE is never positive — the empty (unfolded) structure is always available at 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Mfe_IsNonPositive()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
        {
            double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            return (mfe <= 1e-9).Label($"MFE={mfe} must be ≤ 0");
        });
    }

    /// <summary>
    /// INV-2 (M): A hairpin closed by a GC stem has a lower (more negative) MFE than the same hairpin
    /// closed by an AU stem of equal length — more GC pairs stabilise the structure.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Mfe_GcStem_IsMoreStableThanAuStem()
    {
        double mfeGc = RnaSecondaryStructure.CalculateMinimumFreeEnergy("GGGGGUUUUCCCCC"); // 5-bp G-C stem
        double mfeAu = RnaSecondaryStructure.CalculateMinimumFreeEnergy("AAAAAGGGUUUUU");   // 5-bp A-U stem

        Assert.Multiple(() =>
        {
            Assert.That(mfeGc, Is.LessThan(0.0), "a GC stem should fold to a negative MFE");
            Assert.That(mfeGc, Is.LessThan(mfeAu), "the GC stem must be more stable than the AU stem");
        });
    }

    /// <summary>
    /// INV-3 (D): MFE computation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Mfe_IsDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
            (RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq) == RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq))
                .Label("CalculateMinimumFreeEnergy must be deterministic"));
    }

    #endregion

    #region RNA-PAIR-001: P: only A-U, G-C, G-U pair; S: canPair(a,b)=canPair(b,a); D: deterministic

    // CanPair recognises the six canonical RNA base pairs — Watson-Crick A-U/G-C and the G-U wobble —
    // in either order and case-insensitively. GetBasePairType classifies them as WatsonCrick / Wobble.

    /// <summary>Bases plus a few non-RNA / mixed-case symbols to exercise the negative cases.</summary>
    private static Arbitrary<char> BaseCharArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'U', 'T', 'N', 'a', 'c', 'g', 'u').ToArbitrary();

    /// <summary>
    /// INV-1 (P): Over {A,C,G,U}, CanPair is true exactly for the six canonical pairs (A-U, U-A,
    /// G-C, C-G, G-U, U-G) and false for everything else (e.g. A-A, G-A); DNA 'T' does not pair.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Pair_OnlyCanonicalPairs_AreAccepted()
    {
        var canonical = new HashSet<(char, char)>
        {
            ('A','U'), ('U','A'), ('G','C'), ('C','G'), ('G','U'), ('U','G')
        };
        foreach (char a in "ACGU")
            foreach (char b in "ACGU")
                Assert.That(RnaSecondaryStructure.CanPair(a, b), Is.EqualTo(canonical.Contains((a, b))),
                    $"CanPair({a},{b}) mismatch");

        Assert.That(RnaSecondaryStructure.CanPair('A', 'T'), Is.False, "RNA pairing does not accept DNA T");
    }

    /// <summary>
    /// INV-2 (S): Pairing is symmetric — CanPair(a,b) == CanPair(b,a) for any bases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pair_IsSymmetric()
    {
        return Prop.ForAll(BaseCharArbitrary(), BaseCharArbitrary(), (a, b) =>
            (RnaSecondaryStructure.CanPair(a, b) == RnaSecondaryStructure.CanPair(b, a))
                .Label($"CanPair({a},{b}) ≠ CanPair({b},{a})"));
    }

    /// <summary>
    /// INV-3 (P, case-insensitive): pairing ignores case.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pair_IsCaseInsensitive()
    {
        return Prop.ForAll(BaseCharArbitrary(), BaseCharArbitrary(), (a, b) =>
        {
            char ta = char.IsUpper(a) ? char.ToLowerInvariant(a) : char.ToUpperInvariant(a);
            char tb = char.IsUpper(b) ? char.ToLowerInvariant(b) : char.ToUpperInvariant(b);
            return (RnaSecondaryStructure.CanPair(a, b) == RnaSecondaryStructure.CanPair(ta, tb))
                .Label($"case sensitivity at ({a},{b})");
        });
    }

    /// <summary>
    /// INV-4 (P, classification): GetBasePairType agrees with CanPair and labels Watson-Crick vs
    /// Wobble correctly — A-U/G-C are WatsonCrick, G-U is Wobble, all others null.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pair_TypeClassification_IsConsistent()
    {
        return Prop.ForAll(BaseCharArbitrary(), BaseCharArbitrary(), (a, b) =>
        {
            var type = RnaSecondaryStructure.GetBasePairType(a, b);
            bool can = RnaSecondaryStructure.CanPair(a, b);
            char ua = char.ToUpperInvariant(a), ub = char.ToUpperInvariant(b);
            bool isWobble = (ua == 'G' && ub == 'U') || (ua == 'U' && ub == 'G');
            bool expectedWc = can && !isWobble;

            bool ok = (can == (type != null))
                      && (!can || (isWobble
                            ? type == RnaSecondaryStructure.BasePairType.Wobble
                            : type == RnaSecondaryStructure.BasePairType.WatsonCrick))
                      && (!expectedWc || type == RnaSecondaryStructure.BasePairType.WatsonCrick);
            return ok.Label($"type inconsistent for ({a},{b}): can={can}, type={type}");
        });
    }

    #endregion

    #region RNA-PARTITION-001: R: Z > 0; R: base-pair probability ∈ [0,1]; D: deterministic

    // CalculatePartitionFunction implements McCaskill (1990): the inside partition function Z (≥ 1,
    // since the all-unpaired structure always weighs 1) and equilibrium base-pair probabilities
    // P[i,j] ∈ [0,1] with Σ_j P[i,j] ≤ 1 for each position (a base pairs with at most one partner).

    /// <summary>
    /// INV-1 (R): The partition function is at least 1 (the empty structure is always in the ensemble).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Partition_Z_IsAtLeastOne()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
        {
            double z = RnaSecondaryStructure.CalculatePartitionFunction(seq).PartitionFunction;
            return (z >= 1.0 - 1e-9).Label($"Z={z} must be ≥ 1");
        });
    }

    /// <summary>
    /// INV-2 (R): Every equilibrium base-pair probability lies in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Partition_Probabilities_InUnitInterval()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
        {
            var probs = RnaSecondaryStructure.CalculatePartitionFunction(seq).BasePairProbabilities;
            return probs.Values.All(p => p >= -1e-9 && p <= 1.0 + 1e-9)
                .Label("a base-pair probability fell outside [0,1]");
        });
    }

    /// <summary>
    /// INV-3 (P): For each position the probabilities of all pairs involving it sum to ≤ 1 — a base
    /// can be paired with at most one partner in any structure.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Partition_PerBasePairing_AtMostOne()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
        {
            var probs = RnaSecondaryStructure.CalculatePartitionFunction(seq).BasePairProbabilities;
            var perBase = new Dictionary<int, double>();
            foreach (var kv in probs)
            {
                perBase[kv.Key.I] = perBase.GetValueOrDefault(kv.Key.I) + kv.Value;
                perBase[kv.Key.J] = perBase.GetValueOrDefault(kv.Key.J) + kv.Value;
            }
            return perBase.Values.All(p => p <= 1.0 + 1e-9)
                .Label("a base's total pairing probability exceeded 1");
        });
    }

    /// <summary>
    /// INV-4 (D): The partition function and probabilities are deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Partition_IsDeterministic()
    {
        return Prop.ForAll(RnaArbitrary(8), seq =>
        {
            var a = RnaSecondaryStructure.CalculatePartitionFunction(seq);
            var b = RnaSecondaryStructure.CalculatePartitionFunction(seq);
            bool same = a.PartitionFunction == b.PartitionFunction
                        && a.BasePairProbabilities.Count == b.BasePairProbabilities.Count
                        && a.BasePairProbabilities.All(kv => b.BasePairProbabilities.TryGetValue(kv.Key, out double v) && v == kv.Value);
            return same.Label("CalculatePartitionFunction must be deterministic");
        });
    }

    /// <summary>
    /// INV-5 (boundary): empty sequence has Z=1 and no pairs; a non-positive temperature is rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Partition_Boundaries()
    {
        Assert.Multiple(() =>
        {
            var empty = RnaSecondaryStructure.CalculatePartitionFunction("");
            Assert.That(empty.PartitionFunction, Is.EqualTo(1.0));
            Assert.That(empty.BasePairProbabilities, Is.Empty);
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RnaSecondaryStructure.CalculatePartitionFunction("GGGAAACCC", temperature: 0));
        });
    }

    #endregion

    #region RNA-PSEUDOKNOT-001: P: detects exactly the crossing pairs; R: positions valid; D: deterministic

    // DetectPseudoknots reports each crossing pair-of-pairs: normalising endpoints to (open<close) and
    // ordering by opening position, two pairs cross iff i < k < j < l (Antczak et al. 2018). Nested or
    // disjoint pairs are not pseudoknots.

    private static RnaSecondaryStructure.BasePair Bp(int p1, int p2) =>
        new(p1, p2, 'A', 'U', RnaSecondaryStructure.BasePairType.WatsonCrick);

    private static bool Cross(RnaSecondaryStructure.BasePair x, RnaSecondaryStructure.BasePair y)
    {
        int i = Math.Min(x.Position1, x.Position2), j = Math.Max(x.Position1, x.Position2);
        int k = Math.Min(y.Position1, y.Position2), l = Math.Max(y.Position1, y.Position2);
        if (k < i) (i, j, k, l) = (k, l, i, j);
        return i < k && k < j && j < l;
    }

    /// <summary>Generates 2..6 base pairs over positions 0..19 (distinct endpoints per pair).</summary>
    private static Arbitrary<List<RnaSecondaryStructure.BasePair>> BasePairsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int m = 2 + rng.Next(5);
            var pairs = new List<RnaSecondaryStructure.BasePair>(m);
            for (int t = 0; t < m; t++)
            {
                int a = rng.Next(20), b = rng.Next(20);
                while (a == b) b = rng.Next(20);
                pairs.Add(Bp(a, b));
            }
            return pairs;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P): The set of reported pseudoknots is exactly the set of crossing pair-of-pairs,
    /// verified against an independent crossing predicate, and each carries its two crossing pairs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pseudoknot_DetectsExactlyCrossingPairs()
    {
        return Prop.ForAll(BasePairsArbitrary(), pairs =>
        {
            int expected = 0;
            for (int a = 0; a < pairs.Count; a++)
                for (int b = a + 1; b < pairs.Count; b++)
                    if (Cross(pairs[a], pairs[b])) expected++;

            var detected = RnaSecondaryStructure.DetectPseudoknots(pairs).ToList();
            bool eachHasTwo = detected.All(pk => pk.CrossingPairs.Count == 2);
            return (detected.Count == expected && eachHasTwo)
                .Label($"detected {detected.Count} pseudoknots, expected {expected}");
        });
    }

    /// <summary>
    /// INV-2 (R): Every reported pseudoknot has the crossing layout Start1 &lt; Start2 &lt; End1 &lt; End2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pseudoknot_Positions_AreCrossing()
    {
        return Prop.ForAll(BasePairsArbitrary(), pairs =>
            RnaSecondaryStructure.DetectPseudoknots(pairs)
                .All(pk => pk.Start1 < pk.Start2 && pk.Start2 < pk.End1 && pk.End1 < pk.End2)
                .Label("a pseudoknot did not satisfy i < k < j < l"));
    }

    /// <summary>
    /// INV-3 (D): Pseudoknot detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Pseudoknot_IsDeterministic()
    {
        return Prop.ForAll(BasePairsArbitrary(), pairs =>
            RnaSecondaryStructure.DetectPseudoknots(pairs).Select(pk => (pk.Start1, pk.End1, pk.Start2, pk.End2))
                .SequenceEqual(RnaSecondaryStructure.DetectPseudoknots(pairs).Select(pk => (pk.Start1, pk.End1, pk.Start2, pk.End2)))
                .Label("DetectPseudoknots must be deterministic"));
    }

    /// <summary>
    /// INV-4 (golden): crossing pairs are flagged; nested and disjoint pairs are not.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Pseudoknot_GoldenCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RnaSecondaryStructure.DetectPseudoknots(new[] { Bp(0, 5), Bp(3, 8) }).Count(), Is.EqualTo(1),
                "0<3<5<8 crosses → pseudoknot");
            Assert.That(RnaSecondaryStructure.DetectPseudoknots(new[] { Bp(0, 10), Bp(2, 8) }), Is.Empty,
                "nested pairs are not a pseudoknot");
            Assert.That(RnaSecondaryStructure.DetectPseudoknots(new[] { Bp(0, 3), Bp(5, 8) }), Is.Empty,
                "disjoint pairs are not a pseudoknot");
        });
    }

    #endregion

    #region RNA-PKPREDICT-001: R: ΔG ≤ pseudoknot-free MFE; P: each base paired ≤ once; P: crossings genuine; D: deterministic

    // PredictStructurePseudoknot (H-type pseudoknot search; Rivas & Eddy 1999) takes the plain MFE
    // structure as its baseline/fallback and only accepts a crossing motif when it lowers the energy.

    // The pseudoknot search is O(n^4); cap the generated length so the property suite stays fast while
    // still spanning lengths that can form an H-type knot (≥ 11 nt).
    private static Arbitrary<string> BoundedRnaArbitrary(int minLen = 11, int maxLen = 30) =>
        (from len in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf(len)
         select new string(chars)).ToArbitrary();

    private static bool EachBasePairedAtMostOnce(IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        var seen = new HashSet<int>();
        foreach (var (p1, p2) in pairs)
        {
            if (!seen.Add(p1) || !seen.Add(p2)) return false;
        }
        return true;
    }

    private static bool HasCrossingPair(IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        for (int a = 0; a < pairs.Count; a++)
        {
            int i = Math.Min(pairs[a].Position1, pairs[a].Position2);
            int j = Math.Max(pairs[a].Position1, pairs[a].Position2);
            for (int b = a + 1; b < pairs.Count; b++)
            {
                int k = Math.Min(pairs[b].Position1, pairs[b].Position2);
                int l = Math.Max(pairs[b].Position1, pairs[b].Position2);
                if (i < k && k < j && j < l) return true; // genuine crossing i<k<j<l
            }
        }
        return false;
    }

    /// <summary>
    /// INV-1 (R): the pseudoknot prediction never scores worse than the pseudoknot-free MFE — the plain
    /// MFE structure is its baseline and a crossing motif is accepted only when it lowers ΔG.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknot_Energy_AtMostPseudoknotFreeMfe()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            double pk = RnaSecondaryStructure.PredictStructurePseudoknot(seq).FreeEnergy;
            double mfe = RnaSecondaryStructure.CalculateMfeStructure(seq).FreeEnergy;
            return (pk <= mfe + 1e-9).Label($"pk ΔG={pk} must be ≤ pseudoknot-free MFE={mfe}");
        });
    }

    /// <summary>
    /// INV-2 (P): the predicted structure is a matching — every nucleotide participates in at most one
    /// base pair (no position appears twice across all pairs).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknot_EachBase_PairedAtMostOnce()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            var s = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            return EachBasePairedAtMostOnce(s.BasePairs)
                .Label("a nucleotide participates in more than one base pair");
        });
    }

    /// <summary>
    /// INV-3 (P): the pseudoknot flag is genuine — HasPseudoknot is true if and only if the predicted
    /// base pairs actually contain a crossing pair (i &lt; k &lt; j &lt; l).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknot_Flag_IffGenuineCrossing()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            var s = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            return (s.HasPseudoknot == HasCrossingPair(s.BasePairs))
                .Label($"HasPseudoknot={s.HasPseudoknot} but crossing={HasCrossingPair(s.BasePairs)}");
        });
    }

    /// <summary>INV-4 (D): pseudoknot prediction is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknot_IsDeterministic()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            var a = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            var b = RnaSecondaryStructure.PredictStructurePseudoknot(seq);
            return (a.DotBracket == b.DotBracket && Math.Abs(a.FreeEnergy - b.FreeEnergy) < 1e-12
                    && a.HasPseudoknot == b.HasPseudoknot && a.BasePairs.SequenceEqual(b.BasePairs))
                .Label("PredictStructurePseudoknot must be deterministic");
        });
    }

    #endregion

    #region RNA-PKRECURSIVE-001: R: ΔG ≤ single-knot result; P: valid nested structure; D: deterministic

    // PredictStructurePseudoknotRecursive recursively folds the gap regions inside a knot, so it can
    // only add stabilising pairs — its energy is no higher than the single-knot PredictStructurePseudoknot.

    /// <summary>
    /// INV-1 (R): the recursive pseudoknot fold never scores worse than the single-knot prediction —
    /// recursively folding the loop regions can only lower (or hold) ΔG.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknotRecursive_Energy_AtMostSingleKnot()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            double recursive = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq).FreeEnergy;
            double single = RnaSecondaryStructure.PredictStructurePseudoknot(seq).FreeEnergy;
            return (recursive <= single + 1e-9).Label($"recursive ΔG={recursive} must be ≤ single-knot ΔG={single}");
        });
    }

    /// <summary>
    /// INV-2 (P): the recursive structure is well formed — every nucleotide is paired at most once.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknotRecursive_EachBase_PairedAtMostOnce()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            var s = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq);
            return EachBasePairedAtMostOnce(s.BasePairs)
                .Label("a nucleotide participates in more than one base pair");
        });
    }

    /// <summary>INV-3 (D): recursive pseudoknot prediction is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property PredictPseudoknotRecursive_IsDeterministic()
    {
        return Prop.ForAll(BoundedRnaArbitrary(), seq =>
        {
            var a = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq);
            var b = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq);
            return (a.DotBracket == b.DotBracket && Math.Abs(a.FreeEnergy - b.FreeEnergy) < 1e-12
                    && a.HasPseudoknot == b.HasPseudoknot && a.BasePairs.SequenceEqual(b.BasePairs))
                .Label("PredictStructurePseudoknotRecursive must be deterministic");
        });
    }

    #endregion

    #region RNA-ACCESS-001: R: 0 ≤ P_unpaired ≤ 1; M: longer region → lower P_unpaired; D: deterministic

    // CalculateRegionUnpairedProbability returns the McCaskill equilibrium probability that an entire
    // window [windowEnd−windowLength+1, windowEnd] is unpaired (RNAplfold; Bernhart et al. 2006).

    private const int AccessWindow = 4;

    /// <summary>
    /// INV-1 (R): the region-unpaired probability is a probability in [0,1] for any in-bounds window.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RegionUnpairedProbability_InUnitInterval()
    {
        var gen = (from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf().Where(a => a.Length >= 12)
                   let seq = new string(chars)
                   from windowEnd in Gen.Choose(AccessWindow - 1, seq.Length - 1)
                   select (seq, windowEnd)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double p = RnaSecondaryStructure.CalculateRegionUnpairedProbability(t.seq, t.windowEnd, AccessWindow);
            return (p >= -1e-9 && p <= 1.0 + 1e-9).Label($"P_unpaired={p} outside [0,1]");
        });
    }

    /// <summary>
    /// INV-2 (M): a longer (superset) region is unpaired with no higher probability — extending the
    /// window at a fixed right edge from length L to L+1 cannot raise the all-unpaired probability.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RegionUnpairedProbability_LongerRegion_NotMoreProbable()
    {
        var gen = (from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf().Where(a => a.Length >= 14)
                   let seq = new string(chars)
                   from len in Gen.Choose(2, 6)
                   from windowEnd in Gen.Choose(len, seq.Length - 1) // ensures windowStart ≥ 0 for len+1
                   select (seq, windowEnd, len)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double pShort = RnaSecondaryStructure.CalculateRegionUnpairedProbability(t.seq, t.windowEnd, t.len);
            double pLong = RnaSecondaryStructure.CalculateRegionUnpairedProbability(t.seq, t.windowEnd, t.len + 1);
            return (pLong <= pShort + 1e-9)
                .Label($"longer region more probable: P(len {t.len + 1})={pLong} > P(len {t.len})={pShort}");
        });
    }

    /// <summary>INV-3 (D): the region-unpaired probability is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property RegionUnpairedProbability_IsDeterministic()
    {
        var gen = (from chars in Gen.Elements('A', 'C', 'G', 'U').ArrayOf().Where(a => a.Length >= 12)
                   let seq = new string(chars)
                   from windowEnd in Gen.Choose(AccessWindow - 1, seq.Length - 1)
                   select (seq, windowEnd)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(t.seq, t.windowEnd, AccessWindow);
            double b = RnaSecondaryStructure.CalculateRegionUnpairedProbability(t.seq, t.windowEnd, AccessWindow);
            return (a == b).Label("CalculateRegionUnpairedProbability must be deterministic");
        });
    }

    #endregion
}
