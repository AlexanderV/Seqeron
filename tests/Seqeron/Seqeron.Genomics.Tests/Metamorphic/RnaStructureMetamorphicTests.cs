namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the RnaStructure area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RNA-STRUCT-001 — secondary-structure prediction (RnaStructure).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 71.
///
/// API under test (RnaSecondaryStructure.PredictStructure):
///   Finds complementary stem-loops, greedily selects a non-overlapping set by energy, and
///   reports the resulting base pairs (i,j), the dot-bracket string and the total MFE.
///   Base pairing requires Watson–Crick / wobble complementarity, a minimum stem length and
///   a loop within [minLoopSize, maxLoopSize].
///
/// Relations (derived from complementary pairing, NOT from output):
///   • COMP (empty ⇒ 0 pairs): an empty sequence — and any sequence with no complementary
///          run able to close a stem (e.g. a homopolymer) — yields no base pairs.
///   • MON  (more complementary bases ⇒ ≥ pairs): lengthening the complementary arms of a
///          hairpin can only add base pairs, so the pair count is non-decreasing in arm length.
///   • INV  (non-pairing 3' tail ⇒ existing pairs preserved): appending a tail that cannot
///          pair (a homopolymer) at the 3' end leaves every base pair of the original hairpin
///          intact at the same positions.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class RnaStructureMetamorphicTests
{
    #region Helpers

    /// <summary>The set of base-pair index pairs (i,j) predicted for a sequence.</summary>
    private static HashSet<(int I, int J)> PairSet(string rna) =>
        RnaSecondaryStructure.PredictStructure(rna)
            .BasePairs.Select(bp => (bp.Position1, bp.Position2))
            .ToHashSet();

    /// <summary>A perfect GC hairpin: an n-base G arm, a 4-base poly-A loop, an n-base C arm.</summary>
    private static string GcHairpin(int armLength) =>
        new string('G', armLength) + "AAAA" + new string('C', armLength);

    #endregion

    #region COMP — no complementarity ⇒ no base pairs

    [Test]
    [Description("COMP: an empty sequence yields no base pairs, an empty dot-bracket and zero MFE.")]
    public void PredictStructure_EmptySequence_HasNoPairs()
    {
        var result = RnaSecondaryStructure.PredictStructure("");

        result.BasePairs.Should().BeEmpty(because: "there are no bases to pair");
        result.DotBracket.Should().BeEmpty(because: "an empty sequence has an empty structure string");
        result.MinimumFreeEnergy.Should().Be(0, because: "no pairs contribute any stacking energy");
    }

    [Test]
    [Description("COMP: a homopolymer cannot form any complementary stem, so no base pairs are predicted.")]
    public void PredictStructure_Homopolymer_HasNoPairs()
    {
        foreach (var homopolymer in new[] { "AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC" })
            RnaSecondaryStructure.PredictStructure(homopolymer)
                .BasePairs.Should().BeEmpty(
                    because: $"'{homopolymer[0]}' cannot pair with itself, so no stem can close");
    }

    #endregion

    #region MON — longer complementary arms cannot reduce the pair count

    [Test]
    [Description("MON: extending the complementary arms of a hairpin can only add base pairs — the pair count is non-decreasing in arm length and strictly grows over the range.")]
    public void PredictStructure_LongerComplementaryArms_DoNotReducePairs()
    {
        int previous = -1;
        var counts = new List<int>();

        foreach (int arm in new[] { 3, 4, 5, 6, 7 })
        {
            int pairs = RnaSecondaryStructure.PredictStructure(GcHairpin(arm)).BasePairs.Count;
            counts.Add(pairs);

            pairs.Should().BeGreaterThanOrEqualTo(previous,
                because: $"a longer ({arm}-base) complementary arm keeps every pair the shorter arm could form and may add more");
            previous = pairs;
        }

        counts.Last().Should().BeGreaterThan(counts.First(),
            because: "the 7-base stem forms strictly more base pairs than the 3-base stem");
    }

    #endregion

    #region INV — a non-pairing 3' tail preserves existing pairs

    [Test]
    [Description("INV: appending a homopolymer tail (which cannot pair) at the 3' end leaves every base pair of the original hairpin intact at the same positions.")]
    public void PredictStructure_NonPairing3PrimeTail_PreservesExistingPairs()
    {
        const string core = "GGGGGAAAACCCCC"; // a 5-bp GC hairpin closing a poly-A loop
        var corePairs = PairSet(core);
        corePairs.Should().NotBeEmpty(because: "the core hairpin must actually fold for the test to be meaningful");

        // A poly-A tail at the 3' end cannot pair (no U to pair with) and does not shift the
        // hairpin's positions, so the original pairs must survive unchanged.
        foreach (var tail in new[] { "A", "AAAA", "AAAAAAAA" })
        {
            var extendedPairs = PairSet(core + tail);

            extendedPairs.IsSupersetOf(corePairs).Should().BeTrue(
                because: $"a non-pairing 3' tail of {tail.Length} A's cannot break or move the existing base pairs");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-STEMLOOP-001 — stem-loop (hairpin) detection (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 72.
    //
    // API under test (RnaSecondaryStructure.FindStemLoops):
    //   For each loop window, extends a stem outward as long as the flanking bases are
    //   complementary (Watson–Crick or wobble), returning a StemLoop whose Stem.Length is the
    //   number of consecutive complementary pairs. The stem is built only from the arms; the
    //   loop bases lie strictly between the arms and are never part of the stem.
    //
    // Relations (derived from outward complementary extension, NOT from output):
    //   • COMP (no complement ⇒ no stem): a sequence with no complementary arms (a homopolymer,
    //          or like-base "arms") yields no stem-loops.
    //   • MON  (longer arms ⇒ longer stem): lengthening the complementary arms of a hairpin
    //          lengthens the maximal stem — its length is non-decreasing and strictly grows.
    //   • INV  (loop content irrelevant to stem pairing): with the arms and loop length fixed,
    //          changing the loop bases leaves the hairpin's stem base-pairs identical, since the
    //          stem is built from the arms alone.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-STEMLOOP-001 — Helpers

    private static HashSet<(int I, int J)> StemPairSet(RnaSecondaryStructure.StemLoop sl) =>
        sl.Stem.BasePairs.Select(bp => (bp.Position1, bp.Position2)).ToHashSet();

    #endregion

    #region RNA-STEMLOOP-001 COMP — no complementary arms ⇒ no stem-loop

    [Test]
    [Description("COMP: a sequence with no complementary arms (a homopolymer or like-base arms) produces no stem-loops.")]
    public void FindStemLoops_NoComplementarity_ReturnsNothing()
    {
        foreach (var seq in new[] { "AAAAAAAAAAAA", "GGGGAAAGGGG" })
            RnaSecondaryStructure.FindStemLoops(seq).Should().BeEmpty(
                because: $"'{seq}' has no two complementary arms able to close a stem");
    }

    #endregion

    #region RNA-STEMLOOP-001 MON — longer complementary arms lengthen the maximal stem

    [Test]
    [Description("MON: lengthening a hairpin's complementary arms lengthens the maximal stem — non-decreasing in arm length and strictly larger over the range.")]
    public void FindStemLoops_LongerArms_LengthenMaximalStem()
    {
        int previous = -1;
        var maxima = new List<int>();

        foreach (int arm in new[] { 3, 4, 5, 6, 7 })
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(GcHairpin(arm)).ToList();
            int maxStem = stemLoops.Count == 0 ? 0 : stemLoops.Max(sl => sl.Stem.Length);
            maxima.Add(maxStem);

            maxStem.Should().BeGreaterThanOrEqualTo(previous,
                because: $"a longer ({arm}-base) complementary arm cannot shorten the maximal stem");
            previous = maxStem;
        }

        maxima.Last().Should().BeGreaterThan(maxima.First(),
            because: "the 7-base arms yield a strictly longer maximal stem than the 3-base arms");
    }

    #endregion

    #region RNA-STEMLOOP-001 INV — loop content does not change the stem base-pairs

    [Test]
    [Description("INV: with the arms (5×G / 5×C) and loop length (4) fixed, changing the loop bases leaves the hairpin's stem base-pairs identical, because the stem is built from the arms alone.")]
    public void FindStemLoops_LoopContent_DoesNotChangeStemPairing()
    {
        // Arms G(0..4) / C(9..13) close a length-4 loop at 5..8 → stem pairs (i, 13-i), i=0..4.
        var expected = Enumerable.Range(0, 5).Select(i => (i, 13 - i)).ToHashSet();

        foreach (var loop in new[] { "AAAA", "GAAA", "UUCG", "CUUG" })
        {
            string seq = "GGGGG" + loop + "CCCCC";
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq).ToList();

            stemLoops.Any(sl => StemPairSet(sl).SetEquals(expected)).Should().BeTrue(
                because: $"the 5-bp G:C stem is determined by the arms, so loop '{loop}' yields the same stem base-pairs");
            stemLoops.Max(sl => sl.Stem.Length).Should().Be(5,
                because: $"the arms are 5 bases long, so the maximal stem is length 5 regardless of loop '{loop}'");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-ENERGY-001 — stem free-energy (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 73.
    //
    // API under test (RnaSecondaryStructure.CalculateStemEnergy):
    //   Sums the Turner-2004 nearest-neighbour stacking energy over consecutive base pairs of
    //   a helix, plus a +0.45 terminal AU/GU penalty at each end that closes with an AU/GU
    //   pair. The energy is a pure function of the ordered base-pair list (more negative ⇒
    //   more stable). G:C steps are far more stabilising than A:U steps.
    //
    // Relations (derived from the additive nearest-neighbour model, NOT from output):
    //   • INV  (same structure ⇒ same ΔG): the energy is deterministic and depends only on the
    //          base identities/order, not on the absolute positions of the pairs.
    //   • MON  (more G:C stacks ⇒ lower ΔG): replacing A:U steps by G:C steps strictly lowers
    //          the energy (all-G:C < mixed < all-A:U at equal length).
    //   • COMP (additivity): for an all-G:C helix (no terminal penalty), the total energy equals
    //          the sum of its consecutive two-pair stacking-step energies.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-ENERGY-001 — Helpers

    /// <summary>Builds a helix base-pair list (5'→3') from (Base1,Base2) pairs, positions offset-shifted.</summary>
    private static List<RnaSecondaryStructure.BasePair> Helix(int offset, params (char A, char B)[] pairs)
    {
        int n = pairs.Length;
        var list = new List<RnaSecondaryStructure.BasePair>(n);
        for (int i = 0; i < n; i++)
        {
            var (a, b) = pairs[i];
            var type = (a == 'G' && b == 'U') || (a == 'U' && b == 'G')
                ? RnaSecondaryStructure.BasePairType.Wobble
                : RnaSecondaryStructure.BasePairType.WatsonCrick;
            list.Add(new RnaSecondaryStructure.BasePair(offset + i, offset + 2 * n - 1 - i, a, b, type));
        }
        return list;
    }

    private static (char, char)[] Repeat((char A, char B) pair, int count) =>
        Enumerable.Repeat(pair, count).ToArray();

    private static double StemEnergy(List<RnaSecondaryStructure.BasePair> pairs) =>
        RnaSecondaryStructure.CalculateStemEnergy("", pairs);

    #endregion

    #region RNA-ENERGY-001 INV — energy is deterministic and position-independent

    [Test]
    [Description("INV: the stem energy is a pure function of the ordered base-pair list — identical across repeated calls and unchanged when all pair positions are shifted by a constant.")]
    public void StemEnergy_SameStructure_SameEnergy()
    {
        var pairs = new[] { ('G', 'C'), ('A', 'U'), ('G', 'C'), ('C', 'G') };

        double first = StemEnergy(Helix(0, pairs));
        double again = StemEnergy(Helix(0, pairs));
        double shifted = StemEnergy(Helix(100, pairs)); // same bases/order, different absolute positions

        again.Should().Be(first, because: "the energy function has no hidden state");
        shifted.Should().Be(first, because: "stacking energy depends on the base identities and order, not on absolute positions");
    }

    #endregion

    #region RNA-ENERGY-001 MON — more G:C stacks lower the free energy

    [Test]
    [Description("MON: at equal helix length, replacing A:U steps with G:C steps strictly lowers ΔG — all-G:C < mixed < all-A:U.")]
    public void StemEnergy_MoreGcStacks_LowerEnergy()
    {
        const int n = 6;

        double allGc = StemEnergy(Helix(0, Repeat(('G', 'C'), n)));
        double allAu = StemEnergy(Helix(0, Repeat(('A', 'U'), n)));
        double mixed = StemEnergy(Helix(0,
            Repeat(('G', 'C'), n / 2).Concat(Repeat(('A', 'U'), n / 2)).ToArray()));

        allGc.Should().BeLessThan(mixed,
            because: "an all-G:C helix stacks more strongly (and pays no terminal AU/GU penalty) than a half-A:U helix");
        mixed.Should().BeLessThan(allAu,
            because: "replacing A:U steps with G:C steps lowers the free energy");
    }

    #endregion

    #region RNA-ENERGY-001 COMP — total energy is the sum of its stacking steps

    [Test]
    [Description("COMP: for an all-G:C helix (which carries no terminal penalty), the total stem energy equals the sum of its consecutive two-pair stacking-step energies.")]
    public void StemEnergy_AllGc_IsSumOfStackingSteps()
    {
        var helix = Helix(0, Repeat(('G', 'C'), 6));

        double total = StemEnergy(helix);

        double sumOfSteps = 0;
        for (int i = 0; i < helix.Count - 1; i++)
            sumOfSteps += StemEnergy(new List<RnaSecondaryStructure.BasePair> { helix[i], helix[i + 1] });

        total.Should().BeApproximately(sumOfSteps, 1e-6,
            because: "the nearest-neighbour model is additive over dinucleotide steps when no terminal penalty intervenes");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-DOTBRACKET-001 — dot-bracket notation (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 149.
    //
    // API under test (RnaSecondaryStructure.ParseDotBracket):
    //   Parses dot-bracket notation to base-pair (i,j) index tuples; the canonical formatter places
    //   '(' at the 5' partner, ')' at the 3' partner and '.' elsewhere (ViennaRNA convention).
    //
    // Relations (derived from the bracket-matching definition, NOT from output):
    //   • RT   (parse ∘ format identity): for a canonical nested notation, formatting the parsed
    //          base pairs back reproduces the original string exactly.
    //   • INV  (pairing preserved under reparse): re-parsing the reformatted string yields the same
    //          base-pair set, so the pairing survives a parse/format/reparse cycle.
    // ───────────────────────────────────────────────────────────────────────────

    // Canonical dot-bracket formatter (ViennaRNA: '(' at the 5' partner, ')' at the 3' partner).
    private static string FormatDotBracket(int length, IEnumerable<(int Position1, int Position2)> pairs)
    {
        var notation = Enumerable.Repeat('.', length).ToArray();
        foreach (var (p1, p2) in pairs)
        {
            notation[System.Math.Min(p1, p2)] = '(';
            notation[System.Math.Max(p1, p2)] = ')';
        }
        return new string(notation);
    }

    private static HashSet<(int, int)> CanonicalPairs(IEnumerable<(int Position1, int Position2)> pairs) =>
        pairs.Select(p => (System.Math.Min(p.Position1, p.Position2), System.Math.Max(p.Position1, p.Position2))).ToHashSet();

    #region RNA-DOTBRACKET-001 RT — parse then format is the identity

    [Test]
    [Description("RT: for a canonical nested dot-bracket string, formatting the parsed base pairs back reproduces the original notation exactly.")]
    public void DotBracket_ParseThenFormat_Identity()
    {
        const string notation = "((..))..((...))";
        var pairs = RnaSecondaryStructure.ParseDotBracket(notation).ToList();

        FormatDotBracket(notation.Length, pairs).Should().Be(notation,
            because: "the canonical formatter is the exact inverse of the parser for nested () notation");
    }

    #endregion

    #region RNA-DOTBRACKET-001 INV — pairing is preserved under reparse

    [Test]
    [Description("INV: re-parsing the reformatted string yields the same base-pair set, so the pairing survives a parse/format/reparse cycle.")]
    public void DotBracket_Reparse_PreservesPairing()
    {
        const string notation = "(((...)))..((..)).";
        var firstParse = RnaSecondaryStructure.ParseDotBracket(notation).ToList();

        string reformatted = FormatDotBracket(notation.Length, firstParse);
        var secondParse = RnaSecondaryStructure.ParseDotBracket(reformatted).ToList();

        CanonicalPairs(secondParse).Should().BeEquivalentTo(CanonicalPairs(firstParse),
            because: "formatting and re-parsing the same pairing is lossless");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-HAIRPIN-001 — hairpin loop free energy (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 150.
    //
    // API under test (RnaSecondaryStructure.CalculateHairpinLoopEnergy):
    //   ΔG = loop-initiation(size) + terminal-mismatch(closing pair, first/last loop base) +
    //   sequence bonuses/penalties (Turner 2004 / NNDB).
    //
    // Relations (derived from the Turner loop model, NOT from output):
    //   • MON  (larger loop ⇒ higher, less-stable energy): in the size-extrapolated regime
    //          (size ≥ 10, ΔG(n)=ΔG(9)+1.75·RT·ln(n/9)) the initiation term grows with loop size,
    //          so — holding the closing pair and terminal-mismatch bases fixed — energy increases.
    //   • INV  (closing-pair/terminal context determines energy): for size ≥ 4 the model reads only
    //          the size, the closing pair and the first/last loop bases, so changing the interior
    //          loop bases (no all-C / special-loop trigger) leaves the energy unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-HAIRPIN-001 MON — a larger loop raises the (less stable) energy

    [Test]
    [Description("MON: in the extrapolated size regime (≥10) the loop-initiation term grows with loop size, so with a fixed G·C closing pair and fixed A…A terminal mismatch the hairpin energy increases with loop length.")]
    public void Hairpin_LargerLoop_HigherEnergy()
    {
        double previous = double.MinValue;
        foreach (int size in new[] { 10, 15, 20, 30, 40 })
        {
            string loop = new string('A', size); // first/last base 'A' fixed; no all-C penalty
            double energy = RnaSecondaryStructure.CalculateHairpinLoopEnergy(loop, 'G', 'C');
            energy.Should().BeGreaterThan(previous, because: $"a size-{size} loop has a higher initiation energy than the smaller one");
            previous = energy;
        }
    }

    #endregion

    #region RNA-HAIRPIN-001 INV — interior loop bases do not change the energy

    [Test]
    [Description("INV: for size ≥ 4 the model uses only the size, the closing pair and the first/last loop bases, so changing the interior bases (no all-C/special trigger) leaves the energy unchanged.")]
    public void Hairpin_InteriorLoopBases_DoNotChangeEnergy()
    {
        double reference = RnaSecondaryStructure.CalculateHairpinLoopEnergy("ACCCA", 'G', 'C');

        foreach (var loop in new[] { "AGGGA", "AUUUA", "ACGUA" }) // same size, same first/last 'A'
            RnaSecondaryStructure.CalculateHairpinLoopEnergy(loop, 'G', 'C')
                .Should().BeApproximately(reference, 1e-9,
                    because: $"loop '{loop}' shares the closing pair, size and terminal mismatch, so only its (ignored) interior differs");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-INVERT-001 — inverted repeats (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 151.
    //
    // API under test (RnaSecondaryStructure.FindInvertedRepeats):
    //   Reports (Start1,End1,Start2,End2,Length) where the right arm is the antiparallel reverse
    //   complement of the left arm (the W…W̄ᴿ stem pattern), with a loop in [minSpacing, maxSpacing].
    //
    // Relations (derived from the reverse-complement stem definition, NOT from output):
    //   • SYM  (arms reverse-complementary): for every reported repeat the right-arm substring equals
    //          the reverse complement of the left-arm substring.
    //   • INV  (revcomp preserves count): reverse-complementing the whole sequence maps each inverted
    //          repeat to an inverted repeat, so the number of repeats is unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    private static string RnaRevComp(string s) =>
        new(s.Reverse().Select(RnaSecondaryStructure.GetComplement).ToArray());

    // Two well-separated RNA stem-loops (arm + loop + reverse-complement arm), U-spacer between them.
    private const string InvertedRepeatSeq = "GGGGAAAACCCC" + "UUUUU" + "GGACAAAGUCC";

    #region RNA-INVERT-001 SYM — the right arm is the reverse complement of the left

    [Test]
    [Description("SYM: each reported inverted repeat has a right arm equal to the reverse complement of its left arm.")]
    public void InvertedRepeats_RightArm_IsReverseComplementOfLeft()
    {
        var repeats = RnaSecondaryStructure.FindInvertedRepeats(InvertedRepeatSeq, minLength: 4, minSpacing: 3, maxSpacing: 8).ToList();
        repeats.Should().NotBeEmpty(because: "the sequence contains complementary stem arms");

        foreach (var r in repeats)
        {
            string left = InvertedRepeatSeq.Substring(r.Start1, r.End1 - r.Start1 + 1);
            string right = InvertedRepeatSeq.Substring(r.Start2, r.End2 - r.Start2 + 1);
            right.Should().Be(RnaRevComp(left),
                because: $"the right arm at [{r.Start2},{r.End2}] must be the reverse complement of the left arm '{left}'");
        }
    }

    #endregion

    #region RNA-INVERT-001 INV — reverse complement preserves the repeat count

    [Test]
    [Description("INV: reverse-complementing the whole sequence maps each inverted repeat to an inverted repeat, so the number of reported repeats is unchanged.")]
    public void InvertedRepeats_ReverseComplement_PreservesCount()
    {
        int original = RnaSecondaryStructure.FindInvertedRepeats(InvertedRepeatSeq, minLength: 4, minSpacing: 3, maxSpacing: 8).Count();
        int reversed = RnaSecondaryStructure.FindInvertedRepeats(RnaRevComp(InvertedRepeatSeq), minLength: 4, minSpacing: 3, maxSpacing: 8).Count();

        original.Should().BeGreaterThan(0, because: "the sequence has at least one stem-loop");
        reversed.Should().Be(original, because: "an inverted repeat stays an inverted repeat under whole-sequence reverse complement");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-MFE-001 — minimum free energy (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 152.
    //
    // API under test (RnaSecondaryStructure.CalculateMinimumFreeEnergy):
    //   Zuker-style DP with Turner 2004 nearest-neighbour parameters; returns ΔG in kcal/mol.
    //
    // Relations (derived from the energy model, NOT from output):
    //   • MON  (more GC pairs ⇒ lower MFE): GC stacks are more stabilising than AU, so lengthening a
    //          pure-GC stem (one more GC pair) lowers (makes more negative) the MFE.
    //   • INV  (U/T case-insensitive): T and U are the same base for folding, and the input is
    //          upper-cased, so the spelling (U vs T, upper vs lower) does not change the MFE.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-MFE-001 MON — more GC pairs lower the MFE

    [Test]
    [Description("MON: extending a pure-GC stem by one more GC pair adds a stabilising stack, so the minimum free energy strictly decreases.")]
    public void Mfe_MoreGcPairs_LowerEnergy()
    {
        double previous = double.MaxValue;
        foreach (int stem in new[] { 3, 4, 5, 6 })
        {
            // GC stem of the given length closing an AAA hairpin loop (loop bases cannot pair).
            string seq = new string('G', stem) + "AAA" + new string('C', stem);
            double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            mfe.Should().BeLessThan(previous, because: $"a {stem}-pair GC stem stacks more favourably than the shorter one");
            previous = mfe;
        }
    }

    #endregion

    #region RNA-MFE-001 INV — MFE is U/T- and case-insensitive

    [Test]
    [Description("INV: T and U denote the same base for folding and the input is upper-cased, so the U-spelled, T-spelled and lower-case forms of a sequence all give the same MFE.")]
    public void Mfe_Ut_And_Case_Insensitive()
    {
        const string rna = "GGGGGUUUCCCCC";
        double reference = RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna);

        RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna.Replace('U', 'T'))
            .Should().BeApproximately(reference, 1e-9, because: "T is read as U for folding");
        RnaSecondaryStructure.CalculateMinimumFreeEnergy(rna.ToLowerInvariant())
            .Should().BeApproximately(reference, 1e-9, because: "the input is upper-cased before folding");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-PAIR-001 — base-pair compatibility (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 153.
    //
    // API under test (RnaSecondaryStructure.CanPair):
    //   True when two bases can form a Watson–Crick (A–U, G–C) or wobble (G–U) pair.
    //
    // Relations (derived from the pairing definition, NOT from output):
    //   • SYM  (canPair(a,b)=canPair(b,a)): base pairing is a symmetric relation.
    //   • INV  (case-insensitive): pairing depends on the nucleotide identity, not letter case.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-PAIR-001 SYM — pairing is symmetric

    [Test]
    [Description("SYM: base pairing is a symmetric relation, so CanPair(a,b) equals CanPair(b,a) for every pair of bases.")]
    public void CanPair_Symmetric()
    {
        const string bases = "ACGU";
        foreach (char a in bases)
            foreach (char b in bases)
                RnaSecondaryStructure.CanPair(a, b).Should().Be(RnaSecondaryStructure.CanPair(b, a),
                    because: $"pairing of {a} and {b} cannot depend on the argument order");

        // Guard against a trivially-constant relation: the canonical pairs must actually pair.
        RnaSecondaryStructure.CanPair('A', 'U').Should().BeTrue(because: "A–U is Watson–Crick");
        RnaSecondaryStructure.CanPair('G', 'C').Should().BeTrue(because: "G–C is Watson–Crick");
        RnaSecondaryStructure.CanPair('G', 'U').Should().BeTrue(because: "G–U is a wobble pair");
        RnaSecondaryStructure.CanPair('A', 'G').Should().BeFalse(because: "A–G is not a valid pair");
    }

    #endregion

    #region RNA-PAIR-001 INV — pairing is case-insensitive

    [Test]
    [Description("INV: pairing depends on the nucleotide identity, not letter case, so lower-case bases pair exactly as their upper-case forms do.")]
    public void CanPair_CaseInsensitive()
    {
        const string bases = "ACGU";
        foreach (char a in bases)
            foreach (char b in bases)
                RnaSecondaryStructure.CanPair(char.ToLowerInvariant(a), char.ToLowerInvariant(b))
                    .Should().Be(RnaSecondaryStructure.CanPair(a, b),
                        because: $"the case of {a}/{b} must not affect whether they pair");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-PARTITION-001 — McCaskill partition function (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 154.
    //
    // API under test (RnaSecondaryStructure.CalculatePartitionFunction):
    //   Z = Σ over all pseudoknot-free structures of exp(−βE); with the simplified fixed-per-pair
    //   model every additional admissible pair adds Boltzmann-weighted structures to the ensemble.
    //
    // Relations (derived from the ensemble sum, NOT from output):
    //   • MON  (more pairing options ⇒ higher Z): a sequence with no possible pair has Z = 1 (only
    //          the empty structure); adding complementary content introduces extra structures, each
    //          with positive Boltzmann weight, so Z strictly increases.
    //   • INV  (deterministic): Z is a pure function of the sequence and parameters.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-PARTITION-001 MON — more pairing options raise Z

    [Test]
    [Description("MON: lengthening a hairpin stem adds admissible base pairs, each contributing extra Boltzmann-weighted structures, so the partition function Z strictly increases (from Z=1 for a non-pairing sequence).")]
    public void Partition_MorePairingOptions_HigherZ()
    {
        var sequences = new[]
        {
            "AAAAAA",          // no complementary pair possible ⇒ Z = 1
            "GAAAAC",          // one possible pair
            "GGAAAACC",        // nested pairs
            "GGGAAAACCC",      // more nested pairs
        };

        double previous = 0.0;
        bool first = true;
        foreach (var seq in sequences)
        {
            double z = RnaSecondaryStructure.CalculatePartitionFunction(seq).PartitionFunction;
            if (first) z.Should().Be(1.0, because: "a sequence with no admissible pair has only the empty structure");
            z.Should().BeGreaterThan(previous, because: $"'{seq}' admits more base-pair options, adding ensemble structures");
            previous = z;
            first = false;
        }
    }

    #endregion

    #region RNA-PARTITION-001 INV — the partition function is deterministic

    [Test]
    [Description("INV: Z is a pure function of the sequence, so repeated evaluations give the identical partition function.")]
    public void Partition_SameSequence_SameZ()
    {
        const string seq = "GGGAAAACCC";
        RnaSecondaryStructure.CalculatePartitionFunction(seq).PartitionFunction
            .Should().Be(RnaSecondaryStructure.CalculatePartitionFunction(seq).PartitionFunction,
                because: "the partition function has no hidden state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-PSEUDOKNOT-001 — pseudoknot detection (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 155.
    //
    // API under test (RnaSecondaryStructure.DetectPseudoknots):
    //   Reports a pseudoknot for every CROSSING pair-of-pairs (i < k < j < l). Nested
    //   (i < k < l < j) and disjoint (j < k) pairs are not pseudoknots.
    //
    // Relations (derived from the crossing condition, NOT from output):
    //   • INV  (nested ⇒ no pseudoknot): a purely nested set of base pairs has no crossing, so no
    //          pseudoknot is reported.
    //   • SHIFT (prepend flank shifts positions): adding a constant offset to every position (a
    //          prepended flank) preserves all order relations, so the same pseudoknots are reported
    //          with their coordinates shifted by the flank length.
    // ───────────────────────────────────────────────────────────────────────────

    private static RnaSecondaryStructure.BasePair Bp(int p1, int p2) =>
        new(p1, p2, 'A', 'U', RnaSecondaryStructure.BasePairType.WatsonCrick);

    #region RNA-PSEUDOKNOT-001 INV — a nested structure has no pseudoknot

    [Test]
    [Description("INV: a purely nested set of base pairs contains no crossing pairs, so no pseudoknot is detected.")]
    public void Pseudoknots_NestedStructure_None()
    {
        var nested = new[] { Bp(0, 9), Bp(1, 8), Bp(2, 7), Bp(3, 6) };

        RnaSecondaryStructure.DetectPseudoknots(nested).Should().BeEmpty(
            because: "nested pairs (i<k<l<j) never cross, so they form no pseudoknot");
    }

    #endregion

    #region RNA-PSEUDOKNOT-001 SHIFT — a prepended flank shifts pseudoknot coordinates

    [Test]
    [Description("SHIFT: offsetting every base-pair position by a flank length preserves the crossing relations, so the same pseudoknots are reported with coordinates shifted by the flank.")]
    public void Pseudoknots_PrependFlank_ShiftsCoordinates()
    {
        var crossing = new[] { Bp(0, 5), Bp(2, 8) }; // i=0 < k=2 < j=5 < l=8 ⇒ one pseudoknot
        var original = RnaSecondaryStructure.DetectPseudoknots(crossing)
            .Select(p => (p.Start1, p.End1, p.Start2, p.End2)).ToList();
        original.Should().ContainSingle(because: "the two pairs cross exactly once");

        foreach (int flank in new[] { 3, 100 })
        {
            var shifted = crossing.Select(b => Bp(b.Position1 + flank, b.Position2 + flank)).ToArray();
            var shiftedResult = RnaSecondaryStructure.DetectPseudoknots(shifted)
                .Select(p => (p.Start1, p.End1, p.Start2, p.End2)).ToList();

            var expected = original.Select(p => (p.Start1 + flank, p.End1 + flank, p.Start2 + flank, p.End2 + flank)).ToList();
            shiftedResult.Should().Equal(expected,
                because: $"a {flank}-base prepended flank shifts every pseudoknot coordinate by {flank} without changing the crossing");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-PKPREDICT-001 — canonical H-type pseudoknot prediction (RnaStructure / Analysis)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Reeder & Giegerich 2004, pknotsRG, BMC Bioinformatics 5:104;
    //   docs/algorithms/RnaStructure/Pseudoknot_Prediction.md §2, §6.1):
    //   PredictStructurePseudoknot folds an RNA into a structure that may contain a single
    //   canonical H-type pseudoknot — two crossing helices a·a' and b·b' with three loops —
    //   scored with the Turner 2004 nearest-neighbour model plus pknotsRG penalties (init 9.0,
    //   0.3 per unpaired pseudoknot-loop nt). A candidate knot is accepted ONLY when its ΔG is
    //   strictly below the plain pseudoknot-free MFE (INV-PK-01/04), and the predictor reads
    //   the input case-insensitively with T as U (§3.3).
    //
    // Two metamorphic relations (checklist row 236):
    //   • INV (known H-type knot recovered): a designed two-crossing-helix sequence is recovered
    //     as a pseudoknot, and that recovery is INVARIANT under representation-preserving input
    //     transforms — case folding (upper/lower/mixed) and the T↔U (DNA↔RNA) spelling — yielding
    //     an identical base-pair set, dot-bracket and ΔG. (The recovery itself, not a magic value,
    //     is the oracle-free property; the unit test pins the exact pairs separately.)
    //   • INV (no spurious knot on a plain hairpin): a simple hairpin is never reported as a
    //     pseudoknot — its structure and ΔG equal the plain MFE — and adding non-pairing 5'/3'
    //     context cannot fabricate a crossing helix (the 9 kcal/mol initiation penalty, INV-PK-04).
    //
    // API under test: RnaSecondaryStructure.PredictStructurePseudoknot / .CalculateMfeStructure
    //   / .DetectPseudoknots (PseudoknotStructure).

    #region RNA-PKPREDICT-001 — H-type pseudoknot prediction

    // Designed canonical H-type knot (two crossing 4-bp G·C helices, AA loops) — proven recovered
    // in the unit test; here it is the non-vacuous fixture for the encoding-invariance relations.
    private const string DesignedHTypeKnot = "GGGGAACCCCAACCCCAAGGGG";

    // Same geometry but with U in the (unpaired) loops, so the T↔U spelling transform is non-trivial.
    private const string DesignedHTypeKnotWithU = "GGGGUUCCCCUUCCCCUUGGGG";

    private static HashSet<(int, int)> PkPairSet(IEnumerable<(int Position1, int Position2)> pairs) =>
        pairs.Select(p => p.Position1 < p.Position2 ? (p.Position1, p.Position2) : (p.Position2, p.Position1)).ToHashSet();

    private static int PkCrossingCount(string seq, IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        var bps = pairs.Select(p => new RnaSecondaryStructure.BasePair(
            p.Position1, p.Position2, seq[p.Position1], seq[p.Position2],
            RnaSecondaryStructure.GetBasePairType(seq[p.Position1], seq[p.Position2])
                ?? RnaSecondaryStructure.BasePairType.NonCanonical)).ToList();
        return RnaSecondaryStructure.DetectPseudoknots(bps).Count();
    }

    [Test]
    [Description("INV: the predictor folds case-insensitively, so a known H-type knot is recovered identically (same pairs/dot-bracket/ΔG) regardless of upper/lower/mixed-case spelling.")]
    public void PkPredict_KnownHTypeKnot_RecoveredIdenticallyUnderCaseFolding()
    {
        var reference = RnaSecondaryStructure.PredictStructurePseudoknot(DesignedHTypeKnot);

        // Non-vacuity: this fixture really is a recovered, genuinely-crossing pseudoknot.
        reference.HasPseudoknot.Should().BeTrue(because: "the designed two-crossing-helix sequence must be folded as a pseudoknot");
        PkCrossingCount(reference.Sequence, reference.BasePairs).Should().BeGreaterThanOrEqualTo(1,
            because: "the recovered structure contains a genuine crossing pair (i<k<j<l)");

        foreach (string spelling in new[]
                 {
                     DesignedHTypeKnot.ToLowerInvariant(),
                     // mixed case: lower the even positions only
                     new string(DesignedHTypeKnot.Select((c, i) => i % 2 == 0 ? char.ToLowerInvariant(c) : c).ToArray()),
                 })
        {
            var pk = RnaSecondaryStructure.PredictStructurePseudoknot(spelling);

            pk.HasPseudoknot.Should().Be(reference.HasPseudoknot,
                because: $"case folding ('{spelling}') must not change the knot/no-knot decision");
            PkPairSet(pk.BasePairs).Should().BeEquivalentTo(PkPairSet(reference.BasePairs),
                because: "the predictor upper-cases its input, so the base-pair set is case-invariant");
            pk.DotBracket.Should().Be(reference.DotBracket, because: "the dot-bracket rendering is case-invariant");
            pk.FreeEnergy.Should().BeApproximately(reference.FreeEnergy, 1e-10, because: "ΔG is case-invariant");
        }
    }

    [Test]
    [Description("INV: T is read as U, so a known H-type knot spelled with T (DNA) folds identically to its U (RNA) spelling — same pairs, dot-bracket and ΔG.")]
    public void PkPredict_KnownHTypeKnot_RecoveredIdenticallyUnderTtoUSpelling()
    {
        string rna = DesignedHTypeKnotWithU;
        string dna = rna.Replace('U', 'T');

        var rnaPk = RnaSecondaryStructure.PredictStructurePseudoknot(rna);

        // Non-vacuity: the RNA spelling is a recovered, genuinely-crossing pseudoknot containing U.
        rna.Should().Contain("U", because: "the T↔U transform is only non-trivial if the RNA spelling contains U");
        rnaPk.HasPseudoknot.Should().BeTrue(because: "the U-loop H-type sequence must be folded as a pseudoknot");
        PkCrossingCount(rnaPk.Sequence, rnaPk.BasePairs).Should().BeGreaterThanOrEqualTo(1);

        var dnaPk = RnaSecondaryStructure.PredictStructurePseudoknot(dna);

        dnaPk.HasPseudoknot.Should().Be(rnaPk.HasPseudoknot, because: "T read as U must not change the knot decision");
        PkPairSet(dnaPk.BasePairs).Should().BeEquivalentTo(PkPairSet(rnaPk.BasePairs),
            because: "T and U pair identically (A–U ≡ A–T), so the base-pair set is spelling-invariant");
        dnaPk.DotBracket.Should().Be(rnaPk.DotBracket, because: "the dot-bracket is spelling-invariant");
        dnaPk.FreeEnergy.Should().BeApproximately(rnaPk.FreeEnergy, 1e-10, because: "ΔG is identical under T↔U");
    }

    [Test]
    [Description("INV: a simple hairpin is never reported as a pseudoknot; its returned structure and ΔG equal the plain pseudoknot-free MFE (the 9 kcal/mol initiation penalty forbids spurious knots).")]
    public void PkPredict_PlainHairpins_NoSpuriousKnot_EqualPlainMfe()
    {
        foreach (int arm in new[] { 4, 5, 6 })
        {
            string hairpin = GcHairpin(arm);
            var pk = RnaSecondaryStructure.PredictStructurePseudoknot(hairpin);
            var mfe = RnaSecondaryStructure.CalculateMfeStructure(hairpin);

            // Non-vacuity: the hairpin genuinely folds (it is a real stem-loop, just not a knot).
            mfe.BasePairs.Should().NotBeEmpty(because: $"a {arm}-bp G·C hairpin folds into a stem — the non-vacuity guard");

            pk.HasPseudoknot.Should().BeFalse(
                because: $"a plain hairpin (arm {arm}) cannot cross, and the 9 kcal/mol penalty forbids any spurious pseudoknot");
            pk.DotBracket.Should().Be(mfe.DotBracket,
                because: "with no accepted knot the returned structure is exactly the plain MFE structure");
            pk.FreeEnergy.Should().BeApproximately(mfe.FreeEnergy, 1e-10,
                because: "with no accepted knot the free energy equals the plain MFE");
        }
    }

    [Test]
    [Description("INV: adding non-pairing 5'/3' context to a plain hairpin cannot fabricate a crossing helix — no spurious pseudoknot is introduced.")]
    public void PkPredict_PlainHairpin_NonPairingContext_NeverFabricatesKnot()
    {
        string hairpin = GcHairpin(5); // GGGGG AAAA CCCCC — A pairs only with U/T, of which the flank has none

        RnaSecondaryStructure.PredictStructurePseudoknot(hairpin).HasPseudoknot.Should().BeFalse(
            because: "the bare hairpin is knot-free — the baseline for the context relation");

        foreach (int flank in new[] { 2, 5, 12 })
        {
            string a = new string('A', flank);
            foreach (string ctx in new[] { a + hairpin, hairpin + a, a + hairpin + a })
            {
                RnaSecondaryStructure.PredictStructurePseudoknot(ctx).HasPseudoknot.Should().BeFalse(
                    because: $"a poly-A flank cannot pair with the G/C/A hairpin, so it cannot create a crossing helix (ctx='{ctx}')");
            }
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-PKRECURSIVE-001 — recursive pknotsRG pseudoknot prediction (RnaStructure / Analysis)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Reeder & Giegerich 2004 / Reeder et al. 2007, pknotsRG;
    //   docs/algorithms/RnaStructure/Pseudoknot_Prediction_Recursive.md §2, §4):
    //   PredictStructurePseudoknotRecursive folds the WHOLE sequence with a memoised interval
    //   recurrence in which a pseudoknot "competes with values of unknotted foldings for the
    //   interval (i,j)", and a knot's three loops may themselves contain further knots. Its
    //   decomposition set includes the plain MFE and the single top-level H-type knot, so the
    //   recursive optimum can never be worse than either; it accepts a knot component only when
    //   it lowers ΔG (9 kcal/mol penalty, INV-PKR-01/04).
    //
    // Two metamorphic relations (checklist row 237):
    //   • MON (recursive ΔG ≤ single-knot ΔG): for every sequence the recursive folder's free
    //     energy dominates (≤) the single-knot folder's, which in turn dominates the plain MFE —
    //     because the recursive search strictly contains both alternatives as candidates. The
    //     dominance is strict on inputs the single-knot method cannot fold (over-arching nested
    //     and multiple knots), confirming the recursion genuinely helps.
    //   • INV (separable knots all recovered): a knot recovered in isolation is recovered
    //     UNCHANGED when a second, separable knot is present elsewhere, and the second knot is
    //     additionally recovered — so the multi-knot fold is the union of the local recoveries,
    //     with the crossing count additive. Recovery is local: a distant separable knot does not
    //     perturb it.
    //
    // API under test: RnaSecondaryStructure.PredictStructurePseudoknotRecursive vs
    //   .PredictStructurePseudoknot vs .CalculateMfeStructure (shared PseudoknotStructure).

    #region RNA-PKRECURSIVE-001 — recursive pseudoknot prediction

    // A single canonical H-type knot in an A·U clamp (knot nested inside an over-arching helix).
    private const string RecursiveNestedKnot = "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU";

    // Two A·U-clamped knots side by side (separable): both recovered by the recursive folder.
    private const string RecursiveTwoKnots =
        "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUUAAAAAAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU";

    [Test]
    [Description("MON: the recursive folder's decomposition set contains the single-knot fold and the plain MFE, so its ΔG dominates both (recursive ≤ single ≤ MFE); the dominance is strict where only recursion can fold the structure.")]
    public void PkRecursive_FreeEnergy_DominatesSingleKnotAndMfe()
    {
        string[] family =
        {
            DesignedHTypeKnot,          // single H-type
            RecursiveNestedKnot,        // over-arching nested knot (only recursion folds it)
            RecursiveTwoKnots,          // two separable knots (only recursion folds them)
            "GGGGAAAACCCC",             // plain hairpin
            GcHairpin(5),
            GcHairpin(6),
        };

        foreach (string seq in family)
        {
            double recursive = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq).FreeEnergy;
            double single = RnaSecondaryStructure.PredictStructurePseudoknot(seq).FreeEnergy;
            double mfe = RnaSecondaryStructure.CalculateMfeStructure(seq).FreeEnergy;

            recursive.Should().BeLessThanOrEqualTo(single + 1e-9,
                because: $"the recursive search includes the single-knot fold as a candidate, so it is never worse for '{seq}'");
            single.Should().BeLessThanOrEqualTo(mfe + 1e-9,
                because: $"the single-knot method falls back to the plain MFE, so it is never worse for '{seq}'");
        }

        // Strict dominance where the single-knot method cannot fold the structure (non-vacuity).
        foreach (string seq in new[] { RecursiveNestedKnot, RecursiveTwoKnots })
        {
            double recursive = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(seq).FreeEnergy;
            double single = RnaSecondaryStructure.PredictStructurePseudoknot(seq).FreeEnergy;
            RnaSecondaryStructure.PredictStructurePseudoknot(seq).HasPseudoknot.Should().BeFalse(
                because: $"the single-knot method cannot place a nested/second knot in '{seq}'");
            recursive.Should().BeLessThan(single - 1e-9,
                because: $"recursion recovers a knot the single-knot method cannot, so its ΔG is strictly lower for '{seq}'");
        }
    }

    [Test]
    [Description("INV: a knot recovered in isolation is recovered unchanged when a second separable knot is present, and the second knot is additionally recovered — the multi-knot fold is the union of the local recoveries (additive crossings).")]
    public void PkRecursive_SeparableKnots_AllRecovered_Compositionally()
    {
        var single = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(RecursiveNestedKnot);
        var both = RnaSecondaryStructure.PredictStructurePseudoknotRecursive(RecursiveTwoKnots);

        int singleCross = PkCrossingCount(single.Sequence, single.BasePairs);

        // Non-vacuity: the isolated knot really is a recovered, genuinely-crossing pseudoknot.
        single.HasPseudoknot.Should().BeTrue(because: "the A·U-clamped knot is recovered in isolation");
        singleCross.Should().BeGreaterThanOrEqualTo(1);

        both.HasPseudoknot.Should().BeTrue(because: "two separable knots must both be recovered");

        // (a) The first knot's pairs are recovered UNCHANGED inside the two-knot fold (local recovery).
        PkPairSet(both.BasePairs).IsSupersetOf(PkPairSet(single.BasePairs)).Should().BeTrue(
            because: "a distant separable knot does not perturb the first knot — its exact pairs persist");

        // (b) A second, independent crossing knot is recovered downstream of the first (indices > 37).
        var downstreamPairs = both.BasePairs.Where(p => p.Position1 > 37 && p.Position2 > 37).ToList();
        PkCrossingCount(both.Sequence, downstreamPairs).Should().BeGreaterThanOrEqualTo(1,
            because: "the second separable knot is recovered as its own crossing helix in the downstream region");

        // (c) Crossings are additive: two separable knots cross at least as much as two isolated ones.
        PkCrossingCount(both.Sequence, both.BasePairs).Should().BeGreaterThanOrEqualTo(2 * singleCross,
            because: "recovering two separable knots yields at least the sum of the two isolated knots' crossings");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-ACCESS-001 — McCaskill unpaired (accessibility) probabilities (RnaStructure)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (McCaskill 1990 Biopolymers 29:1105; Bernhart, Hofacker, Stadler 2006, RNAplfold;
    //   tests/TestSpecs/RNA-ACCESS-001.md; docs/Validation/reports/RNA-ACCESS-001.md):
    //   CalculateRegionUnpairedProbability returns the equilibrium probability that a contiguous
    //   window of L nucleotides ending at position e is ENTIRELY unpaired, computed exactly as
    //   Z_open(window) / Z from the Turner-2004 McCaskill partition function. Z is the full
    //   ensemble sum and is fixed for a sequence; Z_open(window) forbids every base in the window
    //   from pairing. CalculateUnpairedProbabilities returns the same Z, the per-base unpaired
    //   probabilities p_unpaired[i] and ΔG_ensemble = −RT·ln Z, from fixed (sequence-independent)
    //   thermodynamic constants (RT at 37 °C; Turner-2004 NN energies).
    //
    // Two metamorphic relations (checklist row 238):
    //   • MON (extending the queried region cannot raise its unpaired probability): enlarging the
    //     window (more bases required unpaired) is a stronger event, so its probability is
    //     non-increasing — Z_open shrinks as more positions are forbidden from pairing while Z is
    //     fixed, hence P(L+1) ≤ P(L) for a window ending at the same position.
    //   • INV (sequence-independent constants reproduce the analytic GAAAC value): because the
    //     constants are fixed and spelling-/case-agnostic, GAAAC (and CAAAAG) reproduce their
    //     hand-derived analytic ensemble (Z, p_unpaired, P, ΔG) exactly — invariantly under case
    //     folding and the T↔U spelling, and deterministically across runs.
    //
    // API under test: RnaSecondaryStructure.CalculateRegionUnpairedProbability /
    //   .CalculateUnpairedProbabilities (UnpairedProbabilityResult).

    #region RNA-ACCESS-001 — McCaskill accessibility probabilities

    [Test]
    [Description("MON: a longer unpaired window is a stronger (subset) event, so accessibility is non-increasing in window length for a window ending at the same position (Z_open shrinks, Z fixed).")]
    public void Access_RegionUnpairedProbability_NonIncreasingInWindowLength()
    {
        // Includes a strongly-pairing hairpin (genuine decrease) and a non-pairing homopolymer (flat at 1).
        string[] sequences = { "GCGCGCAAAAGCGCGC", "GGGGAAAACCCC", "GGGAAACCCAAAGGGAAACCC", "AAAAAAAAAA" };
        bool sawStrictDecrease = false;

        foreach (string seq in sequences)
        {
            foreach (int end in new[] { seq.Length - 1, seq.Length / 2 + 1 })
            {
                double previous = double.PositiveInfinity;
                for (int len = 1; len <= end + 1; len++)
                {
                    double p = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, len);

                    p.Should().BeInRange(-1e-12, 1 + 1e-12, because: "an accessibility is a probability in [0,1]");
                    p.Should().BeLessThanOrEqualTo(previous + 1e-9,
                        because: $"extending the window ending at {end} of '{seq}' to length {len} cannot raise its unpaired probability");
                    if (p < previous - 1e-6) sawStrictDecrease = true;
                    previous = p;
                }
            }
        }

        sawStrictDecrease.Should().BeTrue(
            because: "on a strongly-pairing sequence a longer window is strictly less likely to be fully unpaired — the relation is non-vacuous");
    }

    [Test]
    [Description("CONS: a length-1 region accessibility at i equals the per-base unpaired probability p_unpaired[i] — both are Z_forbid(i)/Z computed from the same partition function.")]
    public void Access_RegionLengthOne_EqualsPerBaseUnpairedProbability()
    {
        foreach (string seq in new[] { "GCGCGCAAAAGCGCGC", "GGGGAAAACCCC", "GGGAAACCCAAAGGGAAACCC" })
        {
            var perBase = RnaSecondaryStructure.CalculateUnpairedProbabilities(seq).UnpairedProbabilities;

            for (int i = 0; i < seq.Length; i++)
            {
                double region1 = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, i, 1);
                region1.Should().BeApproximately(perBase[i], 1e-9,
                    because: $"the length-1 accessibility at {i} of '{seq}' is exactly the per-base unpaired probability");
            }
        }
    }

    [Test]
    [Description("INV: from fixed sequence-independent constants, GAAAC and CAAAAG reproduce their analytic McCaskill ensemble exactly — invariantly under case folding and T↔U spelling, and deterministically.")]
    public void Access_AnalyticPins_ReproducedInvariantlyUnderSpelling()
    {
        // GAAAC: only the open chain or the single G(0)·C(4) hairpin (3-nt loop). Analytic pins.
        var gaaac = RnaSecondaryStructure.CalculateUnpairedProbabilities("GAAAC");
        gaaac.PartitionFunction.Should().BeApproximately(1.0001565052764922, 1e-12, because: "GAAAC analytic Z");
        gaaac.BasePairProbabilities[(0, 4)].Should().BeApproximately(0.00015648078642340854, 1e-12, because: "GAAAC analytic P(0,4)");
        gaaac.UnpairedProbabilities[0].Should().BeApproximately(0.9998435192135765, 1e-12, because: "GAAAC analytic p_unpaired(0)");
        gaaac.EnsembleFreeEnergy.Should().BeApproximately(-9.64416549414892e-05, 1e-12, because: "GAAAC analytic ΔG = −RT·ln Z");

        // The all-unpaired 5-window of GAAAC is the open-chain-only term = 1/Z.
        RnaSecondaryStructure.CalculateRegionUnpairedProbability("GAAAC", windowEnd: 4, windowLength: 5)
            .Should().BeApproximately(1.0 / gaaac.PartitionFunction, 1e-12,
                because: "forbidding all five bases from pairing leaves only the empty structure, weight 1 → 1/Z");

        // CAAAAG: 4-nt loop with a terminal mismatch. Second analytic pin.
        var caaaag = RnaSecondaryStructure.CalculateUnpairedProbabilities("CAAAAG");
        caaaag.PartitionFunction.Should().BeApproximately(1.0012902114608, 1e-12, because: "CAAAAG analytic Z");
        caaaag.BasePairProbabilities[(0, 5)].Should().BeApproximately(0.0012885489601637966, 1e-12, because: "CAAAAG analytic P(0,5)");
        caaaag.UnpairedProbabilities[0].Should().BeApproximately(0.9987114510398362, 1e-12, because: "CAAAAG analytic p_unpaired(0)");

        // INV under case folding: the constants are case-agnostic.
        var lower = RnaSecondaryStructure.CalculateUnpairedProbabilities("gaaac");
        lower.PartitionFunction.Should().BeApproximately(gaaac.PartitionFunction, 1e-15, because: "Z is case-invariant");
        lower.UnpairedProbabilities[0].Should().BeApproximately(gaaac.UnpairedProbabilities[0], 1e-15, because: "p_unpaired is case-invariant");

        // INV under T↔U spelling: a U-containing RNA and its DNA (T) spelling fold identically.
        const string rna = "GGGAAACCCUUU";
        var rnaRes = RnaSecondaryStructure.CalculateUnpairedProbabilities(rna);
        var dnaRes = RnaSecondaryStructure.CalculateUnpairedProbabilities(rna.Replace('U', 'T'));
        dnaRes.PartitionFunction.Should().BeApproximately(rnaRes.PartitionFunction, 1e-12, because: "T read as U → identical Z");
        for (int i = 0; i < rna.Length; i++)
            dnaRes.UnpairedProbabilities[i].Should().BeApproximately(rnaRes.UnpairedProbabilities[i], 1e-12,
                because: $"T read as U → identical p_unpaired at {i}");

        // INV determinism: repeated computation is bit-stable.
        var repeat = RnaSecondaryStructure.CalculateUnpairedProbabilities("GAAAC");
        repeat.PartitionFunction.Should().Be(gaaac.PartitionFunction, because: "the McCaskill computation is deterministic");
    }

    #endregion
}
