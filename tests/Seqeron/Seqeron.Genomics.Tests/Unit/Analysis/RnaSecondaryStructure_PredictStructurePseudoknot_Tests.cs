// RNA-PKPREDICT-001 — Pseudoknot structure prediction, single canonical H-type knot (pknotsRG class csr-PK)
// Evidence: docs/Evidence/RNA-PKPREDICT-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PKPREDICT-001.md
// Sources: Reeder J, Giegerich R (2004) BMC Bioinformatics 5:104 (canonical simple recursive pseudoknots, rules 1-3);
//          pknotsRG reference source (github.com/jensreeder/pknotsRG: Energy.lhs pkinit 9.0 / npp 0.3, Algebras.lhs pkml 6.0,
//          pknotsRG-mfe.lhs `pknot` grammar: bridging bases, v >= 0, rule 3 `getbeta`);
//          ViennaRNA 2.7.2 (E_Hairpin / E_IntLoop / E_MLstem / E_ExtLoop primitives, NNDB-precision parameter file) for the
//          oracle values of PknotsRgOracleData.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_PredictStructurePseudoknot_Tests
{
    // Designed canonical H-type pseudoknot: a=GGGG[0-3] u=AA b=CCCC[6-9] v=AA a'=CCCC[12-15] w=AA b'=GGGG[18-21].
    private const string DesignedHTypeSeq = "GGGGAACCCCAACCCCAAGGGG";

    private static HashSet<(int, int)> PairSet(IEnumerable<(int Position1, int Position2)> pairs)
        => pairs.Select(p => (Math.Min(p.Position1, p.Position2), Math.Max(p.Position1, p.Position2))).ToHashSet();

    private static int CrossingCount(string seq, IReadOnlyList<(int Position1, int Position2)> pairs)
        => DetectPseudoknots(pairs.Select(p => new BasePair(
            p.Position1, p.Position2, seq[p.Position1], seq[p.Position2],
            GetBasePairType(seq[p.Position1], seq[p.Position2]) ?? BasePairType.NonCanonical)).ToList()).Count();

    private static void AssertValidStructure(PseudoknotStructure s)
    {
        int n = s.Sequence.Length;
        var seen = new HashSet<int>();
        foreach (var (a, b) in s.BasePairs)
        {
            Assert.That(a, Is.InRange(0, n - 1));
            Assert.That(b, Is.InRange(0, n - 1));
            Assert.That(a, Is.LessThan(b), "pairs are 5' < 3'");
            Assert.That(b - a - 1, Is.GreaterThanOrEqualTo(3), "hairpin loop >= 3 nt");
            Assert.That(CanPair(s.Sequence[a], s.Sequence[b]), Is.True, $"pair ({a},{b}) must be canonical/wobble");
            Assert.That(seen.Add(a) && seen.Add(b), Is.True, "every position pairs at most once");
        }
        Assert.That(s.BasePairs.Select(p => p.Position1), Is.Ordered, "sorted by 5' position");
        Assert.That(s.DotBracket, Has.Length.EqualTo(n));
        Assert.That(PairSet(ParseDotBracket(s.DotBracket)), Is.EquivalentTo(PairSet(s.BasePairs)),
            "dot-bracket round trip (ToDotBracket / ParseDotBracket)");
        Assert.That(s.HasPseudoknot, Is.EqualTo(CrossingCount(s.Sequence, s.BasePairs) > 0),
            "HasPseudoknot <=> the pair set has a crossing pair-of-pairs");
    }

    #region M — designed knots, energies locked to the independent ViennaRNA-primitive oracle

    // M1 — Designed H-type: both crossing helices recovered exactly; ΔG = oracle (−13.56), the plain MFE is −7.84.
    [Test]
    public void PredictStructurePseudoknot_DesignedHType_RecoversBothHelices_EnergyLocked()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);
        var expected = new HashSet<(int, int)>
        {
            (0, 15), (1, 14), (2, 13), (3, 12), (6, 21), (7, 20), (8, 19), (9, 18),
        };
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.True);
            Assert.That(PairSet(pk.BasePairs), Is.EquivalentTo(expected));
            Assert.That(pk.DotBracket, Is.EqualTo("((((..[[[[..))))..]]]]"));
            Assert.That(pk.FreeEnergy, Is.EqualTo(-13.56).Within(0.005),
                "9.0 + 0.9 (3 bridging nt) + 6×(−3.26) stacks + 0.3 (u') + 0.6 (v) − 4.80 (four helix ends) = −13.56");
            Assert.That(CalculateMfeStructure(DesignedHTypeSeq).FreeEnergy, Is.EqualTo(-7.84).Within(0.005));
        });
    }

    // M2 — a knot is reported only when strictly better than the pseudoknot-free optimum.
    [Test]
    public void PredictStructurePseudoknot_Knot_StrictlyBeatsPlainMfe()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);
        Assert.That(pk.FreeEnergy, Is.LessThan(CalculateMfeStructure(DesignedHTypeSeq).FreeEnergy - 1.0));
    }

    // M3 — DEFECT REPRO (before: flanks ignored): a hairpin 5' of the knot is folded as well.
    // Before the fix `GGGCGAAAGCCC AA <knot>` returned the plain MFE (−15.98, no knot): the knot-only candidate (−8.76, flanks left
    // unfolded and unscored) never beat the whole-sequence MFE. After: hairpin + knot, −21.70 (oracle).
    [Test]
    public void PredictStructurePseudoknot_KnotWithFoldedFlank_FoldsFlanksToo()
    {
        const string seq = "GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG";
        var pk = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.True);
            Assert.That(pk.DotBracket, Is.EqualTo("((((....))))..((((..[[[[..))))..]]]]"));
            Assert.That(pk.FreeEnergy, Is.EqualTo(-21.70).Within(0.005));
            Assert.That(CalculateMfeStructure(seq).FreeEnergy, Is.EqualTo(-15.98).Within(0.005));
            Assert.That(pk.BasePairs, Does.Contain((0, 11)), "the 5' hairpin stem closes the flank");
        });
    }

    // M4 — DEFECT REPRO (before: loops v ≥ 1 and overlap candidates dropped): knot with an EMPTY middle loop v
    // (pknotsRG `emptymid`, rule 3: b is shortened where it would compete with a'). Before: no knot, −4.22. After: −5.93 (oracle).
    [Test]
    public void PredictStructurePseudoknot_EmptyMiddleLoop_KnotFound()
    {
        const string seq = "UCAUGUUAGCUCAUGAAUGGCU";
        var pk = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(pk.DotBracket, Is.EqualTo("(((((..[[[[)))))..]]]]"));
            Assert.That(pk.FreeEnergy, Is.EqualTo(-5.93).Within(0.005));
            Assert.That(CalculateMfeStructure(seq).FreeEnergy, Is.EqualTo(-4.22).Within(0.005));
        });
    }

    // M5 — pknotsRG README example (sequence of the reference distribution) and BWYV frameshift pseudoknot (PDB 437D sequence):
    // the model optimum contains a knot; ΔG locked to the oracle (the crystal structure is tertiary-stabilised and is NOT claimed to match).
    [TestCase("UCAAGUAUUCCGAAGCUCAACGGGAAAAUGAGCUA", -6.54, -12.31, ".......(((((.[[[[[[.)))))...]]]]]].")]
    [TestCase("GGCGCGGCACCGUCCGCGGAACAAACGG", -7.90, -10.51, "..(((((..[[[[)))))......]]]]")]
    public void PredictStructurePseudoknot_PublishedSequences_LockedToOracle(string seq, double mfe, double single, string db)
    {
        var pk = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(CalculateMfeStructure(seq).FreeEnergy, Is.EqualTo(mfe).Within(0.005));
            Assert.That(pk.FreeEnergy, Is.EqualTo(single).Within(0.005));
            Assert.That(pk.DotBracket, Is.EqualTo(db));
        });
    }

    #endregion

    #region Oracle lock table (independent ViennaRNA-primitive implementation of the grammar)

    [TestCaseSource(typeof(PknotsRgOracleData), nameof(PknotsRgOracleData.RowSource))]
    public void PredictStructurePseudoknot_OracleRows_EnergyAndStructureLocked(
        string seq, double vienna, double single, double recursive, string singleDb, string recursiveDb)
    {
        Assert.That(recursive, Is.LessThanOrEqualTo(single + 1e-9), "oracle sanity: the recursive optimum cannot be worse");
        var pk = PredictStructurePseudoknot(seq);
        var mfe = CalculateMfeStructure(seq);
        Assert.Multiple(() =>
        {
            Assert.That(mfe.FreeEnergy, Is.EqualTo(vienna).Within(0.005), "plain MFE = RNA.fold (NNDB-precision parameters)");
            Assert.That(pk.FreeEnergy, Is.EqualTo(single).Within(0.005), "single-knot optimum = oracle");
            Assert.That(pk.DotBracket, Is.EqualTo(singleDb));
            if (!pk.HasPseudoknot)
            {
                Assert.That(pk.DotBracket, Is.EqualTo(mfe.DotBracket), "no knot => exactly the plain MFE structure");
                Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy));
            }
        });
        AssertValidStructure(pk);
    }

    #endregion

    #region Exhaustive reference (table-free)

    // Random H-type-biased sequence of 11–16 nt: a·u·b·v·a'·w·b' with 2–3 nt helices (G/C rich) and short loops.
    internal static string RandomSmallKnotSequence(Random rnd)
    {
        static string Rc(string s) => new(s.Reverse().Select(c => c switch { 'A' => 'U', 'U' => 'A', 'G' => 'C', _ => 'G' }).ToArray());
        string Helix(int len) => new(Enumerable.Range(0, len).Select(_ => "GGCCGCAU"[rnd.Next(8)]).ToArray());
        string Loop(int len) => new(Enumerable.Range(0, len).Select(_ => "AAUAGC"[rnd.Next(6)]).ToArray());
        while (true)
        {
            string a = Helix(2 + rnd.Next(2)), b = Helix(2 + rnd.Next(2));
            string seq = a + Loop(1 + rnd.Next(2)) + b + Loop(rnd.Next(3)) + Rc(a) + Loop(2 + rnd.Next(2)) + Rc(b);
            if (seq.Length <= 16 && !seq.Contains("GGUC")) return seq;
        }
    }

    // The DP equals the exhaustive optimum of the same search space (all flank structures × all canonical knots × all loop structures).
    [Test]
    public void PredictStructurePseudoknot_RandomSequences_EqualsExhaustiveOptimum()
    {
        var rnd = new Random(20260611);
        int knots = 0;
        for (int c = 0; c < 70; c++)
        {
            string seq = c % 3 == 0
                ? GenerateRandomRna(11 + rnd.Next(5), rnd, 0.7)
                : RandomSmallKnotSequence(rnd);
            if (seq.Contains("GGUC")) continue; // the special GGUC/CUGG three-stack of a knot helix is not enumerated here
            double brute = new PknotsRgBruteForce(seq).SingleKnotOptimum();
            var pk = PredictStructurePseudoknot(seq);
            Assert.That(pk.FreeEnergy, Is.EqualTo(brute).Within(0.006), seq);
            if (pk.HasPseudoknot) knots++;
        }
        Assert.That(knots, Is.GreaterThan(5), "the sweep must contain knotted optima");
    }

    // Knot-rich designed inputs (long helices, short loops) against the exhaustive optimum.
    [TestCase("GGGAACCCAACCCAAGGG")]
    [TestCase("GCGAGCCGCAGCGCUCGCCUA")]
    [TestCase("CCGGAUCCGGAUCCGGAUC")]
    [TestCase("GGAGCCUCCAGGCUCC")]
    public void PredictStructurePseudoknot_DesignedSmall_EqualsExhaustiveOptimum(string seq)
    {
        double brute = new PknotsRgBruteForce(seq).SingleKnotOptimum();
        Assert.That(PredictStructurePseudoknot(seq).FreeEnergy, Is.EqualTo(brute).Within(0.006));
    }

    #endregion

    #region Invariants

    // INV: ΔG ≤ plain MFE; no knot ⇒ plain MFE; structure valid; deterministic; recursive ≤ single.
    [Test]
    public void PredictStructurePseudoknot_RandomSweep_Invariants()
    {
        var rnd = new Random(7);
        for (int c = 0; c < 80; c++)
        {
            string seq = GenerateRandomRna(11 + rnd.Next(30), rnd, 0.6);
            var pk = PredictStructurePseudoknot(seq);
            var mfe = CalculateMfeStructure(seq);
            Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9), seq);
            if (!pk.HasPseudoknot) Assert.That(pk.DotBracket, Is.EqualTo(mfe.DotBracket), seq);
            else Assert.That(pk.FreeEnergy, Is.LessThan(mfe.FreeEnergy), seq);
            AssertValidStructure(pk);
            Assert.That(PredictStructurePseudoknotRecursive(seq).FreeEnergy, Is.LessThanOrEqualTo(pk.FreeEnergy + 1e-9), seq);
        }
    }

    // At most ONE knot: the crossing pairs of the result form a single H-type (two helices a/b).
    [Test]
    public void PredictStructurePseudoknot_AtMostOneKnot_EvenWhenTwoWouldFit()
    {
        const string knot = "GGGGAACCCCAACCCCAAGGGG";
        string seq = knot + "AAA" + knot;
        var single = PredictStructurePseudoknot(seq);
        var recursive = PredictStructurePseudoknotRecursive(seq);
        Assert.Multiple(() =>
        {
            Assert.That(single.HasPseudoknot, Is.True);
            Assert.That(single.DotBracket.Count(c => c == '['), Is.EqualTo(4), "one knot: exactly one crossing helix");
            Assert.That(recursive.DotBracket.Count(c => c == '['), Is.EqualTo(8), "two knots side by side");
            Assert.That(recursive.FreeEnergy, Is.LessThan(single.FreeEnergy - 2.0));
        });
    }

    // Spelling: DNA T, lower case and RNA U give identical results.
    [Test]
    public void PredictStructurePseudoknot_DnaAndLowerCase_ParityWithRna()
    {
        var rna = PredictStructurePseudoknot("UCAUGUUAGCUCAUGAAUGGCU");
        foreach (string variant in new[] { "TCATGTTAGCTCATGAATGGCT", "ucauguuagcucaugaauggcu", "tcatgttagctcatgaatggct" })
        {
            var other = PredictStructurePseudoknot(variant);
            Assert.That(other.FreeEnergy, Is.EqualTo(rna.FreeEnergy));
            Assert.That(other.DotBracket, Is.EqualTo(rna.DotBracket));
            Assert.That(other.Sequence, Is.EqualTo("UCAUGUUAGCUCAUGAAUGGCU"));
        }
    }

    // minLoopSize < 3 is clamped to the NNDB minimum 3.
    [Test]
    public void PredictStructurePseudoknot_MinLoopSizeBelowThree_ClampedToThree()
    {
        var a = PredictStructurePseudoknot(DesignedHTypeSeq, 0);
        var b = PredictStructurePseudoknot(DesignedHTypeSeq);
        Assert.That(a.DotBracket, Is.EqualTo(b.DotBracket));
        Assert.That(a.FreeEnergy, Is.EqualTo(b.FreeEnergy));
    }

    // Temperature: the Turner terms follow the shared TurnerModel; 37 °C = default; the knot energy changes with T.
    [Test]
    public void PredictStructurePseudoknot_Temperature_UsesSharedModel()
    {
        var at37 = PredictStructurePseudoknot(DesignedHTypeSeq, 3, 310.15);
        var at25 = PredictStructurePseudoknot(DesignedHTypeSeq, 3, 298.15);
        Assert.Multiple(() =>
        {
            Assert.That(at37.FreeEnergy, Is.EqualTo(PredictStructurePseudoknot(DesignedHTypeSeq).FreeEnergy));
            Assert.That(at25.FreeEnergy, Is.LessThan(at37.FreeEnergy), "colder => more stable stacks");
            Assert.That(at25.HasPseudoknot, Is.True);
        });
    }

    // Thread safety: static pure function; concurrent identical calls agree.
    [Test]
    public void PredictStructurePseudoknot_Concurrent_SameResult()
    {
        var reference = PredictStructurePseudoknot("GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC");
        var results = new PseudoknotStructure[8];
        Parallel.For(0, results.Length, i => results[i] = PredictStructurePseudoknot("GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC"));
        foreach (var r in results)
        {
            Assert.That(r.DotBracket, Is.EqualTo(reference.DotBracket));
            Assert.That(r.FreeEnergy, Is.EqualTo(reference.FreeEnergy));
        }
    }

    // Long input stays tractable (O(n^4) worst case, near n^3 in practice): 150 nt completes.
    [Test, CancelAfter(120000)]
    public void PredictStructurePseudoknot_LongSequence_Completes()
    {
        string seq = GenerateRandomRna(150, new Random(3), 0.55);
        var pk = PredictStructurePseudoknot(seq);
        Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(CalculateMfeStructure(seq).FreeEnergy + 1e-9));
        AssertValidStructure(pk);
    }

    #endregion

    #region Edge cases

    [TestCase(null)]
    [TestCase("")]
    public void PredictStructurePseudoknot_NullOrEmpty_EmptyStructure(string? seq)
    {
        var pk = PredictStructurePseudoknot(seq!);
        Assert.Multiple(() =>
        {
            Assert.That(pk.Sequence, Is.Empty);
            Assert.That(pk.DotBracket, Is.Empty);
            Assert.That(pk.BasePairs, Is.Empty);
            Assert.That(pk.FreeEnergy, Is.EqualTo(0.0));
            Assert.That(pk.HasPseudoknot, Is.False);
        });
    }

    // < 11 nt: no canonical knot fits; the result is the plain MFE.
    [TestCase("GGGCCC")]
    [TestCase("GGGGAAAACC")]
    [TestCase("GGGAAACCCA")]
    public void PredictStructurePseudoknot_TooShortForKnot_PlainMfe(string seq)
    {
        var pk = PredictStructurePseudoknot(seq);
        var mfe = CalculateMfeStructure(seq);
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False);
            Assert.That(pk.DotBracket, Is.EqualTo(mfe.DotBracket));
            Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy));
        });
    }

    // A strong hairpin never becomes a pseudoknot (no spurious knot; 9.0 initiation).
    [Test]
    public void PredictStructurePseudoknot_Hairpin_NoSpuriousKnot()
    {
        var pk = PredictStructurePseudoknot("GGGGAAAACCCC");
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False);
            Assert.That(pk.DotBracket, Is.EqualTo("((((....))))"));
            Assert.That(pk.FreeEnergy, Is.EqualTo(-5.28).Within(0.005));
        });
    }

    // N / IUPAC symbols never pair; the result stays valid and ΔG ≤ MFE.
    [TestCase("GGGGAANCCCCAACCCCAAGGGG")]
    [TestCase("GGGGAACCCCRAACCCCAAGGGG")]
    [TestCase("NNNNNNNNNNNNNNNNNNNN")]
    public void PredictStructurePseudoknot_AmbiguousBases_ValidAndNeverWorseThanMfe(string seq)
    {
        var pk = PredictStructurePseudoknot(seq);
        AssertValidStructure(pk);
        Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(CalculateMfeStructure(seq).FreeEnergy + 1e-9));
    }

    #endregion
}
