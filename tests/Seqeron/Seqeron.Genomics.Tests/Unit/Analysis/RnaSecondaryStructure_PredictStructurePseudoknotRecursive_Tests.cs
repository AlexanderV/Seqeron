// RNA-PKRECURSIVE-001 — Recursive pknotsRG grammar (multiple / nested canonical simple recursive pseudoknots)
// Evidence: docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PKRECURSIVE-001.md
// Sources: Reeder J, Giegerich R (2004) BMC Bioinformatics 5:104 (the loops u, v, w "fold internally in an arbitrary way,
//          including simple recursive pseudoknots"); Reeder, Steffen & Giegerich (2007) NAR 35:W320 (per-interval competition);
//          pknotsRG reference source (pknotsRG-mfe.lhs `block` / `comps` / `knot`; Algebras.lhs `pkml` = 6.0 for a knot used as
//          a multiloop or knot-loop component); ViennaRNA 2.7.2 primitives for the oracle rows of PknotsRgOracleData.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests
{
    private static HashSet<(int, int)> PairSet(IEnumerable<(int Position1, int Position2)> pairs)
        => pairs.Select(p => (Math.Min(p.Position1, p.Position2), Math.Max(p.Position1, p.Position2))).ToHashSet();

    private static void AssertValid(PseudoknotStructure s)
    {
        int n = s.Sequence.Length;
        var seen = new HashSet<int>();
        foreach (var (a, b) in s.BasePairs)
        {
            Assert.That(a, Is.InRange(0, n - 1));
            Assert.That(b, Is.InRange(0, n - 1));
            Assert.That(b - a - 1, Is.GreaterThanOrEqualTo(3));
            Assert.That(CanPair(s.Sequence[a], s.Sequence[b]), Is.True);
            Assert.That(seen.Add(a) && seen.Add(b), Is.True, "each position pairs at most once");
        }
        Assert.That(PairSet(ParseDotBracket(s.DotBracket)), Is.EquivalentTo(PairSet(s.BasePairs)), "dot-bracket round trip");
        var bps = s.BasePairs.Select(p => new BasePair(p.Position1, p.Position2, s.Sequence[p.Position1], s.Sequence[p.Position2],
            GetBasePairType(s.Sequence[p.Position1], s.Sequence[p.Position2]) ?? BasePairType.NonCanonical)).ToList();
        Assert.That(s.HasPseudoknot, Is.EqualTo(DetectPseudoknots(bps).Any()));
    }

    #region Oracle rows

    // Recursive optimum = independent ViennaRNA-primitive implementation of the full grammar (every row, 0.00).
    [TestCaseSource(typeof(PknotsRgOracleData), nameof(PknotsRgOracleData.RowSource))]
    public void Recursive_OracleRows_EnergyAndStructureLocked(
        string seq, double vienna, double single, double recursive, string singleDb, string recursiveDb)
    {
        var r = PredictStructurePseudoknotRecursive(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.FreeEnergy, Is.EqualTo(recursive).Within(0.005));
            Assert.That(r.DotBracket, Is.EqualTo(recursiveDb));
            Assert.That(r.FreeEnergy, Is.LessThanOrEqualTo(single + 1e-9));
            Assert.That(single, Is.LessThanOrEqualTo(vienna + 1e-9));
        });
        AssertValid(r);
    }

    #endregion

    #region Engineered grammar cases (oracle-verified)

    // Two knots side by side in the exterior loop (pknotsRG `struct`: components chained). The single-knot method places one.
    [Test]
    public void Recursive_TwoKnotsSideBySide_BothRecoveredAndBeatSingle()
    {
        const string seq = "CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC";
        var r = PredictStructurePseudoknotRecursive(seq);
        var s = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo("((((..[[[[[))))..]]]]]...((((.[[[[))))..]]]]"));
            Assert.That(r.FreeEnergy, Is.EqualTo(-22.38).Within(0.005));
            Assert.That(s.FreeEnergy, Is.EqualTo(-18.85).Within(0.005));
            Assert.That(r.DotBracket.Count(c => c == '['), Is.EqualTo(9), "5 + 4 crossing-helix pairs: two knots");
        });
    }

    // Designed twin H-type knots: one knot is the single-method optimum (−24.46), two knots the recursive optimum (−27.52).
    [Test]
    public void Recursive_TwinDesignedKnots_BothPlaced()
    {
        const string seq = "GGGGAACCCCAACCCCAAGGGGAAAGGGGAACCCCAACCCCAAGGGG";
        var r = PredictStructurePseudoknotRecursive(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo("((((..[[[[..))))..]]]]...((((..[[[[..))))..]]]]"));
            Assert.That(r.FreeEnergy, Is.EqualTo(-27.52).Within(0.005));
            Assert.That(PredictStructurePseudoknot(seq).FreeEnergy, Is.EqualTo(-24.46).Within(0.005));
        });
    }

    // A knot inside the loop v of another knot (the grammar's recursion: loops fold "including simple recursive pseudoknots").
    // Outer knot a=(0-2)/(a'=34..36 ...) with the inner knot ((((..[[[[..))))...]]]] between its b and a'.
    [Test]
    public void Recursive_KnotInsideKnotLoop_Recovered()
    {
        const string seq = "GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC";
        var r = PredictStructurePseudoknotRecursive(seq);
        var s = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo("(((.[[[((((..[[[[..))))...]]]])))..]]]"));
            Assert.That(r.FreeEnergy, Is.EqualTo(-12.15).Within(0.005));
            Assert.That(s.DotBracket, Is.EqualTo("(((.[[[......((((.((....)))))))))..]]]"), "single-knot: loop v folds knot-free");
            Assert.That(s.FreeEnergy, Is.EqualTo(-11.56).Within(0.005));
        });
        AssertValid(r);
    }

    // A knot as a branch of a multiloop (pknotsRG `pkml`: +6.0 in place of the branch term). Only the full grammar can place it.
    [Test]
    public void Recursive_KnotAsMultiloopBranch_Recovered()
    {
        const string seq = "CGCCCGGGCAAGGCGGGAAACCGCCAAGCGGGAGCCCCCCCCCGCCUGGGGGCAAGCCCGGGCG";
        var r = PredictStructurePseudoknotRecursive(seq);
        var s = PredictStructurePseudoknot(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo("(((((((((..(((((....)))))..(((((.[[[[[[.)))))..]]]]]]..)))))))))"));
            Assert.That(r.FreeEnergy, Is.EqualTo(-46.66).Within(0.005));
            Assert.That(s.HasPseudoknot, Is.False, "single-knot method has no knot in a multiloop");
            Assert.That(s.FreeEnergy, Is.EqualTo(-43.46).Within(0.005));
            Assert.That(CalculateMfeStructure(seq).FreeEnergy, Is.EqualTo(-43.46).Within(0.005));
        });
        AssertValid(r);
    }

    // A strong hairpin never becomes a knot.
    [Test]
    public void Recursive_PlainHairpin_NoSpuriousPseudoknot()
    {
        var r = PredictStructurePseudoknotRecursive("GGGGAAAACCCC");
        Assert.Multiple(() =>
        {
            Assert.That(r.HasPseudoknot, Is.False);
            Assert.That(r.DotBracket, Is.EqualTo("((((....))))"));
        });
    }

    #endregion

    #region Exhaustive equivalence and invariants

    // Below the size where a second / nested knot fits, the recursive optimum IS the single-knot optimum, which the exhaustive
    // reference computes independently (11–15 nt: a nested knot needs ≥ 2·11 nt, a multiloop knot ≥ 20 nt).
    [Test]
    public void Recursive_SmallSequences_EqualExhaustiveSingleKnotOptimum()
    {
        var rnd = new Random(99);
        int n = 0;
        while (n < 40)
        {
            string seq = n % 2 == 0 ? GenerateRandomRna(11 + rnd.Next(5), rnd, 0.7)
                : RnaSecondaryStructure_PredictStructurePseudoknot_Tests.RandomSmallKnotSequence(rnd);
            if (seq.Contains("GGUC")) continue;
            Assert.That(PredictStructurePseudoknotRecursive(seq).FreeEnergy,
                Is.EqualTo(new PknotsRgBruteForce(seq).SingleKnotOptimum()).Within(0.006), seq);
            n++;
        }
    }

    [Test]
    public void Recursive_RandomSweep_Invariants()
    {
        var rnd = new Random(12);
        for (int c = 0; c < 60; c++)
        {
            string seq = GenerateRandomRna(11 + rnd.Next(45), rnd, 0.65);
            var r = PredictStructurePseudoknotRecursive(seq);
            var s = PredictStructurePseudoknot(seq);
            Assert.That(r.FreeEnergy, Is.LessThanOrEqualTo(s.FreeEnergy + 1e-9), seq);
            Assert.That(s.FreeEnergy, Is.LessThanOrEqualTo(CalculateMfeStructure(seq).FreeEnergy + 1e-9), seq);
            AssertValid(r);
            if (!r.HasPseudoknot) Assert.That(r.DotBracket, Is.EqualTo(CalculateMfeStructure(seq).DotBracket), seq);
        }
    }

    // Recursion terminates and stays tractable: 150 nt (the previous interval-memo recursion needed 6.9 s for 80 nt).
    [Test, CancelAfter(120000)]
    public void Recursive_LongSequence_Completes()
    {
        string seq = GenerateRandomRna(150, new Random(5), 0.6);
        var r = PredictStructurePseudoknotRecursive(seq);
        Assert.That(r.FreeEnergy, Is.LessThanOrEqualTo(CalculateMfeStructure(seq).FreeEnergy + 1e-9));
        AssertValid(r);
    }

    [Test]
    public void Recursive_Concurrent_SameResult()
    {
        const string seq = "CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC";
        var reference = PredictStructurePseudoknotRecursive(seq);
        var results = new PseudoknotStructure[6];
        Parallel.For(0, results.Length, i => results[i] = PredictStructurePseudoknotRecursive(seq));
        Assert.That(results.Select(x => x.DotBracket), Is.All.EqualTo(reference.DotBracket));
    }

    #endregion

    #region Edge cases

    [TestCase(null)]
    [TestCase("")]
    public void Recursive_NullOrEmpty_EmptyStructure(string? seq)
    {
        var r = PredictStructurePseudoknotRecursive(seq!);
        Assert.Multiple(() =>
        {
            Assert.That(r.Sequence, Is.Empty);
            Assert.That(r.DotBracket, Is.Empty);
            Assert.That(r.BasePairs, Is.Empty);
            Assert.That(r.FreeEnergy, Is.EqualTo(0.0));
            Assert.That(r.HasPseudoknot, Is.False);
        });
    }

    [TestCase("GGGCCC")]
    [TestCase("GGGGAAAACC")]
    public void Recursive_TooShort_PlainMfe(string seq)
    {
        var r = PredictStructurePseudoknotRecursive(seq);
        Assert.Multiple(() =>
        {
            Assert.That(r.HasPseudoknot, Is.False);
            Assert.That(r.DotBracket, Is.EqualTo(CalculateMfeStructure(seq).DotBracket));
        });
    }

    [Test]
    public void Recursive_DnaAndLowerCase_ParityWithRna()
    {
        const string rna = "CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC";
        var a = PredictStructurePseudoknotRecursive(rna);
        var b = PredictStructurePseudoknotRecursive(rna.Replace('U', 'T').ToLowerInvariant());
        Assert.That(b.DotBracket, Is.EqualTo(a.DotBracket));
        Assert.That(b.FreeEnergy, Is.EqualTo(a.FreeEnergy));
        Assert.That(b.Sequence, Is.EqualTo(rna));
    }

    [Test]
    public void Recursive_MinLoopSizeBelowThree_ClampedToThree()
    {
        const string seq = "GGGGAACCCCAACCCCAAGGGG";
        Assert.That(PredictStructurePseudoknotRecursive(seq, 0).DotBracket, Is.EqualTo(PredictStructurePseudoknotRecursive(seq).DotBracket));
    }

    [Test]
    public void Recursive_Temperature_UsesSharedModel()
    {
        const string seq = "CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC";
        Assert.That(PredictStructurePseudoknotRecursive(seq, 3, 310.15).FreeEnergy, Is.EqualTo(PredictStructurePseudoknotRecursive(seq).FreeEnergy));
        Assert.That(PredictStructurePseudoknotRecursive(seq, 3, 298.15).FreeEnergy, Is.LessThan(PredictStructurePseudoknotRecursive(seq).FreeEnergy));
    }

    #endregion
}
