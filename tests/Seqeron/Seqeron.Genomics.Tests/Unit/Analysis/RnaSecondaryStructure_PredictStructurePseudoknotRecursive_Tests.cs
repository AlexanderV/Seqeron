// RNA-PKRECURSIVE-001 — Recursive pknotsRG pseudoknot prediction (nested / multiple H-type knots)
// Evidence: docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PKRECURSIVE-001.md
// Source: Reeder J, Giegerich R (2004). BMC Bioinformatics 5:104 (loops fold "including simple
//         recursive pseudoknots"; penalties 9.0 / 0.3 / 0.0 kcal/mol); Reeder J, Steffen P,
//         Giegerich R (2007). NAR 35:W320 (pseudoknot value competes with unknotted foldings per
//         interval); Antczak M et al. (2018). Bioinformatics 34(8):1304-1312 (crossing i<k<j<l).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests
{
    // Designed single canonical H-type (RNA-PKPREDICT-001 dataset):
    //   stem1 (0,15)(1,14)(2,13)(3,12); stem2 (6,21)(7,20)(8,19)(9,18); ΔG -8.76.
    private const string SingleHType = "GGGGAACCCCAACCCCAAGGGG";

    // Knot nested inside an over-arching A·U helix (Evidence dataset):
    //   clamp5'=A×8[0-7] · H-type[8-29] · clamp3'=U×8[30-37].
    //   outer helix (0,37)..(7,30); inner stem1 (8,23)..(11,20); inner crossing stem2 (14,29)..(17,26).
    private const string NestedKnot = "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU";

    // Two A·U-clamped knots side by side (Evidence dataset): both knots recovered, crossing-count 32.
    private const string TwoKnots =
        "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUUAAAAAAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU";

    private static HashSet<(int, int)> PairSet(IEnumerable<(int Position1, int Position2)> pairs)
        => pairs.Select(p => (Math.Min(p.Position1, p.Position2), Math.Max(p.Position1, p.Position2))).ToHashSet();

    private static int CrossingCount(string seq, IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        var bps = pairs.Select(p => new BasePair(
            p.Position1, p.Position2, seq[p.Position1], seq[p.Position2],
            GetBasePairType(seq[p.Position1], seq[p.Position2]) ?? BasePairType.NonCanonical)).ToList();
        return DetectPseudoknots(bps).Count();
    }

    #region Recursive recovery — knot nested inside an over-arching helix

    // M1 — Over-arching nested knot: outer A·U helix + inner crossing knot recovered exactly.
    [Test]
    public void Recursive_NestedKnot_RecoversOuterHelixAndInnerCrossingKnot()
    {
        var pk = PredictStructurePseudoknotRecursive(NestedKnot);

        var expected = new HashSet<(int, int)>
        {
            (0, 37), (1, 36), (2, 35), (3, 34), (4, 33), (5, 32), (6, 31), (7, 30), // outer A·U helix
            (8, 23), (9, 22), (10, 21), (11, 20),  // inner stem 1
            (14, 29), (15, 28), (16, 27), (17, 26), // inner crossing stem 2
        };

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.True,
                "A pseudoknot nested inside the outer helix loop must be predicted by the recursive folder.");
            Assert.That(PairSet(pk.BasePairs), Is.EquivalentTo(expected),
                "Outer A·U helix plus the inner crossing H-type knot must be recovered exactly (Evidence dataset).");
            Assert.That(pk.DotBracket, Is.EqualTo("((((((((((((..[[[[..))))..]]]]))))))))"),
                "Two-layer dot-bracket: () outer+stem1, [] crossing stem2.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(-14.37).Within(1e-2),
                "Recursive ΔG of the over-arching nested knot is -14.37 kcal/mol (derived from the energy model).");
        });
    }

    // M2 — The recursive method recovers a knot the single-knot method cannot (and beats its ΔG).
    [Test]
    public void Recursive_NestedKnot_BeatsSingleKnotMethodAndMfe()
    {
        var recursive = PredictStructurePseudoknotRecursive(NestedKnot);
        var single = PredictStructurePseudoknot(NestedKnot);
        var mfe = CalculateMfeStructure(NestedKnot);

        Assert.Multiple(() =>
        {
            Assert.That(single.HasPseudoknot, Is.False,
                "The single-knot method cannot combine the outer helix with an inner knot, so it finds no knot here.");
            Assert.That(single.FreeEnergy, Is.EqualTo(-13.05).Within(1e-2),
                "Single-knot method falls back to the pseudoknot-free fold (-13.05).");
            Assert.That(recursive.FreeEnergy, Is.LessThan(single.FreeEnergy - 1e-9),
                "The recursive method's over-arching knot is strictly more stable than the single-knot result.");
            Assert.That(recursive.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9),
                "INV-PKR-01: recursive ΔG never exceeds the plain MFE.");
        });
    }

    #endregion

    #region Recursive recovery — two separate (non-nested) pseudoknots

    // M3 — Two separate clamped knots are both recovered (two genuine crossings).
    [Test]
    public void Recursive_TwoSeparateKnots_RecoversBothCrossings()
    {
        var pk = PredictStructurePseudoknotRecursive(TwoKnots);

        var expected = new HashSet<(int, int)>
        {
            // first clamped knot
            (0, 37), (1, 36), (2, 35), (3, 34), (4, 33), (5, 32), (6, 31), (7, 30),
            (8, 23), (9, 22), (10, 21), (11, 20), (14, 29), (15, 28), (16, 27), (17, 26),
            // second clamped knot
            (38, 79), (39, 78), (40, 77), (41, 76), (42, 75), (43, 74), (44, 73), (45, 72),
            (50, 65), (51, 64), (52, 63), (53, 62), (56, 71), (57, 70), (58, 69), (59, 68),
        };

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.True, "Two side-by-side knots must be predicted.");
            Assert.That(PairSet(pk.BasePairs), Is.EquivalentTo(expected),
                "Both clamped H-type knots must be recovered exactly (Evidence dataset).");
            Assert.That(CrossingCount(pk.Sequence, pk.BasePairs), Is.EqualTo(32),
                "DetectPseudoknots reports 32 crossing pair-of-pairs (two genuine knots).");
            Assert.That(pk.FreeEnergy, Is.EqualTo(-28.74).Within(1e-2),
                "Recursive ΔG of the two-knot structure is -28.74 kcal/mol.");
        });
    }

    // M4 — The single-knot method recovers NEITHER knot; the recursive method beats it.
    [Test]
    public void Recursive_TwoSeparateKnots_BeatsSingleKnotMethod()
    {
        var recursive = PredictStructurePseudoknotRecursive(TwoKnots);
        var single = PredictStructurePseudoknot(TwoKnots);
        var mfe = CalculateMfeStructure(TwoKnots);

        Assert.Multiple(() =>
        {
            Assert.That(single.HasPseudoknot, Is.False,
                "The single-knot method recovers a pseudoknot-free fold (it cannot place two knots).");
            Assert.That(single.FreeEnergy, Is.EqualTo(-27.14).Within(1e-2),
                "Single-knot method ΔG equals the plain MFE (-27.14).");
            Assert.That(recursive.FreeEnergy, Is.LessThan(single.FreeEnergy - 1e-9),
                "The recursive two-knot structure is strictly more stable than the single-knot result.");
            Assert.That(recursive.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9),
                "INV-PKR-01 holds on the two-knot sequence.");
        });
    }

    #endregion

    #region No spurious pseudoknots + invariants

    // M5 — A plain hairpin must NOT be reported as a pseudoknot; result equals the plain MFE.
    [Test]
    public void Recursive_PlainHairpin_NoSpuriousPseudoknot()
    {
        const string hairpin = "GGGGAAAACCCC";
        var pk = PredictStructurePseudoknotRecursive(hairpin);
        var mfe = CalculateMfeStructure(hairpin);

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False,
                "The 9 kcal/mol initiation penalty prevents a spurious pseudoknot on a simple hairpin (INV-PKR-04).");
            Assert.That(pk.DotBracket, Is.EqualTo("((((....))))"), "Returned structure equals the plain MFE.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy).Within(1e-10), "ΔG equals the plain MFE (-5.28).");
            Assert.That(CrossingCount(pk.Sequence, pk.BasePairs), Is.EqualTo(0), "No crossing pairs.");
        });
    }

    // M6 — Random sweep: recursive ΔG ≤ MFE; any reported knot is valid and genuinely crosses.
    [Test]
    public void Recursive_RandomSweep_NeverWorseThanMfe_NoSpuriousKnots()
    {
        var rng = new Random(20260623); // fixed seed → deterministic
        Assert.Multiple(() =>
        {
            for (int t = 0; t < 150; t++)
            {
                string seq = GenerateRandomRna(rng.Next(12, 38), rng, 0.5);
                var pk = PredictStructurePseudoknotRecursive(seq);
                var mfe = CalculateMfeStructure(seq);

                Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9),
                    $"INV-PKR-01: never worse than the plain MFE. seq={seq}");

                if (pk.HasPseudoknot)
                {
                    Assert.That(pk.FreeEnergy, Is.LessThan(mfe.FreeEnergy),
                        $"a reported pseudoknot strictly beats the MFE. seq={seq}");
                    Assert.That(CrossingCount(pk.Sequence, pk.BasePairs), Is.GreaterThanOrEqualTo(1),
                        $"INV-PKR-02: a reported pseudoknot genuinely crosses. seq={seq}");

                    var seen = new HashSet<int>();
                    foreach (var (a, b) in pk.BasePairs)
                    {
                        Assert.That(seen.Add(a), Is.True, $"INV-PKR-03 position {a} once. seq={seq}");
                        Assert.That(seen.Add(b), Is.True, $"INV-PKR-03 position {b} once. seq={seq}");
                    }
                }
            }
        });
    }

    // M7 — Validity: indices in range, each position paired at most once (INV-PKR-03).
    [Test]
    public void Recursive_RecoveredStructures_AreValid()
    {
        Assert.Multiple(() =>
        {
            foreach (var seq in new[] { NestedKnot, TwoKnots, SingleHType })
            {
                var pk = PredictStructurePseudoknotRecursive(seq);
                int n = pk.Sequence.Length;
                var seen = new HashSet<int>();
                foreach (var (a, b) in pk.BasePairs)
                {
                    Assert.That(a, Is.InRange(0, n - 1), $"5' index in range. seq={seq}");
                    Assert.That(b, Is.InRange(0, n - 1), $"3' index in range. seq={seq}");
                    Assert.That(a, Is.Not.EqualTo(b), $"a base cannot pair itself. seq={seq}");
                    Assert.That(seen.Add(a), Is.True, $"position {a} paired once. seq={seq}");
                    Assert.That(seen.Add(b), Is.True, $"position {b} paired once. seq={seq}");
                }
            }
        });
    }

    // M8 — Genuine crossings: nested case ≥1, two-knot case ≥2 (INV-PKR-02).
    [Test]
    public void Recursive_KnotCases_HaveGenuineCrossings()
    {
        var nested = PredictStructurePseudoknotRecursive(NestedKnot);
        var two = PredictStructurePseudoknotRecursive(TwoKnots);

        Assert.Multiple(() =>
        {
            Assert.That(CrossingCount(nested.Sequence, nested.BasePairs), Is.GreaterThanOrEqualTo(1),
                "The nested case contains at least one genuine crossing.");
            Assert.That(CrossingCount(two.Sequence, two.BasePairs), Is.GreaterThanOrEqualTo(2),
                "The two-knot case contains at least two genuine crossings.");
        });
    }

    #endregion

    #region Edge cases and parity

    // S1 — null input → empty pseudoknot-free structure.
    [Test]
    public void Recursive_NullInput_EmptyStructure()
    {
        var pk = PredictStructurePseudoknotRecursive(null!);
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False, "null input has no pseudoknot.");
            Assert.That(pk.BasePairs, Is.Empty, "null input yields no base pairs.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(0.0).Within(1e-10), "null input ΔG is 0.");
            Assert.That(pk.Sequence, Is.EqualTo(string.Empty), "null input folds the empty sequence.");
        });
    }

    // S2 — empty input → empty structure.
    [Test]
    public void Recursive_EmptyInput_EmptyStructure()
    {
        var pk = PredictStructurePseudoknotRecursive("");
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False, "empty input has no pseudoknot.");
            Assert.That(pk.BasePairs, Is.Empty, "empty input yields no base pairs.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(0.0).Within(1e-10), "empty input ΔG is 0.");
        });
    }

    // S3 — too short for any canonical knot (< 11 nt) → no pseudoknot, equals plain MFE.
    [Test]
    public void Recursive_TooShort_NoPseudoknot()
    {
        const string shortSeq = "GGGCCC";
        var pk = PredictStructurePseudoknotRecursive(shortSeq);
        var mfe = CalculateMfeStructure(shortSeq);

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False, "Below 11 nt no canonical knot can form.");
            Assert.That(pk.DotBracket, Is.EqualTo(mfe.DotBracket), "Returned structure equals the plain MFE.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy).Within(1e-10), "ΔG equals the plain MFE.");
        });
    }

    // S4 — Single canonical H-type: recursion must not regress the single-knot case.
    [Test]
    public void Recursive_SingleHType_MatchesSingleKnotMethod()
    {
        var recursive = PredictStructurePseudoknotRecursive(SingleHType);
        var single = PredictStructurePseudoknot(SingleHType);

        Assert.Multiple(() =>
        {
            Assert.That(recursive.HasPseudoknot, Is.True, "Single H-type is still predicted.");
            Assert.That(recursive.DotBracket, Is.EqualTo(single.DotBracket),
                "Recursive method reproduces the single-knot dot-bracket on the single-knot case.");
            Assert.That(PairSet(recursive.BasePairs), Is.EquivalentTo(PairSet(single.BasePairs)),
                "Recursive method reproduces the single-knot base pairs.");
            Assert.That(recursive.FreeEnergy, Is.EqualTo(single.FreeEnergy).Within(1e-10),
                "Recursive method reproduces the single-knot ΔG (-8.76).");
        });
    }

    // C1 — DNA spelling (T) folds identically to the RNA spelling (T read as U).
    [Test]
    public void Recursive_DnaInput_FoldsIdenticallyToRna()
    {
        string dna = NestedKnot.Replace('U', 'T');
        var rnaPk = PredictStructurePseudoknotRecursive(NestedKnot);
        var dnaPk = PredictStructurePseudoknotRecursive(dna);

        Assert.Multiple(() =>
        {
            Assert.That(dnaPk.HasPseudoknot, Is.EqualTo(rnaPk.HasPseudoknot), "T read as U → same knot decision.");
            Assert.That(PairSet(dnaPk.BasePairs), Is.EquivalentTo(PairSet(rnaPk.BasePairs)),
                "T read as U → identical base pairs.");
            Assert.That(dnaPk.FreeEnergy, Is.EqualTo(rnaPk.FreeEnergy).Within(1e-10),
                "T read as U → identical free energy.");
        });
    }

    // C2 — minLoopSize below NNDB minimum is clamped to 3 (same result as default).
    [Test]
    public void Recursive_MinLoopSizeZero_ClampedToThree()
    {
        var clamped = PredictStructurePseudoknotRecursive(NestedKnot, minLoopSize: 0);
        var def = PredictStructurePseudoknotRecursive(NestedKnot, minLoopSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(clamped.HasPseudoknot, Is.EqualTo(def.HasPseudoknot), "minLoopSize<3 clamps to 3.");
            Assert.That(clamped.DotBracket, Is.EqualTo(def.DotBracket), "clamped result equals the default.");
            Assert.That(clamped.FreeEnergy, Is.EqualTo(def.FreeEnergy).Within(1e-10), "clamped ΔG equals default.");
        });
    }

    #endregion
}
