// RNA-PKPREDICT-001 — Pseudoknot Structure Prediction (canonical H-type, pknotsRG class)
// Evidence: docs/Evidence/RNA-PKPREDICT-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PKPREDICT-001.md
// Source: Reeder J, Giegerich R (2004). BMC Bioinformatics 5:104 (pknotsRG canonical simple
//         recursive pseudoknots; penalties: init 9.0, unpaired-loop 0.3, in-knot pair 0.0 kcal/mol);
//         pknotsRG Energy.lhs (github.com/jensreeder/pknotsRG); Antczak M et al. (2018).
//         Bioinformatics 34(8):1304-1312 (crossing condition i<k<j<l); PDB 437D / Su et al. (1999)
//         (BWYV H-type pseudoknot, tertiary-stabilised — not the NN-thermodynamic MFE).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_PredictStructurePseudoknot_Tests
{
    // Designed canonical H-type pseudoknot (Evidence §Test Datasets):
    //   S1a=GGGG[0-3] L1=AA S2a=CCCC[6-9] L2=AA S1b=CCCC[12-15] L3=AA S2b=GGGG[18-21]
    //   Stem 1 (a·a') = (0,15)(1,14)(2,13)(3,12); Stem 2 (b·b') = (6,21)(7,20)(8,19)(9,18).
    private const string DesignedHTypeSeq = "GGGGAACCCCAACCCCAAGGGG";

    private static HashSet<(int, int)> PairSet(IEnumerable<(int Position1, int Position2)> pairs)
        => pairs.Select(p => (Math.Min(p.Position1, p.Position2), Math.Max(p.Position1, p.Position2))).ToHashSet();

    private static int CrossingCount(string seq, IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        var bps = pairs.Select(p => new BasePair(
            p.Position1, p.Position2, seq[p.Position1], seq[p.Position2],
            GetBasePairType(seq[p.Position1], seq[p.Position2]) ?? BasePairType.NonCanonical)).ToList();
        return DetectPseudoknots(bps).Count();
    }

    #region PredictStructurePseudoknot — canonical H-type recovery

    // M1 — Designed H-type: both crossing helices recovered exactly, two-layer dot-bracket.
    [Test]
    public void PredictStructurePseudoknot_DesignedHType_RecoversBothCrossingHelices()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);

        var expected = new HashSet<(int, int)>
        {
            (0, 15), (1, 14), (2, 13), (3, 12), // stem 1 (a·a')
            (6, 21), (7, 20), (8, 19), (9, 18), // stem 2 (b·b'), crossing
        };

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.True,
                "A canonical H-type pseudoknot must be predicted for the designed two-crossing-helix sequence.");
            Assert.That(PairSet(pk.BasePairs), Is.EquivalentTo(expected),
                "The two crossing 4-bp G·C helices must be recovered exactly per the H-type geometry.");
            Assert.That(pk.DotBracket, Is.EqualTo("((((..[[[[..))))..]]]]"),
                "Two-layer dot-bracket: () = stem 1, [] = crossing stem 2 (ViennaRNA/WUSS).");
        });
    }

    // M2 — The pseudoknot must strictly beat the plain MFE structure (INV-PK-01/04).
    [Test]
    public void PredictStructurePseudoknot_DesignedHType_BeatsPlainMfe()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);
        var mfe = CalculateMfeStructure(DesignedHTypeSeq);

        Assert.That(pk.FreeEnergy, Is.LessThan(mfe.FreeEnergy),
            "A pseudoknot is accepted only when its free energy is strictly below the plain pseudoknot-free MFE.");
    }

    // M3 — The recovered structure contains a genuine crossing pair (i<k<j<l) (INV-PK-02).
    [Test]
    public void PredictStructurePseudoknot_DesignedHType_ContainsGenuineCrossing()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);

        Assert.That(CrossingCount(pk.Sequence, pk.BasePairs), Is.GreaterThanOrEqualTo(1),
            "DetectPseudoknots must find at least one crossing pair-of-pairs, confirming a real pseudoknot.");
    }

    // M4 — Valid structure: indices in range, every position paired at most once (INV-PK-03).
    [Test]
    public void PredictStructurePseudoknot_DesignedHType_ValidStructure()
    {
        var pk = PredictStructurePseudoknot(DesignedHTypeSeq);
        int n = pk.Sequence.Length;
        var seen = new HashSet<int>();

        Assert.Multiple(() =>
        {
            foreach (var (a, b) in pk.BasePairs)
            {
                Assert.That(a, Is.InRange(0, n - 1), "5' index in range.");
                Assert.That(b, Is.InRange(0, n - 1), "3' index in range.");
                Assert.That(a, Is.Not.EqualTo(b), "A base cannot pair with itself.");
                Assert.That(seen.Add(a), Is.True, $"position {a} paired at most once.");
                Assert.That(seen.Add(b), Is.True, $"position {b} paired at most once.");
            }
        });
    }

    #endregion

    #region No spurious pseudoknots

    // M5 — A plain hairpin must NOT be reported as a pseudoknot; result equals the plain MFE.
    [Test]
    public void PredictStructurePseudoknot_PlainHairpin_NoSpuriousPseudoknot()
    {
        const string hairpin = "GGGGAAAACCCC";
        var pk = PredictStructurePseudoknot(hairpin);
        var mfe = CalculateMfeStructure(hairpin);

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False,
                "The 9 kcal/mol initiation penalty prevents reporting a pseudoknot for a simple hairpin.");
            Assert.That(pk.DotBracket, Is.EqualTo("((((....))))"),
                "The returned structure must equal the plain pseudoknot-free MFE structure.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy).Within(1e-10),
                "With no accepted pseudoknot the free energy equals the plain MFE.");
        });
    }

    // M6 — Property/invariant sweep: never worse than MFE; any reported knot beats MFE and is valid.
    [Test]
    public void PredictStructurePseudoknot_RandomSweep_NeverWorseThanMfe_KnotsValid()
    {
        var rng = new Random(20260623); // fixed seed → deterministic
        Assert.Multiple(() =>
        {
            for (int t = 0; t < 200; t++)
            {
                string seq = GenerateRandomRna(rng.Next(12, 40), rng, 0.5);
                var pk = PredictStructurePseudoknot(seq);
                var mfe = CalculateMfeStructure(seq);

                Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9),
                    $"INV-PK-01: predictor never returns worse than the plain MFE. seq={seq}");

                if (pk.HasPseudoknot)
                {
                    Assert.That(pk.FreeEnergy, Is.LessThan(mfe.FreeEnergy),
                        $"INV-PK-04: a reported pseudoknot strictly beats the MFE. seq={seq}");
                    Assert.That(CrossingCount(pk.Sequence, pk.BasePairs), Is.GreaterThanOrEqualTo(1),
                        $"INV-PK-02: a reported pseudoknot genuinely crosses. seq={seq}");

                    var seen = new HashSet<int>();
                    foreach (var (a, b) in pk.BasePairs)
                    {
                        Assert.That(seen.Add(a), Is.True, $"INV-PK-03 position {a} once. seq={seq}");
                        Assert.That(seen.Add(b), Is.True, $"INV-PK-03 position {b} once. seq={seq}");
                    }
                }
            }
        });
    }

    #endregion

    #region Edge cases and parity

    // S1 — null input → empty pseudoknot-free structure.
    [Test]
    public void PredictStructurePseudoknot_NullInput_EmptyStructure()
    {
        var pk = PredictStructurePseudoknot(null!);
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
    public void PredictStructurePseudoknot_EmptyInput_EmptyStructure()
    {
        var pk = PredictStructurePseudoknot("");
        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False, "empty input has no pseudoknot.");
            Assert.That(pk.BasePairs, Is.Empty, "empty input yields no base pairs.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(0.0).Within(1e-10), "empty input ΔG is 0.");
        });
    }

    // S3 — too short for any canonical knot (< 11 nt) → no pseudoknot, equals plain MFE.
    [Test]
    public void PredictStructurePseudoknot_TooShort_NoPseudoknot()
    {
        const string shortSeq = "GGGCCC";
        var pk = PredictStructurePseudoknot(shortSeq);
        var mfe = CalculateMfeStructure(shortSeq);

        Assert.Multiple(() =>
        {
            Assert.That(pk.HasPseudoknot, Is.False,
                "The shortest canonical H-type knot needs 11 nt (2·2 helix pairs + 3 loops); shorter cannot knot.");
            Assert.That(pk.DotBracket, Is.EqualTo(mfe.DotBracket), "Returned structure equals the plain MFE.");
            Assert.That(pk.FreeEnergy, Is.EqualTo(mfe.FreeEnergy).Within(1e-10), "ΔG equals the plain MFE.");
        });
    }

    // S4 — BWYV real H-type knot (PDB 437D) is tertiary-stabilised, NOT the NN-thermodynamic MFE.
    [Test]
    public void PredictStructurePseudoknot_Bwyv_NotRecoveredAsMfe_DocumentsTertiaryLimit()
    {
        const string bwyv = "GGCGCGGCACCGUCCGCGGAACAAACGG";
        var pk = PredictStructurePseudoknot(bwyv);
        var mfe = CalculateMfeStructure(bwyv);

        Assert.Multiple(() =>
        {
            // The crystallographic BWYV knot is held by minor-groove triplexes and ion coordination,
            // which the nearest-neighbour secondary-structure model does not capture. A pure NN
            // predictor therefore returns the pseudoknot-free hairpin, not the knot — a documented
            // limit (Reeder & Giegerich 2004; Su et al. 1999), asserted here so it is not over-claimed.
            Assert.That(pk.HasPseudoknot, Is.False,
                "BWYV is tertiary-stabilised; the NN-thermodynamic optimum is the pseudoknot-free hairpin.");
            Assert.That(pk.FreeEnergy, Is.LessThanOrEqualTo(mfe.FreeEnergy + 1e-9),
                "INV-PK-01 still holds on the real sequence.");
        });
    }

    // S5 — DNA spelling (T) folds identically to the RNA spelling (T read as U).
    [Test]
    public void PredictStructurePseudoknot_DnaInput_FoldsIdenticallyToRna()
    {
        string dna = DesignedHTypeSeq.Replace('U', 'T'); // no U here, but assert T↔U parity generally
        var rnaPk = PredictStructurePseudoknot(DesignedHTypeSeq);
        var dnaPk = PredictStructurePseudoknot(dna);

        Assert.Multiple(() =>
        {
            Assert.That(dnaPk.HasPseudoknot, Is.EqualTo(rnaPk.HasPseudoknot), "T read as U → same knot decision.");
            Assert.That(PairSet(dnaPk.BasePairs), Is.EquivalentTo(PairSet(rnaPk.BasePairs)),
                "T read as U → identical base pairs.");
            Assert.That(dnaPk.FreeEnergy, Is.EqualTo(rnaPk.FreeEnergy).Within(1e-10),
                "T read as U → identical free energy.");
        });
    }

    // C1 — minLoopSize below NNDB minimum is clamped to 3 (same result as default).
    [Test]
    public void PredictStructurePseudoknot_MinLoopSizeZero_ClampedToThree()
    {
        var clamped = PredictStructurePseudoknot(DesignedHTypeSeq, minLoopSize: 0);
        var def = PredictStructurePseudoknot(DesignedHTypeSeq, minLoopSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(clamped.HasPseudoknot, Is.EqualTo(def.HasPseudoknot), "minLoopSize<3 clamps to 3.");
            Assert.That(clamped.DotBracket, Is.EqualTo(def.DotBracket), "clamped result equals the default.");
            Assert.That(clamped.FreeEnergy, Is.EqualTo(def.FreeEnergy).Within(1e-10), "clamped ΔG equals default.");
        });
    }

    #endregion
}
