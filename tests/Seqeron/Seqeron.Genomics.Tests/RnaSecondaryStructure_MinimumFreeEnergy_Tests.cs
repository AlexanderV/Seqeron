// RNA-MFE-001 — Minimum Free Energy (Zuker–Stiegler DP, Turner 2004 parameters)
// Evidence: docs/Evidence/RNA-MFE-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-MFE-001.md
// Source: Zuker M, Stiegler P (1981). Nucleic Acids Res. 9(1):133-148.
//         Mathews DH et al. (2004). PNAS 101:7287-7292 (Turner 2004 / NNDB).

using NUnit.Framework;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class RnaSecondaryStructure_MinimumFreeEnergy_Tests
{
    #region CalculateMinimumFreeEnergy

    // M1 — NNDB Turner 2004 Hairpin Example 1 reconstructed as a full unimolecular hairpin.
    // Sequence CACAAAAAAAUGUG: stem pairs C-G,A-U,C-G,A-U (5'CACA / 3'UGUG), 6-nt loop AAAAAA.
    // ΔG° = 3 stacks (CG/AU -2.11, AU/CG -2.24, CG/AU -2.11) + AU end +0.45
    //       + terminal mismatch AU·AA -0.8 + hairpin initiation(6) +5.4 = -1.41 (NNDB rounds to -1.4).
    [Test]
    public void CalculateMinimumFreeEnergy_NndbHairpinExample1_ReturnsMinus1_41()
    {
        double expectedExact = -2.11 - 2.24 - 2.11 + 0.45 - 0.8 + 5.4; // = -1.41

        double mfe = CalculateMinimumFreeEnergy("CACAAAAAAAUGUG");

        Assert.Multiple(() =>
        {
            Assert.That(mfe, Is.EqualTo(expectedExact).Within(1e-9),
                "MFE must equal the exact sum of NNDB Example 1 per-term parameters (-1.41 kcal/mol).");
            // NNDB tabulates the total to one decimal place; -1.41 rounds to the published -1.4.
            Assert.That(System.Math.Round(mfe, 1), Is.EqualTo(-1.4).Within(1e-9),
                "MFE rounded to NNDB's one-decimal convention must be the tabulated -1.4 kcal/mol.");
        });
    }

    // M2 — NNDB Turner 2004 Hairpin Example 2 (5-nt loop with GG first mismatch).
    // Sequence CACAGAAAGUGUG: same stem, loop GAAAG (first=last=G).
    // ΔG° = 3 stacks + AU end +0.45 + terminal mismatch AU·GG -0.8 + GG first-mismatch -0.8
    //       + hairpin initiation(5) +5.7 = -1.91 (NNDB rounds to -1.9).
    [Test]
    public void CalculateMinimumFreeEnergy_NndbHairpinExample2_ReturnsMinus1_91()
    {
        double expectedExact = -2.11 - 2.24 - 2.11 + 0.45 - 0.8 - 0.8 + 5.7; // = -1.91

        double mfe = CalculateMinimumFreeEnergy("CACAGAAAGUGUG");

        Assert.Multiple(() =>
        {
            Assert.That(mfe, Is.EqualTo(expectedExact).Within(1e-9),
                "MFE must equal the exact sum of NNDB Example 2 per-term parameters (-1.91 kcal/mol).");
            // NNDB tabulates the total to one decimal place; -1.91 rounds to the published -1.9.
            Assert.That(System.Math.Round(mfe, 1), Is.EqualTo(-1.9).Within(1e-9),
                "MFE rounded to NNDB's one-decimal convention must be the tabulated -1.9 kcal/mol.");
        });
    }

    // M3 — A homopolymer has no complementary bases, so no base pair can form; the only
    // available structure is the open chain (ΔG = 0). Per Zuker–Stiegler the optimum
    // over a set that always contains the 0-energy open chain is 0 here.
    [Test]
    public void CalculateMinimumFreeEnergy_HomopolymerNoPairs_ReturnsZero()
    {
        double mfe = CalculateMinimumFreeEnergy("AAAAAAAA");

        Assert.That(mfe, Is.EqualTo(0.0).Within(1e-9),
            "A homopolymer forms no pairs; MFE is the 0-energy open chain.");
    }

    // M4 — Null and empty input return 0 (no sequence to fold).
    [Test]
    public void CalculateMinimumFreeEnergy_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateMinimumFreeEnergy(""), Is.EqualTo(0.0).Within(1e-9),
                "Empty sequence has no structure; MFE = 0.");
            Assert.That(CalculateMinimumFreeEnergy(null!), Is.EqualTo(0.0).Within(1e-9),
                "Null sequence has no structure; MFE = 0.");
        });
    }

    // M5 — A sequence shorter than minLoopSize+2 cannot enclose a 3-nt hairpin loop,
    // so no pair is possible and MFE = 0 (nearest-neighbor minimum-loop rule).
    [Test]
    public void CalculateMinimumFreeEnergy_ShorterThanMinLoop_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateMinimumFreeEnergy("GCGC"), Is.EqualTo(0.0).Within(1e-9),
                "Length 4 < minLoopSize+2 (=5): cannot form a hairpin; MFE = 0.");
            Assert.That(CalculateMinimumFreeEnergy("GC"), Is.EqualTo(0.0).Within(1e-9),
                "Length 2 is too short to enclose a 3-nt loop; MFE = 0.");
        });
    }

    // M6 — INV-01: MFE is never positive (the 0-energy open chain is always available).
    [Test]
    public void CalculateMinimumFreeEnergy_Invariant_NeverPositive()
    {
        string[] sequences =
        {
            "CACAAAAAAAUGUG", "GGGGGAAAUCCCCC", "GCGCGCAAAGCGCGC",
            "AUAUAUAAAUAUAU", "GCAUGCAUAGCAUGCAU", "AAAAAAAAAA",
        };

        Assert.Multiple(() =>
        {
            foreach (var s in sequences)
                Assert.That(CalculateMinimumFreeEnergy(s), Is.LessThanOrEqualTo(0.0),
                    $"INV-01: MFE must be <= 0 for '{s}' (open chain always available).");
        });
    }

    // M7 — INV-02: MFE is non-increasing under suffix extension. Extending a sequence
    // only adds folding options, so the optimum cannot become larger (less negative).
    [Test]
    public void CalculateMinimumFreeEnergy_Invariant_MonotonicUnderExtension()
    {
        string prefix = "CACAAAAAAAUGUG";          // a stable hairpin, MFE < 0
        string extended = prefix + "GCGCAAAGCGC";   // add another foldable region

        double mfePrefix = CalculateMinimumFreeEnergy(prefix);
        double mfeExtended = CalculateMinimumFreeEnergy(extended);

        Assert.That(mfeExtended, Is.LessThanOrEqualTo(mfePrefix + 1e-9),
            "INV-02: extending the sequence must not raise the MFE (more options can only lower it).");
    }

    #endregion

    #region PredictStructure

    // M8 — PredictStructure on the NNDB Example 1 hairpin must recover the 4 base pairs
    // C-G/A-U/C-G/A-U and the dot-bracket ((((......)))).
    [Test]
    public void PredictStructure_NndbHairpinExample1_PairsAndDotBracket()
    {
        var result = PredictStructure("CACAAAAAAAUGUG");

        Assert.Multiple(() =>
        {
            Assert.That(result.DotBracket, Is.EqualTo("((((......))))"),
                "Predicted structure of the Example 1 hairpin is a 4-bp stem over a 6-nt loop.");
            Assert.That(result.BasePairs.Count, Is.EqualTo(4),
                "The hairpin has exactly 4 base pairs.");
            // Verify the exact pairs (outermost first): (0,13) (1,12) (2,11) (3,10).
            Assert.That(
                System.Linq.Enumerable.Select(result.BasePairs, bp => (bp.Position1, bp.Position2)),
                Is.EquivalentTo(new[] { (0, 13), (1, 12), (2, 11), (3, 10) }),
                "Base pairs must be C-G(0,13), A-U(1,12), C-G(2,11), A-U(3,10).");
        });
    }

    // S1 — Empty input yields an empty structure (no pairs, empty dot-bracket).
    [Test]
    public void PredictStructure_EmptySequence_ReturnsEmptyResult()
    {
        var result = PredictStructure("");

        Assert.Multiple(() =>
        {
            Assert.That(result.BasePairs.Count, Is.EqualTo(0), "Empty sequence has no base pairs.");
            Assert.That(result.DotBracket, Is.EqualTo(""), "Empty sequence has an empty dot-bracket.");
            Assert.That(result.MinimumFreeEnergy, Is.EqualTo(0.0).Within(1e-9), "Empty sequence MFE = 0.");
        });
    }

    // S2 — A homopolymer has no structure: all dots, no pairs.
    [Test]
    public void PredictStructure_HomopolymerNoStructure_AllDots()
    {
        var result = PredictStructure("AAAAAAAAAA");

        Assert.Multiple(() =>
        {
            Assert.That(result.BasePairs.Count, Is.EqualTo(0), "Homopolymer forms no base pairs.");
            Assert.That(result.DotBracket, Is.EqualTo(new string('.', 10)),
                "Homopolymer dot-bracket is all unpaired.");
        });
    }

    #endregion

    #region Stability ordering (COULD)

    // C1 — A GC stem is more stable (more negative MFE) than an AU stem of identical
    // geometry, because G-C stacking energies are far stronger than A-U
    // (NNDB GC/CG -3.42 vs AU/UA -1.10 kcal/mol).
    [Test]
    public void CalculateMinimumFreeEnergy_GcStem_MoreStableThanAuStem()
    {
        double gcStem = CalculateMinimumFreeEnergy("GCGCGCAAAGCGCGC");
        double auStem = CalculateMinimumFreeEnergy("AUAUAUAAAUAUAU");

        Assert.That(gcStem, Is.LessThan(auStem),
            "A GC-rich stem must have a lower (more stable) MFE than an AU-rich stem of equal length.");
    }

    #endregion
}
