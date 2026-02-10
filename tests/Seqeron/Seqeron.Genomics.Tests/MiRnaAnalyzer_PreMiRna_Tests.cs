using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// MIRNA-PRECURSOR-001: Pre-miRNA Hairpin Detection
/// Canonical test file for FindPreMiRnaHairpins (and indirectly AnalyzeHairpin).
/// Evidence: Bartel (2004), Ambros (2003), Krol (2004), Bartel (2009), miRBase.
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_PreMiRna_Tests
{
    #region Test Data

    // A well-formed hairpin: 22 bp stem (G-C pairs) + 7 nt loop = 51 nt
    // Extended with flanking pairs to reach ≥55 nt.
    // 25 G's + 7 A's loop + 25 C's = 57 nt total, stem=25, loop=7
    private const string ValidHairpin57 =
        "GGGGGGGGGGGGGGGGGGGGGGGGG" +   // 25 nt 5' stem (G)
        "AAAAAAA" +                       // 7 nt loop
        "CCCCCCCCCCCCCCCCCCCCCCCCC";      // 25 nt 3' stem (C pairs with G)

    // 20 G's + 7 A's loop + 20 C's = 47 nt (too short for default min=55)
    private const string ShortHairpin47 =
        "GGGGGGGGGGGGGGGGGGGG" +          // 20 nt
        "AAAAAAA" +                        // 7 nt
        "CCCCCCCCCCCCCCCCCCCC";            // 20 nt

    // Hairpin with G:U wobble pairs in stem
    // G pairs with U (wobble), and U pairs with G (wobble)
    // 20 G's + 7 A's loop + 20 U's = 47 nt; needs extension
    // 25 G's + 7 A's + 25 U's = 57 nt
    private const string WobbleHairpin57 =
        "GGGGGGGGGGGGGGGGGGGGGGGGG" +     // 25 nt 5' stem (G)
        "AAAAAAA" +                        // 7 nt loop
        "UUUUUUUUUUUUUUUUUUUUUUUUU";      // 25 nt 3' stem (U wobble-pairs with G)

    // Stem too short: only 10 bp (< 18 required)
    // 10 G's + 7 A's + 10 C's = 27 nt
    private const string ShortStemHairpin =
        "GGGGGGGGGG" +
        "AAAAAAA" +
        "CCCCCCCCCC";

    // No complementarity: all A's (A does not pair with A)
    private const string NoComplementarity =
        "AAAAAAAAAAAAAAAAAAAAAAAAA" +      // 25 A's
        "AAAAAAAAAAAAAAAAAAAAAAAAA" +      // 25 A's
        "AAAAAAAAAAAAAAAAAAAA";            // 20 A's = 70 total

    // DNA input with T instead of U
    private const string DnaHairpin57 =
        "GGGGGGGGGGGGGGGGGGGGGGGGG" +     // 25 nt (G)
        "TTTTTTT" +                        // 7 nt loop (T instead of A)
        "CCCCCCCCCCCCCCCCCCCCCCCCC";       // 25 nt (C)

    // Large hairpin to test maxHairpinLength filtering
    // 50 G's + 10 A's + 50 C's = 110 nt
    private const string LargeHairpin110 =
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" + // 50
        "AAAAAAAAAA" +                                          // 10
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";   // 50

    // Two hairpins separated by a spacer
    private static readonly string TwoHairpinSequence =
        ValidHairpin57 + new string('A', 20) + ValidHairpin57;

    #endregion

    #region M1: Null Input — Defensive coding

    [Test]
    public void FindPreMiRnaHairpins_NullInput_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins(null!).ToList();

        Assert.That(results, Is.Empty,
            "Null input must return empty enumerable, not throw");
    }

    #endregion

    #region M2: Empty Input — Defensive coding

    [Test]
    public void FindPreMiRnaHairpins_EmptyInput_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins("").ToList();

        Assert.That(results, Is.Empty,
            "Empty string must return empty enumerable");
    }

    #endregion

    #region M3: Short Sequence — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        // 47 nt hairpin is below default minHairpinLength=55
        var results = FindPreMiRnaHairpins(ShortHairpin47, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Sequence shorter than minHairpinLength (55) must return empty — " +
            "pre-miRNAs are ≥55 nt (Bartel 2004, miRBase)");
    }

    #endregion

    #region M4: Valid Hairpin Detection — Bartel (2004), Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_ValidHairpin_DetectsCandidate()
    {
        // 57 nt hairpin: 25 bp G-C stem + 7 nt loop
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1),
            "A 57 nt sequence with 25 bp perfect stem and 7 nt loop must be detected — " +
            "meets all structural criteria (stem ≥ 18 bp, loop 3–25 nt, length ≥ 55)");
    }

    #endregion

    #region M5: Position Correct — Array bounds

    [Test]
    public void FindPreMiRnaHairpins_Position_WithinInputBounds()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must find at least one hairpin");

        foreach (var pre in results)
        {
            Assert.Multiple(() =>
            {
                Assert.That(pre.Start, Is.GreaterThanOrEqualTo(0),
                    "Start must be ≥ 0 (0-based index)");
                Assert.That(pre.End, Is.LessThan(ValidHairpin57.Length),
                    "End must be < input length");
                Assert.That(pre.End, Is.GreaterThanOrEqualTo(pre.Start),
                    "End must be ≥ Start");
            });
        }
    }

    #endregion

    #region M6: Mature From 5' Arm — Bartel (2009)

    [Test]
    public void FindPreMiRnaHairpins_MatureSequence_ExtractedFrom5PrimeArm()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.MatureSequence, Is.Not.Null.And.Not.Empty,
                "Mature sequence must be extracted");
            Assert.That(pre.MatureSequence.Length, Is.LessThanOrEqualTo(22),
                "Mature miRNA is ~22 nt (Bartel 2009)");
            Assert.That(pre.Sequence.StartsWith(pre.MatureSequence, StringComparison.Ordinal),
                "Mature must come from 5' arm (beginning of hairpin sequence)");
        });
    }

    #endregion

    #region M7: Star From 3' Arm — Bartel (2009)

    [Test]
    public void FindPreMiRnaHairpins_StarSequence_ExtractedFrom3PrimeArm()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.StarSequence, Is.Not.Null.And.Not.Empty,
                "Star sequence must be extracted");
            Assert.That(pre.StarSequence.Length, Is.EqualTo(pre.MatureSequence.Length),
                "Star and mature must have equal length (duplex symmetry)");
            Assert.That(pre.Sequence.EndsWith(pre.StarSequence, StringComparison.Ordinal),
                "Star must come from 3' arm (end of hairpin sequence)");
        });
    }

    #endregion

    #region M8: Dot-Bracket Structure — Standard RNA notation

    [Test]
    public void FindPreMiRnaHairpins_Structure_MatchesStemLoopPattern()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];
        string structure = pre.Structure;

        Assert.Multiple(() =>
        {
            Assert.That(structure.Length, Is.EqualTo(pre.Sequence.Length),
                "Structure length must equal sequence length (INV-5)");
            Assert.That(structure, Does.Match(@"^\(+\.+\)+$"),
                "Structure must follow pattern: '('+ for 5' stem, '.'+ for loop, ')'+ for 3' stem");
        });
    }

    #endregion

    #region M9: Dot-Bracket Balanced — Notation definition

    [Test]
    public void FindPreMiRnaHairpins_Structure_BalancedParentheses()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        foreach (var pre in results)
        {
            int openCount = pre.Structure.Count(c => c == '(');
            int closeCount = pre.Structure.Count(c => c == ')');
            int dotCount = pre.Structure.Count(c => c == '.');

            Assert.Multiple(() =>
            {
                Assert.That(openCount, Is.EqualTo(closeCount),
                    "Count of '(' must equal count of ')' — balanced base pairs (INV-7)");
                Assert.That(dotCount, Is.GreaterThanOrEqualTo(3),
                    "Loop must have ≥ 3 unpaired positions (Bartel 2004)");
                Assert.That(openCount + closeCount + dotCount, Is.EqualTo(pre.Structure.Length),
                    "Structure must contain only '(', ')', and '.' characters (INV-6)");
            });
        }
    }

    #endregion

    #region M10: Free Energy Negative — Thermodynamic stability

    [Test]
    public void FindPreMiRnaHairpins_FreeEnergy_IsNegative()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        foreach (var pre in results)
        {
            Assert.That(pre.FreeEnergy, Is.LessThan(0),
                "Free energy must be negative for stable hairpin structures — " +
                "ASSUMPTION: simplified model (-stemLength*1.5 + loopSize*0.5)");
        }
    }

    #endregion

    #region M11: Free Energy Ordering — Turner (2004) principles

    [Test]
    public void FindPreMiRnaHairpins_FreeEnergy_LongerStemMoreNegative()
    {
        // Both hairpins ≥ 55 nt (AnalyzeHairpin hardcodes n < 55 rejection)
        // Hairpin A: 25G + 7A + 25C = 57 nt → stem=25, loop=7
        //   energy = -25*1.5 + 7*0.5 = -34.0
        // Hairpin B: 20G + 15A + 20C = 55 nt → stem=20, loop=15
        //   energy = -20*1.5 + 15*0.5 = -22.5
        string stem20Hairpin55 = new string('G', 20) + new string('A', 15) + new string('C', 20);

        var resultsLong = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();
        var resultsShort = FindPreMiRnaHairpins(stem20Hairpin55, minHairpinLength: 55).ToList();

        Assert.That(resultsLong, Is.Not.Empty, "25 bp stem hairpin (57 nt) must be detected");
        Assert.That(resultsShort, Is.Not.Empty, "20 bp stem hairpin (55 nt) must be detected");

        double energyLong = resultsLong[0].FreeEnergy;
        double energyShort = resultsShort[0].FreeEnergy;

        Assert.That(energyLong, Is.LessThan(energyShort),
            "Longer stem (25 bp) must produce more negative energy than shorter stem (20 bp) — " +
            "reflects greater thermodynamic stability");
    }

    #endregion

    #region M12: Stem Too Short Rejected — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_StemTooShort_ReturnsEmpty()
    {
        // 10 bp stem (< 18 required) + 7 loop = 27 nt
        // Even with lowered minHairpinLength, stem < 18 should reject
        var results = FindPreMiRnaHairpins(ShortStemHairpin, minHairpinLength: 20).ToList();

        Assert.That(results, Is.Empty,
            "Sequence with only 10 bp stem must be rejected — " +
            "pre-miRNA requires ≥ 18 bp stem (Krol 2004)");
    }

    #endregion

    #region M13: Loop Size Boundary Rejected — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_LoopTooLarge_Rejected()
    {
        // Note: Loop-too-small (< 3) is structurally unreachable because
        // maxStem = n/2 - 5 ensures minimum loop ≥ 10 nt for any valid input.
        // Instead test the reachable arm: loop > 25 → rejected.
        // 18G + 30A + 18C = 66 nt. Stem=18 (G-C pairs from ends), loop=30 > 25
        string largeLoop = new string('G', 18) + new string('A', 30) + new string('C', 18);
        var results = FindPreMiRnaHairpins(largeLoop, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Candidate with 30 nt loop must be rejected — " +
            "loop must be ≤ 25 nt (Bartel 2004, implementation threshold)");
    }

    #endregion

    #region M14: DNA T→U Conversion — RNA biology standard

    [Test]
    public void FindPreMiRnaHairpins_DnaInput_ConvertsToRna()
    {
        // DNA input: G's + T's (loop) + C's
        var results = FindPreMiRnaHairpins(DnaHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty,
            "DNA input with T must be handled (T→U conversion)");

        var pre = results[0];
        Assert.That(pre.Sequence, Does.Not.Contain("T"),
            "Output sequence must use U not T (RNA convention) — INV-10");
    }

    #endregion

    #region M15: G:U Wobble Pairs in Stem — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_GU_WobblePairs_AcceptedInStem()
    {
        // 25 G's + 7 A's loop + 25 U's = 57 nt
        // G-U wobble pairs should be accepted by CanPair
        var results = FindPreMiRnaHairpins(WobbleHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1),
            "G-U wobble pairs must be accepted as valid base pairs in stem — " +
            "common in real pre-miRNAs (Krol 2004)");
    }

    #endregion

    #region M16: Sequence Length In Range — Scanning invariant

    [Test]
    public void FindPreMiRnaHairpins_SequenceLength_WithinMinMaxRange()
    {
        const int min = 55;
        const int max = 120;

        var results = FindPreMiRnaHairpins(LargeHairpin110, minHairpinLength: min, maxHairpinLength: max).ToList();

        foreach (var pre in results)
        {
            Assert.That(pre.Sequence.Length, Is.InRange(min, max),
                $"Returned hairpin length {pre.Sequence.Length} must be within [{min}, {max}] — INV-1");
        }
    }

    #endregion

    #region M17: All Invariants Hold — Multiple sources

    [Test]
    public void FindPreMiRnaHairpins_AllInvariants_Verified()
    {
        const int min = 55;
        const int max = 120;
        const int matureLen = 22;

        var results = FindPreMiRnaHairpins(ValidHairpin57, min, max, matureLen).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        foreach (var pre in results)
        {
            Assert.Multiple(() =>
            {
                // INV-1: Length in range
                Assert.That(pre.Sequence.Length, Is.InRange(min, max),
                    "INV-1: Sequence length must be within [minHairpinLength, maxHairpinLength]");

                // INV-2: Start/End in bounds
                Assert.That(pre.Start, Is.GreaterThanOrEqualTo(0),
                    "INV-2: Start ≥ 0");
                Assert.That(pre.End, Is.LessThan(ValidHairpin57.Length),
                    "INV-2: End < input length");

                // INV-3: MatureSequence length
                Assert.That(pre.MatureSequence.Length, Is.GreaterThan(0)
                    .And.LessThanOrEqualTo(matureLen),
                    "INV-3: Mature length in (0, matureLength]");

                // INV-4: Star length matches mature
                Assert.That(pre.StarSequence.Length, Is.EqualTo(pre.MatureSequence.Length),
                    "INV-4: Star length must equal mature length");

                // INV-5: Structure length matches sequence
                Assert.That(pre.Structure.Length, Is.EqualTo(pre.Sequence.Length),
                    "INV-5: Structure length must equal sequence length");

                // INV-6: Only valid characters
                Assert.That(pre.Structure, Does.Match(@"^[(.)]+$"),
                    "INV-6: Structure must contain only '(', ')', '.' characters");

                // INV-7: Balanced parentheses
                int opens = pre.Structure.Count(c => c == '(');
                int closes = pre.Structure.Count(c => c == ')');
                Assert.That(opens, Is.EqualTo(closes),
                    "INV-7: Count of '(' must equal count of ')'");

                // INV-8: Negative free energy
                Assert.That(pre.FreeEnergy, Is.LessThan(0),
                    "INV-8: Free energy must be negative for stable hairpin");

                // INV-10: RNA uppercase
                Assert.That(pre.Sequence, Does.Match(@"^[AUGC]+$"),
                    "INV-10: Sequence must be uppercase RNA (A, U, G, C only)");
            });
        }
    }

    #endregion

    #region S1: Multiple Hairpins in Long Sequence — Biological miRNA clusters

    [Test]
    public void FindPreMiRnaHairpins_LongSequence_FindsMultipleCandidates()
    {
        // Two hairpins separated by 20 A's spacer
        var results = FindPreMiRnaHairpins(TwoHairpinSequence, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(2),
            "Sequence with two distinct hairpin regions should yield ≥ 2 candidates — " +
            "miRNA clusters can contain multiple precursors (Bartel 2004)");
    }

    #endregion

    #region S2: MaxHairpinLength Respected — Parameter contract

    [Test]
    public void FindPreMiRnaHairpins_MaxHairpinLength_EnforcedAsUpperBound()
    {
        // 110 nt hairpin with max=80 should only return sub-windows ≤ 80
        var results = FindPreMiRnaHairpins(LargeHairpin110, minHairpinLength: 55, maxHairpinLength: 80).ToList();

        foreach (var pre in results)
        {
            Assert.That(pre.Sequence.Length, Is.LessThanOrEqualTo(80),
                "No returned hairpin may exceed maxHairpinLength parameter");
        }
    }

    #endregion

    #region S3: Custom MinHairpinLength — Parameter contract

    [Test]
    public void FindPreMiRnaHairpins_CustomMinHairpinLength_Applied()
    {
        // AnalyzeHairpin hardcodes n < 55 rejection, so minHairpinLength < 55
        // has no practical effect. Test above the 55 threshold instead:
        // ValidHairpin57 is exactly 57 nt.
        // With minHairpinLength=57, windows of 57 nt are scanned → detected
        // With minHairpinLength=58, 57 nt input is too short → empty
        var resultsAt57 = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 57).ToList();
        var resultsAt58 = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 58).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsAt57, Is.Not.Empty,
                "57 nt hairpin must be detected with minHairpinLength=57");
            Assert.That(resultsAt58, Is.Empty,
                "57 nt hairpin must be rejected when minHairpinLength=58 (sequence too short)");
        });
    }

    #endregion

    #region S4: No Complementarity — No stem possible

    [Test]
    public void FindPreMiRnaHairpins_NoComplementarity_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins(NoComplementarity, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "All-A sequence has no self-complementarity → no stem → no hairpin");
    }

    #endregion

    #region Helpers

    /// <summary>
    /// Computes expected free energy using the implementation's simplified model.
    /// Formula: -stemLength * 1.5 + loopSize * 0.5
    /// Source: ASSUMPTION — simplified energy model (not Turner nearest-neighbor).
    /// </summary>
    private static double ComputeExpectedEnergy(int stemLength, int loopSize)
    {
        return -stemLength * 1.5 + loopSize * 0.5;
    }

    #endregion
}
