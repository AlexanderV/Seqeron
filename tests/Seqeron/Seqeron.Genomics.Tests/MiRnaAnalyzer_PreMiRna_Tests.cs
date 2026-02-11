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

    // ── Biologically plausible mixed-base hairpins ─────────────────────────
    // All test hairpins use mixed nucleotide composition (AUGC), matching
    // the diversity found in real pre-miRNAs. Mono-nucleotide sequences
    // (all-G stems, all-C stems) are avoided as they are biologically
    // implausible and can mask implementation bugs.

    // 25 bp Watson-Crick stem + 7 nt loop = 57 nt
    // maxStem = 57/2 - 5 = 23 → effective stem=23, effective loop=11
    private const string ValidHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" +   // 25 nt 5' stem (mixed)
        "GAAAUUU" +                       // 7 nt loop
        "UAGCUAGCUAGCUAGCUAGCUAUGC";     // 25 nt 3' stem (reverse complement)

    // 20 bp Watson-Crick stem + 7 nt loop = 47 nt (below default min=55)
    private const string ShortHairpin47 =
        "GCAUAGCUAGCUAGCUAGCU" +          // 20 nt 5' stem
        "GAAAUUU" +                        // 7 nt loop
        "AGCUAGCUAGCUAGCUAUGC";           // 20 nt 3' stem (reverse complement)

    // 25 bp stem with 4 G:U wobble pairs (at pairing positions 0, 5, 9, 17) + 7 nt loop = 57 nt
    // Wobble pairs are common in real pre-miRNAs — Krol (2004)
    private const string WobbleHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" +     // 25 nt 5' stem (same as ValidHairpin57)
        "GAAAUUU" +                        // 7 nt loop
        "UAGCUAGUUAGCUAGUUAGUUAUGU";      // 25 nt 3' stem (4× C→U = G:U wobble)

    // 15 bp stem + 25 nt loop = 55 nt — tests stem rejection (not n<55 hardcode)
    // AnalyzeHairpin rejects stems < 18 bp; this sequence is ≥55 nt total
    // so the rejection must come from stem length, not the n<55 early-exit.
    private const string ShortStemHairpin55 =
        "GCAUAGCUAGCUAGC" +               // 15 nt 5' stem
        "AUGCAUGCAUGCAUGCAUGCAUGCA" +      // 25 nt loop (mixed, non-self-complementary)
        "GCUAGCUAGCUAUGC";                 // 15 nt 3' stem (reverse complement)

    // No complementarity: alternating purines (A-G never pair with each other)
    private const string NoComplementarity =
        "AGAGAGAGAGAGAGAGAGAGAGAGA" +      // 25 nt
        "GAGAGAGAGAGAGAGAGAGAGAGAG" +      // 25 nt
        "AGAGAGAGAGAGAGAGAGAG";            // 20 nt = 70 total

    // DNA input with T instead of U (same stem structure as ValidHairpin57)
    private const string DnaHairpin57 =
        "GCATAGCTAGCTAGCTAGCTAGCTA" +     // 25 nt 5' stem (DNA: T for U)
        "GAAATTT" +                        // 7 nt loop
        "TAGCTAGCTAGCTAGCTAGCTATGC";       // 25 nt 3' stem

    // Large hairpin: 50 bp mixed stem + 10 nt loop = 110 nt
    private static readonly string LargeHairpin110 = BuildPerfectHairpin(50, 10);

    // Two hairpins separated by non-complementary spacer
    private static readonly string TwoHairpinSequence =
        ValidHairpin57 + "AGAGAGAGAGAGAGAGAGAG" + ValidHairpin57;

    // ── Real miRBase pre-miRNA sequences (known-limitation reference) ──────
    // These sequences are from miRBase v22 and represent real pre-miRNAs.
    // The simplified consecutive-pairing model CANNOT detect them because
    // real pre-miRNAs have internal mismatches and bulges that break
    // consecutive pairing from the ends.

    /// <summary>hsa-mir-21 (MI0000077) — 71 nt, miRBase v22</summary>
    private const string HsaMir21_MI0000077 =
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGA" +
        "CUGUUGAAUCUCAUGGCAACACCAGUCGA" +
        "UGGGCUGUCUGACA";

    /// <summary>hsa-let-7a-1 (MI0000060) — 78 nt, miRBase v22</summary>
    private const string HsaLet7a1_MI0000060 =
        "UGGGAUGAGGUAGUAGGUUGUAUAGUUUU" +
        "AGGGUCACACCCACCACUGGGAGAUAACU" +
        "AUACAAUCUACUGUCUUUCCUA";

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
        // Hairpin A: 57 nt, maxStem=23 → effective stem=23, loop=11
        //   energy = -23*1.5 + 11*0.5 = -29.0
        // Hairpin B: 55 nt, maxStem=22, actual stem=20 (pairing breaks at loop)
        //   energy = -20*1.5 + 15*0.5 = -22.5
        string stem20Hairpin55 =
            "GCAUAGCUAGCUAGCUAGCU" +     // 20 nt 5' stem
            "AUGCAUGCAUGCAUG" +           // 15 nt loop
            "AGCUAGCUAGCUAGCUAUGC";       // 20 nt 3' stem (reverse complement)

        var resultsLong = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();
        var resultsShort = FindPreMiRnaHairpins(stem20Hairpin55, minHairpinLength: 55).ToList();

        Assert.That(resultsLong, Is.Not.Empty, "23 bp effective stem hairpin (57 nt) must be detected");
        Assert.That(resultsShort, Is.Not.Empty, "20 bp stem hairpin (55 nt) must be detected");

        double energyLong = resultsLong[0].FreeEnergy;
        double energyShort = resultsShort[0].FreeEnergy;

        Assert.That(energyLong, Is.LessThan(energyShort),
            "Longer effective stem (23 bp) must produce more negative energy than shorter stem (20 bp) — " +
            "reflects greater thermodynamic stability");
    }

    #endregion

    #region M12: Stem Too Short Rejected — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_StemTooShort_ReturnsEmpty()
    {
        // 15 bp stem + 25 nt loop = 55 nt total (passes n<55 hardcode)
        // Rejection must come from stemLength (15) < 18, not from length filter
        var results = FindPreMiRnaHairpins(ShortStemHairpin55, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Sequence with only 15 bp stem must be rejected — " +
            "pre-miRNA requires ≥ 18 bp stem (Krol 2004). " +
            "Note: sequence is 55 nt, so this tests stem rejection, not the n<55 early-exit");
    }

    #endregion

    #region M13: Loop Size Boundary Rejected — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_LoopTooLarge_Rejected()
    {
        // Note: Loop-too-small (< 3) is structurally unreachable because
        // maxStem = n/2 - 5 ensures minimum loop ≥ 10 nt for any valid input.
        // Instead test the reachable arm: loop > 25 → rejected.
        // 18 bp mixed stem + 30 nt loop = 66 nt. Stem=18, loop=30 > 25
        // Loop uses all-A to prevent self-complementary pairing within the loop
        string largeLoop =
            "GCAUAGCUAGCUAGCUAG" +             // 18 nt 5' stem (mixed)
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + // 30 nt loop (A-A never pairs)
            "CUAGCUAGCUAGCUAUGC";              // 18 nt 3' stem (reverse complement)
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
            "All-purine (AG) sequence has no self-complementarity → no stem → no hairpin");
    }

    #endregion

    #region Known Limitation: Real miRBase pre-miRNAs — consecutive-pairing model

    [Test]
    public void FindPreMiRnaHairpins_RealMiRBase_HsaMir21_NotDetected_KnownLimitation()
    {
        // hsa-mir-21 (MI0000077, miRBase v22) — 71 nt
        // Real pre-miRNA with internal mismatches and bulges.
        // Consecutive-pairing model gets only ~16 pairs from ends (< 18 required).
        // This is a KNOWN LIMITATION of the simplified algorithm.
        var results = FindPreMiRnaHairpins(HsaMir21_MI0000077, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "KNOWN LIMITATION: hsa-mir-21 (MI0000077) is a real pre-miRNA from miRBase " +
            "but is not detected because the simplified consecutive-pairing model " +
            "cannot handle internal mismatches and bulges (Assumption #2)");
    }

    [Test]
    public void FindPreMiRnaHairpins_RealMiRBase_HsaLet7a1_NotDetected_KnownLimitation()
    {
        // hsa-let-7a-1 (MI0000060, miRBase v22) — 78 nt
        // Real pre-miRNA; consecutive pairing from ends yields only ~5 pairs.
        var results = FindPreMiRnaHairpins(HsaLet7a1_MI0000060, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "KNOWN LIMITATION: hsa-let-7a-1 (MI0000060) is a real pre-miRNA from miRBase " +
            "but is not detected because the simplified consecutive-pairing model " +
            "cannot handle internal mismatches and bulges (Assumption #2)");
    }

    #endregion

    #region Helpers

    /// <summary>
    /// Builds a perfect hairpin with mixed-base Watson-Crick stem.
    /// Pattern: GCAU cycling for 5' stem, all-A loop, reverse complement for 3' stem.
    /// </summary>
    private static string BuildPerfectHairpin(int stemLength, int loopSize)
    {
        var bases = new[] { 'G', 'C', 'A', 'U' };
        var stem5 = new char[stemLength];
        for (int i = 0; i < stemLength; i++)
            stem5[i] = bases[i % 4];

        string loop = new string('A', loopSize);

        var stem3 = new char[stemLength];
        for (int i = 0; i < stemLength; i++)
        {
            stem3[i] = stem5[stemLength - 1 - i] switch
            {
                'A' => 'U', 'U' => 'A',
                'G' => 'C', 'C' => 'G',
                _ => stem5[stemLength - 1 - i]
            };
        }

        return new string(stem5) + loop + new string(stem3);
    }

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
