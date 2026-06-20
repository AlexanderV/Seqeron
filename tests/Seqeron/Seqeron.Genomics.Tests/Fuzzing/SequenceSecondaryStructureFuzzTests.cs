using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area protein secondary-structure unit — the
/// Chou-Fasman windowed mean-propensity profile.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (DivideByZeroException,
/// KeyNotFoundException, NullReferenceException, OverflowException, NaN result, …).
/// Every input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation outcome (here: an empty enumerable). A raw
/// runtime exception, a NaN, or a hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-SECSTRUCT-001 — Secondary Structure Prediction (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 126.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — empty / null sequence, single residue (n = 1),
///          windowSize boundaries (0, 1, n, n+1), and the degenerate all-unknown
///          window where the per-window known-residue count = 0.
///   • MC = Malformed Content — unknown residues (IUPAC ambiguity B/Z/X/J, gaps,
///          stop codons, digits, unicode / control junk) that have no Chou-Fasman
///          table entry and must be silently excluded from the window mean.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Secondary_Structure_Prediction.md, SEQ-SECSTRUCT-001)
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.PredictSecondaryStructure(string, int windowSize = 7)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 840–870), a lazy iterator yielding (double Helix, double Sheet, double Turn).
///
/// Documented behaviour:
///   • §2.2 / INV-02: for a window, each component is the arithmetic mean of the
///     member residues' propensity (Pα helix, Pβ sheet, Pt turn), over KNOWN residues
///     only:  mean_c(w) = Σ_{r∈w, known} P_c(r) / (#known in w).
///   • §2.4 INV-01: a window of ONE residue returns that residue's (Pα, Pβ, Pt) tuple.
///   • §2.4 INV-03: for an all-known length-n sequence, #emitted windows = max(0, n − w + 1).
///   • §2.4 INV-04: case-insensitive (input ToUpperInvariant'd before lookup).
///   • §2.4 INV-05 / §6.1: UNKNOWN residues (anything outside the 20 standard AAs) are
///     EXCLUDED from a window's count and mean — no KeyNotFound (source uses TryGetValue).
///     A window of ONLY unknown residues emits NOTHING for that position (count = 0,
///     guarded by `count > 0`, so NO DivideByZero).
///   • §2.4 INV-06 / §6.1: null/empty input, windowSize > n, or windowSize < 1 →
///     EMPTY enumerable, no exception.
///
/// Chou-Fasman propensity table pinned in this fixture (Secondary_Structure_Prediction.md
/// §4.2 — the materially behaviour-defining lookup table; propensity = published int ÷ 100):
///   A 1.42/0.83/0.66  R 0.98/0.93/0.95  N 0.67/0.89/1.56  D 1.01/0.54/1.46
///   C 0.70/1.19/1.19  E 1.51/0.37/0.74  Q 1.11/1.10/0.98  G 0.57/0.75/1.56
///   H 1.00/0.87/0.95  I 1.08/1.60/0.47  L 1.21/1.30/0.59  K 1.14/0.74/1.01
///   M 1.45/1.05/0.60  F 1.13/1.38/0.60  P 0.57/0.55/1.52  S 0.77/0.75/1.43
///   T 0.83/1.19/0.96  W 1.08/1.37/0.96  Y 0.69/1.47/1.14  V 1.06/1.70/0.50
/// (each entry is Pα / Pβ / Pt).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceSecondaryStructureFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The 20 standard one-letter residues — the ONLY characters with a Chou-Fasman entry.</summary>
    private const string StandardResidues = "ARNDCEQGHILKMFPSTWYV";

    /// <summary>Independent oracle copy of the documented Chou-Fasman table
    /// (Secondary_Structure_Prediction.md §4.2). Used to cross-check the profile.</summary>
    private static readonly IReadOnlyDictionary<char, (double Helix, double Sheet, double Turn)> Cf =
        new Dictionary<char, (double, double, double)>
        {
            { 'A', (1.42, 0.83, 0.66) }, { 'R', (0.98, 0.93, 0.95) },
            { 'N', (0.67, 0.89, 1.56) }, { 'D', (1.01, 0.54, 1.46) },
            { 'C', (0.70, 1.19, 1.19) }, { 'E', (1.51, 0.37, 0.74) },
            { 'Q', (1.11, 1.10, 0.98) }, { 'G', (0.57, 0.75, 1.56) },
            { 'H', (1.00, 0.87, 0.95) }, { 'I', (1.08, 1.60, 0.47) },
            { 'L', (1.21, 1.30, 0.59) }, { 'K', (1.14, 0.74, 1.01) },
            { 'M', (1.45, 1.05, 0.60) }, { 'F', (1.13, 1.38, 0.60) },
            { 'P', (0.57, 0.55, 1.52) }, { 'S', (0.77, 0.75, 1.43) },
            { 'T', (0.83, 1.19, 0.96) }, { 'W', (1.08, 1.37, 0.96) },
            { 'Y', (0.69, 1.47, 1.14) }, { 'V', (1.06, 1.70, 0.50) }
        };

    /// <summary>The minimum and maximum value present in ANY component of the table — the
    /// envelope that any window mean (an average of in-table values) must stay within.</summary>
    private const double TableMin = 0.37; // E sheet
    private const double TableMax = 1.70; // V sheet

    /// <summary>Independent oracle: the documented windowed mean-propensity profile.
    /// Mirrors the contract exactly — empty when null/empty, w &lt; 1, or w &gt; n; per-window
    /// mean over KNOWN residues; window with no known residue is skipped (emits nothing).</summary>
    private static List<(double Helix, double Sheet, double Turn)> ProfileOracle(string seq, int windowSize)
    {
        var result = new List<(double, double, double)>();
        if (string.IsNullOrEmpty(seq) || windowSize < 1 || windowSize > seq.Length)
            return result;

        string upper = seq.ToUpperInvariant();
        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double h = 0, s = 0, t = 0;
            int count = 0;
            for (int j = 0; j < windowSize; j++)
            {
                if (Cf.TryGetValue(upper[i + j], out var p))
                {
                    h += p.Helix; s += p.Sheet; t += p.Turn; count++;
                }
            }
            if (count > 0)
                result.Add((h / count, s / count, t / count));
        }
        return result;
    }

    /// <summary>The universal well-formedness contract for ANY single emitted tuple: every
    /// component is a finite number (never NaN / ±Infinity — no DivideByZero) and lies within
    /// the table envelope [TableMin, TableMax], because a mean of in-table values stays inside it.</summary>
    private static void AssertWellFormed((double Helix, double Sheet, double Turn) tup)
    {
        foreach (double v in new[] { tup.Helix, tup.Sheet, tup.Turn })
        {
            double.IsFinite(v).Should().BeTrue(
                "every propensity component must be finite — no NaN, no ±Infinity, no DivideByZero (INV-05/INV-06)");
            v.Should().BeInRange(TableMin - Tolerance, TableMax + Tolerance,
                "a window mean of Chou-Fasman propensities stays within the table envelope [0.37, 1.70]");
        }
    }

    private static void AssertWellFormed(IEnumerable<(double Helix, double Sheet, double Turn)> profile)
    {
        foreach (var tup in profile) AssertWellFormed(tup);
    }

    private static void AssertMatchesOracle(IReadOnlyList<(double Helix, double Sheet, double Turn)> actual, string seq, int w)
    {
        var oracle = ProfileOracle(seq, w);
        actual.Count.Should().Be(oracle.Count, "the number of emitted windows must match the documented profile length");
        for (int i = 0; i < oracle.Count; i++)
        {
            actual[i].Helix.Should().BeApproximately(oracle[i].Helix, Tolerance, "window {0} helix mean", i);
            actual[i].Sheet.Should().BeApproximately(oracle[i].Sheet, Tolerance, "window {0} sheet mean", i);
            actual[i].Turn.Should().BeApproximately(oracle[i].Turn, Tolerance, "window {0} turn mean", i);
        }
    }

    /// <summary>Random string of arbitrary BMP code points (0x0000–0xFFFF): control chars,
    /// the null byte, lone surrogate halves, unicode letters/digits — none of which (except
    /// the 20 standard residues) has a Chou-Fasman entry. Fuzz fodder.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-SECSTRUCT-001 — Chou-Fasman propensity profile : fuzz targets (BE, MC)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact propensity / dominant conformation

    /// <summary>
    /// Positive baseline (not a boundary): the documented worked examples must reproduce
    /// EXACTLY. "AE" with window 2 → ((1.42+1.51)/2, (0.83+0.37)/2, (0.66+0.74)/2) =
    /// (1.465, 0.60, 0.70). "AEV" with window 2 emits two windows: [A,E] = (1.465,0.60,0.70)
    /// and [E,V] = (1.285, 1.035, 0.62). Confirms the suite asserts the real BUSINESS contract
    /// (per-window mean of the Chou-Fasman table), not just non-throwing.
    /// — Secondary_Structure_Prediction.md §7.1 worked example + numerical walk-through.
    /// </summary>
    [Test]
    public void Profile_DocumentedWorkedExamples_MatchHandComputedExactly()
    {
        var ae = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 2).Single();
        ae.Helix.Should().BeApproximately(1.465, Tolerance, "(1.42+1.51)/2 = 1.465 (§7.1)");
        ae.Sheet.Should().BeApproximately(0.60, Tolerance, "(0.83+0.37)/2 = 0.60 (§7.1)");
        ae.Turn.Should().BeApproximately(0.70, Tolerance, "(0.66+0.74)/2 = 0.70 (§7.1)");

        var aev = SequenceStatistics.PredictSecondaryStructure("AEV", windowSize: 2).ToList();
        aev.Should().HaveCount(2, "n−w+1 = 3−2+1 = 2 windows (INV-03)");
        aev[0].Helix.Should().BeApproximately(1.465, Tolerance);
        aev[0].Sheet.Should().BeApproximately(0.60, Tolerance);
        aev[0].Turn.Should().BeApproximately(0.70, Tolerance);
        aev[1].Helix.Should().BeApproximately(1.285, Tolerance, "(1.51+1.06)/2 = 1.285");
        aev[1].Sheet.Should().BeApproximately(1.035, Tolerance, "(0.37+1.70)/2 = 1.035");
        aev[1].Turn.Should().BeApproximately(0.62, Tolerance, "(0.74+0.50)/2 = 0.62");
    }

    /// <summary>
    /// Positive baseline: a strongly HELIX-forming homopeptide must be predicted helix-dominant.
    /// poly-Glu (Pα 1.51, the strongest helix former) and poly-Ala (Pα 1.42) yield window
    /// means equal to their own propensities, with Helix &gt; Sheet and Helix &gt; Turn, and a
    /// helix mean &gt; 1.0 ("former"). Conversely a strongly SHEET-forming homopeptide — poly-Val
    /// (Pβ 1.70, the strongest sheet former) and poly-Ile (Pβ 1.60) — must be sheet-dominant
    /// (Sheet &gt; Helix and Sheet &gt; Turn, sheet mean &gt; 1.0). Pins the documented propensity
    /// convention (mean &gt; 1.0 ⇒ favours that conformation).
    /// — Secondary_Structure_Prediction.md §2.1 / §4.2 (E,A helix formers; V,I sheet formers).
    /// </summary>
    [TestCase("EEEEEEE", 1.51, true)]   // poly-Glu, helix-dominant
    [TestCase("AAAAAAA", 1.42, true)]   // poly-Ala, helix-dominant
    [TestCase("VVVVVVV", 1.70, false)]  // poly-Val, sheet-dominant
    [TestCase("IIIIIII", 1.60, false)]  // poly-Ile, sheet-dominant
    public void Profile_StrongFormerHomopeptide_IsDominantInExpectedConformation(
        string seq, double dominantValue, bool helixDominant)
    {
        var profile = SequenceStatistics.PredictSecondaryStructure(seq, windowSize: 5).ToList();
        profile.Should().NotBeEmpty("a length-7 sequence with window 5 emits 3 windows");

        foreach (var tup in profile)
        {
            if (helixDominant)
            {
                tup.Helix.Should().BeApproximately(dominantValue, Tolerance,
                    "a homopeptide window mean equals the residue's Pα");
                tup.Helix.Should().BeGreaterThan(tup.Sheet, "helix former: Pα > Pβ");
                tup.Helix.Should().BeGreaterThan(tup.Turn, "helix former: Pα > Pt");
                tup.Helix.Should().BeGreaterThan(1.0, "a helix former has Pα > 1.0");
            }
            else
            {
                tup.Sheet.Should().BeApproximately(dominantValue, Tolerance,
                    "a homopeptide window mean equals the residue's Pβ");
                tup.Sheet.Should().BeGreaterThan(tup.Helix, "sheet former: Pβ > Pα");
                tup.Sheet.Should().BeGreaterThan(tup.Turn, "sheet former: Pβ > Pt");
                tup.Sheet.Should().BeGreaterThan(1.0, "a sheet former has Pβ > 1.0");
            }
            AssertWellFormed(tup);
        }
    }

    /// <summary>
    /// Positive baseline: a strongly TURN-forming homopeptide is turn-dominant. poly-Gly
    /// (Pt 1.56) and poly-Pro (Pt 1.52) have the highest turn propensities; each window
    /// mean equals its Pt with Turn &gt; Helix and Turn &gt; Sheet and Turn &gt; 1.0.
    /// — Secondary_Structure_Prediction.md §4.2 (G,N,P turn formers).
    /// </summary>
    [TestCase("GGGGGG", 1.56)]
    [TestCase("PPPPPP", 1.52)]
    public void Profile_TurnFormerHomopeptide_IsTurnDominant(string seq, double pt)
    {
        foreach (var tup in SequenceStatistics.PredictSecondaryStructure(seq, windowSize: 4))
        {
            tup.Turn.Should().BeApproximately(pt, Tolerance, "homopeptide turn mean equals Pt");
            tup.Turn.Should().BeGreaterThan(tup.Helix, "turn former: Pt > Pα");
            tup.Turn.Should().BeGreaterThan(tup.Sheet, "turn former: Pt > Pβ");
            tup.Turn.Should().BeGreaterThan(1.0, "a turn former has Pt > 1.0");
            AssertWellFormed(tup);
        }
    }

    /// <summary>
    /// Positive baseline: INV-04 case-insensitivity — the lowercase form of a mixed peptide
    /// yields the exact same profile as the uppercase form. Guards a missing ToUpperInvariant
    /// that would treat lowercase residues as unknown and corrupt the window means.
    /// — Secondary_Structure_Prediction.md §2.4 INV-04 / §3.3.
    /// </summary>
    [Test]
    public void Profile_LowercaseInput_EqualsUppercaseProfile()
    {
        var upper = SequenceStatistics.PredictSecondaryStructure("ACDEFGHIKLMNPQRSTVWY", 7).ToList();
        var lower = SequenceStatistics.PredictSecondaryStructure("acdefghiklmnpqrstvwy", 7).ToList();

        lower.Should().HaveCount(upper.Count, "INV-04: same window count regardless of case");
        for (int i = 0; i < upper.Count; i++)
        {
            lower[i].Helix.Should().BeApproximately(upper[i].Helix, Tolerance);
            lower[i].Sheet.Should().BeApproximately(upper[i].Sheet, Tolerance);
            lower[i].Turn.Should().BeApproximately(upper[i].Turn, Tolerance);
        }
        AssertWellFormed(upper);
    }

    #endregion

    #region BE — Boundary: empty / null (empty enumerable, no exception)

    /// <summary>
    /// BE: the empty string is the lower size boundary — zero residues, no window fits.
    /// The documented result is an EMPTY enumerable (INV-06), with no DivideByZero and no
    /// NullReference. Materialized to force the lazy iterator to run.
    /// — Secondary_Structure_Prediction.md §6.1 / INV-06 (empty → empty result).
    /// </summary>
    [Test]
    public void Profile_EmptyString_IsEmpty_NoThrow()
    {
        var act = () => SequenceStatistics.PredictSecondaryStructure(string.Empty).ToList();

        act.Should().NotThrow("the empty string is a defined boundary, not an error (INV-06)");
        act().Should().BeEmpty("no residues ⇒ no window fits ⇒ empty result");
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 844) — empty enumerable, no NullReferenceException.
    /// — Secondary_Structure_Prediction.md §3.3 (null/empty → empty result).
    /// </summary>
    [Test]
    public void Profile_Null_IsEmpty_NoThrow()
    {
        var act = () => SequenceStatistics.PredictSecondaryStructure(null!).ToList();

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().BeEmpty();
    }

    #endregion

    #region BE — Boundary: windowSize edges (0, 1, n, n+1)

    /// <summary>
    /// BE: windowSize boundaries on a fixed sequence "AEVL" (n=4). w &lt; 1 (0 and the
    /// negative −1, INV-06) → empty; w &gt; n (5) → empty (no window fits); w = n (4) →
    /// exactly one window = the whole-sequence mean; w = 1 → exactly n single-residue
    /// windows (each = that residue's tuple, INV-01). Pins the documented window-count
    /// rule max(0, n−w+1) at every boundary, with NO out-of-range/empty-sequence crash.
    /// — Secondary_Structure_Prediction.md §3.1 / INV-03 / INV-06.
    /// </summary>
    [TestCase(-1, 0)]
    [TestCase(0, 0)]
    [TestCase(1, 4)]
    [TestCase(3, 2)]
    [TestCase(4, 1)]
    [TestCase(5, 0)]
    [TestCase(99, 0)]
    public void Profile_WindowSizeBoundaries_EmitDocumentedWindowCount(int windowSize, int expectedCount)
    {
        const string seq = "AEVL";
        var profile = SequenceStatistics.PredictSecondaryStructure(seq, windowSize).ToList();

        profile.Should().HaveCount(expectedCount,
            "window count is max(0, n−w+1) for n=4, w={0} (INV-03/INV-06)", windowSize);
        AssertMatchesOracle(profile, seq, windowSize);
        AssertWellFormed(profile);
    }

    #endregion

    #region BE — Boundary: single residue (n = 1, INV-01)

    /// <summary>
    /// BE/INV-01: a single recognized residue with window 1 is the n=1 boundary — the lone
    /// window's mean is computed over count = 1, so the emitted tuple must equal that
    /// residue's (Pα, Pβ, Pt) EXACTLY. Verified for ALL 20 standard residues against the
    /// oracle table; guards an off-by-one denominator or a wrong table value, and confirms
    /// no variance/window crash at length 1.
    /// — Secondary_Structure_Prediction.md §2.4 INV-01 (window of one ⇒ its tuple).
    /// </summary>
    [Test]
    public void Profile_SingleStandardResidue_EqualsItsPropensityTuple()
    {
        foreach (char r in StandardResidues)
        {
            var profile = SequenceStatistics.PredictSecondaryStructure(r.ToString(), windowSize: 1).ToList();

            profile.Should().HaveCount(1, $"single residue '{r}', window 1 ⇒ one tuple (INV-01)");
            var (eh, es, et) = Cf[r];
            profile[0].Helix.Should().BeApproximately(eh, Tolerance, $"Pα of '{r}'");
            profile[0].Sheet.Should().BeApproximately(es, Tolerance, $"Pβ of '{r}'");
            profile[0].Turn.Should().BeApproximately(et, Tolerance, $"Pt of '{r}'");
            AssertWellFormed(profile[0]);
        }
    }

    /// <summary>
    /// BE/INV-06: a single residue with the DEFAULT window 7 — the window (7) exceeds the
    /// sequence length (1), so the documented result is empty, NOT a crash. Confirms the
    /// w &gt; n guard at the n=1 boundary.
    /// — Secondary_Structure_Prediction.md §6.1 (windowSize > length → empty).
    /// </summary>
    [Test]
    public void Profile_SingleResidue_DefaultWindow_IsEmpty()
    {
        SequenceStatistics.PredictSecondaryStructure("A").Should().BeEmpty(
            "default window 7 > length 1 ⇒ no window fits ⇒ empty (INV-06)");
    }

    /// <summary>
    /// BE/INV-01/INV-04: a single LOWERCASE recognized residue with window 1 still yields its
    /// tuple (uppercased before lookup) — the n=1 boundary combined with case-insensitivity.
    /// </summary>
    [TestCase('v')]
    [TestCase('e')]
    [TestCase('g')]
    public void Profile_SingleLowercaseResidue_EqualsItsPropensityTuple(char r)
    {
        char up = char.ToUpperInvariant(r);
        var tup = SequenceStatistics.PredictSecondaryStructure(r.ToString(), windowSize: 1).Single();
        var (eh, es, et) = Cf[up];
        tup.Helix.Should().BeApproximately(eh, Tolerance);
        tup.Sheet.Should().BeApproximately(es, Tolerance);
        tup.Turn.Should().BeApproximately(et, Tolerance);
    }

    /// <summary>
    /// BE/MC/INV-05: a single UNKNOWN residue with window 1 — the only residue is not in the
    /// Chou-Fasman table, so the lone window has count 0 and emits NOTHING (empty result),
    /// reached via TryGetValue (no KeyNotFound) and the `count > 0` guard (no DivideByZero).
    /// — Secondary_Structure_Prediction.md §6.1 (window of only unknown residues → no tuple).
    /// </summary>
    [TestCase("X")]
    [TestCase("B")]
    [TestCase("Z")]
    [TestCase("J")]
    [TestCase("-")]
    [TestCase("*")]
    [TestCase("1")]
    [TestCase("?")]
    public void Profile_SingleUnknownResidue_Window1_IsEmpty_NoThrow(string seq)
    {
        var act = () => SequenceStatistics.PredictSecondaryStructure(seq, windowSize: 1).ToList();

        act.Should().NotThrow($"'{seq}' has no Chou-Fasman entry but must not throw (TryGetValue)");
        act().Should().BeEmpty("the only residue is unknown ⇒ count 0 ⇒ no tuple emitted (INV-05)");
    }

    #endregion

    #region MC — Malformed Content: unknown residues excluded from the window mean

    /// <summary>
    /// MC/INV-05: a window's mean is computed over KNOWN residues ONLY. Interleaving unknown
    /// characters among known residues must NOT change the propensity tuple — the unknowns are
    /// excluded from both the sum and the count. With window = full length, "AE" and the
    /// junk-spliced "AX-E*" must produce the SAME single tuple (1.465, 0.60, 0.70), because
    /// only A and E are counted. Guards a length-based denominator that would dilute the mean.
    /// — Secondary_Structure_Prediction.md §2.2 / INV-05 (unknown excluded from count and mean).
    /// </summary>
    [Test]
    public void Profile_UnknownInterleaved_DenominatorExcludesUnknown()
    {
        var clean = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 2).Single();
        var dirty = SequenceStatistics.PredictSecondaryStructure("AX-E*", windowSize: 5).Single();

        clean.Helix.Should().BeApproximately(1.465, Tolerance);
        dirty.Helix.Should().BeApproximately(clean.Helix, Tolerance, "unknowns excluded ⇒ mean over A,E only");
        dirty.Sheet.Should().BeApproximately(clean.Sheet, Tolerance);
        dirty.Turn.Should().BeApproximately(clean.Turn, Tolerance);
        AssertWellFormed(dirty);
    }

    /// <summary>
    /// MC/INV-05: a window of ONLY unknown residues emits nothing, but a profile over a
    /// MIXED sequence still emits a tuple for every window that contains ≥ 1 known residue.
    /// "AXXXXA" with window 6 (= full length) has 2 known residues (both A) ⇒ one tuple equal
    /// to A's own propensities (mean of two identical entries). No DivideByZero, no KeyNotFound.
    /// — Secondary_Structure_Prediction.md §2.2 / §6.1 (skip-and-exclude unknowns).
    /// </summary>
    [Test]
    public void Profile_WindowWithMostlyUnknown_AveragesKnownOnly()
    {
        var tup = SequenceStatistics.PredictSecondaryStructure("AXXXXA", windowSize: 6).Single();
        var (ah, ash, at) = Cf['A'];
        tup.Helix.Should().BeApproximately(ah, Tolerance, "mean of the two known A residues = A's Pα");
        tup.Sheet.Should().BeApproximately(ash, Tolerance);
        tup.Turn.Should().BeApproximately(at, Tolerance);
    }

    /// <summary>
    /// MC: an all-unknown sequence — EVERY window has count 0, so NO tuple is emitted at any
    /// position: the result is empty. Holds at the n=1 boundary, at scale, case-insensitively,
    /// and across the common non-standard codes. No DivideByZero, no KeyNotFound.
    /// — Secondary_Structure_Prediction.md §6.1 (all-unknown window emits nothing).
    /// </summary>
    [TestCase("X", 1)]
    [TestCase("XX", 1)]
    [TestCase("XXXXXXXX", 3)]
    [TestCase("xxxxx", 2)]
    [TestCase("BZJ-*0123", 4)]
    public void Profile_AllUnknown_IsEmpty_NoThrow(string seq, int windowSize)
    {
        var act = () => SequenceStatistics.PredictSecondaryStructure(seq, windowSize).ToList();

        act.Should().NotThrow($"all-unknown '{seq}' must never throw (TryGetValue + count>0 guard)");
        act().Should().BeEmpty("no window has a known residue ⇒ no tuple emitted");
    }

    /// <summary>
    /// MC: a long all-X homopolymer — the all-unknown boundary at scale. Still empty, still
    /// no DivideByZero / KeyNotFound, and bounded in time (hang guard).
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Profile_LongAllX_IsEmpty_NoCrash()
    {
        SequenceStatistics.PredictSecondaryStructure(new string('X', 5000), windowSize: 7)
            .Should().BeEmpty("no X is recognized ⇒ empty profile at scale");
    }

    #endregion

    #region MC / BE — Random garbage: never throws, always matches the oracle

    /// <summary>
    /// MC/BE: a large batch of arbitrary BMP strings (control chars, the null byte, lone
    /// surrogate halves, unicode letters/digits, occasionally seeded with real residues),
    /// fuzzed over random window sizes (incl. 0, negatives, and oversize), must NEVER throw
    /// and must ALWAYS produce a profile that matches the independent oracle exactly and is
    /// well-formed (finite, within the table envelope). Core fuzz guarantee: no DivideByZero,
    /// no KeyNotFound, no NaN, no overflow, no hang on garbage of any shape/length (incl. 0).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Profile_RandomGarbageStrings_NeverThrow_MatchOracle()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);
            int windowSize = rng.Next(-2, 210); // include w<1, w in range, w>n

            List<(double Helix, double Sheet, double Turn)> profile = null!;
            var act = () => profile = SequenceStatistics.PredictSecondaryStructure(input, windowSize).ToList();

            act.Should().NotThrow($"garbage input (len {len}, w {windowSize}) must never crash");
            AssertMatchesOracle(profile, input, windowSize);
            AssertWellFormed(profile);
        }
    }

    /// <summary>
    /// MC/BE: randomly built sequences over the FULL standard residue alphabet (n ≥ 1) with
    /// random valid window sizes must equal the independent oracle profile exactly over many
    /// shapes — cross-checks the per-window summation, the known-residue denominator, and the
    /// window-count rule, and confirms every emitted tuple is well-formed.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Profile_RandomStandardPeptides_MatchOracleProfile()
    {
        var rng = new Random(13572468);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(1, 300);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = StandardResidues[rng.Next(StandardResidues.Length)];
            string seq = new string(chars);
            int windowSize = rng.Next(1, len + 1); // valid window in [1, n]

            var profile = SequenceStatistics.PredictSecondaryStructure(seq, windowSize).ToList();

            profile.Should().HaveCount(len - windowSize + 1, "INV-03: window count = n−w+1 for all-known input");
            AssertMatchesOracle(profile, seq, windowSize);
            AssertWellFormed(profile);
        }
    }

    /// <summary>
    /// MC/BE: standard residues randomly polluted with unknown junk, over random window sizes
    /// — the full mixed fuzz surface. Must never throw and must always match the oracle (which
    /// excludes unknowns from each window mean and drops all-unknown windows) and stay
    /// well-formed. Exercises the documented skip-and-exclude path at random scale.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Profile_StandardResiduesPollutedWithJunk_MatchOracle()
    {
        var rng = new Random(99887766);
        const string junk = "XBZJ-*0123 .?\t\n";

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(1, 120);
            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < len; i++)
            {
                // ~60% known residue, ~40% junk — guarantees a mix of known/unknown windows.
                if (rng.NextDouble() < 0.6)
                    sb.Append(StandardResidues[rng.Next(StandardResidues.Length)]);
                else
                    sb.Append(junk[rng.Next(junk.Length)]);
            }
            string seq = sb.ToString();
            int windowSize = rng.Next(1, seq.Length + 1);

            List<(double Helix, double Sheet, double Turn)> profile = null!;
            var act = () => profile = SequenceStatistics.PredictSecondaryStructure(seq, windowSize).ToList();

            act.Should().NotThrow($"polluted input (len {seq.Length}, w {windowSize}) must never crash");
            AssertMatchesOracle(profile, seq, windowSize);
            AssertWellFormed(profile);
        }
    }

    #endregion
}
