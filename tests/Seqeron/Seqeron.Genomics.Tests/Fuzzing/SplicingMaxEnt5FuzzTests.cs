using System;
using System.Collections.Generic;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Splicing area — MaxEntScan <c>score5ss</c> maximum-entropy
/// 5' (donor) splice-site model (SPLICE-MAXENT5-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception. Every input must resolve to
/// EITHER a well-defined, theory-correct result, OR a *documented, intentional*
/// validation outcome (here a typed Argument… exception). — docs/ADVANCED_TESTING_CHECKLIST.md
/// §8 "Fuzzing"; docs/checklists/03_FUZZING.md §Description (MC = Malformed Content,
/// BE = Boundary Exploitation).
///
/// For MaxEntScan score5ss the headline hazards a fuzzer must rule out are:
///   • IndexOutOfRange — from a window whose length ≠ 9 reaching the fixed GT-position
///     (3..4) reads or the 7-nt "rest" Span slices. The contract REJECTS any wrong
///     length up front with an ArgumentException (NOT an IndexOutOfRange).
///   • KeyNotFoundException — from a non-A/C/G/T(/U) character producing an out-of-table
///     7-mer rest key, or from the GT-position consensus/background dictionary lookups.
///     The contract REJECTS any non-A/C/G/T(/U) character with an ArgumentException
///     (via NormalizeNucleotide) before any table/dictionary lookup.
///   • NaN / ±Infinity — from log2(0) when the maxent rest probability is zero. The score
///     is log2(P_maxent/P_background); for ANY all-A/C/G/T(/U) 9-mer the embedded rest-7mer
///     probability is strictly positive and the consensus/background GT term is strictly
///     positive, so the product is positive and the log2 is FINITE. A NaN/Inf on a valid
///     9-mer would be a real bug.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The score5ss contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The 5' splice DONOR site marks the START of an intron; its consensus is
/// <c>MAG|GTRAGT</c>, with the almost-invariant intronic <c>GT</c> immediately after the
/// exon/intron junction. MaxEntScan score5ss (Yeo &amp; Burge 2004) models a fixed 9-nt
/// window — 3 exonic + 6 intronic nt, with the conserved GT at 0-based window positions
/// 3–4 — as
///   score = log2( P_maxent(window) / P_background(window) )   (in bits),
/// computed by removing the GT (scored by a consensus/background dinucleotide model) and
/// looking up the maximum-entropy probability of the remaining 7 ("rest") positions
/// (<c>window[0:3] + window[5:9]</c>) DIRECTLY in a single embedded table — unlike score3,
/// score5 is single-matrix (no overlapping sub-windows). score5ss is CALIBRATED so true 5'
/// sites score positive: the canonical documented example
/// score5('cagGTAAGT') == 10.858313 bits ("10.86"); a stronger site 'gagGTAAGT' == 11.078494
/// ("11.08"); a weak non-GT site 'taaATAAGT' == −0.116791 ("−0.12").
/// — docs/algorithms/Splicing/Donor_Site_Detection.md §5.3 (MaxEntScan paragraph), §6.1;
///   tests/TestSpecs/SPLICE-MAXENT5-001.md (validated full-precision oracle);
///   src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs (#region MaxEntScan score5ss).
///
/// Method under test (src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs ~L1347):
///   double ScoreDonorMaxEnt(string window)
///     → the score5ss in bits.
///
/// Documented input handling (XML on ScoreDonorMaxEnt; §5.3; §6.1 "Lowercase input"):
///   • null window               → ArgumentNullException.
///   • window.Length ≠ 9         → ArgumentException (NOT IndexOutOfRange).
///   • non-A/C/G/T(/U) character → ArgumentException (NOT KeyNotFound / IndexOutOfRange),
///                                  via NormalizeNucleotide.
///   • lowercase / mixed case    → folded via ToUpperInvariant(), then scored — NOT rejected.
///   • T and U are equivalent (U is normalised to T in both the rest key and the GT model).
///   • The GT at 3–4 is NOT *required* by the math: a non-GT window is SCORED (the consensus
///     GT term heavily penalises it → score drops far below the canonical GT window), never
///     rejected. The GT is only the biological consensus, not an input precondition — the
///     documented weak example 'taaATAAGT' (dinucleotide "AT", not "GT") is itself scored.
///
/// Theory-correct invariants asserted (independently derived from §5.3 / Yeo &amp; Burge 2004
/// and the maxentpy score5 reference recorded in tests/TestSpecs/SPLICE-MAXENT5-001.md):
///   • FINITE — every all-A/C/G/T(/U) 9-mer yields a finite score (no NaN, no ±Inf): the
///     score is log2 of (strictly-positive rest-7mer probability × strictly-positive GT
///     consensus/background term). Pinned over hundreds of LOCALLY-seeded random 9-mers AND
///     every A/C/G/T homopolymer.
///   • DOCUMENTED ANCHORS — the published maxentpy score5 reference values (the independent
///     oracle, recorded verbatim in the spec — NOT read off the code's own tables):
///       score5('cagGTAAGT') ≈ 10.858313  (canonical, "10.86 bits")
///       score5('gagGTAAGT') ≈ 11.078494  (stronger site, "11.08")
///       score5('taaATAAGT') ≈ −0.116791  (weak non-GT site, "−0.12")
///   • CONSENSUS &gt; RANDOM — a strong consensus donor (CAG|GTAAGT) scores ABOVE a random /
///     GC-rich 9-mer; score5ss is calibrated so true sites score high.
///   • CASE / T==U INVARIANCE — lowercase, uppercase and mixed-case spellings of the same
///     window, and the T-vs-U spelling, all produce the identical score.
///   • GT-NOT-REQUIRED — replacing the 3–4 GT with a non-GT dinucleotide does NOT throw and
///     drops the score far below the canonical GT window (the model scores it, never rejects).
///
/// Discipline: tests encode the DOCUMENTED theory, derived independently; a test that
/// passes against a wrong implementation is invalid. If a source-derived test and the code
/// disagree the CODE is wrong (fixed minimally per the doc). A NaN/Inf score or an
/// IndexOutOfRange/KeyNotFound on a malformed window would be a REAL bug.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SplicingMaxEnt5FuzzTests
{
    #region SPLICE-MAXENT5-001 — MaxEntScan score5ss 5' donor maximum-entropy model

    // The canonical documented MaxEntScan score5 examples and their published full-precision
    // maxentpy reference values (the independent oracle; Donor_Site_Detection.md §5.3 and
    // tests/TestSpecs/SPLICE-MAXENT5-001.md — NOT read off the embedded tables).
    private const string CanonicalWindow = "cagGTAAGT"; // 9 nt, GT at 3-4
    private const double CanonicalScore = 10.858313;    // "10.86 bits"
    private const string StrongWindow = "gagGTAAGT";
    private const double StrongScore = 11.078494;       // "11.08"
    private const string WeakWindow = "taaATAAGT";      // dinucleotide at 3-4 is "AT", NOT "GT"
    private const double WeakScore = -0.116791;         // "−0.12"

    private const int WindowLength = 9;
    private const int GtStart = 3; // 0-based position of the conserved GT dinucleotide

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz windows are reproducible.</summary>
    private static string RandomWindow(int length, int seed, string alphabet = "ACGT")
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    /// <summary>A random 9-mer whose 3–4 dinucleotide is forced to the consensus GT.</summary>
    private static string RandomGtWindow(int seed)
    {
        char[] w = RandomWindow(WindowLength, seed).ToCharArray();
        w[GtStart] = 'G';
        w[GtStart + 1] = 'T';
        return new string(w);
    }

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity — the documented anchors (proves we score a real model)
    // ───────────────────────────────────────────────────────────────────

    #region Anchor — documented maxentpy reference values are reproduced exactly

    /// <summary>
    /// Positive sanity / anchor: the three published MaxEntScan score5 reference windows
    /// reproduce their documented full-precision values exactly (10.858313 / 11.078494 /
    /// −0.116791; Donor_Site_Detection.md §5.3, tests/TestSpecs/SPLICE-MAXENT5-001.md,
    /// Yeo &amp; Burge 2004 / maxentpy). This pins the fuzz harness to a detector that actually
    /// computes the real maximum-entropy score — not a no-op — so the boundary/malformed
    /// assertions below are meaningful. The values are the INDEPENDENT oracle (the maxentpy
    /// reference), never read off the code's tables.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_DocumentedReferenceWindows_ReproduceExactValues()
    {
        ScoreDonorMaxEnt(CanonicalWindow).Should().BeApproximately(CanonicalScore, 1e-6,
            "the canonical documented window scores 10.858313 bits (maxentpy score5 oracle, \"10.86\")");
        ScoreDonorMaxEnt(StrongWindow).Should().BeApproximately(StrongScore, 1e-6,
            "the stronger-site documented window scores 11.078494 bits (\"11.08\")");
        ScoreDonorMaxEnt(WeakWindow).Should().BeApproximately(WeakScore, 1e-6,
            "the weak non-GT documented window scores −0.116791 bits (\"−0.12\", a calibrated near-zero site)");
    }

    /// <summary>
    /// Ordering anchor (Donor_Site_Detection.md §5.3 — score5ss is calibrated so true 5'
    /// sites score positive). The two consensus GT donors score strictly above zero, the
    /// stronger ranks above the canonical, and the weak non-GT site ranks far below both.
    /// This pins the calibration ordering of the documented examples, independent of the
    /// exact bit values.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_DocumentedExamples_RankByCalibration()
    {
        double canonical = ScoreDonorMaxEnt(CanonicalWindow);
        double strong = ScoreDonorMaxEnt(StrongWindow);
        double weak = ScoreDonorMaxEnt(WeakWindow);

        canonical.Should().BeGreaterThan(0.0, "a CAG|GTAAGT consensus donor is a calibrated true 5' site → positive bits");
        strong.Should().BeGreaterThan(canonical, "GAG|GTAAGT is the stronger documented site (11.08 > 10.86)");
        weak.Should().BeLessThan(canonical, "the weak non-GT site ranks far below the canonical GT consensus");
    }

    /// <summary>
    /// CONSENSUS &gt; RANDOM anchor (Donor_Site_Detection.md §5.3 — score5ss favours the donor
    /// consensus). A strong consensus donor (the exonic CAG followed by the intronic GTAAGT)
    /// must score strictly ABOVE a random and a GC-rich 9-mer of the same length. We sweep
    /// LOCALLY-seeded random windows (each forced to carry the GT at 3–4 so only the
    /// background varies) so the ordering is pinned against arbitrary backgrounds, not one
    /// hand-picked decoy.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_ConsensusDonor_RanksAboveRandomAndGcRich()
    {
        double consensusScore = ScoreDonorMaxEnt(CanonicalWindow); // CAG|GTAAGT

        // GC-rich decoy carrying the same terminal GT consensus: a GC-rich context is not the
        // donor consensus, so the maximum-entropy model scores it far below the consensus.
        double gcRich = ScoreDonorMaxEnt("gcgGTcgcg"); // GT at 3-4, GC-rich flanks
        consensusScore.Should().BeGreaterThan(gcRich,
            "a true CAG|GTAAGT donor outscores a GC-rich 9-mer (score5ss favours the donor consensus)");

        // Random A/C/G/T 9-mers (always with the GT at 3–4 so only the background varies):
        // the consensus donor must outrank every random background.
        foreach (int seed in new[] { 7, 41, 137, 911, 2026 })
        {
            string random = RandomGtWindow(seed);
            double randomScore = ScoreDonorMaxEnt(random);
            consensusScore.Should().BeGreaterThan(randomScore,
                $"the CAG|GTAAGT consensus must outscore the random background (seed {seed}, '{random}')");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  FINITE — every valid 9-mer yields a finite score (no NaN / ±Inf)
    // ───────────────────────────────────────────────────────────────────

    #region BE/MC — finite score on any A/C/G/T(/U) 9-mer (no log2(0) → no NaN/Inf)

    /// <summary>
    /// The headline numeric invariant: for ANY all-A/C/G/T(/U) 9-mer the score is FINITE —
    /// never NaN, never ±Infinity. The score is log2(P_maxent/P_background); the rest-7mer
    /// probability is a strictly-positive embedded table entry and the GT consensus/background
    /// term is strictly positive, so the argument of log2 is strictly positive and the result
    /// is finite (Yeo &amp; Burge 2004; SpliceSitePredictor.ScoreDonorMaxEnt). We sweep hundreds
    /// of LOCALLY-seeded random 9-mers across A/C/G/T spellings; a NaN/Inf on any valid window
    /// would be a real bug.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_AnyAcgtWindow_ProducesFiniteScore()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            string window = RandomWindow(WindowLength, seed);
            window.Length.Should().Be(WindowLength, "the fuzz fixture builds an exactly-9-nt window");

            double score = ScoreDonorMaxEnt(window);

            double.IsNaN(score).Should().BeFalse($"score5 of an A/C/G/T 9-mer must not be NaN (seed {seed}, '{window}')");
            double.IsInfinity(score).Should().BeFalse($"score5 of an A/C/G/T 9-mer must not be ±Inf (seed {seed}, '{window}')");
            double.IsFinite(score).Should().BeTrue($"score5 of an A/C/G/T 9-mer is finite (seed {seed}, '{window}')");
        }
    }

    /// <summary>
    /// The four homopolymer corners (AAA…/CCC…/GGG…/TTT…) are the extreme low-information
    /// windows — none carries the GT consensus, all four are valid A/C/G/T(/U) 9-mers, and
    /// each must still produce a FINITE score (no zero-probability rest-7mer anywhere). They
    /// are scored, never rejected — the GT is not an input precondition. The U-spelling
    /// homopolymer must match its T-spelling exactly (T==U).
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_Homopolymers_ScoredAndFinite()
    {
        foreach (char b in new[] { 'A', 'C', 'G', 'T' })
        {
            string window = new string(b, WindowLength);
            double score = ScoreDonorMaxEnt(window);
            double.IsFinite(score).Should().BeTrue($"the {b}-homopolymer 9-mer scores a finite value, not NaN/Inf");
        }

        // T- and U-homopolymers are the same window after U→T normalisation → identical score.
        ScoreDonorMaxEnt(new string('U', WindowLength))
            .Should().Be(ScoreDonorMaxEnt(new string('T', WindowLength)),
                "T==U: the U-homopolymer scores identically to the T-homopolymer");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — window length ≠ 9 → documented ArgumentException, NOT IndexOutOfRange
    // ───────────────────────────────────────────────────────────────────

    #region BE — wrong-length and empty windows throw ArgumentException (not IndexOutOfRange)

    /// <summary>
    /// Target "window ≠ 9 nt" / "empty seq" (KEY boundary). The length is validated up front:
    /// ANY window whose length ≠ 9 — including the empty string, lengths just below/above 9,
    /// and a window that merely PROVISIONALLY looks donor-like — throws an ArgumentException
    /// (a DOCUMENTED validation outcome), NOT an IndexOutOfRangeException off the fixed
    /// GT-position reads (3..4) or the 7-nt rest Span slices. We sweep lengths 0..40 except 9,
    /// building each from a real A/C/G/T alphabet so the rejection is on LENGTH, not alphabet.
    /// — ScoreDonorMaxEnt XML (&lt;exception ArgumentException&gt; "not exactly 9 nt").
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_WindowLengthNot9_ThrowsArgumentException()
    {
        for (int len = 0; len <= 40; len++)
        {
            if (len == WindowLength)
                continue;

            string window = RandomWindow(len, seed: 1000 + len);
            window.Length.Should().Be(len);

            Action act = () => ScoreDonorMaxEnt(window);
            act.Should().Throw<ArgumentException>(
                    $"a {len}-nt window (≠ 9) is rejected up front, never crashing with IndexOutOfRange")
                .Which.Should().NotBeOfType<ArgumentNullException>(
                    "a non-null wrong-length window is an ArgumentException, not a null-argument exception");
        }
    }

    /// <summary>
    /// The empty and whitespace-padded boundaries explicitly: "" (length 0) and a window that
    /// is 9 visible nt but padded with whitespace to a wrong length both throw ArgumentException
    /// on length. An empty window must NOT slip through to read window[3] (IndexOutOfRange);
    /// padding makes a *would-be* valid window the wrong length, which length-validation catches
    /// BEFORE the alphabet check ever sees the whitespace.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_EmptyAndPaddedWindows_ThrowArgumentException()
    {
        foreach (string bad in new[]
                 {
                     "",                 // empty → length 0
                     "   ",              // whitespace only → length 3
                     " cagGTAAGT ",      // 9-nt window padded to length 11
                     "cagGTAAGT ",       // trailing space → length 10
                     "cagGTAAG",         // one short → length 8
                     "cagGTAAGTA",       // one long → length 10
                 })
        {
            Action act = () => ScoreDonorMaxEnt(bad);
            act.Should().Throw<ArgumentException>(
                $"a length-{bad.Length} window (≠ 9) is rejected on length, never IndexOutOfRange ('{bad}')");
        }
    }

    /// <summary>
    /// null throws ArgumentNullException (the documented null contract), distinct from the
    /// wrong-length ArgumentException above — and never a NullReferenceException leaking from
    /// an unchecked .Length / .ToUpperInvariant() access.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_NullWindow_ThrowsArgumentNullException()
    {
        Action act = () => ScoreDonorMaxEnt(null!);
        act.Should().Throw<ArgumentNullException>("a null window is the documented null contract");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — non-A/C/G/T(/U) characters → documented ArgumentException
    // ───────────────────────────────────────────────────────────────────

    #region MC — non-ACGT characters throw ArgumentException (not KeyNotFound/IndexOutOfRange)

    /// <summary>
    /// Target "non-ACGT" (Malformed Content). A length-correct (9-nt) window that contains a
    /// character outside {A,C,G,T,U} — IUPAC ambiguity codes (N, R, Y…), digits, punctuation —
    /// throws an ArgumentException (the DOCUMENTED handling, via NormalizeNucleotide), NOT a
    /// KeyNotFoundException from an out-of-table rest 7-mer key and NOT an IndexOutOfRange. The
    /// invalid character may sit in the GT positions (consensus/background dictionary lookup) OR
    /// anywhere in the 7-nt rest (the embedded table key) — both paths must reject cleanly.
    /// — ScoreDonorMaxEnt XML (&lt;exception ArgumentException&gt; "contains a non-A/C/G/T(/U) character").
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_NonAcgtCharacter_ThrowsArgumentException()
    {
        foreach (string bad in new[]
                 {
                     "cagGTAAGN",   // N in the rest (exon-distal) — rest 7-mer key
                     "Nagctagct",   // leading junk — rest 7-mer key
                     "cagNTAAGT",   // N at GT position 3 (the G) — consensus lookup
                     "cagGNAAGT",   // N at GT position 4 (the T) — consensus lookup
                     "cagGTAAG!",   // punctuation
                     "cagGTAAG5",   // a digit
                     "cagRYAAGT",   // IUPAC ambiguity codes at the GT positions
                     "rrrGTrrrr",   // an ambiguity-code tract in the rest
                 })
        {
            bad.Length.Should().Be(WindowLength, $"the malformed-content fixture stays length 9 ('{bad}')");

            Action act = () => ScoreDonorMaxEnt(bad);
            act.Should().Throw<ArgumentException>(
                $"a non-A/C/G/T(/U) character is rejected, never KeyNotFound/IndexOutOfRange ('{bad}')");
        }
    }

    /// <summary>
    /// An all-ambiguity ("NNN…") 9-nt window — length-correct but with no concrete nucleotide
    /// anywhere — throws ArgumentException, never wedging or leaking a KeyNotFound from the rest
    /// 7-mer lookup. This is the malformed-content corner where every position is junk.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_AllAmbiguityWindow_ThrowsArgumentException()
    {
        string window = new string('N', WindowLength);
        Action act = () => ScoreDonorMaxEnt(window);
        act.Should().Throw<ArgumentException>(
            "an all-N 9-mer has no valid nucleotide and is rejected, not crashed");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — lowercase / mixed case is folded (documented), not rejected
    // ───────────────────────────────────────────────────────────────────

    #region MC — lowercase/mixed-case folding and T==U equivalence

    /// <summary>
    /// Target "lowercase" (documented: case-fold to ACGT, not reject — ScoreDonorMaxEnt
    /// uppercases via ToUpperInvariant; Donor_Site_Detection.md §6.1 "Lowercase input").
    /// The fully-lowercase, fully-uppercase and a deliberately mixed-case spelling of the SAME
    /// 9-mer all produce the IDENTICAL score — case carries no information. We pin this on the
    /// three documented reference windows AND on LOCALLY-seeded random GT windows so
    /// case-folding holds for arbitrary content, not just the anchors.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_CaseFolding_AllSpellingsScoreIdentically()
    {
        var windows = new List<string> { CanonicalWindow, StrongWindow, WeakWindow };
        for (int seed = 50; seed < 60; seed++)
            windows.Add(RandomGtWindow(seed));

        foreach (string window in windows)
        {
            string lower = window.ToLowerInvariant();
            string upper = window.ToUpperInvariant();
            // Alternating-case spelling: same letters, scrambled case.
            char[] mixedChars = lower.ToCharArray();
            for (int i = 0; i < mixedChars.Length; i += 2)
                mixedChars[i] = char.ToUpperInvariant(mixedChars[i]);
            string mixed = new string(mixedChars);

            double lowerScore = ScoreDonorMaxEnt(lower);
            double upperScore = ScoreDonorMaxEnt(upper);
            double mixedScore = ScoreDonorMaxEnt(mixed);

            upperScore.Should().Be(lowerScore,
                $"case is folded: '{upper}' must score identically to '{lower}'");
            mixedScore.Should().Be(lowerScore,
                $"mixed case is folded: '{mixed}' must score identically to '{lower}'");
            double.IsFinite(lowerScore).Should().BeTrue("a folded valid window still scores finite");
        }
    }

    /// <summary>
    /// T==U equivalence (the model normalises U→T in both the rest 7-mer key and the GT
    /// consensus model, so DNA T and RNA U are identical). The T-spelling and the U-spelling
    /// of the same windows — in any case — score identically, across the documented anchors and
    /// random windows. This is the alphabet-equivalence boundary: a fuzzer must not see U
    /// treated as out-of-alphabet junk.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_TandUSpellings_ScoreIdentically()
    {
        var windows = new List<string> { CanonicalWindow, StrongWindow, WeakWindow };
        for (int seed = 70; seed < 80; seed++)
            windows.Add(RandomGtWindow(seed));

        foreach (string window in windows)
        {
            string tForm = window.ToUpperInvariant();                   // DNA spelling
            string uForm = tForm.Replace('T', 'U');                     // RNA spelling
            string uLower = tForm.ToLowerInvariant().Replace('t', 'u'); // lowercase RNA

            double tScore = ScoreDonorMaxEnt(tForm);
            ScoreDonorMaxEnt(uForm).Should().Be(tScore,
                $"T==U: the U-form '{uForm}' scores identically to the T-form '{tForm}'");
            ScoreDonorMaxEnt(uLower).Should().Be(tScore,
                $"T==U and case-folding compose: lowercase RNA '{uLower}' scores identically");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — the GT at 3-4 is NOT required by the math (scored, not rejected)
    // ───────────────────────────────────────────────────────────────────

    #region MC — non-GT dinucleotide is scored (not rejected) and ranks below the GT window

    /// <summary>
    /// Target "no GT": the conserved GT at window positions 3–4 is the biological consensus,
    /// NOT an input precondition — maxentpy score5 performs no GT validation (the documented
    /// weak example 'taaATAAGT' has dinucleotide "AT", not "GT", yet is itself scored to
    /// −0.116791). A window whose 3–4 dinucleotide is any other valid A/C/G/T(/U) pair is
    /// SCORED (not rejected, not crashed): the consensus GT term heavily penalises the
    /// mismatch, so the score drops far BELOW the canonical GT window
    /// (Donor_Site_Detection.md §5.3). We sweep every non-GT dinucleotide over the canonical
    /// flanks; each must score finite and strictly less than the 10.858313-bit GT window. This
    /// is the "assert the documented no-GT behavior" requirement.
    /// </summary>
    [Test]
    public void ScoreDonorMaxEnt_NonGtDinucleotide_ScoredFiniteAndBelowGtWindow()
    {
        const string prefix = "cag"; // 3 nt exon (positions 0..2)
        const string suffix = "AAGT"; // 4 nt intron tail (positions 5..8)

        foreach (char g in "ACGT")
        foreach (char t in "ACGT")
        {
            if (g == 'G' && t == 'T')
                continue; // skip the real GT — that is the consensus, scored above

            string window = prefix + g + t + suffix;
            window.Length.Should().Be(WindowLength);

            double score = ScoreDonorMaxEnt(window);

            double.IsFinite(score).Should().BeTrue(
                $"a non-GT dinucleotide ('{g}{t}') window is scored to a finite value, never NaN/Inf");
            score.Should().BeLessThan(CanonicalScore,
                $"a non-GT ('{g}{t}') donor scores below the canonical GT window (the GT consensus is penalised)");
        }
    }

    #endregion

    #endregion
}
