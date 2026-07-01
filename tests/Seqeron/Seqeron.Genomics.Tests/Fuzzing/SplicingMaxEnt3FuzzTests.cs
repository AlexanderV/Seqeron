using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Splicing area — MaxEntScan <c>score3ss</c> maximum-entropy
/// 3' (acceptor) splice-site model (SPLICE-MAXENT3-001).
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
/// For MaxEntScan score3ss the headline hazards a fuzzer must rule out are:
///   • IndexOutOfRange — from a window whose length ≠ 23 reaching the fixed AG-position
///     (18..19) reads or the 21-nt "rest" Span slices. The contract REJECTS any wrong
///     length up front with an ArgumentException (NOT an IndexOutOfRange).
///   • KeyNotFoundException — from a non-A/C/G/T(/U) character producing an out-of-table
///     hash, or from the AG-position consensus/background dictionary lookups. The contract
///     REJECTS any non-A/C/G/T(/U) character with an ArgumentException before any lookup.
///   • NaN / ±Infinity — from log2(0) when a maxent sub-sequence probability is zero. The
///     score is log2(P_maxent/P_background); for ANY all-A/C/G/T(/U) 23-mer every embedded
///     table entry is strictly positive, so the product is positive and the log2 is FINITE.
///     A NaN/Inf on a valid 23-mer would be a real bug.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The score3ss contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The 3' splice ACCEPTOR site marks the END of an intron; its consensus is a
/// polypyrimidine tract followed by the almost-invariant terminal AG. MaxEntScan score3ss
/// (Yeo &amp; Burge 2004) models a fixed 23-nt window — 20 intronic + 3 exonic nt, with the
/// conserved AG at 0-based window positions 18–19 — as
///   score = log2( P_maxent(window) / P_background(window) )   (in bits),
/// computed by removing the AG (scored by a consensus/background dinucleotide model) and
/// factorising the remaining 21 ("rest") positions over nine overlapping sub-sequences
/// (five multiplied, four divided — inclusion/exclusion) via embedded probability tables.
/// score3ss is CALIBRATED so true 3' sites score positive: the canonical documented example
/// score3('ttccaaacgaacttttgtAGgga') == 2.89 bits; a strong site 'tgtctttttctgtgtggcAGtgg'
/// == 8.19; a weak site 'ttctctcttcagacttatAGcaa' == −0.08.
/// — docs/algorithms/Splicing/Acceptor_Site_Detection.md §5.3 (MaxEntScan paragraph), §6.2;
///   src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs (#region MaxEntScan score3ss).
///
/// Method under test (src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs ~L1176):
///   double ScoreAcceptorMaxEnt(string window)
///     → the score3ss in bits.
///
/// Documented input handling (XML on ScoreAcceptorMaxEnt; §5.3; §6.1 "Lowercase input"):
///   • null window               → ArgumentNullException.
///   • window.Length ≠ 23        → ArgumentException (NOT IndexOutOfRange).
///   • non-A/C/G/T(/U) character → ArgumentException (NOT KeyNotFound / IndexOutOfRange).
///   • lowercase / mixed case    → folded via ToUpperInvariant(), then scored — NOT rejected.
///   • T and U are equivalent (T==U in both the hash encoding and the AG model).
///   • The AG at 18–19 is NOT *required* by the math: a non-AG window is SCORED (the
///     consensus AG term heavily penalises it → strongly negative), never rejected. The AG
///     is only the biological consensus, not an input precondition.
///
/// Theory-correct invariants asserted (independently derived from §5.3 / Yeo &amp; Burge 2004):
///   • FINITE — every all-A/C/G/T(/U) 23-mer yields a finite score (no NaN, no ±Inf): the
///     maxent factorisation is a product of strictly-positive table probabilities, so its
///     log2 is finite. Pinned over hundreds of LOCALLY-seeded random 23-mers AND every
///     A/C/G/T homopolymer.
///   • DOCUMENTED ANCHORS — the published maxentpy reference values (the independent oracle,
///     recorded verbatim in the source/spec — NOT read off the code's own tables):
///       score3('ttccaaacgaacttttgtAGgga') ≈ 2.886773  (canonical, "2.89 bits")
///       score3('tgtctttttctgtgtggcAGtgg') ≈ 8.190965  (strong site)
///       score3('ttctctcttcagacttatAGcaa') ≈ −0.080278 (weak site)
///   • CONSENSUS &gt; RANDOM — a strong consensus acceptor (poly-pyrimidine tract + AG) scores
///     ABOVE a random / GC-rich 23-mer; score3ss is calibrated so true sites score high.
///   • CASE / T==U INVARIANCE — lowercase, uppercase and mixed-case spellings of the same
///     window, and the T-vs-U spelling, all produce the identical score.
///   • AG-NOT-REQUIRED — replacing the 18–19 AG with a non-AG dinucleotide does NOT throw and
///     drops the score far below the canonical AG window (the model scores it, never rejects).
///
/// Discipline: tests encode the DOCUMENTED theory, derived independently; a test that
/// passes against a wrong implementation is invalid. If a source-derived test and the code
/// disagree the CODE is wrong (fixed minimally per the doc). A NaN/Inf score or an
/// IndexOutOfRange/KeyNotFound on a malformed window would be a REAL bug.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SplicingMaxEnt3FuzzTests
{
    #region SPLICE-MAXENT3-001 — MaxEntScan score3ss 3' acceptor maximum-entropy model

    // The canonical documented MaxEntScan score3 examples and their published full-precision
    // maxentpy reference values (the independent oracle; Acceptor_Site_Detection.md §5.3 and
    // the SpliceSitePredictor source — NOT read off the embedded tables).
    private const string CanonicalWindow = "ttccaaacgaacttttgtAGgga"; // 23 nt, AG at 18-19
    private const double CanonicalScore = 2.886773;                    // "2.89 bits"
    private const string StrongWindow = "tgtctttttctgtgtggcAGtgg";
    private const double StrongScore = 8.190965;
    private const string WeakWindow = "ttctctcttcagacttatAGcaa";
    private const double WeakScore = -0.080278;

    private const int WindowLength = 23;

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz windows are reproducible.</summary>
    private static string RandomWindow(int length, int seed, string alphabet = "ACGT")
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity — the documented anchors (proves we score a real model)
    // ───────────────────────────────────────────────────────────────────

    #region Anchor — documented maxentpy reference values are reproduced exactly

    /// <summary>
    /// Positive sanity / anchor: the three published MaxEntScan score3 reference windows
    /// reproduce their documented full-precision values exactly (2.886773 / 8.190965 /
    /// −0.080278; Acceptor_Site_Detection.md §5.3, Yeo &amp; Burge 2004 / maxentpy). This pins
    /// the fuzz harness to a detector that actually computes the real maximum-entropy score
    /// — not a no-op — so the boundary/malformed assertions below are meaningful. The values
    /// are the INDEPENDENT oracle (the maxentpy reference), never read off the code's tables.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_DocumentedReferenceWindows_ReproduceExactValues()
    {
        ScoreAcceptorMaxEnt(CanonicalWindow).Should().BeApproximately(CanonicalScore, 1e-6,
            "the canonical documented window scores 2.886773 bits (maxentpy score3 oracle)");
        ScoreAcceptorMaxEnt(StrongWindow).Should().BeApproximately(StrongScore, 1e-6,
            "the strong-site documented window scores 8.190965 bits");
        ScoreAcceptorMaxEnt(WeakWindow).Should().BeApproximately(WeakScore, 1e-6,
            "the weak-site documented window scores −0.080278 bits (a calibrated near-zero site)");
    }

    /// <summary>
    /// CONSENSUS &gt; RANDOM anchor (Acceptor_Site_Detection.md §5.3 — score3ss is calibrated
    /// so true 3' sites score positive). A strong consensus acceptor (a long polypyrimidine
    /// tract immediately upstream of the conserved AG) must score strictly ABOVE a random and
    /// a GC-rich 23-mer of the same length. We sweep LOCALLY-seeded random windows so the
    /// ordering is pinned against arbitrary backgrounds, not one hand-picked decoy.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_ConsensusAcceptor_RanksAboveRandomAndGcRich()
    {
        // Strong consensus 3' site: 18-nt polypyrimidine tract + AG + 3-nt exon.
        const string consensus = "uuuuucuuuuccuuuucuAGgua"; // 23 nt, AG at 18-19
        double consensusScore = ScoreAcceptorMaxEnt(consensus);

        // A consensus PPT+AG site is a strong, positive-scoring acceptor.
        consensusScore.Should().BeGreaterThan(0.0,
            "a polypyrimidine-tract + AG consensus is a calibrated true 3' site → positive bits");

        // GC-rich decoy with the same terminal AG: a GC tract is NOT a polypyrimidine tract,
        // so the maximum-entropy model scores it far below the consensus.
        double gcRich = ScoreAcceptorMaxEnt("gcgcgcgcgcgcgcgcgcAGgcg");
        consensusScore.Should().BeGreaterThan(gcRich,
            "a true PPT+AG acceptor outscores a GC-rich 23-mer (score3ss favours pyrimidine tracts)");

        // Random A/C/G/T 23-mers (always with the AG at 18-19 so only the background varies):
        // the consensus PPT site must outrank every random background.
        foreach (int seed in new[] { 7, 41, 137, 911, 2026 })
        {
            char[] w = RandomWindow(WindowLength, seed).ToCharArray();
            w[18] = 'A';
            w[19] = 'G';
            string random = new string(w);
            double randomScore = ScoreAcceptorMaxEnt(random);
            consensusScore.Should().BeGreaterThan(randomScore,
                $"the PPT+AG consensus must outscore the random background (seed {seed}, '{random}')");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  FINITE — every valid 23-mer yields a finite score (no NaN / ±Inf)
    // ───────────────────────────────────────────────────────────────────

    #region BE/MC — finite score on any A/C/G/T(/U) 23-mer (no log2(0) → no NaN/Inf)

    /// <summary>
    /// The headline numeric invariant: for ANY all-A/C/G/T(/U) 23-mer the score is FINITE —
    /// never NaN, never ±Infinity. The score is log2(P_maxent/P_background); the maxent
    /// factorisation is a product of strictly-positive embedded table probabilities (five
    /// multiplied, four divided) times a strictly-positive consensus/background AG term, so
    /// the argument of log2 is strictly positive and the result is finite (Yeo &amp; Burge 2004;
    /// SpliceSitePredictor.ScoreAcceptorMaxEnt). We sweep hundreds of LOCALLY-seeded random
    /// 23-mers across A/C/G/T spellings; a NaN/Inf on any valid window would be a real bug.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_AnyAcgtWindow_ProducesFiniteScore()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            string window = RandomWindow(WindowLength, seed);
            window.Length.Should().Be(WindowLength, "the fuzz fixture builds an exactly-23-nt window");

            double score = ScoreAcceptorMaxEnt(window);

            double.IsNaN(score).Should().BeFalse($"score3 of an A/C/G/T 23-mer must not be NaN (seed {seed}, '{window}')");
            double.IsInfinity(score).Should().BeFalse($"score3 of an A/C/G/T 23-mer must not be ±Inf (seed {seed}, '{window}')");
            double.IsFinite(score).Should().BeTrue($"score3 of an A/C/G/T 23-mer is finite (seed {seed}, '{window}')");
        }
    }

    /// <summary>
    /// The four homopolymer corners (AAA…/CCC…/GGG…/TTT…) are the extreme low-information
    /// windows — none contains the AG consensus, all four are valid A/C/G/T(/U) 23-mers, and
    /// each must still produce a FINITE score (no zero-probability sub-sequence anywhere in
    /// the factorisation). They are scored, never rejected — the AG is not an input
    /// precondition. The U-spelling homopolymer must match its T-spelling exactly (T==U).
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_Homopolymers_ScoredAndFinite()
    {
        foreach (char b in new[] { 'A', 'C', 'G', 'T' })
        {
            string window = new string(b, WindowLength);
            double score = ScoreAcceptorMaxEnt(window);
            double.IsFinite(score).Should().BeTrue($"the {b}-homopolymer 23-mer scores a finite value, not NaN/Inf");
        }

        // T- and U-homopolymers are the same window after T→U normalisation → identical score.
        ScoreAcceptorMaxEnt(new string('U', WindowLength))
            .Should().Be(ScoreAcceptorMaxEnt(new string('T', WindowLength)),
                "T==U: the U-homopolymer scores identically to the T-homopolymer");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — window length ≠ 23 → documented ArgumentException, NOT IndexOutOfRange
    // ───────────────────────────────────────────────────────────────────

    #region BE — wrong-length and empty windows throw ArgumentException (not IndexOutOfRange)

    /// <summary>
    /// Target "window ≠ 23 nt" / "empty seq" (KEY boundary). The length is validated up front:
    /// ANY window whose length ≠ 23 — including the empty string, lengths just below/above 23,
    /// and a window that merely PROVISIONALLY looks acceptor-like — throws an ArgumentException
    /// (a DOCUMENTED validation outcome), NOT an IndexOutOfRangeException off the fixed
    /// AG-position reads (18..19) or the 21-nt rest Span slices. We sweep lengths 0..50 except
    /// 23, building each from a real A/C/G/T alphabet so the rejection is on LENGTH, not
    /// alphabet. — ScoreAcceptorMaxEnt XML (&lt;exception ArgumentException&gt; "not exactly 23 nt").
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_WindowLengthNot23_ThrowsArgumentException()
    {
        for (int len = 0; len <= 50; len++)
        {
            if (len == WindowLength)
                continue;

            string window = RandomWindow(len, seed: 1000 + len);
            window.Length.Should().Be(len);

            Action act = () => ScoreAcceptorMaxEnt(window);
            act.Should().Throw<ArgumentException>(
                    $"a {len}-nt window (≠ 23) is rejected up front, never crashing with IndexOutOfRange")
                .Which.Should().NotBeOfType<ArgumentNullException>(
                    "a non-null wrong-length window is an ArgumentException, not a null-argument exception");
        }
    }

    /// <summary>
    /// The empty and whitespace-padded boundaries explicitly: "" (length 0) and a window that
    /// is 23 visible nt but padded with whitespace to a wrong length both throw ArgumentException
    /// on length. An empty window must NOT slip through to read window[18] (IndexOutOfRange);
    /// padding makes a *would-be* valid window the wrong length, which length-validation catches
    /// BEFORE the alphabet check ever sees the whitespace.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_EmptyAndPaddedWindows_ThrowArgumentException()
    {
        foreach (string bad in new[]
                 {
                     "",                              // empty → length 0
                     "   ",                           // whitespace only → length 3
                     " ttccaaacgaacttttgtAGgga ",     // 23-nt window padded to length 25
                     "ttccaaacgaacttttgtAGgga ",      // trailing space → length 24
                 })
        {
            Action act = () => ScoreAcceptorMaxEnt(bad);
            act.Should().Throw<ArgumentException>(
                $"a length-{bad.Length} window (≠ 23) is rejected on length, never IndexOutOfRange ('{bad}')");
        }
    }

    /// <summary>
    /// null throws ArgumentNullException (the documented null contract), distinct from the
    /// wrong-length ArgumentException above — and never a NullReferenceException leaking from
    /// an unchecked .Length / .ToUpperInvariant() access.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_NullWindow_ThrowsArgumentNullException()
    {
        Action act = () => ScoreAcceptorMaxEnt(null!);
        act.Should().Throw<ArgumentNullException>("a null window is the documented null contract");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — non-A/C/G/T(/U) characters → documented ArgumentException
    // ───────────────────────────────────────────────────────────────────

    #region MC — non-ACGT characters throw ArgumentException (not KeyNotFound/IndexOutOfRange)

    /// <summary>
    /// Target "non-ACGT" (Malformed Content). A length-correct (23-nt) window that contains a
    /// character outside {A,C,G,T,U} — IUPAC ambiguity codes (N, R, Y…), digits, punctuation —
    /// throws an ArgumentException (the DOCUMENTED handling), NOT a KeyNotFoundException from an
    /// out-of-table hash and NOT an IndexOutOfRange. The invalid character may sit in the AG
    /// positions (consensus/background dictionary lookup) OR anywhere in the 21-nt rest
    /// (HashMaxEntSubsequence) — both paths must reject cleanly. — ScoreAcceptorMaxEnt XML
    /// (&lt;exception ArgumentException&gt; "contains a non-A/C/G/T(/U) character").
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_NonAcgtCharacter_ThrowsArgumentException()
    {
        foreach (string bad in new[]
                 {
                     "ttccaaacgaacttttgtAGggN",   // N in the exon — rest hash
                     "Nttccaaacgaacttttgt AGg".Replace(" ", "N"), // keep length 23 with leading junk
                     "ttccaaacgaacttttgtNGgga",   // N at AG position 18 (the A) — consensus lookup
                     "ttccaaacgaacttttgtANgga",   // N at AG position 19 (the G) — consensus lookup
                     "ttccaaacgaacttttgtAGgg!",   // punctuation
                     "ttccaaacgaacttttgtAGgg5",   // a digit
                     "ttccaaacgaacttttgtRYgga",   // IUPAC ambiguity codes at the AG positions
                     "rrrrrrrrrrrrrrrrrrAGgga",    // an ambiguity-code tract in the intron
                 })
        {
            bad.Length.Should().Be(WindowLength, $"the malformed-content fixture stays length 23 ('{bad}')");

            Action act = () => ScoreAcceptorMaxEnt(bad);
            act.Should().Throw<ArgumentException>(
                $"a non-A/C/G/T(/U) character is rejected, never KeyNotFound/IndexOutOfRange ('{bad}')");
        }
    }

    /// <summary>
    /// An all-ambiguity ("NNN…") 23-nt window — length-correct but with no concrete nucleotide
    /// anywhere — throws ArgumentException, never wedging or leaking a KeyNotFound from the very
    /// first sub-sequence hash. This is the malformed-content corner where every position is junk.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_AllAmbiguityWindow_ThrowsArgumentException()
    {
        string window = new string('N', WindowLength);
        Action act = () => ScoreAcceptorMaxEnt(window);
        act.Should().Throw<ArgumentException>(
            "an all-N 23-mer has no valid nucleotide and is rejected, not crashed");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — lowercase / mixed case is folded (documented), not rejected
    // ───────────────────────────────────────────────────────────────────

    #region MC — lowercase/mixed-case folding and T==U equivalence

    /// <summary>
    /// Target "lowercase" (documented: case-fold to ACGT, not reject — ScoreAcceptorMaxEnt
    /// uppercases via ToUpperInvariant; Acceptor_Site_Detection.md §6.1 "Lowercase input").
    /// The fully-lowercase, fully-uppercase and a deliberately mixed-case spelling of the SAME
    /// 23-mer all produce the IDENTICAL score — case carries no information. We pin this on the
    /// three documented reference windows AND on LOCALLY-seeded random windows so case-folding
    /// holds for arbitrary content, not just the anchors.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_CaseFolding_AllSpellingsScoreIdentically()
    {
        var windows = new List<string> { CanonicalWindow, StrongWindow, WeakWindow };
        for (int seed = 50; seed < 60; seed++)
        {
            char[] w = RandomWindow(WindowLength, seed).ToCharArray();
            w[18] = 'A';
            w[19] = 'G';
            windows.Add(new string(w));
        }

        foreach (string window in windows)
        {
            string lower = window.ToLowerInvariant();
            string upper = window.ToUpperInvariant();
            // Alternating-case spelling: same letters, scrambled case.
            char[] mixedChars = lower.ToCharArray();
            for (int i = 0; i < mixedChars.Length; i += 2)
                mixedChars[i] = char.ToUpperInvariant(mixedChars[i]);
            string mixed = new string(mixedChars);

            double lowerScore = ScoreAcceptorMaxEnt(lower);
            double upperScore = ScoreAcceptorMaxEnt(upper);
            double mixedScore = ScoreAcceptorMaxEnt(mixed);

            upperScore.Should().Be(lowerScore,
                $"case is folded: '{upper}' must score identically to '{lower}'");
            mixedScore.Should().Be(lowerScore,
                $"mixed case is folded: '{mixed}' must score identically to '{lower}'");
            double.IsFinite(lowerScore).Should().BeTrue("a folded valid window still scores finite");
        }
    }

    /// <summary>
    /// T==U equivalence (the model treats DNA T and RNA U identically in both the hash encoding
    /// and the AG consensus model). The T-spelling and the U-spelling of the same windows — in
    /// any case — score identically, across the documented anchors and random windows. This is
    /// the alphabet-equivalence boundary: a fuzzer must not see U treated as out-of-alphabet junk.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_TandUSpellings_ScoreIdentically()
    {
        var windows = new List<string> { CanonicalWindow, StrongWindow, WeakWindow };
        for (int seed = 70; seed < 80; seed++)
        {
            char[] w = RandomWindow(WindowLength, seed).ToCharArray();
            w[18] = 'A';
            w[19] = 'G';
            windows.Add(new string(w));
        }

        foreach (string window in windows)
        {
            string tForm = window.ToUpperInvariant();                 // DNA spelling
            string uForm = tForm.Replace('T', 'U');                   // RNA spelling
            string uLower = tForm.ToLowerInvariant().Replace('t', 'u'); // lowercase RNA

            double tScore = ScoreAcceptorMaxEnt(tForm);
            ScoreAcceptorMaxEnt(uForm).Should().Be(tScore,
                $"T==U: the U-form '{uForm}' scores identically to the T-form '{tForm}'");
            ScoreAcceptorMaxEnt(uLower).Should().Be(tScore,
                $"T==U and case-folding compose: lowercase RNA '{uLower}' scores identically");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  MC — the AG at 18-19 is NOT required by the math (scored, not rejected)
    // ───────────────────────────────────────────────────────────────────

    #region MC — non-AG dinucleotide is scored (not rejected) and ranks below the AG window

    /// <summary>
    /// Target "no AG": the conserved AG at window positions 18–19 is the biological consensus,
    /// NOT an input precondition — maxentpy score3 performs no AG validation. A window whose
    /// 18–19 dinucleotide is any other valid A/C/G/T(/U) pair is SCORED (not rejected, not
    /// crashed): the consensus AG term heavily penalises the mismatch, so the score drops far
    /// BELOW the canonical AG window (Acceptor_Site_Detection.md §5.3). We sweep every non-AG
    /// dinucleotide over the canonical flanks; each must score finite and strictly less than the
    /// 2.886773-bit AG window. This is the "assert the documented no-AG behavior" requirement.
    /// </summary>
    [Test]
    public void ScoreAcceptorMaxEnt_NonAgDinucleotide_ScoredFiniteAndBelowAgWindow()
    {
        const string prefix = "ttccaaacgaacttttgt"; // 18 nt intron (positions 0..17)
        const string suffix = "gga";                 // 3 nt exon (positions 20..22)

        foreach (char a in "ACGT")
        foreach (char g in "ACGT")
        {
            if (a == 'A' && g == 'G')
                continue; // skip the real AG — that is the consensus, scored above

            string window = prefix + a + g + suffix;
            window.Length.Should().Be(WindowLength);

            double score = ScoreAcceptorMaxEnt(window);

            double.IsFinite(score).Should().BeTrue(
                $"a non-AG dinucleotide ('{a}{g}') window is scored to a finite value, never NaN/Inf");
            score.Should().BeLessThan(CanonicalScore,
                $"a non-AG ('{a}{g}') acceptor scores below the canonical AG window (the AG consensus is penalised)");
        }
    }

    #endregion

    #endregion
}
