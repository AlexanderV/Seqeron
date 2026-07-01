using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;
using Enc = Seqeron.Genomics.IO.QualityScoreAnalyzer.QualityEncoding;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Quality area — QUALITY-PHRED-001 (Phred score handling).
/// The unit under test is the Phred ASCII ↔ score codec
/// <see cref="QualityScoreAnalyzer.ParseQualityString"/> (decode: Q = ord(char) −
/// offset, with range validation), its inverse
/// <see cref="QualityScoreAnalyzer.ToQualityString"/> (encode: char = chr(Q +
/// offset)), the cross-encoding re-offset
/// <see cref="QualityScoreAnalyzer.ConvertEncoding"/>, and the score→error-
/// probability map <see cref="QualityScoreAnalyzer.PhredToErrorProbability"/>
/// (P = 10^(−Q/10)); implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no NaN/Infinity, no nonsense output, and no *unhandled* runtime
/// exception. Every input must resolve to EITHER a well-defined, theory-correct
/// value OR a *documented, intentional* outcome — here, an
/// ArgumentOutOfRangeException for a character that decodes outside the
/// encoding's valid Phred range, and an ArgumentNullException for null input.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: QUALITY-PHRED-001 — Phred Score Handling (Quality)
/// Checklist: docs/checklists/03_FUZZING.md, row 219.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///   • MC = Malformed Content — невалідний / out-of-alphabet content.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// MAPPING of the checklist's "empty, out-of-range char, wrong offset" targets
/// onto THIS unit's documented contract
/// (docs/algorithms/Quality/Phred_Score_Handling.md):
///   • "empty" (BE) → empty quality string ⇒ empty score array; empty score
///       list ⇒ empty string. Identity boundary (§6.1, §3.2).
///   • "out-of-range char" (MC) → a character whose decoded Q falls outside the
///       encoding's valid range: a char BELOW the offset gives a negative Q
///       (Phred Q ≥ 0), or a char above ASCII 126 gives Q above the max. The doc
///       requires this be REJECTED with ArgumentOutOfRangeException, never a
///       silent negative/oversized score and never a crash (§3.3, §6.1, INV-01/02).
///   • "wrong offset" (BE) → decoding/re-encoding under a DIFFERENT encoding than
///       the data was written in. The model defines this exactly: convert is a
///       pure score-preserving re-offset (a ±31 byte shift); a Phred+33 string
///       carrying Q > 62 has NO Phred+64 representation and must raise
///       ArgumentOutOfRangeException, whereas every Phred+64 score (0–62) fits in
///       Phred+33 (0–93) and always succeeds (§2.2, §6.1, INV-04).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (docs/algorithms/Quality/Phred_Score_Handling.md)
/// Independently re-derived from Cock et al. (2010), NOT read off the code:
/// ───────────────────────────────────────────────────────────────────────────
///   • Phred score Q = −10·log₁₀(P), P = base-call error probability (§2.2).
///   • Decode: Q = ord(char) − offset; Encode: char = chr(Q + offset)  (§2.2).
///   • Phred+33 (Sanger/Illumina 1.8+): offset 33, ASCII 33–126 ⇒ Q 0–93 (INV-01).
///   • Phred+64 (Illumina 1.3+):        offset 64, ASCII 64–126 ⇒ Q 0–62 (INV-02).
///   • ToQualityString(ParseQualityString(s,e),e) == s on the valid range (INV-03).
///   • ConvertEncoding preserves the Phred score across variants (INV-04).
///   • A char decoding to Q < 0 or Q > max ⇒ ArgumentOutOfRangeException (§3.3).
///   • Null ⇒ ArgumentNullException; empty ⇒ empty (§6.1).
///   • P_error = 10^(−Q/10); finite and within (0, 1] for any valid Q ≥ 0.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class QualityPhredFuzzTests
{
    // ── Encoding constants, independently fixed from Cock et al. (2010) — NOT
    //    echoed from the implementation's private constants. ───────────────────
    private const int Phred33Offset = 33;   // ASCII '!'
    private const int Phred64Offset = 64;   // ASCII '@'
    private const int Phred33Max = 93;      // ASCII 126 '~' − 33
    private const int Phred64Max = 62;      // ASCII 126 '~' − 64

    // Independent oracle for the decode relation Q = ord(char) − offset.
    private static int ExpectedQ(char c, int offset) => c - offset;

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-PHRED-001 — Phred Score Handling (positive sanity)
    // ═════════════════════════════════════════════════════════════════════════

    // ── POSITIVE sanity: the hand-checkable worked example from the doc (§7.1).
    //    "!5?I~" under Phred+33 decodes to [0, 20, 30, 40, 93], each value
    //    re-derived here as ord(char) − 33:
    //      '!' = 33 → 0 ;  '5' = 53 → 20 ;  '?' = 63 → 30 ;
    //      'I' = 73 → 40 ;  '~' = 126 → 93 (the Phred+33 maximum).
    //    The codec round-trips back to the exact same string (INV-03). ──
    [Test]
    public void ParseQualityString_DocWorkedExample_Phred33_HandComputed()
    {
        int[] q = QualityScoreAnalyzer.ParseQualityString("!5?I~", Enc.Phred33);

        q.Should().Equal(new[] { 0, 20, 30, 40, 93 },
            "Q = ord(char) − 33 per Cock et al. (2010): !→0, 5→20, ?→30, I→40, ~→93 (§7.1)");

        // Round-trip is the identity on the valid range (INV-03).
        QualityScoreAnalyzer.ToQualityString(q, Enc.Phred33)
            .Should().Be("!5?I~", "decode∘encode is the identity on valid scores (INV-03)");
    }

    // ── POSITIVE sanity: the doc's ConvertEncoding worked example (§7.1).
    //    "@h~" written Phred+64 carries scores @=64→0, h=104→40, ~=126→62.
    //    Re-encoding those SAME scores under Phred+33 yields 0+33='!',
    //    40+33='I', 62+33='_'  ⇒  "!I_". A pure score-preserving re-offset
    //    (INV-04). ──
    [Test]
    public void ConvertEncoding_DocWorkedExample_Phred64ToPhred33_PreservesScore()
    {
        // Source scores under Phred+64, derived independently.
        QualityScoreAnalyzer.ParseQualityString("@h~", Enc.Phred64)
            .Should().Equal(new[] { 0, 40, 62 }, "Q = ord − 64: @→0, h→40, ~→62");

        string p33 = QualityScoreAnalyzer.ConvertEncoding("@h~", Enc.Phred64, Enc.Phred33);

        p33.Should().Be("!I_", "re-encoding scores 0/40/62 at offset 33 ⇒ '!','I','_' (§7.1, INV-04)");

        // The conversion must preserve the underlying Phred scores exactly.
        QualityScoreAnalyzer.ParseQualityString(p33, Enc.Phred33)
            .Should().Equal(new[] { 0, 40, 62 }, "ConvertEncoding is variant-invariant on the score (INV-04)");
    }

    // ── POSITIVE sanity: P_error = 10^(−Q/10), exact at hand-checkable points.
    //    Q0  → 10^0    = 1.0     (a base call certain to be wrong)
    //    Q10 → 10^-1   = 0.1
    //    Q20 → 10^-2   = 0.01
    //    Q30 → 10^-3   = 0.001
    //    Q40 → 10^-4   = 0.0001  (the Illumina Q40 benchmark). ──
    [TestCase(0, 1.0)]
    [TestCase(10, 0.1)]
    [TestCase(20, 0.01)]
    [TestCase(30, 0.001)]
    [TestCase(40, 0.0001)]
    public void PhredToErrorProbability_KnownPoints_MatchesFormula(int q, double expected)
    {
        QualityScoreAnalyzer.PhredToErrorProbability(q)
            .Should().BeApproximately(expected, 1e-12,
                "P_error = 10^(−Q/10) per Cock et al. (2010) §2.2");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-PHRED-001 — BE: empty input (identity boundary)
    // ═════════════════════════════════════════════════════════════════════════

    // ── BE empty: empty string ⇒ empty score array (§3.2, §6.1). ──
    [TestCase(Enc.Phred33)]
    [TestCase(Enc.Phred64)]
    [TestCase(Enc.Auto)]
    public void ParseQualityString_Empty_ReturnsEmpty(Enc enc)
    {
        QualityScoreAnalyzer.ParseQualityString(string.Empty, enc)
            .Should().BeEmpty("empty input is the identity boundary ⇒ empty output (§6.1)");
    }

    // ── BE empty: empty score list ⇒ empty string; round-trip stays empty. ──
    [TestCase(Enc.Phred33)]
    [TestCase(Enc.Phred64)]
    public void ToQualityString_Empty_ReturnsEmpty(Enc enc)
    {
        QualityScoreAnalyzer.ToQualityString(Array.Empty<int>(), enc)
            .Should().BeEmpty("encoding zero scores yields the empty quality line (§3.2)");
    }

    // ── BE empty: converting an empty string yields an empty string. ──
    [Test]
    public void ConvertEncoding_Empty_ReturnsEmpty()
    {
        QualityScoreAnalyzer.ConvertEncoding(string.Empty, Enc.Phred33, Enc.Phred64)
            .Should().BeEmpty("an empty quality line re-offsets to empty (§6.1)");
    }

    // ── BE null: null input is a documented ArgumentNullException, not a crash. ──
    [Test]
    public void ParseQualityString_Null_Throws_ArgumentNull()
    {
        Action act = () => QualityScoreAnalyzer.ParseQualityString(null!, Enc.Phred33);
        act.Should().Throw<ArgumentNullException>("null is a contract violation (§3.3, §6.1)");
    }

    [Test]
    public void ToQualityString_Null_Throws_ArgumentNull()
    {
        Action act = () => QualityScoreAnalyzer.ToQualityString(null!, Enc.Phred33);
        act.Should().Throw<ArgumentNullException>("null scores is a contract violation (§3.3)");
    }

    [Test]
    public void ConvertEncoding_Null_Throws_ArgumentNull()
    {
        Action act = () => QualityScoreAnalyzer.ConvertEncoding(null!, Enc.Phred33, Enc.Phred64);
        act.Should().Throw<ArgumentNullException>("null is a contract violation (§3.3)");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-PHRED-001 — MC: out-of-range character (malformed content)
    // ═════════════════════════════════════════════════════════════════════════

    // ── MC out-of-range, NEGATIVE Q: any char below the offset decodes to Q < 0,
    //    which is forbidden (Phred Q ≥ 0). The doc requires
    //    ArgumentOutOfRangeException, NOT a silent negative score (§6.1, INV-01).
    //    Phred+33: ' ' (ASCII 32) → −1, the just-below-offset boundary; NUL (0)
    //    and a control char are deeper out of range. ──
    [TestCase(' ',  Enc.Phred33)]   // 32 → Q −1  (boundary just below offset)
    [TestCase('\0', Enc.Phred33)]   // 0  → Q −33
    [TestCase('\t', Enc.Phred33)]   // 9  → Q −24
    [TestCase('?',  Enc.Phred64)]   // 63 → Q −1  (Phred+64 offset is 64)
    [TestCase('!',  Enc.Phred64)]   // 33 → Q −31
    public void ParseQualityString_CharBelowOffset_Throws_OutOfRange(char c, Enc enc)
    {
        // Independently confirm the char really decodes below 0 for this offset.
        int offset = enc == Enc.Phred64 ? Phred64Offset : Phred33Offset;
        ExpectedQ(c, offset).Should().BeLessThan(0, "test fixture sanity: char is below the offset");

        Action act = () => QualityScoreAnalyzer.ParseQualityString(c.ToString(), enc);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "a negative Phred score is malformed and must be rejected, not crash (§3.3, INV-01/02)");
    }

    // ── MC out-of-range, ABOVE max: a char above ASCII 126 decodes above the
    //    Phred+33 maximum (93) and must be rejected. Under Phred+64 the ceiling is
    //    lower (62) — ASCII 127+ overshoots even harder. ──
    [TestCase((char)127, Enc.Phred33)]   // 127 → Q 94  (> 93)
    [TestCase((char)200, Enc.Phred33)]   // 200 → Q 167
    [TestCase((char)127, Enc.Phred64)]   // 127 → Q 63  (> 62)
    public void ParseQualityString_CharAboveMax_Throws_OutOfRange(char c, Enc enc)
    {
        int offset = enc == Enc.Phred64 ? Phred64Offset : Phred33Offset;
        int max = enc == Enc.Phred64 ? Phred64Max : Phred33Max;
        ExpectedQ(c, offset).Should().BeGreaterThan(max, "test fixture sanity: char overshoots the max");

        Action act = () => QualityScoreAnalyzer.ParseQualityString(c.ToString(), enc);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "a Phred score above the encoding ceiling is malformed (§3.3, INV-01/02)");
    }

    // ── MC out-of-range on ENCODE: a negative or over-ceiling score is rejected
    //    symmetrically (§3.3). −1 and MaxInt are the BE extremes. ──
    [TestCase(-1, Enc.Phred33)]
    [TestCase(94, Enc.Phred33)]            // one above the Phred+33 ceiling
    [TestCase(int.MaxValue, Enc.Phred33)]  // BE MaxInt
    [TestCase(int.MinValue, Enc.Phred33)]  // BE most-negative
    [TestCase(63, Enc.Phred64)]            // one above the Phred+64 ceiling
    public void ToQualityString_ScoreOutOfRange_Throws_OutOfRange(int score, Enc enc)
    {
        Action act = () => QualityScoreAnalyzer.ToQualityString(new[] { score }, enc);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "scores outside [0, max] have no ASCII representation and must be rejected (§3.3)");
    }

    // ── MC out-of-range, one bad char in an otherwise valid run: the malformed
    //    char must still be caught (no partial-corruption / silent skip). ──
    [Test]
    public void ParseQualityString_OneBadCharAmongValid_Throws_OutOfRange()
    {
        // "II II" — the space (ASCII 32) decodes to Q −1 under Phred+33.
        Action act = () => QualityScoreAnalyzer.ParseQualityString("II II", Enc.Phred33);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "a single malformed char in the line still fails the whole parse (§3.3)");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-PHRED-001 — BE: wrong offset (cross-encoding behaviour)
    // ═════════════════════════════════════════════════════════════════════════

    // ── Wrong offset, defined SUCCESS: Phred+64 → Phred+33 always succeeds because
    //    Phred+33 range (0–93) ⊇ Phred+64 range (0–62) (§6.1, INV-04). Each
    //    converted char must carry the SAME Phred score, i.e. shift DOWN by 31. ──
    [TestCase("@",   "!")]    // Q0  : 64 → 33
    [TestCase("h",   "I")]    // Q40 : 104 → 73
    [TestCase("~",   "_")]    // Q62 : 126 → 95
    [TestCase("@h~", "!I_")]
    public void ConvertEncoding_Phred64ToPhred33_AlwaysSucceeds_ShiftsBy31(string src, string expected)
    {
        string got = QualityScoreAnalyzer.ConvertEncoding(src, Enc.Phred64, Enc.Phred33);
        got.Should().Be(expected, "Phred+64→+33 is a −31 re-offset preserving Q (INV-04)");

        // Score-preservation: decode both ways and compare.
        QualityScoreAnalyzer.ParseQualityString(got, Enc.Phred33)
            .Should().Equal(QualityScoreAnalyzer.ParseQualityString(src, Enc.Phred64),
                "the Phred score is variant-invariant (INV-04)");
    }

    // ── Wrong offset, defined OVERFLOW: a Phred+33 string carrying Q > 62 has NO
    //    Phred+64 representation and must raise ArgumentOutOfRangeException — a
    //    DEFINED detection, not a crash and not silent truncation (§6.1, INV-04).
    //    '~' = Q93 under Phred+33; 93 > 62 ⇒ unrepresentable in Phred+64. ──
    [Test]
    public void ConvertEncoding_Phred33ToPhred64_ScoreAbove62_Throws_OutOfRange()
    {
        // Independently: '~' under Phred+33 is Q 93, beyond the Phred+64 ceiling 62.
        ExpectedQ('~', Phred33Offset).Should().Be(93);
        93.Should().BeGreaterThan(Phred64Max);

        Action act = () => QualityScoreAnalyzer.ConvertEncoding("~", Enc.Phred33, Enc.Phred64);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "Q93 has no Phred+64 representation (max 62) — defined detection, not corruption (§6.1, INV-04)");
    }

    // ── Wrong offset, MISREAD: decoding a Phred+33-written string with the
    //    Phred+64 offset yields the SAME bytes shifted by −31 in score space.
    //    "II" (ASCII 73) is Q40 under +33 but Q9 under +64 (73−64). The misread is
    //    a well-defined wrong-but-not-crashing value; we pin that the off-by-31 is
    //    exactly the offset difference (no silent re-interpretation/crash). ──
    [Test]
    public void ParseQualityString_Phred33DataReadAsPhred64_OffsetDifferenceIs31()
    {
        int[] as33 = QualityScoreAnalyzer.ParseQualityString("II", Enc.Phred33);
        int[] as64 = QualityScoreAnalyzer.ParseQualityString("II", Enc.Phred64);

        as33.Should().Equal(new[] { 40, 40 }, "I = 73 − 33 = Q40");
        as64.Should().Equal(new[] { 9, 9 }, "I = 73 − 64 = Q9 — the wrong-offset misread is defined, not a crash");

        for (int i = 0; i < as33.Length; i++)
            (as33[i] - as64[i]).Should().Be(Phred64Offset - Phred33Offset,
                "the two offsets differ by exactly 31");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-PHRED-001 — Randomized boundary sweep
    // ═════════════════════════════════════════════════════════════════════════

    // ── Randomized round-trip sweep: for many random VALID quality lines under a
    //    random encoding, decode then re-encode must be the identity (INV-03),
    //    every decoded score must sit in [0, max], and every P_error must be a
    //    finite value in (0, 1]. No crash / hang / NaN / Infinity over the sweep.
    //    Locally seeded; bounded by [CancelAfter]. ──
    [Test]
    [CancelAfter(20_000)]
    public void RandomizedSweep_ValidLines_RoundTripAndFiniteProbabilities()
    {
        var rng = new Random(219_001);
        for (int iter = 0; iter < 4000; iter++)
        {
            bool phred64 = rng.Next(2) == 0;
            Enc enc = phred64 ? Enc.Phred64 : Enc.Phred33;
            int offset = phred64 ? Phred64Offset : Phred33Offset;
            int max = phred64 ? Phred64Max : Phred33Max;

            int len = rng.Next(0, 40);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = (char)(offset + rng.Next(0, max + 1));   // always in-range
            string line = new string(chars);

            int[] scores = QualityScoreAnalyzer.ParseQualityString(line, enc);
            scores.Length.Should().Be(len, "one score per character");

            foreach (int s in scores)
            {
                s.Should().BeInRange(0, max, "every decoded score is within the encoding range (INV-01/02)");
                double p = QualityScoreAnalyzer.PhredToErrorProbability(s);
                double.IsNaN(p).Should().BeFalse("P_error must never be NaN");
                double.IsInfinity(p).Should().BeFalse("P_error must never be Infinity");
                p.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(1.0,
                    "P_error = 10^(−Q/10) ∈ (0, 1] for Q ≥ 0");
            }

            QualityScoreAnalyzer.ToQualityString(scores, enc)
                .Should().Be(line, "decode∘encode is the identity on valid input (INV-03)");
        }
    }

    // ── Randomized arbitrary-char sweep: feed RANDOM 16-bit chars (incl. unicode,
    //    control, out-of-range) under a random encoding. The codec must either
    //    return a correct, in-range score array OR throw ArgumentOutOfRangeException
    //    — never a different/unhandled exception, never a silent out-of-range score,
    //    never a hang. ──
    [Test]
    [CancelAfter(20_000)]
    public void RandomizedSweep_ArbitraryChars_EitherValidOrDocumentedThrow()
    {
        var rng = new Random(219_002);
        for (int iter = 0; iter < 4000; iter++)
        {
            bool phred64 = rng.Next(2) == 0;
            Enc enc = phred64 ? Enc.Phred64 : Enc.Phred33;
            int offset = phred64 ? Phred64Offset : Phred33Offset;
            int max = phred64 ? Phred64Max : Phred33Max;

            int len = rng.Next(0, 24);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = (char)rng.Next(0, 0x10000);   // any UTF-16 code unit
            string line = new string(chars);

            bool allInRange = line.All(c => { int q = c - offset; return q >= 0 && q <= max; });

            try
            {
                int[] scores = QualityScoreAnalyzer.ParseQualityString(line, enc);
                // No throw ⇒ every char MUST have been in range, and the result correct.
                allInRange.Should().BeTrue("a successful parse implies every char decoded in range");
                scores.Length.Should().Be(len);
                for (int i = 0; i < len; i++)
                {
                    int expected = ExpectedQ(line[i], offset);
                    scores[i].Should().Be(expected, "Q = ord(char) − offset (§2.2)");
                    scores[i].Should().BeInRange(0, max);
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                // Documented outcome: at least one char was out of range.
                allInRange.Should().BeFalse("an out-of-range throw implies some char decoded out of range");
            }
        }
    }

    #endregion
}
