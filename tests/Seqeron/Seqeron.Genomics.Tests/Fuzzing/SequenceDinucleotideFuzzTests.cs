using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Statistics-area dinucleotide-frequency unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception. A raw runtime exception or a
/// hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-DINUC-001 — Dinucleotide analysis (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 122.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — length&lt;2 (empty / single base), all-N,
///          homopolymer (plus the lower-size and degenerate boundaries around
///          the N−1 dinucleotide count and the zero-total denominator).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The dinucleotide contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateDinucleotideFrequencies(string)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 602–631), returning IReadOnlyDictionary&lt;string,double&gt;.
///
/// Documented behaviour (Dinucleotide_Analysis.md, Test Unit ID SEQ-DINUC-001):
///   • §2.2 / §5.1: f_XY = count(XY) / (N − 1) — normalized frequency over the
///     N−1 ADJACENT (overlapping) dinucleotide positions of a length-N sequence
///     (Karlin genomic-signature N−1 convention).
///   • §3.3 / §6.1: null / empty / length &lt; 2 → EMPTY dictionary. With fewer
///     than 2 bases there are N−1 ≤ 0 dinucleotide positions, so no dinucleotide
///     exists and there is no zero-total denominator to divide by.
///   • §3.3 / §4.2: only dinucleotides over the alphabet {A,T,G,C,U} are counted;
///     any pair containing a non-alphabet base (N, degenerate IUPAC, junk) is
///     EXCLUDED from the counts and the N−1 total. Hence an all-N input yields the
///     EMPTY dictionary ("NN" is not an alphabet dinucleotide), with no crash and
///     no DivideByZero (the frequency loop never runs when total = 0).
///   • §3.3: input is upper-cased (ToUpperInvariant) before counting, so counting
///     is case-insensitive.
///
/// Invariants pinned (Dinucleotide_Analysis.md §2.4):
///   • INV-01: dinucleotide frequencies sum to 1.0 whenever ≥ 1 valid dinucleotide
///     was counted (each counted position contributes 1/total).
///   • INV-04: every returned frequency is finite and ≥ 0 (count ≥ 1, total ≥ 1).
///   • N−1 count identity: when EVERY adjacent pair is alphabet-valid, the counts
///     sum to exactly Length − 1 (no off-by-one).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceDinucleotideFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The dinucleotide alphabet: a pair counts only if BOTH bases are here.</summary>
    private const string DinucAlphabet = "ATGCU";

    /// <summary>Generates a random string of arbitrary BMP code points (0x0000–0xFFFF) —
    /// control chars, the null byte, lone surrogate halves, unicode letters/digits:
    /// random-byte fuzz fodder for the dinucleotide counter.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal well-formedness contract that must hold for ANY input:
    /// every key is a length-2 alphabet dinucleotide, every frequency is finite,
    /// strictly positive and ≤ 1, and — per INV-01 — the frequencies sum to 1.0
    /// whenever the dictionary is non-empty (i.e. ≥ 1 valid dinucleotide).
    /// </summary>
    private static void AssertWellFormed(IReadOnlyDictionary<string, double> freq)
    {
        freq.Should().NotBeNull("the method must always return a dictionary, never null");

        foreach (var (dinuc, f) in freq)
        {
            dinuc.Should().HaveLength(2, "a dinucleotide key is exactly two bases");
            dinuc.All(c => DinucAlphabet.Contains(c))
                .Should().BeTrue($"'{dinuc}' must be over the {{A,T,G,C,U}} alphabet");

            double.IsFinite(f).Should().BeTrue($"frequency of '{dinuc}' must be finite (INV-04)");
            f.Should().BeGreaterThan(0.0, $"a counted dinucleotide '{dinuc}' has count ≥ 1 (INV-04)");
            f.Should().BeLessThanOrEqualTo(1.0 + Tolerance,
                $"a frequency '{dinuc}' is count/total ≤ 1");
        }

        if (freq.Count > 0)
            freq.Values.Sum().Should().BeApproximately(1.0, Tolerance,
                "INV-01: frequencies over the counted set sum to 1");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-DINUC-001 — dinucleotide frequency : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact result

    /// <summary>
    /// Positive baseline (not a boundary): a known sequence must yield the documented
    /// dinucleotide counts EXACTLY, with the total equal to N−1 and fractions summing
    /// to 1. "ATGCGCGT" (N=8) has 7 adjacent pairs: AT, TG, GC, CG, GC, CG, GT ⇒
    /// GC×2, CG×2, AT×1, TG×1, GT×1 (total 7). So f_GC = f_CG = 2/7 and
    /// f_AT = f_TG = f_GT = 1/7. Confirms the suite asserts the BUSINESS contract.
    /// — Dinucleotide_Analysis.md §7.1 worked example (ATGCGCGT, 7 dinuc positions).
    /// </summary>
    [Test]
    public void Dinuc_KnownSequence_MatchesHandComputedCountsAndFractions()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        // Five distinct dinucleotides over 7 positions (N−1 = 8−1).
        freq.Should().HaveCount(5);
        freq["GC"].Should().BeApproximately(2.0 / 7.0, Tolerance);
        freq["CG"].Should().BeApproximately(2.0 / 7.0, Tolerance);
        freq["AT"].Should().BeApproximately(1.0 / 7.0, Tolerance);
        freq["TG"].Should().BeApproximately(1.0 / 7.0, Tolerance);
        freq["GT"].Should().BeApproximately(1.0 / 7.0, Tolerance);

        // N−1 count identity and INV-01.
        freq.Values.Sum().Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: length < 2 (empty / null / single base)

    /// <summary>
    /// BE: the empty string is the lower size boundary. With N=0 there are N−1 = −1
    /// (≤ 0) dinucleotide positions, so the documented result is the EMPTY dictionary —
    /// NO DivideByZero on a zero total, no exception.
    /// — Dinucleotide_Analysis.md §3.3 / §6.1 (null/empty/length&lt;2 → empty dictionary).
    /// </summary>
    [Test]
    public void Dinuc_EmptyString_IsEmptyAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateDinucleotideFrequencies(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        var freq = act();
        freq.Should().BeEmpty("N−1 ≤ 0 ⇒ no dinucleotide positions");
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 607) — empty dictionary, no NullReferenceException.
    /// — Dinucleotide_Analysis.md §3.3 (null/empty → empty dictionary).
    /// </summary>
    [Test]
    public void Dinuc_Null_IsEmptyAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateDinucleotideFrequencies(null!);

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().BeEmpty();
    }

    /// <summary>
    /// BE: a one-base sequence (N=1) is the critical length&lt;2 boundary — N−1 = 0
    /// dinucleotide positions, so the result is EMPTY and the zero-total denominator
    /// is never reached (guards against an off-by-one that would form a 1-char
    /// "dinucleotide" or divide by a zero total). Every canonical base behaves so.
    /// — Dinucleotide_Analysis.md §6.1 (length&lt;2 → empty dictionary).
    /// </summary>
    [TestCase("A")]
    [TestCase("C")]
    [TestCase("G")]
    [TestCase("T")]
    [TestCase("U")]
    [TestCase("N")]
    [TestCase("a")]
    public void Dinuc_SingleBase_IsEmpty_NoDinucleotideExists(string seq)
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

        freq.Should().BeEmpty("a length-1 sequence has N−1 = 0 dinucleotide positions");
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: all-N

    /// <summary>
    /// BE: every base is 'N'. The only adjacent pair is "NN", which is NOT over the
    /// {A,T,G,C,U} alphabet, so EVERY pair is excluded: the total is 0 and the
    /// documented result is the EMPTY dictionary — the frequency-conversion loop never
    /// runs, so there is NO DivideByZero by the zero total. Holds at the length=2
    /// boundary and at scale, and case-insensitively for lowercase 'n'.
    /// — Dinucleotide_Analysis.md §3.3 / §4.2 (non-alphabet dinucleotides excluded).
    /// </summary>
    [TestCase("NN")]
    [TestCase("NNN")]
    [TestCase("nnnn")]
    public void Dinuc_AllN_IsEmpty_NoDivideByZero(string seq)
    {
        var act = () => SequenceStatistics.CalculateDinucleotideFrequencies(seq);

        act.Should().NotThrow("an all-N input must never throw (no zero-total division)");
        act().Should().BeEmpty("'NN' is not an {A,T,G,C,U} dinucleotide ⇒ excluded");
    }

    /// <summary>
    /// BE: a long all-N homopolymer — the all-N boundary at scale. Still empty, still
    /// no DivideByZero, and bounded in time.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Dinuc_LongAllN_IsEmpty_NoCrash()
    {
        string seq = new string('N', 5000);

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

        freq.Should().BeEmpty("no NN pair is over the dinucleotide alphabet");
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: homopolymer

    /// <summary>
    /// BE: a homopolymer (e.g. "AAAA…") is the degenerate single-dinucleotide boundary.
    /// A length-L run of base X has L−1 adjacent "XX" pairs and NO other pair, so the
    /// result is a SINGLE entry {"XX": 1.0}: f_XX = (L−1)/(L−1) = 1.0 and 0 for every
    /// other dinucleotide. Pins the N−1 count (off-by-one guard) and INV-01 at the
    /// extreme. Verified for each canonical base and several lengths from the minimum.
    /// — Dinucleotide_Analysis.md §2.2 (f_XY = count/(N−1)); §2.4 INV-01.
    /// </summary>
    [TestCase('A', 2)]
    [TestCase('A', 5)]
    [TestCase('C', 4)]
    [TestCase('G', 10)]
    [TestCase('T', 3)]
    [TestCase('U', 7)]
    public void Dinuc_Homopolymer_IsSingleDinucleotideWithFractionOne(char baseChar, int length)
    {
        string seq = new string(baseChar, length);
        string expectedKey = $"{baseChar}{baseChar}";

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

        freq.Should().HaveCount(1, "a homopolymer has exactly one distinct dinucleotide");
        freq.Should().ContainKey(expectedKey);
        freq[expectedKey].Should().BeApproximately(1.0, Tolerance,
            "the only dinucleotide accounts for all N−1 positions");

        // INV-01 / fraction over the single counted dinucleotide.
        freq.Values.Sum().Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: lowercase homopolymer must produce the SAME upper-cased dinucleotide key —
    /// guards against a case-sensitivity bug (ToUpperInvariant, line 610).
    /// </summary>
    [Test]
    public void Dinuc_LowercaseHomopolymer_KeyIsUpperCased()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("gggg");

        freq.Should().HaveCount(1);
        freq.Should().ContainKey("GG", "input is upper-cased before counting");
        freq["GG"].Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — degenerate / junk interleaving

    /// <summary>
    /// BE: a pair is counted ONLY if BOTH bases are over {A,T,G,C,U}. A canonical base
    /// adjacent to N (or junk) forms a non-alphabet pair that is excluded, which also
    /// breaks the adjacency chain. "ANA" (N=3) has pairs AN, NA — both excluded — so
    /// the result is EMPTY even though two canonical bases are present. "AANAA" yields
    /// only the AA pairs (positions 0 and 3): {"AA": 1.0}. No crash, INV-01 holds.
    /// — Dinucleotide_Analysis.md §3.3 (dinucleotides outside the alphabet excluded).
    /// </summary>
    [Test]
    public void Dinuc_CanonicalSeparatedByN_PairsExcluded()
    {
        SequenceStatistics.CalculateDinucleotideFrequencies("ANA")
            .Should().BeEmpty("AN and NA are not {A,T,G,C,U} dinucleotides");

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("AANAA");
        freq.Should().HaveCount(1);
        freq["AA"].Should().BeApproximately(1.0, Tolerance,
            "only the two AA pairs are alphabet-valid; AN/NA excluded");
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Random / RB fuzz: never throws, always well-formed

    /// <summary>
    /// BE/RB: a large batch of arbitrary BMP strings (control chars, null byte, lone
    /// surrogate halves, unicode letters/digits) must NEVER throw and must ALWAYS
    /// produce a well-formed dictionary satisfying every invariant. Core fuzz
    /// guarantee: no DivideByZero, no IndexOutOfRange, no overflow on garbage.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Dinuc_RandomGarbageStrings_NeverThrow_AlwaysWellFormed()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            IReadOnlyDictionary<string, double> freq = null!;
            var act = () => freq = SequenceStatistics.CalculateDinucleotideFrequencies(input);

            act.Should().NotThrow($"garbage input (len {len}) must never crash dinuc counting");
            AssertWellFormed(freq);
        }
    }

    /// <summary>
    /// BE: a randomly built canonical-only sequence (N ≥ 2, alphabet {A,T,G,C,U}) has
    /// EVERY adjacent pair alphabet-valid, so the counts must sum to exactly N−1 (the
    /// N−1 identity / off-by-one guard) and the frequencies must equal an independent
    /// re-count over all N−1 sliding pairs, summing to 1. Cross-checks the counting
    /// logic against a simple oracle over many shapes.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Dinuc_RandomCanonicalSequences_CountsMatchSlidingOracle()
    {
        var rng = new Random(424242);

        for (int iteration = 0; iteration < 1000; iteration++)
        {
            int len = rng.Next(2, 300);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = DinucAlphabet[rng.Next(DinucAlphabet.Length)];
            string seq = new string(chars);

            // Oracle: every one of the N−1 overlapping pairs is counted (all canonical).
            var oracle = new Dictionary<string, int>();
            for (int i = 0; i < len - 1; i++)
            {
                string d = seq.Substring(i, 2);
                oracle[d] = oracle.GetValueOrDefault(d) + 1;
            }
            int total = len - 1;

            var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);

            freq.Should().HaveCount(oracle.Count, "same set of distinct dinucleotides");
            foreach (var (d, cnt) in oracle)
                freq[d].Should().BeApproximately((double)cnt / total, Tolerance,
                    $"f_{d} = count/(N−1)");

            // N−1 identity: counts (= freq×total) sum to exactly N−1.
            freq.Values.Sum().Should().BeApproximately(1.0, Tolerance, "INV-01");
            oracle.Values.Sum().Should().Be(total, "N−1 overlapping pairs, no off-by-one");
            AssertWellFormed(freq);
        }
    }

    #endregion
}
