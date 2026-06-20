using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Complexity area — Lempel–Ziv (compression-based) complexity.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, OverflowException, …). Every input must result in
/// EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation exception (ArgumentException / ArgumentNullException). A raw
/// runtime exception, a hang, or a NaN/Infinity leak on garbage input is a bug,
/// not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPLEX-COMPRESS-001 — Lempel–Ziv complexity (Complexity)
/// Checklist: docs/checklists/03_FUZZING.md, row 228.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, homopolymer
///           (maximally compressible), random sequence (near-incompressible),
///           and a very long sequence (O(n²) worst-case guard under CancelAfter).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Lempel–Ziv-complexity contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-COMPLEX-COMPRESS-001 is the Lempel–Ziv (1976) complexity: parse the
/// sequence left-to-right by an exhaustive-history production process, counting
/// the number of distinct components c(S) needed to reconstruct it. Repetitive
/// input yields FEW components (low complexity); diverse/random input yields MANY
/// (high complexity). A length-normalized variant removes the length dependence:
///     LZ_norm = c / (n / log_b(n)),   b = number of distinct symbols present,
/// with LZ_norm → 1 for a maximally complex (random) sequence and smaller values
/// for more compressible input.
///   — docs/algorithms/Complexity/Lempel_Ziv_Complexity.md
///     §2.2 (core model + normalization), §2.4 (invariants INV-01..INV-05),
///     §3 (contract), §6.1 (edge cases), §7.1 (worked example).
///     Sources: Lempel & Ziv (1976) [1]; Naereen reference parse [3];
///     entropy/antropy lziv_complexity normalization [4]; Zhang et al. (2009) [5].
///
/// Every expected value below is derived INDEPENDENTLY from the doc and the
/// primary-source parse rule (Lempel_Ziv_Complexity.md §2.2 / §7.1), NOT read off
/// the code's arrays. The raw-count walk-throughs were reproduced by hand:
///   • "1001111011000010" → 1 / 0 / 01 / 11 / 10 / 110 / 00 / 010 → c = 8;
///     n=16, b=2, log₂16=4, b(n)=4, LZ_norm = 8/4 = 2.0 (§7.1).
///   • homopolymer "0"×16 → 0 / 00 / 000 / 0000 / 00000 → c = 5 (§6.1).
///   • single base "A" → c = 1 (INV-02).
/// A test that would still pass against an implementation that, say, dropped the
/// normalization or mis-counted the trailing partial component is invalid.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Surfaces and their documented contracts
/// ───────────────────────────────────────────────────────────────────────────
/// The unit exposes raw, normalized, and a `EstimateCompressionRatio` delegate,
/// each over BOTH a strict typed `DnaSequence` overload and a lenient `string`
/// overload (SequenceComplexity.cs lines 453–577). Fuzzing pins both surfaces and
/// the boundary between them so neither drifts:
///
/// (1) The TYPED overloads — Calculate*…(DnaSequence) (lines 460–491). The
///     DnaSequence argument has ALREADY passed the strict ctor validation gate, so
///     only A/C/G/T (upper-cased) ever reach the metric. Documented validation
///     (§3.3): null DnaSequence → ArgumentNullException (explicit guard, never a
///     NullReferenceException). The parse is alphabet-agnostic and total over the
///     valid input, so the result is always a finite count / finite ratio.
///
/// (2) The LENIENT string overloads — Calculate*…(string) (lines 469–522).
///     LENIENT by documented design (§3.3, §6.1): null/empty short-circuits to 0,
///     the input is upper-cased (ToUpperInvariant), and the parse is alphabet-
///     agnostic — ANY character (digits, gaps, the binary {0,1} alphabet of the
///     worked example, '\0', unicode, surrogate halves) is parsed as an opaque
///     symbol and NEVER throws. We pin that this surface stays finite and never
///     throws on pure random-byte garbage.
///
/// Documented degenerate handling pinned here (§2.2, §3.3, §6.1, deviation #2):
///   • empty / null  → 0 (raw) and 0 (normalized) — INV-01, no division.
///   • single-symbol input (alphabet b<2): the log base is undefined, so the
///     normalizer clamps b := max(b, 2); for the length-1 degenerate case
///     (log_b(1)=0) it returns the RAW count. So a homopolymer's normalized value
///     is c / (n / log₂ n) and is NOT bounded by 1 (it can exceed 1 — e.g. "AAAA"
///     gives 2 / (4/2) = 1.0, and "0"×16 gives 5 / (16/4) = 1.25). This is a
///     *defined* consequence of the b<2 clamp, NOT a bug — we pin it explicitly so
///     the homopolymer-vs-random ORDERING is asserted via the RAW count (which is
///     monotone and model-clean), reserving the exact normalized homopolymer value
///     for the hand-checked worked examples.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComplexityFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;
    private const int Seed = 20260620;

    /// <summary>Generates a random valid DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(Random rng, int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Generates a random string of arbitrary BMP code points (0x0000–0xFFFF),
    /// deliberately spanning control characters, the null byte, and lone surrogate
    /// halves — pure random-byte fuzz fodder for the lenient string surface.
    /// </summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>
    /// Independent reference parse of the raw Lempel–Ziv (1976) complexity, written
    /// straight from Lempel_Ziv_Complexity.md §2.2 (exhaustive-history rule) and the
    /// Naereen reference [3] — deliberately NOT calling the production code, so it
    /// can cross-check the implementation rather than echo it.
    /// </summary>
    private static int ReferenceLempelZiv(string seq)
    {
        var components = new HashSet<string>();
        int ind = 0, inc = 1;
        while (ind + inc <= seq.Length)
        {
            string sub = seq.Substring(ind, inc);
            if (components.Contains(sub))
            {
                inc++;
            }
            else
            {
                components.Add(sub);
                ind += inc;
                inc = 1;
            }
        }
        return components.Count;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  Positive sanity — hand-checkable worked examples (doc §7.1, §6.1)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — worked examples

    /// <summary>
    /// The §7.1 worked example, derived independently: "1001111011000010" parses as
    /// 1 / 0 / 01 / 11 / 10 / 110 / 00 / 010 → c = 8, and normalized with n=16,
    /// b=2, log₂16=4, b(n)=4 gives LZ_norm = 8/4 = 2.0. EstimateCompressionRatio is
    /// a thin delegate to the normalized value (INV-05), so it must return the same
    /// 2.0. This is the alphabet-agnostic binary example; it runs on the lenient
    /// string surface because it is not DNA.
    /// </summary>
    [Test]
    public void LempelZiv_BinaryWorkedExample_MatchesDocExactly()
    {
        const string s = "1001111011000010";

        int raw = SequenceComplexity.CalculateLempelZivComplexity(s);
        raw.Should().Be(8, "the §7.1 walk-through 1/0/01/11/10/110/00/010 produces 8 components");

        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(s);
        norm.Should().BeApproximately(2.0, Tolerance,
            "n=16, b=2, b(n)=16/log₂16=4 ⇒ LZ_norm = 8/4 = 2.0 (§7.1)");

        SequenceComplexity.EstimateCompressionRatio(s).Should().BeApproximately(2.0, Tolerance,
            "EstimateCompressionRatio delegates to the normalized value (INV-05)");
    }

    /// <summary>
    /// The homopolymer edge case from §6.1, derived independently: "0"×16 parses as
    /// 0 / 00 / 000 / 0000 / 00000 → c = 5. Pinned exactly to guard the
    /// productivity-buildup behaviour (INV-04) at the raw-count level.
    /// </summary>
    [Test]
    public void LempelZiv_Homopolymer16_RawIsFive()
    {
        int raw = SequenceComplexity.CalculateLempelZivComplexity(new string('0', 16));

        raw.Should().Be(5, "0/00/000/0000/00000 is the exhaustive-history parse of \"0\"×16 (§6.1)");
    }

    /// <summary>
    /// INV-04, the core compression contract, on DNA over the strict typed surface:
    /// a homopolymer (maximally compressible) has STRICTLY LOWER raw complexity than
    /// an all-distinct string of equal length, which in turn is ≤ a random sequence.
    /// Asserted on the RAW count (monotone, model-clean) rather than the normalized
    /// value, because the homopolymer's b<2 normalization clamp gives it a value not
    /// directly comparable to the b=4 random case (see class header).
    /// </summary>
    [Test]
    public void LempelZiv_HomopolymerStrictlyLessComplexThanDiverse_OnDna()
    {
        var rng = new Random(Seed);
        const int n = 256;

        int homopolymer = SequenceComplexity.CalculateLempelZivComplexity(
            new DnaSequence(new string('A', n)));
        int cyclic = SequenceComplexity.CalculateLempelZivComplexity(
            new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", n / 4))));
        int random = SequenceComplexity.CalculateLempelZivComplexity(
            new DnaSequence(RandomDna(rng, n)));

        homopolymer.Should().BeLessThan(cyclic,
            "a homopolymer reuses components; a diverse periodic string starts more (INV-04)");
        cyclic.Should().BeLessThanOrEqualTo(random,
            "a random sequence is near-incompressible, so it yields at least as many components");
        random.Should().BeGreaterThan(homopolymer,
            "random DNA is far more complex than a homopolymer of the same length");
    }

    /// <summary>
    /// Normalized ordering on equal-length DNA where BOTH operands have b≥2 (so the
    /// normalization is fully defined and comparable): a low-diversity, highly
    /// periodic sequence has a strictly smaller normalized complexity than random
    /// DNA, and random DNA sits near the theoretical maximum of 1 (§2.2: LZ_norm → 1
    /// for a maximally complex sequence). This pins the documented "low for
    /// repetitive, high for random" behaviour at the normalized level.
    /// </summary>
    [Test]
    public void NormalizedLempelZiv_RepetitiveLowerThanRandom_AndRandomNearOne()
    {
        var rng = new Random(Seed + 1);
        const int n = 2048;

        // Periodic dinucleotide: two distinct symbols (b=2) but maximally repetitive.
        double repetitive = SequenceComplexity.CalculateNormalizedLempelZivComplexity(
            new DnaSequence(string.Concat(Enumerable.Repeat("AC", n / 2))));
        double random = SequenceComplexity.CalculateNormalizedLempelZivComplexity(
            new DnaSequence(RandomDna(rng, n)));

        repetitive.Should().BeLessThan(random,
            "a maximally periodic sequence is more compressible ⇒ smaller normalized LZ (§2.2)");
        random.Should().BeInRange(0.5, 1.5,
            "random DNA approaches the LZ_norm → 1 asymptote (§2.2); pinned with generous slack");
        repetitive.Should().BeGreaterThan(0.0,
            "any non-empty input yields at least one component (INV-02)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  BE — Boundary: empty input
    // ═══════════════════════════════════════════════════════════════════

    #region BE — Boundary: empty / null

    /// <summary>
    /// BE: the empty string is the lower size boundary. INV-01: c(empty) = 0, and
    /// the normalized value is 0 by the explicit n==0 guard — NO division, NO NaN,
    /// NO exception (§6.1). Pinned across raw, normalized, and the delegate on the
    /// lenient string surface (which short-circuits null/empty to 0).
    /// </summary>
    [Test]
    public void EmptyAndNull_String_AreZeroAndDoNotThrow()
    {
        foreach (var input in new[] { string.Empty, (string?)null })
        {
            var act = () =>
            {
                SequenceComplexity.CalculateLempelZivComplexity(input!).Should().Be(0,
                    "c(empty) = 0 (INV-01)");
                SequenceComplexity.CalculateNormalizedLempelZivComplexity(input!).Should().Be(0.0,
                    "normalized complexity of empty/null is the defined 0, no division (§6.1)");
                SequenceComplexity.EstimateCompressionRatio(input!).Should().Be(0.0,
                    "the delegate inherits the empty/null → 0 contract");
            };

            act.Should().NotThrow("empty/null is a defined boundary input on the lenient string surface, not an error");
        }
    }

    /// <summary>
    /// BE: the empty DnaSequence (built from "" via the ctor short-circuit) is a
    /// defined input whose raw and normalized complexity are both 0 — no exception.
    /// </summary>
    [Test]
    public void EmptyDnaSequence_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            var empty = new DnaSequence(string.Empty);
            empty.Length.Should().Be(0);
            SequenceComplexity.CalculateLempelZivComplexity(empty).Should().Be(0);
            SequenceComplexity.CalculateNormalizedLempelZivComplexity(empty).Should().Be(0.0);
            SequenceComplexity.EstimateCompressionRatio(empty).Should().Be(0.0);
        };

        act.Should().NotThrow("an empty sequence is a defined boundary, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the documented ArgumentNullException boundary on the
    /// typed surface (§3.3) — an *intentional* validation exception, never a raw
    /// NullReferenceException. Pinned for all three typed entry points.
    /// </summary>
    [Test]
    public void NullDnaSequence_ThrowsArgumentNullException()
    {
        var raw = () => SequenceComplexity.CalculateLempelZivComplexity((DnaSequence)null!);
        var norm = () => SequenceComplexity.CalculateNormalizedLempelZivComplexity((DnaSequence)null!);
        var ratio = () => SequenceComplexity.EstimateCompressionRatio((DnaSequence)null!);

        raw.Should().Throw<ArgumentNullException>("the typed overload guards null explicitly (§3.3)");
        norm.Should().Throw<ArgumentNullException>("the typed overload guards null explicitly (§3.3)");
        ratio.Should().Throw<ArgumentNullException>("the typed delegate guards null explicitly (§3.3)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  BE — Boundary: single character
    // ═══════════════════════════════════════════════════════════════════

    #region BE — Boundary: single char

    /// <summary>
    /// BE: a one-symbol sequence is the minimal non-empty input. INV-02: the first
    /// symbol is always a new component, so raw c = 1 for every single character.
    /// For the single-symbol input the alphabet is b=1 (clamped to 2) but n=1 makes
    /// log_b(1)=0, so the documented degenerate guard returns the RAW count 1 as the
    /// normalized value (§2.2, §3.3, deviation #2). Pinned across DNA bases and a
    /// non-DNA symbol on the lenient surface — never a 0/0 NaN.
    /// </summary>
    [TestCase("A")]
    [TestCase("C")]
    [TestCase("G")]
    [TestCase("T")]
    [TestCase("0")]
    [TestCase("Z")]
    public void SingleCharacter_RawIsOne_AndNormalizedIsDegenerateRawCount(string single)
    {
        int raw = SequenceComplexity.CalculateLempelZivComplexity(single);
        raw.Should().Be(1, "the first symbol is always a new component (INV-02)");

        double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(single);
        norm.Should().Be(1.0,
            "for a length-1 single-symbol input log_b(1)=0, so the degenerate guard returns the raw count 1 (§2.2)");
        double.IsNaN(norm).Should().BeFalse("the degenerate guard avoids a 0/0 NaN");
    }

    /// <summary>
    /// BE: single DNA base over the strict typed surface — c = 1 (INV-02), and the
    /// normalized degenerate guard returns the raw 1. Mirrors the string case to pin
    /// the strict/lenient boundary at the single-char minimum.
    /// </summary>
    [TestCase('A')]
    [TestCase('C')]
    [TestCase('G')]
    [TestCase('T')]
    public void SingleDnaBase_RawIsOne(char baseChar)
    {
        var seq = new DnaSequence(baseChar.ToString());

        SequenceComplexity.CalculateLempelZivComplexity(seq).Should().Be(1, "INV-02");
        SequenceComplexity.CalculateNormalizedLempelZivComplexity(seq).Should().Be(1.0,
            "length-1 degenerate guard returns the raw count (§2.2)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  BE — Boundary: homopolymer (maximally compressible)
    // ═══════════════════════════════════════════════════════════════════

    #region BE — Boundary: homopolymer

    /// <summary>
    /// BE: a homopolymer is maximally compressible ⇒ minimal complexity. For "X"×n
    /// the exhaustive-history parse yields components X / XX / XXX / … whose total
    /// length is 1+2+…+m = m(m+1)/2 ≤ n, so c(S) is the largest m with m(m+1)/2 ≤ n
    /// (a triangular bound) — independently of WHICH symbol repeats. We assert this
    /// closed form against the implementation across a range of lengths and several
    /// symbols (DNA and non-DNA), and confirm c grows only ~√(2n), far below n. This
    /// pins INV-03 (c ≤ n) and INV-04 (homopolymer minimal) at the extreme.
    /// </summary>
    [Test]
    public void Homopolymer_RawComplexity_IsTriangularBound()
    {
        foreach (char sym in new[] { 'A', 'G', '0', 'Z' })
        {
            foreach (int n in new[] { 1, 2, 3, 4, 9, 10, 16, 17, 100, 1000 })
            {
                // Largest m with m(m+1)/2 <= n: this is exactly the number of
                // components 1/2/.../m the homopolymer parse emits (derived from §6.1).
                int expected = 0;
                while ((long)(expected + 1) * (expected + 2) / 2 <= n) expected++;

                string s = new string(sym, n);
                int raw = SequenceComplexity.CalculateLempelZivComplexity(s);

                raw.Should().Be(expected,
                    $"\"{sym}\"×{n} parses into the triangular number of distinct runs (§6.1)");
                raw.Should().BeLessThanOrEqualTo(n, "c(S) ≤ n (INV-03)");
            }
        }
    }

    /// <summary>
    /// BE: the §6.1 anchor "0"×16 → c=5 must equal the triangular closed form, and a
    /// same-length all-distinct-ish DNA sequence must be strictly more complex —
    /// pinning the homopolymer as the minimal-complexity extreme (INV-04) at a
    /// hand-checked length.
    /// </summary>
    [Test]
    public void Homopolymer_IsMinimalComplexity_VersusDiverseSameLength()
    {
        int homopolymer = SequenceComplexity.CalculateLempelZivComplexity(new string('0', 16));
        homopolymer.Should().Be(5, "0/00/000/0000/00000 (§6.1)");

        int diverse = SequenceComplexity.CalculateLempelZivComplexity("ACGTACGTACGTACGT");
        diverse.Should().BeGreaterThan(homopolymer, "a more diverse string of equal length is more complex (INV-04)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  BE / RB — random sequence (near-incompressible) + boundary sweep
    // ═══════════════════════════════════════════════════════════════════

    #region BE/RB — random sequence + randomized sweep

    /// <summary>
    /// BE/RB: a randomized boundary sweep. Over many fixed-seed random DNA sequences
    /// of varied length AND deliberately included degenerate boundaries (length 0,
    /// 1, 2, homopolymers), the raw and normalized complexities must ALWAYS be
    /// well-formed: raw an integer in [0, n], normalized finite (no NaN/Infinity),
    /// EstimateCompressionRatio identical to the normalized value (INV-05), the raw
    /// count equal to the independent reference parse, and the random sequence
    /// strictly more complex than a homopolymer of the same length (INV-04). Never a
    /// throw, hang, NaN, or overflow.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void RandomizedSweep_AlwaysWellFormed()
    {
        var rng = new Random(Seed + 2);

        for (int iter = 0; iter < 400; iter++)
        {
            int n = rng.Next(0, 600);                 // includes 0, 1, 2 …
            int mode = rng.Next(3);
            string s = mode switch
            {
                0 => RandomDna(rng, n),                                // random DNA (b up to 4)
                1 => new string("ACGT"[rng.Next(4)], n),               // homopolymer
                _ => string.Concat(Enumerable.Repeat("AC", n / 2)),    // periodic (b=2)
            };

            int raw = SequenceComplexity.CalculateLempelZivComplexity(s);
            raw.Should().Be(ReferenceLempelZiv(s),
                "the production raw count must match the independent reference parse");
            raw.Should().BeInRange(0, s.Length, "0 ≤ c(S) ≤ n (INV-01..INV-03)");
            if (s.Length > 0) raw.Should().BeGreaterThan(0, "non-empty input yields ≥ 1 component (INV-02)");

            double norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(s);
            double.IsNaN(norm).Should().BeFalse("normalized complexity must never be NaN");
            double.IsInfinity(norm).Should().BeFalse("normalized complexity must never be Infinity");
            norm.Should().BeGreaterThanOrEqualTo(0.0, "complexity is non-negative");

            SequenceComplexity.EstimateCompressionRatio(s).Should().Be(norm,
                "EstimateCompressionRatio is a thin delegate to the normalized value (INV-05)");

            // INV-04 at every sampled length: random/periodic ≥ homopolymer.
            if (s.Length >= 4)
            {
                int homo = SequenceComplexity.CalculateLempelZivComplexity(new string('A', s.Length));
                raw.Should().BeGreaterThanOrEqualTo(homo,
                    "no sequence is more compressible than a homopolymer of equal length (INV-04)");
            }
        }
    }

    /// <summary>
    /// RB: the lenient string surface must NEVER throw on pure random-byte garbage —
    /// arbitrary BMP code points including control chars, the null byte, and lone
    /// surrogate halves (§3.3, §6.1: alphabet-agnostic parse). The parse treats each
    /// char as an opaque symbol, so the result must be a finite non-negative value
    /// matching the independent reference parse, with no encoding surprise.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void LenientStringSurface_RandomBmpGarbage_NeverThrows()
    {
        var rng = new Random(Seed + 3);

        for (int iter = 0; iter < 300; iter++)
        {
            string s = RandomBmpChars(rng, rng.Next(0, 300));

            int raw = 0;
            double norm = 0;
            var act = () =>
            {
                raw = SequenceComplexity.CalculateLempelZivComplexity(s);
                norm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(s);
                _ = SequenceComplexity.EstimateCompressionRatio(s);
            };

            act.Should().NotThrow("the lenient string parse is alphabet-agnostic and total over char (§3.3)");

            // Compare against the reference on the UPPER-CASED input, because the
            // surface upper-cases (ToUpperInvariant) before parsing (§3.3).
            raw.Should().Be(ReferenceLempelZiv(s.ToUpperInvariant()),
                "raw count must equal the reference parse of the upper-cased input");
            double.IsNaN(norm).Should().BeFalse("never NaN");
            double.IsInfinity(norm).Should().BeFalse("never Infinity");
        }
    }

    /// <summary>
    /// BE/OVF: a very long sequence (200,000 bases) must compute without overflow or
    /// hang under a CancelAfter guard (the parse is O(n²) worst case on highly
    /// repetitive input — §4.3 — so the long homopolymer is the stress case). Both a
    /// long homopolymer and a long random sequence must yield a finite raw count in
    /// [0, n] and a finite normalized value, and random must be far more complex than
    /// the homopolymer at scale (INV-04).
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void VeryLong_NoOverflowNoHang_AndContractHolds()
    {
        var rng = new Random(Seed + 4);
        const int n = 200_000;

        string homopolymer = new string('A', n);
        string random = RandomDna(rng, n);

        int homoRaw = SequenceComplexity.CalculateLempelZivComplexity(homopolymer);
        int randomRaw = SequenceComplexity.CalculateLempelZivComplexity(random);

        homoRaw.Should().BeInRange(0, n, "c(S) ≤ n even at scale (INV-03)");
        randomRaw.Should().BeInRange(0, n, "c(S) ≤ n even at scale (INV-03)");
        randomRaw.Should().BeGreaterThan(homoRaw,
            "random DNA is near-incompressible; a homopolymer is maximally compressible (INV-04)");

        double homoNorm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(homopolymer);
        double randomNorm = SequenceComplexity.CalculateNormalizedLempelZivComplexity(random);
        double.IsNaN(homoNorm).Should().BeFalse();
        double.IsNaN(randomNorm).Should().BeFalse();
        double.IsInfinity(homoNorm).Should().BeFalse();
        double.IsInfinity(randomNorm).Should().BeFalse();
        randomNorm.Should().BeInRange(0.5, 1.5,
            "at scale random DNA sits near the LZ_norm → 1 asymptote (§2.2)");
    }

    #endregion
}
