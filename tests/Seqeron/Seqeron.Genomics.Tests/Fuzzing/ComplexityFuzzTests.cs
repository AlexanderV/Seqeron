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

    /// <summary>
    /// Independent reference DUST score, written straight from DUST_Score.md §2.2 and the
    /// Morgulis et al. (2006) / Li (2025) restatement — deliberately NOT calling the
    /// production code, so it cross-checks the implementation rather than echoes it:
    ///   S(x) = ( Σ_t c_t·(c_t − 1)/2 ) / (L − wordSize + 1),
    /// where c_t is the occurrence count of each overlapping word t. Returns 0 for the
    /// documented degenerate L &lt; wordSize case (no word exists, §3.3 / §6.1). Inputs
    /// are treated alphabet-agnostically (each char is an opaque symbol), mirroring the
    /// lenient string surface AFTER upper-casing.
    /// </summary>
    private static double ReferenceDustScore(string seq, int wordSize = 3)
    {
        if (wordSize < 1) throw new ArgumentOutOfRangeException(nameof(wordSize));
        if (string.IsNullOrEmpty(seq) || seq.Length < wordSize) return 0.0;

        int wordCount = seq.Length - wordSize + 1;
        var counts = new Dictionary<string, int>();
        for (int i = 0; i < wordCount; i++)
        {
            string w = seq.Substring(i, wordSize);
            counts[w] = counts.TryGetValue(w, out int c) ? c + 1 : 1;
        }

        double numerator = 0.0;
        foreach (int c in counts.Values)
            numerator += c * (c - 1) / 2.0;

        return numerator / wordCount;
    }

    /// <summary>
    /// Independent reference K-mer entropy, written straight from K-mer_Entropy.md
    /// §2.2 (overlapping k-mer decomposition, N = L − k + 1) and the Shannon (1948)
    /// formula H = −Σ p_i·log₂(p_i), p_i = n_i / N — deliberately NOT calling the
    /// production code, so it cross-checks the implementation rather than echoes it.
    /// Returns the documented degenerate 0 for L &lt; k (§3.3 / §6.1: no k-mer exists).
    /// Inputs are treated alphabet-agnostically (each length-k substring is an opaque
    /// symbol), mirroring the lenient string surface AFTER upper-casing.
    /// </summary>
    private static double ReferenceKmerEntropy(string seq, int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));
        if (string.IsNullOrEmpty(seq) || seq.Length < k) return 0.0;

        int n = seq.Length - k + 1;
        var counts = new Dictionary<string, int>();
        for (int i = 0; i < n; i++)
        {
            string kmer = seq.Substring(i, k);
            counts[kmer] = counts.TryGetValue(kmer, out int c) ? c + 1 : 1;
        }

        double entropy = 0.0;
        foreach (int c in counts.Values)
        {
            double p = (double)c / n;
            entropy -= p * Math.Log2(p);
        }

        return entropy;
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

    // ═══════════════════════════════════════════════════════════════════════════
    //
    //  SEQ-COMPLEX-DUST-001 — DUST low-complexity score (CalculateDustScore)
    //
    // ═══════════════════════════════════════════════════════════════════════════

    #region SEQ-COMPLEX-DUST-001 — DUST low-complexity score

    /*  ─────────────────────────────────────────────────────────────────────────
     *  Unit: SEQ-COMPLEX-DUST-001 — DUST score (SequenceComplexity.CalculateDustScore)
     *  Checklist: docs/checklists/03_FUZZING.md, row 229.
     *  Fuzz strategies exercised for THIS unit:
     *    • BE = Boundary Exploitation — empty input, input SHORTER THAN the word
     *      (triplet) length, homopolymer, non-ACGT characters.
     *    — docs/checklists/03_FUZZING.md §Description (strategy codes; BE = 0/-1/MaxInt/empty),
     *      and docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
     *
     *  ─────────────────────────────────────────────────────────────────────────
     *  The DUST contract under test (DUST_Score.md)
     *  ─────────────────────────────────────────────────────────────────────────
     *  DUST (Morgulis et al. 2006 [1]; longdust restatement, Li 2025 [2]) is a
     *  triplet-frequency heuristic for low-complexity DNA. For a length-L sequence,
     *  let c_t be the number of occurrences of each overlapping word t among the
     *  L − wordSize + 1 windows (wordSize = 3 by default). The score is (§2.2):
     *
     *      S(x) = ( Σ_t c_t·(c_t − 1)/2 ) / (L − wordSize + 1)
     *
     *  The numerator counts pairs of IDENTICAL triplets; the denominator is the
     *  number of triplets. A HIGHER score ⇒ LOWER complexity (more repeated words);
     *  fully distinct triplets give exactly 0 (§3.2, INV-02, INV-04). Documented
     *  invariants pinned here:
     *    • INV-01 S(x) ≥ 0 (each c(c−1)/2 ≥ 0, divisor > 0 for L ≥ 3).
     *    • INV-02 all-distinct triplets ⇒ S = 0.
     *    • INV-03 homopolymer of length L ⇒ S = (L−3)/2 (one triplet repeated
     *             L−2 times: (L−2)(L−3)/2 / (L−2)).
     *    • INV-04 higher S ⇒ lower complexity (repeats raise Σ c(c−1)/2).
     *  Mask threshold (§4.2): 2.0 (reference level T = 20); S > 2.0 ⇒ masked.
     *
     *  Every expected value below is derived INDEPENDENTLY from DUST_Score.md and
     *  the Morgulis/Li formula (see ReferenceDustScore + the by-hand worked
     *  examples), NOT read off the code's arrays. A test that would still pass
     *  against an implementation that, say, divided by (words − 1) instead of the
     *  word count (the §5.2 historical bug) is invalid.
     *
     *  ─────────────────────────────────────────────────────────────────────────
     *  Surfaces and their documented validation (§3.1, §3.3, §6.1)
     *  ─────────────────────────────────────────────────────────────────────────
     *  (1) TYPED  CalculateDustScore(DnaSequence, wordSize=3): the DnaSequence has
     *      already passed strict ctor validation, so only A/C/G/T reach the core.
     *      null DnaSequence ⇒ ArgumentNullException (explicit guard, never an NRE).
     *  (2) LENIENT CalculateDustScore(string, wordSize=3): null/empty ⇒ 0; the input
     *      is upper-cased (ToUpperInvariant) and the word tally is alphabet-agnostic,
     *      so ANY character (digits, gaps, '\0', unicode, surrogate halves) is parsed
     *      as an opaque symbol and NEVER throws.
     *  Both surfaces: wordSize < 1 ⇒ ArgumentOutOfRangeException; L < wordSize ⇒ 0
     *  (no word exists — defined-output convention, §3.3 / deviation #1).
     *  ───────────────────────────────────────────────────────────────────────── */

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity — hand-checkable worked examples (DUST §7.1, INV-03)
    // ───────────────────────────────────────────────────────────────────

    #region DUST — Positive sanity — worked examples

    /// <summary>
    /// The §7.1 worked example, derived independently: "AAAAAA" (L=6) has triplet
    /// AAA at positions 0..3 ⇒ c_AAA = 4, numerator = 4·3/2 = 6, divisor = L−2 = 4,
    /// so S = 6/4 = 1.5. This also equals INV-03's (L−3)/2 = 3/2. Pinned exactly to
    /// guard the triplet-count formula AND the §5.2-corrected divisor (number of
    /// words, NOT words−1, which would give 6/3 = 2.0).
    /// </summary>
    [Test]
    public void Dust_HexamerA_WorkedExample_IsOnePointFive()
    {
        double score = SequenceComplexity.CalculateDustScore("AAAAAA");

        score.Should().BeApproximately(1.5, Tolerance,
            "AAA repeats 4×: numerator 4·3/2=6, divisor L−2=4 ⇒ 6/4 = 1.5 (§7.1, INV-03)");
        score.Should().BeApproximately(ReferenceDustScore("AAAAAA"), Tolerance,
            "must match the independent triplet-frequency reference");
    }

    /// <summary>
    /// The second §7.1 worked example, derived independently: "ACGTACGT" (L=8) has
    /// triplets ACG=2, CGT=2, GTA=1, TAC=1 ⇒ numerator = 1 + 1 + 0 + 0 = 2, divisor
    /// = L−2 = 6, so S = 2/6 = 1/3 ≈ 0.3333. Run on BOTH the lenient string surface
    /// and the strict typed DnaSequence surface (the input is valid DNA), pinning the
    /// strict/lenient agreement on a non-degenerate case.
    /// </summary>
    [Test]
    public void Dust_AcgtAcgt_WorkedExample_IsOneThird()
    {
        const string s = "ACGTACGT";

        double viaString = SequenceComplexity.CalculateDustScore(s);
        viaString.Should().BeApproximately(1.0 / 3.0, Tolerance,
            "ACG=2,CGT=2,GTA=1,TAC=1 ⇒ numerator 1+1=2, divisor 6 ⇒ 1/3 (§7.1)");

        double viaDna = SequenceComplexity.CalculateDustScore(new DnaSequence(s));
        viaDna.Should().BeApproximately(viaString, Tolerance,
            "the typed and lenient surfaces agree on valid DNA");
        viaString.Should().BeApproximately(ReferenceDustScore(s), Tolerance,
            "must match the independent triplet-frequency reference");
    }

    /// <summary>
    /// INV-02, the zero-floor contract: a sequence whose every triplet is distinct
    /// scores exactly 0. "ACGTACG" reuses words, so instead we use a hand-built
    /// all-distinct-triplet string. The shortest non-trivial all-distinct case is two
    /// triplets that differ: "ACGA" (L=4) has triplets ACG, CGA — both unique — so
    /// numerator = 0, divisor = 2, S = 0. Pinned across a few all-distinct strings.
    /// </summary>
    [Test]
    public void Dust_AllDistinctTriplets_IsZero()
    {
        foreach (string s in new[] { "ACGA", "ACGT", "ACGTG", "ACGTGCA" })
        {
            double score = SequenceComplexity.CalculateDustScore(s);
            score.Should().BeApproximately(0.0, Tolerance,
                $"\"{s}\" has all-distinct triplets ⇒ every c_t=1 ⇒ S=0 (INV-02)");
            score.Should().BeApproximately(ReferenceDustScore(s), Tolerance,
                "must match the independent reference");
        }
    }

    /// <summary>
    /// INV-03 / INV-04, the core "higher score ⇒ lower complexity" contract: a
    /// homopolymer (one repeated triplet) attains the documented (L−3)/2 score, which
    /// is STRICTLY higher than the score of a diverse same-length sequence (which
    /// approaches 0). This is the discriminating fuzz assertion — a degenerate
    /// implementation returning a constant or echoing the count would fail it.
    /// </summary>
    [Test]
    public void Dust_HomopolymerScoresHigherThanDiverse_SameLength()
    {
        const int n = 64;

        double homopolymer = SequenceComplexity.CalculateDustScore(new DnaSequence(new string('A', n)));
        double diverse = SequenceComplexity.CalculateDustScore(
            new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", n / 4))));

        homopolymer.Should().BeApproximately((n - 3) / 2.0, Tolerance,
            "homopolymer of length L scores (L−3)/2 (INV-03)");
        homopolymer.Should().BeGreaterThan(diverse,
            "a homopolymer is far lower complexity ⇒ far higher DUST score (INV-04)");
        diverse.Should().BeGreaterThanOrEqualTo(0.0, "S ≥ 0 (INV-01)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: empty / null
    // ───────────────────────────────────────────────────────────────────

    #region DUST — BE — Boundary: empty / null

    /// <summary>
    /// BE: empty/null on the lenient string surface short-circuits to a defined 0 —
    /// NO division by L−2 (which would be −2, then 0/0 on the count), NO NaN, NO
    /// exception (§3.3, §6.1: "no words ⇒ minimal complexity").
    /// </summary>
    [Test]
    public void Dust_EmptyAndNullString_AreZeroAndDoNotThrow()
    {
        foreach (var input in new[] { string.Empty, (string?)null })
        {
            var act = () => SequenceComplexity.CalculateDustScore(input!);
            act.Should().NotThrow("empty/null is a defined boundary on the lenient surface (§6.1)");
            SequenceComplexity.CalculateDustScore(input!).Should().Be(0.0,
                "no words ⇒ S = 0, no division (§6.1)");
        }
    }

    /// <summary>
    /// BE: the empty DnaSequence (built from "" via the ctor short-circuit) scores a
    /// defined 0 — no division-by-zero on the L−2 = −2 word count, no NaN.
    /// </summary>
    [Test]
    public void Dust_EmptyDnaSequence_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            var empty = new DnaSequence(string.Empty);
            empty.Length.Should().Be(0);
            double score = SequenceComplexity.CalculateDustScore(empty);
            score.Should().Be(0.0, "no triplet exists ⇒ S = 0 (§6.1)");
            double.IsNaN(score).Should().BeFalse("the L < wordSize guard avoids a 0/0 NaN");
        };

        act.Should().NotThrow("an empty sequence is a defined boundary, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the documented ArgumentNullException boundary on the
    /// typed surface (§3.3, §6.1) — an INTENTIONAL validation exception, never a raw
    /// NullReferenceException.
    /// </summary>
    [Test]
    public void Dust_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => SequenceComplexity.CalculateDustScore((DnaSequence)null!);
        act.Should().Throw<ArgumentNullException>("the typed overload guards null explicitly (§3.3)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: sequence shorter than the word (triplet) length
    // ───────────────────────────────────────────────────────────────────

    #region DUST — BE — Boundary: shorter than the word

    /// <summary>
    /// BE: a sequence SHORTER THAN one word forms no triplet, so the documented
    /// convention (§3.3, §6.1, deviation #1) returns a defined 0 — never a divide by
    /// the negative/zero word count L − wordSize + 1, never a NaN, never an
    /// IndexOutOfRange from Substring(i, wordSize). Pinned for L = 0,1,2 with the
    /// default triplet word, on both the lenient string and strict typed surfaces.
    /// </summary>
    [TestCase("")]
    [TestCase("A")]
    [TestCase("AC")]
    [TestCase("AG")]
    public void Dust_ShorterThanWord_IsZero(string s)
    {
        double viaString = SequenceComplexity.CalculateDustScore(s);   // wordSize = 3
        viaString.Should().Be(0.0, "L < wordSize ⇒ no triplet ⇒ S = 0 (§3.3, §6.1)");
        double.IsNaN(viaString).Should().BeFalse("no 0/0 on the negative/zero word count");
        viaString.Should().Be(ReferenceDustScore(s), "matches the independent reference");

        if (s.Length > 0)
        {
            double viaDna = SequenceComplexity.CalculateDustScore(new DnaSequence(s));
            viaDna.Should().Be(0.0, "strict surface obeys the same L < wordSize ⇒ 0 convention");
        }
    }

    /// <summary>
    /// BE: the "shorter than the WORD" boundary is parameterised by wordSize. For any
    /// wordSize > L the score is the defined 0; the exact boundary L == wordSize gives
    /// exactly ONE word (S = 0, a single distinct word ⇒ c=1 ⇒ c(c−1)/2 = 0, divided
    /// by 1). This pins the off-by-one around the window edge — the first length at
    /// which a word exists yields a finite, non-NaN 0, not a crash.
    /// </summary>
    [TestCase("ACGTAC", 4)]   // L=6 > wordSize=4 ⇒ words exist
    [TestCase("ACG", 3)]      // L == wordSize ⇒ exactly one word, S = 0
    [TestCase("AC", 3)]       // L < wordSize ⇒ 0
    [TestCase("A", 5)]        // L << wordSize ⇒ 0
    [TestCase("ACGTACGT", 8)] // whole sequence is one word, S = 0
    public void Dust_WordSizeBoundary_IsDefinedAndFinite(string s, int wordSize)
    {
        double score = SequenceComplexity.CalculateDustScore(s, wordSize);

        double.IsNaN(score).Should().BeFalse("no NaN at the window edge");
        double.IsInfinity(score).Should().BeFalse("no Infinity at the window edge");
        score.Should().BeGreaterThanOrEqualTo(0.0, "S ≥ 0 (INV-01)");
        score.Should().BeApproximately(ReferenceDustScore(s, wordSize), Tolerance,
            "must match the independent generalized reference S = Σ c(c−1)/2 / (L − wordSize + 1)");
    }

    /// <summary>
    /// BE: wordSize &lt; 1 is the documented ArgumentOutOfRangeException boundary
    /// (§3.3) on both surfaces — including the BE archetype 0 and −1. An intentional
    /// validation throw, never a DivideByZero or empty-loop silent 0.
    /// </summary>
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void Dust_WordSizeBelowOne_ThrowsArgumentOutOfRange(int wordSize)
    {
        var viaString = () => SequenceComplexity.CalculateDustScore("ACGTACGT", wordSize);
        var viaDna = () => SequenceComplexity.CalculateDustScore(new DnaSequence("ACGTACGT"), wordSize);

        viaString.Should().Throw<ArgumentOutOfRangeException>("wordSize < 1 is invalid (§3.3)");
        viaDna.Should().Throw<ArgumentOutOfRangeException>("wordSize < 1 is invalid (§3.3)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: homopolymer (maximal DUST score / lowest complexity)
    // ───────────────────────────────────────────────────────────────────

    #region DUST — BE — Boundary: homopolymer

    /// <summary>
    /// BE: a homopolymer is the lowest-complexity extreme ⇒ MAXIMAL DUST score. For
    /// "X"×L (L ≥ 3) there is exactly ONE triplet repeated L−2 times, so numerator =
    /// (L−2)(L−3)/2 and S = (L−3)/2 (INV-03), independently of WHICH base/symbol
    /// repeats. Verified against the closed form across many lengths and several
    /// symbols (DNA bases and a non-DNA symbol on the lenient surface). This pins both
    /// the triplet-count formula and the symbol-agnostic core.
    /// </summary>
    [Test]
    public void Dust_Homopolymer_ScoreIsLMinus3Over2()
    {
        foreach (char sym in new[] { 'A', 'C', 'G', 'T', '0', 'Z' })
        {
            foreach (int n in new[] { 3, 4, 5, 6, 10, 16, 100, 1000 })
            {
                string s = new string(sym, n);
                double score = SequenceComplexity.CalculateDustScore(s);

                score.Should().BeApproximately((n - 3) / 2.0, Tolerance,
                    $"\"{sym}\"×{n}: one triplet repeated {n - 2}× ⇒ S = (L−3)/2 (INV-03)");
                score.Should().BeApproximately(ReferenceDustScore(s), Tolerance,
                    "matches the independent reference");
            }
        }
    }

    /// <summary>
    /// BE: at the homopolymer extreme the DUST score grows without bound with length
    /// ((L−3)/2 → ∞), so a long homopolymer must EXCEED the masking threshold 2.0
    /// (§4.2) while a short one (L ≤ 7 ⇒ (L−3)/2 ≤ 2) does not. This pins the
    /// documented masking semantics at the boundary length L = 7 (score exactly 2.0,
    /// NOT strictly above threshold) and L = 8 (score 2.5, above).
    /// </summary>
    [Test]
    public void Dust_Homopolymer_MaskThresholdBoundary()
    {
        // L = 7 ⇒ (7−3)/2 = 2.0, exactly the threshold (not strictly above).
        SequenceComplexity.CalculateDustScore(new string('A', 7))
            .Should().BeApproximately(2.0, Tolerance, "(L−3)/2 = 2.0 at L=7 (§4.2 threshold)");
        // L = 8 ⇒ 2.5, strictly above ⇒ would be masked.
        double l8 = SequenceComplexity.CalculateDustScore(new string('A', 8));
        l8.Should().BeApproximately(2.5, Tolerance, "(L−3)/2 = 2.5 at L=8");
        l8.Should().BeGreaterThan(2.0, "a longer homopolymer exceeds the mask threshold (§4.2)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE / INJ — non-ACGT characters
    // ───────────────────────────────────────────────────────────────────

    #region DUST — BE/INJ — non-ACGT characters

    /// <summary>
    /// BE/INJ: non-ACGT characters. The STRICT typed surface rejects them at the
    /// DnaSequence ctor with ArgumentException (so the metric never sees garbage),
    /// while the LENIENT string surface is alphabet-agnostic: it upper-cases and
    /// tallies each char as an opaque triplet symbol, never throwing (§3.3). We pin
    /// both: a non-DNA homopolymer "NNNNNN" scores the same (L−3)/2 = 1.5 as "AAAAAA"
    /// (the symbol identity is irrelevant), and gaps/digits form valid distinct
    /// triplets that score 0 when all-distinct.
    /// </summary>
    [Test]
    public void Dust_NonAcgt_LenientStringScoresStructurally_StrictRejects()
    {
        // Lenient surface: alphabet-agnostic, non-DNA homopolymer behaves like a DNA one.
        SequenceComplexity.CalculateDustScore("NNNNNN")
            .Should().BeApproximately(1.5, Tolerance,
                "the core is symbol-agnostic: \"N\"×6 scores (L−3)/2 = 1.5 like \"A\"×6 (§3.3)");
        SequenceComplexity.CalculateDustScore("------")
            .Should().BeApproximately(1.5, Tolerance, "gap homopolymer scores identically (§3.3)");

        // Lenient: a non-DNA all-distinct-triplet string still scores 0 (INV-02).
        double digits = SequenceComplexity.CalculateDustScore("12345");
        digits.Should().BeApproximately(ReferenceDustScore("12345"), Tolerance,
            "matches the independent reference");
        digits.Should().Be(0.0, "digits 123/234/345 are all distinct ⇒ S = 0 (INV-02)");

        // Strict surface: non-ACGT is rejected at the ctor, never reaching the metric.
        var build = () => new DnaSequence("ACGTNACGT");
        build.Should().Throw<ArgumentException>("the DnaSequence ctor rejects non-ACGT (strict gate)");
    }

    /// <summary>
    /// BE/INJ: the lenient string surface must NEVER throw on pure random-byte garbage
    /// — arbitrary BMP code points including control chars, the null byte, and lone
    /// surrogate halves. Each char is parsed as an opaque triplet symbol after upper-
    /// casing (§3.3), so the result must be a finite non-negative value matching the
    /// independent reference parse of the UPPER-CASED input, with no encoding surprise.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Dust_LenientStringSurface_RandomBmpGarbage_NeverThrows()
    {
        var rng = new Random(Seed + 5);

        for (int iter = 0; iter < 300; iter++)
        {
            string s = RandomBmpChars(rng, rng.Next(0, 300));

            double score = 0;
            var act = () => { score = SequenceComplexity.CalculateDustScore(s); };
            act.Should().NotThrow("the lenient string parse is alphabet-agnostic and total over char (§3.3)");

            double.IsNaN(score).Should().BeFalse("never NaN on garbage");
            double.IsInfinity(score).Should().BeFalse("never Infinity on garbage");
            score.Should().BeGreaterThanOrEqualTo(0.0, "S ≥ 0 (INV-01)");
            score.Should().BeApproximately(ReferenceDustScore(s.ToUpperInvariant()), Tolerance,
                "must equal the reference score of the upper-cased input (§3.3)");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE / RB — randomized boundary sweep + scale guard
    // ───────────────────────────────────────────────────────────────────

    #region DUST — BE/RB — randomized sweep + scale guard

    /// <summary>
    /// BE/RB: a randomized boundary sweep. Over many fixed-seed inputs of varied
    /// length AND deliberately included degenerate boundaries (length 0,1,2 shorter-
    /// than-word, homopolymers, periodic strings), the DUST score must ALWAYS be
    /// well-formed: finite (no NaN/Infinity), ≥ 0 (INV-01), exactly equal to the
    /// independent triplet-frequency reference, never bounded above the homopolymer
    /// score of the same length (the maximal-score extreme, INV-04), and 0 whenever
    /// the input is shorter than a triplet (§6.1). Never a throw, hang, or overflow.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Dust_RandomizedSweep_AlwaysWellFormed()
    {
        var rng = new Random(Seed + 6);

        for (int iter = 0; iter < 400; iter++)
        {
            int n = rng.Next(0, 600);                 // includes 0, 1, 2 …
            int mode = rng.Next(3);
            string s = mode switch
            {
                0 => RandomDna(rng, n),                                 // random DNA
                1 => new string("ACGT"[rng.Next(4)], n),                // homopolymer (max score)
                _ => string.Concat(Enumerable.Repeat("AC", n / 2)),     // periodic dinucleotide
            };

            double score = SequenceComplexity.CalculateDustScore(s);

            double.IsNaN(score).Should().BeFalse("DUST score must never be NaN");
            double.IsInfinity(score).Should().BeFalse("DUST score must never be Infinity");
            score.Should().BeGreaterThanOrEqualTo(0.0, "S ≥ 0 (INV-01)");
            score.Should().BeApproximately(ReferenceDustScore(s), Tolerance,
                "the production score must match the independent triplet-frequency reference");

            if (s.Length < 3)
                score.Should().Be(0.0, "shorter than a triplet ⇒ S = 0 (§6.1)");

            // INV-04: no sequence scores higher than a homopolymer of equal length.
            if (s.Length >= 3)
            {
                double homo = SequenceComplexity.CalculateDustScore(new string('A', s.Length));
                score.Should().BeLessThanOrEqualTo(homo + Tolerance,
                    "a homopolymer is the lowest-complexity ⇒ highest-score input at each length (INV-04)");
            }
        }
    }

    /// <summary>
    /// BE/OVF: a very long sequence (200,000 bases) must score without overflow or
    /// hang under a CancelAfter guard. Both a long homopolymer (the maximal-score
    /// extreme, S = (L−3)/2 ≈ 1e5) and a long random sequence must yield a finite,
    /// non-negative score matching the closed-form/reference, with the homopolymer far
    /// exceeding the random score (INV-04) and the masking threshold (§4.2).
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Dust_VeryLong_NoOverflowNoHang_AndContractHolds()
    {
        var rng = new Random(Seed + 7);
        const int n = 200_000;

        string homopolymer = new string('A', n);
        string random = RandomDna(rng, n);

        double homoScore = SequenceComplexity.CalculateDustScore(homopolymer);
        double randomScore = SequenceComplexity.CalculateDustScore(random);

        double.IsNaN(homoScore).Should().BeFalse();
        double.IsNaN(randomScore).Should().BeFalse();
        double.IsInfinity(homoScore).Should().BeFalse();
        double.IsInfinity(randomScore).Should().BeFalse();

        homoScore.Should().BeApproximately((n - 3) / 2.0, Tolerance,
            "homopolymer score (L−3)/2 holds at scale (INV-03)");
        randomScore.Should().BeGreaterThanOrEqualTo(0.0, "S ≥ 0 (INV-01)");
        homoScore.Should().BeGreaterThan(randomScore,
            "the homopolymer is the lowest-complexity ⇒ highest-score input (INV-04)");
        homoScore.Should().BeGreaterThan(2.0, "a long homopolymer is far above the mask threshold (§4.2)");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════════════
    //
    //  SEQ-COMPLEX-KMER-001 — K-mer entropy (CalculateKmerEntropy)
    //
    // ═══════════════════════════════════════════════════════════════════════════

    #region SEQ-COMPLEX-KMER-001 — K-mer entropy

    /*  ─────────────────────────────────────────────────────────────────────────
     *  Unit: SEQ-COMPLEX-KMER-001 — K-mer entropy (SequenceComplexity.CalculateKmerEntropy)
     *  Checklist: docs/checklists/03_FUZZING.md, row 230.
     *  Fuzz strategies exercised for THIS unit:
     *    • BE = Boundary Exploitation — empty input, len < k (no k-mer fits), the
     *      degenerate k=0 parameter (the BE archetype 0, also −1 / MinInt), a
     *      homopolymer (single distinct k-mer ⇒ minimal entropy 0), and a very long
     *      sequence (O(N·k) scale guard under CancelAfter).
     *    — docs/checklists/03_FUZZING.md §Description (strategy codes; BE = 0/-1/MaxInt/empty),
     *      and docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
     *
     *  ─────────────────────────────────────────────────────────────────────────
     *  The K-mer-entropy contract under test (K-mer_Entropy.md)
     *  ─────────────────────────────────────────────────────────────────────────
     *  K-mer entropy (longdust, Li 2025 [1]; Shannon 1948 [3]) measures complexity
     *  as the Shannon entropy (in BITS) of the frequency distribution of a sequence's
     *  overlapping k-mers. Decompose a length-L sequence into N = L − k + 1 overlapping
     *  k-mers (sliding window, step 1); with n_i the count of distinct k-mer i and
     *  p_i = n_i / N, the score is (§2.2):
     *
     *      H = −Σ_i p_i · log₂(p_i)      (bits)
     *
     *  A few dominant k-mers (repeats, homopolymers) ⇒ LOW H; a near-uniform k-mer
     *  distribution (random) ⇒ HIGH H. Documented invariants pinned here:
     *    • INV-01 0 ≤ H ≤ log₂(N), N = L − k + 1 (Shannon bound [3]).
     *    • INV-02 a single distinct k-mer (deterministic, p=1) ⇒ H = 0 — homopolymer.
     *    • INV-03 all k-mers distinct (uniform) ⇒ H = log₂(N) — the maximum.
     *    • INV-04 result is case-invariant (input upper-cased before counting).
     *
     *  NOTE ON UNIT IDENTITY: the doc carrying Test Unit ID SEQ-COMPLEX-KMER-001 is
     *  docs/algorithms/Complexity/K-mer_Entropy.md, whose implementation is
     *  CalculateKmerEntropy (Shannon entropy of the k-mer distribution, range
     *  [0, log₂ N] BITS). It is NOT the observed/possible-subword "linguistic
     *  complexity" ratio (CalculateLinguisticComplexity = SEQ-COMPLEX-001, row 5).
     *  The checklist row-230 targets ("len < k", "k=0") match the k-mer-entropy
     *  contract exactly; we therefore test the entropy contract per the cited doc.
     *
     *  Every expected value below is derived INDEPENDENTLY from K-mer_Entropy.md and
     *  the Shannon formula (see ReferenceKmerEntropy + the by-hand §7.1 walk-through),
     *  NOT read off the code's arrays. A test that would still pass against an
     *  implementation that, say, used a natural-log base, or divided by L instead of
     *  N = L − k + 1, or omitted the L < k guard, is invalid.
     *
     *  ─────────────────────────────────────────────────────────────────────────
     *  Surfaces and their documented validation (§3.1, §3.3, §6.1)
     *  ─────────────────────────────────────────────────────────────────────────
     *  (1) TYPED  CalculateKmerEntropy(DnaSequence, k=2): the DnaSequence has already
     *      passed strict ctor validation, so only A/C/G/T reach the core.
     *      null DnaSequence ⇒ ArgumentNullException (explicit guard, never an NRE).
     *  (2) LENIENT CalculateKmerEntropy(string, k=2): null/empty ⇒ 0; the input is
     *      upper-cased (ToUpperInvariant) and the k-mer tally is alphabet-agnostic,
     *      so ANY character (digits, gaps, '\0', unicode, surrogate halves) is parsed
     *      as an opaque k-mer symbol and NEVER throws.
     *  Both surfaces: k < 1 ⇒ ArgumentOutOfRangeException (the degenerate-k guard,
     *  including the BE archetype 0); L < k ⇒ 0 (no k-mer exists, §3.3 / §6.1).
     *  ───────────────────────────────────────────────────────────────────────── */

    // ───────────────────────────────────────────────────────────────────
    //  Positive sanity — hand-checkable worked examples (KMER §7.1, INV-02/03)
    // ───────────────────────────────────────────────────────────────────

    #region KMER — Positive sanity — worked examples

    /// <summary>
    /// The §7.1 worked example, derived independently: "ATATAT", k=2 → dimers
    /// AT,TA,AT,TA,AT ⇒ AT=3, TA=2, N=5; p = 0.6, 0.4; H = −0.6·log₂0.6 − 0.4·log₂0.4
    /// = 0.4421793565 + 0.5287712380 = 0.9709505944546686 bits. The distinct-k-mer
    /// count (2) and the denominator N = L−k+1 = 5 are each verified independently of
    /// the code. Run on BOTH the strict typed DnaSequence surface and the lenient
    /// string surface (the input is valid DNA), pinning their agreement.
    /// </summary>
    [Test]
    public void Kmer_AtatatK2_WorkedExample_MatchesDocExactly()
    {
        const string s = "ATATAT";
        const double expected = 0.9709505944546686;

        double viaDna = SequenceComplexity.CalculateKmerEntropy(new DnaSequence(s), 2);
        viaDna.Should().BeApproximately(expected, Tolerance,
            "AT=3,TA=2 among N=5 dimers ⇒ H = −0.6log₂0.6 − 0.4log₂0.4 = 0.97095059… (§7.1)");

        double viaString = SequenceComplexity.CalculateKmerEntropy(s, 2);
        viaString.Should().BeApproximately(viaDna, Tolerance,
            "the typed and lenient surfaces agree on valid DNA");
        viaString.Should().BeApproximately(ReferenceKmerEntropy(s, 2), Tolerance,
            "must match the independent Shannon k-mer-entropy reference");
    }

    /// <summary>
    /// INV-03, the maximum: when every k-mer is distinct the distribution is uniform
    /// and H = log₂(N). "ACGT", k=2 → AC,CG,GT (all distinct), N = 3, so H = log₂3 =
    /// 1.5849625007211562 bits; and "ACGT", k=1 → A,C,G,T all distinct, N = 4, so
    /// H = log₂4 = 2.0. The distinct count equals N in each case, independently checked.
    /// </summary>
    [Test]
    public void Kmer_AllDistinct_IsLog2N()
    {
        // ACGT, k=2: AC,CG,GT all distinct ⇒ N=3 ⇒ H = log₂3.
        SequenceComplexity.CalculateKmerEntropy(new DnaSequence("ACGT"), 2)
            .Should().BeApproximately(Math.Log2(3), Tolerance,
                "3 distinct dimers over N=3 ⇒ uniform ⇒ H = log₂3 (INV-03)");

        // ACGT, k=1: A,C,G,T all distinct over N=4 ⇒ H = log₂4 = 2.0.
        SequenceComplexity.CalculateKmerEntropy(new DnaSequence("ACGT"), 1)
            .Should().BeApproximately(2.0, Tolerance,
                "A,C,G,T all distinct over N=4 ⇒ H = log₂4 = 2.0 (INV-03)");
    }

    /// <summary>
    /// INV-02, the minimum: a homopolymer produces ONE distinct k-mer repeated N
    /// times (deterministic distribution, p=1), so H = 0 exactly — the lowest
    /// complexity. Pinned across several lengths/symbols (DNA + a non-DNA symbol on
    /// the lenient surface) so the value is symbol-agnostic, never a 0/0 NaN.
    /// </summary>
    [Test]
    public void Kmer_Homopolymer_EntropyIsZero()
    {
        foreach (char sym in new[] { 'A', 'C', 'G', 'T', '0', 'Z' })
        {
            foreach (int n in new[] { 2, 3, 5, 10, 64, 1000 })
            {
                string s = new string(sym, n);
                double h = SequenceComplexity.CalculateKmerEntropy(s, 2);

                h.Should().BeApproximately(0.0, Tolerance,
                    $"\"{sym}\"×{n}, k=2: a single distinct dimer (p=1) ⇒ H = 0 (INV-02)");
                double.IsNaN(h).Should().BeFalse("a single deterministic k-mer never yields a 0/0 NaN");
                h.Should().BeApproximately(ReferenceKmerEntropy(s, 2), Tolerance,
                    "matches the independent reference");
            }
        }
    }

    /// <summary>
    /// INV-02 vs the general case, the discriminating contract: a homopolymer (one
    /// distinct k-mer, H=0) has STRICTLY LOWER k-mer entropy than a diverse
    /// same-length sequence, which in turn must satisfy the Shannon upper bound
    /// H ≤ log₂(N). A degenerate implementation returning a constant or echoing the
    /// count would fail this.
    /// </summary>
    [Test]
    public void Kmer_HomopolymerStrictlyLessEntropyThanDiverse_SameLength()
    {
        const int n = 256;

        double homopolymer = SequenceComplexity.CalculateKmerEntropy(
            new DnaSequence(new string('A', n)), 3);
        double diverse = SequenceComplexity.CalculateKmerEntropy(
            new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", n / 4))), 3);

        homopolymer.Should().BeApproximately(0.0, Tolerance, "single distinct 3-mer ⇒ H = 0 (INV-02)");
        diverse.Should().BeGreaterThan(homopolymer,
            "a diverse sequence spreads probability over many k-mers ⇒ higher entropy");
        diverse.Should().BeLessThanOrEqualTo(Math.Log2(n - 3 + 1) + Tolerance,
            "H ≤ log₂(N), N = L − k + 1 (INV-01)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: empty / null
    // ───────────────────────────────────────────────────────────────────

    #region KMER — BE — Boundary: empty / null

    /// <summary>
    /// BE: empty/null on the lenient string surface short-circuits to a defined 0 —
    /// NO division by N, NO NaN, NO exception (§3.3, §6.1: null/empty string ⇒ 0).
    /// </summary>
    [Test]
    public void Kmer_EmptyAndNullString_AreZeroAndDoNotThrow()
    {
        foreach (var input in new[] { string.Empty, (string?)null })
        {
            var act = () => SequenceComplexity.CalculateKmerEntropy(input!, 2);
            act.Should().NotThrow("empty/null is a defined boundary on the lenient surface (§6.1)");
            SequenceComplexity.CalculateKmerEntropy(input!, 2).Should().Be(0.0,
                "null/empty string ⇒ H = 0, no division (§3.3, §6.1)");
        }
    }

    /// <summary>
    /// BE: the empty DnaSequence (built from "" via the ctor short-circuit) yields a
    /// defined 0 (L = 0 &lt; k ⇒ no k-mer) — no division-by-zero on N, no NaN.
    /// </summary>
    [Test]
    public void Kmer_EmptyDnaSequence_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            var empty = new DnaSequence(string.Empty);
            empty.Length.Should().Be(0);
            double h = SequenceComplexity.CalculateKmerEntropy(empty, 2);
            h.Should().Be(0.0, "L = 0 < k ⇒ no k-mer ⇒ H = 0 (§6.1)");
            double.IsNaN(h).Should().BeFalse("the L < k guard avoids a 0/0 NaN");
        };

        act.Should().NotThrow("an empty sequence is a defined boundary, not an error");
    }

    /// <summary>
    /// BE: a null DnaSequence is the documented ArgumentNullException boundary on the
    /// typed surface (§3.3, §6.1) — an INTENTIONAL validation exception, never a raw
    /// NullReferenceException.
    /// </summary>
    [Test]
    public void Kmer_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => SequenceComplexity.CalculateKmerEntropy((DnaSequence)null!, 2);
        act.Should().Throw<ArgumentNullException>("the typed overload guards null explicitly (§3.3)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: len < k  (no k-mer fits)
    // ───────────────────────────────────────────────────────────────────

    #region KMER — BE — Boundary: len < k

    /// <summary>
    /// BE: a sequence SHORTER THAN k forms no k-mer, so the documented convention
    /// (§3.3, §6.1: "k > L ⇒ 0") returns a defined 0 — never a divide by the
    /// negative/zero count N = L − k + 1, never a NaN, never an IndexOutOfRange from
    /// Substring(i, k). Pinned across L = 0,1,2,3,4 with k = 5, and the exact boundary
    /// L == k (one k-mer ⇒ deterministic ⇒ H = 0), on both surfaces.
    /// </summary>
    [TestCase("", 5)]
    [TestCase("A", 5)]
    [TestCase("AC", 5)]
    [TestCase("ACG", 5)]
    [TestCase("ACGT", 5)]    // L=4 < k=5 ⇒ 0
    [TestCase("ACGTA", 5)]   // L == k ⇒ exactly one k-mer, deterministic ⇒ H = 0
    public void Kmer_LenLessThanK_IsZero(string s, int k)
    {
        double viaString = SequenceComplexity.CalculateKmerEntropy(s, k);
        viaString.Should().Be(0.0,
            "L < k ⇒ no k-mer ⇒ H = 0; L == k ⇒ one k-mer (p=1) ⇒ H = 0 (§3.3, §6.1)");
        double.IsNaN(viaString).Should().BeFalse("no 0/0 on the negative/zero k-mer count N");
        viaString.Should().Be(ReferenceKmerEntropy(s, k), "matches the independent reference");

        if (s.Length > 0)
        {
            double viaDna = SequenceComplexity.CalculateKmerEntropy(new DnaSequence(s), k);
            viaDna.Should().Be(0.0, "strict surface obeys the same k > L ⇒ 0 convention");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE — Boundary: degenerate k  (k = 0, −1, MinInt)
    // ───────────────────────────────────────────────────────────────────

    #region KMER — BE — Boundary: degenerate k

    /// <summary>
    /// BE: k &lt; 1 is the documented ArgumentOutOfRangeException boundary (§3.3) on
    /// both surfaces — including the BE archetype 0 and −1 and the MinInt extreme. An
    /// INTENTIONAL validation throw guarding the would-be divide-by-zero / empty-window
    /// degenerate, never a DivideByZero, never an IndexOutOfRange, never a silent 0.
    /// The k &lt; 1 guard fires even on empty/valid input (it is checked before the
    /// null/empty short-circuit on the string surface — source-verified), so a garbage
    /// k can never slip a NaN through on any input.
    /// </summary>
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void Kmer_KBelowOne_ThrowsArgumentOutOfRange(int k)
    {
        var viaString = () => SequenceComplexity.CalculateKmerEntropy("ACGTACGT", k);
        var viaDna = () => SequenceComplexity.CalculateKmerEntropy(new DnaSequence("ACGTACGT"), k);
        // Guard precedes the empty/null short-circuit on the string surface.
        var viaEmptyString = () => SequenceComplexity.CalculateKmerEntropy(string.Empty, k);

        viaString.Should().Throw<ArgumentOutOfRangeException>("k < 1 is invalid (§3.3)");
        viaDna.Should().Throw<ArgumentOutOfRangeException>("k < 1 is invalid (§3.3)");
        viaEmptyString.Should().Throw<ArgumentOutOfRangeException>(
            "the k < 1 guard fires before the empty short-circuit (source-verified)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE / INJ — non-ACGT characters
    // ───────────────────────────────────────────────────────────────────

    #region KMER — BE/INJ — non-ACGT characters

    /// <summary>
    /// BE/INJ: non-ACGT characters. The STRICT typed surface rejects them at the
    /// DnaSequence ctor with ArgumentException (so the metric never sees garbage),
    /// while the LENIENT string surface is alphabet-agnostic: it upper-cases and tallies
    /// each char as an opaque k-mer symbol, never throwing (§3.3, INV-04). We pin both:
    /// a non-DNA homopolymer "NNNNNN" scores H = 0 like "AAAAAA" (symbol identity is
    /// irrelevant), case-folding gives an identical result (INV-04), and a non-DNA
    /// all-distinct-k-mer string scores its uniform maximum log₂(N).
    /// </summary>
    [Test]
    public void Kmer_NonAcgt_LenientStringScoresStructurally_StrictRejects()
    {
        // Lenient surface: alphabet-agnostic, non-DNA homopolymer behaves like a DNA one.
        SequenceComplexity.CalculateKmerEntropy("NNNNNN", 2)
            .Should().BeApproximately(0.0, Tolerance,
                "the core is symbol-agnostic: \"N\"×6 has one distinct dimer ⇒ H = 0 (INV-02)");
        SequenceComplexity.CalculateKmerEntropy("------", 2)
            .Should().BeApproximately(0.0, Tolerance, "gap homopolymer scores identically (§3.3)");

        // INV-04: case-invariance — lower-case input is upper-cased before counting.
        SequenceComplexity.CalculateKmerEntropy("atatat", 2)
            .Should().BeApproximately(SequenceComplexity.CalculateKmerEntropy("ATATAT", 2), Tolerance,
                "result is invariant to letter case (INV-04)");

        // Lenient: a non-DNA all-distinct-k-mer string scores its uniform max log₂(N).
        double digits = SequenceComplexity.CalculateKmerEntropy("123456", 2);
        digits.Should().BeApproximately(ReferenceKmerEntropy("123456", 2), Tolerance,
            "matches the independent reference");
        digits.Should().BeApproximately(Math.Log2(5), Tolerance,
            "12,23,34,45,56 all distinct over N=5 ⇒ H = log₂5 (INV-03)");

        // Strict surface: non-ACGT is rejected at the ctor, never reaching the metric.
        var build = () => new DnaSequence("ACGTNACGT");
        build.Should().Throw<ArgumentException>("the DnaSequence ctor rejects non-ACGT (strict gate)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  BE / RB — randomized boundary sweep + scale guard
    // ───────────────────────────────────────────────────────────────────

    #region KMER — BE/RB — randomized sweep + scale guard

    /// <summary>
    /// BE/RB: a randomized boundary sweep. Over many fixed-seed inputs of varied length
    /// AND varied k, deliberately including degenerate boundaries (len &lt; k,
    /// homopolymers, periodic strings), the k-mer entropy must ALWAYS be well-formed:
    /// finite (no NaN/Infinity), within the Shannon bounds 0 ≤ H ≤ log₂(N) (INV-01),
    /// exactly equal to the independent Shannon reference, 0 whenever len &lt; k or the
    /// sequence is a homopolymer (INV-02). Never a throw, hang, or overflow.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Kmer_RandomizedSweep_AlwaysWellFormed()
    {
        var rng = new Random(Seed + 8);

        for (int iter = 0; iter < 400; iter++)
        {
            int n = rng.Next(0, 600);                 // includes 0, 1, 2 …
            int k = rng.Next(1, 8);                    // valid k ≥ 1
            int mode = rng.Next(3);
            string s = mode switch
            {
                0 => RandomDna(rng, n),                                 // random DNA
                1 => new string("ACGT"[rng.Next(4)], n),                // homopolymer (H = 0)
                _ => string.Concat(Enumerable.Repeat("AC", n / 2)),     // periodic dinucleotide
            };

            double h = SequenceComplexity.CalculateKmerEntropy(s, k);

            double.IsNaN(h).Should().BeFalse("k-mer entropy must never be NaN");
            double.IsInfinity(h).Should().BeFalse("k-mer entropy must never be Infinity");
            h.Should().BeGreaterThanOrEqualTo(0.0, "H ≥ 0 (INV-01)");
            h.Should().BeApproximately(ReferenceKmerEntropy(s, k), Tolerance,
                "the production entropy must match the independent Shannon reference");

            int nKmers = s.Length - k + 1;
            if (nKmers >= 1)
                h.Should().BeLessThanOrEqualTo(Math.Log2(nKmers) + Tolerance,
                    "H ≤ log₂(N), N = L − k + 1 (INV-01)");

            if (s.Length < k)
                h.Should().Be(0.0, "len < k ⇒ no k-mer ⇒ H = 0 (§6.1)");

            // INV-02: a homopolymer is the minimal-entropy (0) extreme at every length/k.
            if (mode == 1 && s.Length >= k)
                h.Should().BeApproximately(0.0, Tolerance,
                    "a homopolymer has one distinct k-mer ⇒ H = 0 (INV-02)");
        }
    }

    /// <summary>
    /// BE/INJ: the lenient string surface must NEVER throw on pure random-byte garbage —
    /// arbitrary BMP code points including control chars, the null byte, and lone
    /// surrogate halves. Each char is parsed as an opaque k-mer symbol after upper-casing
    /// (§3.3), so the result must be a finite non-negative value matching the independent
    /// reference parse of the UPPER-CASED input, with no encoding surprise.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Kmer_LenientStringSurface_RandomBmpGarbage_NeverThrows()
    {
        var rng = new Random(Seed + 9);

        for (int iter = 0; iter < 300; iter++)
        {
            string s = RandomBmpChars(rng, rng.Next(0, 300));
            int k = rng.Next(1, 6);

            double h = 0;
            var act = () => { h = SequenceComplexity.CalculateKmerEntropy(s, k); };
            act.Should().NotThrow("the lenient string parse is alphabet-agnostic and total over char (§3.3)");

            double.IsNaN(h).Should().BeFalse("never NaN on garbage");
            double.IsInfinity(h).Should().BeFalse("never Infinity on garbage");
            h.Should().BeGreaterThanOrEqualTo(0.0, "H ≥ 0 (INV-01)");
            h.Should().BeApproximately(ReferenceKmerEntropy(s.ToUpperInvariant(), k), Tolerance,
                "must equal the reference entropy of the upper-cased input (§3.3)");
        }
    }

    /// <summary>
    /// BE/OVF: a very long sequence (200,000 bases) must compute without overflow or
    /// hang under a CancelAfter guard (the scan is O(N·k), §4.3). A long homopolymer
    /// (the minimal-entropy extreme, H = 0) and a long random sequence must each yield a
    /// finite entropy within [0, log₂ N]; random must be far more complex than the
    /// homopolymer (INV-02) and sit high within its log₂(N) ceiling, with the homopolymer
    /// pinned to exactly 0.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Kmer_VeryLong_NoOverflowNoHang_AndContractHolds()
    {
        var rng = new Random(Seed + 10);
        const int n = 200_000;
        const int k = 6;

        string homopolymer = new string('A', n);
        string random = RandomDna(rng, n);

        double homoH = SequenceComplexity.CalculateKmerEntropy(homopolymer, k);
        double randomH = SequenceComplexity.CalculateKmerEntropy(random, k);

        double.IsNaN(homoH).Should().BeFalse();
        double.IsNaN(randomH).Should().BeFalse();
        double.IsInfinity(homoH).Should().BeFalse();
        double.IsInfinity(randomH).Should().BeFalse();

        homoH.Should().BeApproximately(0.0, Tolerance,
            "a homopolymer has one distinct k-mer at scale ⇒ H = 0 (INV-02)");
        randomH.Should().BeGreaterThan(homoH,
            "random DNA spreads probability over many k-mers ⇒ far higher entropy (INV-02)");
        randomH.Should().BeLessThanOrEqualTo(Math.Log2(n - k + 1) + Tolerance,
            "H ≤ log₂(N) even at scale (INV-01)");
        randomH.Should().BeGreaterThan(2.0,
            "random DNA at scale sits high within [0, log₂(N≈2e5)≈17.6 bits]");
    }

    #endregion

    #endregion
}
