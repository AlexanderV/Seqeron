using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area codon-frequency unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …) and no
/// NaN / Infinity leaking into a frequency. Every input must result in EITHER a
/// well-defined, theory-correct value, OR a *documented, intentional* validation
/// outcome. A raw runtime exception or a hang on garbage input is a bug, not a
/// passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-CODON-FREQ-001 — Codon frequencies (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 227.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — empty / length&lt;3, length not a multiple of
///          3 (trailing 1–2 bases), non-ACGT characters, lowercase, very long
///          input (plus the all-ambiguous zero-total denominator boundary).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The codon-frequency contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateCodonFrequencies(string, int readingFrame = 0)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 688–723), returning IReadOnlyDictionary&lt;string,double&gt;.
///
/// This is the STATISTICS-module codon-frequency (count / total counted codons)
/// of the Kazusa CUTG convention — distinct from Codon-area CODON-* and the
/// Annotation-area codon-usage units.
///
/// Documented behaviour (docs/algorithms/Statistics/Codon_Frequencies.md,
/// Test Unit ID SEQ-CODON-FREQ-001):
///   • §2.2 / §4.1: read NON-OVERLAPPING in-frame triplets c_k over
///     [f+3(k−1), f+3k); frequency(x) = count(x) / total, where total is the
///     number of triplets composed solely of {A,C,G,T} (§2.2). The denominator
///     is the count of COMPLETE, ALL-ACGT codons — NOT the sequence length.
///   • §3.3 / §6.1: null / empty / length &lt; 3 → EMPTY dictionary (no full codon
///     to count, no zero-total division).
///   • §3.3 / §6.1: only COMPLETE non-overlapping triplets are read; trailing 1–2
///     bases (length not a multiple of 3 from the frame) are IGNORED.
///   • §2.2 / §3.3 / §6.1 (INV-03): a triplet containing ANY non-ACGT base is
///     excluded from BOTH count and total (ambiguous codons excluded). U is a
///     non-ACGT base here — no T↔U conversion (§3.3). When EVERY triplet is
///     ambiguous (total = 0) the result is the EMPTY dictionary — no DivideByZero.
///   • §3.3 / §6.1 (INV-04): input is upper-cased (ToUpperInvariant) before
///     counting, so counting is case-insensitive.
///
/// Invariants pinned (Codon_Frequencies.md §2.4):
///   • INV-01: every frequency ∈ (0, 1]; only OBSERVED codons are keys
///     (count(x) ≥ 1, count(x) ≤ total).
///   • INV-02: Σ frequency(x) = 1 over all keys whenever total ≥ 1.
///   • INV-03: codons with any non-ACGT base never appear and never change total.
///   • INV-04: output is independent of input letter case.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceCodonFrequencyFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The codon alphabet: a triplet is counted only if ALL three bases are here.</summary>
    private const string CodonAlphabet = "ACGT";

    /// <summary>Generates a random string of arbitrary BMP code points (0x0000–0xFFFF) —
    /// control chars, the null byte, lone surrogate halves, unicode letters/digits:
    /// random-byte fuzz fodder for the codon counter.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>
    /// Independent oracle: re-derives the documented count/total table by reading
    /// non-overlapping in-frame triplets and excluding any triplet with a non-ACGT
    /// base (Codon_Frequencies.md §2.2). The case-fold is applied with the SAME
    /// ToUpperInvariant the doc cites (§3.3), but the counting/exclusion logic is
    /// written here from the spec, NOT copied from the implementation.
    /// </summary>
    private static Dictionary<string, double> OracleFrequencies(string seq, int frame)
    {
        var counts = new Dictionary<string, int>();
        int total = 0;
        if (string.IsNullOrEmpty(seq) || seq.Length < 3)
            return new Dictionary<string, double>();

        string upper = seq.ToUpperInvariant();
        for (int i = frame; i <= upper.Length - 3; i += 3)
        {
            string codon = upper.Substring(i, 3);
            if (codon.All(c => CodonAlphabet.Contains(c)))
            {
                counts[codon] = counts.GetValueOrDefault(codon) + 1;
                total++;
            }
        }

        var freq = new Dictionary<string, double>();
        foreach (var (codon, c) in counts)
            freq[codon] = (double)c / total;
        return freq;
    }

    /// <summary>
    /// Asserts the universal well-formedness contract that must hold for ANY input:
    /// every key is a length-3 ACGT codon, every frequency is finite (never NaN /
    /// Infinity), strictly positive and ≤ 1 (INV-01), and — per INV-02 — the
    /// frequencies sum to 1.0 whenever the dictionary is non-empty (total ≥ 1).
    /// </summary>
    private static void AssertWellFormed(IReadOnlyDictionary<string, double> freq)
    {
        freq.Should().NotBeNull("the method must always return a dictionary, never null");

        foreach (var (codon, f) in freq)
        {
            codon.Should().HaveLength(3, "a codon key is exactly three bases");
            codon.All(c => CodonAlphabet.Contains(c))
                .Should().BeTrue($"'{codon}' must be over the {{A,C,G,T}} alphabet (INV-03)");

            double.IsFinite(f).Should().BeTrue(
                $"frequency of '{codon}' must be finite — never NaN/Infinity");
            f.Should().BeGreaterThan(0.0,
                $"a counted codon '{codon}' has count ≥ 1, so frequency > 0 (INV-01)");
            f.Should().BeLessThanOrEqualTo(1.0 + Tolerance,
                $"a frequency '{codon}' is count/total ≤ 1 (INV-01)");
        }

        if (freq.Count > 0)
            freq.Values.Sum().Should().BeApproximately(1.0, Tolerance,
                "INV-02: frequencies over the counted codons sum to 1");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-CODON-FREQ-001 — codon frequency : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact result

    /// <summary>
    /// Positive baseline (not a boundary): the doc's worked example must yield the
    /// documented codon frequencies EXACTLY. "ATGATGAAA" from frame 0 reads the
    /// non-overlapping triplets ATG, ATG, AAA ⇒ total = 3, count(ATG) = 2,
    /// count(AAA) = 1, so f_ATG = 2/3 and f_AAA = 1/3 (summing to 1). Confirms the
    /// suite asserts the BUSINESS contract, not just absence of crashes.
    /// — Codon_Frequencies.md §7.1 worked example.
    /// </summary>
    [Test]
    public void Codon_WorkedExample_MatchesHandComputedFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA");

        freq.Should().HaveCount(2, "two distinct codons: ATG and AAA");
        freq["ATG"].Should().BeApproximately(2.0 / 3.0, Tolerance);
        freq["AAA"].Should().BeApproximately(1.0 / 3.0, Tolerance);

        freq.Values.Sum().Should().BeApproximately(1.0, Tolerance, "INV-02");
        AssertWellFormed(freq);
    }

    /// <summary>
    /// Positive: INV-05 numeric anchor independent of the worked example. A sequence
    /// built so one codon occurs 22 times out of 386 counted codons must give that
    /// codon frequency 22/386 = 0.0569948…, i.e. the CUTG per-thousand value 56.995
    /// divided by 1000 (Codon_Frequencies.md §2.4 INV-05). 386 codons = 1158 bases.
    /// </summary>
    [Test]
    public void Codon_CutgPerThousandAnchor_Frequency_IsCountOverTotal()
    {
        // 22 × GCT, then 364 × AAA → 386 codons total.
        var sb = string.Concat(Enumerable.Repeat("GCT", 22)) +
                 string.Concat(Enumerable.Repeat("AAA", 364));

        var freq = SequenceStatistics.CalculateCodonFrequencies(sb);

        freq["GCT"].Should().BeApproximately(22.0 / 386.0, Tolerance,
            "INV-05: count/total = CUTG per-thousand (56.995) ÷ 1000");
        freq["AAA"].Should().BeApproximately(364.0 / 386.0, Tolerance);
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: empty / null / length < 3

    /// <summary>
    /// BE: the empty string is the lower size boundary. With no full codon there is
    /// nothing to count, so the documented result is the EMPTY dictionary — NO
    /// DivideByZero on a zero total, no exception.
    /// — Codon_Frequencies.md §3.3 / §6.1 (null/empty/length&lt;3 → empty dictionary).
    /// </summary>
    [Test]
    public void Codon_EmptyString_IsEmptyAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateCodonFrequencies(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        var freq = act();
        freq.Should().BeEmpty("no complete codon to count ⇒ empty table");
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 698) — empty dictionary, no NullReferenceException.
    /// — Codon_Frequencies.md §3.3 (null → empty dictionary).
    /// </summary>
    [Test]
    public void Codon_Null_IsEmptyAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateCodonFrequencies(null!);

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().BeEmpty();
    }

    /// <summary>
    /// BE: any sequence shorter than one codon (length 1 or 2) is the critical
    /// length&lt;3 boundary — no complete triplet exists, so the result is EMPTY and
    /// the zero-total denominator is never reached (guards against forming a partial
    /// "codon" or dividing by a zero total).
    /// — Codon_Frequencies.md §6.1 (length&lt;3 → empty dictionary).
    /// </summary>
    [TestCase("A")]
    [TestCase("AT")]
    [TestCase("g")]
    [TestCase("gc")]
    [TestCase("NN")]
    public void Codon_ShorterThanOneCodon_IsEmpty(string seq)
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies(seq);

        freq.Should().BeEmpty("length < 3 has no complete codon");
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: length not a multiple of 3 (trailing 1–2 bases ignored)

    /// <summary>
    /// BE: a sequence whose length is not a multiple of 3 must IGNORE the trailing
    /// 1–2 bases (only complete non-overlapping triplets are read). "ATGATGAA"
    /// (length 8) reads ATG, ATG and discards the trailing "AA": total = 2,
    /// f_ATG = 1.0. "ATGAAAC" (length 7) reads ATG, AAA, discards trailing "C":
    /// total = 2, f_ATG = f_AAA = 1/2. The trailing bases never inflate the total.
    /// — Codon_Frequencies.md §3.3 / §6.1 (trailing 1–2 bases ignored).
    /// </summary>
    [Test]
    public void Codon_LengthNotMultipleOfThree_IgnoresTrailingBases()
    {
        // length 8 ≡ 2 (mod 3): trailing "AA" discarded.
        var f8 = SequenceStatistics.CalculateCodonFrequencies("ATGATGAA");
        f8.Should().HaveCount(1);
        f8["ATG"].Should().BeApproximately(1.0, Tolerance, "two ATG over total 2; AA discarded");
        AssertWellFormed(f8);

        // length 7 ≡ 1 (mod 3): trailing "C" discarded.
        var f7 = SequenceStatistics.CalculateCodonFrequencies("ATGAAAC");
        f7.Should().HaveCount(2);
        f7["ATG"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        f7["AAA"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        f7.Values.Sum().Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(f7);
    }

    /// <summary>
    /// BE: the trailing-base rule under a non-zero reading frame. From frame 1,
    /// "XATGAAAyz" semantics: triplets start at offset 1. "GATGAAATG" from frame 1
    /// reads bases [1..] = "ATGAAATG" as ATG, AAA, then trailing "TG" (only 2 bases)
    /// is discarded ⇒ total = 2, f_ATG = f_AAA = 1/2. Confirms the frame offset and
    /// the trailing-base discard compose correctly (no IndexOutOfRange at the tail).
    /// — Codon_Frequencies.md §2.2 (frame offset f) / §3.3 (trailing bases ignored).
    /// </summary>
    [Test]
    public void Codon_ReadingFrame_OffsetsStartAndDiscardsTrailing()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("GATGAAATG", readingFrame: 1);

        freq.Should().HaveCount(2, "from frame 1: ATG, AAA; trailing 'TG' discarded");
        freq["ATG"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        freq["AAA"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        AssertWellFormed(freq);
    }

    #endregion

    #region BE — Boundary: non-ACGT characters (ambiguous codons excluded)

    /// <summary>
    /// BE: a triplet containing ANY non-ACGT base is excluded from BOTH count and
    /// total (INV-03). "ATGNNNAAA" reads ATG, NNN, AAA: NNN is excluded ⇒ total = 2,
    /// f_ATG = f_AAA = 1/2 (NOT 1/3). The ambiguous codon never appears as a key and
    /// never enlarges the denominator.
    /// — Codon_Frequencies.md §2.2 / §6.1 / INV-03 (ambiguous codons excluded).
    /// </summary>
    [Test]
    public void Codon_AmbiguousTriplet_ExcludedFromCountAndTotal()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGNNNAAA");

        freq.Should().HaveCount(2, "NNN is excluded; only ATG and AAA counted");
        freq.Should().NotContainKey("NNN", "INV-03: ambiguous codons never appear");
        freq["ATG"].Should().BeApproximately(1.0 / 2.0, Tolerance,
            "total = 2 (NNN excluded), not 3");
        freq["AAA"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: a SINGLE non-ACGT base anywhere in a triplet voids the WHOLE triplet —
    /// "ANG" / "ATN" / "NTG" are all excluded, not partially counted. "AAAANGTTT"
    /// reads AAA, ANG, TTT ⇒ ANG (has N) excluded, total = 2, f_AAA = f_TTT = 1/2.
    /// Also covers RNA 'U' being treated as non-ACGT (no T↔U conversion, §3.3):
    /// "AUGAAA" reads AUG (has U → excluded) and AAA ⇒ {AAA: 1.0}.
    /// — Codon_Frequencies.md §2.2 / §3.3 (U is a non-ACGT base, no U→T).
    /// </summary>
    [Test]
    public void Codon_NonAcgtBaseVoidsWholeTriplet_AndUIsAmbiguous()
    {
        var f1 = SequenceStatistics.CalculateCodonFrequencies("AAAANGTTT");
        f1.Should().HaveCount(2, "ANG (contains N) is excluded entirely");
        f1["AAA"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        f1["TTT"].Should().BeApproximately(1.0 / 2.0, Tolerance);
        AssertWellFormed(f1);

        // RNA U is non-ACGT here: AUG excluded, only AAA counted.
        var f2 = SequenceStatistics.CalculateCodonFrequencies("AUGAAA");
        f2.Should().HaveCount(1);
        f2.Should().NotContainKey("AUG", "U is a non-ACGT base; no T↔U conversion (§3.3)");
        f2["AAA"].Should().BeApproximately(1.0, Tolerance);
        AssertWellFormed(f2);
    }

    /// <summary>
    /// BE: when EVERY triplet is ambiguous the total is 0 and the documented result
    /// is the EMPTY dictionary — the frequency-conversion loop never runs, so there
    /// is NO DivideByZero by the zero total. Holds at the length=3 boundary and at
    /// scale, and case-insensitively for lowercase 'n'.
    /// — Codon_Frequencies.md §6.1 (all triplets ambiguous, total = 0 → empty).
    /// </summary>
    [TestCase("NNN")]
    [TestCase("NNNNNN")]
    [TestCase("nnnnnnnnn")]
    [TestCase("XYZxyz")]
    public void Codon_AllAmbiguous_IsEmpty_NoDivideByZero(string seq)
    {
        var act = () => SequenceStatistics.CalculateCodonFrequencies(seq);

        act.Should().NotThrow("an all-ambiguous input must never throw (no zero-total division)");
        act().Should().BeEmpty("every triplet has a non-ACGT base ⇒ total = 0 ⇒ empty");
    }

    #endregion

    #region BE — Boundary: lowercase (case-insensitive, INV-04)

    /// <summary>
    /// BE: lowercase input must produce the SAME upper-cased codon keys and identical
    /// frequencies as the uppercase form (INV-04: output independent of case, via
    /// ToUpperInvariant, SequenceStatistics.cs line 701). "atgatgaaa" must equal the
    /// worked example: f_ATG = 2/3, f_AAA = 1/3, with keys upper-cased.
    /// — Codon_Frequencies.md §3.3 / §6.1 / INV-04.
    /// </summary>
    [Test]
    public void Codon_LowercaseInput_SameUpperCasedKeysAndFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("atgatgaaa");

        freq.Should().HaveCount(2);
        freq.Should().ContainKey("ATG", "input is upper-cased before counting (INV-04)");
        freq.Should().ContainKey("AAA");
        freq["ATG"].Should().BeApproximately(2.0 / 3.0, Tolerance);
        freq["AAA"].Should().BeApproximately(1.0 / 3.0, Tolerance);
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: mixed-case must be folded too — "AtGaTgAaA" equals the all-upper worked
    /// example exactly. Guards against a half-applied case-fold.
    /// — Codon_Frequencies.md INV-04.
    /// </summary>
    [Test]
    public void Codon_MixedCase_EqualsUpperCaseResult()
    {
        var mixed = SequenceStatistics.CalculateCodonFrequencies("AtGaTgAaA");
        var upper = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA");

        mixed.Should().HaveCount(upper.Count);
        foreach (var (codon, f) in upper)
            mixed[codon].Should().BeApproximately(f, Tolerance);
        AssertWellFormed(mixed);
    }

    #endregion

    #region BE — Boundary: very long sequence (O(n), no overflow)

    /// <summary>
    /// BE: a very long sequence must be processed in O(n) without overflow, hang or
    /// precision corruption. 300 000 bases = 100 000 codons, all "AAA": total =
    /// 100 000, f_AAA = 1.0 exactly. Pins the large-total denominator and bounds the
    /// runtime under [CancelAfter].
    /// — Codon_Frequencies.md §4.3 (O(n)); INV-02.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Codon_VeryLongHomopolymer_SingleCodonFrequencyOne()
    {
        const int codonCount = 100_000;
        string seq = string.Concat(Enumerable.Repeat("AAA", codonCount)); // 300 000 bases

        var freq = SequenceStatistics.CalculateCodonFrequencies(seq);

        freq.Should().HaveCount(1, "a homopolymer of AAA has exactly one distinct codon");
        freq["AAA"].Should().BeApproximately(1.0, Tolerance,
            $"all {codonCount} codons are AAA ⇒ count/total = 1");
        AssertWellFormed(freq);
    }

    /// <summary>
    /// BE: a very long MIXED sequence whose length is not a multiple of 3 and that
    /// carries interspersed ambiguous bases — the heaviest realistic fuzz shape. The
    /// result must match an independent count/total oracle (§2.2) exactly and remain
    /// finite, summing to 1, with no overflow or hang.
    /// — Codon_Frequencies.md §2.2 / §4.3.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Codon_VeryLongMixedSequence_MatchesOracle()
    {
        var rng = new Random(2026_06_21);
        const string alphabet = "ACGTNacgtn-"; // canonical, lowercase, N and junk
        int len = 299_999; // ≡ 2 (mod 3): trailing bases must be ignored
        var chars = new char[len];
        for (int i = 0; i < len; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        string seq = new string(chars);

        var freq = SequenceStatistics.CalculateCodonFrequencies(seq);
        var oracle = OracleFrequencies(seq, 0);

        freq.Should().HaveCount(oracle.Count, "same set of distinct counted codons");
        foreach (var (codon, f) in oracle)
            freq[codon].Should().BeApproximately(f, Tolerance, $"f_{codon} = count/total");
        AssertWellFormed(freq);
    }

    #endregion

    #region BE / RB — Random fuzz: never throws, always well-formed, matches oracle

    /// <summary>
    /// BE/RB: a large batch of arbitrary BMP strings (control chars, null byte, lone
    /// surrogate halves, unicode letters/digits) must NEVER throw and must ALWAYS
    /// produce a well-formed dictionary satisfying every invariant. Core fuzz
    /// guarantee: no DivideByZero, no IndexOutOfRange, no overflow, no NaN on garbage.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Codon_RandomGarbageStrings_NeverThrow_AlwaysWellFormed()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);
            int frame = rng.Next(0, 3);

            IReadOnlyDictionary<string, double> freq = null!;
            var act = () => freq = SequenceStatistics.CalculateCodonFrequencies(input, frame);

            act.Should().NotThrow($"garbage input (len {len}, frame {frame}) must never crash");
            AssertWellFormed(freq);
        }
    }

    /// <summary>
    /// BE: a randomized sweep over realistic sequences (canonical + lowercase + N +
    /// junk, random length and reading frame) must match the independent count/total
    /// oracle (§2.2) EXACTLY for every shape — the spec-faithful cross-check that the
    /// counts, the case-fold, the ambiguous-codon exclusion, the frame offset and the
    /// trailing-base discard all compose correctly, with frequencies summing to 1.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Codon_RandomRealisticSequences_MatchOracleAcrossFrames()
    {
        var rng = new Random(424242);
        const string alphabet = "ACGTacgtN-x";

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 120);
            int frame = rng.Next(0, 3);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            IReadOnlyDictionary<string, double> freq = null!;
            var act = () => freq = SequenceStatistics.CalculateCodonFrequencies(seq, frame);
            act.Should().NotThrow($"len {len}, frame {frame} must not crash");

            var oracle = OracleFrequencies(seq, frame);
            freq.Should().HaveCount(oracle.Count,
                $"same distinct codons (len {len}, frame {frame})");
            foreach (var (codon, f) in oracle)
                freq[codon].Should().BeApproximately(f, Tolerance,
                    $"f_{codon} = count/total (len {len}, frame {frame})");
            AssertWellFormed(freq);
        }
    }

    #endregion
}
