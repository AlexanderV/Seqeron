using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Epigenetics area — in-silico bisulfite conversion
/// (EPIGEN-BISULF-001). The single public entry point under test lives in
/// <see cref="EpigeneticsAnalyzer"/>:
///   • <see cref="EpigeneticsAnalyzer.SimulateBisulfiteConversion(string, IReadOnlySet{int})"/>
///     — simulates sodium-bisulfite treatment of ONE DNA strand: every
///     unprotected cytosine is converted to thymine (uracil, read as T), every
///     5-methyl (protected) cytosine whose 0-based index is in
///     <c>methylatedPositions</c> stays a cytosine, and every non-cytosine base
///     is returned unchanged.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (IndexOutOfRange /
/// NullReference / ArgumentOutOfRange). Every input must resolve to EITHER a
/// well-defined, theory-correct value OR a *documented, intentional* outcome.
/// For the bisulfite simulator the headline hazards are:
///   • an output whose length differs from the input (a dropped/duplicated base
///     — violates INV-01);
///   • the WRONG base being converted (a non-C altered, or a methylated C
///     converted — violates INV-02/INV-03);
///   • a NullReferenceException when <c>methylatedPositions</c> is null
///     (it is documented as "null = none protected");
///   • a crash on the empty/null sequence instead of the documented "".
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-BISULF-001 — Bisulfite Sequencing Analysis (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 182.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///     Targets (checklist row 182): "no C, all-C, all-methylated, empty".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The bisulfite-conversion contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Conversion rule (Frommer et al. 1992), per base b at position i of one strand:
///   • b is an unprotected cytosine        → thymine (C→T, c→t; case preserved);
///   • b is a protected (5-methyl) cytosine (i ∈ methylatedPositions) → unchanged C;
///   • otherwise                           → b unchanged.
///   — docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md §2.2, §4.1.
///
/// Documented invariants exercised here:
///   • INV-01: output length == input length (per-base substitution only).
///   • INV-02: protected cytosines stay C; unprotected C→T (case preserved).
///   • INV-03: non-cytosine bases are unchanged.
///   — docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md §2.4.
///
/// Boundary / degenerate handling (docs §3.1, §3.3, §6.1):
///   • null / empty sequence              → "" (no crash).
///   • no cytosines                       → returned unchanged (nothing to convert).
///   • methylatedPositions == null        → none protected (every C converts).
///   • lowercase c                        → t (case preserved).
///   • indices in methylatedPositions that do NOT point at a C are inert
///     (only cytosines react at all).
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticBisulfiteFuzzTests
{
    // ── Well-formed-result assertion helper ─────────────────────────────────
    // Pins INV-01..INV-03 against an independent re-derivation of the documented
    // rule. For every position, the output is checked base-by-base:
    //   • length must equal the input length (INV-01);
    //   • a non-cytosine input base must appear unchanged (INV-03);
    //   • a cytosine whose index is protected must remain a cytosine (INV-02);
    //   • an unprotected cytosine must be T (upper C) or t (lower c) (INV-02);
    //   • NO output position may differ from the input except an unprotected
    //     cytosine, and the only legal change is C→T / c→t.
    // This is what stops a test from rubber-stamping a corrupted strand green.
    private static void AssertWellFormedConversion(
        string input,
        string output,
        IReadOnlySet<int>? methylatedPositions)
    {
        // INV-01: per-base substitution preserves length exactly.
        output.Length.Should().Be(input.Length,
            because: "INV-01: conversion is a per-base substitution");

        for (int i = 0; i < input.Length; i++)
        {
            char inBase = input[i];
            char outBase = output[i];
            bool isUpperC = inBase == 'C';
            bool isLowerC = inBase == 'c';
            bool isProtected = methylatedPositions?.Contains(i) == true;

            if (!isUpperC && !isLowerC)
            {
                // INV-03: non-cytosine bases are never touched.
                outBase.Should().Be(inBase,
                    because: $"INV-03: non-cytosine base at {i} must be unchanged");
            }
            else if (isProtected)
            {
                // INV-02: a protected cytosine stays exactly the same cytosine.
                outBase.Should().Be(inBase,
                    because: $"INV-02: protected cytosine at {i} stays C");
            }
            else
            {
                // INV-02: an unprotected cytosine converts C→T / c→t (case kept).
                char expected = isUpperC ? 'T' : 't';
                outBase.Should().Be(expected,
                    because: $"INV-02: unprotected cytosine at {i} converts (case preserved)");
            }
        }
    }

    #region EPIGEN-BISULF-001 — Bisulfite conversion (SimulateBisulfiteConversion)

    // ── BE: empty / null sequence → documented "" ───────────────────────────

    [Test]
    public void Convert_NullSequence_ReturnsEmptyString_NoCrash()
    {
        // Docs §3.1 / §6.1: "Null or empty sequence → ''".
        Action act = () => SimulateBisulfiteConversion(null!);
        act.Should().NotThrow();
        SimulateBisulfiteConversion(null!).Should().Be("");
    }

    [Test]
    public void Convert_EmptySequence_ReturnsEmptyString()
    {
        SimulateBisulfiteConversion("").Should().Be("");
    }

    [Test]
    public void Convert_NullSequence_WithMethylationSet_ReturnsEmptyString()
    {
        // An empty/null sequence short-circuits before the methylation set is
        // ever consulted; a populated set must not change the "" result.
        SimulateBisulfiteConversion(null!, new HashSet<int> { 0, 1, 2 }).Should().Be("");
        SimulateBisulfiteConversion("", new HashSet<int> { 0, 1, 2 }).Should().Be("");
    }

    // ── BE: no cytosine → returned unchanged (nothing to convert) ────────────

    [Test]
    public void Convert_NoCytosine_ReturnsInputUnchanged()
    {
        // Docs §6.1: "No cytosines → returned unchanged". Only C reacts with
        // bisulfite, so an A/G/T-only strand passes through verbatim.
        string seq = "ATGATGATG";
        string converted = SimulateBisulfiteConversion(seq);
        converted.Should().Be(seq);
        AssertWellFormedConversion(seq, converted, null);
    }

    [Test]
    public void Convert_AllAdenine_ReturnsInputUnchanged()
    {
        string seq = new string('A', 50);
        SimulateBisulfiteConversion(seq).Should().Be(seq);
    }

    [Test]
    public void Convert_NoCytosineLowercase_ReturnsInputUnchanged()
    {
        // Lowercase non-C bases are equally inert (case preserved, not altered).
        string seq = "atgatgatg";
        SimulateBisulfiteConversion(seq).Should().Be(seq);
    }

    [Test]
    public void Convert_NoCytosine_WithMethylationSet_StillUnchanged()
    {
        // No cytosine exists at the "protected" indices, so a populated set is
        // inert: only cytosines are ever affected at all (INV-03).
        string seq = "ATGATGATG";
        var meth = new HashSet<int> { 0, 1, 2, 5, 8 };
        SimulateBisulfiteConversion(seq, meth).Should().Be(seq);
    }

    // ── BE: all-C, all unmethylated → every C converts to T ──────────────────

    [Test]
    public void Convert_AllC_AllUnmethylated_AllConvertToT()
    {
        // Frommer rule: with nothing protected, every cytosine converts.
        string seq = new string('C', 40);
        string converted = SimulateBisulfiteConversion(seq);
        converted.Should().Be(new string('T', 40));
        AssertWellFormedConversion(seq, converted, null);
    }

    [Test]
    public void Convert_AllLowercaseC_AllUnmethylated_AllConvertToLowercaseT()
    {
        // Docs §6.1: "Lowercase c → t (case preserved)".
        string seq = new string('c', 40);
        SimulateBisulfiteConversion(seq).Should().Be(new string('t', 40));
    }

    [Test]
    public void Convert_AllC_EmptyMethylationSet_AllConvertToT()
    {
        // An explicitly empty set is equivalent to null = nothing protected.
        string seq = new string('C', 25);
        SimulateBisulfiteConversion(seq, new HashSet<int>()).Should().Be(new string('T', 25));
    }

    // ── BE: all-C, all methylated → fully protected, output == input ─────────

    [Test]
    public void Convert_AllC_AllMethylated_OutputEqualsInput()
    {
        // INV-02 protected branch: when EVERY cytosine is methylated, none react,
        // so the converted strand is identical to the input (the protected case).
        string seq = new string('C', 40);
        var meth = new HashSet<int>(Enumerable.Range(0, seq.Length));
        string converted = SimulateBisulfiteConversion(seq, meth);
        converted.Should().Be(seq);
        AssertWellFormedConversion(seq, converted, meth);
    }

    [Test]
    public void Convert_AllLowercaseC_AllMethylated_OutputEqualsInput()
    {
        string seq = new string('c', 30);
        var meth = new HashSet<int>(Enumerable.Range(0, seq.Length));
        SimulateBisulfiteConversion(seq, meth).Should().Be(seq);
    }

    [Test]
    public void Convert_MixedSequence_AllCytosinesMethylated_OutputEqualsInput()
    {
        // Protecting every cytosine position in a mixed strand must yield the
        // input verbatim — no base of any kind is altered.
        string seq = "ACGTCGAACCGGTTAACC";
        var methCols = Enumerable.Range(0, seq.Length)
            .Where(i => seq[i] is 'C' or 'c')
            .ToHashSet();
        string converted = SimulateBisulfiteConversion(seq, methCols);
        converted.Should().Be(seq);
        AssertWellFormedConversion(seq, converted, methCols);
    }

    // ── POSITIVE sanity: hand-computed mixed methylation pattern ─────────────

    [Test]
    public void Convert_DocumentedWorkedExample_MatchesExactly()
    {
        // Docs §7.1 worked example:
        //   SimulateBisulfiteConversion("ACGTCGAA", {1}) == "ACGTTGAA".
        // Index map:   A C G T C G A A
        //              0 1 2 3 4 5 6 7
        // C@1 is protected (stays C); C@4 is unprotected (→ T); rest unchanged.
        string seq = "ACGTCGAA";
        var meth = new HashSet<int> { 1 };
        string converted = SimulateBisulfiteConversion(seq, meth);
        converted.Should().Be("ACGTTGAA");
        AssertWellFormedConversion(seq, converted, meth);
    }

    [Test]
    public void Convert_KnownMix_ConvertsExactlyPerRule()
    {
        // Hand-computed: "CCGGCC" with C@0 and C@5 protected.
        // index:   C  C  G  G  C  C
        //          0  1  2  3  4  5
        // protected {0,5}: C@0 stays C, C@1 → T, G@2 G@3 unchanged, C@4 → T, C@5 stays C.
        // expected: C  T  G  G  T  C  = "CTGGTC".
        string seq = "CCGGCC";
        var meth = new HashSet<int> { 0, 5 };
        string converted = SimulateBisulfiteConversion(seq, meth);
        converted.Should().Be("CTGGTC");
        AssertWellFormedConversion(seq, converted, meth);
    }

    [Test]
    public void Convert_MixedCase_PreservesCaseWhileConverting()
    {
        // "AcGCgT": positions of cytosines are 1 (c) and 3 (C) and 4 (g is NOT c).
        // index:   A  c  G  C  g  T
        //          0  1  2  3  4  5
        // nothing protected: c@1 → t, C@3 → T; A,G,g,T unchanged.
        // expected: A  t  G  T  g  T  = "AtGTgT".
        string seq = "AcGCgT";
        string converted = SimulateBisulfiteConversion(seq);
        converted.Should().Be("AtGTgT");
        AssertWellFormedConversion(seq, converted, null);
    }

    [Test]
    public void Convert_NonCBaseInMethylationSet_IsInert()
    {
        // A protected index that does NOT point at a cytosine cannot shield a
        // non-existent C and cannot alter the (non-C) base sitting there.
        // "ATGCAT": only C is at index 3; protect index 0 (A) and 3 (C).
        // index:   A  T  G  C  A  T
        //          0  1  2  3  4  5
        // protect {0,3}: A@0 unchanged (non-C, set inert), C@3 protected → stays C.
        // expected: A  T  G  C  A  T  = unchanged.
        string seq = "ATGCAT";
        var meth = new HashSet<int> { 0, 3 };
        string converted = SimulateBisulfiteConversion(seq, meth);
        converted.Should().Be("ATGCAT");
        AssertWellFormedConversion(seq, converted, meth);
    }

    [Test]
    public void Convert_NonCanonicalBases_PassThroughUnchanged()
    {
        // N / gaps / IUPAC codes are not cytosine, so INV-03 leaves them intact;
        // only the real cytosine converts.
        string seq = "NNCNN-RYC";
        // cytosines at indices 2 and 7, nothing protected → both → T, rest as-is.
        string converted = SimulateBisulfiteConversion(seq);
        converted.Should().Be("NNTNN-RYT");
        AssertWellFormedConversion(seq, converted, null);
    }

    // ── BE: single-base boundaries ───────────────────────────────────────────

    [Test]
    public void Convert_SingleUnmethylatedC_ConvertsToT()
    {
        SimulateBisulfiteConversion("C").Should().Be("T");
        SimulateBisulfiteConversion("c").Should().Be("t");
    }

    [Test]
    public void Convert_SingleMethylatedC_StaysC()
    {
        SimulateBisulfiteConversion("C", new HashSet<int> { 0 }).Should().Be("C");
        SimulateBisulfiteConversion("c", new HashSet<int> { 0 }).Should().Be("c");
    }

    [Test]
    public void Convert_SingleNonCytosine_Unchanged()
    {
        SimulateBisulfiteConversion("A").Should().Be("A");
        SimulateBisulfiteConversion("G").Should().Be("G");
        SimulateBisulfiteConversion("T").Should().Be("T");
    }

    [Test]
    public void Convert_OutOfRangeMethylationIndices_AreIgnored_NoCrash()
    {
        // Indices outside [0, len) (including negative and huge values) cannot
        // protect any in-range cytosine and must not crash. Every in-range C is
        // therefore unprotected and converts.
        string seq = "CCCC";
        var meth = new HashSet<int> { -1, 100, int.MaxValue, int.MinValue };
        Action act = () => SimulateBisulfiteConversion(seq, meth);
        act.Should().NotThrow();
        SimulateBisulfiteConversion(seq, meth).Should().Be("TTTT");
    }

    // ── BE/robustness: random fuzz — never crash, always well-formed ─────────

    [Test]
    [CancelAfter(30000)]
    public void Convert_RandomSequences_NeverCrash_AlwaysWellFormed()
    {
        // Random strands over an alphabet that mixes cytosine (both cases) with
        // non-cytosine and non-canonical bases, plus a random protected subset.
        // Every result must satisfy INV-01..INV-03 against an independent
        // re-derivation of the documented rule.
        const string alphabet = "ACGTacgtNn-RY";
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 400);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            // Random protected subset: include some out-of-range indices too.
            var meth = new HashSet<int>();
            int protectedCount = len == 0 ? 0 : rng.Next(0, len + 1);
            for (int k = 0; k < protectedCount; k++)
                meth.Add(rng.Next(-5, len + 5));

            // Roughly half the time pass null instead of an empty/built set.
            IReadOnlySet<int>? methArg = (meth.Count == 0 && rng.Next(2) == 0) ? null : meth;

            string converted = null!;
            Action act = () => converted = SimulateBisulfiteConversion(seq, methArg);
            act.Should().NotThrow($"seed={seed}, len={len}");

            AssertWellFormedConversion(seq, converted, methArg);
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void Convert_NullVsEmptyMethylationSet_AreEquivalent()
    {
        // "null = none protected" must be EXACTLY equivalent to an empty set for
        // any input, over a fuzzed corpus.
        const string alphabet = "ACGTacgtN";
        var empty = new HashSet<int>();
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 200);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            string viaNull = SimulateBisulfiteConversion(seq, null);
            string viaEmpty = SimulateBisulfiteConversion(seq, empty);
            viaNull.Should().Be(viaEmpty, $"seed={seed}: null must equal empty set");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void Convert_OutputLengthAlwaysEqualsInput_FuzzedLengths()
    {
        // INV-01 stress: a wide range of lengths, including very long, must
        // always preserve length exactly (no buffer drift on the StringBuilder).
        const string alphabet = "ACGTacgtNn";
        int[] lengths = { 0, 1, 2, 3, 7, 16, 63, 64, 65, 1000, 100_000 };
        var rng = new Random(12345);
        foreach (int len in lengths)
        {
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            SimulateBisulfiteConversion(seq).Length.Should().Be(len, $"len={len}");
        }
    }

    #endregion
}
