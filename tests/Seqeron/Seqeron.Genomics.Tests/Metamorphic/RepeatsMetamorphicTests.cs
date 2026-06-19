using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Repeats area.
///
/// Each test encodes a metamorphic relation (MR) — a property that relates the
/// outputs of multiple executions under input transformations, without requiring
/// a hardcoded test oracle. The relations are derived from the *definition* of a
/// DNA palindrome (a segment equal to its own reverse complement), not from the
/// current implementation's output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-PALIN-001 — DNA palindrome detection (Repeats)
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 17.
/// Source: docs/algorithms/Repeat_Analysis/Palindrome_Detection.md;
///         RepeatFinder.FindPalindromes(DnaSequence, int minLength=4, int maxLength=12)
///         returning PalindromeResult(int Position, string Sequence, int Length).
///
/// API semantics confirmed from the algorithm doc + RepeatFinder.cs:
///   - A DNA palindrome S satisfies S == ReverseComplement(S) (NOT textual symmetry).
///   - Only EVEN candidate lengths are scanned (len steps by 2); minLength must be
///     even and >= 4, maxLength must be >= minLength.
///   - Overlapping palindromes of different even lengths can both be reported.
///   - No gapped/loop palindromes (loop = 0; a palindrome is the zero-loop special
///     case of an inverted repeat — inverted repeats are a SEPARATE unit, REP-INV-001).
///
/// Relations verified:
///   - INV: every reported segment equals its own reverse complement; known
///          restriction-site palindromes (EcoRI/BamHI) are detected.
///   - MON: widening [minLength, maxLength] yields a superset of results.
///   - SHIFT: prepending a flank shifts every position by |F| and preserves the set.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class RepeatsMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20240617);

    /// <summary>Generates a random DNA string of the given length (fixed seed).</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(4)];
        return new string(chars);
    }

    private static List<PalindromeResult> Palindromes(string seq, int minLen, int maxLen) =>
        RepeatFinder.FindPalindromes(new DnaSequence(seq), minLen, maxLen).ToList();

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-PALIN-001 — DNA palindrome detection
    // ═══════════════════════════════════════════════════════════════════

    #region MR-INV: every palindrome equals its own reverse complement

    /// <summary>
    /// INV (defining invariant): for EVERY detected palindrome, the reported segment
    /// equals its own reverse complement. This is the consistency invariant that
    /// makes the result a DNA palindrome at all, and it must hold on every result
    /// across hand-built and random inputs.
    /// </summary>
    [Test]
    public void FindPalindromes_EveryResult_EqualsOwnReverseComplement()
    {
        var sequences = new[]
        {
            "GAATTCGGATCCAAGCTT",                 // EcoRI + BamHI + HindIII back to back
            "ACGTAAATTTACGT",                     // mixed
            "GGGGCCCCATATATGCGC",                 // several even windows
            RandomDna(60),
            RandomDna(120),
            RandomDna(200),
        };

        foreach (var seq in sequences)
        {
            var results = Palindromes(seq, 4, 12);

            foreach (var p in results)
            {
                // The reported segment is genuinely self-reverse-complementary.
                string revComp = DnaSequence.GetReverseComplementString(p.Sequence);
                p.Sequence.Should().Be(revComp,
                    because: $"a DNA palindrome must equal its own reverse complement; "
                           + $"'{p.Sequence}' at {p.Position} violates the definition");

                // Cross-check against the source text at the reported position/length.
                seq.Substring(p.Position, p.Length).Should().Be(p.Sequence,
                    because: $"reported Sequence must match the window at Position={p.Position}, Length={p.Length}");

                // INV-02 from the doc: every reported length is even.
                (p.Length % 2).Should().Be(0,
                    because: "DNA palindromes require even length so every base pairs across the axis");
            }
        }
    }

    /// <summary>
    /// INV: a hand-built sequence containing a known restriction-site palindrome
    /// (EcoRI = GAATTC, BamHI = GGATCC) is detected at the expected position with
    /// the expected segment. Theory pins both the position and the segment here.
    /// </summary>
    [Test]
    public void FindPalindromes_KnownRestrictionSites_AreDetected()
    {
        // GAATTC (EcoRI) starts at index 4; GGATCC (BamHI) starts at index 14.
        const string seq = "TTTTGAATTCGGGGGGATCCAAAA";

        var results = Palindromes(seq, 4, 12);

        results.Should().Contain(p => p.Sequence == "GAATTC" && p.Position == 4,
            because: "EcoRI site GAATTC is self-reverse-complementary and present at index 4");
        results.Should().Contain(p => p.Sequence == "GGATCC" && p.Position == 14,
            because: "BamHI site GGATCC is self-reverse-complementary and present at index 14");

        // Sanity: GAATTC and GGATCC really are their own reverse complements.
        DnaSequence.GetReverseComplementString("GAATTC").Should().Be("GAATTC");
        DnaSequence.GetReverseComplementString("GGATCC").Should().Be("GGATCC");
    }

    #endregion

    #region MR-MON: wider [minLength, maxLength] range → ≥ palindromes (superset)

    /// <summary>
    /// MON: widening the length window [minLength, maxLength] can only add results.
    /// Results(narrow) ⊆ Results(wide) by (Position, Sequence), and count is
    /// non-decreasing. The narrow window [4,6] is contained in the wide window [4,12].
    /// </summary>
    [Test]
    public void FindPalindromes_WiderLengthRange_SupersetOfNarrower()
    {
        var sequences = new[]
        {
            "GAATTCGGATCCGCGGCCGCAAGCTT",   // includes NotI (GCGGCCGC, len 8) among len-4/6 sites
            "ACGTGGGGCCCCATATGCGCGC",
            RandomDna(80),
            RandomDna(150),
        };

        foreach (var seq in sequences)
        {
            var narrow = Palindromes(seq, 4, 6);
            var wide = Palindromes(seq, 4, 12);

            wide.Count.Should().BeGreaterThanOrEqualTo(narrow.Count,
                because: "a wider length window scans every candidate length the narrow one did, plus more");

            var wideSet = wide.Select(p => (p.Position, p.Sequence)).ToHashSet();
            foreach (var p in narrow)
            {
                wideSet.Should().Contain((p.Position, p.Sequence),
                    because: $"palindrome '{p.Sequence}' at {p.Position} from [4,6] must persist in [4,12]");
            }
        }
    }

    /// <summary>
    /// MON: lengthening the upper bound from a value that excludes an 8-bp palindrome
    /// to one that includes it strictly grows the set with that palindrome.
    /// </summary>
    [Test]
    public void FindPalindromes_RaisingMaxLength_AddsLongerPalindrome()
    {
        // GCGGCCGC (NotI) is an 8-bp palindrome.
        const string seq = "TTTTGCGGCCGCAAAA";

        var upTo6 = Palindromes(seq, 4, 6);
        var upTo8 = Palindromes(seq, 4, 8);

        upTo8.Should().Contain(p => p.Sequence == "GCGGCCGC" && p.Position == 4,
            because: "raising maxLength to 8 must surface the 8-bp NotI palindrome");
        upTo6.Should().NotContain(p => p.Length == 8,
            because: "maxLength=6 cannot report any 8-bp window");

        // Superset still holds: everything from [4,6] survives in [4,8].
        var widerSet = upTo8.Select(p => (p.Position, p.Sequence)).ToHashSet();
        foreach (var p in upTo6)
            widerSet.Should().Contain((p.Position, p.Sequence),
                because: "wider maxLength is a superset of the narrower scan");
    }

    #endregion

    #region MR-SHIFT: prepending a flank shifts positions by |F|, preserves the set

    /// <summary>
    /// SHIFT: prepending a flank F that creates no palindrome at the junction shifts
    /// every palindrome position by exactly |F| and preserves the count and segments.
    /// The flank is a homopolymer run ("AAAAAA"): a run of identical bases can never
    /// be self-reverse-complementary (revcomp of A…A is T…T), and the A|G junction
    /// with the interior introduces no new even self-complementary window here.
    /// </summary>
    [Test]
    public void FindPalindromes_PrependFlank_ShiftsPositionsByFlankLength()
    {
        const string flank = "AAAAAA"; // homopolymer: never palindromic, neutral junction
        var sequences = new[]
        {
            "GAATTCGGGGGGATCC",
            "GCGGCCGCATATCCCGGG",
            RandomDna(70),
            RandomDna(130),
        };

        foreach (var seq in sequences)
        {
            var original = Palindromes(seq, 4, 12);
            var shifted = Palindromes(flank + seq, 4, 12);

            // Count preserved: flank adds no palindromes and destroys none.
            shifted.Count.Should().Be(original.Count,
                because: $"a non-palindrome-creating flank of length {flank.Length} preserves the palindrome set");

            // Exact +|F| shift and same segments.
            var shiftedSet = shifted.Select(p => (p.Position, p.Sequence)).ToHashSet();
            foreach (var p in original)
            {
                shiftedSet.Should().Contain((p.Position + flank.Length, p.Sequence),
                    because: $"palindrome '{p.Sequence}' at {p.Position} must move to {p.Position + flank.Length}");
            }

            // And the shifted set introduces nothing that isn't an original shifted by |F|.
            var expected = original.Select(p => (p.Position + flank.Length, p.Sequence)).ToHashSet();
            shiftedSet.Should().BeEquivalentTo(expected,
                because: "the only change from a neutral flank is a uniform +|F| positional shift");
        }
    }

    #endregion
}
