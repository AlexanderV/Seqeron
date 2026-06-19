using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Composition area.
///
/// Each test encodes a metamorphic relation (MR) — a property that relates the
/// outputs of multiple executions under input transformations, without requiring
/// a hardcoded test oracle. The relations are derived from the GC-content
/// *definition*, not from the current implementation's output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-GC-001 — GC content (Composition)
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 1.
/// Relations: INV complement preserves GC%; INV shuffle preserves GC%;
///            INV case-insensitive (+ derived INV reverse-complement,
///            ADD concatenation-additivity of the GC count).
///
/// Source (GC% definition):
///   GC% = (G + C) / (A + T + G + C) × 100
///   — docs/algorithms/Statistics/GC_Content_Profile.md §2.2 [Wikipedia GC-content];
///     docs/algorithms/Sequence_Composition/Sequence_Composition.md §2.2.
///   API: DnaSequence.GcContent() returns a percentage in [0, 100]
///        (src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs);
///        backed by SequenceExtensions.CalculateGcContent (case-insensitive,
///        counts G/C over A/T/G/C/U).
///
/// Why these relations hold (from the definition above):
///   • Complement maps A↔T and C↔G. A 'G' becomes a 'C' and a 'C' becomes a 'G',
///     so the (G + C) count is invariant; A↔T leaves (A + T) invariant. With both
///     numerator and denominator unchanged, GC% is invariant.
///   • GC% depends only on base *counts* (composition), not on order, so any
///     permutation (shuffle) of the bases leaves GC% unchanged.
///   • Counting is case-folded, so GC% is identical for lower-, upper-, and
///     mixed-case spellings of the same sequence.
///   • Reverse-complement = reverse (a permutation, GC%-invariant) ∘ complement
///     (GC%-invariant), so it preserves GC% too.
///   • The GC count is additive over concatenation: count(a + b) = count(a) +
///     count(b), and likewise for the length, which pins GC%(a+b) exactly.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class CompositionMetamorphicTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed so shuffles/random inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260619);

    private const double Tolerance = 1e-9;

    /// <summary>Generates a random DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Deterministic Fisher–Yates shuffle of a string (permutes bases, preserves multiset).</summary>
    private static string Shuffle(string s)
    {
        var chars = s.ToCharArray();
        for (int i = chars.Length - 1; i > 0; i--)
        {
            int j = Rng.Next(i + 1);
            (chars[i], chars[j]) = (chars[j], chars[i]);
        }
        return new string(chars);
    }

    /// <summary>Representative DNA sequences plus a few fixed-seed random ones.</summary>
    private static string[] SampleSequences() => new[]
    {
        "ACGT",
        "GGGGCCCC",
        "AAAATTTT",
        "ATGCATGCATGC",
        "GATTACA",
        "CGCGCGCGAT",
        "TTTTTTTTTT",
        RandomDna(37),
        RandomDna(64),
        RandomDna(101),
    };

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-001 — GC content
    // ═══════════════════════════════════════════════════════════════════

    #region MR1: INV — complement preserves GC%

    /// <summary>
    /// MR1: GC%(seq) == GC%(complement(seq)).
    /// Complement maps A↔T and C↔G; a G becomes a C and a C becomes a G, so the
    /// (G + C) count — and thus GC% — is invariant. Verified on fixed and random
    /// (fixed-seed) sequences.
    /// </summary>
    [Test]
    public void GcContent_Complement_PreservesGcPercent()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            var complement = seq.Complement();

            complement.GcContent().Should().BeApproximately(seq.GcContent(), Tolerance,
                because: $"complement maps C↔G (G+C count invariant), so GC% of '{s}' must equal GC% of its complement");
        }
    }

    #endregion

    #region MR2: INV — shuffle (permutation) preserves GC%

    /// <summary>
    /// MR2: GC%(seq) == GC%(shuffle(seq)).
    /// GC% is a function of base *counts* only, not of order, so any permutation of
    /// the bases leaves it unchanged. Uses a deterministic fixed-seed Fisher–Yates
    /// shuffle so the test is reproducible.
    /// </summary>
    [Test]
    public void GcContent_Shuffle_PreservesGcPercent()
    {
        foreach (var s in SampleSequences())
        {
            var original = new DnaSequence(s).GcContent();

            // A few independent shuffles to exercise different permutations.
            for (int t = 0; t < 5; t++)
            {
                var shuffled = new DnaSequence(Shuffle(s)).GcContent();
                shuffled.Should().BeApproximately(original, Tolerance,
                    because: $"GC% depends only on composition, not order, so shuffling '{s}' must not change it");
            }
        }
    }

    #endregion

    #region MR3: INV — case-insensitive

    /// <summary>
    /// MR3: GC%(seq) == GC%(seq.ToLower()) == GC%(seq.ToUpper()).
    /// Counting is case-folded. We exercise the underlying case-sensitive string API
    /// (CalculateGcContentFast) directly — DnaSequence upper-cases at construction,
    /// so testing only DnaSequence would not exercise the lower-case code path. The
    /// DnaSequence facade is also asserted to confirm the contract end-to-end.
    /// </summary>
    [Test]
    public void GcContent_CaseFolding_IsCaseInsensitive()
    {
        foreach (var s in SampleSequences())
        {
            string upper = s.ToUpperInvariant();
            string lower = s.ToLowerInvariant();

            double gcUpper = upper.CalculateGcContentFast();
            double gcLower = lower.CalculateGcContentFast();

            gcLower.Should().BeApproximately(gcUpper, Tolerance,
                because: $"GC counting is case-insensitive, so '{lower}' and '{upper}' must give the same GC%");

            // End-to-end via the DnaSequence facade (constructs from a lower-case spelling).
            new DnaSequence(lower).GcContent().Should().BeApproximately(new DnaSequence(upper).GcContent(), Tolerance,
                because: "DnaSequence.GcContent() must be case-insensitive regardless of input spelling");
        }
    }

    #endregion

    #region MR4: INV — reverse-complement preserves GC% (derived)

    /// <summary>
    /// MR4 (derived): GC%(seq) == GC%(reverseComplement(seq)).
    /// Reverse-complement = reverse (a permutation — GC%-invariant by MR2) composed
    /// with complement (GC%-invariant by MR1). The composition is therefore also
    /// GC%-invariant. This is theoretically guaranteed, not a weakened heuristic.
    /// </summary>
    [Test]
    public void GcContent_ReverseComplement_PreservesGcPercent()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            var revComp = seq.ReverseComplement();

            revComp.GcContent().Should().BeApproximately(seq.GcContent(), Tolerance,
                because: $"reverse-complement = permutation ∘ complement, both GC%-invariant, so GC% of '{s}' is preserved");
        }
    }

    #endregion

    #region MR5: ADD — GC count is additive over concatenation (derived)

    /// <summary>
    /// MR5 (derived): the GC *count* is additive over concatenation, so
    ///   GC%(a + b) × len(a + b) == GC%(a) × len(a) + GC%(b) × len(b).
    /// Because (G + C) and length are both additive when two sequences are joined,
    /// the GC percentage of the concatenation is pinned exactly by the two parts.
    /// We compare GC *counts* (percentage × length / 100) to avoid percentage
    /// renormalisation ambiguity. All sequences are non-empty over {A,C,G,T}.
    /// </summary>
    [Test]
    public void GcContent_Concatenation_GcCountIsAdditive()
    {
        var samples = SampleSequences();
        foreach (var a in samples)
        {
            foreach (var b in samples)
            {
                var seqA = new DnaSequence(a);
                var seqB = new DnaSequence(b);
                var seqAb = new DnaSequence(a + b);

                // GC count = GC% × length / 100. Lengths are exact (DnaSequence over A/C/G/T).
                double gcCountA = seqA.GcContent() * seqA.Length / 100.0;
                double gcCountB = seqB.GcContent() * seqB.Length / 100.0;
                double gcCountAb = seqAb.GcContent() * seqAb.Length / 100.0;

                gcCountAb.Should().BeApproximately(gcCountA + gcCountB, 1e-6,
                    because: $"(G+C) count is additive over concatenation: count('{a}'+'{b}') = count('{a}') + count('{b}')");
            }
        }
    }

    #endregion
}
