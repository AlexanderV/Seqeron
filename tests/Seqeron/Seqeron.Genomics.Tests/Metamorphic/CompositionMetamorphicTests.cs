using System;
using System.Linq;
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
/// Units: SEQ-GC-001 — GC content (Composition); SEQ-COMP-001 — DNA complement (Composition);
///        SEQ-REVCOMP-001 — reverse complement (Composition)
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, rows 1–3.
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

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMP-001 — DNA complement
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Sequence_Composition/Sequence_Composition.md;
    //   SequenceExtensions.GetComplementBase — "Source: Wikipedia Nucleic acid
    //   notation — IUPAC complement table"):
    //   The complement operation maps each base to its Watson–Crick partner —
    //   A↔T, C↔G for canonical DNA — and, for ambiguity codes, to the complement
    //   of the *set* it denotes: R↔Y, K↔M, B↔V, D↔H, with the self-complementary
    //   codes S↔S, W↔W, N↔N. Every entry of this table is paired and self-inverse,
    //   so complement is an INVOLUTION (complement∘complement = identity), and it
    //   rewrites positions one-for-one, so it is LENGTH-PRESERVING. These relations
    //   follow from the table alone, independent of any particular input.
    //
    // API surface under test:
    //   • DnaSequence.Complement() — canonical {A,C,G,T} path. DnaSequence validates
    //     strictly to A/C/G/T at construction, so it cannot carry IUPAC ambiguity
    //     codes; the IUPAC alphabet is therefore exercised through the lower-level
    //     SequenceExtensions.GetComplementBase / TryGetComplement (the same table
    //     DnaSequence.Complement() is built on).

    #region MR6: INV — complement is an involution over {A,C,G,T} (DnaSequence)

    /// <summary>
    /// MR6: complement(complement(x)) == x for canonical DNA via the DnaSequence facade.
    /// The complement table pairs A↔T and C↔G; applying it twice returns each base to
    /// itself, so the double complement reproduces the original sequence exactly.
    /// Verified on fixed and fixed-seed random sequences.
    /// </summary>
    [Test]
    public void Complement_AppliedTwice_IsIdentity_Dna()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);

            var doubleComplement = seq.Complement().Complement();

            doubleComplement.Sequence.Should().Be(seq.Sequence,
                because: $"complement pairs A↔T and C↔G (each mapping self-inverse), so complement∘complement must return '{s}' unchanged");
        }
    }

    #endregion

    #region MR7: INV — complement is an involution over the full IUPAC alphabet (GetComplementBase)

    /// <summary>
    /// MR7: GetComplementBase(GetComplementBase(c)) == c for every IUPAC code.
    /// The IUPAC complement table is a fixed-point-free / self-paired permutation:
    /// A↔T, C↔G, U→A→T (note: U complements to A, A to T — so U is NOT a fixed point
    /// of the double map and is excluded below), R↔Y, K↔M, B↔V, D↔H, and the
    /// self-complementary S↔S, W↔W, N↔N. Over the DNA-emitting alphabet (no U input)
    /// applying the map twice is therefore the identity.
    /// </summary>
    [Test]
    public void Complement_AppliedTwice_IsIdentity_AllIupacCodes()
    {
        // DNA IUPAC alphabet (uppercase). U is excluded: GetComplementBase('U')='A'
        // and 'A'→'T', so U is intentionally not an involution fixed point in the
        // DNA-emitting table; the involution is over the {A,C,G,T}+ambiguity set.
        const string iupac = "ACGTRYSWKMBDHVN";

        foreach (char c in iupac)
        {
            char twice = SequenceExtensions.GetComplementBase(SequenceExtensions.GetComplementBase(c));

            twice.Should().Be(c,
                because: $"the IUPAC complement table is self-paired, so complementing '{c}' twice must yield '{c}'");
        }
    }

    /// <summary>
    /// MR7-b: case-insensitive involution — lower-case IUPAC input complemented twice
    /// returns its upper-case canonical form (the table upper-cases recognised bases).
    /// Confirms the involution holds regardless of input spelling.
    /// </summary>
    [Test]
    public void Complement_AppliedTwice_IsIdentity_LowerCaseIupac()
    {
        const string iupac = "acgtryswkmbdhvn";

        foreach (char c in iupac)
        {
            char twice = SequenceExtensions.GetComplementBase(SequenceExtensions.GetComplementBase(c));

            twice.Should().Be(char.ToUpperInvariant(c),
                because: $"complement is case-folding and self-paired, so '{c}' complemented twice must yield '{char.ToUpperInvariant(c)}'");
        }
    }

    /// <summary>
    /// MR7-c: whole-string IUPAC involution via the span-based TryGetComplement.
    /// Complementing an IUPAC string twice reproduces it (in upper case). This
    /// exercises the same table over the path DnaSequence cannot reach (it rejects
    /// ambiguity codes at construction).
    /// </summary>
    [Test]
    public void Complement_AppliedTwice_IsIdentity_IupacString()
    {
        var rng = new Random(20260619);
        const string iupac = "ACGTRYSWKMBDHVN";

        for (int t = 0; t < 20; t++)
        {
            int len = 1 + rng.Next(40);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = iupac[rng.Next(iupac.Length)];
            string original = new string(chars);

            var once = new char[len];
            original.AsSpan().TryGetComplement(once).Should().BeTrue();
            var twice = new char[len];
            ((ReadOnlySpan<char>)once).TryGetComplement(twice).Should().BeTrue();

            new string(twice).Should().Be(original,
                because: $"the IUPAC complement table is self-paired, so complementing '{original}' twice must reproduce it");
        }
    }

    #endregion

    #region MR8: INV — complement preserves length

    /// <summary>
    /// MR8: complement(x).Length == x.Length.
    /// Complement rewrites each position to exactly one output character — it never
    /// inserts or deletes — so the output length always equals the input length.
    /// </summary>
    [Test]
    public void Complement_PreservesLength()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);

            seq.Complement().Length.Should().Be(seq.Length,
                because: $"complement is a per-position rewrite (no insert/delete), so length of '{s}' is unchanged");
        }
    }

    #endregion

    #region MR9: P — complement swaps complementary base counts (derived from theory)

    /// <summary>
    /// MR9 (derived): complement maps A↔T and C↔G, so the count of A in x equals the
    /// count of T in complement(x), and count of C in x equals count of G in
    /// complement(x) (and symmetrically). This is a direct, theory-guaranteed
    /// consequence of the canonical complement table over {A,C,G,T}.
    /// </summary>
    [Test]
    public void Complement_SwapsComplementaryBaseCounts()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            string comp = seq.Complement().Sequence;

            int aIn = seq.Sequence.Count(ch => ch == 'A');
            int tIn = seq.Sequence.Count(ch => ch == 'T');
            int cIn = seq.Sequence.Count(ch => ch == 'C');
            int gIn = seq.Sequence.Count(ch => ch == 'G');

            comp.Count(ch => ch == 'T').Should().Be(aIn,
                because: $"A↔T: every A in '{s}' becomes a T in its complement");
            comp.Count(ch => ch == 'A').Should().Be(tIn,
                because: $"A↔T: every T in '{s}' becomes an A in its complement");
            comp.Count(ch => ch == 'G').Should().Be(cIn,
                because: $"C↔G: every C in '{s}' becomes a G in its complement");
            comp.Count(ch => ch == 'C').Should().Be(gIn,
                because: $"C↔G: every G in '{s}' becomes a C in its complement");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-REVCOMP-001 — reverse complement
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Sequence_Composition/Sequence_Composition.md;
    //   Watson–Crick antiparallel base pairing):
    //   The reverse complement is the sequence read 3'→5' on the opposing strand,
    //   i.e. reverse-complement = reverse ∘ complement = complement ∘ reverse
    //   (the two operations commute — reversing positions and rewriting each base
    //   are independent). Because complement is an INVOLUTION (A↔T, C↔G, each
    //   self-inverse — SEQ-COMP-001) and reverse is its own inverse, their
    //   composition is also an INVOLUTION: revcomp(revcomp(x)) = x. Each step is a
    //   per-position rewrite or reordering that never inserts/deletes, so revcomp is
    //   LENGTH-PRESERVING. These relations follow from the definition alone,
    //   independent of any particular input. We do NOT re-test the plain-complement
    //   involution here (that is SEQ-COMP-001, above).
    //
    // API surface under test:
    //   • DnaSequence.ReverseComplement() — canonical {A,C,G,T} facade.
    //   • DnaSequence.GetReverseComplementString(string) — static string helper.

    #region MR10: INV — reverse-complement is an involution over {A,C,G,T}

    /// <summary>
    /// MR10: revcomp(revcomp(x)) == x.
    /// Reverse-complement = reverse ∘ complement, and both reverse and complement are
    /// self-inverse, so applying the composition twice returns the original sequence
    /// exactly. Verified on fixed and fixed-seed random sequences via the DnaSequence
    /// facade and the static string helper.
    /// </summary>
    [Test]
    public void ReverseComplement_AppliedTwice_IsIdentity_Dna()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);

            var doubleRevComp = seq.ReverseComplement().ReverseComplement();
            doubleRevComp.Sequence.Should().Be(seq.Sequence,
                because: $"revcomp = reverse ∘ complement, both self-inverse, so revcomp∘revcomp must return '{s}' unchanged");

            // Same involution via the static string helper.
            string twiceString = DnaSequence.GetReverseComplementString(
                DnaSequence.GetReverseComplementString(s));
            twiceString.Should().Be(s.ToUpperInvariant(),
                because: $"GetReverseComplementString is the same involution, so applying it twice must reproduce '{s}'");
        }
    }

    #endregion

    #region MR11: INV — reverse-complement preserves length

    /// <summary>
    /// MR11: revcomp(x).Length == x.Length.
    /// Reverse reorders positions and complement rewrites each base one-for-one;
    /// neither inserts nor deletes, so the output length always equals the input
    /// length.
    /// </summary>
    [Test]
    public void ReverseComplement_PreservesLength()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);

            seq.ReverseComplement().Length.Should().Be(seq.Length,
                because: $"revcomp is a reorder + per-position rewrite (no insert/delete), so length of '{s}' is unchanged");
        }
    }

    #endregion

    #region MR12: COMP — revcomp = reverse ∘ complement = complement ∘ reverse (derived)

    /// <summary>
    /// MR12 (derived): revcomp(x) == reverse(complement(x)) == complement(reverse(x)).
    /// By definition the reverse complement is the complement read in reverse, and
    /// reversing positions commutes with the per-position complement rewrite, so both
    /// orders of composition agree with ReverseComplement(). This pins the exact
    /// output, not just an invariant.
    /// </summary>
    [Test]
    public void ReverseComplement_EqualsReverseComposedWithComplement()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            string revComp = seq.ReverseComplement().Sequence;

            // reverse(complement(x))
            string reverseOfComplement = new string(seq.Complement().Sequence.Reverse().ToArray());
            // complement(reverse(x))  — build a reversed DnaSequence, then complement it
            var reversed = new DnaSequence(new string(seq.Sequence.Reverse().ToArray()));
            string complementOfReverse = reversed.Complement().Sequence;

            revComp.Should().Be(reverseOfComplement,
                because: $"revcomp is the complement read in reverse, so it must equal reverse(complement('{s}'))");
            revComp.Should().Be(complementOfReverse,
                because: $"reverse and complement commute, so revcomp('{s}') must also equal complement(reverse('{s}'))");
        }
    }

    #endregion

    #region MR13: INV — reverse-complement preserves GC% (derived)

    /// <summary>
    /// MR13 (derived): GC%(revcomp(x)) == GC%(x).
    /// Reverse is a permutation (GC%-invariant) and complement maps C↔G leaving the
    /// (G + C) count unchanged (GC%-invariant), so their composition preserves GC%.
    /// This restates the revcomp = permutation ∘ complement decomposition from a
    /// composition angle (distinct from SEQ-COMP-001's plain-complement relation).
    /// </summary>
    [Test]
    public void ReverseComplement_PreservesGcPercent()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);

            seq.ReverseComplement().GcContent().Should().BeApproximately(seq.GcContent(), Tolerance,
                because: $"revcomp = permutation ∘ complement, both GC%-invariant, so GC% of '{s}' is preserved");
        }
    }

    #endregion
}
