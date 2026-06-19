using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

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
///        SEQ-REVCOMP-001 — reverse complement (Composition);
///        SEQ-VALID-001 — sequence validation (Composition);
///        SEQ-COMPLEX-001 — sequence complexity measure (Composition);
///        SEQ-ENTROPY-001 — Shannon entropy of base composition (Composition);
///        SEQ-GCSKEW-001 — GC skew (Composition)
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, rows 1–7.
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

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-VALID-001 — sequence validation
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Sequence_Composition; IUPAC-IUB 1970 / NC-IUB 1984
    //   nucleotide notation):
    //   • A *valid DNA* string contains only the four canonical bases A/C/G/T,
    //     case-insensitively. Validity is a property of the SET of characters
    //     present — neither order nor multiplicity matter.
    //   • The *IUPAC nucleotide alphabet* is a strict SUPERSET of the DNA alphabet:
    //     it adds the ambiguity codes R,Y,S,W,K,M,B,D,H,V,N (and, in this codebase,
    //     U and the gap symbols '-','.'). Hence every DNA-valid string is also
    //     IUPAC-valid, but not conversely — there exist IUPAC-valid strings
    //     (e.g. containing R/Y/N) that are NOT DNA-valid. This makes the
    //     subset relation a PROPER subset, so the COMP relation is non-vacuous.
    //
    // Three metamorphic relations (checklist row 4):
    //   INV  case conversion preserves validity:  IsValid(x)=IsValid(lower)=IsValid(upper)
    //   COMP valid DNA ⊂ valid IUPAC:             DNA-valid ⇒ IUPAC-valid (proper subset)
    //   INV  repeat seq → same result:            IsValid(x)=IsValid(x repeated k)
    //
    // API surface under test (confirmed in src/.../Seqeron.Genomics.Core):
    //   • DNA validity:  SequenceExtensions.IsValidDna(this ReadOnlySpan<char>)
    //                    — case-insensitive (char.ToUpperInvariant per char), accepts
    //                      ONLY A/C/G/T (SequenceExtensions.cs).
    //   • IUPAC validity: IupacDnaSequence.IsValid() (inherited SequenceBase.IsValid,
    //                    ISequence.cs) — tests each char against IupacDnaSequence.Alphabet,
    //                    the IUPAC superset {A,C,G,T,U,N,R,Y,W,S,K,M,B,D,H,V,-,.}.
    //   ADAPTATION (documented): there is no separate case-insensitive IsValidIupac
    //   helper, and IupacDnaSequence.Alphabet is UPPERCASE-only, so IupacDnaSequence
    //   .IsValid() is case-SENSITIVE. The COMP subset relation is therefore evaluated
    //   on canonical UPPERCASE spellings (the relation "valid DNA ⊂ valid IUPAC" is a
    //   statement over the alphabets, independent of casing). The case-INVARIANCE
    //   relation is tested separately against the case-insensitive IsValidDna API.

    #region SEQ-VALID-001 — sequence validation

    #region MR14: INV — case conversion preserves DNA validity

    /// <summary>
    /// MR14: IsValidDna(x) == IsValidDna(x.ToLower()) == IsValidDna(x.ToUpper()).
    /// DNA validity is case-folded (IsValidDna upper-cases each char before checking),
    /// so the lower-, upper-, and mixed-case spellings of any string agree on validity.
    /// Exercised on DNA-valid samples AND on strings carrying out-of-alphabet
    /// characters, so the invariance is shown for BOTH the true and the false verdict.
    /// </summary>
    [Test]
    public void IsValidDna_CaseConversion_PreservesValidity()
    {
        var inputs = SampleSequences()
            .Concat(new[] { "acgt", "AcGtAcGt", "GATTACA", "XYZ", "ACGTN", "ACGU", "ACG-T", "" })
            .ToArray();

        foreach (var s in inputs)
        {
            bool asIs = s.AsSpan().IsValidDna();
            bool lower = s.ToLowerInvariant().AsSpan().IsValidDna();
            bool upper = s.ToUpperInvariant().AsSpan().IsValidDna();

            lower.Should().Be(asIs,
                because: $"IsValidDna is case-insensitive, so lower-casing '{s}' must not change its validity verdict");
            upper.Should().Be(asIs,
                because: $"IsValidDna is case-insensitive, so upper-casing '{s}' must not change its validity verdict");
        }
    }

    #endregion

    #region MR15: COMP — valid DNA ⊂ valid IUPAC (proper subset)

    /// <summary>
    /// MR15: every DNA-valid string is also IUPAC-valid (subset direction).
    /// The IUPAC alphabet is a superset of {A,C,G,T}, so a string that validates as
    /// DNA must validate as IUPAC. Verified across many fixed and fixed-seed-random
    /// DNA strings. Validity is evaluated on canonical UPPERCASE spellings because
    /// IupacDnaSequence.IsValid() is case-sensitive (see file note); the subset
    /// relation is a statement over alphabets and holds regardless of casing.
    /// </summary>
    [Test]
    public void Validity_EveryDnaValidString_IsAlsoIupacValid_Subset()
    {
        var dnaStrings = SampleSequences()
            .Concat(Enumerable.Range(0, 30).Select(_ => RandomDna(1 + Rng.Next(60))))
            .ToArray();

        foreach (var s in dnaStrings)
        {
            string canonical = s.ToUpperInvariant();

            // Precondition: this sample really is DNA-valid (otherwise the implication is vacuous here).
            canonical.AsSpan().IsValidDna().Should().BeTrue(
                because: $"'{s}' is drawn from the {{A,C,G,T}} alphabet, so it must be DNA-valid");

            new IupacDnaSequence(canonical).IsValid().Should().BeTrue(
                because: $"the IUPAC alphabet is a superset of {{A,C,G,T}}, so DNA-valid '{s}' must also be IUPAC-valid");
        }
    }

    /// <summary>
    /// MR15-b: the subset is PROPER — there EXIST IUPAC-valid strings that are NOT
    /// DNA-valid (they contain ambiguity codes R/Y/S/W/K/M/B/D/H/V/N absent from the
    /// DNA alphabet). This confirms the COMP relation is non-vacuous: valid DNA is a
    /// strict, not an equal, subset of valid IUPAC.
    /// </summary>
    [Test]
    public void Validity_IupacAmbiguityCodes_AreIupacValidButNotDnaValid_ProperSubset()
    {
        // Each string is built from the IUPAC superset and contains ≥1 ambiguity code,
        // so it lies in (valid IUPAC) \ (valid DNA).
        string[] iupacOnly =
        {
            "ACGTN", "RYSWKM", "BDHV", "ACGTRYN", "NNNN", "ACGTRYSWKMBDHVN",
        };

        foreach (var s in iupacOnly)
        {
            new IupacDnaSequence(s).IsValid().Should().BeTrue(
                because: $"'{s}' is composed entirely of IUPAC codes, so it must be IUPAC-valid");

            s.AsSpan().IsValidDna().Should().BeFalse(
                because: $"'{s}' contains an IUPAC ambiguity code outside {{A,C,G,T}}, so it must NOT be DNA-valid — proving the subset is proper");
        }
    }

    #endregion

    #region MR16: INV — repeating the sequence preserves validity

    /// <summary>
    /// MR16: IsValidDna(x) == IsValidDna(x repeated k times).
    /// Validity depends only on the SET of characters present, not on length or order,
    /// so concatenating a sequence with itself any number of times cannot change the
    /// verdict: a valid sequence stays valid (its character set is unchanged) and an
    /// invalid one stays invalid (the offending character recurs in every copy).
    /// Verified for valid and invalid inputs across several repeat counts.
    /// </summary>
    [Test]
    public void IsValidDna_RepeatedSequence_PreservesValidity()
    {
        var inputs = SampleSequences()
            .Concat(new[] { "GATTACA", "XYZ", "ACGTN", "ACG-T" })
            .ToArray();

        foreach (var s in inputs)
        {
            bool original = s.AsSpan().IsValidDna();

            foreach (int k in new[] { 2, 3, 5 })
            {
                string repeated = string.Concat(Enumerable.Repeat(s, k));
                repeated.AsSpan().IsValidDna().Should().Be(original,
                    because: $"validity is a property of the character set only, so repeating '{s}' {k}× must not change its verdict");
            }
        }
    }

    /// <summary>
    /// MR16-b: concatenation/duplication corollaries.
    /// Joining two DNA-valid sequences stays valid (the union of two subsets of
    /// {A,C,G,T} is still a subset), while appending a single out-of-alphabet
    /// character to ANY sequence makes it invalid (the character set now escapes the
    /// DNA alphabet). Both follow directly from validity being a set-membership test.
    /// </summary>
    [Test]
    public void IsValidDna_Concatenation_FollowsCharacterSetSemantics()
    {
        var dnaSamples = SampleSequences();

        foreach (var a in dnaSamples)
        {
            foreach (var b in dnaSamples)
            {
                (a + b).AsSpan().IsValidDna().Should().BeTrue(
                    because: $"concatenating DNA-valid '{a}' and '{b}' keeps the character set within {{A,C,G,T}}, so the result stays DNA-valid");
            }

            // Appending an out-of-alphabet character escapes the DNA alphabet.
            (a + "N").AsSpan().IsValidDna().Should().BeFalse(
                because: $"appending the ambiguity code 'N' to '{a}' introduces a character outside {{A,C,G,T}}, so the result must be invalid");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMPLEX-001 — sequence complexity measure
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Statistics/Entropy_Profile.md §"Description"; Shannon 1948;
    //   tests/TestSpecs/SEQ-COMPLEX-001.md §3.2):
    //   The checklist row-5 relations — INV permutation preserves complexity; MON homopolymer
    //   → MIN complexity; MON random → higher — only hold for a COMPOSITION-based measure, i.e.
    //   one that is a function of the per-symbol (mononucleotide) frequency distribution alone.
    //   Of the measures exposed by SequenceComplexity, the order-dependent ones — linguistic
    //   complexity (counts distinct subwords of length ≥ 2), DUST, k-mer entropy and Lempel–Ziv
    //   — are NOT permutation-invariant (a shuffle changes their subword/triplet/parse content),
    //   and those order-dependent measures are the SEPARATE units SEQ-COMPLEX-DUST/KMER/
    //   COMPRESS/WINDOW-001 (checklist rows 228–231). The composition-based complexity measure
    //   the class exposes is the SHANNON ENTROPY of the mononucleotide distribution,
    //   SequenceComplexity.CalculateShannonEntropy:
    //
    //       H = −Σ_{b∈{A,C,G,T}} p_b · log₂ p_b ,   p_b = count(b) / N
    //
    //   It depends only on the four base COUNTS, so it is permutation-invariant; a homopolymer
    //   (one base, p=1) gives H = 0 — the documented MINIMUM — and any sequence with ≥ 2
    //   distinct bases has H > 0; a uniform 4-base sequence gives the 2-bit maximum (log₂4).
    //   These are exactly the three row-5 relations, and they are theory-guaranteed, not weakened.
    //
    //   RECONCILIATION: row 5's label is the generic "sequence complexity measure". The only
    //   API in SequenceComplexity for which all three stated relations genuinely hold is the
    //   composition-based Shannon entropy; the order-dependent complexity APIs cannot satisfy
    //   "permutation preserves complexity" and are covered by their own units. We therefore
    //   encode SEQ-COMPLEX-001 against CalculateShannonEntropy (the documented minimum H = 0 is
    //   used as the exact endpoint; monotonicity uses orderings, never magic absolute values).
    //
    // API surface under test:
    //   • SequenceComplexity.CalculateShannonEntropy(string)  — composition-based complexity.
    //   • SequenceComplexity.CalculateShannonEntropy(DnaSequence) — facade (same core).

    #region SEQ-COMPLEX-001 — sequence complexity

    #region MR17: INV — permutation preserves complexity

    /// <summary>
    /// MR17: complexity(x) == complexity(permute(x)).
    /// The composition-based complexity (Shannon entropy of the mononucleotide distribution)
    /// is a function of base COUNTS only, so any permutation of the bases — here a fixed-seed
    /// Fisher–Yates shuffle — leaves it unchanged. Verified on fixed and fixed-seed random
    /// sequences, with several independent shuffles per input.
    /// </summary>
    [Test]
    public void Complexity_Permutation_PreservesComplexity()
    {
        foreach (var s in SampleSequences())
        {
            double original = SequenceComplexity.CalculateShannonEntropy(s);

            for (int t = 0; t < 5; t++)
            {
                double shuffled = SequenceComplexity.CalculateShannonEntropy(Shuffle(s));
                shuffled.Should().BeApproximately(original, Tolerance,
                    because: $"composition-based complexity depends only on base counts, not order, so permuting '{s}' must not change it");
            }
        }
    }

    #endregion

    #region MR18: MON — homopolymer → minimum complexity

    /// <summary>
    /// MR18: a single-symbol run yields the documented MINIMUM complexity (Shannon entropy
    /// H = 0 exactly) and is ≤ the complexity of any other sequence; any sequence with ≥ 2
    /// distinct bases has strictly greater complexity (H > 0). The minimum is pinned exactly
    /// (theory endpoint); the comparison against other sequences is an ordering, not a magic
    /// absolute value.
    /// </summary>
    [Test]
    public void Complexity_Homopolymer_IsMinimumAndStrictlyBelowDiverse()
    {
        // Homopolymers over each base reach the exact documented minimum, H = 0.
        foreach (char b in "ACGT")
        {
            foreach (int len in new[] { 1, 4, 16, 64 })
            {
                string homo = new string(b, len);
                SequenceComplexity.CalculateShannonEntropy(homo).Should().BeApproximately(0.0, Tolerance,
                    because: $"a homopolymer '{homo}' has a single symbol (p=1), so its Shannon entropy is exactly the documented minimum 0");
            }
        }

        // Homopolymer ≤ any other sequence, with strict inequality once ≥ 2 distinct bases appear.
        string homopolymer = new string('A', 32);
        double homoComplexity = SequenceComplexity.CalculateShannonEntropy(homopolymer);

        foreach (var s in SampleSequences())
        {
            double sComplexity = SequenceComplexity.CalculateShannonEntropy(s);

            sComplexity.Should().BeGreaterThanOrEqualTo(homoComplexity - Tolerance,
                because: $"the homopolymer attains the minimum complexity, so no sequence ('{s}') can fall below it");

            if (s.Distinct().Count(c => "ACGT".Contains(c)) >= 2)
            {
                sComplexity.Should().BeGreaterThan(homoComplexity + Tolerance,
                    because: $"'{s}' has ≥ 2 distinct bases, so its composition-based complexity must strictly exceed the homopolymer minimum");
            }
        }
    }

    #endregion

    #region MR19: MON — random (high-diversity) → higher complexity than low-diversity

    /// <summary>
    /// MR19: a high-diversity sequence (fixed-seed random over all four bases, roughly equal)
    /// has higher complexity than a low-diversity one (homopolymer or near-homopolymer), and a
    /// perfectly uniform 4-base sequence attains the documented MAXIMUM, H = log₂4 = 2 bits.
    /// The ordering is asserted (low &lt; high), and only the theory-pinned uniform maximum uses
    /// an exact endpoint — no magic intermediate constants.
    /// </summary>
    [Test]
    public void Complexity_Random_HigherThanLowDiversity()
    {
        const int length = 200;

        // Low-diversity references.
        double homopolymer = SequenceComplexity.CalculateShannonEntropy(new string('G', length));
        // Near-homopolymer: overwhelmingly one base with a single foreign base.
        double nearHomopolymer = SequenceComplexity.CalculateShannonEntropy(new string('A', length - 1) + "C");

        // High-diversity reference: fixed-seed random over {A,C,G,T}.
        double random = SequenceComplexity.CalculateShannonEntropy(RandomDna(length));

        random.Should().BeGreaterThan(homopolymer + Tolerance,
            because: "a random sequence using all four bases is compositionally more diverse than a homopolymer, so it has higher complexity");
        random.Should().BeGreaterThan(nearHomopolymer + Tolerance,
            because: "a balanced random sequence has a flatter base distribution than a near-homopolymer, so its Shannon complexity is higher");
        nearHomopolymer.Should().BeGreaterThan(homopolymer + Tolerance,
            because: "introducing a second base raises the distribution above the single-symbol minimum, so a near-homopolymer exceeds a homopolymer");

        // Theory endpoint: a perfectly uniform 4-base sequence attains the documented maximum, 2 bits.
        string uniform = string.Concat(Enumerable.Repeat("ACGT", length / 4));
        SequenceComplexity.CalculateShannonEntropy(uniform).Should().BeApproximately(2.0, Tolerance,
            because: "an equal A/C/G/T composition is the maximum-entropy distribution for DNA, so H = log₂4 = 2 bits");

        // And the random sample never exceeds that documented maximum.
        random.Should().BeLessThanOrEqualTo(2.0 + Tolerance,
            because: "Shannon entropy over four symbols is bounded above by log₂4 = 2 bits");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-ENTROPY-001 — Shannon entropy of base composition
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Sequence_Composition/Shannon_Entropy.md; Shannon 1948;
    //   Cover & Thomas 1991):
    //       H(X) = −Σ_{b} p_b · log₂ p_b ,   p_b = count(b) / N
    //   over the canonical DNA base distribution {A,C,G,T} (the implementation counts
    //   ONLY A/T/G/C; base-2 log ⇒ bits; empty/null ⇒ 0 — §2.4, §6.1). Entropy is a
    //   functional of the FREQUENCY DISTRIBUTION alone. Three entropy-specific theorems
    //   pin this unit (distinct from the SEQ-COMPLEX-001 "complexity-measure" framing of
    //   the same method above):
    //
    //   • INV (permutation): H depends only on the multiset of base counts {n_A,n_C,n_G,
    //     n_T}, not on their arrangement, so H(x) = H(any permutation of x). This is the
    //     information-theoretic invariance of entropy under relabelling of positions.
    //
    //   • MON (uniform → MAX): for a distribution supported on k symbols, H ≤ log₂k with
    //     EQUALITY iff the distribution is uniform (maximum-entropy characterisation; a
    //     corollary of Gibbs' inequality −Σ p log₂ p ≤ −Σ p log₂ q with q ≡ 1/k). So an
    //     equiprobable k-symbol sequence attains the EXACT endpoint H = log₂k (k=2 → 1
    //     bit, k=3 → log₂3, k=4 → 2 bits), and ANY non-uniform distribution on the same
    //     support has STRICTLY smaller H. We pin the endpoint exactly at several k and
    //     assert the strict Gibbs ordering for non-uniform composition.
    //
    //   • MON (single symbol → 0): a one-symbol distribution has p=1 (and 0·log0≡0 for
    //     the rest), so H = 0 EXACTLY — the global minimum (H ≥ 0 always). Introducing a
    //     second distinct base lifts the support to k=2 and forces H > 0.
    //
    //   DIFFERENTIATION FROM SEQ-COMPLEX-001: that unit frames CalculateShannonEntropy as
    //   a generic complexity measure (permutation-invariance, homopolymer = minimum, and a
    //   random > low-diversity ordering). THIS unit encodes the entropy THEORY proper: the
    //   maximum-entropy law H_max = log₂(k) verified at multiple alphabet sizes k∈{2,3,4},
    //   the strict Gibbs inequality (uniform > any non-uniform of equal support), and the
    //   exact-0 single-symbol endpoint with the "second symbol ⇒ H>0" lift. The test
    //   bodies, sequences and assertions are deliberately non-duplicative.
    //
    // API surface under test:
    //   • SequenceComplexity.CalculateShannonEntropy(string)      — canonical A/T/G/C path.
    //   • SequenceComplexity.CalculateShannonEntropy(DnaSequence) — facade (same core).

    #region SEQ-ENTROPY-001 — Shannon entropy

    #region MR20: INV — permutation preserves Shannon entropy

    /// <summary>
    /// MR20: H(x) == H(permute(x)).
    /// Shannon entropy is a functional of the base-frequency distribution {n_A,n_C,n_G,n_T}
    /// alone; relabelling positions (a fixed-seed Fisher–Yates shuffle) leaves the multiset of
    /// counts — and therefore H — invariant. Asserted via both the string and DnaSequence APIs,
    /// with several independent shuffles per input.
    /// </summary>
    [Test]
    public void ShannonEntropy_Permutation_PreservesEntropy()
    {
        foreach (var s in SampleSequences())
        {
            double original = SequenceComplexity.CalculateShannonEntropy(s);

            for (int t = 0; t < 5; t++)
            {
                string permuted = Shuffle(s);

                SequenceComplexity.CalculateShannonEntropy(permuted).Should().BeApproximately(original, Tolerance,
                    because: $"entropy is a functional of the base-frequency distribution only, so permuting '{s}' leaves H unchanged");

                // Same invariance through the DnaSequence facade (constructs over A/C/G/T).
                SequenceComplexity.CalculateShannonEntropy(new DnaSequence(permuted)).Should()
                    .BeApproximately(original, Tolerance,
                        because: $"the DnaSequence entropy facade is the same frequency functional, so permuting '{s}' must not change H");
            }
        }
    }

    #endregion

    #region MR21: MON — uniform → maximum entropy (H = log₂ k), strictly above any non-uniform

    /// <summary>
    /// MR21-a: a k-symbol EQUIPROBABLE sequence attains the EXACT maximum H = log₂(k).
    /// Verified at several alphabet sizes — k=2 (→1 bit), k=3 (→log₂3≈1.585 bits) and k=4
    /// (→2 bits) — by repeating a balanced unit so every counted base occurs equally often.
    /// The maximum-entropy law H_max = log₂(k) is pinned exactly (theory endpoint), not by a
    /// magic constant.
    /// </summary>
    [Test]
    public void ShannonEntropy_UniformDistribution_EqualsLog2OfAlphabetSize()
    {
        // (unit, k) pairs: each unit uses k distinct bases exactly once → equiprobable when repeated.
        var cases = new[]
        {
            (unit: "AT",   k: 2),
            (unit: "GC",   k: 2),
            (unit: "ACG",  k: 3),
            (unit: "ACGT", k: 4),
        };

        foreach (var (unit, k) in cases)
        {
            double expectedMax = Math.Log2(k);

            // Repeat the balanced unit so counts stay perfectly equal across many copies.
            foreach (int reps in new[] { 1, 4, 25 })
            {
                string uniform = string.Concat(Enumerable.Repeat(unit, reps));

                SequenceComplexity.CalculateShannonEntropy(uniform).Should().BeApproximately(expectedMax, Tolerance,
                    because: $"an equiprobable {k}-symbol sequence ('{unit}'×{reps}) is the maximum-entropy distribution, so H = log₂{k} exactly");
            }
        }
    }

    /// <summary>
    /// MR21-b: Gibbs' inequality (strict) — for a fixed support of k symbols, the UNIFORM
    /// distribution strictly maximises entropy, so any NON-uniform composition on the same
    /// k bases has H strictly below log₂(k). Encoded by skewing a balanced sequence (adding
    /// extra copies of one base) and asserting H(non-uniform) &lt; H(uniform) = log₂(k).
    /// </summary>
    [Test]
    public void ShannonEntropy_NonUniform_StrictlyBelowUniformMaximum()
    {
        // Support of all four bases, but skewed so the distribution is non-uniform.
        var skewedSequences = new[]
        {
            "AAAACGT",            // A over-represented
            "ACGTGGGG",          // G over-represented
            "AACCGTTTTTTTT",      // T over-represented
            string.Concat(Enumerable.Repeat("ACGT", 10)) + "AAAA", // near-uniform but tilted toward A
        };

        foreach (var s in skewedSequences)
        {
            int support = s.Distinct().Count(c => "ACGT".Contains(c));
            double uniformMax = Math.Log2(support);
            double h = SequenceComplexity.CalculateShannonEntropy(s);

            h.Should().BeLessThan(uniformMax - Tolerance,
                because: $"'{s}' uses {support} bases but non-uniformly, so by Gibbs' inequality H < log₂{support} (the uniform maximum)");
            h.Should().BeGreaterThan(0.0 + Tolerance,
                because: $"'{s}' has ≥ 2 distinct bases, so its entropy is strictly positive");
        }

        // Direct uniform-vs-skewed comparison on the SAME 4-base support: uniform > skewed.
        double hUniform = SequenceComplexity.CalculateShannonEntropy(string.Concat(Enumerable.Repeat("ACGT", 8)));
        double hSkewed = SequenceComplexity.CalculateShannonEntropy(string.Concat(Enumerable.Repeat("ACGT", 8)) + new string('A', 16));

        hUniform.Should().BeGreaterThan(hSkewed + Tolerance,
            because: "with the same {A,C,G,T} support, the uniform distribution strictly maximises entropy, so skewing toward A lowers H");
    }

    #endregion

    #region MR22: MON — single symbol → 0 exactly; second distinct symbol lifts H above 0

    /// <summary>
    /// MR22-a: a single-symbol sequence has H = 0 EXACTLY — the global minimum.
    /// One base carries probability 1 (and 0·log₂0 ≡ 0 for the rest), so there is no
    /// uncertainty and entropy collapses to 0. Verified for every base and several lengths,
    /// through both the string and DnaSequence APIs.
    /// </summary>
    [Test]
    public void ShannonEntropy_SingleSymbol_IsExactlyZero()
    {
        foreach (char b in "ACGT")
        {
            foreach (int len in new[] { 1, 2, 8, 50, 256 })
            {
                string homopolymer = new string(b, len);

                SequenceComplexity.CalculateShannonEntropy(homopolymer).Should().BeApproximately(0.0, Tolerance,
                    because: $"the single-symbol distribution '{homopolymer}' has p=1, so H = 0 exactly (the global minimum)");

                SequenceComplexity.CalculateShannonEntropy(new DnaSequence(homopolymer)).Should()
                    .BeApproximately(0.0, Tolerance,
                        because: $"the DnaSequence facade must also report H = 0 for the homopolymer '{homopolymer}'");
            }
        }
    }

    /// <summary>
    /// MR22-b: introducing a SECOND distinct base lifts entropy strictly above the
    /// single-symbol minimum (H &gt; 0), because the support grows from k=1 to k=2 and a
    /// two-symbol distribution can no longer be deterministic. The 50/50 two-symbol case is
    /// additionally pinned to its exact value H = log₂2 = 1 bit (maximum entropy at k=2).
    /// </summary>
    [Test]
    public void ShannonEntropy_AddingSecondSymbol_LiftsEntropyAboveZero()
    {
        // A single foreign base among many copies of one base: H must exceed the 0 minimum.
        foreach (int len in new[] { 8, 50, 256 })
        {
            string nearHomopolymer = new string('A', len - 1) + "C";

            SequenceComplexity.CalculateShannonEntropy(nearHomopolymer).Should().BeGreaterThan(0.0 + Tolerance,
                because: $"adding a single 'C' to {len - 1} 'A's grows the support to k=2, so H > 0 (no longer deterministic)");
        }

        // Exact 50/50 two-symbol endpoint: maximum entropy at k=2 is log₂2 = 1 bit.
        foreach (var pair in new[] { "AT", "GC", "AG", "CT" })
        {
            string balancedTwoSymbol = string.Concat(Enumerable.Repeat(pair, 32));

            SequenceComplexity.CalculateShannonEntropy(balancedTwoSymbol).Should().BeApproximately(1.0, Tolerance,
                because: $"an equiprobable two-symbol sequence ('{pair}'×32) attains H = log₂2 = 1 bit exactly");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GCSKEW-001 — GC skew
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (docs/algorithms/Sequence_Composition/GC_Skew.md §2.2; Lobry 1996;
    //   Grigoriev 1998 "cumulative skew diagrams"):
    //
    //       GC skew = (G − C) / (G + C)        (∈ [-1, 1]; 0 when G + C = 0)
    //
    //   The scalar skew is a function of only the G and C COUNTS. Three theory-pinned
    //   metamorphic relations follow from this definition alone:
    //
    //   • SYM (complement flips sign): the DNA complement maps G↔C (and A↔T), so it
    //     swaps the G and C counts. (G − C) negates while (G + C) is invariant, hence
    //     skew(complement(x)) = −skew(x) EXACTLY. The same swap negates every window
    //     skew, so each entry of the cumulative-skew profile negates as well. (Here we
    //     swap only G↔C — the relation depends solely on the G/C subcomposition, and
    //     A↔T does not touch the GC skew — but the canonical DnaSequence.Complement()
    //     does the full A↔T,C↔G map and is asserted too.)
    //
    //   • INV (reverse relates the cumulative profile): the implemented cumulative
    //     profile (CalculateCumulativeGcSkew) sums consecutive non-overlapping window
    //     skews (stepSize = windowSize internally). With windowSize = 1 every base is
    //     its own window, so the per-window skew sequence is exactly +1 (G), −1 (C),
    //     0 (A/T) and the CumulativeGcSkew column is the running prefix sum (#G − #C).
    //     Reversing the sequence reverses the ORDER of accumulation, so:
    //       – the per-window skew sequence is the exact order-reversal of the original
    //         (entry i of the reverse equals entry n−1−i of the original), and
    //       – the FINAL cumulative value is UNCHANGED, because the total (#G − #C) over
    //         the whole sequence is order-independent.
    //     We encode THAT precise pair of guarantees (reversed per-window profile; equal
    //     final cumulative total) — not a naive "profile negates", which the cumulative
    //     definition does not satisfy under reversal.
    //
    //   • INV (all-G → max positive): an all-G sequence has C = 0, so skew = G/G = +1
    //     EXACTLY — the documented maximum (INV-01: −1 ≤ skew ≤ 1). Symmetrically an
    //     all-C sequence gives skew = −1 EXACTLY (the minimum).
    //
    // API surface under test (src/.../Seqeron.Genomics.Analysis/GcSkewCalculator.cs):
    //   • GcSkewCalculator.CalculateGcSkew(string)      — scalar (G−C)/(G+C).
    //   • GcSkewCalculator.CalculateGcSkew(DnaSequence) — facade (same core).
    //   • GcSkewCalculator.CalculateCumulativeGcSkew(string, windowSize) — IEnumerable
    //     of CumulativeGcSkewPoint(Position, GcSkew, CumulativeGcSkew); non-overlapping
    //     windows (internal stepSize = windowSize); trailing partial window is dropped.

    #region SEQ-GCSKEW-001 — GC skew

    #region MR23: SYM — complement flips the skew sign

    /// <summary>
    /// MR23-a: skew(swapGC(x)) == −skew(x).
    /// Swapping G↔C exchanges the G and C counts, so the numerator (G − C) negates while
    /// the denominator (G + C) is invariant; the scalar GC skew therefore flips sign
    /// exactly. Verified on fixed and fixed-seed random sequences.
    /// </summary>
    [Test]
    public void GcSkew_GcSwap_FlipsSign()
    {
        foreach (var s in SampleSequences())
        {
            string swapped = new(s.Select(c => c switch { 'G' => 'C', 'C' => 'G', _ => c }).ToArray());

            double original = GcSkewCalculator.CalculateGcSkew(s);
            double flipped = GcSkewCalculator.CalculateGcSkew(swapped);

            flipped.Should().BeApproximately(-original, Tolerance,
                because: $"swapping G↔C in '{s}' negates (G−C) and leaves (G+C), so the skew must flip sign");
        }
    }

    /// <summary>
    /// MR23-b: skew(complement(x)) == −skew(x) via the canonical DnaSequence.Complement()
    /// (full A↔T, C↔G map). The A↔T half does not touch G/C, and the C↔G half swaps the
    /// G and C counts, so the GC skew negates exactly. Confirms the sign-flip end-to-end
    /// through the production complement, not just a hand-rolled G↔C swap.
    /// </summary>
    [Test]
    public void GcSkew_Complement_FlipsSign()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            double original = GcSkewCalculator.CalculateGcSkew(seq);
            double complemented = GcSkewCalculator.CalculateGcSkew(seq.Complement());

            complemented.Should().BeApproximately(-original, Tolerance,
                because: $"complement maps C↔G (swapping the G and C counts) and A↔T (no GC effect), so skew of '{s}' negates");
        }
    }

    /// <summary>
    /// MR23-c: each entry of the cumulative-skew profile negates under complement.
    /// Because every window's scalar skew flips sign, both the per-window GcSkew column
    /// and the running CumulativeGcSkew column negate position-for-position. Uses
    /// windowSize = 1 (each base is its own window) so the profile is the full per-base
    /// trace.
    /// </summary>
    [Test]
    public void GcSkew_CumulativeProfile_Complement_NegatesEachEntry()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            var original = GcSkewCalculator.CalculateCumulativeGcSkew(seq.Sequence, windowSize: 1).ToList();
            var complemented = GcSkewCalculator.CalculateCumulativeGcSkew(seq.Complement().Sequence, windowSize: 1).ToList();

            complemented.Should().HaveCount(original.Count,
                because: $"complement is a per-position rewrite, so '{s}' and its complement yield the same number of windows");

            for (int i = 0; i < original.Count; i++)
            {
                complemented[i].GcSkew.Should().BeApproximately(-original[i].GcSkew, Tolerance,
                    because: $"window {i} of complement('{s}') swaps G↔C, so its skew negates");
                complemented[i].CumulativeGcSkew.Should().BeApproximately(-original[i].CumulativeGcSkew, Tolerance,
                    because: $"each window skew negates, so the running cumulative skew at {i} negates too for '{s}'");
            }
        }
    }

    #endregion

    #region MR24: INV — reverse relates the cumulative skew profile

    /// <summary>
    /// MR24-a: with windowSize = 1 the per-window skew sequence of reverse(x) is the exact
    /// order-reversal of that of x. Each base is its own window, so reversing the sequence
    /// reverses the order in which the +1/−1/0 window skews are emitted — entry i of the
    /// reversed profile equals entry n−1−i of the original. This is the precise INV the
    /// cumulative definition guarantees under reversal (the accumulation order is reversed).
    /// </summary>
    [Test]
    public void GcSkew_CumulativeProfile_Reverse_ReversesPerWindowSkews()
    {
        foreach (var s in SampleSequences())
        {
            string reversed = new(s.Reverse().ToArray());

            var forward = GcSkewCalculator.CalculateCumulativeGcSkew(s, windowSize: 1)
                .Select(p => p.GcSkew).ToList();
            var backward = GcSkewCalculator.CalculateCumulativeGcSkew(reversed, windowSize: 1)
                .Select(p => p.GcSkew).ToList();

            backward.Should().HaveCount(forward.Count,
                because: $"windowSize = 1 emits one window per base, so '{s}' and its reverse have equal window counts");

            for (int i = 0; i < forward.Count; i++)
            {
                backward[i].Should().BeApproximately(forward[forward.Count - 1 - i], Tolerance,
                    because: $"reversing '{s}' reverses the accumulation order, so per-window skew {i} of the reverse equals window {forward.Count - 1 - i} of the original");
            }
        }
    }

    /// <summary>
    /// MR24-b: the FINAL cumulative skew value is invariant under reversal. The last
    /// CumulativeGcSkew (with windowSize = 1) is the total #G − #C over the whole sequence,
    /// which is order-independent, so reversing the sequence leaves it unchanged even though
    /// the intermediate trace is reversed (MR24-a). Empty profiles (no full window) are
    /// skipped — there is no final value to compare.
    /// </summary>
    [Test]
    public void GcSkew_CumulativeProfile_Reverse_PreservesFinalCumulative()
    {
        foreach (var s in SampleSequences())
        {
            string reversed = new(s.Reverse().ToArray());

            var forward = GcSkewCalculator.CalculateCumulativeGcSkew(s, windowSize: 1).ToList();
            var backward = GcSkewCalculator.CalculateCumulativeGcSkew(reversed, windowSize: 1).ToList();

            if (forward.Count == 0)
                continue;

            backward[^1].CumulativeGcSkew.Should().BeApproximately(forward[^1].CumulativeGcSkew, Tolerance,
                because: $"the final cumulative skew is the total #G − #C of '{s}', which is order-independent, so reversal preserves it");
        }
    }

    #endregion

    #region MR25: INV — all-G → max positive (+1); all-C → min (−1)

    /// <summary>
    /// MR25: an all-G sequence has C = 0, so skew = G/G = +1 EXACTLY — the documented
    /// maximum (INV-01: −1 ≤ skew ≤ 1). Symmetrically an all-C sequence gives −1 EXACTLY.
    /// These are exact theory endpoints, verified across several lengths through both the
    /// string and DnaSequence APIs.
    /// </summary>
    [Test]
    public void GcSkew_Homopolymer_ReachesExactEndpoints()
    {
        foreach (int len in new[] { 1, 4, 16, 100 })
        {
            string allG = new('G', len);
            string allC = new('C', len);

            GcSkewCalculator.CalculateGcSkew(allG).Should().BeApproximately(1.0, Tolerance,
                because: $"an all-G sequence of length {len} has C = 0, so skew = G/G = +1 (the documented maximum)");
            GcSkewCalculator.CalculateGcSkew(new DnaSequence(allG)).Should().BeApproximately(1.0, Tolerance,
                because: $"the DnaSequence facade must also report +1 for an all-G sequence of length {len}");

            GcSkewCalculator.CalculateGcSkew(allC).Should().BeApproximately(-1.0, Tolerance,
                because: $"an all-C sequence of length {len} has G = 0, so skew = −C/C = −1 (the documented minimum)");
            GcSkewCalculator.CalculateGcSkew(new DnaSequence(allC)).Should().BeApproximately(-1.0, Tolerance,
                because: $"the DnaSequence facade must also report −1 for an all-C sequence of length {len}");
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-ATSKEW-001 — AT skew
    // ═══════════════════════════════════════════════════════════════════
    //
    // API under test (src/.../Seqeron.Genomics.Analysis/GcSkewCalculator.cs):
    //   • GcSkewCalculator.CalculateAtSkew(string)      — scalar (A − T) / (A + T).
    //   • GcSkewCalculator.CalculateAtSkew(DnaSequence) — facade (same core).
    //
    // Relations (derived from the (A − T)/(A + T) definition, NOT from output):
    //   • SYM (complement reverses sign): complement maps A↔T, swapping the A and T counts, so the
    //         numerator (A − T) negates while the denominator (A + T) is invariant — the skew flips
    //         sign exactly. (Sequences with no A and no T have skew 0 = −0, consistent.)
    //   • INV (cumulative length = seq length): AT skew is a composition statistic. Partition the
    //         sequence into consecutive non-overlapping chunks that TILE it — the chunk lengths sum
    //         to the sequence length — and the per-chunk A and T counts accumulate to the whole's, so
    //         the skew recomposed from the cumulative counts equals the whole-sequence skew. (As a
    //         corollary the statistic is order-independent: any permutation preserves it.)

    #region SEQ-ATSKEW-001 — AT skew

    #region SYM — complement flips the AT-skew sign

    [Test]
    public void AtSkew_Complement_FlipsSign()
    {
        foreach (var s in SampleSequences())
        {
            var seq = new DnaSequence(s);
            double original = GcSkewCalculator.CalculateAtSkew(seq);
            double complemented = GcSkewCalculator.CalculateAtSkew(seq.Complement());

            complemented.Should().BeApproximately(-original, Tolerance,
                because: $"complement maps A↔T (swapping the A and T counts) and C↔G (no AT effect), so the AT skew of '{s}' negates");
        }
    }

    #endregion

    #region INV — permutation leaves the AT skew unchanged (composition statistic)

    [Test]
    public void AtSkew_Shuffle_PreservesSkew()
    {
        foreach (var s in SampleSequences())
        {
            double original = GcSkewCalculator.CalculateAtSkew(s);

            for (int t = 0; t < 5; t++)
                GcSkewCalculator.CalculateAtSkew(Shuffle(s)).Should().BeApproximately(original, Tolerance,
                    because: $"AT skew depends only on the A and T counts, not their order, so shuffling '{s}' must not change it");
        }
    }

    #endregion

    #region INV — cumulative tiling: chunk lengths sum to seq length and conserve the skew

    [Test]
    public void AtSkew_PartitionTiling_ConservesSkewAndLength()
    {
        foreach (var s in SampleSequences())
        {
            string upper = s.ToUpperInvariant();

            // Split into consecutive non-overlapping chunks that tile the whole sequence.
            var chunks = new System.Collections.Generic.List<string>();
            for (int i = 0; i < upper.Length; i += 3)
                chunks.Add(upper.Substring(i, Math.Min(3, upper.Length - i)));

            // "cumulative length = seq length": the chunk lengths exhaust the sequence.
            chunks.Sum(c => c.Length).Should().Be(upper.Length,
                because: $"the partition of '{s}' tiles it, so the cumulative chunk length equals the sequence length");

            // A and T counts accumulate over the tiling to the whole's counts.
            int aWhole = upper.Count(c => c == 'A');
            int tWhole = upper.Count(c => c == 'T');
            chunks.Sum(c => c.Count(ch => ch == 'A')).Should().Be(aWhole, because: "A counts are additive over a tiling partition");
            chunks.Sum(c => c.Count(ch => ch == 'T')).Should().Be(tWhole, because: "T counts are additive over a tiling partition");

            // The skew recomposed from the cumulative counts equals the whole-sequence skew.
            int denom = aWhole + tWhole;
            double recomposed = denom > 0 ? (double)(aWhole - tWhole) / denom : 0;
            recomposed.Should().BeApproximately(GcSkewCalculator.CalculateAtSkew(s), Tolerance,
                because: $"the cumulative A and T counts over the tiling reproduce the AT skew of '{s}'");
        }
    }

    #endregion

    #endregion
}
