using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Composition area.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What algebraic testing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Algebraic testing pins the formal laws an operation must obey for EVERY input,
/// expressed as equations rather than point examples: identity (ID), idempotence
/// (IDEMP), commutativity (COMM), involution (INV), distributivity / conservation
/// (DIST), round-trip / isomorphism (RT), triangle inequality (TRI). A law that
/// holds is a structural guarantee about the algorithm, not a single green dot.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description;
///   docs/ADVANCED_TESTING_CHECKLIST.md §4 "Algebraic Testing".
///
/// Each law below is derived from the algorithm's published model and is asserted
/// over FsCheck-generated DNA (for the ∀-quantified laws) and over explicit
/// boundary witnesses (for the identity / neutral-element laws).
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Composition")]
public class CompositionAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-GC-001 — GC content (Composition)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 1.
    //
    // Model: GC% = (nG + nC) / (nA + nT + nG + nC) × 100, case-insensitive;
    //        empty / no-canonical-base input → 0 (defined, not an exception).
    //   — docs/algorithms/Sequence_Composition/Sequence_Composition.md §2.2, §3.3;
    //     SequenceExtensions.CalculateGcContent / DnaSequence.GcContent().
    //
    // Laws under test (checklist row 1):
    //   • ID    — GC("") = 0            (neutral input is the additive identity).
    //   • IDEMP — GC(seq) is stable on recomputation (pure, side-effect-free).
    //   • DIST  — GC(seq) = GC(complement(seq))  (complement swaps A↔T and C↔G;
    //             the C↔G swap preserves the GC count, so GC% is invariant).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: GC of the empty sequence is the additive identity 0.
    /// Evidence: empty input returns an all-zero composition, no exception
    /// (Sequence_Composition.md §3.3).
    /// </summary>
    [Test]
    public void Gc_Identity_EmptySequenceIsZero()
    {
        new DnaSequence(string.Empty).GcContent().Should().Be(0.0);
        string.Empty.AsSpan().CalculateGcContent().Should().Be(0.0);
    }

    /// <summary>
    /// IDEMP: GC content is a pure function — recomputing on the same sequence
    /// yields the identical value (no hidden state, no drift between calls).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Gc_Idempotent_StableOnRecompute()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            double first = dna.GcContent();
            double second = dna.GcContent();
            return (first == second)
                .Label($"GC recompute drifted: {first} vs {second} for \"{seq}\"");
        });
    }

    /// <summary>
    /// DIST: GC content is invariant under base complement: GC(seq) = GC(complement(seq)).
    /// Evidence: complement maps A↔T (no GC contribution either way) and C↔G
    /// (one GC base for one GC base), so the (G+C) count and the canonical total
    /// are both preserved — Watson–Crick pairing, Sequence_Composition.md §2.2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Gc_Distributive_InvariantUnderComplement()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            double original = dna.GcContent();
            double complemented = dna.Complement().GcContent();
            return (original == complemented)
                .Label($"GC not complement-invariant: {original} vs {complemented} for \"{seq}\"");
        });
    }

    /// <summary>
    /// DIST witness: a worked example pinning the law on a known sequence so the
    /// equation cannot be satisfied vacuously. "GGCCATAT" has GC=4/8=50%; its
    /// complement "CCGGTATA" also has GC=4/8=50%.
    /// </summary>
    [Test]
    public void Gc_Distributive_WorkedExample()
    {
        var dna = new DnaSequence("GGCCATAT");
        dna.GcContent().Should().Be(50.0);
        dna.Complement().Sequence.Should().Be("CCGGTATA");
        dna.Complement().GcContent().Should().Be(50.0);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-COMP-001 — DNA complement (Composition)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 2.
    //
    // Model: complement maps A↔T and C↔G (Watson–Crick base pairing), applied
    //        base-by-base with no reordering.
    //   — docs/algorithms/Sequence_Composition/RNA_Complement.md §2.1–2.2;
    //     DnaSequence.Complement() / SequenceExtensions.GetComplementBase.
    //
    // Laws under test (checklist row 2):
    //   • INV — complement(complement(x)) = x. The base map is its own inverse
    //           (A→T→A, C→G→C), and complement preserves position, so the
    //           composition is the identity on the whole sequence.
    //   • ID  — complement preserves length: |complement(x)| = |x|
    //           (a length-preserving structural identity).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// INV: complement is an involution — complement(complement(x)) = x for all DNA.
    /// Evidence: the base permutation {A↔T, C↔G} is an involution; applied
    /// position-wise twice it returns the original sequence (RNA_Complement.md §2.2).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_Involution_TwiceIsIdentity()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            string twice = dna.Complement().Complement().Sequence;
            return (twice == dna.Sequence)
                .Label($"complement∘complement(\"{seq}\") = \"{twice}\"");
        });
    }

    /// <summary>
    /// ID: complement preserves length — |complement(x)| = |x| for all DNA.
    /// Evidence: the map is total and position-wise (one output base per input base).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_PreservesLength()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            return (dna.Complement().Length == dna.Length)
                .Label($"|complement| = {dna.Complement().Length}, |x| = {dna.Length}");
        });
    }

    /// <summary>
    /// INV witness: a worked example so the involution cannot pass vacuously.
    /// complement("ACGT") = "TGCA"; complement("TGCA") = "ACGT".
    /// </summary>
    [Test]
    public void Complement_Involution_WorkedExample()
    {
        var dna = new DnaSequence("ACGT");
        dna.Complement().Sequence.Should().Be("TGCA");
        dna.Complement().Complement().Sequence.Should().Be("ACGT");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-REVCOMP-001 — DNA reverse complement (Composition)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 3.
    //
    // Model: reverse complement = reverse ∘ complement (read the opposite strand
    //        5'→3'); base map A↔T, C↔G with position reversal.
    //   — docs/algorithms/Sequence_Composition/RNA_Complement.md §2.1–2.2;
    //     DnaSequence.ReverseComplement().
    //
    // Laws under test (checklist row 3):
    //   • INV — revcomp(revcomp(x)) = x. Reversal is an involution and the base
    //           complement is an involution, and the two operations commute, so
    //           their composition applied twice is the identity.
    //   • ID  — revcomp preserves length: |revcomp(x)| = |x|.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// INV: reverse complement is an involution — revcomp(revcomp(x)) = x.
    /// Evidence: reverse∘reverse = id and complement∘complement = id, and reverse
    /// commutes with the position-wise complement (RNA_Complement.md §2.2).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReverseComplement_Involution_TwiceIsIdentity()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            string twice = dna.ReverseComplement().ReverseComplement().Sequence;
            return (twice == dna.Sequence)
                .Label($"revcomp∘revcomp(\"{seq}\") = \"{twice}\"");
        });
    }

    /// <summary>
    /// ID: reverse complement preserves length — |revcomp(x)| = |x|.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReverseComplement_PreservesLength()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            return (dna.ReverseComplement().Length == dna.Length)
                .Label($"|revcomp| = {dna.ReverseComplement().Length}, |x| = {dna.Length}");
        });
    }

    /// <summary>
    /// INV witness: revcomp("AACG") = "CGTT" (complement "TTGC" reversed);
    /// applying revcomp again returns "AACG".
    /// </summary>
    [Test]
    public void ReverseComplement_Involution_WorkedExample()
    {
        var dna = new DnaSequence("AACG");
        dna.ReverseComplement().Sequence.Should().Be("CGTT");
        dna.ReverseComplement().ReverseComplement().Sequence.Should().Be("AACG");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-GCSKEW-001 — GC skew (Composition)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 7.
    //
    // Model: GC skew = (nG − nC) / (nG + nC); zero-denominator guard returns 0
    //        (empty / no G or C). Range [−1, 1].
    //   — docs/algorithms/Sequence_Composition/Sequence_Composition.md §2.2 [2][3];
    //     GcSkewCalculator.CalculateGcSkew.
    //
    // Laws under test (checklist row 7):
    //   • ID   — skew("") = 0 (zero-denominator guard → neutral value).
    //   • DIST — skew(G-only) = +1, skew(C-only) = −1 (the conservation/extremal
    //            values of (G−C)/(G+C) when one of the two GC bases is absent).
    //   • INV  — complement negates skew: skew(complement(x)) = −skew(x), because
    //            complement swaps G↔C so the numerator (G−C) flips sign while the
    //            denominator (G+C) is preserved.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: GC skew of the empty sequence is the neutral value 0 (zero-denominator
    /// guard), not an exception (GcSkewCalculator string overload).
    /// </summary>
    [Test]
    public void GcSkew_Identity_EmptyIsZero()
    {
        GcSkewCalculator.CalculateGcSkew(string.Empty).Should().Be(0.0);
        GcSkewCalculator.CalculateGcSkew(new DnaSequence(string.Empty)).Should().Be(0.0);
    }

    /// <summary>
    /// DIST: extremal conservation values — an all-G sequence has skew +1 and an
    /// all-C sequence has skew −1, the boundaries of (G−C)/(G+C).
    /// </summary>
    [Test]
    public void GcSkew_Distributive_ExtremalValues()
    {
        GcSkewCalculator.CalculateGcSkew(new DnaSequence("GGGGGG")).Should().Be(1.0);
        GcSkewCalculator.CalculateGcSkew(new DnaSequence("CCCCCC")).Should().Be(-1.0);
    }

    /// <summary>
    /// INV: complement negates GC skew — skew(complement(x)) = −skew(x).
    /// Evidence: complement swaps G↔C, flipping the sign of (G−C) while leaving
    /// (G+C) unchanged. AT-only sequences (skew 0) satisfy −0 == 0 trivially.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcSkew_Involution_ComplementNegatesSkew()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            double skew = GcSkewCalculator.CalculateGcSkew(dna);
            double compSkew = GcSkewCalculator.CalculateGcSkew(dna.Complement());
            return (Math.Abs(compSkew - (-skew)) < 1e-12)
                .Label($"skew={skew}, skew(complement)={compSkew} for \"{seq}\"");
        });
    }

    /// <summary>
    /// Range invariant supporting the laws above: GC skew ∈ [−1, 1] for all DNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcSkew_AlwaysInRange()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double skew = GcSkewCalculator.CalculateGcSkew(new DnaSequence(seq));
            return (skew >= -1.0 && skew <= 1.0).Label($"skew={skew} out of [-1,1] for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-GC-ANALYSIS-001 — Comprehensive GC analysis (Composition)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 233.
    //
    // Model: a comprehensive GC analysis reports the overall GC content (plus
    //        windowed skew/content); the overall GC content is the same (G+C)/total
    //        measure as SEQ-GC-001 and shares its laws.
    //   — docs/algorithms/Sequence_Composition; GcSkewCalculator.AnalyzeGcContent.
    //
    // Laws (row 233): ID — analysis("").OverallGcContent = 0.
    //                 IDEMP — deterministic on recompute.
    //                 DIST — OverallGcContent(seq) = OverallGcContent(complement(seq)).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ID: the empty sequence yields zero overall GC content.</summary>
    [Test]
    public void GcAnalysis_Identity_EmptyIsZero()
    {
        GcSkewCalculator.AnalyzeGcContent(string.Empty).OverallGcContent.Should().Be(0.0);
    }

    /// <summary>IDEMP: the overall GC content is stable on recomputation.</summary>
    [FsCheck.NUnit.Property]
    public Property GcAnalysis_Idempotent_StableOnRecompute()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double first = GcSkewCalculator.AnalyzeGcContent(seq).OverallGcContent;
            double second = GcSkewCalculator.AnalyzeGcContent(seq).OverallGcContent;
            return (first == second).Label($"GC analysis drifted: {first} vs {second}");
        });
    }

    /// <summary>DIST: overall GC content is invariant under base complement.</summary>
    [FsCheck.NUnit.Property]
    public Property GcAnalysis_Distributive_InvariantUnderComplement()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double original = GcSkewCalculator.AnalyzeGcContent(seq).OverallGcContent;
            double complemented = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(seq).Complement().Sequence)
                .OverallGcContent;
            return (System.Math.Abs(original - complemented) < 1e-12)
                .Label($"GC analysis not complement-invariant: {original} vs {complemented}");
        });
    }
}
