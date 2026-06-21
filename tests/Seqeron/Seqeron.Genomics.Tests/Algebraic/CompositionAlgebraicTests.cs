using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

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
}
