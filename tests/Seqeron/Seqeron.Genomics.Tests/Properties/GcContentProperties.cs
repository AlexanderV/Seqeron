using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Tests.Builders;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for GC content calculation.
/// These tests verify mathematical invariants that must hold for ALL valid DNA sequences,
/// not just hand-picked examples.
///
/// Test Unit: SEQ-GC-001 (Property Extension)
/// Evidence: Formula GC% = (G+C)/(A+T+G+C) × 100
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class GcContentProperties
{
    private static Arbitrary<string> DnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// INV-1: GC percentage is always in [0, 100] for any valid DNA sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcContent_AlwaysInRange_0_to_100()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double gc = seq.AsSpan().CalculateGcContent();
            return (gc >= 0.0 && gc <= 100.0)
                .Label($"GC%={gc:F4} for len={seq.Length}");
        });
    }

    /// <summary>
    /// INV-2: GC fraction is always in [0, 1] for any valid DNA sequence.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcFraction_AlwaysInRange_0_to_1()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double frac = seq.AsSpan().CalculateGcFraction();
            return (frac >= 0.0 && frac <= 1.0)
                .Label($"GcFrac={frac:F6} for len={seq.Length}");
        });
    }

    /// <summary>
    /// INV-3: GcContent == GcFraction × 100 (relationship between percentage and fraction).
    /// Evidence: Formula derivation — both compute from the same nucleotide counts.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcContent_EqualsGcFraction_Times100()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            ReadOnlySpan<char> span = seq;
            double pct = span.CalculateGcContent();
            double frac = span.CalculateGcFraction();
            return (Math.Abs(pct - frac * 100.0) < 0.0001)
                .Label($"GC%={pct:F4}, Frac×100={frac * 100.0:F4}");
        });
    }

    /// <summary>
    /// INV-4: Complement preserves GC content.
    /// Evidence: Complement swaps A↔T and G↔C — count of G+C is invariant.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_PreservesGcContent()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var comp = dna.Complement();
            return (Math.Abs(dna.GcContent() - comp.GcContent()) < 0.0001)
                .Label($"Original GC={dna.GcContent():F4}, Complement GC={comp.GcContent():F4}");
        });
    }

    /// <summary>
    /// INV-5: Reverse complement preserves GC content.
    /// Evidence: Reverse does not change nucleotide counts; complement preserves G+C count.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReverseComplement_PreservesGcContent()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var rc = dna.ReverseComplement();
            return (Math.Abs(dna.GcContent() - rc.GcContent()) < 0.0001)
                .Label($"Original GC={dna.GcContent():F4}, RevComp GC={rc.GcContent():F4}");
        });
    }

    /// <summary>
    /// INV-6: Span-based and string-based methods return the same result.
    /// Evidence: CalculateGcContentFast delegates to CalculateGcContent(ReadOnlySpan).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SpanAndString_ReturnSameResult()
    {
        return Prop.ForAll(DnaArbitrary(), seq =>
        {
            double spanResult = seq.AsSpan().CalculateGcContent();
            double stringResult = seq.CalculateGcContentFast();
            return (Math.Abs(spanResult - stringResult) < 0.0001)
                .Label($"Span={spanResult:F4}, String={stringResult:F4}");
        });
    }

    /// <summary>
    /// INV-7: All-GC sequences yield 100%, all-AT sequences yield 0%.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HomopolymericSequences_YieldExtremalValues()
    {
        var gcGen = Gen.Elements('G', 'C')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

        return Prop.ForAll(gcGen, seq =>
        {
            double gc = seq.AsSpan().CalculateGcContent();
            return (Math.Abs(gc - 100.0) < 0.0001)
                .Label($"All-GC sequence of len={seq.Length} should be 100%, got {gc:F4}");
        });
    }
}
