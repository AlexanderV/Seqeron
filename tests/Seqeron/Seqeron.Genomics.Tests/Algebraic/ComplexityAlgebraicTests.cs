using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Complexity area (windowed complexity profile).
///
/// Algebraic testing pins the window-count conservation and determinism of the
/// sliding-window complexity profile.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 231.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Complexity")]
public class ComplexityAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-COMPLEX-WINDOW-001 — Windowed complexity profile (Complexity)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 231.
    //
    // Model: a complexity profile slides a fixed window across the sequence and
    //        emits one complexity point per fully contained window; with step 1 the
    //        number of windows is exactly len − w + 1.
    //   — docs/algorithms/Complexity; SequenceComplexity.CalculateWindowedComplexity.
    //
    // Laws under test (checklist row 231):
    //   • DIST  — window count = len − w + 1 (with step 1).
    //   • IDEMP — deterministic: re-running yields the identical profile.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property WindowedComplexity_Distributive_WindowCount()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 16), seq =>
        {
            const int w = 10;
            int count = SequenceComplexity
                .CalculateWindowedComplexity(new DnaSequence(seq), windowSize: w, stepSize: 1).Count();
            return (count == seq.Length - w + 1).Label($"count={count}, expected={seq.Length - w + 1}");
        });
    }

    [FsCheck.NUnit.Property]
    public Property WindowedComplexity_Idempotent_Deterministic()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 16), seq =>
        {
            const int w = 10;
            var a = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), w, 1)
                .Select(p => p.ShannonEntropy).ToList();
            var b = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), w, 1)
                .Select(p => p.ShannonEntropy).ToList();
            return a.SequenceEqual(b).ToProperty();
        });
    }
}
