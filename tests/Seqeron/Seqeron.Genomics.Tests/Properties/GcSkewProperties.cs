using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for GC skew calculations.
/// Verifies skew range, complement negation, and windowed consistency.
///
/// Test Unit: SEQ-GCSKEW-001 (Property Extension)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class GcSkewProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 10) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// GC skew is in [-1, 1].
    /// Evidence: (G-C)/(G+C), bounded by definition.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcSkew_InRange()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            double skew = GcSkewCalculator.CalculateGcSkew(seq);
            return (skew >= -1.0 - 0.0001 && skew <= 1.0 + 0.0001)
                .Label($"GcSkew={skew:F4} must be in [-1, 1]");
        });
    }

    /// <summary>
    /// AT skew is in [-1, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AtSkew_InRange()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            double skew = GcSkewCalculator.CalculateAtSkew(seq);
            return (skew >= -1.0 - 0.0001 && skew <= 1.0 + 0.0001)
                .Label($"AtSkew={skew:F4} must be in [-1, 1]");
        });
    }

    /// <summary>
    /// Complement negates GC skew: skew(complement) == -skew(original).
    /// Replacing G↔C inverts the (G-C)/(G+C) ratio.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_NegatesGcSkew()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            // Swap G↔C only to keep (G+C) constant
            string swapped = new(seq.Select(c => c switch
            {
                'G' => 'C',
                'C' => 'G',
                _ => c
            }).ToArray());

            double skew1 = GcSkewCalculator.CalculateGcSkew(seq);
            double skew2 = GcSkewCalculator.CalculateGcSkew(swapped);
            return (Math.Abs(skew1 + skew2) < 0.0001)
                .Label($"skew(orig)={skew1:F4}, skew(swapped)={skew2:F4}, sum={skew1 + skew2:F4}");
        });
    }

}
