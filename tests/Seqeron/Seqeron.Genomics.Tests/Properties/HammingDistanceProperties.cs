using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for Hamming Distance.
/// Verifies metric space axioms that must hold for ALL equal-length string pairs.
///
/// Test Unit: PAT-APPROX-001 (Property Extension)
/// Evidence: Wikipedia — Hamming distance is a metric on strings of equal length.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Alignment")]
public class HammingDistanceProperties
{
    private static Arbitrary<string> DnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length > 0)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Metric axiom: d(a, a) == 0 — distance from a string to itself is always zero.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Identity_DistanceToSelfIsZero()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            return (ApproximateMatcher.HammingDistance(s, s) == 0)
                .Label($"d(\"{s[..Math.Min(s.Length, 10)]}\", self) should be 0");
        });
    }

    /// <summary>
    /// Metric axiom: d(a, b) == d(b, a) — Hamming distance is symmetric.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Symmetry()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            // Create second string of same length with some mutations
            var rng = new Random(s.GetHashCode());
            var chars = s.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
                if (rng.NextDouble() < 0.3)
                    chars[i] = "ACGT"[rng.Next(4)];
            string t = new(chars);

            int dAB = ApproximateMatcher.HammingDistance(s, t);
            int dBA = ApproximateMatcher.HammingDistance(t, s);
            return (dAB == dBA).Label($"d(a,b)={dAB}, d(b,a)={dBA}");
        });
    }

    /// <summary>
    /// Non-negativity: d(a, b) >= 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NonNegativity()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            var t = new string(s.Select(c => c == 'A' ? 'T' : c).ToArray());
            return (ApproximateMatcher.HammingDistance(s, t) >= 0)
                .Label("Distance must be non-negative");
        });
    }

    /// <summary>
    /// Upper bound: d(a, b) &lt;= |a| — at most every position differs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property UpperBound()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            string complement = new(s.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', _ => 'G' }).ToArray());
            int d = ApproximateMatcher.HammingDistance(s, complement);
            return (d <= s.Length).Label($"d={d}, len={s.Length}");
        });
    }

    /// <summary>
    /// Triangle inequality: d(a, c) &lt;= d(a, b) + d(b, c).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TriangleInequality()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            if (s.Length < 2) return true.ToProperty();
            var rng = new Random(s.GetHashCode());
            var b = s.ToCharArray();
            var c = s.ToCharArray();
            for (int i = 0; i < s.Length; i++)
            {
                if (rng.NextDouble() < 0.2) b[i] = "ACGT"[rng.Next(4)];
                if (rng.NextDouble() < 0.3) c[i] = "ACGT"[rng.Next(4)];
            }
            int ab = ApproximateMatcher.HammingDistance(s, new string(b));
            int bc = ApproximateMatcher.HammingDistance(new string(b), new string(c));
            int ac = ApproximateMatcher.HammingDistance(s, new string(c));
            return (ac <= ab + bc).Label($"d(a,c)={ac} <= d(a,b)={ab} + d(b,c)={bc}");
        });
    }

    /// <summary>
    /// Complement has maximal distance: d(seq, complement) == |seq| for pure ACGT sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_HasMaximalDistance()
    {
        return Prop.ForAll(DnaArbitrary(), s =>
        {
            var dna = new DnaSequence(s);
            string comp = dna.Complement().Sequence;
            int d = ApproximateMatcher.HammingDistance(s, comp);
            return (d == s.Length).Label($"d(seq, complement)={d}, expected={s.Length}");
        });
    }
}
