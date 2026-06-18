using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for comparative genomics algorithms.
/// Verifies invariants drawn from the literature each algorithm implements.
///
/// Test Units: COMPGEN-ANI-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("ComparativeGenomics")]
public class ComparativeGenomicsProperties
{
    /// <summary>Generates a DNA string of exactly <paramref name="len"/> bases over {A,C,G,T}.</summary>
    private static Arbitrary<string> DnaArbitrary(int len) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= len)
            .Select(a => new string(a, 0, len))
            .ToArbitrary();

    #region COMPGEN-ANI-001: R: ANI ∈ [0,1]; S: ANI(A,B)=ANI(B,A); I: ANI(A,A)=1; D: deterministic

    // ANIb (Goris et al. 2007, Int J Syst Evol Microbiol 57:81-91) cuts the query genome into
    // consecutive fragments, places each ungapped on the reference, keeps matches passing the
    // identity / alignable-fraction cut-offs, and averages their identities. The implementation
    // returns ANI as a *fraction in [0,1]* (Goris' percentage / 100), so the perfect-identity
    // value is 1.0 rather than 100.

    // A fragment length that fits inside the short sequences used by the property generators, so
    // that whole-sequence fragmentation actually produces qualifying fragments.
    private const int FragLen = 20;

    /// <summary>
    /// INV-1 (R): ANI is always a fraction in [0,1] for arbitrary query/reference DNA.
    /// Evidence: each fragment's contribution is min(matches/fragLen, 1.0) ∈ [0,1]; ANI is the
    /// mean of such values, hence bounded by [0,1] (Goris et al. 2007).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_InUnitInterval()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), DnaArbitrary(3 * FragLen), (query, reference) =>
        {
            double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            return (ani >= 0.0 && ani <= 1.0)
                .Label($"ANI={ani} outside [0,1]");
        });
    }

    /// <summary>
    /// INV-2 (I): Self-identity. ANI(A,A) == 1.0 whenever A admits at least one fragment.
    /// Evidence: every fragment of A is a substring of A, so its best ungapped placement is exact
    /// (identity 1.0) and qualifies; the mean of 1.0 values is 1.0. This is the biological anchor
    /// that a genome is 100% identical to itself.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_SelfIdentity_IsOne()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), seq =>
        {
            double ani = ComparativeGenomics.CalculateANI(seq, seq, fragmentLength: FragLen);
            return (Math.Abs(ani - 1.0) < 1e-12)
                .Label($"ANI(A,A)={ani}, expected 1.0");
        });
    }

    /// <summary>
    /// INV-3 (S): Symmetry. When both genomes have length equal to the fragment length there is a
    /// single fragment placed at the single possible offset, so the match count A↔B equals B↔A and
    /// ANI(A,B) == ANI(B,A) exactly. (In general ANIb is only approximately symmetric because the
    /// query alone is fragmented; this configuration isolates the exact-symmetry regime.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_IsSymmetric_ForEqualLengthGenomes()
    {
        return Prop.ForAll(DnaArbitrary(FragLen), DnaArbitrary(FragLen), (a, b) =>
        {
            double ab = ComparativeGenomics.CalculateANI(a, b, fragmentLength: FragLen, minIdentity: 0.0);
            double ba = ComparativeGenomics.CalculateANI(b, a, fragmentLength: FragLen, minIdentity: 0.0);
            return (Math.Abs(ab - ba) < 1e-12)
                .Label($"ANI(A,B)={ab} ≠ ANI(B,A)={ba}");
        });
    }

    /// <summary>
    /// INV-4 (M): Divergence does not increase identity. Mutating more positions of a copy of A
    /// yields ANI no higher than mutating fewer positions: ANI(A, more-mutated) ≤ ANI(A, less-mutated).
    /// Evidence: with one full-length fragment the qualifying identity is matches/length; extra
    /// substitutions can only reduce the match count (Goris et al. 2007 identity definition).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Ani_MoreMutations_DoNotIncreaseIdentity()
    {
        const string baseSeq = "ACGTACGTACGTACGTACGT"; // length 20 == FragLen → single fragment
        string fewMutations = "TCGTACGTACGTACGTACGT"; // 1 substitution
        string moreMutations = "TCGAACGTTCGTACGTACGT"; // superset: positions 0,3,8 changed

        double aniFew = ComparativeGenomics.CalculateANI(baseSeq, fewMutations, fragmentLength: FragLen, minIdentity: 0.0);
        double aniMore = ComparativeGenomics.CalculateANI(baseSeq, moreMutations, fragmentLength: FragLen, minIdentity: 0.0);

        Assert.That(aniMore, Is.LessThanOrEqualTo(aniFew + 1e-12),
            $"More divergence raised ANI: ANI(more)={aniMore} > ANI(few)={aniFew}");
    }

    /// <summary>
    /// INV-5 (D): Determinism. The same query/reference pair yields the same ANI on repeated calls.
    /// Evidence: CalculateANI is a pure function of its inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Ani_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(3 * FragLen), DnaArbitrary(3 * FragLen), (query, reference) =>
        {
            double a1 = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            double a2 = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: FragLen);
            return (a1 == a2).Label("CalculateANI must be deterministic");
        });
    }

    /// <summary>
    /// INV-6 (boundary): ANI is 0 when no full-length fragment fits — empty inputs or a reference
    /// shorter than the fragment leave no qualifying match, so the mean is defined as 0
    /// (Goris et al. 2007: only fragments with a placeable, qualifying match contribute).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Ani_ReturnsZero_WhenNoFragmentQualifies()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ComparativeGenomics.CalculateANI("", "ACGTACGTACGT", fragmentLength: FragLen),
                Is.EqualTo(0.0), "empty query");
            Assert.That(ComparativeGenomics.CalculateANI("ACGTACGTACGT", "", fragmentLength: FragLen),
                Is.EqualTo(0.0), "empty reference");
            // Reference shorter than the fragment: no full-length placement exists.
            Assert.That(ComparativeGenomics.CalculateANI(new string('A', FragLen), "ACGT", fragmentLength: FragLen),
                Is.EqualTo(0.0), "reference shorter than fragment");
        });
    }

    #endregion
}
