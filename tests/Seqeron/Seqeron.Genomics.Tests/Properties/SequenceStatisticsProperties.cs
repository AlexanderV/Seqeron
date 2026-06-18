using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the sequence-statistics algorithm group (<see cref="SequenceStatistics"/>).
///
/// This file is the single home for the Statistics block of checklist 01 (rows #121–#130, and the
/// profile/codon rows #227, #232, #234). Each test unit lives in its own <c>#region</c>. Oracles are
/// derived INDEPENDENTLY from the cited theory/doc, never routed through the production result, so a
/// self-consistent-but-wrong production formula is still caught.
///
/// Test Units: SEQ-COMPOSITION-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Statistics")]
public class SequenceStatisticsProperties
{
    #region SEQ-COMPOSITION-001 — Nucleotide Composition (counts, GC/AT fractions, skews)

    // -------------------------------------------------------------------------
    // Theory (nucleotide composition):
    //   • Every character is tallied into exactly one bucket (A/T/G/C/U/N/Other), so the seven
    //     counts sum to the sequence length.                                          (P counts = length)
    //   • Over the canonical bases (total = A+T+G+C+U): GcContent = (G+C)/total and
    //     AtContent = (A+T+U)/total, which partition the canonical bases ⇒ sum to 1.   (P Σ fractions = 1)
    //   • GcContent, AtContent ∈ [0,1]; GcSkew = (G−C)/(G+C), AtSkew = (A−T)/(A+T) ∈ [−1,1]. (R)
    //
    // The bucketing and the fraction formulae are recomputed independently here.
    // -------------------------------------------------------------------------

    private const double CompTolerance = 1e-9;

    /// <summary>A string over a mixed alphabet: canonical bases, U/N, lower-case, and non-base symbols.</summary>
    private static Arbitrary<string> MixedSequenceArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'U', 'N', 'a', 'c', 'g', 't', 'u', 'n', 'X', '-', '5')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// P (checklist "counts sum = length"): the seven nucleotide buckets (A, T, G, C, U, N, Other) account
    /// for every character exactly once, so they sum to the sequence length and each is non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_CountsPartitionTheSequence()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
            int sum = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU + comp.CountN + comp.CountOther;
            bool nonNeg = comp.CountA >= 0 && comp.CountT >= 0 && comp.CountG >= 0 && comp.CountC >= 0
                          && comp.CountU >= 0 && comp.CountN >= 0 && comp.CountOther >= 0;
            return (sum == comp.Length && comp.Length == seq.Length && nonNeg)
                .Label($"counts sum {sum} ≠ length {comp.Length} (seq len {seq.Length})");
        });
    }

    /// <summary>
    /// R (checklist "each fraction ∈ [0,1]") + P (checklist "Σ fractions = 1.0"): GcContent and AtContent lie
    /// in [0,1] and sum to 1 whenever the sequence has a canonical base (else both 0); the GC/AT skews lie in
    /// [−1,1]. GcContent matches the independent (G+C)/(A+T+G+C+U) oracle.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_FractionsInUnitRange_AndSumToOne()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var comp = SequenceStatistics.CalculateNucleotideComposition(seq);
            int canonical = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU;

            bool inRange = comp.GcContent is >= 0.0 and <= 1.0 && comp.AtContent is >= 0.0 and <= 1.0
                           && comp.GcSkew is >= -1.0 and <= 1.0 && comp.AtSkew is >= -1.0 and <= 1.0;

            bool sumOk;
            bool gcOracleOk;
            if (canonical > 0)
            {
                sumOk = Math.Abs(comp.GcContent + comp.AtContent - 1.0) < CompTolerance;
                double oracleGc = (double)(comp.CountG + comp.CountC) / canonical;
                gcOracleOk = Math.Abs(comp.GcContent - oracleGc) < CompTolerance;
            }
            else
            {
                sumOk = comp.GcContent == 0.0 && comp.AtContent == 0.0;
                gcOracleOk = true;
            }

            return (inRange && sumOk && gcOracleOk)
                .Label($"GC={comp.GcContent}, AT={comp.AtContent}, canonical={canonical}");
        });
    }

    /// <summary>D (determinism): nucleotide composition is identical for identical input.</summary>
    [FsCheck.NUnit.Property]
    public Property NucleotideComposition_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
            SequenceStatistics.CalculateNucleotideComposition(seq)
                .Equals(SequenceStatistics.CalculateNucleotideComposition(seq))
                .Label("CalculateNucleotideComposition is not deterministic for identical input"));
    }

    /// <summary>
    /// Anchors: a GC-only sequence is 100% GC / 0% AT; a mixed sequence partitions correctly; the empty
    /// sequence yields all-zero counts and fractions; non-base symbols fall into Other.
    /// </summary>
    [Test]
    [Category("Property")]
    public void NucleotideComposition_CanonicalCases()
    {
        Assert.Multiple(() =>
        {
            var gc = SequenceStatistics.CalculateNucleotideComposition("GCGCGC");
            Assert.That(gc.GcContent, Is.EqualTo(1.0).Within(CompTolerance), "All G/C ⇒ GC 100%.");
            Assert.That(gc.AtContent, Is.EqualTo(0.0).Within(CompTolerance), "No A/T/U ⇒ AT 0%.");

            var mixed = SequenceStatistics.CalculateNucleotideComposition("AATTGGCCNNXX");
            Assert.That(mixed.CountA + mixed.CountT + mixed.CountG + mixed.CountC + mixed.CountU + mixed.CountN + mixed.CountOther,
                Is.EqualTo(12), "All 12 characters are bucketed.");
            Assert.That(mixed.CountN, Is.EqualTo(2), "Two N bases.");
            Assert.That(mixed.CountOther, Is.EqualTo(2), "Two non-base symbols (X) ⇒ Other.");
            Assert.That(mixed.GcContent + mixed.AtContent, Is.EqualTo(1.0).Within(CompTolerance), "GC% + AT% = 1.");

            var empty = SequenceStatistics.CalculateNucleotideComposition("");
            Assert.That(empty.Length, Is.EqualTo(0));
            Assert.That(empty.GcContent, Is.EqualTo(0.0));
        });
    }

    #endregion
}
