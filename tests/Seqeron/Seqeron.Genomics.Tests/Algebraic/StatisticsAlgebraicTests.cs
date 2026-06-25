using System.Collections.Generic;
using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Statistics area (Wallace Tm, codon frequencies,
/// entropy profile, GC-content profile).
///
/// Algebraic testing pins the length-linearity of the Wallace Tm, the per-amino /
/// per-codon normalization of codon frequencies, and the complement-invariance and
/// window-count conservation of the sliding-window profiles.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 130, 227, 232, 234.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Statistics")]
public class StatisticsAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    private static string Complement(string seq) => new DnaSequence(seq).Complement().Sequence;

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-TM-001 — Melting temperature (Statistics), row 130.
    // ID — Tm("") = 0.  HOMO — homopolymer Wallace Tm scales linearly with length.
    //   — SequenceStatistics.CalculateMeltingTemperature (Wallace 2(A+T)+4(G+C)).
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Tm_Identity_EmptyIsZero()
    {
        SequenceStatistics.CalculateMeltingTemperature("").Should().Be(0.0);
    }

    /// <summary>
    /// HOMO: within the Wallace regime (length &lt; 14) a homopolymer's Tm is exactly
    /// linear in length — poly-A gives 2·L, poly-G gives 4·L.
    /// </summary>
    [Test]
    public void Tm_Homogeneous_HomopolymerLinearInLength()
    {
        for (int len = 1; len <= 13; len++)
        {
            SequenceStatistics.CalculateMeltingTemperature(new string('A', len)).Should().Be(2.0 * len);
            SequenceStatistics.CalculateMeltingTemperature(new string('G', len)).Should().Be(4.0 * len);
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-CODON-FREQ-001 — Codon frequencies (Statistics), row 227.
    // ID — freq("") = ∅.  IDEMP — deterministic.  DIST — Σ counts = ⌊len/3⌋.
    //   — SequenceStatistics.CalculateCodonFrequencies.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void CodonFreq_Identity_EmptyIsEmpty()
    {
        SequenceStatistics.CalculateCodonFrequencies("").Should().BeEmpty();
    }

    /// <summary>DIST: the frequencies sum to 1, and recovering integer counts via
    /// the codon total ⌊len/3⌋ accounts for exactly ⌊len/3⌋ codons.</summary>
    [FsCheck.NUnit.Property]
    public Property CodonFreq_Distributive_CountsSumToCodonTotal()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 3), seq =>
        {
            var freqs = SequenceStatistics.CalculateCodonFrequencies(seq);
            int totalCodons = seq.Length / 3;
            double freqSum = freqs.Values.Sum();
            // Pure A/C/G/T input has no ambiguous codons, so every codon is counted.
            int recoveredCounts = freqs.Values.Sum(f => (int)System.Math.Round(f * totalCodons));
            return (System.Math.Abs(freqSum - 1.0) < 1e-9 && recoveredCounts == totalCodons)
                .Label($"freqSum={freqSum}, recovered={recoveredCounts}, total={totalCodons}");
        });
    }

    [FsCheck.NUnit.Property]
    public Property CodonFreq_Idempotent_Deterministic()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 3), seq =>
        {
            var a = SequenceStatistics.CalculateCodonFrequencies(seq);
            var b = SequenceStatistics.CalculateCodonFrequencies(seq);
            return a.OrderBy(k => k.Key).SequenceEqual(b.OrderBy(k => k.Key)).ToProperty();
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-ENTROPY-PROFILE-001 — Sliding Shannon entropy profile (Statistics), row 232.
    // IDEMP — deterministic.  INVAR — complement-invariant.  DIST — length = len−w+1.
    //   — SequenceStatistics.CalculateEntropyProfile.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property EntropyProfile_Distributive_LengthIsWindowCount()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 12), seq =>
        {
            const int w = 8;
            int count = SequenceStatistics.CalculateEntropyProfile(seq, windowSize: w, stepSize: 1).Count();
            return (count == seq.Length - w + 1).Label($"count={count}, expected={seq.Length - w + 1}");
        });
    }

    [FsCheck.NUnit.Property]
    public Property EntropyProfile_Invariant_ComplementInvariant()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 12), seq =>
        {
            const int w = 8;
            var orig = SequenceStatistics.CalculateEntropyProfile(seq, w, 1).ToList();
            var comp = SequenceStatistics.CalculateEntropyProfile(Complement(seq), w, 1).ToList();
            bool ok = orig.Count == comp.Count
                && orig.Zip(comp, (a, b) => System.Math.Abs(a - b) < 1e-9).All(x => x);
            return ok.Label($"entropy profile not complement-invariant for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-GC-PROFILE-001 — Sliding GC-content profile (Statistics), row 234.
    // IDEMP — deterministic.  INVAR — complement-invariant.  DIST — length = len−w+1.
    //   — SequenceStatistics.CalculateGcContentProfile.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property GcProfile_Distributive_LengthIsWindowCount()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 12), seq =>
        {
            const int w = 8;
            int count = SequenceStatistics.CalculateGcContentProfile(seq, windowSize: w, stepSize: 1).Count();
            return (count == seq.Length - w + 1).Label($"count={count}, expected={seq.Length - w + 1}");
        });
    }

    [FsCheck.NUnit.Property]
    public Property GcProfile_Invariant_ComplementInvariant()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 12), seq =>
        {
            const int w = 8;
            var orig = SequenceStatistics.CalculateGcContentProfile(seq, w, 1).ToList();
            var comp = SequenceStatistics.CalculateGcContentProfile(Complement(seq), w, 1).ToList();
            bool ok = orig.Count == comp.Count
                && orig.Zip(comp, (a, b) => System.Math.Abs(a - b) < 1e-9).All(x => x);
            return ok.Label($"GC profile not complement-invariant for \"{seq}\"");
        });
    }
}
