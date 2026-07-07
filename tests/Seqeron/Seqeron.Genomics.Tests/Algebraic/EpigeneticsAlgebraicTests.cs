using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

using MethylationSite = EpigeneticsAnalyzer.MethylationSite;
using MethylationType = EpigeneticsAnalyzer.MethylationType;

/// <summary>
/// Algebraic-law tests for the Epigenetics area (CpG O/E ratio, DMR detection).
///
/// Algebraic testing pins the closed-form O/E ratio identity and its zero on a
/// CpG-free sequence, plus the identity/symmetry of differential-methylation
/// calling.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 85, 184.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Epigenetics")]
public class EpigeneticsAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-CPG-001 — CpG observed/expected ratio (Epigenetics)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 85.
    //
    // Model: Gardiner-Garden & Frommer O/E = observed_CpG / ((nC × nG) / length).
    //        With no CpG dinucleotide the observed count is 0, so the ratio is 0.
    //   — docs/algorithms/Epigenetics; EpigeneticsAnalyzer.CalculateCpGObservedExpected.
    //
    // Laws (row 85): ID — no CG dinucleotide → ratio 0.
    //                DIST — the O/E formula: ratio = obsCpG / ((C×G)/length).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ID: a sequence with no C-then-G has CpG O/E ratio 0.</summary>
    [Test]
    public void Cpg_Identity_NoCpGIsZero()
    {
        // "GGGCCC" contains C and G but no "CG" dinucleotide (only "GC").
        EpigeneticsAnalyzer.CalculateCpGObservedExpected("GGGCCC").Should().Be(0.0);
        EpigeneticsAnalyzer.CalculateCpGObservedExpected("AAAATTTT").Should().Be(0.0);
    }

    /// <summary>DIST: the returned ratio equals the Gardiner-Garden &amp; Frommer
    /// closed form obsCpG / ((C×G)/length).</summary>
    [FsCheck.NUnit.Property]
    public Property Cpg_Distributive_ObservedExpectedFormula()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 2), seq =>
        {
            var s = seq.ToUpperInvariant();
            int c = s.Count(ch => ch == 'C');
            int g = s.Count(ch => ch == 'G');
            int cpg = Enumerable.Range(0, s.Length - 1).Count(i => s[i] == 'C' && s[i + 1] == 'G');
            double expected = (double)c * g / s.Length;
            double reference = expected > 0 ? cpg / expected : 0;
            double actual = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            return (System.Math.Abs(actual - reference) < 1e-9)
                .Label($"ratio {actual} != {reference} for \"{seq}\"");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-DMR-001 — Differentially methylated regions (Epigenetics)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 184.
    //
    // Model: a DMR is a window whose mean methylation difference (sample2 − sample1)
    //        exceeds a cutoff. Identical methylomes have zero difference everywhere
    //        (no DMR); swapping the two samples negates each difference but leaves
    //        its magnitude — hence the called regions — unchanged.
    //   — docs/algorithms/Epigenetics; EpigeneticsAnalyzer.FindDMRs.
    //
    // Laws (row 184): ID — identical methylomes → no DMR.
    //                 COMM — |Δ| symmetric: the DMR regions of (s1,s2) and (s2,s1)
    //                 coincide, with opposite-signed but equal-magnitude differences.
    // ═══════════════════════════════════════════════════════════════════════

    private static List<MethylationSite> Methylome(double level) =>
        new[] { 100, 200, 300, 400 }
            .Select(p => new MethylationSite(p, MethylationType.CpG, "CpG", level, 20))
            .ToList();

    [Test]
    public void Dmr_Identity_IdenticalMethylomesHaveNoDmr()
    {
        var sample = Methylome(0.5);
        EpigeneticsAnalyzer.FindDMRs(sample, sample).Should().BeEmpty();
    }

    [Test]
    public void Dmr_Commutative_AbsoluteDifferenceSymmetric()
    {
        var low = Methylome(0.1);
        var high = Methylome(0.9);

        var forward = EpigeneticsAnalyzer.FindDMRs(low, high).ToList();
        var reverse = EpigeneticsAnalyzer.FindDMRs(high, low).ToList();

        forward.Should().NotBeEmpty();
        reverse.Select(d => (d.Start, d.End)).Should()
            .BeEquivalentTo(forward.Select(d => (d.Start, d.End)));

        // Same regions, opposite-signed but equal-magnitude mean differences.
        for (int i = 0; i < forward.Count; i++)
        {
            System.Math.Abs(forward[i].MeanDifference).Should()
                .BeApproximately(System.Math.Abs(reverse[i].MeanDifference), 1e-9);
            (forward[i].MeanDifference * reverse[i].MeanDifference).Should().BeLessThan(0);
        }
    }
}
