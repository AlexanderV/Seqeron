// ONCO-ASCAT-001 — Battenberg callSubclones copy-number driver (B24 F61)
// Evidence: docs/Evidence/ONCO-ASCAT-001-Evidence.md (§ F61)
// TestSpec: tests/TestSpecs/ONCO-ASCAT-001.md (§20)
// Source: Wedge-lab/battenberg master (57a8f7e) R/fitcopynumber.R callSubclones: set.seed(seed) →
//         determine_copynumber → merge_segments → determine_copynumber → mask_high_cn_segments (one RNG stream across
//         both determine_copynumber calls); R 4.3.3 RNG.c (Mersenne-Twister, R_unif_index).
// Inputs and R expectations: BattenbergCallSubclonesData (G1–G4). G1 and G3 have sub-clonal segments in both
// determine_copynumber calls, so the final bootstrap columns (SDfrac_BS, frac1_0.025, frac1_0.975) depend on the
// draws consumed by the first call.

using static Seqeron.Genomics.Tests.Unit.Oncology.BattenbergCallSubclonesData;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_BattenbergCallSubclones_Tests
{
    private static OncologyAnalyzer.BattenbergSubcloneCalls Run(string name)
    {
        Genome g = Genomes[name];
        var (snps, logR) = Generate(g);
        return OncologyAnalyzer.CallBattenbergSubclones(
            snps, logR, g.Rho, g.Psit, g.RSeed, maxAllowedState: g.MaxState, bafOption: g.Option);
    }

    [TestCase("G1")]
    [TestCase("G2")]
    [TestCase("G3")]
    [TestCase("G4")]
    public void CallBattenbergSubclones_MatchesR(string name)
    {
        var result = Run(name);
        Expectation e = Expected[name];
        Assert.Multiple(() =>
        {
            AssertRows(result.InitialCalls, e.Initial, $"{name} determine_copynumber #1 (_1.txt)");
            Assert.That(result.Merge.Segments.Select(s => (s.Chromosome, s.Start, s.End)),
                Is.EqualTo(e.Merged.Select(m => (m.Chr, m.Start, m.End))), $"{name}: merged extents.");
            AssertRows(result.Calls, e.Final, $"{name} final (determine_copynumber #2 + mask)");
            Assert.That((result.Mask.MaskedCount, result.Mask.MaskedSize), Is.EqualTo((e.MaskedCount, e.MaskedSize)), $"{name}: masking details.");
        });
    }

    [Test]
    public void CallBattenbergSubclones_SecondCallContinuesTheSeededStream()
    {
        // G1: the final table is not the one FitSubclonalCopyNumberWithBootstrap gives with a fresh set.seed(2024) on the
        // merged rows (the F41 caveat) — the driver keeps R's stream, so its bootstrap columns differ from that and match R.
        Genome g = Genomes["G1"];
        var (_, logR) = Generate(g);
        var result = Run("G1");
        var fresh = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(
            OncologyAnalyzer.BuildBattenbergSegments(result.Merge.SegmentedSnps, logR), g.Rho, g.Psit, g.RSeed);
        int sub = result.Calls.Select((c, i) => (c, i)).First(t => t.c.Fit.IsSubclonal).i;
        Assert.Multiple(() =>
        {
            Assert.That(fresh[sub].Solutions[0].Fraction1, Is.EqualTo(result.Calls[sub].Solutions[0].Fraction1), "τ does not depend on the RNG.");
            Assert.That(fresh[sub].Solutions[0].FractionBootstrapSd, Is.Not.EqualTo(result.Calls[sub].Solutions[0].FractionBootstrapSd),
                "a re-seeded second call draws different resamples.");
            Assert.That(result.Calls[sub].Solutions[0].FractionBootstrapSd, Is.EqualTo(Expected["G1"].Final[sub].Row[11]).Within(1e-12),
                "the continued stream reproduces R's SDfrac_A_BS.");
        });
    }

    [Test]
    public void CallBattenbergSubclones_InvalidInput_Throws()
    {
        var (snps, logR) = Generate(Genomes["G2"]);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CallBattenbergSubclones(null!, logR, 0.85, 3.1, 1));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CallBattenbergSubclones(snps, null!, 0.85, 3.1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CallBattenbergSubclones(snps, logR, 0.85, 3.1, 1, permutations: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CallBattenbergSubclones(snps, logR, 1.5, 3.1, 1));
        });
    }
}
