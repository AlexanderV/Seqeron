// POP-ROH-001 — Runs of Homozygosity (consecutive-runs detection + F_ROH)
// Evidence: docs/Evidence/POP-ROH-001-Evidence.md
// TestSpec: tests/TestSpecs/POP-ROH-001.md
// Source: McQuillan R, et al. (2008). Am J Hum Genet 83(3):359-372 (F_ROH = ΣL_roh / L_auto).
//         Marras G, et al. (2015). Anim Genet 46(2):110-121 (consecutive-runs method).
//         Chang CC, et al. (2015). GigaScience 4:7 — PLINK 1.9 --homozyg defaults.

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Population.PopulationGeneticsAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Population;

[TestFixture]
public class PopulationGeneticsAnalyzer_FindROH_Tests
{
    private const double Tol = 1e-10;

    // Builds n consecutive SNPs at the given spacing (bp) with a homozygous genotype,
    // optionally inserting heterozygous (genotype 1) calls at the listed indices.
    private static List<(int Position, int Genotype)> Snps(
        int n, int spacing, params int[] hetIndices)
    {
        var het = new HashSet<int>(hetIndices);
        var list = new List<(int, int)>(n);
        for (int i = 0; i < n; i++)
            list.Add((i * spacing, het.Contains(i) ? 1 : 0));
        return list;
    }

    #region FindROH

    // M1 — A single uninterrupted homozygous run is reported once with exact bounds and count.
    // 100 SNPs spaced 20 kb apart span [0, 1,980,000] bp = 1.98 Mb ≥ 1 Mb (PLINK --homozyg-kb),
    // 100 SNPs ≥ 100 (PLINK --homozyg-snp); zero opposite genotypes (Marras et al. 2015).
    [Test]
    public void FindROH_SingleHomozygousRun_ReportsExactBoundsAndCount()
    {
        var snps = Snps(100, 20_000);

        var roh = FindROH(snps).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "One uninterrupted homozygous run must yield exactly one segment.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "Run starts at the first homozygous SNP position.");
            Assert.That(roh[0].End, Is.EqualTo(1_980_000), "Run ends at the last homozygous SNP position (99 × 20 kb).");
            Assert.That(roh[0].SnpCount, Is.EqualTo(100), "All 100 homozygous SNPs are in the run.");
        });
    }

    // M2 — One tolerated interior heterozygote (maxHeterozygotes = 1) does NOT break the run.
    // Marras et al. (2015): up to maxOppRun opposite genotypes are allowed inside a run.
    [Test]
    public void FindROH_OneToleratedHeterozygote_KeepsSingleRun()
    {
        // 101 SNPs, the SNP at index 50 is heterozygous; one het ≤ maxHeterozygotes = 1.
        var snps = Snps(101, 20_000, 50);

        var roh = FindROH(snps, maxHeterozygotes: 1).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "A single tolerated heterozygote must not split the run.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "Run still starts at the first SNP.");
            Assert.That(roh[0].End, Is.EqualTo(2_000_000), "Run ends at the last homozygous SNP (index 100 × 20 kb).");
            Assert.That(roh[0].SnpCount, Is.EqualTo(101), "Run spans all 101 SNPs including the one tolerated het.");
        });
    }

    // M3 — A second heterozygote beyond the tolerance splits the run into two segments,
    // each closed at its last homozygous SNP and started after the breaking het.
    [Test]
    public void FindROH_HeterozygoteBeyondTolerance_SplitsIntoTwoRuns()
    {
        // 201 SNPs; hets at indices 50 (tolerated) and 100 (the breaker). With maxHet = 1 the
        // first run is [0..99] (the het at 100 breaks it; last hom = index 99), second run
        // restarts at index 101 and ends at index 200.
        var snps = Snps(201, 20_000, 50, 100);

        var roh = FindROH(snps, minSnps: 50, minLength: 500_000, maxHeterozygotes: 1).ToList();

        Assert.That(roh, Has.Count.EqualTo(2), "A heterozygote beyond the tolerance must split the run.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "First run starts at SNP 0.");
            Assert.That(roh[0].End, Is.EqualTo(99 * 20_000), "First run ends at the last homozygous SNP before the breaker (index 99).");
            Assert.That(roh[0].SnpCount, Is.EqualTo(100), "First run holds SNPs 0..99 (incl. tolerated het at 50).");
            Assert.That(roh[1].Start, Is.EqualTo(101 * 20_000), "Second run restarts at the first homozygous SNP after the breaker.");
            Assert.That(roh[1].End, Is.EqualTo(200 * 20_000), "Second run ends at the final SNP.");
            Assert.That(roh[1].SnpCount, Is.EqualTo(100), "Second run holds SNPs 101..200.");
        });
    }

    // M4 — A physical gap larger than maxGap breaks the run even with all-homozygous SNPs
    // (PLINK --homozyg-gap default 1000 kb; Chang et al. 2015).
    [Test]
    public void FindROH_GapExceedsMaxGap_BreaksRun()
    {
        // Two blocks of 60 homozygous SNPs at 20 kb spacing, separated by a 2 Mb jump > maxGap.
        var block1 = Enumerable.Range(0, 60).Select(i => (i * 20_000, 0));
        var block2 = Enumerable.Range(0, 60).Select(i => (60 * 20_000 + 2_000_000 + i * 20_000, 0));
        var snps = block1.Concat(block2).ToList();

        var roh = FindROH(snps, minSnps: 50, minLength: 500_000, maxGap: 1_000_000).ToList();

        Assert.That(roh, Has.Count.EqualTo(2), "A gap larger than maxGap must terminate the run and start a new one.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "First run starts at SNP 0.");
            Assert.That(roh[0].End, Is.EqualTo(59 * 20_000), "First run ends before the gap.");
            Assert.That(roh[0].SnpCount, Is.EqualTo(60), "First block holds 60 SNPs.");
            Assert.That(roh[1].SnpCount, Is.EqualTo(60), "Second block holds 60 SNPs.");
        });
    }

    // M5 — A run shorter than minSnps is discarded (PLINK --homozyg-snp; Marras minSNP).
    [Test]
    public void FindROH_FewerThanMinSnps_NotReported()
    {
        var snps = Snps(20, 20_000); // 20 homozygous SNPs, < default minSnps = 100

        var roh = FindROH(snps).ToList();

        Assert.That(roh, Is.Empty, "A run with fewer than minSnps SNPs must not be reported.");
    }

    // M6 — A run shorter than minLength is discarded even when SNP count passes
    // (PLINK --homozyg-kb 1000; Chang et al. 2015).
    [Test]
    public void FindROH_ShorterThanMinLength_NotReported()
    {
        // 120 SNPs at 1 kb spacing span only 119 kb < 1 Mb default minLength.
        var snps = Snps(120, 1_000);

        var roh = FindROH(snps, minSnps: 100, minLength: 1_000_000).ToList();

        Assert.That(roh, Is.Empty, "A run shorter than minLength must not be reported even if SNP count passes.");
    }

    // S1 — Inputs need not be pre-sorted; positions are ordered ascending internally.
    [Test]
    public void FindROH_UnsortedInput_OrdersByPosition()
    {
        var sorted = Snps(100, 20_000);
        var shuffled = sorted.AsEnumerable().Reverse().ToList();

        var roh = FindROH(shuffled).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "Unsorted input must yield the same single run.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "Start is the minimum position after internal sort.");
            Assert.That(roh[0].End, Is.EqualTo(1_980_000), "End is the maximum position after internal sort.");
        });
    }

    // S2 — Empty input yields no runs.
    [Test]
    public void FindROH_EmptyInput_ReturnsEmpty()
    {
        var roh = FindROH(new List<(int, int)>()).ToList();

        Assert.That(roh, Is.Empty, "No genotypes means no runs.");
    }

    // S3 — Homozygous-alternate genotype (2) is homozygous and counts toward a run, like 0.
    [Test]
    public void FindROH_HomozygousAlternateGenotype_CountsAsHomozygous()
    {
        // All genotype 2 (homozygous alt): a valid ROH per the same rule as genotype 0.
        var snps = Enumerable.Range(0, 100).Select(i => (i * 20_000, 2)).ToList();

        var roh = FindROH(snps).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "Homozygous-alternate SNPs form a run just like homozygous-reference SNPs.");
        Assert.That(roh[0].SnpCount, Is.EqualTo(100), "All 100 homozygous-alt SNPs are in the run.");
    }

    // S4 — Leading heterozygotes are skipped; the run starts at the first homozygous SNP.
    [Test]
    public void FindROH_LeadingHeterozygotes_RunStartsAtFirstHomozygous()
    {
        var snps = Snps(102, 20_000, 0, 1); // first two SNPs heterozygous

        var roh = FindROH(snps, maxHeterozygotes: 1, minSnps: 100, minLength: 1_000_000).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "Run forms from the first homozygous SNP onward.");
        Assert.That(roh[0].Start, Is.EqualTo(2 * 20_000), "Leading heterozygotes are skipped; run starts at index 2.");
    }

    // S6 — A run ends at its last homozygous SNP; a trailing tolerated heterozygote is NOT part
    // of the emitted run. A ROH is by definition bounded by homozygous markers (Marras et al.
    // 2015: a run is a stretch of homozygous SNPs; opposite genotypes are only tolerated to
    // bridge interior error, not to extend the run past its last homozygous locus). With a het
    // at the FINAL index (100), the run must close at index 99 with SnpCount = 100, End =
    // 1,980,000 — i.e. the trailing het is excluded from both End and SnpCount.
    [Test]
    public void FindROH_TrailingHeterozygote_RunEndsAtLastHomozygousSnp()
    {
        var snps = Snps(101, 20_000, 100); // 100 hom SNPs followed by one trailing het

        var roh = FindROH(snps, maxHeterozygotes: 1, minSnps: 100, minLength: 1_000_000).ToList();

        Assert.That(roh, Has.Count.EqualTo(1), "The trailing het does not start a new run.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "Run starts at the first homozygous SNP.");
            Assert.That(roh[0].End, Is.EqualTo(1_980_000), "Run ends at the last homozygous SNP (index 99), excluding the trailing het.");
            Assert.That(roh[0].SnpCount, Is.EqualTo(100), "The trailing het is excluded from the run's SNP count.");
        });
    }

    // S7 — With zero tolerance (maxHeterozygotes = 0) any heterozygote breaks the run; this is
    // the strict consecutive-homozygous-run case. A het at index 50 splits 100 homozygous SNPs
    // into [0..49] (50 SNPs) and [51..99] (49 SNPs) — Marras et al. (2015) maxOppRun = 0.
    [Test]
    public void FindROH_ZeroToleranceHeterozygote_BreaksRunImmediately()
    {
        var snps = Snps(100, 20_000, 50); // het at the midpoint

        var roh = FindROH(snps, minSnps: 40, minLength: 500_000, maxHeterozygotes: 0).ToList();

        Assert.That(roh, Has.Count.EqualTo(2), "With zero tolerance, a single het splits the run.");
        Assert.Multiple(() =>
        {
            Assert.That(roh[0].Start, Is.EqualTo(0), "First run starts at SNP 0.");
            Assert.That(roh[0].End, Is.EqualTo(49 * 20_000), "First run ends at the last homozygous SNP before the het.");
            Assert.That(roh[0].SnpCount, Is.EqualTo(50), "First run holds the 50 homozygous SNPs 0..49.");
            Assert.That(roh[1].Start, Is.EqualTo(51 * 20_000), "Second run restarts after the breaking het.");
            Assert.That(roh[1].End, Is.EqualTo(99 * 20_000), "Second run ends at the final SNP.");
            Assert.That(roh[1].SnpCount, Is.EqualTo(49), "Second run holds the 49 homozygous SNPs 51..99.");
        });
    }

    // C1 — Invalid arguments raise ArgumentOutOfRangeException / ArgumentNullException.
    [Test]
    public void FindROH_InvalidArguments_Throw()
    {
        var snps = Snps(10, 20_000);
        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentNullException>(() => FindROH(null!).ToList(),
                "Null genotypes must throw ArgumentNullException.");
            Assert.Throws<System.ArgumentOutOfRangeException>(() => FindROH(snps, minSnps: 0).ToList(),
                "minSnps < 1 must throw.");
            Assert.Throws<System.ArgumentOutOfRangeException>(() => FindROH(snps, minLength: -1).ToList(),
                "Negative minLength must throw.");
            Assert.Throws<System.ArgumentOutOfRangeException>(() => FindROH(snps, maxHeterozygotes: -1).ToList(),
                "Negative maxHeterozygotes must throw.");
            Assert.Throws<System.ArgumentOutOfRangeException>(() => FindROH(snps, maxGap: -1).ToList(),
                "Negative maxGap must throw.");
        });
    }

    #endregion

    #region CalculateInbreedingFromROH

    // M7 — F_ROH = ΣL_roh / L_auto (McQuillan et al. 2008). Two segments of 10 Mb each over a
    // 100 Mb genome: (10M + 10M) / 100M = 0.20.
    [Test]
    public void CalculateInbreedingFromROH_TwoSegments_MatchesFRohFormula()
    {
        var segments = new List<(int Start, int End)>
        {
            (0, 10_000_000),
            (50_000_000, 60_000_000),
        };

        double f = CalculateInbreedingFromROH(segments, 100_000_000);

        Assert.That(f, Is.EqualTo(0.20).Within(Tol),
            "F_ROH = (10M + 10M) / 100M = 0.20 per McQuillan et al. (2008) ΣL_roh / L_auto.");
    }

    // M8 — A single ROH covering the whole genome gives F_ROH = 1 (full autozygosity bound).
    [Test]
    public void CalculateInbreedingFromROH_WholeGenomeROH_ReturnsOne()
    {
        var segments = new List<(int Start, int End)> { (0, 2_673_768) };

        double f = CalculateInbreedingFromROH(segments, 2_673_768);

        Assert.That(f, Is.EqualTo(1.0).Within(Tol), "A ROH spanning all of L_auto yields F_ROH = 1.");
    }

    // S5 — No ROH segments gives F_ROH = 0.
    [Test]
    public void CalculateInbreedingFromROH_NoSegments_ReturnsZero()
    {
        double f = CalculateInbreedingFromROH(new List<(int, int)>(), 300_000_000);

        Assert.That(f, Is.EqualTo(0.0).Within(Tol), "No ROH means ΣL_roh = 0, so F_ROH = 0.");
    }

    // C2 — Non-positive genome length returns 0 (no defined denominator); null segments throw.
    [Test]
    public void CalculateInbreedingFromROH_InvalidInput_HandledPerContract()
    {
        var segments = new List<(int, int)> { (0, 1_000_000) };
        Assert.Multiple(() =>
        {
            Assert.That(CalculateInbreedingFromROH(segments, 0), Is.EqualTo(0.0).Within(Tol),
                "genomeLength = 0 has no denominator; result is 0.");
            Assert.That(CalculateInbreedingFromROH(segments, -1), Is.EqualTo(0.0).Within(Tol),
                "Negative genomeLength has no denominator; result is 0.");
            Assert.Throws<System.ArgumentNullException>(() => CalculateInbreedingFromROH(null!, 100),
                "Null segments must throw ArgumentNullException.");
        });
    }

    #endregion
}
