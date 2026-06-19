using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Chromosome;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Chromosome area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-TELO-001 — telomere analysis (Chromosome).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 48.
///
/// API under test (ChromosomeAnalyzer.AnalyzeTelomeres):
///   Scans the 3' end for tandem TTAGGG repeats (≥70 % per-hexamer identity) counting inward,
///   and the 5' start for the reverse-complement CCCTAA; reports the terminal repeat-tract
///   LENGTHS and purities at each end.
///
/// Relations (derived from that definition, NOT from output):
///   • MON (more repeats ⇒ longer telomere): a 3' tract of k TTAGGG units (with a
///          non-telomeric interior) is measured as length 6k, strictly increasing in k; the
///          5' CCCTAA tract behaves symmetrically.
///   • INV (interior flank doesn't affect the core): inserting non-telomeric sequence BETWEEN
///          the two terminal tracts leaves both telomere lengths and purities unchanged — the
///          ends (the biological "core" of the measurement) are untouched.
///   • SHIFT/strand-duality: this API reports terminal LENGTHS, not interior coordinates, so
///          the checklist's positional-shift relation is realised as the strand duality that
///          underlies telomere calling — reverse-complementing the chromosome maps the 3'
///          TTAGGG tract onto a 5' CCCTAA tract, swapping the 5' and 3' measurements exactly.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ChromosomeMetamorphicTests
{
    #region Helpers

    private const string TeloRepeat = "TTAGGG";              // 3' tract
    private const string TeloRepeatRc = "CCCTAA";            // 5' tract (reverse complement)

    private static string Repeat(string unit, int count) => string.Concat(Enumerable.Repeat(unit, count));

    /// <summary>Non-telomeric filler whose every hexamer falls well below the 70 % match threshold for both tracts.</summary>
    private static string Filler(int length)
    {
        var sb = new System.Text.StringBuilder(length);
        while (sb.Length < length) sb.Append("GC");
        return sb.ToString(0, length);
    }

    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    private static readonly Random Rng = new(20260619);

    /// <summary>Pseudo-random non-repetitive DNA (15-mers essentially unique ⇒ ~0 repeat content).</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    #region MON — more terminal repeats give a longer telomere

    [Test]
    [Description("MON: a 3' tract of k TTAGGG units is measured as telomere length 6k, strictly increasing in k; the 5' CCCTAA tract behaves symmetrically.")]
    public void AnalyzeTelomeres_MoreRepeats_IncreasesLength()
    {
        int prev3 = int.MinValue, prev5 = int.MinValue;

        foreach (int k in new[] { 1, 2, 5, 10, 20, 50 })
        {
            var three = ChromosomeAnalyzer.AnalyzeTelomeres("chr", Filler(60) + Repeat(TeloRepeat, k));
            three.TelomereLength3Prime.Should().Be(6 * k,
                because: $"the 3' end carries {k} contiguous TTAGGG units, each contributing 6 bp to the measured telomere");
            three.TelomereLength3Prime.Should().BeGreaterThan(prev3,
                because: "adding TTAGGG units to the 3' end strictly lengthens the measured telomere");
            prev3 = three.TelomereLength3Prime;

            var five = ChromosomeAnalyzer.AnalyzeTelomeres("chr", Repeat(TeloRepeatRc, k) + Filler(60));
            five.TelomereLength5Prime.Should().Be(6 * k,
                because: $"the 5' end carries {k} contiguous CCCTAA units (the reverse complement of TTAGGG)");
            five.TelomereLength5Prime.Should().BeGreaterThan(prev5,
                because: "adding CCCTAA units to the 5' end strictly lengthens the measured telomere");
            prev5 = five.TelomereLength5Prime;
        }
    }

    #endregion

    #region INV — interior non-telomeric sequence does not change the telomeres

    [Test]
    [Description("INV: inserting non-telomeric sequence between the two terminal tracts leaves both telomere lengths and purities unchanged.")]
    public void AnalyzeTelomeres_InteriorFlank_PreservesTelomeres()
    {
        const int a = 30, b = 40;   // 5' and 3' tract repeat counts
        string fivePrime = Repeat(TeloRepeatRc, a);
        string threePrime = Repeat(TeloRepeat, b);

        var baseline = ChromosomeAnalyzer.AnalyzeTelomeres("chr", fivePrime + Filler(100) + threePrime);

        foreach (int extra in new[] { 200, 1000, 5000 })
        {
            var grown = ChromosomeAnalyzer.AnalyzeTelomeres("chr", fivePrime + Filler(100 + extra) + threePrime);

            grown.TelomereLength5Prime.Should().Be(baseline.TelomereLength5Prime,
                because: "the 5' terminal tract is untouched by interior insertions, so its measured length is unchanged");
            grown.TelomereLength3Prime.Should().Be(baseline.TelomereLength3Prime,
                because: "the 3' terminal tract is untouched by interior insertions, so its measured length is unchanged");
            grown.RepeatPurity5Prime.Should().BeApproximately(baseline.RepeatPurity5Prime, 1e-12,
                because: "the 5' tract content is unchanged, so its repeat purity is unchanged");
            grown.RepeatPurity3Prime.Should().BeApproximately(baseline.RepeatPurity3Prime, 1e-12,
                because: "the 3' tract content is unchanged, so its repeat purity is unchanged");
        }
    }

    #endregion

    #region SHIFT/strand-duality — reverse complement swaps the 5' and 3' telomeres

    [Test]
    [Description("SHIFT/strand-duality: reverse-complementing the chromosome maps the 3' TTAGGG tract onto a 5' CCCTAA tract, swapping the 5' and 3' telomere measurements.")]
    public void AnalyzeTelomeres_ReverseComplement_SwapsEnds()
    {
        const int a = 25, b = 60;   // asymmetric tracts so the swap is observable
        string seq = Repeat(TeloRepeatRc, a) + Filler(120) + Repeat(TeloRepeat, b);

        var forward = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
        var reversed = ChromosomeAnalyzer.AnalyzeTelomeres("chr", RevComp(seq));

        reversed.TelomereLength5Prime.Should().Be(forward.TelomereLength3Prime,
            because: "reverse-complementing turns the 3' TTAGGG tract into a 5' CCCTAA tract of equal length");
        reversed.TelomereLength3Prime.Should().Be(forward.TelomereLength5Prime,
            because: "reverse-complementing turns the 5' CCCTAA tract into a 3' TTAGGG tract of equal length");
        reversed.RepeatPurity5Prime.Should().BeApproximately(forward.RepeatPurity3Prime, 1e-12,
            because: "the swapped tracts carry the same residues, so their purities swap too");
        reversed.RepeatPurity3Prime.Should().BeApproximately(forward.RepeatPurity5Prime, 1e-12,
            because: "the swapped tracts carry the same residues, so their purities swap too");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CHROM-CENT-001 — centromere analysis (Chromosome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 49.
    //
    // API under test (ChromosomeAnalyzer.AnalyzeCentromere):
    //   Scans fixed-width windows for the centromeric signal score = repeatContent·(1−gcVariability)
    //   and reports the best-scoring region (Start/End) and its score as AlphaSatelliteContent.
    //   repeatContent = fraction of window positions whose 15-mer recurs; for windows below the
    //   1000-bp GC sub-window the gcVariability term is 0, so the score reduces to repeatContent.
    //
    // Relations (derived from the actual scoring, NOT from output):
    //   • MON (more satellite-repeat content ⇒ higher score): embedding a longer tandem array
    //          raises the window's repeat content, so AlphaSatelliteContent is non-decreasing in
    //          the array length, and a long array is detected while pure-random sequence is not.
    //          ── Reconciliation: the checklist says "more AT-rich → higher score", but the
    //          implementation scores TANDEM-REPEAT content (and GC uniformity), NOT base
    //          composition — a GC-only tandem array scores just as high as an AT-only one. The
    //          rigorous monotone driver (repeat content, the alpha-satellite array signal) is
    //          therefore what is tested, rather than a composition relation the code does not compute.
    //   • INV (non-centromeric flank append ⇒ same position): appending non-repetitive sequence
    //          AFTER the satellite array leaves the detected Start/End unchanged.
    //   • SHIFT (prepended flank shifts the position): prepending p bases of non-repetitive
    //          sequence (p on the scan grid) shifts both Start and End by exactly p.
    // ───────────────────────────────────────────────────────────────────────────

    #region Centromere helpers

    private const int CentWindow = 400;   // < 1000 ⇒ gcVariability term is 0 ⇒ score = repeatContent
    private static string Satellite(int length) => Repeat("AT", length / 2);

    #endregion

    #region MON — a longer satellite array gives a non-decreasing centromere score

    [Test]
    [Description("MON: embedding a longer tandem (alpha-satellite-like) array raises the repeat-content score, so AlphaSatelliteContent is non-decreasing in array length; a long array is detected while random sequence is not.")]
    public void AnalyzeCentromere_MoreRepeatContent_IncreasesScore()
    {
        double previous = double.NegativeInfinity;
        ChromosomeAnalyzer.CentromereResult? smallest = null, largest = null;

        foreach (int arrayLen in new[] { 0, 100, 300, 800 })
        {
            string seq = RandomDna(400) + Satellite(arrayLen) + RandomDna(400);
            var result = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize: CentWindow);

            result.AlphaSatelliteContent.Should().BeGreaterThanOrEqualTo(previous - 1e-12,
                because: "a longer tandem array raises the window's repeat content, so the centromere score does not decrease");
            previous = result.AlphaSatelliteContent;

            smallest ??= result;
            largest = result;
        }

        largest!.Value.AlphaSatelliteContent.Should().BeGreaterThan(smallest!.Value.AlphaSatelliteContent,
            because: "an 800-bp satellite array carries far more repeat content than pure-random sequence");
        largest.Value.Start.Should().NotBeNull(
            because: "a long satellite array exceeds the alpha-satellite content threshold and is detected as a centromere");
    }

    #endregion

    #region INV — appending non-centromeric sequence preserves the centromere position

    [Test]
    [Description("INV: appending non-repetitive sequence after the satellite array leaves the detected centromere Start/End unchanged.")]
    public void AnalyzeCentromere_AppendNonCentromericFlank_PreservesPosition()
    {
        string core = Satellite(800);
        var baseline = ChromosomeAnalyzer.AnalyzeCentromere("chr", core + RandomDna(800), windowSize: CentWindow);
        baseline.Start.Should().NotBeNull(because: "the 800-bp satellite array is detected as a centromere");

        foreach (int extra in new[] { 400, 1200, 3000 })
        {
            var grown = ChromosomeAnalyzer.AnalyzeCentromere("chr", core + RandomDna(800 + extra), windowSize: CentWindow);

            grown.Start.Should().Be(baseline.Start,
                because: "the satellite array stays at the front, so appending non-repetitive sequence cannot move the centromere start");
            grown.End.Should().Be(baseline.End,
                because: "right-extension stops at the array/flank boundary regardless of how much non-repetitive sequence follows");
        }
    }

    #endregion

    #region SHIFT — prepending non-centromeric sequence shifts the centromere position

    [Test]
    [Description("SHIFT: prepending p bases of non-repetitive sequence (p on the scan grid) shifts both the centromere Start and End by exactly p.")]
    public void AnalyzeCentromere_PrependFlank_ShiftsPosition()
    {
        string core = Satellite(800);
        var baseline = ChromosomeAnalyzer.AnalyzeCentromere("chr", core + RandomDna(800), windowSize: CentWindow);
        baseline.Start.Should().NotBeNull();

        foreach (int p in new[] { 100, 200, 400 })   // multiples of the scan step (windowSize / 4 = 100)
        {
            var shifted = ChromosomeAnalyzer.AnalyzeCentromere("chr", RandomDna(p) + core + RandomDna(800), windowSize: CentWindow);

            shifted.Start.Should().Be(baseline.Start + p,
                because: $"prepending {p} non-repetitive bases moves the satellite array — and the detected start — right by {p}");
            shifted.End.Should().Be(baseline.End + p,
                because: $"the whole detected region is translated right by the {p}-base prefix");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CHROM-KARYO-001 — karyotype analysis (Chromosome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 50.
    //
    // API under test (ChromosomeAnalyzer.AnalyzeKaryotype):
    //   Aggregates (Name, Length, IsSexChromosome) records into a Karyotype: total/autosome/
    //   sex counts, total and mean length, and per-base-name copy-number aneuploidy calls.
    //
    // Relations (derived from the aggregation, NOT from output):
    //   • COMP (N chromosomes ⇒ N entries): TotalChromosomes = N, autosome + sex counts = N,
    //          total genome size = Σ lengths, mean = total / N.
    //   • INV (order independence): every output is an order-insensitive aggregate (count, sum,
    //          group-by), so permuting the input chromosomes leaves the karyotype unchanged —
    //          the list fields compared as sets.
    // ───────────────────────────────────────────────────────────────────────────

    #region Karyotype helpers

    private static List<(string Name, long Length, bool IsSexChromosome)> DiploidSet() => new()
    {
        ("chr1_1", 248_000_000, false), ("chr1_2", 248_000_000, false),
        ("chr2_1", 242_000_000, false), ("chr2_2", 242_000_000, false),
        ("chr3_1", 198_000_000, false), ("chr3_2", 198_000_000, false),
        ("X", 156_000_000, true), ("Y", 57_000_000, true),
    };

    private static List<(string Name, long Length, bool IsSexChromosome)> AneuploidSet() => new()
    {
        ("chr1_1", 248_000_000, false), ("chr1_2", 248_000_000, false), ("chr1_3", 248_000_000, false), // trisomy
        ("chr2_1", 242_000_000, false),                                                                  // monosomy
        ("chr3_1", 198_000_000, false), ("chr3_2", 198_000_000, false),
        ("X", 156_000_000, true), ("X", 156_000_000, true),
    };

    private List<(string Name, long Length, bool IsSexChromosome)> Shuffle(
        List<(string, long, bool)> items)
    {
        var list = items.ToList();
        for (int i = list.Count - 1; i > 0; i--)
        {
            int j = Rng.Next(i + 1);
            (list[i], list[j]) = (list[j], list[i]);
        }
        return list;
    }

    #endregion

    #region COMP — N input chromosomes give a karyotype of N entries

    [Test]
    [Description("COMP: the karyotype accounts for exactly the N input chromosomes — total = N, autosome + sex counts = N, total size = Σ lengths, mean = total / N.")]
    public void AnalyzeKaryotype_AccountsForAllChromosomes()
    {
        foreach (var set in new[] { DiploidSet(), AneuploidSet() })
        {
            int n = set.Count;
            long expectedTotal = set.Sum(c => c.Length);
            int expectedSex = set.Count(c => c.IsSexChromosome);

            var k = ChromosomeAnalyzer.AnalyzeKaryotype(set);

            k.TotalChromosomes.Should().Be(n, because: "every input chromosome is counted exactly once");
            (k.AutosomeCount + k.SexChromosomes.Count).Should().Be(n,
                because: "each chromosome is classified as either an autosome or a sex chromosome");
            k.SexChromosomes.Count.Should().Be(expectedSex, because: "the sex-chromosome count matches the flagged inputs");
            k.TotalGenomeSize.Should().Be(expectedTotal, because: "the genome size is the sum of the chromosome lengths");
            k.MeanChromosomeLength.Should().BeApproximately(expectedTotal / (double)n, 1e-6,
                because: "the mean length is the total size divided by the chromosome count");
        }
    }

    #endregion

    #region INV — chromosome input order does not affect the karyotype

    [Test]
    [Description("INV: permuting the input chromosomes leaves every karyotype field unchanged, as all outputs are order-insensitive aggregates.")]
    public void AnalyzeKaryotype_InputOrder_DoesNotAffectResult()
    {
        foreach (var set in new[] { DiploidSet(), AneuploidSet() })
        {
            var baseline = ChromosomeAnalyzer.AnalyzeKaryotype(set);

            for (int trial = 0; trial < 5; trial++)
            {
                var shuffled = ChromosomeAnalyzer.AnalyzeKaryotype(Shuffle(set));

                shuffled.TotalChromosomes.Should().Be(baseline.TotalChromosomes);
                shuffled.AutosomeCount.Should().Be(baseline.AutosomeCount);
                shuffled.TotalGenomeSize.Should().Be(baseline.TotalGenomeSize);
                shuffled.MeanChromosomeLength.Should().BeApproximately(baseline.MeanChromosomeLength, 1e-6);
                shuffled.PloidyLevel.Should().Be(baseline.PloidyLevel);
                shuffled.HasAneuploidy.Should().Be(baseline.HasAneuploidy);
                shuffled.SexChromosomes.Should().BeEquivalentTo(baseline.SexChromosomes,
                    because: "the set of sex chromosomes is independent of input order");
                shuffled.Abnormalities.Should().BeEquivalentTo(baseline.Abnormalities,
                    because: "aneuploidy calls come from order-insensitive group-by counts");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CHROM-ANEU-001 — ploidy / copy-number from read depth (Chromosome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 51.
    //
    // API under test (ChromosomeAnalyzer.DetectPloidy):
    //   ploidy = round(2 · median(normalizedDepths) / expectedDiploidDepth), clamped to [1,8].
    //
    // Relations (derived from that formula, NOT from output):
    //   • MON (doubled depth ⇒ doubled CN): doubling every depth doubles the median and hence
    //          the ratio, so the rounded copy number doubles (within the [1,8] clamp); and CN is
    //          non-decreasing in depth.
    //   • INV (neighbouring region doesn't affect local CN): the estimate is a function of the
    //          MEDIAN, which is robust — appending a minority of contaminating depths from a
    //          neighbouring amplified/deleted region (fewer than the local count) does not move
    //          the median, so the local copy number is unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region MON — doubling all depths doubles the copy-number estimate

    [Test]
    [Description("MON: doubling every normalized depth doubles the median and the rounded copy number (within the [1,8] clamp).")]
    public void DetectPloidy_DoubledDepth_DoublesCopyNumber()
    {
        foreach (double c in new[] { 0.5, 1.0, 1.5, 2.0 })   // base ploidy 1,2,3,4 ⇒ doubled 2,4,6,8 (in range)
        {
            var baseDepths = Enumerable.Repeat(c, 25).ToList();
            var doubled = baseDepths.Select(d => 2 * d).ToList();

            int basePloidy = ChromosomeAnalyzer.DetectPloidy(baseDepths).PloidyLevel;
            int doubledPloidy = ChromosomeAnalyzer.DetectPloidy(doubled).PloidyLevel;

            doubledPloidy.Should().Be(2 * basePloidy,
                because: $"the median doubles from {c} to {2 * c}, so round(2·median) doubles the copy number");
        }
    }

    [Test]
    [Description("MON: the copy-number estimate is non-decreasing as the (uniform) read depth increases.")]
    public void DetectPloidy_HigherDepth_DoesNotDecreaseCopyNumber()
    {
        int previous = int.MinValue;
        foreach (double c in new[] { 0.4, 0.6, 1.0, 1.4, 1.8, 2.5, 3.5 })
        {
            int ploidy = ChromosomeAnalyzer.DetectPloidy(Enumerable.Repeat(c, 25).ToList()).PloidyLevel;
            ploidy.Should().BeGreaterThanOrEqualTo(previous,
                because: "a higher median depth maps to a higher-or-equal rounded copy number");
            previous = ploidy;
        }
    }

    #endregion

    #region INV — a minority of neighbouring-region depths does not change the local CN

    [Test]
    [Description("INV: appending a minority of contaminating depths from a neighbouring amplified or deleted region leaves the median-based local copy number unchanged.")]
    public void DetectPloidy_MinorityNeighbourContamination_PreservesLocalCopyNumber()
    {
        const int localCount = 25;
        const double localDepth = 1.0;   // diploid local region
        int localPloidy = ChromosomeAnalyzer.DetectPloidy(
            Enumerable.Repeat(localDepth, localCount).ToList()).PloidyLevel;

        foreach (double neighbourDepth in new[] { 4.0, 0.05 })   // a neighbouring amplification / deletion
        {
            foreach (int k in new[] { 1, 5, 12 })   // strictly fewer than the local count ⇒ median stays local
            {
                var contaminated = Enumerable.Repeat(localDepth, localCount)
                    .Concat(Enumerable.Repeat(neighbourDepth, k))
                    .ToList();

                ChromosomeAnalyzer.DetectPloidy(contaminated).PloidyLevel.Should().Be(localPloidy,
                    because: $"a minority ({k} < {localCount}) of neighbouring-region depths cannot move the median off the local depth, so the local copy number is unchanged");
            }
        }
    }

    #endregion
}
