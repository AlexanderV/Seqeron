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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: CHROM-SYNT-001 — synteny block identification (Chromosome).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 52.
    //
    // API under test (ChromosomeAnalyzer.FindSyntenyBlocks):
    //   Groups ortholog pairs by (chr1, chr2), sorts by genome-1 position, and emits collinear
    //   runs of ≥ minGenes as SyntenyBlocks (Species1/Species2 spans, strand, gene count).
    //
    // Relations (derived from that definition, NOT from output):
    //   • SYM (species swap mirrors the blocks): swapping the two species in every ortholog pair
    //          produces the same forward-collinear blocks with the Species1 and Species2
    //          coordinate fields exchanged and the gene count and strand preserved.
    //   • INV (non-syntenic insert ⇒ same blocks): adding ortholog pairs on a DIFFERENT
    //          chromosome pair with fewer than minGenes forms no block and lies in a separate
    //          group, so the original synteny blocks are unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region Synteny helpers

    // (Chr1, Start1, End1, Gene1, Chr2, Start2, End2, Gene2)
    private static List<(string, int, int, string, string, int, int, string)> CollinearOrthologs(
        string chr1, string chr2, int n)
    {
        var pairs = new List<(string, int, int, string, string, int, int, string)>();
        for (int i = 0; i < n; i++)
            pairs.Add((chr1, i * 1000, i * 1000 + 500, $"{chr1}_g{i}",
                       chr2, i * 1000, i * 1000 + 500, $"{chr2}_g{i}"));
        return pairs;
    }

    private static (string, int, int, string, int, int, char, int) Key(ChromosomeAnalyzer.SyntenyBlock b) =>
        (b.Species1Chromosome, b.Species1Start, b.Species1End,
         b.Species2Chromosome, b.Species2Start, b.Species2End, b.Strand, b.GeneCount);

    #endregion

    #region SYM — swapping the two species mirrors the blocks

    [Test]
    [Description("SYM: swapping the two species in every ortholog pair yields the same forward-collinear block with the Species1/Species2 coordinate fields exchanged and the gene count/strand preserved.")]
    public void FindSyntenyBlocks_SpeciesSwap_MirrorsBlocks()
    {
        var ab = CollinearOrthologs("chrA", "chrB", 5);
        var ba = ab.Select(p => (p.Item5, p.Item6, p.Item7, p.Item8, p.Item1, p.Item2, p.Item3, p.Item4)).ToList();

        var blocksAb = ChromosomeAnalyzer.FindSyntenyBlocks(ab).ToList();
        var blocksBa = ChromosomeAnalyzer.FindSyntenyBlocks(ba).ToList();

        blocksAb.Should().HaveCount(1, because: "the five collinear orthologs form a single synteny block");
        blocksBa.Should().HaveCount(blocksAb.Count, because: "swapping the species cannot change how many collinear blocks exist");

        var x = blocksAb[0];
        var y = blocksBa[0];
        y.Species1Chromosome.Should().Be(x.Species2Chromosome, because: "the second species becomes the first under the swap");
        y.Species2Chromosome.Should().Be(x.Species1Chromosome, because: "the first species becomes the second under the swap");
        y.Species1Start.Should().Be(x.Species2Start);
        y.Species1End.Should().Be(x.Species2End);
        y.Species2Start.Should().Be(x.Species1Start);
        y.Species2End.Should().Be(x.Species1End);
        y.GeneCount.Should().Be(x.GeneCount, because: "the same orthologs are grouped, just with the two genomes exchanged");
        y.Strand.Should().Be(x.Strand, because: "co-increasing coordinates remain a forward (+) block after the swap");
    }

    #endregion

    #region INV — a non-syntenic insert on a separate chromosome pair leaves blocks unchanged

    [Test]
    [Description("INV: adding ortholog pairs on a different chromosome pair with fewer than minGenes forms no block and lies in a separate group, so the original synteny blocks are unchanged.")]
    public void FindSyntenyBlocks_NonSyntenicInsert_PreservesBlocks()
    {
        var baseline = CollinearOrthologs("chrA", "chrB", 5);
        var baseBlocks = ChromosomeAnalyzer.FindSyntenyBlocks(baseline).Select(Key).ToHashSet();

        // Two orthologs on an unrelated chromosome pair — below minGenes (3), so no block forms.
        var withInsert = baseline
            .Concat(new[]
            {
                ("chrC", 5000, 5500, "chrC_g0", "chrD", 80000, 80500, "chrD_g0"),
                ("chrC", 9000, 9500, "chrC_g1", "chrD", 12000, 12500, "chrD_g1"),
            })
            .ToList();

        var withInsertBlocks = ChromosomeAnalyzer.FindSyntenyBlocks(withInsert).Select(Key).ToHashSet();

        withInsertBlocks.Should().BeEquivalentTo(baseBlocks,
            because: "the inserted orthologs form their own sub-minGenes group, contributing no block and leaving the chrA/chrB synteny untouched");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-ALPHASAT-001 / CHROM-HOR-001 — alpha-satellite & HOR detection (Chromosome)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Willard 1985; Waye & Willard 1987; Masumoto et al. 1989 CENP-B box;
    //   McNulty & Sullivan 2018; Miga et al. 2014/T2T; tests/TestSpecs/CHROM-ALPHASAT/HOR-001.md):
    //   Alpha-satellite DNA is a tandem array of ~171-bp AT-rich monomers, many carrying a 17-bp
    //   CENP-B box (consensus YTTCGTTGGAARCGGGA). A higher-order repeat (HOR) is a block of k DISTINCT
    //   monomers tandemly repeated, so copies of the unit (inter-HOR) are far more identical than the
    //   monomers within a unit (intra-HOR). Four metamorphic relations (checklist rows 257–258):
    //
    //   • INV (a tandem 171-bp array is detected as alpha-satellite): an AT-rich 171-bp monomer tandem
    //     is called alpha-satellite with a detected period ≈ 171 bp.
    //   • MON (more CENP-B boxes → stronger call): adding CENP-B boxes to more monomers raises the
    //     reported CENP-B box count while the array stays alpha-satellite.
    //   • MON (more HOR copies → stronger periodicity): repeating a k-monomer HOR unit more times
    //     raises the detected HOR copy number (and keeps the higher-order structure).
    //   • INV (a pure monomeric array has no HOR): an array of one identical monomer has no
    //     higher-order structure (period 1).
    //
    // API under test: ChromosomeAnalyzer.DetectAlphaSatellite / .DetectHigherOrderRepeat.

    #region CHROM-ALPHASAT-001 / CHROM-HOR-001 — alpha-satellite & HOR

    private const int AlphaMonomerLen = ChromosomeAnalyzer.AlphaSatelliteMonomerLength; // 171
    private const string CenpBBox = "CTTCGTTGGAAACGGGA"; // a valid CENP-B box (Y=C, R=A), 17 bp

    private static string AtRichFiller(int length, int seed)
    {
        var rng = new Random(seed);
        char[] atRich = { 'A', 'A', 'T', 'T', 'A', 'T', 'G', 'C', 'A', 'T' }; // 70% AT
        var chars = new char[length];
        for (int i = 0; i < length; i++) chars[i] = atRich[rng.Next(atRich.Length)];
        return new string(chars);
    }

    [Test]
    [Description("INV: an AT-rich 171-bp monomer tandem is detected as alpha-satellite with a monomer period ≈ 171 bp.")]
    public void AlphaSat_Tandem171bpArray_DetectedAsAlphaSatellite()
    {
        string monomer = CenpBBox + AtRichFiller(AlphaMonomerLen - CenpBBox.Length, seed: 7);
        string array = string.Concat(Enumerable.Repeat(monomer, 5));

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

        result.IsAlphaSatellite.Should().BeTrue(because: "a tandem AT-rich 171-bp monomer array is alpha-satellite");
        result.BestPeriod.Should().BeInRange(AlphaMonomerLen - 5, AlphaMonomerLen + 5,
            because: "the detected monomer period is ≈ 171 bp");
    }

    [Test]
    [Description("MON: adding a CENP-B box to more monomers raises the reported CENP-B box count while the AT-rich array stays alpha-satellite.")]
    public void AlphaSat_MoreCenpBBoxes_StrengthensCall()
    {
        string tail = AtRichFiller(AlphaMonomerLen - CenpBBox.Length, seed: 11); // shared 154-bp tail
        string plainHead = AtRichFiller(CenpBBox.Length, seed: 23);              // AT-rich, non-box head
        string boxMonomer = CenpBBox + tail;
        string plainMonomer = plainHead + tail;
        const int total = 6;

        int previous = -1;
        bool sawStrictIncrease = false;
        for (int boxes = 0; boxes <= total; boxes++)
        {
            string array = string.Concat(Enumerable.Repeat(boxMonomer, boxes))
                           + string.Concat(Enumerable.Repeat(plainMonomer, total - boxes));
            var result = ChromosomeAnalyzer.DetectAlphaSatellite(array);

            result.IsAlphaSatellite.Should().BeTrue(because: $"the AT-rich 171-bp array stays alpha-satellite with {boxes} CENP-B boxes");
            result.CenpBBoxCount.Should().BeGreaterThanOrEqualTo(previous,
                because: $"adding a CENP-B box to one more monomer cannot lower the box count ({boxes} boxes)");
            if (result.CenpBBoxCount > previous) sawStrictIncrease = true;
            previous = result.CenpBBoxCount;
        }

        sawStrictIncrease.Should().BeTrue(because: "more CENP-B boxes genuinely raise the count — the relation is non-vacuous");
    }

    /// <summary>Builds a HOR of k distinct random 171-bp monomers, the unit repeated c times.</summary>
    private static string BuildHor(int seed, int k, int copies)
    {
        var rng = new Random(seed);
        char[] bases = { 'A', 'C', 'G', 'T' };
        var monomers = new string[k];
        for (int j = 0; j < k; j++)
        {
            var m = new char[AlphaMonomerLen];
            for (int i = 0; i < AlphaMonomerLen; i++) m[i] = bases[rng.Next(4)];
            monomers[j] = new string(m);
        }

        return string.Concat(Enumerable.Repeat(string.Concat(monomers), copies));
    }

    [Test]
    [Description("MON: repeating a 2-monomer HOR unit more times raises the detected HOR copy number and keeps the higher-order structure.")]
    public void Hor_MoreCopies_StrengthensPeriodicity()
    {
        int previous = -1;
        foreach (int copies in new[] { 2, 3, 4, 5 })
        {
            string array = BuildHor(seed: 2026, k: 2, copies: copies);
            var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array, AlphaMonomerLen);

            result.HasHigherOrderStructure.Should().BeTrue(because: $"a 2-monomer unit repeated {copies}× is a HOR");
            result.MonomersPerUnit.Should().Be(2, because: "the HOR period is the 2-monomer unit");
            result.HorCopyNumber.Should().BeGreaterThan(previous,
                because: $"more unit copies ({copies}) raise the detected HOR copy number");
            result.HorUnitLengthBp.Should().Be(2 * AlphaMonomerLen, because: "the HOR unit is 2 × 171 bp");
            previous = result.HorCopyNumber;
        }
    }

    [Test]
    [Description("INV: an array of one identical 171-bp monomer has no higher-order structure (period 1).")]
    public void Hor_PureMonomericArray_HasNoHor()
    {
        foreach (int copies in new[] { 4, 6, 10 })
        {
            string monomer = BuildHor(seed: 99, k: 1, copies: 1); // one 171-bp monomer
            string array = string.Concat(Enumerable.Repeat(monomer, copies));

            var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array, AlphaMonomerLen);

            result.HasHigherOrderStructure.Should().BeFalse(
                because: $"{copies} identical monomers form a monomeric (period-1) array, not a HOR");
            result.MonomersPerUnit.Should().Be(1, because: "a purely monomeric array has a single-monomer period");
        }
    }

    #endregion
}
