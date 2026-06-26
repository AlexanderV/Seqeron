using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Chromosome area (synteny block detection).
///
/// Algebraic testing pins the role-swap symmetry of synteny (comparing genome A
/// to B versus B to A yields the same blocks with the two species' coordinates
/// exchanged) and the self-comparison identity (a genome compared to itself is
/// fully syntenic — one collinear block spanning every ortholog).
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 52.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Chromosome")]
public class ChromosomeAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-SYNT-001 — Synteny block detection (Chromosome)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 52.
    //
    // Model: FindSyntenyBlocks groups collinear ortholog pairs into conserved
    //        blocks. The two genomes enter symmetrically: swapping their roles in
    //        each ortholog pair produces the same blocks with Species1 and Species2
    //        coordinates exchanged. A genome compared to itself maps every gene to
    //        its own location, so a forward-collinear run forms one block over all.
    //   — docs/algorithms/Chromosome_Analysis; ChromosomeAnalyzer.FindSyntenyBlocks.
    //
    // Laws under test (checklist row 52):
    //   • ID   — self-comparison → full synteny (one '+' block spanning all genes).
    //   • COMM — blocks(A, B) = role-swapped blocks(B, A): each block's Species1
    //            and Species2 coordinates exchange when the inputs are swapped.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Builds the 8-field ortholog-pair tuple FindSyntenyBlocks consumes
    /// (Chr1, Start1, End1, Gene1, Chr2, Start2, End2, Gene2).
    /// </summary>
    private static (string, int, int, string, string, int, int, string) Pair(
        string chr1, int s1, int e1, string g1, string chr2, int s2, int e2, string g2)
        => (chr1, s1, e1, g1, chr2, s2, e2, g2);

    /// <summary>
    /// ID: a genome compared to itself — every gene maps to its own coordinates —
    /// yields a single forward block covering all genes (full synteny).
    /// </summary>
    [Test]
    public void Synteny_Identity_SelfComparisonIsFullSynteny()
    {
        var self = new[]
        {
            Pair("chr1", 0, 100, "g1", "chr1", 0, 100, "g1"),
            Pair("chr1", 200, 300, "g2", "chr1", 200, 300, "g2"),
            Pair("chr1", 400, 500, "g3", "chr1", 400, 500, "g3"),
            Pair("chr1", 600, 700, "g4", "chr1", 600, 700, "g4"),
            Pair("chr1", 800, 900, "g5", "chr1", 800, 900, "g5"),
        };

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(self, minGenes: 3).ToList();

        blocks.Should().HaveCount(1);
        var b = blocks[0];
        b.GeneCount.Should().Be(5);
        b.Strand.Should().Be('+');
        b.Species1Start.Should().Be(0);
        b.Species1End.Should().Be(900);
        b.Species2Start.Should().Be(0);
        b.Species2End.Should().Be(900);
    }

    /// <summary>
    /// COMM: comparing A→B and B→A on the same collinear ortholog set yields the
    /// same block with the two species' coordinates exchanged.
    /// </summary>
    [Test]
    public void Synteny_Commutative_RoleSwapExchangesCoordinates()
    {
        var ab = new[]
        {
            Pair("chrA", 0, 100, "a1", "chrB", 1000, 1100, "b1"),
            Pair("chrA", 200, 300, "a2", "chrB", 1200, 1300, "b2"),
            Pair("chrA", 400, 500, "a3", "chrB", 1400, 1500, "b3"),
            Pair("chrA", 600, 700, "a4", "chrB", 1600, 1700, "b4"),
        };
        // Swap the two genomes' roles in every pair.
        var ba = ab.Select(p => Pair(p.Item5, p.Item6, p.Item7, p.Item8,
                                     p.Item1, p.Item2, p.Item3, p.Item4)).ToArray();

        var blocksAb = ChromosomeAnalyzer.FindSyntenyBlocks(ab, minGenes: 3).ToList();
        var blocksBa = ChromosomeAnalyzer.FindSyntenyBlocks(ba, minGenes: 3).ToList();

        blocksAb.Should().HaveCount(1);
        blocksBa.Should().HaveCount(1);

        var x = blocksAb[0];
        var y = blocksBa[0];

        // Species roles are exchanged between the two directions.
        y.Species1Chromosome.Should().Be(x.Species2Chromosome);
        y.Species1Start.Should().Be(x.Species2Start);
        y.Species1End.Should().Be(x.Species2End);
        y.Species2Chromosome.Should().Be(x.Species1Chromosome);
        y.Species2Start.Should().Be(x.Species1Start);
        y.Species2End.Should().Be(x.Species1End);
        y.GeneCount.Should().Be(x.GeneCount);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-ALPHASAT-001 — Alpha-satellite detection (Chromosome), row 257.
    //
    // Model: a centromeric alpha-satellite array is a tandem array of ~171-bp AT-rich monomers,
    //        many carrying the 17-bp CENP-B box (Willard 1985; Masumoto et al. 1989). Detection
    //        requires all three signals: 171-bp periodicity, AT-richness, and the CENP-B motif.
    //   — ChromosomeAnalyzer.DetectAlphaSatellite; TestSpec CHROM-ALPHASAT-001.
    //
    // Laws (row 257): ID — a 171-bp-periodic + AT-rich + CENP-B-box tandem array is detected
    //                 (BestPeriod = 171, one CENP-B box per monomer).  IDEMP — deterministic.
    // ═══════════════════════════════════════════════════════════════════════

    // 171-bp AT-rich monomer carrying one CENP-B box: 77 'A' + box(17) + 77 'T'.
    // box = TTTCGTTGGAAGCGGGA is a Y→T, R→G instance of the consensus YTTCGTTGGAARCGGGA.
    private const string AlphaSatMonomer =
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        + "TTTCGTTGGAAGCGGGA"
        + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

    private static string AlphaSatArray(int copies) => string.Concat(Enumerable.Repeat(AlphaSatMonomer, copies));

    [Test]
    public void AlphaSatellite_Identity_PeriodicAtRichCenpBArrayIsDetected()
    {
        AlphaSatMonomer.Length.Should().Be(171, "fixture sanity: monomer is one 171-bp alpha-satellite unit");
        const int copies = 10;

        var result = ChromosomeAnalyzer.DetectAlphaSatellite(AlphaSatArray(copies));

        result.IsAlphaSatellite.Should().BeTrue("171-bp-periodic + AT-rich + CENP-B array is alpha-satellite");
        result.BestPeriod.Should().Be(ChromosomeAnalyzer.AlphaSatelliteMonomerLength); // 171
        result.AtContent.Should().BeGreaterThan(0.50);
        result.CenpBBoxCount.Should().Be(copies, "exactly one CENP-B box per monomer");
    }

    [Test]
    public void AlphaSatellite_Idempotent_Deterministic()
    {
        var a = ChromosomeAnalyzer.DetectAlphaSatellite(AlphaSatArray(10));
        var b = ChromosomeAnalyzer.DetectAlphaSatellite(AlphaSatArray(10));
        a.IsAlphaSatellite.Should().Be(b.IsAlphaSatellite);
        a.BestPeriod.Should().Be(b.BestPeriod);
        a.PeriodicityScore.Should().Be(b.PeriodicityScore);
        a.CenpBBoxCount.Should().Be(b.CenpBBoxCount);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CHROM-HOR-001 — Higher-order repeat detection (Chromosome), row 258.
    //
    // Model: an alpha-satellite higher-order repeat (HOR) is a block of k distinct ~171-bp
    //        monomers (50–70% mutually identical) that is itself tandemly repeated with high
    //        (≥95%) inter-copy identity. The HOR period is k monomers; the unit length is k×171 bp.
    //   — ChromosomeAnalyzer.DetectHigherOrderRepeat; McNulty & Sullivan (2018); TestSpec CHROM-HOR-001.
    //
    // Laws (row 258): ID — a k-monomer HOR unit repeated n× → period k, unit length k×171.
    //                 IDEMP — deterministic.
    // ═══════════════════════════════════════════════════════════════════════

    private const int HorMonomerLength = 171; // ChromosomeAnalyzer.AlphaSatelliteMonomerLength

    // A fixed high-complexity 171-bp background (deterministic LCG, no homopolymer runs).
    private static readonly char[] HorBackground = BuildHorBackground();

    private static char[] BuildHorBackground()
    {
        var bases = new[] { 'A', 'C', 'G', 'T' };
        var chars = new char[HorMonomerLength];
        long state = 123456789L;
        char prev = '\0';
        for (int i = 0; i < HorMonomerLength; i++)
        {
            state = (1103515245L * state + 12345L) & 0x7fffffff;
            char b = bases[(int)(state % 4)];
            if (b == prev) b = bases[(int)((state + 1) % 4)];
            chars[i] = b;
            prev = b;
        }
        return chars;
    }

    private static char HorNextBase(char b) => b switch { 'A' => 'C', 'C' => 'G', 'G' => 'T', _ => 'A' };

    // A monomer = background with the given scattered positions overwritten by a different base.
    private static string HorMonomer(int[] substituted)
    {
        var chars = (char[])HorBackground.Clone();
        foreach (int p in substituted) chars[p] = HorNextBase(HorBackground[p]);
        return new string(chars);
    }

    // Positions spaced 2 apart: three disjoint sets ⇒ pairwise ≈57.9% identity (50–70% intra-HOR band).
    private static int[] HorScattered(int start, int count) =>
        Enumerable.Range(0, count).Select(k => start + 2 * k).ToArray();

    [Test]
    public void HigherOrderRepeat_Identity_KMonomerUnitHasPeriodKTimes171()
    {
        // 3 distinct monomers (A,B,C) tandemly repeated 5× with EXACT inter-HOR copies.
        string a = HorMonomer(HorScattered(0, 36));
        string b = HorMonomer(HorScattered(1, 36));   // odd positions, disjoint from A
        string c = HorMonomer(HorScattered(72, 36));  // disjoint from A and B
        string unit = a + b + c;                       // k = 3 monomers
        string array = string.Concat(Enumerable.Repeat(unit, 5)); // 15 monomers, 2565 bp

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        result.HasHigherOrderStructure.Should().BeTrue("a 3-monomer unit repeated 5× is a HOR");
        result.MonomersPerUnit.Should().Be(3, "the HOR period is the k = 3 monomers per unit");
        result.HorUnitLengthBp.Should().Be(3 * HorMonomerLength, "unit length = k × 171 = 513 bp");
        result.HorCopyNumber.Should().Be(5, "15 monomers / period 3 = 5 HOR copies");
        result.MonomerCount.Should().Be(15);
        result.MeanInterHorIdentity.Should().BeGreaterThan(result.MeanIntraHorIdentity,
            "HOR hallmark: inter-copy identity ≫ intra-unit identity");
    }

    [Test]
    public void HigherOrderRepeat_Idempotent_Deterministic()
    {
        string unit = HorMonomer(HorScattered(0, 36)) + HorMonomer(HorScattered(1, 36)) + HorMonomer(HorScattered(72, 36));
        string array = string.Concat(Enumerable.Repeat(unit, 5));

        var a = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);
        var b = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);
        a.MonomersPerUnit.Should().Be(b.MonomersPerUnit);
        a.HorUnitLengthBp.Should().Be(b.HorUnitLengthBp);
        a.HorCopyNumber.Should().Be(b.HorCopyNumber);
        a.MeanInterHorIdentity.Should().Be(b.MeanInterHorIdentity);
    }
}
