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
}
