using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CHROM-SYNT-001: Synteny Analysis
/// Covers FindSyntenyBlocks and DetectRearrangements methods.
/// </summary>
/// <remarks>
/// Sources:
/// - Wikipedia (Synteny): https://en.wikipedia.org/wiki/Synteny
/// - Wikipedia (Comparative genomics): https://en.wikipedia.org/wiki/Comparative_genomics
/// - Wikipedia (Chromosomal rearrangement): https://en.wikipedia.org/wiki/Chromosomal_rearrangement
/// - Wang et al. (2012) MCScanX, Nucleic Acids Res. 40(7):e49
/// - Goel et al. (2019) SyRI, Genome Biology 20:277
/// </remarks>
[TestFixture]
public class ChromosomeAnalyzer_Synteny_Tests
{
    #region FindSyntenyBlocks - MUST Tests

    /// <summary>
    /// M1: Core functionality - detect forward collinear blocks.
    /// Source: Wikipedia (Synteny) - collinearity definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_CollinearForward_ReturnsBlockWithPlusStrand()
    {
        // Arrange: 4 genes in same order in both genomes (forward collinearity)
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert: single block with exact properties
        // Hand-calculated: 4 collinear forward genes → 1 block, strand '+'
        // Coordinates span first.Start1=1000 to last.End1=8000
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "4 collinear genes should form exactly 1 block");
            Assert.That(blocks[0].Strand, Is.EqualTo('+'), "Forward collinearity → '+' strand");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(4), "All 4 genes should be in the block");
            Assert.That(blocks[0].Species1Chromosome, Is.EqualTo("chr1"));
            Assert.That(blocks[0].Species2Chromosome, Is.EqualTo("chrA"));
            Assert.That(blocks[0].Species1Start, Is.EqualTo(1000), "Block starts at first gene");
            Assert.That(blocks[0].Species1End, Is.EqualTo(8000), "Block ends at last gene");
            Assert.That(blocks[0].Species2Start, Is.EqualTo(1000), "Target starts at first gene");
            Assert.That(blocks[0].Species2End, Is.EqualTo(8000), "Target ends at last gene");
        });
    }

    /// <summary>
    /// M2: Core functionality - detect inverted (reverse) collinear blocks.
    /// Source: Wikipedia (Synteny) - inverted synteny blocks
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_CollinearReverse_ReturnsBlockWithMinusStrand()
    {
        // Arrange: genes in reverse order in target genome
        var orthologPairs = CreateReverseCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert: single block with exact properties
        // Hand-calculated: 4 genes with decreasing species2 positions → strand '-'
        // Species2: Min(8000,2000)=2000, Max(9000,3000)=9000
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "4 reverse-collinear genes → exactly 1 block");
            Assert.That(blocks[0].Strand, Is.EqualTo('-'), "Reverse collinearity → '-' strand");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(4), "All 4 genes should be in the block");
            Assert.That(blocks[0].Species1Start, Is.EqualTo(1000));
            Assert.That(blocks[0].Species1End, Is.EqualTo(8000));
            Assert.That(blocks[0].Species2Start, Is.EqualTo(2000), "Min of first/last species2 start");
            Assert.That(blocks[0].Species2End, Is.EqualTo(9000), "Max of first/last species2 end");
        });
    }

    /// <summary>
    /// M3: minGenes threshold enforcement.
    /// Source: Definition - blocks require minimum gene count
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_TooFewGenes_ReturnsEmpty()
    {
        // Arrange: only 1 gene (below minGenes=3)
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
        };

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3).ToList();

        // Assert
        Assert.That(blocks, Is.Empty, "Should return empty when genes < minGenes");
    }

    /// <summary>
    /// M4: Edge case - empty input returns empty.
    /// Source: Implementation
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_EmptyInput_ReturnsEmpty()
    {
        // Arrange
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3).ToList();

        // Assert
        Assert.That(blocks, Is.Empty, "Empty input should return empty result");
    }

    /// <summary>
    /// M5: Boundary - exactly minGenes should return a block.
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_ExactlyMinGenes_ReturnsBlock()
    {
        // Arrange: exactly 3 genes with minGenes=3
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
        };

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert: hand-calculated: exactly minGenes forward-collinear genes → 1 block
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(1), "Exactly minGenes should form 1 block");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(3), "GeneCount should equal minGenes");
            Assert.That(blocks[0].Strand, Is.EqualTo('+'), "Forward collinearity");
            Assert.That(blocks[0].Species1Start, Is.EqualTo(1000));
            Assert.That(blocks[0].Species1End, Is.EqualTo(6000));
        });
    }

    /// <summary>
    /// M6: Multiple chromosome pairs produce separate blocks.
    /// Source: Definition - blocks are grouped by chromosome pairs
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_MultipleChromosomePairs_SeparateBlocks()
    {
        // Arrange: genes on two different chromosome pairs
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            // First pair: chr1 -> chrA
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
            // Second pair: chr2 -> chrB
            ("chr2", 1000, 2000, "gene4", "chrB", 1000, 2000, "geneD"),
            ("chr2", 3000, 4000, "gene5", "chrB", 3000, 4000, "geneE"),
            ("chr2", 5000, 6000, "gene6", "chrB", 5000, 6000, "geneF"),
        };

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert: hand-calculated: 2 chromosome pairs → 2 separate blocks, 3 genes each
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(2), "2 chromosome pairs → 2 blocks");
            var block1 = blocks.First(b => b.Species1Chromosome == "chr1");
            var block2 = blocks.First(b => b.Species1Chromosome == "chr2");
            Assert.That(block1.Species2Chromosome, Is.EqualTo("chrA"));
            Assert.That(block1.GeneCount, Is.EqualTo(3));
            Assert.That(block1.Species1Start, Is.EqualTo(1000));
            Assert.That(block1.Species1End, Is.EqualTo(6000));
            Assert.That(block2.Species2Chromosome, Is.EqualTo("chrB"));
            Assert.That(block2.GeneCount, Is.EqualTo(3));
            Assert.That(block2.Species1Start, Is.EqualTo(1000));
            Assert.That(block2.Species1End, Is.EqualTo(6000));
        });
    }

    /// <summary>
    /// M16: maxGap parameter splits blocks when gap exceeds threshold.
    /// Source: MCScanX (Wang et al. 2012) - gap tolerance in synteny detection
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_GapExceedsMaxGap_SplitsIntoSeparateBlocks()
    {
        // Arrange: 6 collinear genes, 3 close + 3MB gap + 3 close
        var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
            // Gap: 3000000 - 6000 = 2994000 bp > maxGap=2 (2MB)
            ("chr1", 3000000, 4000000, "gene4", "chrA", 3000000, 4000000, "geneD"),
            ("chr1", 5000000, 6000000, "gene5", "chrA", 5000000, 6000000, "geneE"),
            ("chr1", 7000000, 8000000, "gene6", "chrA", 7000000, 8000000, "geneF"),
        };

        // Act: maxGap=2 means 2MB threshold; gap between genes 3-4 is ~3MB -> split
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 2).ToList();

        // Assert: hand-calculated: 2 blocks of 3 genes each
        // Block 1: genes 1-3, chr1:1000-6000 -> chrA:1000-6000
        // Block 2: genes 4-6, chr1:3000000-8000000 -> chrA:3000000-8000000
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.EqualTo(2), "Gap > maxGap should split into 2 blocks");
            Assert.That(blocks[0].GeneCount, Is.EqualTo(3));
            Assert.That(blocks[0].Species1Start, Is.EqualTo(1000));
            Assert.That(blocks[0].Species1End, Is.EqualTo(6000));
            Assert.That(blocks[1].GeneCount, Is.EqualTo(3));
            Assert.That(blocks[1].Species1Start, Is.EqualTo(3000000));
            Assert.That(blocks[1].Species1End, Is.EqualTo(8000000));
        });
    }

    #endregion

    #region DetectRearrangements - MUST Tests

    /// <summary>
    /// M9: Core functionality - detect inversions.
    /// Source: Wikipedia (Chromosomal rearrangement) - inversion definition
    /// </summary>
    [Test]
    public void DetectRearrangements_Inversion_DetectsInversion()
    {
        // Arrange: adjacent blocks with strand change (same target chromosome)
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrA", 60000, 100000, '-', 8, 0.93)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert: hand-calculated from code logic
        // Same chr2 + different strand → Inversion
        // Position1 = current.Species1End = 50000
        // Position2 = next.Species1Start = 60000
        // Size = 60000 - 50000 = 10000
        Assert.Multiple(() =>
        {
            Assert.That(rearrangements, Has.Count.EqualTo(1), "Exactly 1 inversion");
            Assert.That(rearrangements[0].Type, Is.EqualTo("Inversion"));
            Assert.That(rearrangements[0].Chromosome1, Is.EqualTo("chr1"));
            Assert.That(rearrangements[0].Position1, Is.EqualTo(50000));
            Assert.That(rearrangements[0].Position2, Is.EqualTo(60000));
            Assert.That(rearrangements[0].Size, Is.EqualTo(10000));
        });
    }

    /// <summary>
    /// M10: Core functionality - detect translocations.
    /// Source: Wikipedia (Chromosomal rearrangement) - translocation definition
    /// </summary>
    [Test]
    public void DetectRearrangements_Translocation_DetectsTranslocation()
    {
        // Arrange: adjacent blocks mapping to different chromosomes
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrB", 1000, 40000, '+', 8, 0.93)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert: hand-calculated from code logic
        // Different chr2 → Translocation
        // Position1 = current.Species1End = 50000
        // Chromosome2 = next.Species2Chromosome = "chrB"
        // Position2 = next.Species2Start = 1000
        Assert.Multiple(() =>
        {
            Assert.That(rearrangements, Has.Count.EqualTo(1), "Exactly 1 translocation");
            Assert.That(rearrangements[0].Type, Is.EqualTo("Translocation"));
            Assert.That(rearrangements[0].Chromosome1, Is.EqualTo("chr1"));
            Assert.That(rearrangements[0].Position1, Is.EqualTo(50000));
            Assert.That(rearrangements[0].Chromosome2, Is.EqualTo("chrB"));
            Assert.That(rearrangements[0].Position2, Is.EqualTo(1000));
            Assert.That(rearrangements[0].Size, Is.Null, "Translocation has no size");
        });
    }

    /// <summary>
    /// M11: Edge case - empty input returns empty.
    /// Source: Implementation
    /// </summary>
    [Test]
    public void DetectRearrangements_EmptyInput_ReturnsEmpty()
    {
        // Arrange
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>();

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.That(rearrangements, Is.Empty, "Empty input should return empty result");
    }

    /// <summary>
    /// M12: Edge case - single block has no adjacent pairs.
    /// Source: Definition
    /// </summary>
    [Test]
    public void DetectRearrangements_SingleBlock_ReturnsEmpty()
    {
        // Arrange
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.That(rearrangements, Is.Empty,
            "Single block has no adjacent pairs for comparison");
    }

    /// <summary>
    /// M13: Collinear genome with no rearrangements.
    /// Source: Definition
    /// </summary>
    [Test]
    public void DetectRearrangements_NoRearrangements_ReturnsEmpty()
    {
        // Arrange: perfectly collinear blocks (same chromosome, same strand)
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrA", 60000, 100000, '+', 8, 0.93),
            new("chr1", 110000, 150000, "chrA", 110000, 150000, '+', 6, 0.91)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.That(rearrangements, Is.Empty,
            "Perfectly collinear blocks should produce no rearrangements");
    }

    #endregion

    #region SHOULD Tests - Invariants

    /// <summary>
    /// S1: Strand invariant - must be '+' or '-'.
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_StrandIsValidChar()
    {
        // Arrange
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.That(blocks.All(b => b.Strand == '+' || b.Strand == '-'), Is.True,
            "All blocks must have strand '+' or '-'");
    }

    /// <summary>
    /// S2: SequenceIdentity is NaN when not computable from coordinate-only input.
    /// Source: MCScanX (Wang et al. 2012) - identity requires BLAST alignment data
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_SequenceIdentityIsNaN()
    {
        // Arrange
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert: Identity not computable from coordinate data alone (no sequences provided)
        Assert.That(blocks.All(b => double.IsNaN(b.SequenceIdentity)), Is.True,
            "SequenceIdentity should be NaN when not computable from coordinate-only input");
    }

    /// <summary>
    /// S3: Coordinates invariant - Start <= End.
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_CoordinatesValid()
    {
        // Arrange
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks.All(b => b.Species1Start <= b.Species1End), Is.True,
                "Species1Start must be <= Species1End");
            Assert.That(blocks.All(b => b.Species2Start <= b.Species2End), Is.True,
                "Species2Start must be <= Species2End");
        });
    }

    /// <summary>
    /// S4: Rearrangement Type invariant - recognized values only.
    /// Source: Definition
    /// </summary>
    [Test]
    public void DetectRearrangements_TypeIsRecognizedValue()
    {
        // Arrange: create blocks that will produce rearrangements
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrA", 60000, 100000, '-', 8, 0.93),
            new("chr1", 110000, 150000, "chrB", 1000, 40000, '+', 6, 0.91)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        var validTypes = new[] { "Inversion", "Translocation", "Deletion", "Duplication" };
        Assert.That(rearrangements.All(r => validTypes.Contains(r.Type)), Is.True,
            "All rearrangement types must be recognized values");
    }

    /// <summary>
    /// S5: Position1 invariant - always set for rearrangements.
    /// Source: Definition
    /// </summary>
    [Test]
    public void DetectRearrangements_Position1AlwaysSet()
    {
        // Arrange
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 60000, 100000, "chrB", 1000, 40000, '+', 8, 0.93)
        };

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.That(rearrangements.All(r => r.Position1 > 0), Is.True,
            "Position1 must be set for all rearrangements");
    }

    /// <summary>
    /// M14: Core functionality - detect deletions.
    /// Source: Wikipedia (Chromosomal rearrangement) - deletion = segment is removed
    /// </summary>
    [Test]
    public void DetectRearrangements_Deletion_DetectsDeletion()
    {
        // Arrange: adjacent blocks on same chromosome pair, same strand,
        // with asymmetric gap (large gap in species 1, small gap in species 2)
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 150000, 200000, "chrA", 55000, 100000, '+', 8, 0.93)
        };
        // Hand-calculated:
        // gap1 = 150000 - 50000 = 100000
        // gap2 = 55000 - 50000 = 5000 (strand='+')
        // gap1 > gap2*2: 100000 > 10000 → true → Deletion
        // Size = gap1 - gap2 = 100000 - 5000 = 95000

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(rearrangements, Has.Count.EqualTo(1), "Exactly 1 deletion");
            Assert.That(rearrangements[0].Type, Is.EqualTo("Deletion"));
            Assert.That(rearrangements[0].Chromosome1, Is.EqualTo("chr1"));
            Assert.That(rearrangements[0].Position1, Is.EqualTo(50000));
            Assert.That(rearrangements[0].Position2, Is.EqualTo(150000));
            Assert.That(rearrangements[0].Size, Is.EqualTo(95000),
                "Deletion size = gap1 - gap2 = 100000 - 5000");
        });
    }

    /// <summary>
    /// M15: Core functionality - detect duplications.
    /// Source: Wikipedia (Chromosomal rearrangement) - duplication = segment is copied
    /// </summary>
    [Test]
    public void DetectRearrangements_Duplication_DetectsDuplication()
    {
        // Arrange: overlapping species 1 coordinates mapping to different species 2 locations
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            new("chr1", 1000, 50000, "chrA", 1000, 50000, '+', 10, 0.95),
            new("chr1", 20000, 70000, "chrA", 200000, 250000, '+', 8, 0.93)
        };
        // Hand-calculated:
        // Overlap in species 1: max(1000,20000)=20000 to min(50000,70000)=50000
        // Different species 2 locations → Duplication
        // Size = 50000 - 20000 = 30000

        // Act
        var rearrangements = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(rearrangements, Has.Count.EqualTo(1), "Exactly 1 duplication");
            Assert.That(rearrangements[0].Type, Is.EqualTo("Duplication"));
            Assert.That(rearrangements[0].Chromosome1, Is.EqualTo("chr1"));
            Assert.That(rearrangements[0].Position1, Is.EqualTo(20000), "Overlap start");
            Assert.That(rearrangements[0].Chromosome2, Is.EqualTo("chrA"));
            Assert.That(rearrangements[0].Position2, Is.EqualTo(200000), "Target block start");
            Assert.That(rearrangements[0].Size, Is.EqualTo(30000),
                "Duplication size = overlap extent = 50000 - 20000");
        });
    }

    #endregion

    #region Helper Methods

    private static List<(string, int, int, string, string, int, int, string)> CreateForwardCollinearPairs()
    {
        return new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
            ("chr1", 7000, 8000, "gene4", "chrA", 7000, 8000, "geneD"),
        };
    }

    private static List<(string, int, int, string, string, int, int, string)> CreateReverseCollinearPairs()
    {
        return new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 8000, 9000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 6000, 7000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 4000, 5000, "geneC"),
            ("chr1", 7000, 8000, "gene4", "chrA", 2000, 3000, "geneD"),
        };
    }

    #endregion
}
