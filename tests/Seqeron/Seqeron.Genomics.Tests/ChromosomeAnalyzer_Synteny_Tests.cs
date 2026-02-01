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
        // Arrange: genes in same order in both genomes (forward collinearity)
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.GreaterThanOrEqualTo(1), "Should find at least one block");
            Assert.That(blocks[0].Strand, Is.EqualTo('+'), "Forward collinearity should have '+' strand");
            Assert.That(blocks[0].GeneCount, Is.GreaterThanOrEqualTo(3), "Block should contain at least minGenes");
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

        // Assert
        Assert.That(blocks.Any(b => b.Strand == '-'), Is.True,
            "Reverse collinearity should produce block with '-' strand");
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

        // Assert
        Assert.That(blocks, Has.Count.GreaterThanOrEqualTo(1),
            "Exactly minGenes should form a valid block");
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

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Count.GreaterThanOrEqualTo(2),
                "Should produce at least 2 blocks for 2 chromosome pairs");
            Assert.That(blocks.Select(b => b.Species1Chromosome).Distinct().Count(),
                Is.GreaterThanOrEqualTo(2), "Blocks should span multiple source chromosomes");
        });
    }

    /// <summary>
    /// M7: GeneCount accurately reflects input genes in block.
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_GeneCountMatchesInput()
    {
        // Arrange: 4 genes
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.That(blocks.Sum(b => b.GeneCount), Is.GreaterThanOrEqualTo(4),
            "Total gene count should reflect input genes");
    }

    /// <summary>
    /// M8: Block coordinates span from first to last gene.
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_CoordinatesSpanAllGenes()
    {
        // Arrange
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(blocks[0].Species1Start, Is.EqualTo(1000),
                "Block should start at first gene");
            Assert.That(blocks[0].Species1End, Is.EqualTo(8000),
                "Block should end at last gene");
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

        // Assert
        Assert.That(rearrangements.Any(r => r.Type == "Inversion"), Is.True,
            "Strand change on same chromosome should detect inversion");
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

        // Assert
        Assert.That(rearrangements.Any(r => r.Type == "Translocation"), Is.True,
            "Chromosome change should detect translocation");
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
    /// S2: SequenceIdentity invariant - must be in [0, 1].
    /// Source: Definition
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_SequenceIdentityInRange()
    {
        // Arrange
        var orthologPairs = CreateForwardCollinearPairs();

        // Act
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            orthologPairs, minGenes: 3, maxGap: 10).ToList();

        // Assert
        Assert.That(blocks.All(b => b.SequenceIdentity >= 0 && b.SequenceIdentity <= 1), Is.True,
            "SequenceIdentity must be in range [0, 1]");
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
        var validTypes = new[] { "Inversion", "Translocation" };
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
