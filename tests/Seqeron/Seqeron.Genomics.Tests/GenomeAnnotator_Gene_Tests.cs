using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Evidence-based tests for GenomeAnnotator gene prediction methods.
/// Test Unit: ANNOT-GENE-001
/// 
/// Canonical methods:
///   - GenomeAnnotator.PredictGenes(dna, minOrfLength, prefix)
///   - GenomeAnnotator.FindRibosomeBindingSites(dna, window, minDist, maxDist)
/// 
/// Evidence:
///   - Wikipedia: Gene prediction, Shine-Dalgarno sequence, Ribosome-binding site
///   - Shine &amp; Dalgarno (1975): SD consensus AGGAGG
///   - Chen et al. (1994): Optimal aligned spacing 5 nt
///   - Laursen et al. (2005): Bacterial translation initiation
/// </summary>
[TestFixture]
[Category("Annotation")]
[Category("GenePrediction")]
public class GenomeAnnotator_Gene_Tests
{
    #region Test Data

    /// <summary>
    /// Minimal valid ORF: ATG + 100 codons (300 bp) + stop codon.
    /// Creates exactly 100 amino acid ORF.
    /// </summary>
    private static string CreateMinimalOrf(int aminoAcidCount = 100)
    {
        // ATG start + (aminoAcidCount * 3 - 3 for stop) nucleotides + TAA stop
        // Total ORF length = aminoAcidCount * 3 + 3 (including stop)
        return "ATG" + new string('A', (aminoAcidCount - 1) * 3) + "TAA";
    }

    /// <summary>
    /// Creates a sequence with Shine-Dalgarno at specific aligned spacing from start codon.
    /// The aligned spacing is measured from SD 3'-end to the first nucleotide of AUG,
    /// per Chen et al. (1994) definition.
    /// </summary>
    private static string CreateSequenceWithSd(string sdMotif, int distanceToStart, int orfLength = 100)
    {
        // Structure: padding + SD + spacer + ORF
        string padding = new string('C', 10); // Avoid false positives
        string spacer = new string('C', distanceToStart);
        string orf = CreateMinimalOrf(orfLength);
        return padding + sdMotif + spacer + orf;
    }

    #endregion

    #region PredictGenes - Must Tests

    /// <summary>
    /// M1: All predicted genes have Type = "CDS"
    /// Source: Gene annotation standard
    /// </summary>
    [Test]
    public void PredictGenes_AllGenesHaveCdsType()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Not.Empty);
        Assert.That(genes.All(g => g.Type == "CDS"), Is.True,
            "All predicted genes should have Type = 'CDS'");
    }

    /// <summary>
    /// M2: Gene IDs follow pattern "{prefix}_{number:D4}" with sequential numbering
    /// Source: Implementation contract
    /// </summary>
    [Test]
    public void PredictGenes_AssignsSequentialGeneIds()
    {
        string sequence = CreateMinimalOrf(100) + new string('C', 20) + CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50, prefix: "test").ToList();

        Assert.That(genes.Count, Is.GreaterThanOrEqualTo(2));
        for (int i = 0; i < genes.Count; i++)
        {
            Assert.That(genes[i].GeneId, Is.EqualTo($"test_{(i + 1):D4}"),
                $"Gene at index {i} should have sequential ID test_{(i + 1):D4}");
        }
    }

    /// <summary>
    /// M3: All genes have strand info ('+' or '-')
    /// Source: Wikipedia - genes exist on both strands
    /// </summary>
    [Test]
    public void PredictGenes_IncludesStrandInformation()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            foreach (var gene in genes)
            {
                Assert.That(gene.Strand, Is.EqualTo('+').Or.EqualTo('-'),
                    $"Gene {gene.GeneId} should have strand '+' or '-'");
            }
        });
    }

    /// <summary>
    /// M4: ORFs shorter than minOrfLength are filtered out
    /// Source: Parameter contract
    /// </summary>
    [Test]
    public void PredictGenes_FiltersOrfsByMinLength()
    {
        // Create ORF with 50 amino acids
        string shortOrf = "ATG" + new string('A', 147) + "TAA"; // 50 codons including start

        var genesStrict = GenomeAnnotator.PredictGenes(shortOrf, minOrfLength: 100).ToList();
        var genesLoose = GenomeAnnotator.PredictGenes(shortOrf, minOrfLength: 30).ToList();

        Assert.That(genesStrict, Is.Empty, "50 aa ORF should be filtered with minOrfLength=100");
        Assert.That(genesLoose.Count, Is.GreaterThan(0), "50 aa ORF should pass with minOrfLength=30");
    }

    /// <summary>
    /// M5: Genes found on both forward and reverse strands
    /// Source: Wikipedia - genes can be on either strand
    /// </summary>
    [Test]
    public void PredictGenes_FindsGenesOnBothStrands()
    {
        // Forward strand ORF: ATG + poly-A + TAA
        string forwardOrf = CreateMinimalOrf(100);
        string spacer = new string('C', 50);
        // Reverse strand ORF: place revcomp of an ORF so the reverse strand reads ATG...TAA
        string reverseTarget = "ATG" + new string('G', 297) + "TAA";
        string reverseOrfOnFwd = DnaSequence.GetReverseComplementString(reverseTarget);
        string sequence = forwardOrf + spacer + reverseOrfOnFwd;

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes.Any(g => g.Strand == '+'), Is.True,
            "Should find gene on forward strand");
        Assert.That(genes.Any(g => g.Strand == '-'), Is.True,
            "Should find gene on reverse strand");
    }

    /// <summary>
    /// M8: Empty sequence returns empty result
    /// Source: Edge case definition
    /// </summary>
    [Test]
    public void PredictGenes_EmptySequence_ReturnsEmpty()
    {
        var genes = GenomeAnnotator.PredictGenes("", minOrfLength: 50).ToList();

        Assert.That(genes, Is.Empty);
    }

    /// <summary>
    /// M9: Sequence without valid ORFs returns empty
    /// Source: No start/stop pattern = no gene
    /// </summary>
    [Test]
    public void PredictGenes_NoValidOrfs_ReturnsEmpty()
    {
        // Sequence with no ATG start codon
        string noStart = new string('A', 500) + "TAA";

        var genes = GenomeAnnotator.PredictGenes(noStart, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Empty);
    }

    /// <summary>
    /// M10: Protein length in attributes matches (End-Start)/3 - 1
    /// Source: Biological definition — protein length excludes stop codon.
    ///         ORF span (End-Start) includes stop codon triplet, so amino acids = codons - 1.
    /// </summary>
    [Test]
    public void PredictGenes_ProteinLengthAttributeIsAccurate()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Not.Empty);
        foreach (var gene in genes)
        {
            Assert.That(gene.Attributes.ContainsKey("protein_length"), Is.True,
                $"Gene {gene.GeneId} should have protein_length attribute");
            int proteinLength = int.Parse(gene.Attributes["protein_length"]);
            int expected = (gene.End - gene.Start) / 3 - 1;
            Assert.That(proteinLength, Is.EqualTo(expected),
                $"Gene {gene.GeneId}: protein_length should be (End-Start)/3 - 1 = {expected} (excludes stop codon)");
        }
    }

    #endregion

    #region PredictGenes - Should Tests

    /// <summary>
    /// S1: Multiple genes in sequence are all detected with correct coordinates
    /// Source: Multi-gene operons are common in prokaryotes
    /// </summary>
    [Test]
    public void PredictGenes_MultipleGenes_AllDetected()
    {
        // Two distinct ORFs separated by intergenic region
        string orf1 = CreateMinimalOrf(100); // 303 bp: ATG + 297 A's + TAA
        string spacer = new string('C', 20);
        string orf2 = CreateMinimalOrf(100); // 303 bp
        string sequence = orf1 + spacer + orf2;

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();
        var forwardGenes = genes.Where(g => g.Strand == '+').OrderBy(g => g.Start).ToList();

        Assert.That(forwardGenes.Count, Is.GreaterThanOrEqualTo(2),
            "Should detect at least 2 forward-strand genes");
        Assert.That(forwardGenes[0].Start, Is.EqualTo(0), "First gene starts at 0");
        Assert.That(forwardGenes[0].End, Is.EqualTo(303), "First gene ends at 303");
        Assert.That(forwardGenes[1].Start, Is.EqualTo(323), "Second gene starts at 323");
        Assert.That(forwardGenes[1].End, Is.EqualTo(626), "Second gene ends at 626");
    }

    /// <summary>
    /// S3: Alternative start codons (GTG, TTG) should be recognized
    /// Source: Wikipedia: Gene prediction — prokaryotic start codons
    /// </summary>
    [Test]
    [TestCase("GTG")]
    [TestCase("TTG")]
    public void PredictGenes_AlternativeStartCodons_Recognized(string startCodon)
    {
        string orf = startCodon + new string('A', 297) + "TAA";

        var genes = GenomeAnnotator.PredictGenes(orf, minOrfLength: 50).ToList();
        var forwardGenes = genes.Where(g => g.Strand == '+').ToList();

        Assert.That(forwardGenes, Is.Not.Empty,
            $"{startCodon} start codon should be recognized");
    }

    #endregion

    #region FindRibosomeBindingSites - Must Tests

    /// <summary>
    /// M6: Detects canonical AGGAGG Shine-Dalgarno sequence
    /// Source: Shine &amp; Dalgarno (1975) - AGGAGG is the consensus
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_DetectsConsensusAggagg()
    {
        // Place AGGAGG 8bp upstream of ATG (within functional range 4-15 bp)
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 8, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.sequence == "AGGAGG"), Is.True,
            "Should detect exact AGGAGG Shine-Dalgarno consensus");
    }

    /// <summary>
    /// RBS too close to start codon should not be detected.
    /// Source: Wikipedia: Shine-Dalgarno sequence — functional range 4-15 bp;
    ///         below 4 bp, ribosome cannot bind properly.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_TooClose_NotDetected()
    {
        // SD at only 2bp aligned spacing from start (below minDistance=4)
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 2, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15).ToList();

        // AGGAGG at aligned spacing 2 should NOT be returned (below minDistance=4).
        // Note: shorter sub-motifs (e.g., AGGA) at the same genomic position
        // have larger aligned spacing due to shorter length, so they may still be valid.
        var fullConsensusAtInvalidDist = sites.Where(s =>
            s.sequence == "AGGAGG" &&
            s.position == 10); // padding(10) is where SD starts in CreateSequenceWithSd
        Assert.That(fullConsensusAtInvalidDist.Count(), Is.EqualTo(0),
            "Full AGGAGG at aligned spacing 2 (below minDistance=4) should be filtered out");
    }

    /// <summary>
    /// Optimal aligned spacing is 5 nt per Chen et al. (1994).
    /// Source: Chen H et al. (1994) Nucleic Acids Research 22(23):4953-4957.
    ///         "Measurements of protein expression demonstrated an optimal aligned
    ///         spacing of 5 nt for both series."
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_OptimalSpacing_5nt_ChenEtAl1994()
    {
        // Chen et al. (1994) found 5 nt is the optimal aligned spacing
        // for SD-to-AUG in E. coli
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 5, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15).ToList();

        var aggaggSites = sites.Where(s => s.sequence == "AGGAGG").ToList();
        Assert.That(aggaggSites, Is.Not.Empty,
            "AGGAGG at optimal aligned spacing of 5 nt (Chen et al. 1994) must be detected");
    }

    /// <summary>
    /// Boundary: SD at exactly minDistance (4 bp) should be detected.
    /// Source: Wikipedia: Shine-Dalgarno sequence — 4 bp is the lower boundary
    ///         of the functional range.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_AtMinDistance_Detected()
    {
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 4, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.sequence == "AGGAGG"), Is.True,
            "AGGAGG at exactly minDistance (4 bp) must be detected");
    }

    /// <summary>
    /// Boundary: SD at exactly maxDistance (15 bp) should be detected.
    /// Source: Wikipedia: Shine-Dalgarno sequence — 15 bp is the upper boundary
    ///         of the functional range.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_AtMaxDistance_Detected()
    {
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 15, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 30, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.sequence == "AGGAGG"), Is.True,
            "AGGAGG at exactly maxDistance (15 bp) must be detected");
    }

    /// <summary>
    /// Boundary: SD beyond maxDistance (16 bp) should NOT be detected.
    /// Source: Wikipedia: Shine-Dalgarno sequence — beyond 15 bp is outside functional range.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_BeyondMaxDistance_NotDetected()
    {
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 16, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 30, minDistance: 4, maxDistance: 15).ToList();

        // AGGAGG at aligned spacing 16 should NOT be detected (above maxDistance=15)
        var aggaggSites = sites.Where(s => s.sequence == "AGGAGG").ToList();
        Assert.That(aggaggSites, Is.Empty,
            "AGGAGG at aligned spacing 16 (above maxDistance=15) should be filtered out");
    }

    #endregion

    #region FindRibosomeBindingSites - Should Tests

    /// <summary>
    /// S4: Shorter SD variants (GGAGG, AGGAG, GAGG, AGGA) are detected
    /// Source: Wikipedia - variant SD sequences exist
    /// </summary>
    [Test]
    [TestCase("GGAGG")]
    [TestCase("AGGAG")]
    [TestCase("GAGG")]
    [TestCase("AGGA")]
    public void FindRibosomeBindingSites_ShorterMotifs_Detected(string sdMotif)
    {
        string sequence = CreateSequenceWithSd(sdMotif, distanceToStart: 8, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.sequence == sdMotif), Is.True,
            $"Should detect exact shorter SD variant: {sdMotif}");
    }

    /// <summary>
    /// S5: Score reflects motif length (longer = higher score)
    /// Source: Implementation - quality metric
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_ScoreReflectsMotifLength()
    {
        // CreateSequenceWithSd places "AGGAGG" upstream of an ORF
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 8, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15).ToList();

        // Sequence with AGGAGG and ORF MUST find RBS sites
        Assert.That(sites, Is.Not.Empty,
            "Sequence with AGGAGG SD motif and ORF must produce RBS sites");

        // Score should be normalized by consensus length (6)
        var aggaggSite = sites.FirstOrDefault(s => s.sequence == "AGGAGG");
        Assert.That(aggaggSite, Is.Not.EqualTo(default((int, string, double))),
            "AGGAGG motif must be found");
        Assert.That(aggaggSite.score, Is.EqualTo(1.0).Within(0.01),
            "Full AGGAGG should have score of 1.0 (6/6)");
    }

    /// <summary>
    /// No ORFs = no RBS sites
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_NoOrfs_ReturnsEmpty()
    {
        // Sequence with SD but no ORF
        string sequence = "CCCCAGGAGGCCCC" + new string('C', 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(sequence).ToList();

        // Should be empty because no ORF to associate with
        Assert.That(sites, Is.Empty);
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Case-insensitive sequence handling
    /// </summary>
    [Test]
    public void PredictGenes_MixedCase_HandledCorrectly()
    {
        string lowerCaseOrf = "atg" + new string('a', 297) + "taa";

        var genes = GenomeAnnotator.PredictGenes(lowerCaseOrf, minOrfLength: 50).ToList();

        Assert.That(genes.Count, Is.GreaterThan(0),
            "Should handle lowercase sequences");
    }

    /// <summary>
    /// Null handling
    /// </summary>
    [Test]
    public void PredictGenes_NullSequence_HandledGracefully()
    {
        var genes = GenomeAnnotator.PredictGenes(null!, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Empty);
    }

    /// <summary>
    /// Default prefix is "gene"
    /// </summary>
    [Test]
    public void PredictGenes_DefaultPrefix_IsGene()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes.Count, Is.GreaterThan(0));
        Assert.That(genes[0].GeneId, Does.StartWith("gene_"));
    }

    /// <summary>
    /// Gene coordinates are valid (Start < End)
    /// </summary>
    [Test]
    public void PredictGenes_CoordinatesAreValid()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.Multiple(() =>
        {
            foreach (var gene in genes)
            {
                Assert.That(gene.Start, Is.LessThan(gene.End),
                    $"Gene {gene.GeneId}: Start ({gene.Start}) should be < End ({gene.End})");
                Assert.That(gene.Start, Is.GreaterThanOrEqualTo(0),
                    $"Gene {gene.GeneId}: Start should be >= 0");
            }
        });
    }

    /// <summary>
    /// Frame attribute is present and valid (1, 2, or 3)
    /// </summary>
    [Test]
    public void PredictGenes_FrameAttributeIsValid()
    {
        string sequence = CreateMinimalOrf(100);

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();

        Assert.That(genes, Is.Not.Empty);
        foreach (var gene in genes)
        {
            Assert.That(gene.Attributes.ContainsKey("frame"), Is.True,
                $"Gene {gene.GeneId} should have frame attribute");
            int frame = int.Parse(gene.Attributes["frame"]);
            Assert.That(frame, Is.InRange(1, 3),
                $"Frame should be 1, 2, or 3; got {frame}");
        }
    }

    /// <summary>
    /// S2: Overlapping genes in different reading frames are both reported
    /// Source: Prokaryotic genomes can have overlapping genes in different frames
    /// </summary>
    [Test]
    public void PredictGenes_OverlappingGenes_BothReported()
    {
        // Frame 0: ATG at pos 0, GGG fill, TAA stop at pos 300–302 (100 aa)
        // Frame 1: ATG at pos 4, GGG fill, TAA stop at pos 304–306 (100 aa)
        // Overlap: positions [4, 303) — 299 bp overlap
        string sequence = "ATG" + "G" + "ATG" + new string('G', 293) + "TAA" + "G" + "TAA";

        var genes = GenomeAnnotator.PredictGenes(sequence, minOrfLength: 50).ToList();
        var forwardGenes = genes.Where(g => g.Strand == '+').OrderBy(g => g.Start).ToList();

        Assert.That(forwardGenes.Count, Is.GreaterThanOrEqualTo(2),
            "Overlapping genes in different frames should both be reported");
        Assert.That(forwardGenes[1].Start, Is.LessThan(forwardGenes[0].End),
            "Genes should overlap (gene2.Start < gene1.End)");
    }

    /// <summary>
    /// C1: Very long ORF (>1000 amino acids) is handled correctly
    /// Source: Stress test — large genes exist in prokaryotic genomes
    /// </summary>
    [Test]
    public void PredictGenes_VeryLongOrf_Handled()
    {
        string longOrf = CreateMinimalOrf(1500); // 1500 aa = 4503 bp

        var genes = GenomeAnnotator.PredictGenes(longOrf, minOrfLength: 100).ToList();
        var forwardGenes = genes.Where(g => g.Strand == '+').ToList();

        Assert.That(forwardGenes, Is.Not.Empty,
            "Very long ORF (>1000 aa) should be handled");
        Assert.That(forwardGenes[0].End - forwardGenes[0].Start, Is.GreaterThan(3000),
            "ORF should span >3000 nucleotides");
    }

    /// <summary>
    /// Edge: ATG without stop codon produces no valid gene (too short / no stop)
    /// </summary>
    [Test]
    public void PredictGenes_StartCodonOnly_NoStop_ReturnsEmpty()
    {
        var genes = GenomeAnnotator.PredictGenes("ATG", minOrfLength: 50).ToList();

        Assert.That(genes, Is.Empty,
            "ATG alone (no stop codon) should produce no genes");
    }

    /// <summary>
    /// C3: Multiple RBS motifs upstream of the same ORF are all detected
    /// Source: Multiple SD-like sequences can exist in upstream region
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_MultipleUpstreamOfSameOrf()
    {
        // AGGAGG at position 10, aligned spacing = 29 - 10 - 6 = 13
        // GAGG at position 20, aligned spacing = 29 - 20 - 4 = 5
        // ATG starts at position 29
        string sequence = new string('C', 10) + "AGGAGG" + new string('C', 4)
            + "GAGG" + new string('C', 5) + CreateMinimalOrf(100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 30, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.sequence == "AGGAGG" && s.position == 10), Is.True,
            "Should detect AGGAGG at position 10 (aligned spacing 13)");
        Assert.That(sites.Any(s => s.sequence == "GAGG" && s.position == 20), Is.True,
            "Should detect standalone GAGG at position 20 (aligned spacing 5)");
    }

    #endregion
}
