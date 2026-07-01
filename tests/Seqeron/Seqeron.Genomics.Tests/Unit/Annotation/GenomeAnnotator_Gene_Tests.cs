using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

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

    /// <summary>
    /// Reverse-complement of a nucleotide string (A↔T, C↔G, reversed).
    /// </summary>
    private static string ReverseComplement(string s)
    {
        var map = new Dictionary<char, char> { ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C' };
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            chars[i] = map[s[s.Length - 1 - i]];
        return new string(chars);
    }

    /// <summary>
    /// Builds a forward-strand genomic sequence that harbours a single REVERSE-strand gene
    /// with a Shine-Dalgarno motif at the given aligned spacing.
    /// </summary>
    /// <remarks>
    /// Constructed as the reverse complement of <see cref="CreateSequenceWithSd"/>. The
    /// reverse complement of that forward SD+ORF construct, read as a forward genomic
    /// sequence, contains exactly one reverse-strand ORF whose SD motif reads 5'→3' on the
    /// reverse strand. On the forward strand the SD region appears as the anti-SD complement
    /// (e.g. AGGAGG → CCTCCT). The reverse SD's forward-strand 5'-base coordinate equals
    /// (length − 10 − sdMotif.Length), since the SD sits at offset 10 in the source construct.
    /// </remarks>
    private static string CreateSequenceWithReverseStrandSd(
        string sdMotif, int distanceToStart, int orfLength = 100)
    {
        return ReverseComplement(CreateSequenceWithSd(sdMotif, distanceToStart, orfLength));
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

    #region FindRibosomeBindingSitesBothStrands - Reverse-Strand Tests (C6)

    /// <summary>
    /// R1: Reverse-strand Shine-Dalgarno is reported with the correct forward coordinate,
    /// motif (read 5'→3' on the reverse strand) and strand label.
    /// Source: Shine &amp; Dalgarno (1975) Nature 254:34-38 — consensus AGGAGG;
    ///         the mRNA of a reverse-strand gene is the reverse complement of the forward
    ///         genomic strand, so the SD motif lies on the reverse strand upstream of that
    ///         gene's start codon (https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence).
    /// Hand computation: source construct places AGGAGG at offset 10, aligned spacing 8.
    ///         Reverse complement length = 10+6+8+303 = 327. Reverse SD 5'-base in forward
    ///         coordinates = 327 − 10 − 6 = 311. Forward bases there are the anti-SD CCTCCT.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_ReverseStrandConsensus_Detected()
    {
        string sequence = CreateSequenceWithReverseStrandSd("AGGAGG", distanceToStart: 8, orfLength: 100);
        // Guard the hand computation: the forward strand carries the anti-SD, not AGGAGG.
        Assert.That(sequence.Substring(311, 6), Is.EqualTo("CCTCCT"),
            "Forward strand at the reverse SD locus must read the anti-SD complement CCTCCT");

        var sites = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15).ToList();

        var reverseHit = sites.SingleOrDefault(s => s.strand == '-' && s.sequence == "AGGAGG");
        Assert.Multiple(() =>
        {
            Assert.That(reverseHit, Is.Not.EqualTo(default((int, string, double, char))),
                "Reverse-strand AGGAGG must be reported");
            Assert.That(reverseHit.position, Is.EqualTo(311),
                "Reverse SD 5'-base forward coordinate must be 327 − 10 − 6 = 311");
            Assert.That(reverseHit.strand, Is.EqualTo('-'),
                "Hit must be labelled reverse strand");
            Assert.That(reverseHit.score, Is.EqualTo(1.0).Within(1e-12),
                "Full AGGAGG consensus scores 6/6 = 1.0");
        });
    }

    /// <summary>
    /// R2: The forward-only overload must NOT report a reverse-strand SD (behaviour preserved).
    /// The reverse-only construct has zero forward ORFs, so the legacy method yields nothing.
    /// This pins down that wrong strand handling (e.g. scanning the forward strand for AGGAGG)
    /// would change the legacy output.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSites_ForwardOnly_IgnoresReverseStrandSd()
    {
        string sequence = CreateSequenceWithReverseStrandSd("AGGAGG", distanceToStart: 8, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSites(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites, Is.Empty,
            "Forward-only helper must not report reverse-strand SD (no forward ORF present)");
    }

    /// <summary>
    /// R3: Both-strands overload reports BOTH a forward and a reverse SD when both genes exist.
    /// The construct concatenates a forward SD+ORF and a reverse SD+ORF (reverse part offset by
    /// the forward part's length). Exact coordinates derived by hand.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_ForwardAndReverse_BothReported()
    {
        string forwardPart = CreateSequenceWithSd("AGGAGG", distanceToStart: 8, orfLength: 100); // len 327, fwd SD at 10
        string reversePart = CreateSequenceWithReverseStrandSd("AGGAGG", distanceToStart: 8, orfLength: 100); // len 327, rev SD at 311 within itself
        string sequence = forwardPart + reversePart;
        int offset = forwardPart.Length; // 327

        var sites = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Any(s => s.strand == '+' && s.sequence == "AGGAGG" && s.position == 10), Is.True,
                "Forward AGGAGG at position 10 must be reported on '+' strand");
            Assert.That(sites.Any(s => s.strand == '-' && s.sequence == "AGGAGG" && s.position == offset + 311), Is.True,
                "Reverse AGGAGG at forward position 327+311=638 must be reported on '-' strand");
        });
    }

    /// <summary>
    /// R4: Reverse-strand SD at exactly maxDistance (15 bp aligned spacing) is detected.
    /// Source: Wikipedia Shine-Dalgarno — functional aligned spacing window upper bound.
    /// Reverse construct length = 10+6+15+303 = 334; reverse SD 5'-base forward coord = 334−10−6 = 318.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_ReverseAtMaxDistance_Detected()
    {
        string sequence = CreateSequenceWithReverseStrandSd("AGGAGG", distanceToStart: 15, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(
            sequence, upstreamWindow: 30, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.strand == '-' && s.sequence == "AGGAGG" && s.position == 318), Is.True,
            "Reverse AGGAGG at aligned spacing 15 (forward pos 318) must be detected");
    }

    /// <summary>
    /// R5: Reverse-strand SD beyond maxDistance (16 bp) is NOT detected.
    /// Source: Wikipedia Shine-Dalgarno — beyond 15 bp is outside the functional window.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_ReverseBeyondMaxDistance_NotDetected()
    {
        string sequence = CreateSequenceWithReverseStrandSd("AGGAGG", distanceToStart: 16, orfLength: 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(
            sequence, upstreamWindow: 30, minDistance: 4, maxDistance: 15).ToList();

        Assert.That(sites.Any(s => s.strand == '-' && s.sequence == "AGGAGG"), Is.False,
            "Reverse AGGAGG at aligned spacing 16 (above maxDistance) must be filtered out");
    }

    /// <summary>
    /// R6: No ORF on either strand ⇒ no hits from the both-strands overload.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_NoOrfs_ReturnsEmpty()
    {
        string sequence = "CCCCAGGAGGCCCC" + new string('C', 100);

        var sites = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(sequence).ToList();

        Assert.That(sites, Is.Empty);
    }

    /// <summary>
    /// R7: Both-strands overload reproduces the legacy forward-only result for a pure
    /// forward construct (same positions, sequences and scores), with every hit on '+'.
    /// Guards against the both-strands path altering forward behaviour.
    /// </summary>
    [Test]
    public void FindRibosomeBindingSitesBothStrands_ForwardConstruct_MatchesLegacyForwardHits()
    {
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 8, orfLength: 100);

        var legacy = GenomeAnnotator.FindRibosomeBindingSites(
                sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15)
            .OrderBy(s => s.position).ThenBy(s => s.sequence).ToList();
        var both = GenomeAnnotator.FindRibosomeBindingSitesBothStrands(
                sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15)
            .OrderBy(s => s.position).ThenBy(s => s.sequence).ToList();

        Assert.That(both.All(s => s.strand == '+'), Is.True,
            "A pure forward construct must yield only '+' hits");
        Assert.That(both.Select(s => (s.position, s.sequence, s.score)).ToList(),
            Is.EqualTo(legacy), "Both-strands forward hits must equal the legacy forward result");
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
