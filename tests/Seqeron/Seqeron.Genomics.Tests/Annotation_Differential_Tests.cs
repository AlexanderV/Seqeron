// 08_DIFFERENTIAL_TESTING rows 28-31 (Annotation). Independent oracles: ORF spec-reconstruction +
// hand-derived count, PredictGenes-vs-FindOrfs mapping, an independent motif scan for promoter boxes,
// and a hand-built GFF3 line for the serializer.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Annotation_Differential_Tests
{
    private static readonly HashSet<string> Starts = new() { "ATG", "GTG", "TTG" };
    private static readonly HashSet<string> Stops = new() { "TAA", "TAG", "TGA" };

    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
    private static string RevComp(string s)
    {
        var a = s.Select(c => Comp[c]).ToArray();
        Array.Reverse(a);
        return new string(a);
    }

    // ---- Row 28: ANNOT-ORF-001 — FindOrfs vs spec reconstruction + hand-derived count ----

    [Test]
    [Category("ANNOT-ORF-001")]
    public void FindOrfs_ForwardOnly_HandDerivedCount()
    {
        // "ATG AAA TAA": exactly one forward ORF in frame 1, [0,9), 2 aa.
        var orfs = GenomeAnnotator.FindOrfs("ATGAAATAA", minLength: 1, searchBothStrands: false, requireStartCodon: true).ToList();
        Assert.That(orfs.Count, Is.EqualTo(1));
        Assert.That((orfs[0].Start, orfs[0].End, orfs[0].Frame, orfs[0].IsReverseComplement), Is.EqualTo((0, 9, 1, false)));
    }

    [Test]
    [Category("ANNOT-ORF-001")]
    [TestCase("ATGAAATAAGGGCCCTTATTTCATGGG")]
    [TestCase("ATGCATGCATGCTAAGGGATGTTTTGAaaa")]
    public void FindOrfs_EveryOrfSatisfiesSpec(string seq)
    {
        seq = seq.ToUpperInvariant();
        foreach (var orf in GenomeAnnotator.FindOrfs(seq, minLength: 1, searchBothStrands: true, requireStartCodon: true))
        {
            // Coordinate + length invariants.
            Assert.That(orf.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(orf.End, Is.GreaterThan(orf.Start).And.LessThanOrEqualTo(seq.Length));
            Assert.That((orf.End - orf.Start) % 3, Is.EqualTo(0), "ORF length is a codon multiple");

            // The strand sequence is the forward slice (forward strand) or its revcomp (reverse strand).
            string strandSeq = orf.IsReverseComplement
                ? RevComp(seq.Substring(orf.Start, orf.End - orf.Start))
                : seq.Substring(orf.Start, orf.End - orf.Start);
            Assert.That(orf.Sequence, Is.EqualTo(strandSeq), "ORF sequence matches strand slice");

            // Starts with a start codon, ends with a stop codon (requireStartCodon=true path).
            Assert.That(Starts.Contains(strandSeq.Substring(0, 3)), Is.True, "starts with start codon");
            Assert.That(Stops.Contains(strandSeq.Substring(strandSeq.Length - 3, 3)), Is.True, "ends with stop codon");

            // No internal in-frame stop before the terminal one (maximal from start to first stop).
            for (int k = 0; k < strandSeq.Length - 3; k += 3)
                Assert.That(Stops.Contains(strandSeq.Substring(k, 3)), Is.False, $"no internal stop at {k}");
        }
    }

    // ---- Row 29: ANNOT-GENE-001 — PredictGenes vs FindOrfs (ordered) + mapping ----

    [Test]
    [Category("ANNOT-GENE-001")]
    [TestCase("ATGAAATAAGGGCCCTTATTTCATGGGATGCCCTGA")]
    public void PredictGenes_MapsOrderedOrfs(string seq)
    {
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 1).ToList();
        var orfs = GenomeAnnotator.FindOrfs(seq, 1, searchBothStrands: true, requireStartCodon: true)
            .OrderBy(o => o.Start).ToList();

        Assert.That(genes.Count, Is.EqualTo(orfs.Count));
        for (int k = 0; k < genes.Count; k++)
        {
            var g = genes[k];
            var o = orfs[k];
            Assert.That(g.Start, Is.EqualTo(o.Start));
            Assert.That(g.End, Is.EqualTo(o.End));
            Assert.That(g.Strand, Is.EqualTo(o.IsReverseComplement ? '-' : '+'));
            Assert.That(g.GeneId, Is.EqualTo($"gene_{k + 1:D4}"));
            Assert.That(g.Type, Is.EqualTo("CDS"));
            Assert.That(g.Attributes["frame"], Is.EqualTo(o.Frame.ToString()));
            Assert.That(g.Attributes["protein_length"], Is.EqualTo(o.ProteinSequence.TrimEnd('*').Length.ToString()));
        }
    }

    // ---- Row 30: ANNOT-PROM-001 — FindPromoterMotifs vs independent motif scan ----

    private static readonly (string motif, double score)[] Minus35 =
        { ("TTGACA", 1.000), ("TTGAC", 0.855), ("TGACA", 0.815), ("TTGA", 0.710) };
    private static readonly (string motif, double score)[] Minus10 =
        { ("TATAAT", 1.000), ("TATAA", 0.801), ("ATAAT", 0.813), ("TATA", 0.665) };

    [Test]
    [Category("ANNOT-PROM-001")]
    [TestCase("GGTTGACATTTTTTATAATGG")]
    [TestCase("AAAATTGACAGGGGTATAATCCCTTGACA")]
    public void FindPromoterMotifs_MatchesIndependentScan(string seq)
    {
        var actual = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        string s = seq.ToUpperInvariant();
        var expected = new List<(int, string, string, double)>();
        foreach (var (motif, score) in Minus35)
            for (int i = 0; i + motif.Length <= s.Length; i++)
                if (s.Substring(i, motif.Length) == motif) expected.Add((i, "-35 box", motif, score));
        foreach (var (motif, score) in Minus10)
            for (int i = 0; i + motif.Length <= s.Length; i++)
                if (s.Substring(i, motif.Length) == motif) expected.Add((i, "-10 box", motif, score));

        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 31: ANNOT-GFF-001 — ToGff3 serializer vs hand-built GFF3 lines ----

    [Test]
    [Category("ANNOT-GFF-001")]
    public void ToGff3_SingleCds_MatchesHandBuiltLine()
    {
        var attrs = new Dictionary<string, string> { ["frame"] = "1", ["protein_length"] = "2", ["translation"] = "MK" };
        var ann = new GenomeAnnotator.GeneAnnotation("gene_0001", 0, 9, '+', "CDS", "hypothetical protein", attrs);

        var lines = GenomeAnnotator.ToGff3(new[] { ann }, "seq1").ToList();

        Assert.That(lines[0], Is.EqualTo("##gff-version 3"));
        // seqId  source(.)  type  start+1  end  score(.)  strand  phase(0 for sole CDS)  attrs(translation skipped)
        Assert.That(lines[1], Is.EqualTo("seq1\t.\tCDS\t1\t9\t.\t+\t0\tID=gene_0001;product=hypothetical protein;frame=1;protein_length=2"));
    }

    [Test]
    [Category("ANNOT-GFF-001")]
    public void ToGff3_EncodesReservedCharsAndNonCdsPhase()
    {
        var attrs = new Dictionary<string, string>();
        // Non-CDS feature => phase "."; product with reserved ';' => %3B encoding.
        var ann = new GenomeAnnotator.GeneAnnotation("g1", 4, 10, '-', "gene", "a;b", attrs);

        var lines = GenomeAnnotator.ToGff3(new[] { ann }, "chr1").ToList();
        Assert.That(lines[1], Is.EqualTo("chr1\t.\tgene\t5\t10\t.\t-\t.\tID=g1;product=a%3Bb"));
    }
}
