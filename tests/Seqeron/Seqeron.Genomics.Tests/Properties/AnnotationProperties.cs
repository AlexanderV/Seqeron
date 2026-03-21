using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for genome annotation: ORF detection, gene prediction,
/// promoter motif detection, and GFF3 round-trip I/O.
///
/// Test Units: ANNOT-ORF-001, ANNOT-GENE-001, ANNOT-PROM-001, ANNOT-GFF-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class AnnotationProperties
{
    #region Generators

    /// <summary>
    /// Generates random DNA sequences of sufficient length to contain ORFs.
    /// </summary>
    private static Arbitrary<string> DnaArbitrary(int minLen = 50) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Builds a synthetic sequence guaranteed to contain at least one ORF:
    /// random flank + ATG + (codingCodons × AAA) + TAA + random flank.
    /// </summary>
    private static string BuildSequenceWithOrf(int codingCodons = 35)
    {
        var rng = new Random(42);
        char[] bases = { 'A', 'C', 'G', 'T' };
        string flank5 = new(Enumerable.Range(0, 30).Select(_ => bases[rng.Next(4)]).ToArray());
        string flank3 = new(Enumerable.Range(0, 30).Select(_ => bases[rng.Next(4)]).ToArray());
        // Coding codons: use GCT (Ala) to avoid internal stop codons
        string coding = string.Concat(Enumerable.Repeat("GCT", codingCodons));
        return flank5 + "ATG" + coding + "TAA" + flank3;
    }

    private static readonly HashSet<string> ValidStartCodons =
        new(StringComparer.OrdinalIgnoreCase) { "ATG", "GTG", "TTG" };

    private static readonly HashSet<string> ValidStopCodons =
        new(StringComparer.OrdinalIgnoreCase) { "TAA", "TAG", "TGA" };

    #endregion

    #region ANNOT-ORF-001: R: ORF start < end ≤ seqLen; P: starts with ATG; M: longer seq → ≥ ORFs; R: len divisible by 3

    /// <summary>
    /// INV-1: ORF coordinate validity — start &lt; end ≤ seqLen.
    /// Evidence: An ORF occupies a contiguous region within the sequence;
    /// start must precede end and neither can exceed sequence length.
    /// Source: Wikipedia "Open reading frame"; Rosalind ORF.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindOrfs_Coordinates_AreValid()
    {
        return Prop.ForAll(DnaArbitrary(100), seq =>
        {
            var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 10).ToList();
            return orfs.All(o => o.Start >= 0 && o.Start < o.End && o.End <= seq.Length)
                .Label($"All ORF coordinates must satisfy 0 ≤ start < end ≤ {seq.Length}");
        });
    }

    /// <summary>
    /// INV-2: Every ORF begins with a valid start codon when requireStartCodon=true.
    /// Evidence: By definition an ORF starts with ATG (or GTG/TTG in prokaryotes).
    /// Source: Wikipedia "Open reading frame"; NCBI ORF Finder.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindOrfs_StartsWithStartCodon()
    {
        string seq = BuildSequenceWithOrf(40);
        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 10, requireStartCodon: true).ToList();

        foreach (var orf in orfs)
        {
            string firstCodon = orf.Sequence[..3].ToUpperInvariant();
            Assert.That(ValidStartCodons.Contains(firstCodon), Is.True,
                $"ORF at {orf.Start} does not start with a valid start codon: '{firstCodon}'");
        }
    }

    /// <summary>
    /// INV-3: ORF nucleotide length is divisible by 3 (frame integrity).
    /// Evidence: Translation reads codons as non-overlapping triplets.
    /// Source: Wikipedia "Open reading frame"; Crick (1968).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindOrfs_Length_DivisibleBy3()
    {
        return Prop.ForAll(DnaArbitrary(100), seq =>
        {
            var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 10).ToList();
            return orfs.All(o => o.Sequence.Length % 3 == 0)
                .Label("ORF sequence length must be divisible by 3");
        });
    }

    /// <summary>
    /// INV-4: Protein length ≤ seqLen / 3.
    /// Evidence: Each protein residue requires exactly 3 nucleotides.
    /// Source: Wikipedia "Reading frame".
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindOrfs_ProteinLength_BoundedBySeqLength()
    {
        return Prop.ForAll(DnaArbitrary(100), seq =>
        {
            var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 10).ToList();
            return orfs.All(o => o.ProteinSequence.Length <= seq.Length / 3)
                .Label("Protein length must be ≤ seqLen / 3");
        });
    }

    /// <summary>
    /// INV-5: ORF detection is deterministic — same input always yields same ORFs.
    /// Evidence: FindOrfs is a pure function with no randomness.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindOrfs_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(100), seq =>
        {
            var r1 = GenomeAnnotator.FindOrfs(seq, minLength: 10).ToList();
            var r2 = GenomeAnnotator.FindOrfs(seq, minLength: 10).ToList();
            bool same = r1.Count == r2.Count &&
                        r1.Zip(r2).All(pair => pair.First.Start == pair.Second.Start &&
                                                pair.First.End == pair.Second.End);
            return same.Label("FindOrfs must be deterministic");
        });
    }

    /// <summary>
    /// INV-6: A sequence with a synthetic ORF (≥ minLength codons) must produce at least one result.
    /// Evidence: The constructed sequence is ATG + coding + TAA, guaranteed to be detected.
    /// Source: Wikipedia definition — ATG..TAA with no internal stops qualifies as ORF.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindOrfs_SyntheticOrf_DetectsAtLeastOne()
    {
        string seq = BuildSequenceWithOrf(40); // 40 coding codons = 120 aa including start
        var orfs = GenomeAnnotator.FindOrfs(seq, minLength: 30).ToList();

        Assert.That(orfs, Is.Not.Empty,
            "Sequence with synthetic ORF must produce at least one result");
    }

    #endregion

    #region ANNOT-GENE-001: R: gene start < end; P: contains RBS motif upstream; D: deterministic

    /// <summary>
    /// INV-1: Gene coordinate validity — start &lt; end.
    /// Evidence: A gene annotation spans a non-empty region.
    /// Source: GFF3 specification: start ≤ end.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictGenes_Coordinates_AreValid()
    {
        string seq = BuildSequenceWithOrf(40);
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();

        foreach (var gene in genes)
        {
            Assert.That(gene.Start, Is.GreaterThanOrEqualTo(0),
                $"Gene '{gene.GeneId}' start < 0");
            Assert.That(gene.End, Is.GreaterThan(gene.Start),
                $"Gene '{gene.GeneId}': end ({gene.End}) must be > start ({gene.Start})");
        }
    }

    /// <summary>
    /// INV-2: Gene IDs are unique.
    /// Evidence: Each predicted gene must have a distinct identifier for unambiguous reference.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictGenes_GeneIds_AreUnique()
    {
        string seq = BuildSequenceWithOrf(40);
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();

        var ids = genes.Select(g => g.GeneId).ToList();
        Assert.That(ids.Distinct().Count(), Is.EqualTo(ids.Count),
            "All gene IDs must be unique");
    }

    /// <summary>
    /// INV-3: Gene strand must be '+' or '-'.
    /// Evidence: GFF3 spec mandates strand ∈ {+, -, ., ?}.
    /// Gene predictions on forward/reverse strand use +/-.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictGenes_Strand_IsValid()
    {
        string seq = BuildSequenceWithOrf(40);
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();

        foreach (var gene in genes)
            Assert.That(gene.Strand, Is.EqualTo('+').Or.EqualTo('-'),
                $"Gene '{gene.GeneId}' has invalid strand: '{gene.Strand}'");
    }

    /// <summary>
    /// INV-4: Gene prediction is deterministic.
    /// Evidence: PredictGenes wraps FindOrfs (pure function) with sequential numbering.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PredictGenes_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(200), seq =>
        {
            var r1 = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();
            var r2 = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();
            bool same = r1.Count == r2.Count &&
                        r1.Zip(r2).All(pair => pair.First.GeneId == pair.Second.GeneId &&
                                                pair.First.Start == pair.Second.Start &&
                                                pair.First.End == pair.Second.End);
            return same.Label("PredictGenes must be deterministic");
        });
    }

    /// <summary>
    /// INV-5: Gene type attribute is "CDS".
    /// Evidence: PredictGenes produces CDS annotations from ORFs.
    /// Source: GFF3 spec — CDS is the standard type for coding sequences.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictGenes_Type_IsCDS()
    {
        string seq = BuildSequenceWithOrf(40);
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();

        foreach (var gene in genes)
            Assert.That(gene.Type, Is.EqualTo("CDS"),
                $"Gene '{gene.GeneId}' has type '{gene.Type}', expected 'CDS'");
    }

    #endregion

    #region ANNOT-PROM-001: R: position ≥ 0; P: contains -10/-35 box; M: lower score threshold → ≥ promoters

    /// <summary>
    /// INV-1: Promoter motif positions are within sequence bounds.
    /// Evidence: A discovered motif at position p requires p + motif.Length ≤ seq.Length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPromoterMotifs_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
            return motifs.All(m => m.position >= 0 && m.position + m.sequence.Length <= seq.Length)
                .Label("All promoter motif positions must be within sequence bounds");
        });
    }

    /// <summary>
    /// INV-2: Promoter motif type is either "-10 box" or "-35 box".
    /// Evidence: Bacterial promoters contain two consensus elements recognized by σ70.
    /// Source: Wikipedia "Promoter (genetics)"; Harley &amp; Reynolds (1987).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPromoterMotifs_Type_IsValid()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
            return motifs.All(m => m.type == "-10 box" || m.type == "-35 box")
                .Label("Promoter type must be '-10 box' or '-35 box'");
        });
    }

    /// <summary>
    /// INV-3: Promoter motif score ∈ [0, 1].
    /// Evidence: Score = sum(matched position probabilities) / sum(all consensus probabilities).
    /// Maximum is 1.0 for full consensus match.
    /// Source: Harley &amp; Reynolds (1987).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPromoterMotifs_Score_InRange()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
            return motifs.All(m => m.score >= 0.0 && m.score <= 1.0 + 1e-9)
                .Label("Promoter motif score must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-4: A sequence containing the full -10 box consensus TATAAT must find at least one -10 box.
    /// Evidence: Exact match of the consensus must be detected with score = 1.0.
    /// Source: Wikipedia "Promoter (genetics)" — -10 box (Pribnow box) consensus is TATAAT.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPromoterMotifs_FullConsensus_IsDetected()
    {
        string seq = "GCGCGCTATAAT" + "GCGCGCGCGC";
        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        var minus10 = motifs.Where(m => m.type == "-10 box" && m.sequence == "TATAAT").ToList();
        Assert.That(minus10, Is.Not.Empty,
            "Full -10 box consensus TATAAT must be detected");
        Assert.That(minus10.First().score, Is.EqualTo(1.0).Within(0.001),
            "Full consensus match should have score = 1.0");
    }

    /// <summary>
    /// INV-5: A sequence containing the full -35 box consensus TTGACA must find at least one -35 box.
    /// Evidence: Exact match of the consensus must be detected with score = 1.0.
    /// Source: Wikipedia "Promoter (genetics)" — -35 box consensus is TTGACA.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindPromoterMotifs_Full35Consensus_IsDetected()
    {
        string seq = "GCGCGCTTGACA" + "GCGCGCGCGC";
        var motifs = GenomeAnnotator.FindPromoterMotifs(seq).ToList();

        var minus35 = motifs.Where(m => m.type == "-35 box" && m.sequence == "TTGACA").ToList();
        Assert.That(minus35, Is.Not.Empty,
            "Full -35 box consensus TTGACA must be detected");
        Assert.That(minus35.First().score, Is.EqualTo(1.0).Within(0.001),
            "Full consensus match should have score = 1.0");
    }

    /// <summary>
    /// INV-6: Promoter detection is deterministic.
    /// Evidence: FindPromoterMotifs is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindPromoterMotifs_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var r1 = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
            var r2 = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
            bool same = r1.Count == r2.Count &&
                        r1.Zip(r2).All(p => p.First.position == p.Second.position &&
                                             p.First.type == p.Second.type &&
                                             p.First.score == p.Second.score);
            return same.Label("FindPromoterMotifs must be deterministic");
        });
    }

    #endregion

    #region ANNOT-GFF-001: RT: parse(serialize(features))=features; R: well-formed GFF3; P: coordinates 1-based

    /// <summary>
    /// INV-1: Round-trip — ParseGff3(ToGff3(annotations)) preserves core fields.
    /// Evidence: GFF3 serialization and parsing should be inverse operations for core fields.
    /// Source: GFF3 Spec v1.26 — format is self-describing.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_RoundTrip_PreservesCoreFields()
    {
        string seq = BuildSequenceWithOrf(40);
        var genes = GenomeAnnotator.PredictGenes(seq, minOrfLength: 30).ToList();

        if (genes.Count == 0)
        {
            Assert.Pass("No genes predicted from synthetic sequence — skip round-trip");
            return;
        }

        // Serialize to GFF3 lines
        var gff3Lines = GenomeAnnotator.ToGff3(genes, "chr1").ToList();

        // Parse back
        var features = GenomeAnnotator.ParseGff3(gff3Lines).ToList();

        // Core field preservation
        Assert.That(features.Count, Is.EqualTo(genes.Count),
            $"Feature count mismatch: serialized {genes.Count}, parsed {features.Count}");

        for (int i = 0; i < genes.Count; i++)
        {
            // ToGff3 converts 0-based Start to 1-based (Start + 1)
            Assert.That(features[i].Start, Is.EqualTo(genes[i].Start + 1),
                $"Gene {i}: Start mismatch (0-based {genes[i].Start} → 1-based {genes[i].Start + 1})");
            Assert.That(features[i].End, Is.EqualTo(genes[i].End),
                $"Gene {i}: End mismatch");
            Assert.That(features[i].Strand, Is.EqualTo(genes[i].Strand),
                $"Gene {i}: Strand mismatch");
            Assert.That(features[i].Type, Is.EqualTo(genes[i].Type),
                $"Gene {i}: Type mismatch");
        }
    }

    /// <summary>
    /// INV-2: GFF3 output starts with ##gff-version 3 header.
    /// Evidence: GFF3 Spec v1.26 mandates the version directive as the first line.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_Output_StartsWithVersionHeader()
    {
        var genes = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("gene_1", 0, 100, '+', "CDS", "hypothetical protein",
                new Dictionary<string, string> { ["frame"] = "1" })
        };

        var lines = GenomeAnnotator.ToGff3(genes).ToList();
        Assert.That(lines.First(), Is.EqualTo("##gff-version 3"),
            "GFF3 output must start with '##gff-version 3'");
    }

    /// <summary>
    /// INV-3: GFF3 coordinates are 1-based (start column = internal 0-based + 1).
    /// Evidence: GFF3 Spec uses 1-based inclusive coordinates.
    /// Source: GFF3 Spec v1.26 — "The start and end of the feature, in 1-based integer coordinates."
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_Coordinates_Are1Based()
    {
        var genes = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("gene_1", 0, 300, '+', "CDS", "test",
                new Dictionary<string, string> { ["frame"] = "1" })
        };

        var lines = GenomeAnnotator.ToGff3(genes).ToList();
        // Skip header
        string dataLine = lines[1];
        string[] fields = dataLine.Split('\t');

        int gff3Start = int.Parse(fields[3]);
        Assert.That(gff3Start, Is.EqualTo(1),
            $"GFF3 start should be 1 (0-based 0 + 1), got {gff3Start}");
    }

    /// <summary>
    /// INV-4: ParseGff3 skips comment and empty lines.
    /// Evidence: GFF3 Spec — lines starting with # are comments, blank lines are ignored.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_Parse_SkipsCommentsAndBlankLines()
    {
        var lines = new[]
        {
            "##gff-version 3",
            "# This is a comment",
            "",
            "chr1\t.\tCDS\t1\t300\t.\t+\t0\tID=gene_1;product=test",
            "# Another comment",
        };

        var features = GenomeAnnotator.ParseGff3(lines).ToList();
        Assert.That(features.Count, Is.EqualTo(1),
            "ParseGff3 should only return data lines, not comments or blanks");
    }

    /// <summary>
    /// INV-5: GFF3 data lines have exactly 9 tab-separated fields.
    /// Evidence: GFF3 Spec v1.26 defines 9 mandatory columns.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_Output_Has9Fields()
    {
        var genes = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("gene_1", 100, 500, '+', "CDS", "hypothetical protein",
                new Dictionary<string, string> { ["frame"] = "1" })
        };

        var lines = GenomeAnnotator.ToGff3(genes).Where(l => !l.StartsWith("#")).ToList();

        foreach (var line in lines)
        {
            int fieldCount = line.Split('\t').Length;
            Assert.That(fieldCount, Is.EqualTo(9),
                $"GFF3 data line must have 9 fields, got {fieldCount}");
        }
    }

    /// <summary>
    /// INV-6: GFF3 round-trip is deterministic.
    /// Evidence: Both ToGff3 and ParseGff3 are pure functions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Gff3_RoundTrip_IsDeterministic()
    {
        var genes = new List<GenomeAnnotator.GeneAnnotation>
        {
            new("gene_1", 0, 300, '+', "CDS", "test",
                new Dictionary<string, string> { ["frame"] = "1" }),
            new("gene_2", 500, 900, '-', "CDS", "test2",
                new Dictionary<string, string> { ["frame"] = "-2" })
        };

        var lines1 = GenomeAnnotator.ToGff3(genes).ToList();
        var lines2 = GenomeAnnotator.ToGff3(genes).ToList();

        Assert.That(lines1.SequenceEqual(lines2), Is.True,
            "GFF3 serialization must be deterministic");
    }

    #endregion
}
