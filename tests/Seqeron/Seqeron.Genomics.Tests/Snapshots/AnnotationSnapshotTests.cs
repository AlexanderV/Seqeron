namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for genome annotation algorithms.
/// Verifies ORF finding, gene prediction, and promoter motif detection output stability.
///
/// Test Units: ANNOT-ORF-001, ANNOT-GENE-001, ANNOT-PROM-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Annotation")]
public class AnnotationSnapshotTests
{
    // A synthetic gene-like sequence: promoter region + start codon + coding + stop codon
    private const string TestSequence =
        "TTGACAAAAAATTTTTTATAATAGCACGTACGATGATGAAAGCGTTCAAGCGTACTGCGTGA" +
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" +
        "ATGATGCCCAAAGGGTTTTAAACGTACGT";

    [Test]
    public Task FindOrfs_KnownSequence_MatchesSnapshot()
    {
        var orfs = GenomeAnnotator.FindOrfs(TestSequence, minLength: 9, searchBothStrands: false)
            .Select(o => new { o.Start, o.End, o.Frame, Length = o.End - o.Start })
            .OrderBy(o => o.Start)
            .Take(5)
            .ToList();

        return Verify(new { OrfCount = orfs.Count, TopOrfs = orfs });
    }

    [Test]
    public Task PredictGenes_KnownSequence_MatchesSnapshot()
    {
        var genes = GenomeAnnotator.PredictGenes(TestSequence, minOrfLength: 9, prefix: "gene")
            .Select(g => new { g.GeneId, g.Start, g.End, g.Strand })
            .OrderBy(g => g.Start)
            .Take(5)
            .ToList();

        return Verify(new { GeneCount = genes.Count, TopGenes = genes });
    }

    [Test]
    public Task FindPromoterMotifs_KnownSequence_MatchesSnapshot()
    {
        var motifs = GenomeAnnotator.FindPromoterMotifs(TestSequence)
            .Select(m => new { m.position, m.type, m.score })
            .OrderByDescending(m => m.score)
            .Take(5)
            .ToList();

        return Verify(new { MotifCount = motifs.Count, TopMotifs = motifs });
    }

    // NOTE: CodingPotential snapshot removed in ANNOT-CODING-001 — CalculateCodingPotential now
    // requires coding/non-coding hexamer frequency tables (CPAT, Wang et al. 2013); it is covered
    // by exact-value tests in GenomeAnnotator_CalculateCodingPotential_Tests.cs.
}
