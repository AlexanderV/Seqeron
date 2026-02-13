namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for molecular tools.
/// Verifies CRISPR guide design, primer design, restriction enzyme analysis.
///
/// Test Units: MOL-PRIMER-001, MOL-RESTR-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("MolTools")]
public class MolToolsSnapshotTests
{
    [Test]
    public Task CodonUsageComparison_TwoSequences_MatchesSnapshot()
    {
        string seq1 = "ATGAAAGCGTTCAAGCGTACTGCGTGA";
        string seq2 = "ATGAAGGCATTCAAACGCACAGCCTGA";

        double similarity = CodonOptimizer.CompareCodonUsage(seq1, seq2);
        var usage1 = CodonOptimizer.CalculateCodonUsage(seq1).OrderBy(kv => kv.Key).ToDictionary(kv => kv.Key, kv => kv.Value);
        var usage2 = CodonOptimizer.CalculateCodonUsage(seq2).OrderBy(kv => kv.Key).ToDictionary(kv => kv.Key, kv => kv.Value);

        return Verify(new
        {
            Similarity = Math.Round(similarity, 4),
            Usage1 = usage1,
            Usage2 = usage2
        });
    }

    [Test]
    public Task CreateCodonTable_FromSequence_MatchesSnapshot()
    {
        string refSeq = "ATGAAAGCGTTCAAGCGTACTGCGATGCCCAAAGGGTTTTAA";
        var table = CodonOptimizer.CreateCodonTableFromSequence(refSeq, "TestOrganism");

        var topCodons = table.CodonFrequencies
            .OrderByDescending(kv => kv.Value)
            .Take(10)
            .ToDictionary(kv => kv.Key, kv => Math.Round(kv.Value, 4));

        return Verify(new
        {
            table.OrganismName,
            TotalCodons = table.CodonFrequencies.Count,
            TopCodons = topCodons
        });
    }
}
