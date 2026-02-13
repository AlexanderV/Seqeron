namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for epigenetics analysis.
/// Verifies CpG island detection, methylation sites, and O/E ratio.
///
/// Test Unit: EPIGEN-CPG-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Annotation")]
public class EpigeneticsSnapshotTests
{
    // CpG-rich island region embedded in AT-rich flanking
    private const string TestSequence =
        "ATATATATATATATATATATATATATATAT" +
        "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" +
        "ATATATATATATATATATATATATATATAT";

    [Test]
    public Task FindCpGIslands_KnownSequence_MatchesSnapshot()
    {
        var islands = EpigeneticsAnalyzer.FindCpGIslands(TestSequence, minLength: 50)
            .Select(i => new
            {
                i.Start,
                i.End,
                Length = i.End - i.Start,
                GcContent = Math.Round(i.GcContent, 4),
                CpGRatio = Math.Round(i.CpGRatio, 4)
            })
            .ToList();

        return Verify(new { IslandCount = islands.Count, Islands = islands });
    }

    [Test]
    public Task FindCpGSites_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACGTCGACGCGTACG";
        var sites = EpigeneticsAnalyzer.FindCpGSites(seq).ToList();
        return Verify(new { SiteCount = sites.Count, Positions = sites });
    }

    [Test]
    public Task CpGObservedExpected_KnownSequence_MatchesSnapshot()
    {
        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(TestSequence);
        return Verify(new { OERatio = Math.Round(ratio, 4) });
    }

    [Test]
    public Task BisulfiteConversion_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACGTCGATCG";
        // No methylated positions — all C→U (shown as T in DNA context)
        string converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq);
        return Verify(new { Original = seq, Converted = converted });
    }
}
