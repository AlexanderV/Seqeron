namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for restriction enzyme analysis.
///
/// Test Units: RESTR-FIND-001, RESTR-DIGEST-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("MolTools")]
public class RestrictionSnapshotTests
{
    private const string TestSequence =
        "AAAGAATTCAAAGGATCCAAAGAATTCAAAGGATCCAAA";

    [Test]
    public Task FindSites_EcoRI_Snapshot()
    {
        var dna = new DnaSequence(TestSequence);
        var sites = RestrictionAnalyzer.FindSites(dna, "EcoRI")
            .Select(s => new { s.Position, Enzyme = s.Enzyme.Name, s.RecognizedSequence, s.IsForwardStrand })
            .ToList();

        return Verify(new { Sites = sites });
    }

    [Test]
    public Task Digest_EcoRI_Snapshot()
    {
        var dna = new DnaSequence(TestSequence);
        var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI")
            .Select(f => new { f.FragmentNumber, f.Length, f.StartPosition, f.LeftEnzyme, f.RightEnzyme })
            .ToList();

        return Verify(new { Fragments = fragments });
    }

    [Test]
    public Task DigestSummary_MultiEnzyme_Snapshot()
    {
        var dna = new DnaSequence(TestSequence);
        var summary = RestrictionAnalyzer.GetDigestSummary(dna, "EcoRI", "BamHI");

        return Verify(new
        {
            summary.TotalFragments,
            summary.LargestFragment,
            summary.SmallestFragment,
            summary.AverageFragmentSize,
            summary.EnzymesUsed
        });
    }
}
