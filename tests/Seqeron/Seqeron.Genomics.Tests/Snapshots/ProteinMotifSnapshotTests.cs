namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for protein motif finding.
///
/// Test Units: PROTMOTIF-FIND-001, PROTMOTIF-PROSITE-001, PROTMOTIF-DOMAIN-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Analysis")]
public class ProteinMotifSnapshotTests
{
    private const string TestProtein = "MKTLLLTLVVVTLVLSSQPVLSRELRECPRGSGKSCQACPAG" +
                                       "NISTYQCQSYVMSHLCSYQCNQRCFQSLENQCQTFHCRGFQF" +
                                       "NSTRTMPLHCRGFQFNSTRTMPLHCRG";

    [Test]
    public Task FindCommonMotifs_Snapshot()
    {
        var motifs = ProteinMotifFinder.FindCommonMotifs(TestProtein)
            .Select(m => new { m.Start, m.End, m.Sequence, m.MotifName })
            .ToList();

        return Verify(new { Motifs = motifs });
    }

    [Test]
    public Task FindMotifByProsite_NGlycosylation_Snapshot()
    {
        string prositePattern = "N-{P}-[ST]-{P}";
        var motifs = ProteinMotifFinder.FindMotifByProsite(TestProtein, prositePattern)
            .Select(m => new { m.Start, m.End, m.Sequence, m.MotifName })
            .ToList();

        return Verify(new { Motifs = motifs });
    }

    [Test]
    public Task PredictSignalPeptide_Snapshot()
    {
        var sp = ProteinMotifFinder.PredictSignalPeptide(TestProtein);

        return Verify(new
        {
            HasSignalPeptide = sp.HasValue,
            sp?.CleavagePosition,
            // Round to 6 dp: the score's last (~15th) double digit is not meaningful precision and
            // varies with FMA/vectorization across platforms, which would make the snapshot flaky.
            Score = sp.HasValue ? Math.Round(sp.Value.Score, 6) : (double?)null,
            sp?.IsLikelySignalPeptide
        });
    }
}
