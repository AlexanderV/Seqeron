namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot (golden-master) tests for pattern matching algorithms.
///
/// Test Units: PAT-EXACT-001, PAT-IUPAC-001, PAT-PWM-001 (Snapshot Extensions)
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Analysis")]
public class PatternMatchingSnapshotTests
{
    [Test]
    public Task FindExactMotif_Snapshot()
    {
        var dna = new DnaSequence("ACGTACGTACGTACGTACGT");
        var positions = MotifFinder.FindExactMotif(dna, "ACGT").ToList();

        return Verify(new { Positions = positions });
    }

    [Test]
    public Task FindDegenerateMotif_Snapshot()
    {
        var dna = new DnaSequence("ACGTACGTACGTACGTACGT");
        // R = A or G
        var matches = MotifFinder.FindDegenerateMotif(dna, "RCGT")
            .Select(m => new { m.Position, m.MatchedSequence })
            .ToList();

        return Verify(new { Matches = matches });
    }

    [Test]
    public Task CreatePwm_Snapshot()
    {
        var sequences = new[] { "ACGT", "ACGC", "ACGA", "ACGT" };
        var pwm = MotifFinder.CreatePwm(sequences);

        return Verify(new { PwmLength = pwm.Length });
    }
}
