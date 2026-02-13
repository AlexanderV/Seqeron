namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for repeat finding algorithms.
/// Verifies palindrome, inverted repeat, direct repeat, and microsatellite outputs.
///
/// Test Units: REP-PALIN-001, REP-INV-001, REP-DIRECT-001, REP-STR-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Analysis")]
public class RepeatSnapshotTests
{
    [Test]
    public Task FindPalindromes_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACGTGAATTCACGTACGTGATATCACGT";
        var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12)
            .Select(p => new { p.Position, p.Sequence, p.Length })
            .OrderBy(p => p.Position)
            .ToList();

        return Verify(new { Count = palindromes.Count, Palindromes = palindromes });
    }

    [Test]
    public Task FindDirectRepeats_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACGTACGTTTTTTTTTACGTACGT";
        var repeats = RepeatFinder.FindDirectRepeats(seq, minLength: 4)
            .Select(r => new { r.FirstPosition, r.SecondPosition, r.RepeatSequence, r.Length, r.Spacing })
            .OrderBy(r => r.FirstPosition)
            .ToList();

        return Verify(new { Count = repeats.Count, Repeats = repeats });
    }

    [Test]
    public Task FindMicrosatellites_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACACACACACACGTGTGTGTGTAAAAAAAAAA";
        var microsats = RepeatFinder.FindMicrosatellites(seq, minUnitLength: 1, maxUnitLength: 4, minRepeats: 3)
            .Select(m => new { m.Position, m.RepeatUnit, m.RepeatCount, m.TotalLength })
            .OrderBy(m => m.Position)
            .ToList();

        return Verify(new { Count = microsats.Count, Microsatellites = microsats });
    }

    [Test]
    public Task FindInvertedRepeats_KnownSequence_MatchesSnapshot()
    {
        string seq = "ACGTACGTAAAAAAATGCATGCA";
        var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4)
            .Select(r => new
            {
                r.LeftArmStart,
                r.RightArmStart,
                r.ArmLength,
                r.LoopLength,
                r.LeftArm,
                r.RightArm,
                r.CanFormHairpin
            })
            .ToList();

        return Verify(new { Count = repeats.Count, Repeats = repeats });
    }
}
