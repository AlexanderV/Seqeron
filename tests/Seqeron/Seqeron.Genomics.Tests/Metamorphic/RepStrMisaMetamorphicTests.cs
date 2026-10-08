namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the REP-STR-001 / REP-TANDEM-001 MISA additions (review B04, audit WP3):
/// strand symmetry of the MISA repeat-type class table, position shift of compound SSRs under an N prefix,
/// refinement of compounds when the maximal interruption grows, case invariance. Fixed seeds; every relation is exact.
/// Test Units: REP-STR-001, REP-TANDEM-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Repeats")]
public class RepStrMisaMetamorphicTests
{
    private static string RandomSsrRich(Random rng, int parts)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < parts; i++)
        {
            if (rng.Next(2) == 0)
            {
                int p = rng.Next(1, 7);
                string u = new(Enumerable.Range(0, p).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                sb.Append(string.Concat(Enumerable.Repeat(u, rng.Next(2, 13))));
            }
            else
            {
                sb.Append(new string(Enumerable.Range(0, rng.Next(0, 40)).Select(_ => "ACGTN"[rng.Next(rng.Next(3) == 0 ? 5 : 4)]).ToArray()));
            }
        }
        return sb.ToString();
    }

    private static string ReverseComplement(string s) =>
        new(s.Reverse().Select(c => "ACGTN"["TGCAN".IndexOf(c)]).ToArray());

    /// <summary>MR: maximal primitive runs of revcomp(S) are the mirrored runs of S (same period, same copies), and the
    /// MISA class of a run is strand-independent, so the class table (and the Krait level-2 table) is identical.</summary>
    [Test]
    public void CanonicalMotifTables_InvariantUnderReverseComplement()
    {
        var rng = new Random(3001);
        for (int t = 0; t < 200; t++)
        {
            string s = RandomSsrRich(rng, rng.Next(1, 12));
            var map = t % 2 == 0 ? RepeatFinder.MisaDefaultMinRepeats : new Dictionary<int, int> { [1] = 3, [2] = 2, [3] = 2, [5] = 2 };
            var a = RepeatFinder.FindMicrosatellites(s, map).ToList();
            var b = RepeatFinder.FindMicrosatellites(ReverseComplement(s), map).ToList();

            Assert.That(RepeatFinder.GetCanonicalMotifFrequencies(b), Is.EquivalentTo(RepeatFinder.GetCanonicalMotifFrequencies(a)), s);
            Assert.That(RepeatFinder.GetStandardMotifFrequencies(b), Is.EquivalentTo(RepeatFinder.GetStandardMotifFrequencies(a)), s);
            Assert.That(b.Select(m => (m.RepeatUnit.Length, m.RepeatCount)), Is.EquivalentTo(a.Select(m => (m.RepeatUnit.Length, m.RepeatCount))), s);
        }
    }

    /// <summary>MR: prefixing k N's (never part of an SSR) shifts every compound by k and leaves notation/type unchanged.</summary>
    [Test]
    public void CompoundMicrosatellites_NPrefix_ShiftsPositions()
    {
        var rng = new Random(3002);
        for (int t = 0; t < 200; t++)
        {
            string s = RandomSsrRich(rng, rng.Next(1, 12));
            int k = rng.Next(1, 30);
            int amb = new[] { 0, 5, 20, 100 }[t % 4];
            var a = RepeatFinder.FindCompoundMicrosatellites(s, 3, amb);
            var b = RepeatFinder.FindCompoundMicrosatellites(new string('N', k) + s, 3, amb);

            Assert.That(b.Select(c => (c.Start - k, c.End - k, c.MisaType, c.Notation)),
                Is.EqualTo(a.Select(c => (c.Start, c.End, c.MisaType, c.Notation))), s);
        }
    }

    /// <summary>MR: raising maxInterruption only adds joins between the same consecutive SSRs, so every compound (and
    /// every single SSR) at a smaller threshold lies inside one compound at the larger threshold.</summary>
    [Test]
    public void CompoundMicrosatellites_LargerInterruption_CoarsensPartition()
    {
        var rng = new Random(3003);
        for (int t = 0; t < 200; t++)
        {
            string s = RandomSsrRich(rng, rng.Next(2, 14));
            var small = RepeatFinder.FindCompoundMicrosatellites(s, 3, 5);
            var large = RepeatFinder.FindCompoundMicrosatellites(s, 3, 60);

            foreach (var c in small)
            {
                var keys = c.Components.Select(m => (m.Position, m.RepeatUnit)).ToHashSet();
                Assert.That(large.Count(l => keys.IsSubsetOf(l.Components.Select(m => (m.Position, m.RepeatUnit)))), Is.EqualTo(1), s);
            }
            Assert.That(large.Sum(c => c.Components.Count), Is.GreaterThanOrEqualTo(small.Sum(c => c.Components.Count)), s);
        }
    }

    /// <summary>MR: lower-casing the input changes nothing (interruptions are lower-cased as in misa.pl).</summary>
    [Test]
    public void CompoundMicrosatellites_CaseInvariant()
    {
        var rng = new Random(3004);
        for (int t = 0; t < 200; t++)
        {
            string s = RandomSsrRich(rng, rng.Next(1, 12));
            var a = RepeatFinder.FindCompoundMicrosatellites(s);
            var b = RepeatFinder.FindCompoundMicrosatellites(s.ToLowerInvariant());
            Assert.That(b.Select(c => (c.Start, c.End, c.Notation)), Is.EqualTo(a.Select(c => (c.Start, c.End, c.Notation))), s);
        }
    }
}
