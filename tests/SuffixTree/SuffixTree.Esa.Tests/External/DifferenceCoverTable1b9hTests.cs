using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The Singer cover of the default period (993, the nearest Singer period to the nominal 1024) is part of the code since work package 1b.9h (debt G(5): the search for a primitive cubic over GF(31) cost 0.27 s in every build, 11 %
/// of a build of 10^7 symbols): the table is the result of the search, it is checked every time it is used, and a new cover is made from it every time. What is pinned here: the table is what the search finds; a damaged table is
/// refused and the cover is then the one of the search; the covers of two calls share no array; parallel callers get the same cover; another period never gets this table.
/// </summary>
[TestFixture]
public class DifferenceCoverTable1b9hTests
{
    private static readonly int[] DefaultResidues = DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod).Residues;

    private static void AssertSameCover(DifferenceCover actual, DifferenceCover expected, string label)
    {
        Assert.That(actual.Period, Is.EqualTo(expected.Period), label);
        Assert.That(actual.Kind, Is.EqualTo(expected.Kind), label);
        Assert.That(actual.Residues, Is.EqualTo(expected.Residues), label);
        Assert.That(actual.IndexOf, Is.EqualTo(expected.IndexOf), label);
        Assert.That(actual.Witness, Is.EqualTo(expected.Witness), label);
    }

    [Test]
    public void TheCoverOfTheDefaultPeriod_IsWhatTheSearchFinds_ResidueForResidue_WithTheTablesOfTheComparator()
    {
        DifferenceCover table = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        DifferenceCover search = DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod);
        AssertSameCover(table, search, "993");
        Assert.That(table.Period, Is.EqualTo(993));
        Assert.That(table.Kind, Is.EqualTo(CoverKind.Singer));
        Assert.That(table.Residues, Has.Length.EqualTo(32), "q + 1 residues for q = 31");
        Assert.That(table.Verify(), Is.True);
        Assert.That(table.IsPerfect(), Is.True);
        Assert.That(table.Witness, Has.None.EqualTo(-1), "every difference has a witness");
        Assert.That(DifferenceCover.DefaultSingerPeriod, Is.EqualTo((31 * 31) + 31 + 1));
        Assert.That(DifferenceCover.NearestSingerPeriod(1024), Is.EqualTo(DifferenceCover.DefaultSingerPeriod), "993 is the Singer period that the builder takes for its default of 1024");
    }

    [Test]
    public void TheCoverThatTheBuilderTakesForItsDefaults_IsTheCoverOfTheDefaultPeriod()
    {
        var defaults = new SaBuildOptions { ScratchDir = "unused" };
        DifferenceCover cover = SaBuilder.CoverFor(defaults.Period, defaults.ExactPeriod);
        AssertSameCover(cover, DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod), "SaBuildOptions defaults");
        Assert.That(new SesaBuildOptions().Period, Is.EqualTo(defaults.Period), "the container builder has the same nominal period");
    }

    [Test]
    public void EveryCallGetsAFreshCover_AndNoCallerCanChangeTheCoverOfAnother()
    {
        DifferenceCover first = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        DifferenceCover second = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        Assert.That(second, Is.Not.SameAs(first));
        Assert.That(second.Residues, Is.Not.SameAs(first.Residues));
        Assert.That(second.IndexOf, Is.Not.SameAs(first.IndexOf));
        Assert.That(second.Witness, Is.Not.SameAs(first.Witness));
        Array.Fill(first.Residues, 7);
        Array.Fill(first.IndexOf, 7);
        Array.Fill(first.Witness, 7);
        DifferenceCover third = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        AssertSameCover(second, DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod), "the second cover after the first was overwritten");
        AssertSameCover(third, DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod), "a cover made after the first was overwritten");
        Assert.That(third.Verify(), Is.True);
    }

    [Test]
    public void ParallelCallers_GetTheSameCover_AndNeverShareAnObject()
    {
        const int callers = 16;
        const int rounds = 40;
        var covers = new DifferenceCover[callers * rounds];
        var errors = new List<Exception>();
        using var start = new ManualResetEventSlim();
        var threads = new Thread[callers];
        for (int t = 0; t < callers; t++)
        {
            int id = t;
            threads[t] = new Thread(() =>
            {
                try
                {
                    start.Wait();
                    for (int r = 0; r < rounds; r++)
                    {
                        DifferenceCover c = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
                        if (r % 3 == 0) Array.Clear(c.Residues); // a caller that damages its own cover, while the others read the table
                        covers[(id * rounds) + r] = c;
                    }
                }
#pragma warning disable CA1031 // collected and asserted below
                catch (Exception e)
                {
                    lock (errors) errors.Add(e);
                }
#pragma warning restore CA1031
            })
            { IsBackground = true };
            threads[t].Start();
        }
        start.Set();
        foreach (Thread th in threads) Assert.That(th.Join(TimeSpan.FromMinutes(2)), Is.True, "a caller did not return");
        Assert.That(errors, Is.Empty);
        Assert.That(covers.Select(c => c.Residues).Distinct().Count(), Is.EqualTo(callers * rounds), "no two callers share an array");
        int damaged = 0;
        foreach (DifferenceCover c in covers)
        {
            if (c.Residues.All(r => r == 0)) damaged++;
            else Assert.That(c.Residues, Is.EqualTo(DefaultResidues));
        }
        Assert.That(damaged, Is.EqualTo(callers * ((rounds + 2) / 3)), "only the covers whose callers damaged them are damaged");
        Assert.That(DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod).Residues, Is.EqualTo(DefaultResidues), "the table is untouched");
    }

    [Test]
    public void ADifferentPeriod_GivesADifferentCover_AndNeverTheTableOfTheDefault()
    {
        DifferenceCover defaults = DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        foreach (int period in new[] { 7, 13, 21, 31, 57, 553, 992, 994, 1024, 1893, 2000 })
        {
            DifferenceCover other = DifferenceCover.Create(period);
            Assert.That(other.Period, Is.EqualTo(period));
            Assert.That(other.Residues, Is.Not.EqualTo(defaults.Residues), $"X = {period}");
            Assert.That(other.Verify(), Is.True, $"X = {period}");
            AssertSameCover(other, DifferenceCover.CreateBySearch(period), $"X = {period} against the search");
        }
        Assert.That(DifferenceCover.Create(553).Kind, Is.EqualTo(CoverKind.Singer), "553 = 23^2 + 23 + 1");
        Assert.That(DifferenceCover.Create(553).Residues, Has.Length.EqualTo(24));
        Assert.That(DifferenceCover.Create(992).Kind, Is.EqualTo(CoverKind.Pruned), "the neighbours of 993 are not Singer periods and never get its table");
        Assert.That(DifferenceCover.Create(994).Kind, Is.EqualTo(CoverKind.Pruned));
    }

    [Test]
    public void TheDefaultPeriod_IsServedByTheTable_WithoutTheSearch_AndEveryOtherSingerPeriodByTheSearch()
    {
        _ = DefaultResidues;
        long before = DifferenceCover.Searches;
        for (int i = 0; i < 5; i++) DifferenceCover.Create(DifferenceCover.DefaultSingerPeriod);
        _ = SaBuilder.CoverFor(1024, false);
        _ = SaBuilder.CoverFor(DifferenceCover.DefaultSingerPeriod, true);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(before), "the default period is made without the search for a primitive cubic");
        _ = DifferenceCover.Create(553);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(before + 1), "another Singer period is searched");
        _ = DifferenceCover.Create(992);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(before + 1), "a period that is no Singer period is pruned, not searched");
        _ = DifferenceCover.CreateBySearch(DifferenceCover.DefaultSingerPeriod);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(before + 2), "CreateBySearch searches also for the default period");
    }

    [Test]
    public void ADamagedTableOfTheDefaultPeriod_IsRefused_TheSearchRuns_AndTheCoverIsTheRightOne()
    {
        int p = DifferenceCover.DefaultSingerPeriod;
        DifferenceCover expected = DifferenceCover.CreateBySearch(p);
        foreach (int[] damaged in new[] { Replace(DefaultResidues, 3, DefaultResidues[3] + 1), DefaultResidues[..^1], Array.Empty<int>(), Replace(DefaultResidues, 0, -4) })
        {
            long before = DifferenceCover.Searches;
            DifferenceCover cover = DifferenceCover.Create(p, damaged);
            Assert.That(DifferenceCover.Searches, Is.EqualTo(before + 1), "the damaged table is not used: the search runs");
            AssertSameCover(cover, expected, "a damaged table");
        }
        long untouched = DifferenceCover.Searches;
        DifferenceCover good = DifferenceCover.Create(p, DefaultResidues);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(untouched), "the right table is used");
        AssertSameCover(good, expected, "the table");
        long other = DifferenceCover.Searches;
        DifferenceCover seven = DifferenceCover.Create(7, DefaultResidues);   // the table is for the default period only
        Assert.That(seven.Period, Is.EqualTo(7));
        Assert.That(seven.Verify(), Is.True);
        Assert.That(DifferenceCover.Searches, Is.EqualTo(other + 1));
    }

    [Test]
    public void TheTableIsCheckedWhenItIsUsed_AndADamagedTableIsRefused()
    {
        int p = DifferenceCover.DefaultSingerPeriod;
        Assert.That(DifferenceCover.FromTable(p, DefaultResidues), Is.Not.Null, "the real table is accepted");
        var cases = new List<(string Name, int[] Table)>
        {
            ("one residue changed to a non-residue of a perfect set", Replace(DefaultResidues, 5, DefaultResidues[5] + 1)),
            ("one residue changed to another residue", Replace(DefaultResidues, 5, DefaultResidues[6])),
            ("a residue below zero", Replace(DefaultResidues, 0, -1)),
            ("a residue at the period", Replace(DefaultResidues, 31, p)),
            ("a residue above the period", Replace(DefaultResidues, 31, p + 5)),
            ("one residue short", DefaultResidues[..^1]),
            ("one residue more", [.. DefaultResidues, p - 1]),
            ("not sorted", Replace(Replace(DefaultResidues, 2, DefaultResidues[3]), 3, DefaultResidues[2])),
            ("the first two swapped (a perfect set, not in order)", Replace(Replace(DefaultResidues, 0, DefaultResidues[1]), 1, DefaultResidues[0])),
            ("the whole set reversed (a perfect set, not in order)", [.. DefaultResidues.Reverse()]),
            ("empty", []),
            ("the last residue replaced by one that is no residue: the right shape, not a perfect set", Replace(DefaultResidues, 31, Enumerable.Range(DefaultResidues[30] + 1, p - DefaultResidues[30] - 1).First(v => v != DefaultResidues[31]))),
        };
        foreach ((string name, int[] table) in cases) Assert.That(DifferenceCover.FromTable(p, table), Is.Null, name);
        Assert.That(DifferenceCover.FromTable(p + 1, DefaultResidues), Is.Null, "a period that is no Singer period takes no table");
        Assert.That(DifferenceCover.FromTable(553, DefaultResidues), Is.Null, "the table of one period is not the cover of another");
    }

    [Test]
    public void AShiftedPerfectSet_IsStillAPerfectSet_ButOnlyTheSearchsOwnSetIsTheTable()
    {
        // any translate of a perfect difference set is one: the check of the table cannot tell it from the real one, the comparison with the search can (the first test), and it is the set whose ranks the prototype's files hold
        int p = DifferenceCover.DefaultSingerPeriod;
        int[] shifted = [.. DefaultResidues.Select(r => (r + 1) % p).Order()];
        DifferenceCover? cover = DifferenceCover.FromTable(p, shifted);
        Assert.That(cover, Is.Not.Null);
        Assert.That(cover!.IsPerfect(), Is.True);
        Assert.That(cover.Residues, Is.Not.EqualTo(DefaultResidues));
        Assert.That(DifferenceCover.Create(p).Residues, Is.EqualTo(DefaultResidues), "Create never takes a shifted set");
    }

    private static int[] Replace(int[] source, int index, int value)
    {
        int[] copy = [.. source];
        copy[index] = value;
        return copy;
    }
}
