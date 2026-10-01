namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-DIRECT-001 Vmatch output modes (B04 completeness audit WP15):
/// <see cref="RepeatFinder.FindDegenerateRepeats(string,int,int,ApproximateRepeatDistance,bool,int,int,DegenerateRepeatReporting,bool)"/>
/// with <see cref="DegenerateRepeatReporting.BestPerSeed"/> (Vmatch's default output without <c>-allmax</c>: one
/// best extension per exact seed by E-value, identity, length) and <c>vmatchCompatible</c> (stock Vmatch 2.3.1
/// left-extension seed shortcut). Expected lists are real Vmatch output: the Debian <c>vmatch</c> binary (stock) and
/// the same release built from source with the shortcut switched off (env <c>VM_NOPRUNE</c> / <c>VM_NOPRUNE_H</c>);
/// tuples (FirstPosition, SecondPosition, FirstLength, SecondLength, Distance) = Vmatch's (i, j, l, r, |d|), sorted.
/// </summary>
[TestFixture]
public class RepeatFinder_VmatchReporting_Tests
{
    // 121 bp: a 30-mer copied with one substitution, then with a deletion (Evidence REP-DIRECT-001 §WP15).
    private const string Copies =
        "CCGGCCCCTGAGTCCGAGGAGGATCACAGTCTACACTGCTCACTCCAACCGAGGATCACAGTTTACACTGCTCACTCCAACCGAGGGTGCTTGGATCACAGTCTACATGCTCACTCCAACC";

    private const string ShortcutCase = "ATCTGGTGTACTCTGCCCACGACTATCGGTGTACTCTGC";

    private static string Run(string seq, int minLength, int k, ApproximateRepeatDistance distance, bool palindromic,
        DegenerateRepeatReporting reporting, bool compatible) =>
        string.Join(";", RepeatFinder.FindDegenerateRepeats(seq, minLength, k, distance, palindromic, int.MaxValue, int.MinValue,
                reporting, compatible)
            .Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.FirstLength},{x.SecondLength},{x.Distance}"));

    private const ApproximateRepeatDistance H = ApproximateRepeatDistance.Hamming;
    private const ApproximateRepeatDistance E = ApproximateRepeatDistance.Edit;
    private const DegenerateRepeatReporting Best = DegenerateRepeatReporting.BestPerSeed;
    private const DegenerateRepeatReporting All = DegenerateRepeatReporting.AllMaximal;

    #region BestPerSeed (vmatch without -allmax)

    /// <summary>vmatch -l 12 -h 2: stock keeps exact seeds whose left scan meets another seed; complete extension
    /// lets every seed reach the E-value-best window (so the same row appears once per seed).</summary>
    [Test]
    public void BestPerSeed_HammingDirect_MatchesVmatch()
    {
        Assert.That(Run(Copies, 12, 2, H, false, Best, compatible: true),
            Is.EqualTo("18,50,36,36,1;20,92,15,15,0;31,63,23,23,0;36,107,14,14,0;52,92,15,15,1;68,107,14,14,0"));
        Assert.That(Run(Copies, 12, 2, H, false, Best, compatible: false),
            Is.EqualTo("18,50,36,36,1;18,50,36,36,1;20,92,15,15,0;36,107,14,14,0;52,92,15,15,1;52,92,15,15,1;68,107,14,14,0"));
    }

    /// <summary>vmatch -l 12 -e 2 (stock / shortcut off).</summary>
    [Test]
    public void BestPerSeed_EditDirect_MatchesVmatch()
    {
        Assert.That(Run(Copies, 12, 2, E, false, Best, compatible: true),
            Is.EqualTo("18,50,36,36,1;20,92,30,29,1;31,63,23,23,0;36,107,14,14,0;52,92,30,29,2;63,103,19,18,1;68,107,14,14,0"));
        Assert.That(Run(Copies, 12, 2, E, false, Best, compatible: false),
            Is.EqualTo("18,50,36,36,1;18,50,36,36,1;20,92,30,29,1;20,92,30,29,1;52,92,30,29,2;52,92,30,29,2;52,92,30,29,2"));
    }

    /// <summary>vmatch -p -l 12 -h 2 and vmatch -p -l 5 -e 3 (palindromic; rows with i &gt; j dropped as in Vmatch).</summary>
    [Test]
    public void BestPerSeed_Palindromic_MatchesVmatch()
    {
        Assert.That(Run(Copies, 12, 2, H, true, Best, compatible: true), Is.EqualTo("26,26,12,12,2;26,58,12,12,2;58,58,12,12,2"));
        Assert.That(Run(Copies, 12, 2, H, true, Best, compatible: false),
            Is.EqualTo("26,26,12,12,2;26,26,12,12,2;26,58,12,12,2;26,58,12,12,2;58,58,12,12,2;58,58,12,12,2"));
        Assert.That(Run("GGACCATGAAGG", 5, 3, E, true, Best, compatible: true),
            Is.EqualTo("0,0,5,5,1;0,0,5,5,3;0,0,5,6,3;2,2,5,5,3;2,4,6,5,2;3,3,9,9,3;3,4,5,5,1;3,4,5,7,2;4,6,5,6,3"));
        Assert.That(Run("GGACCATGAAGG", 5, 3, E, true, Best, compatible: false),
            Is.EqualTo("0,0,5,5,1;0,0,5,5,1;0,0,5,5,3;0,0,5,5,3;0,0,5,6,3;0,0,6,5,3;2,2,5,5,3;2,2,5,5,3;2,4,6,5,2;3,3,5,6,2;" +
                       "3,3,9,9,3;3,3,9,9,3;3,3,9,9,3;3,3,9,9,3;3,4,5,5,1;3,4,5,7,2;3,6,7,6,3;3,6,7,6,3;4,6,5,6,3"));
    }

    /// <summary>
    /// Vmatch cmpmatches prefers the smaller E-value even over a longer window: vmatch -l 15 -e 3 reports the
    /// distance-1 window (0,24,16,15,1) once per seed, not the maximal distance-3 matches of -allmax; -h 3 gives
    /// (1,24,15,15,3).
    /// </summary>
    [Test]
    public void BestPerSeed_PrefersSmallerEvalueOverLongerMatch()
    {
        Assert.That(Run(ShortcutCase, 15, 3, E, false, Best, compatible: true), Is.EqualTo("0,24,16,15,1"));
        Assert.That(Run(ShortcutCase, 15, 3, E, false, Best, compatible: false), Is.EqualTo("0,24,16,15,1;0,24,16,15,1"));
        Assert.That(Run(ShortcutCase, 15, 3, H, false, Best, compatible: true), Is.EqualTo("1,24,15,15,3"));
        Assert.That(Run(ShortcutCase, 15, 3, H, false, Best, compatible: false), Is.EqualTo("1,24,15,15,3"));
    }

    /// <summary>Every best-per-seed row is a real match: length ≥ minLength, Hamming distance = Distance, copies of the row.</summary>
    [Test]
    public void BestPerSeed_RowsAreValidMatches()
    {
        var rng = new Random(15);
        for (int t = 0; t < 40; t++)
        {
            string s = RandomRepeatRich(rng, rng.Next(30, 120));
            foreach (bool compatible in new[] { false, true })
            foreach (var r in RepeatFinder.FindDegenerateRepeats(s, 10, 2, H, false, int.MaxValue, int.MinValue, Best, compatible))
            {
                Assert.That(r.FirstLength, Is.EqualTo(r.SecondLength));
                Assert.That(r.FirstLength, Is.GreaterThanOrEqualTo(10));
                Assert.That(r.FirstCopy, Is.EqualTo(s.Substring(r.FirstPosition, r.FirstLength)));
                int mm = Enumerable.Range(0, r.FirstLength).Count(x => r.FirstCopy[x] != r.SecondCopy[x]);
                Assert.That(mm, Is.EqualTo(r.Distance));
            }
        }
    }

    #endregion

    #region vmatchCompatible with -allmax

    /// <summary>
    /// Stock vmatch -l 15 -e 3 -allmax: the shortcut loses (0,21,16,18,3) and reports the contained (0,22,16,17,3);
    /// the default stays the complete set (= shortcut off = brute force, REP-DIRECT-001 §WP8).
    /// </summary>
    [Test]
    public void AllMaximal_Compatible_ReproducesStockShortcutLoss()
    {
        Assert.That(Run(ShortcutCase, 15, 3, E, false, All, compatible: true), Is.EqualTo("0,22,16,17,3;0,23,17,16,3;0,24,18,15,3"));
        Assert.That(Run(ShortcutCase, 15, 3, E, false, All, compatible: false), Is.EqualTo("0,21,16,18,3;0,23,17,16,3;0,24,18,15,3"));
        Assert.That(Run("GGACCATGAAGG", 5, 3, E, true, All, compatible: true),
            Is.EqualTo("0,0,5,7,3;0,0,6,6,3;0,0,7,5,3;1,4,7,6,3;2,2,5,5,3;2,3,7,7,3;2,4,6,8,3;2,6,7,6,3;3,3,9,9,3"));
        Assert.That(Run("GGACCATGAAGG", 5, 3, E, true, All, compatible: false),
            Is.EqualTo("0,0,5,7,3;0,0,6,6,3;0,0,7,5,3;1,4,7,7,3;2,2,5,5,3;2,3,7,7,3;2,4,6,8,3;2,6,7,6,3;3,3,9,9,3"));
    }

    /// <summary>
    /// Vmatch's matchcontainer keeps the distance of the first seed that reaches a match (not the edit distance):
    /// vmatch -p -l 8 -e 3 -allmax prints (0,0,8,8,3) for CAAATTTT (d_E(CAAATTTT, AAAATTTG) = 2), and
    /// vmatch -l 7 -e 3 -allmax prints (0,2,11,11,3) for ACAAAAAAAACAC (d_E = 2). The default reports the distance.
    /// </summary>
    [Test]
    public void AllMaximal_Compatible_ReproducesFirstSeedDistanceLabel()
    {
        Assert.That(Run("CAAATTTT", 8, 3, E, true, All, compatible: true), Is.EqualTo("0,0,8,8,3"));
        Assert.That(Run("CAAATTTT", 8, 3, E, true, All, compatible: false), Is.EqualTo("0,0,8,8,2"));
        Assert.That(Run("ACAAAAAAAACAC", 7, 3, E, false, All, compatible: true), Is.EqualTo("0,1,9,11,3;0,2,11,11,3;0,3,12,10,3"));
        Assert.That(Run("ACAAAAAAAACAC", 7, 3, E, false, All, compatible: false), Is.EqualTo("0,1,9,11,3;0,2,11,11,2;0,3,12,10,3"));
    }

    /// <summary>For Hamming distance the shortcut never changes the maximal set (stock = definition, §WP8).</summary>
    [Test]
    public void AllMaximal_Hamming_CompatibleEqualsDefault()
    {
        var rng = new Random(7);
        for (int t = 0; t < 40; t++)
        {
            string s = RandomRepeatRich(rng, rng.Next(20, 100));
            foreach (bool pal in new[] { false, true })
                Assert.That(Run(s, 8, 2, H, pal, All, compatible: true), Is.EqualTo(Run(s, 8, 2, H, pal, All, compatible: false)), s);
        }
    }

    #endregion

    #region API

    [Test]
    public void ShortOverload_EqualsAllMaximalDefault_AndDnaOverloadEqualsString()
    {
        string legacy = string.Join(";", RepeatFinder.FindDegenerateRepeats(ShortcutCase, 15, 3, E, false, int.MaxValue, int.MinValue)
            .Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.FirstLength},{x.SecondLength},{x.Distance}"));
        Assert.That(Run(ShortcutCase, 15, 3, E, false, All, compatible: false), Is.EqualTo(legacy));

        var viaDna = RepeatFinder.FindDegenerateRepeats(new DnaSequence(Copies), 12, 2, E, false, int.MaxValue, int.MinValue, Best, true);
        var viaString = RepeatFinder.FindDegenerateRepeats(Copies, 12, 2, E, false, int.MaxValue, int.MinValue, Best, true);
        Assert.That(viaDna, Is.EqualTo(viaString));
        Assert.That(RepeatFinder.FindDegenerateRepeats("", 12, 2, E, false, int.MaxValue, 1, Best), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindDegenerateRepeats(Copies, 12, 2, E, false, int.MaxValue, 1, (DegenerateRepeatReporting)7));
    }

    [Test]
    public void ApproximateDirectRepeats_BestPerSeed_IsHammingDirectBestPerSeed()
    {
        var a = RepeatFinder.FindApproximateDirectRepeats(Copies, 12, 2, int.MaxValue, int.MinValue, false, Best, vmatchCompatible: true).ToList();
        Assert.That(string.Join(";", a.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.Length},{x.Mismatches}")),
            Is.EqualTo("18,50,36,1;20,92,15,0;31,63,23,0;36,107,14,0;52,92,15,1;68,107,14,0"));
        Assert.That(a[0].Spacing, Is.EqualTo(50 - 18 - 36));
        Assert.That(RepeatFinder.FindApproximateDirectRepeats(Copies, 12, 2, int.MaxValue, 1, true, All),
            Is.EqualTo(RepeatFinder.FindApproximateDirectRepeats(Copies, 12, 2, int.MaxValue, 1, true)));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindApproximateDirectRepeats(Copies, 12, 0, int.MaxValue, 1, false, Best));
    }

    #endregion

    private static string RandomRepeatRich(Random rng, int n)
    {
        string motif = new(Enumerable.Range(0, rng.Next(6, 16)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var sb = new System.Text.StringBuilder();
        while (sb.Length < n)
        {
            foreach (char c in motif) sb.Append(rng.Next(10) == 0 ? "ACGT"[rng.Next(4)] : c);
            for (int x = rng.Next(0, 5); x > 0; x--) sb.Append("ACGT"[rng.Next(4)]);
        }
        return sb.ToString(0, n);
    }
}
