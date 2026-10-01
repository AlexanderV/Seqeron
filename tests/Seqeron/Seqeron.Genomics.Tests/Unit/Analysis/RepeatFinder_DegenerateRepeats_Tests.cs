namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-DIRECT-001 degenerate repeats (B04 completeness audit WP8, L11):
/// <see cref="RepeatFinder.FindDegenerateRepeats(string,int,int,ApproximateRepeatDistance,bool,int,int)"/> —
/// maximal k-differences (unit-cost edit distance) and k-mismatch repeats, direct and palindromic, with the
/// Vmatch manual App. A definitions (Kurtz et al. 2001 REPuter; <c>vmatch [-p] -l m (-e|-h) k -allmax</c>).
/// Expected lists: Vmatch 2.3.1 (Debian binary and the same release built from source with the left-extension
/// seed shortcut disabled, <c>VM_NOPRUNE</c>) and an independent brute force of the definition; tuples are
/// (FirstPosition, SecondPosition, FirstLength, SecondLength, Distance), i.e. Vmatch's (i, j, l, r, distance).
/// </summary>
[TestFixture]
public class RepeatFinder_DegenerateRepeats_Tests
{
    private static List<DegenerateRepeatResult> Find(
        string seq, int k, int minLength, ApproximateRepeatDistance distance, bool palindromic) =>
        RepeatFinder.FindDegenerateRepeats(seq, minLength, k, distance, palindromic, int.MaxValue, int.MinValue).ToList();

    private static string Tuples(IEnumerable<DegenerateRepeatResult> r) =>
        string.Join(";", r.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.FirstLength},{x.SecondLength},{x.Distance}"));

    private static string RevComp(string s) =>
        new(s.Reverse().Select(SequenceExtensions.GetComplementBase).ToArray());

    #region Vmatch-locked lists

    /// <summary>vmatch -l 8 -e 1 -allmax (stock = shortcut disabled = brute force).</summary>
    [Test]
    public void EditDirect_MatchesVmatchAllmax()
    {
        var r = Find("ACGTTGCATGCAAACGTAGCATGCAGGGTTTACGTTGCTTGCAAACG", 1, 8, ApproximateRepeatDistance.Edit, false);
        Assert.That(Tuples(r), Is.EqualTo("0,13,12,12,1;0,31,16,16,1;5,18,8,8,1;8,39,9,8,1"));
        Assert.That(r[3].FirstCopy, Is.EqualTo("TGCAAACGT"));
        Assert.That(r[3].SecondCopy, Is.EqualTo("TGCAAACG"));
        Assert.That(r[3].Spacing, Is.EqualTo(39 - 8 - 9));
        Assert.That(r.All(x => !x.IsReverseComplement));
    }

    /// <summary>Insertions and deletions: GATTACA … GATTCAGATTACA (vmatch -l 6 -e 1 -allmax).</summary>
    [Test]
    public void EditDirect_Indels_MatchesVmatchAllmax()
    {
        var r = Find("GATTACAGATTACATTTTTGATTCAGATTACA", 1, 6, ApproximateRepeatDistance.Edit, false);
        Assert.That(Tuples(r), Is.EqualTo("0,6,7,8,1;0,7,8,8,1;0,19,14,13,1;5,23,10,9,1;19,25,6,7,1"));
    }

    /// <summary>
    /// The Vmatch left-extension shortcut (stop where another exact match ≥ the seed length starts) loses
    /// (0, 21, 16, 18, 3) here and reports the contained (0, 22, 16, 17, 3) instead; with the shortcut disabled Vmatch
    /// and the brute force give the complete set returned by this method.
    /// </summary>
    [Test]
    public void EditDirect_CompleteWhereVmatchShortcutMissesAMatch()
    {
        var r = Find("ATCTGGTGTACTCTGCCCACGACTATCGGTGTACTCTGC", 3, 15, ApproximateRepeatDistance.Edit, false);
        Assert.That(Tuples(r), Is.EqualTo("0,21,16,18,3;0,23,17,16,3;0,24,18,15,3"));
    }

    /// <summary>
    /// Vmatch acceptmatch: trivial self-alignments of a tandem array (non-overlapping part ≤ distance, or the right
    /// instance embedded in the left one) are not matches; (AC)9 -l 4 -e 1 gives exactly two shifted copies.
    /// </summary>
    [Test]
    public void EditDirect_TandemArray_AcceptRuleAsVmatch()
    {
        var r = Find("ACACACACACACACACAC", 1, 4, ApproximateRepeatDistance.Edit, false);
        Assert.That(Tuples(r), Is.EqualTo("0,1,16,17,1;0,2,17,16,1"));
    }

    /// <summary>A wildcard never matches but can be inserted/deleted (vmatch -l 6 -e 1 -allmax).</summary>
    [Test]
    public void EditDirect_WildcardIsAnEditedSymbol()
    {
        var r = Find("GATTACANGATTACA", 1, 6, ApproximateRepeatDistance.Edit, false);
        Assert.That(Tuples(r), Is.EqualTo("0,7,7,8,1;0,8,8,7,1"));
        Assert.That(Find("gattacangattaca", 1, 6, ApproximateRepeatDistance.Edit, false).Select(x => (x.FirstPosition, x.SecondPosition)),
            Is.EqualTo(r.Select(x => (x.FirstPosition, x.SecondPosition))));
    }

    /// <summary>
    /// vmatch -p -l 5 -e 3 -allmax with the shortcut disabled (= brute force): palindromic matches incl. i = j with
    /// l ≠ r in both orientations; stock Vmatch reports the contained (1, 4, 7, 6, 3) instead of (1, 4, 7, 7, 3).
    /// </summary>
    [Test]
    public void EditPalindromic_MatchesVmatchAllmax()
    {
        var r = Find("GGACCATGAAGG", 3, 5, ApproximateRepeatDistance.Edit, true);
        Assert.That(Tuples(r), Is.EqualTo(
            "0,0,5,7,3;0,0,6,6,3;0,0,7,5,3;1,4,7,7,3;2,2,5,5,3;2,3,7,7,3;2,4,6,8,3;2,6,7,6,3;3,3,9,9,3"));
        Assert.That(r.All(x => x.IsReverseComplement && x.FirstPosition <= x.SecondPosition));
    }

    /// <summary>A hairpin with one inserted base: vmatch -p -l 8 -e 1 vs -p -l 8 -h 1.</summary>
    [Test]
    public void Palindromic_EditAndHamming_MatchVmatchAllmax()
    {
        const string seq = "TTGACCGTAACCCCCGTTACGGTCAACC";
        Assert.That(Tuples(Find(seq, 1, 8, ApproximateRepeatDistance.Edit, true)),
            Is.EqualTo("0,14,12,12,1;0,15,11,12,1;13,13,9,9,1"));
        Assert.That(Tuples(Find(seq, 1, 8, ApproximateRepeatDistance.Hamming, true)),
            Is.EqualTo("0,14,12,12,1;13,13,9,9,1"));
    }

    /// <summary>Direct Hamming via this method = FindApproximateDirectRepeats(excludeContained: true) (vmatch -h k -allmax).</summary>
    [Test]
    public void HammingDirect_EqualsApproximateDirectRepeatsExcludeContained()
    {
        var rng = new Random(8);
        for (int t = 0; t < 40; t++)
        {
            string seq = new(Enumerable.Range(0, 120).Select(_ => "ACGT"[rng.Next(t % 2 == 0 ? 4 : 2)]).ToArray());
            int k = 1 + t % 3;
            var expected = RepeatFinder.FindApproximateDirectRepeats(seq, 8, k, int.MaxValue, int.MinValue, excludeContained: true)
                .Select(x => (x.FirstPosition, x.SecondPosition, x.Length, x.Length, x.Mismatches));
            var actual = Find(seq, k, 8, ApproximateRepeatDistance.Hamming, false)
                .Select(x => (x.FirstPosition, x.SecondPosition, x.FirstLength, x.SecondLength, x.Distance));
            Assert.That(actual, Is.EqualTo(expected), seq);
        }
    }

    #endregion

    #region Definition (brute force)

    /// <summary>
    /// Brute force of the Vmatch App. A definitions on 160 random / repeat-rich inputs (all four modes): every
    /// k-differences (or k-mismatch) match, acceptmatch for direct edit matches, palindromic containment over both
    /// orientations, maximal = not contained, both instances ≥ minLength.
    /// </summary>
    [Test]
    public void AllModes_EqualBruteForceOfDefinition()
    {
        var rng = new Random(2026);
        for (int t = 0; t < 160; t++)
        {
            string seq = RandomRepeatRich(rng, rng.Next(8, 26));
            int k = 1 + rng.Next(2);
            int minLength = k + 1 + rng.Next(5);
            var distance = t % 2 == 0 ? ApproximateRepeatDistance.Edit : ApproximateRepeatDistance.Hamming;
            bool palindromic = t % 4 >= 2;
            var expected = BruteForce(seq, k, minLength, distance == ApproximateRepeatDistance.Edit, palindromic);
            var actual = Find(seq, k, minLength, distance, palindromic)
                .Select(x => (x.FirstPosition, x.SecondPosition, x.FirstLength, x.SecondLength, x.Distance)).ToList();
            Assert.That(actual, Is.EquivalentTo(expected), $"{seq} k={k} m={minLength} {distance} pal={palindromic}");
        }
    }

    /// <summary>Every reported repeat has the stated distance between its copies (canonical edit / Hamming distance).</summary>
    [Test]
    public void Results_DistanceAndCopiesAreConsistent()
    {
        var rng = new Random(77);
        string seq = RandomRepeatRich(rng, 400).Replace('N', 'A'); // canonical distances treat N = N
        foreach (bool pal in new[] { false, true })
        {
            foreach (var x in Find(seq, 2, 10, ApproximateRepeatDistance.Edit, pal))
            {
                string second = pal ? RevComp(x.SecondCopy) : x.SecondCopy;
                Assert.That(ApproximateMatcher.EditDistance(x.FirstCopy, second), Is.EqualTo(x.Distance), $"{seq} {x}");
                Assert.That(x.FirstCopy, Is.EqualTo(seq.Substring(x.FirstPosition, x.FirstLength)));
                Assert.That(x.SecondCopy, Is.EqualTo(seq.Substring(x.SecondPosition, x.SecondLength)));
            }
            foreach (var x in Find(seq, 2, 10, ApproximateRepeatDistance.Hamming, pal))
            {
                string second = pal ? RevComp(x.SecondCopy) : x.SecondCopy;
                Assert.That(x.FirstLength, Is.EqualTo(x.SecondLength));
                Assert.That(ApproximateMatcher.HammingDistance(x.FirstCopy, second), Is.EqualTo(x.Distance));
            }
        }
    }

    #endregion

    #region Filters, overloads, validation

    [Test]
    public void Filters_AppliedAfterMaximality()
    {
        const string seq = "GATTACAGATTACATTTTTGATTCAGATTACA";
        var all = Find(seq, 1, 6, ApproximateRepeatDistance.Edit, false);
        var spaced = RepeatFinder.FindDegenerateRepeats(seq, 6, 1).ToList(); // defaults: Edit, direct, minSpacing 1
        var shortOnly = RepeatFinder.FindDegenerateRepeats(seq, 6, 1, ApproximateRepeatDistance.Edit, false, 9, int.MinValue).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(spaced, Is.EqualTo(all.Where(x => x.Spacing >= 1)));
            Assert.That(Tuples(spaced), Is.EqualTo("0,19,14,13,1;5,23,10,9,1"));
            Assert.That(shortOnly, Is.EqualTo(all.Where(x => x.FirstLength <= 9 && x.SecondLength <= 9)));
        });
    }

    [Test]
    public void DnaOverloadEqualsString_EmptyInput_AndValidation()
    {
        const string seq = "ACGTTGCATGCAAACGTAGCATGCAGGGTTTACGTTGCTTGCAAACG";
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindDegenerateRepeats(new DnaSequence(seq), 8, 1, ApproximateRepeatDistance.Edit, true),
                Is.EqualTo(RepeatFinder.FindDegenerateRepeats(seq, 8, 1, ApproximateRepeatDistance.Edit, true)));
            Assert.That(RepeatFinder.FindDegenerateRepeats(string.Empty), Is.Empty);
            Assert.That(RepeatFinder.FindDegenerateRepeats((string)null!), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindDegenerateRepeats((DnaSequence)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindDegenerateRepeats(seq, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindDegenerateRepeats(seq, 8, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindDegenerateRepeats(seq, 8, 8));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindDegenerateRepeats(seq, 8, 1, maxLength: 7));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindDegenerateRepeats(seq, 8, 1, (ApproximateRepeatDistance)7));
        });
    }

    #endregion

    #region Brute force

    private static string RandomRepeatRich(Random rng, int n)
    {
        string alphabet = rng.Next(3) switch { 0 => "ACGT", 1 => "AC", _ => "ACGTN" };
        if (rng.Next(2) == 0)
            return new string(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
        string motif = new(Enumerable.Range(0, rng.Next(4, 9)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var sb = new System.Text.StringBuilder();
        while (sb.Length < n)
        {
            foreach (char c in motif)
            {
                int r = rng.Next(20);
                if (r == 0) continue;                                  // deletion
                sb.Append(r == 1 ? "ACGT"[rng.Next(4)] : c);           // substitution / copy
                if (r == 2) sb.Append("ACGT"[rng.Next(4)]);            // insertion
            }
            if (rng.Next(2) == 0) sb.Append("ACGT"[rng.Next(4)]);
        }
        return sb.ToString(0, n);
    }

    private static int Code(char c) => "ACGT".IndexOf(c);

    private static List<(int, int, int, int, int)> BruteForce(string s, int k, int minLength, bool edit, bool pal)
    {
        int n = s.Length;
        var good = new int[n + 1, n + 1, n + 1, n + 1]; // [i, j, e1, e2] = distance + 1 (0 = not a match)
        var dp = new int[n + 1, n + 1];
        // Direct: u = S[i..], w = S[j..]. Palindromic: u = S[i..], w read as revcomp ending at e2.
        for (int i = 0; i < n; i++)
        {
            for (int other = pal ? 1 : i + 1; other <= (pal ? n : n - 1); other++)
            {
                int ml = n - i, mr = pal ? other : n - other;
                for (int a = 0; a <= ml; a++)
                {
                    for (int b = 0; b <= mr; b++)
                    {
                        int v;
                        if (a == 0) v = b;
                        else if (b == 0) v = a;
                        else
                        {
                            char x = s[i + a - 1];
                            int y = pal ? 3 - Code(s[other - b]) : Code(s[other + b - 1]);
                            bool match = Code(x) >= 0 && (pal ? Code(s[other - b]) >= 0 && Code(x) == y : Code(x) == y);
                            v = dp[a - 1, b - 1] + (match ? 0 : 1);
                            if (edit) v = Math.Min(v, Math.Min(dp[a - 1, b], dp[a, b - 1]) + 1);
                        }
                        if (!edit && a != b) v = int.MaxValue / 2;
                        dp[a, b] = v;
                        if (a == 0 || b == 0 || v > k) continue;
                        int j = pal ? other - b : other, e2 = pal ? other : other + b;
                        if (!pal && edit && !Accept(v, a, i, b, j)) continue;
                        good[i, j, i + a, e2] = v + 1;
                    }
                }
            }
        }

        // dom[i, j, e1, e2]: some match (i′ ≤ i, j′ ≤ j, e1′ ≥ e1, e2′ ≥ e2) exists — containment by dominance.
        var dom = new bool[n + 1, n + 1, n + 2, n + 2];
        for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
        for (int e1 = n; e1 >= 0; e1--)
        for (int e2 = n; e2 >= 0; e2--)
        {
            dom[i, j, e1, e2] = good[i, j, e1, e2] != 0
                || (i > 0 && dom[i - 1, j, e1, e2]) || (j > 0 && dom[i, j - 1, e1, e2])
                || dom[i, j, e1 + 1, e2] || dom[i, j, e1, e2 + 1];
        }

        var result = new List<(int, int, int, int, int)>();
        for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        for (int e1 = i + 1; e1 <= n; e1++)
        for (int e2 = j + 1; e2 <= n; e2++)
        {
            if (good[i, j, e1, e2] == 0 || (pal && j < i) || e1 - i < minLength || e2 - j < minLength) continue;
            bool contained = (i > 0 && dom[i - 1, j, e1, e2]) || (j > 0 && dom[i, j - 1, e1, e2])
                || dom[i, j, e1 + 1, e2] || dom[i, j, e1, e2 + 1];
            if (!contained) result.Add((i, j, e1 - i, e2 - j, good[i, j, e1, e2] - 1));
        }
        return result;
    }

    /// <summary>Vmatch <c>acceptmatch</c> (kurtz/extendED.c), positions p1 &lt; p2.</summary>
    private static bool Accept(int dist, int l1, int p1, int l2, int p2)
    {
        if (p1 >= p2) return false;
        if (p1 + l1 - 1 < p2) return true;
        if (p1 + l1 >= p2 + l2) return false;
        return (p2 - p1) + (p2 + l2) - (p1 + l1) > dist;
    }

    #endregion
}
