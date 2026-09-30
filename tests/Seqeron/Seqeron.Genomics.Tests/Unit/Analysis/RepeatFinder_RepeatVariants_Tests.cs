namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-DIRECT-001 repeat-enumeration variants (B04 completeness audit WP2):
/// <list type="bullet">
/// <item><see cref="RepeatFinder.FindReverseComplementRepeats(string,int,int,int)"/> — maximal reverse-complement
/// pairs = MUMmer <c>repeat-match</c> (without <c>-f</c>) <c>r</c> lines = Vmatch <c>-p</c>.</item>
/// <item><see cref="RepeatFinder.FindApproximateDirectRepeats(string,int,int,int,int,bool)"/> — maximal k-mismatch
/// (Hamming) repeats, REPuter / Vmatch <c>-h k</c> (<c>excludeContained</c> = Vmatch <c>-h k -allmax</c>).</item>
/// <item><see cref="RepeatFinder.FindSupermaximalRepeats(string,int)"/> — Gusfield 1997 §7.12.1 = Vmatch <c>-supermax</c>.</item>
/// </list>
/// Every expected list below was produced by the reference binaries (mummer4 <c>repeat-match</c> compiled from
/// source; Vmatch 2.3.1 Debian build: <c>mkvtree -dna -pl -allout</c> + <c>vmatch</c>) and agrees with an
/// independent brute force of the published definition.
/// </summary>
[TestFixture]
public class RepeatFinder_RepeatVariants_Tests
{
    private static string Triples(IEnumerable<ReverseComplementRepeatResult> r) =>
        string.Join(";", r.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.Length}"));

    private static string Quads(IEnumerable<ApproximateDirectRepeatResult> r) =>
        string.Join(";", r.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.Length},{x.Mismatches}"));

    private static string Supers(IEnumerable<SupermaximalRepeatResult> r) =>
        string.Join(";", r.Select(x => $"{x.Length}:{string.Join(",", x.Positions)}"));

    private static string RandomDna(int n, int seed, string alphabet = "ACGT")
    {
        var rng = new Random(seed);
        return new string(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
    }

    private static string RevComp(string s) =>
        new(s.Reverse().Select(SequenceExtensions.GetComplementBase).ToArray());

    #region Reverse-complement maximal pairs

    /// <summary>
    /// repeat-match (without -f, sequence + N sentinel) prints 8 19r 12 / 7 12r 6 / 16 19r 4; Vmatch -p -l 3
    /// prints (i, k, L) = (7,7,12) (15,15,4) (6,6,6). Start2 of repeat-match = k + L (last base of copy 2, 1-based).
    /// </summary>
    [Test]
    public void FindReverseComplementRepeats_MatchesRepeatMatchAndVmatch()
    {
        const string seq = "AAAAAAAACGTTGCAACGTAAAA";
        var results = RepeatFinder.FindReverseComplementRepeats(seq, 3, int.MaxValue, int.MinValue).ToList();

        Assert.That(Triples(results), Is.EqualTo("6,6,6;7,7,12;15,15,4"));
        var repeatMatchLines = results.Select(r => (r.FirstPosition + 1, r.SecondPosition + r.Length, r.Length))
            .OrderBy(t => t).ToList();
        Assert.That(repeatMatchLines, Is.EqualTo(new[] { (7, 12, 6), (8, 19, 12), (16, 19, 4) }));
        Assert.That(results[1].RepeatSequence, Is.EqualTo("ACGTTGCAACGT"));
        Assert.That(results[1].Spacing, Is.EqualTo(-12));
    }

    /// <summary>Hairpin-rich sequence: all 14 repeat-match reverse lines (min 4), Vmatch -p identical.</summary>
    [Test]
    public void FindReverseComplementRepeats_Hairpins_FullRepeatMatchList()
    {
        var results = RepeatFinder.FindReverseComplementRepeats(
            "TTGCATGCAAAAAATTTTTTTGCATGCAA", 4, int.MaxValue, int.MinValue);

        Assert.That(Triples(results), Is.EqualTo(
            "0,0,10;0,15,14;1,1,4;1,20,4;5,5,4;5,24,4;8,8,12;8,14,4;8,14,5;9,16,5;10,17,4;19,19,10;20,20,4;24,24,4"));
    }

    /// <summary>Defaults (min 5, max 50, spacing ≥ 1) keep only separated copies: (0,15,14) hairpin stem etc.</summary>
    [Test]
    public void FindReverseComplementRepeats_Defaults_SeparatedCopiesOnly()
    {
        var results = RepeatFinder.FindReverseComplementRepeats("TTGCATGCAAAAAATTTTTTTGCATGCAA").ToList();

        Assert.That(results.Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Spacing)),
            Is.EqualTo(new[] { (0, 15, 14, 1), (8, 14, 5, 1), (9, 16, 5, 2) }));
        Assert.That(results[0].SecondSequence, Is.EqualTo(RevComp(results[0].RepeatSequence)));
    }

    /// <summary>T₇ facing A₇: one (i, k) with several maximal lengths on different anti-diagonals (Vmatch -p).</summary>
    [Test]
    public void FindReverseComplementRepeats_SamePairSeveralLengths_OrderedByLength()
    {
        var results = RepeatFinder.FindReverseComplementRepeats(
            "GGATCCTTTTTTTGGATCCAAAAAAA", 4, int.MaxValue, int.MinValue);

        Assert.That(Triples(results), Is.EqualTo("0,0,6;0,13,6;6,6,20;6,19,4;6,19,5;6,19,6;7,20,6;8,21,5;9,22,4"));
    }

    [Test]
    public void FindReverseComplementRepeats_NonAcgtNeverPairs_CaseInsensitive()
    {
        Assert.That(Triples(RepeatFinder.FindReverseComplementRepeats("GAATTCNGANTTC", 3, int.MaxValue, int.MinValue)),
            Is.EqualTo("0,0,6;0,10,3"));
        Assert.That(RepeatFinder.FindReverseComplementRepeats(new string('N', 30), 2, 50, int.MinValue), Is.Empty);
        var lower = RepeatFinder.FindReverseComplementRepeats("gaattcAAAAAgaattc", 6).ToList();
        Assert.That(lower.Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Spacing, r.RepeatSequence)),
            Is.EqualTo(new[] { (0, 11, 6, 5, "GAATTC") }));
    }

    [Test]
    public void FindReverseComplementRepeats_RandomSequence_PairsAreMaximalReverseComplements()
    {
        string seq = RandomDna(2000, 11);
        var results = RepeatFinder.FindReverseComplementRepeats(seq, 6, int.MaxValue, int.MinValue).ToList();

        Assert.That(results, Is.Not.Empty);
        Assert.That(results.Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).Distinct().Count(), Is.EqualTo(results.Count));
        Assert.That(results, Is.Ordered.By(nameof(ReverseComplementRepeatResult.FirstPosition))
            .Then.By(nameof(ReverseComplementRepeatResult.SecondPosition)).Then.By(nameof(ReverseComplementRepeatResult.Length)));
        int n = seq.Length;
        foreach (var r in results)
        {
            int i = r.FirstPosition, k = r.SecondPosition, len = r.Length;
            Assert.That(i, Is.LessThanOrEqualTo(k));
            Assert.That(seq.Substring(i, len), Is.EqualTo(RevComp(seq.Substring(k, len))));
            Assert.That(r.Spacing, Is.EqualTo(k - i - len));
            bool outward = i > 0 && k + len < n && SequenceExtensions.GetComplementBase(seq[i - 1]) == seq[k + len];
            bool inward = i + len < n && k > 0 && SequenceExtensions.GetComplementBase(seq[i + len]) == seq[k - 1];
            Assert.That(outward || inward, Is.False, $"not maximal: ({i},{k},{len})");
        }
    }

    [Test]
    public void FindReverseComplementRepeats_DnaOverloadEqualsString_AndValidatesEagerly()
    {
        string seq = RandomDna(500, 3);
        Assert.That(RepeatFinder.FindReverseComplementRepeats(new DnaSequence(seq), 4, 30, 0),
            Is.EqualTo(RepeatFinder.FindReverseComplementRepeats(seq, 4, 30, 0)));
        Assert.That(RepeatFinder.FindReverseComplementRepeats((string)null!, 5), Is.Empty);
        Assert.That(RepeatFinder.FindReverseComplementRepeats("", 5), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindReverseComplementRepeats((DnaSequence)null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindReverseComplementRepeats("ACGT", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindReverseComplementRepeats("ACGT", 5, 4));
    }

    #endregion

    #region Degenerate (k-mismatch) direct repeats

    /// <summary>Vmatch -l 10 -h 2 (and -allmax): one repeat (0, 21, 17) with 2 mismatches (T/A at 3 and 11).</summary>
    [Test]
    public void FindApproximateDirectRepeats_TwoMismatches_MatchesVmatch()
    {
        var results = RepeatFinder.FindApproximateDirectRepeats(
            "ACGTTGCAAGCTTACGGGGGGACGATGCAAGCATACGG", 10, 2, int.MaxValue, int.MinValue).ToList();

        Assert.That(Quads(results), Is.EqualTo("0,21,17,2"));
        Assert.Multiple(() =>
        {
            Assert.That(results[0].FirstCopy, Is.EqualTo("ACGTTGCAAGCTTACGG"));
            Assert.That(results[0].SecondCopy, Is.EqualTo("ACGATGCAAGCATACGG"));
            Assert.That(results[0].Spacing, Is.EqualTo(4));
        });
    }

    /// <summary>A boundary-limited exact repeat is maximal with 0 ≤ k mismatches (Vmatch -h 1: 0,7,7 distance 0).</summary>
    [Test]
    public void FindApproximateDirectRepeats_ExactRepeatAtBoundaries_ZeroMismatches()
    {
        Assert.That(Quads(RepeatFinder.FindApproximateDirectRepeats("GATTACAGATTACA", 6, 1, int.MaxValue, int.MinValue)),
            Is.EqualTo("0,7,7,0"));
    }

    /// <summary>
    /// Per-diagonal maximality (default) vs the literal Vmatch containment rule: Vmatch -l 5 -h 1 -allmax prints
    /// exactly the excludeContained list; (0,2,7) (0,3,6) (0,4,5) lie inside (0,1,8) on another diagonal.
    /// </summary>
    [Test]
    public void FindApproximateDirectRepeats_ExcludeContained_EqualsVmatchAllmax()
    {
        const string seq = "AAAAAAAACGTTGCAACGTAAAA";
        Assert.That(Quads(RepeatFinder.FindApproximateDirectRepeats(seq, 5, 1, int.MaxValue, int.MinValue)),
            Is.EqualTo("0,1,8,1;0,2,7,1;0,3,6,1;0,4,5,1;0,18,5,1;1,18,5,1;2,18,5,1;3,18,5,1;5,13,6,1;6,14,6,1"));
        Assert.That(Quads(RepeatFinder.FindApproximateDirectRepeats(seq, 5, 1, int.MaxValue, int.MinValue, excludeContained: true)),
            Is.EqualTo("0,1,8,1;0,18,5,1;1,18,5,1;2,18,5,1;3,18,5,1;5,13,6,1;6,14,6,1"));
    }

    /// <summary>N is a mismatch (Vmatch wildcard rule): -l 8 -h 1 → (0,4,8,1) and (0,14,10,1) with A/N.</summary>
    [Test]
    public void FindApproximateDirectRepeats_NonAcgtCountsAsMismatch()
    {
        const string seq = "ACGTACGTACTTTTACGTNCGTAC";
        var results = RepeatFinder.FindApproximateDirectRepeats(seq, 8, 1, int.MaxValue, int.MinValue).ToList();

        Assert.That(Quads(results), Is.EqualTo("0,4,8,1;0,14,10,1"));
        Assert.That(results[1].SecondCopy, Is.EqualTo("ACGTNCGTAC"));
        Assert.That(RepeatFinder.FindApproximateDirectRepeats(new string('N', 40), 5, 2, int.MaxValue, int.MinValue), Is.Empty);
    }

    [Test]
    public void FindApproximateDirectRepeats_ZeroMismatches_EqualsFindDirectRepeats()
    {
        foreach (int seed in new[] { 1, 2, 3 })
        {
            string seq = RandomDna(1500, seed, seed == 3 ? "ACGTN" : "ACGT");
            foreach (var (min, max, sp) in new[] { (5, 50, 1), (4, int.MaxValue, int.MinValue), (6, 8, 0) })
            {
                var direct = RepeatFinder.FindDirectRepeats(seq, min, max, sp)
                    .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
                var approx = RepeatFinder.FindApproximateDirectRepeats(seq, min, 0, max, sp)
                    .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
                Assert.That(approx, Is.EqualTo(direct), $"seed {seed} ({min},{max},{sp})");
            }
        }
    }

    [Test]
    public void FindApproximateDirectRepeats_Random_EachRepeatMaximalWithinK()
    {
        string seq = RandomDna(1200, 21);
        const int k = 2;
        var results = RepeatFinder.FindApproximateDirectRepeats(seq, 12, k, int.MaxValue, int.MinValue).ToList();

        Assert.That(results, Is.Not.Empty);
        Assert.That(results, Is.Ordered.By(nameof(ApproximateDirectRepeatResult.FirstPosition))
            .Then.By(nameof(ApproximateDirectRepeatResult.SecondPosition)).Then.By(nameof(ApproximateDirectRepeatResult.Length)));
        int n = seq.Length;
        foreach (var r in results)
        {
            int d = r.SecondPosition - r.FirstPosition;
            int hamming = Enumerable.Range(0, r.Length).Count(t => r.FirstCopy[t] != r.SecondCopy[t]);
            Assert.That(hamming, Is.EqualTo(r.Mismatches));
            Assert.That(r.Mismatches, Is.LessThanOrEqualTo(k));
            Assert.That(d, Is.GreaterThan(0));
            bool leftExt = r.FirstPosition > 0 && r.Mismatches + (seq[r.FirstPosition - 1] != seq[r.SecondPosition - 1] ? 1 : 0) <= k;
            bool rightExt = r.SecondPosition + r.Length < n
                && r.Mismatches + (seq[r.FirstPosition + r.Length] != seq[r.SecondPosition + r.Length] ? 1 : 0) <= k;
            Assert.That(leftExt || rightExt, Is.False, $"not maximal: {r}");
        }
    }

    [Test]
    public void FindApproximateDirectRepeats_Validation()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 5, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 5, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 5, 1, 4));
        Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindApproximateDirectRepeats((DnaSequence)null!));
        Assert.That(RepeatFinder.FindApproximateDirectRepeats((string)null!), Is.Empty);
        string seq = RandomDna(400, 5);
        Assert.That(RepeatFinder.FindApproximateDirectRepeats(new DnaSequence(seq), 8, 1, 30, 0),
            Is.EqualTo(RepeatFinder.FindApproximateDirectRepeats(seq, 8, 1, 30, 0)));
    }

    #endregion

    #region Supermaximal repeats

    /// <summary>Values = Vmatch -supermax -l m (which prints every position pair of each repeat).</summary>
    [TestCase("AAAAAAAACGTTGCAACGTAAAA", 2, "7:0,1;5:6,14")]
    [TestCase("GATTACAGATTACA", 3, "7:0,7")]
    [TestCase("ACGTACGTTTTTTTTTACGTACGT", 3, "8:0,16;8:7,8")]
    [TestCase("CAGCAGCAGTTTCAGCAG", 3, "6:0,3,12")]
    public void FindSupermaximalRepeats_MatchesVmatchSupermax(string seq, int minLength, string expected)
    {
        Assert.That(Supers(RepeatFinder.FindSupermaximalRepeats(seq, minLength)), Is.EqualTo(expected));
    }

    [Test]
    public void FindSupermaximalRepeats_SequenceAndOccurrences()
    {
        var r = RepeatFinder.FindSupermaximalRepeats("CAGCAGCAGTTTCAGCAG", 3).Single();
        Assert.That(r.Sequence, Is.EqualTo("CAGCAG"));
        Assert.That(r.Length, Is.EqualTo(6));

        // N never matches and a suffix preceded by N has an undefined (distinct) left character.
        Assert.That(Supers(RepeatFinder.FindSupermaximalRepeats("acgtnacgtn", 3)), Is.EqualTo("4:0,5"));
        Assert.That(RepeatFinder.FindSupermaximalRepeats(new string('N', 20), 1), Is.Empty);
    }

    [Test]
    public void FindSupermaximalRepeats_Random_AreMaximalRepeatsNotContainedInOthers()
    {
        string seq = RandomDna(3000, 17);
        var supers = RepeatFinder.FindSupermaximalRepeats(seq, 1).ToList();
        var maximalRepeats = RepeatFinder.FindDirectRepeats(seq, 2, int.MaxValue, int.MinValue)
            .Select(x => x.RepeatSequence).ToHashSet();

        Assert.That(supers, Is.Not.Empty);
        foreach (var s in supers.Where(s => s.Length >= 2))
        {
            Assert.That(maximalRepeats, Does.Contain(s.Sequence));
            Assert.That(maximalRepeats.Any(o => o.Length > s.Length && o.Contains(s.Sequence, StringComparison.Ordinal)), Is.False);
            Assert.That(s.Positions, Is.Ordered);
            Assert.That(s.Positions.All(p => seq.Substring(p, s.Length) == s.Sequence));
            int occurrences = Enumerable.Range(0, seq.Length - s.Length + 1)
                .Count(p => string.CompareOrdinal(seq, p, s.Sequence, 0, s.Length) == 0);
            Assert.That(s.Positions, Has.Count.EqualTo(occurrences));
        }
    }

    [Test]
    public void FindSupermaximalRepeats_Validation()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindSupermaximalRepeats("ACGT", 0));
        Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindSupermaximalRepeats((DnaSequence)null!));
        Assert.That(RepeatFinder.FindSupermaximalRepeats((string)null!), Is.Empty);
        Assert.That(RepeatFinder.FindSupermaximalRepeats("A"), Is.Empty);
        string seq = RandomDna(300, 9);
        Assert.That(Supers(RepeatFinder.FindSupermaximalRepeats(new DnaSequence(seq), 3)),
            Is.EqualTo(Supers(RepeatFinder.FindSupermaximalRepeats(seq, 3))));
    }

    #endregion
}
