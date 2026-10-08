// PAT-APPROX-002 — weighted Damerau distances (Lowrance & Wagner 1975) and transposition-aware
// tracebacks (GetDamerauLevenshteinAlignment / GetOptimalStringAlignment).
// TestSpec: tests/TestSpecs/PAT-APPROX-002.md
// Review 2026-09 B05 audit round 2 group G3 (docs/Validation/review-2026-09/B05.md, F34).
// Sources: Lowrance & Wagner (1975) J. ACM 22(2):177–183 (costs W_I, W_D, W_C, W_S; exact when
//   2·W_S ≥ W_I + W_D) — paper not reachable (ACM DL / Springer / BMC proxy-blocked), recurrence and
//   condition taken from the WebSearch records of Zhao & Sahni (2019) BMC Bioinformatics 20:277;
//   R stringdist 0.9.12 src/osa.c, src/dl.c (raw.githubusercontent.com/markvanderloo/stringdist).
// References (computed 2026-10-01): R stringdist 0.9.12 stringdist(a, b, method = "osa" | "dl",
//   weight = c(ins, del, sub, trans) / max) * max — note stringdist's "d" weight is charged for
//   characters of b (our insertion); osa agrees on 5000/5000 random cases, dl on 1723/1723 cases
//   with del = ins = trans (it charges in-between characters at the transposition weight otherwise);
//   exact minimum over all edit sequences = Python Dijkstra over the edit graph (ins/del/sub/adjacent
//   swap), which equals the Lowrance–Wagner recurrence on 3225/3225 random cases with
//   2·W_S ≥ W_I + W_D and differs on 105 cases violating it.
namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class ApproximateMatcher_WeightedDamerau_Tests
{
    // (s1, s2, ins, del, sub, trans, OSA, DL): OSA = stringdist "osa"; DL = Dijkstra exact minimum
    // (also = stringdist "dl" except where noted in the trailing comment).
    [TestCase("CA", "ABC", 1, 1, 1, 1, 3, 2)] // stringdist osa=3 dl=2 dijkstra=2
    [TestCase("CA", "ABC", 2, 3, 4, 3, 7, 5)] // stringdist osa=7 dl=6 dijkstra=5
    [TestCase("ABC", "CA", 2, 3, 4, 3, 8, 6)] // stringdist osa=8 dl=6 dijkstra=6
    [TestCase("ab", "ba", 1, 1, 1, 1, 1, 1)] // stringdist osa=1 dl=1 dijkstra=1
    [TestCase("ab", "ba", 1, 1, 5, 1, 1, 1)] // stringdist osa=1 dl=1 dijkstra=1
    [TestCase("ab", "ba", 1, 1, 5, 3, 2, 2)] // stringdist osa=2 dl=2 dijkstra=2
    [TestCase("ABCD", "DABC", 1, 1, 1, 1, 2, 2)] // stringdist osa=2 dl=2 dijkstra=2
    [TestCase("ABCD", "DABC", 2, 2, 3, 2, 4, 4)] // stringdist osa=4 dl=4 dijkstra=4
    [TestCase("ACGT", "TGCA", 1, 2, 3, 2, 8, 8)] // stringdist osa=8 dl=8 dijkstra=8
    [TestCase("kitten", "sitting", 1, 1, 2, 1, 5, 5)] // stringdist osa=5 dl=5 dijkstra=5
    [TestCase("a cat", "an act", 2, 2, 3, 2, 4, 4)] // stringdist osa=4 dl=4 dijkstra=4
    [TestCase("CAAAACA", "CA", 1, 7, 5, 5, 35, 35)] // stringdist osa=35 dl=32 dijkstra=35
    [TestCase("AGGG", "GAGCACA", 1, 3, 5, 2, 9, 9)] // stringdist osa=9 dl=8 dijkstra=9
    [TestCase("ABC", "BCA", 2, 2, 4, 2, 4, 4)] // stringdist osa=4 dl=4 dijkstra=4
    [TestCase("ACGTACGT", "CAGTCAGT", 1, 1, 3, 1, 2, 2)] // stringdist osa=2 dl=2 dijkstra=2
    [TestCase("ACGTACGT", "CAGTCAGT", 3, 2, 4, 3, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("", "ACG", 2, 3, 1, 3, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("ACG", "", 2, 3, 1, 3, 9, 9)] // stringdist osa=9 dl=9 dijkstra=9
    [TestCase("ACGT", "ACGT", 5, 5, 5, 5, 0, 0)] // stringdist osa=0 dl=0 dijkstra=0
    [TestCase("AXB", "BYA", 1, 1, 1, 1, 3, 3)] // stringdist osa=3 dl=3 dijkstra=3
    [TestCase("AXB", "BYA", 2, 1, 3, 2, 6, 5)] // stringdist osa=6 dl=6 dijkstra=5
    [TestCase("ADCB", "ABDC", 1, 1, 1, 1, 2, 2)] // stringdist osa=2 dl=2 dijkstra=2
    [TestCase("ADCB", "ABDC", 4, 2, 5, 3, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("BBCCAA", "BCCAAB", 3, 2, 1, 5, 3, 3)] // stringdist osa=3 dl=3 dijkstra=3
    [TestCase("TTCAAA", "CC", 5, 1, 4, 3, 8, 8)] // stringdist osa=8 dl=8 dijkstra=8
    [TestCase("TT", "CTC", 6, 2, 4, 6, 10, 10)] // stringdist osa=10 dl=10 dijkstra=10
    [TestCase("BB", "ABC", 3, 1, 2, 2, 5, 5)] // stringdist osa=5 dl=5 dijkstra=5
    [TestCase("CAGGTC", "TCTCCG", 2, 5, 2, 6, 12, 12)] // stringdist osa=12 dl=12 dijkstra=12
    [TestCase("CTT", "ATAAAA", 5, 3, 2, 6, 19, 19)] // stringdist osa=19 dl=19 dijkstra=19
    [TestCase("CCAA", "TTCA", 1, 2, 3, 3, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("AG", "GTTTT", 5, 1, 5, 5, 21, 21)] // stringdist osa=21 dl=21 dijkstra=21
    [TestCase("ATAC", "TT", 3, 1, 3, 3, 5, 5)] // stringdist osa=5 dl=5 dijkstra=5
    [TestCase("CCGTTC", "TGCT", 2, 2, 4, 2, 8, 8)] // stringdist osa=8 dl=8 dijkstra=8
    [TestCase("AAB", "ACAB", 2, 1, 4, 4, 2, 2)] // stringdist osa=2 dl=2 dijkstra=2
    [TestCase("BBC", "BCBCCB", 3, 1, 1, 3, 9, 9)] // stringdist osa=9 dl=9 dijkstra=9
    [TestCase("CC", "ACAA", 4, 3, 2, 6, 10, 10)] // stringdist osa=10 dl=10 dijkstra=10
    [TestCase("GTGCA", "CC", 4, 1, 3, 5, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("AA", "AACBC", 2, 4, 2, 5, 6, 6)] // stringdist osa=6 dl=6 dijkstra=6
    [TestCase("CABAC", "BCBABC", 1, 5, 1, 5, 3, 3)] // stringdist osa=3 dl=3 dijkstra=3
    [TestCase("CA", "BCCABA", 6, 2, 1, 5, 24, 24)] // stringdist osa=24 dl=24 dijkstra=24 
    public void WeightedDistances_AndAlignments_EqualReferences(
        string s1, string s2, int ins, int del, int sub, int trans, int osa, int dl)
    {
        var costs = new DamerauCosts(ins, del, sub, trans);
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(s1, s2, costs), Is.EqualTo(osa));
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(s1, s2, new EditCosts(ins, del, sub), trans), Is.EqualTo(osa));
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(s1, s2, costs), Is.EqualTo(dl));
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(s1, s2, new EditCosts(ins, del, sub), trans), Is.EqualTo(dl));

            var o = ApproximateMatcher.GetOptimalStringAlignment(s1, s2, costs);
            Assert.That(o.Distance, Is.EqualTo(osa));
            AssertReplay(s1, s2, costs, o, restricted: true);

            var d = ApproximateMatcher.GetDamerauLevenshteinAlignment(s1, s2, costs);
            Assert.That(d.Distance, Is.EqualTo(dl));
            AssertReplay(s1, s2, costs, d, restricted: false);
        });
    }

    [Test]
    public void WeightedDamerauLevenshtein_RandomSmall_EqualsExhaustiveDijkstra()
    {
        var rng = new Random(20261001);
        int checkedCases = 0;
        while (checkedCases < 600)
        {
            string a = RandomString(rng, "ABC", rng.Next(0, 5));
            string b = RandomString(rng, "ABC", rng.Next(0, 5));
            var c = new DamerauCosts(rng.Next(0, 7), rng.Next(0, 7), rng.Next(0, 7), rng.Next(0, 7));
            if (2 * c.Transposition < c.Insertion + c.Deletion)
                continue;
            checkedCases++;

            long expected = Dijkstra(a, b, c);
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, c), Is.EqualTo(expected), $"{a} → {b} {c}");
            var al = ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b, c);
            Assert.That(al.Distance, Is.EqualTo(expected), $"{a} → {b} {c}");
            AssertReplay(a, b, c, al, restricted: false);

            var osa = ApproximateMatcher.GetOptimalStringAlignment(a, b, c);
            Assert.That(osa.Distance, Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, c)));
            AssertReplay(a, b, c, osa, restricted: true);
        }
    }

    [Test]
    public void UnitAndUniformCosts_AreBitIdenticalToUnitEngines_Exhaustive()
    {
        var strings = AllStrings("ABC", 4);
        foreach (string a in strings)
        {
            foreach (string b in strings)
            {
                int dl = ApproximateMatcher.DamerauLevenshteinDistance(a, b);
                int osa = ApproximateMatcher.OptimalStringAlignmentDistance(a, b);
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, DamerauCosts.Unit), Is.EqualTo(dl));
                Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, DamerauCosts.Unit), Is.EqualTo(osa));
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, new DamerauCosts(3, 3, 3, 3)), Is.EqualTo(3 * dl));
                Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, new DamerauCosts(3, 3, 3, 3)), Is.EqualTo(3 * osa));

                // The weighted full-matrix kernels with unit costs (reached through the alignments).
                var da = ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b);
                var oa = ApproximateMatcher.GetOptimalStringAlignment(a, b);
                Assert.That(da.Distance, Is.EqualTo(dl), $"{a} → {b}");
                Assert.That(oa.Distance, Is.EqualTo(osa), $"{a} → {b}");
                AssertReplay(a, b, DamerauCosts.Unit, da, restricted: false);
                AssertReplay(a, b, DamerauCosts.Unit, oa, restricted: true);
            }
        }
    }

    [Test]
    public void Alignment_CaToAbc_LowranceWagnerBlockVersusOsa()
    {
        // DL(CA, ABC) = 2: swap CA → AC with B inserted between (one block); OSA = 3.
        var dl = ApproximateMatcher.GetDamerauLevenshteinAlignment("CA", "ABC");
        Assert.Multiple(() =>
        {
            Assert.That(dl.Distance, Is.EqualTo(2));
            Assert.That(dl.Script, Is.EqualTo("Td"));
            Assert.That(dl.TranspositionCount, Is.EqualTo(1));
            Assert.That(dl.Operations, Is.EqualTo(new[]
            {
                new DamerauEditOperation(DamerauOperationKind.Transposition, 0, 2, 0, 3, 2),
            }));
        });

        var ab = ApproximateMatcher.GetOptimalStringAlignment("ab", "ba");
        Assert.That(ab.Script, Is.EqualTo("T"));
        Assert.That(ab.Operations[0], Is.EqualTo(new DamerauEditOperation(DamerauOperationKind.Transposition, 0, 2, 0, 2, 1)));

        var osa = ApproximateMatcher.GetOptimalStringAlignment("CA", "ABC");
        Assert.That(osa.Distance, Is.EqualTo(3));
        Assert.That(osa.TranspositionCount, Is.EqualTo(0));
    }

    [Test]
    public void Alignment_TranspositionWithDeletionBetween()
    {
        // Block a_k … a_i = "AXB" → b_l … b_j = "BA": delete X between the pair, then swap (W_D + W_S).
        var al = ApproximateMatcher.GetDamerauLevenshteinAlignment("AXB", "BA");
        Assert.Multiple(() =>
        {
            Assert.That(al.Distance, Is.EqualTo(2));
            Assert.That(al.Script, Is.EqualTo("Ti"));
            Assert.That(al.Operations.Single(), Is.EqualTo(new DamerauEditOperation(DamerauOperationKind.Transposition, 0, 3, 0, 2, 2)));
        });
    }

    [Test]
    public void LowranceWagnerConditionViolated_Throws_OsaStillDefined()
    {
        // W_I = W_D = 2, W_S = 1: ABC → BCA costs 2 by two swaps; the recurrence would give 4.
        var c = new DamerauCosts(2, 2, 4, 1);
        Assert.Throws<ArgumentException>(() => ApproximateMatcher.DamerauLevenshteinDistance("ABC", "BCA", c));
        Assert.Throws<ArgumentException>(() => ApproximateMatcher.GetDamerauLevenshteinAlignment("ABC", "BCA", c));
        Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance("ABC", "BCA", c), Is.EqualTo(4)); // stringdist osa
        Assert.That(ApproximateMatcher.GetOptimalStringAlignment("ABC", "BCA", c).Distance, Is.EqualTo(4));
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.DamerauLevenshteinDistance(null!, "A", DamerauCosts.Unit));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.OptimalStringAlignmentDistance("A", null!, DamerauCosts.Unit));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetDamerauLevenshteinAlignment(null!, "A"));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetOptimalStringAlignment("A", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.DamerauLevenshteinDistance("A", "B", new DamerauCosts(1, 1, 1, -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.OptimalStringAlignmentDistance("A", "B", new DamerauCosts(1, -1, 1, 1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.GetOptimalStringAlignment("A", "B", new DamerauCosts(1, 1, -1, 1)));
    }

    [Test]
    public void LargeCosts_OverflowChecked()
    {
        var c = new DamerauCosts(int.MaxValue, int.MaxValue, int.MaxValue, int.MaxValue);
        Assert.Throws<OverflowException>(() => ApproximateMatcher.DamerauLevenshteinDistance("AA", "", c));
        var w = new DamerauCosts(int.MaxValue / 2, int.MaxValue, 1, int.MaxValue);
        Assert.Throws<OverflowException>(() => ApproximateMatcher.OptimalStringAlignmentDistance("AA", "", w));
        Assert.That(ApproximateMatcher.DamerauLevenshteinDistance("A", "", w), Is.EqualTo(int.MaxValue));
    }

    /// <summary>
    /// Replays the operations on s1: positions are contiguous, matches / transposition ends agree
    /// with both strings, the produced string is s2 and the summed cost equals Distance.
    /// </summary>
    private static void AssertReplay(string s1, string s2, DamerauCosts c, DamerauAlignment al, bool restricted)
    {
        var sb = new System.Text.StringBuilder();
        int p = 0, q = 0;
        long cost = 0;
        foreach (var op in al.Operations)
        {
            Assert.That(op.SourcePosition, Is.EqualTo(p), "source contiguity");
            Assert.That(op.TargetPosition, Is.EqualTo(q), "target contiguity");
            long expectedCost;
            switch (op.Kind)
            {
                case DamerauOperationKind.Match:
                    Assert.That(s1[p], Is.EqualTo(s2[q]));
                    sb.Append(s1[p]);
                    expectedCost = 0;
                    break;
                case DamerauOperationKind.Substitution:
                    Assert.That(s1[p], Is.Not.EqualTo(s2[q]));
                    sb.Append(s2[q]);
                    expectedCost = c.Substitution;
                    break;
                case DamerauOperationKind.Insertion:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((0, 1)));
                    sb.Append(s2[q]);
                    expectedCost = c.Insertion;
                    break;
                case DamerauOperationKind.Deletion:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((1, 0)));
                    expectedCost = c.Deletion;
                    break;
                default:
                    int sl = op.SourceLength, tl = op.TargetLength;
                    Assert.That(sl, Is.GreaterThanOrEqualTo(2));
                    Assert.That(tl, Is.GreaterThanOrEqualTo(2));
                    if (restricted)
                        Assert.That((sl, tl), Is.EqualTo((2, 2)));
                    char first = s1[p], last = s1[p + sl - 1];
                    Assert.That(s2[q], Is.EqualTo(last), "a_i = b_l");
                    Assert.That(s2[q + tl - 1], Is.EqualTo(first), "a_k = b_j");
                    // delete the in-between source characters, swap, insert b_{l+1..j−1} between.
                    sb.Append(last).Append(s2, q + 1, tl - 2).Append(first);
                    expectedCost = (long)(sl - 2) * c.Deletion + c.Transposition + (long)(tl - 2) * c.Insertion;
                    break;
            }
            Assert.That(op.Cost, Is.EqualTo(expectedCost));
            cost += op.Cost;
            p += op.SourceLength;
            q += op.TargetLength;
        }
        Assert.That(p, Is.EqualTo(s1.Length));
        Assert.That(sb.ToString(), Is.EqualTo(s2), $"replay of {al.Script}");
        Assert.That(cost, Is.EqualTo(al.Distance), "summed cost = distance");
    }

    /// <summary>Exact minimum cost over all edit sequences (ins / del / sub / adjacent swap).</summary>
    private static long Dijkstra(string a, string b, DamerauCosts c)
    {
        var alphabet = (a + b).Distinct().ToArray();
        int maxLen = Math.Max(a.Length, b.Length) + 1;
        var dist = new Dictionary<string, long> { [a] = 0 };
        var pq = new PriorityQueue<string, long>();
        pq.Enqueue(a, 0);
        while (pq.TryDequeue(out var s, out long d))
        {
            if (s == b)
                return d;
            if (d > dist[s])
                continue;
            void Relax(string t, long w)
            {
                long nd = d + w;
                if (!dist.TryGetValue(t, out long old) || nd < old)
                {
                    dist[t] = nd;
                    pq.Enqueue(t, nd);
                }
            }
            for (int i = 0; i < s.Length; i++)
            {
                Relax(s.Remove(i, 1), c.Deletion);
                foreach (char ch in alphabet)
                    if (ch != s[i])
                        Relax(s.Substring(0, i) + ch + s.Substring(i + 1), c.Substitution);
                if (i + 1 < s.Length && s[i] != s[i + 1])
                    Relax(s.Substring(0, i) + s[i + 1] + s[i] + s.Substring(i + 2), c.Transposition);
            }
            if (s.Length < maxLen)
                for (int i = 0; i <= s.Length; i++)
                    foreach (char ch in alphabet)
                        Relax(s.Insert(i, ch.ToString()), c.Insertion);
        }
        throw new InvalidOperationException("unreachable");
    }

    private static string RandomString(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static List<string> AllStrings(string alphabet, int maxLength)
    {
        var result = new List<string> { string.Empty };
        var frontier = new List<string> { string.Empty };
        for (int len = 1; len <= maxLength; len++)
        {
            frontier = frontier.SelectMany(s => alphabet.Select(ch => s + ch)).ToList();
            result.AddRange(frontier);
        }
        return result;
    }
}
