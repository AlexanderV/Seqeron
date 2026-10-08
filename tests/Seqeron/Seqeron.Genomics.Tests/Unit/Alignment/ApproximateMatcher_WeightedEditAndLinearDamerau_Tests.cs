// PAT-APPROX-002 — Linear-space unrestricted Damerau–Levenshtein (Zhao & Sahni) and weighted
// (non-unit) Levenshtein distance / alignment / Sellers search.
// TestSpec: tests/TestSpecs/PAT-APPROX-002.md
// Review 2026-09 B05 audit group B (docs/Validation/review-2026-09/B05.md, F27).
// Sources: Zhao & Sahni (2020) BMC Bioinformatics 21(Suppl 1), doi:10.1186/s12859-019-3184-8, as
//   transcribed in rapidfuzz-cpp rapidfuzz/distance/DamerauLevenshtein_impl.hpp
//   (damerau_levenshtein_distance_zhao); rapidfuzz-cpp Levenshtein_impl.hpp
//   (generalized_levenshtein_wagner_fischer, LevenshteinWeightTable {insert, delete, replace}).
// References (values computed 2026-10-01): rapidfuzz 3.14.6
//   DamerauLevenshtein.distance / Levenshtein.distance(s1, s2, weights=(ins, del, sub)),
//   jellyfish damerau_levenshtein_distance; weighted Sellers sets = Python brute force
//   min_i Levenshtein.distance(p, t[i..j], weights).
namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class ApproximateMatcher_WeightedEditAndLinearDamerau_Tests
{
    #region DamerauLevenshteinDistance — Zhao & Sahni linear-space engine

    [TestCase("CA", "ABC", 2)]
    [TestCase("ABCDEF", "BADCFE", 3)]
    [TestCase("kitten", "sitting", 3)]
    [TestCase("ab", "ba", 1)]
    [TestCase("abcdef", "fedcba", 5)]
    [TestCase("ACGTACGT", "CAGTACTG", 2)]
    [TestCase("αβγ", "βαγ", 1)]
    [TestCase("", "ACG", 3)]
    [TestCase("GATTACA", "", 7)]
    [TestCase("", "", 0)]
    public void DamerauLevenshtein_RapidfuzzJellyfishValues(string s1, string s2, int expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(s1, s2), Is.EqualTo(expected));
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistanceFullMatrix(s1, s2), Is.EqualTo(expected));
        });
    }

    [Test]
    public void DamerauLevenshtein_ExhaustiveSmall_EqualsLowranceWagnerFullMatrix()
    {
        // Every pair of strings of length ≤ 5 over {A, B, C}: 364 × 364 pairs.
        var words = new List<string> { "" };
        for (int len = 1; len <= 5; len++)
            words.AddRange(Enumerate("ABC", len));

        int mismatches = 0;
        foreach (var a in words)
            foreach (var b in words)
                if (ApproximateMatcher.DamerauLevenshteinDistance(a, b) != ApproximateMatcher.DamerauLevenshteinDistanceFullMatrix(a, b))
                    mismatches++;

        Assert.That(mismatches, Is.Zero);
    }

    [Test]
    public void DamerauLevenshtein_Random_EqualsFullMatrix_IncludingNonAscii()
    {
        var rng = new Random(20261001);
        string[] alphabets = { "AC", "ACGT", "abcdefghij", "αβγδ€жA", "ACGTN-ü中" };
        for (int iter = 0; iter < 1500; iter++)
        {
            string al = alphabets[rng.Next(alphabets.Length)];
            string a = RandomString(rng, al, rng.Next(0, 60));
            string b = rng.Next(2) == 0 ? Mutate(rng, a, al) : RandomString(rng, al, rng.Next(0, 60));
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b),
                Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistanceFullMatrix(a, b)), $"{a} / {b}");
        }
    }

    [Test]
    public void DamerauLevenshtein_Long_EqualsFullMatrix()
    {
        var rng = new Random(7);
        string a = RandomString(rng, "ACGT", 1500);
        string b = Mutate(rng, a, "ACGT");
        Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b),
            Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistanceFullMatrix(a, b)));
    }

    #endregion

    #region EditDistance — weighted costs (rapidfuzz Levenshtein weights=(ins, del, sub))

    [TestCase("kitten", "sitting", 1, 1, 2, 5)]
    [TestCase("kitten", "sitting", 2, 1, 1, 4)]
    [TestCase("kitten", "sitting", 1, 3, 5, 9)]
    [TestCase("ACGT", "AGT", 1, 5, 1, 5)]
    [TestCase("AGT", "ACGT", 1, 5, 1, 1)]
    [TestCase("ACGT", "TGCA", 0, 1, 1, 3)]
    [TestCase("ACGT", "TGCA", 1, 0, 1, 3)]
    [TestCase("ACGT", "TGCA", 3, 3, 0, 0)]
    [TestCase("GATTACA", "GCATGCU", 2, 3, 4, 13)]
    [TestCase("", "ACG", 4, 1, 1, 12)]
    [TestCase("ACG", "", 4, 1, 1, 3)]
    [TestCase("ABC", "XYZ", 1, 1, 5, 6)]
    [TestCase("abc", "abc", 7, 7, 7, 0)]
    [TestCase("kitten", "sitting", 3, 3, 3, 9)]
    [TestCase("ACGTACGT", "ACGTTACG", 1000, 2, 999, 1002)]
    public void EditDistance_Weighted_RapidfuzzValues(string s1, string s2, int ins, int del, int sub, int expected)
    {
        var costs = new EditCosts(ins, del, sub);
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.EditDistance(s1, s2, ins, del, sub), Is.EqualTo(expected));
            Assert.That(ApproximateMatcher.EditDistance(s1, s2, costs), Is.EqualTo(expected));

            var a = ApproximateMatcher.GetEditAlignment(s1, s2, costs);
            Assert.That(a.Distance, Is.EqualTo(expected));
            Assert.That(WeightedReplayCost(a.Operations, s1, s2, costs), Is.EqualTo(expected));

            var l = ApproximateMatcher.GetEditAlignmentLinearSpace(s1, s2, costs);
            Assert.That(l.Distance, Is.EqualTo(expected));
            Assert.That(WeightedReplayCost(l.Operations, s1, s2, costs), Is.EqualTo(expected));
        });
    }

    [Test]
    public void EditDistance_Weighted_InsertionDeletionDirection_FollowsS1ToS2()
    {
        // rapidfuzz: deleting from s1 costs deletion, inserting characters of s2 costs insertion.
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.EditDistance("ACGT", "AGT", 1, 5, 9), Is.EqualTo(5));
            Assert.That(ApproximateMatcher.EditDistance("AGT", "ACGT", 1, 5, 9), Is.EqualTo(1));
            Assert.That(ApproximateMatcher.GetEditAlignment("ACGT", "AGT", new EditCosts(1, 5, 9)).Cigar, Is.EqualTo("1=1I2="));
            Assert.That(ApproximateMatcher.GetEditAlignment("AGT", "ACGT", new EditCosts(1, 5, 9)).Cigar, Is.EqualTo("1=1D2="));
        });
    }

    [Test]
    public void EditDistance_Weighted_SubstitutionAboveIndelPair_UsesIndels()
    {
        // sub (5) > ins + del (2): replacing is never chosen.
        var a = ApproximateMatcher.GetEditAlignment("ABC", "XYZ", new EditCosts(1, 1, 5));
        Assert.Multiple(() =>
        {
            Assert.That(a.Distance, Is.EqualTo(6));
            Assert.That(a.Operations, Does.Not.Contain("X"));
            Assert.That(ApproximateMatcher.GetEditAlignmentLinearSpace("ABC", "XYZ", new EditCosts(1, 1, 5)).Operations, Does.Not.Contain("X"));
        });
    }

    [Test]
    public void EditDistance_UnitCosts_EqualsUnitMethods_PathsIdentical()
    {
        var rng = new Random(424242);
        string[] alphabets = { "AC", "ACGT", "αβγδ" };
        for (int iter = 0; iter < 600; iter++)
        {
            string al = alphabets[rng.Next(alphabets.Length)];
            string a = RandomString(rng, al, rng.Next(0, 40));
            string b = rng.Next(2) == 0 ? Mutate(rng, a, al) : RandomString(rng, al, rng.Next(0, 40));
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance(a, b, EditCosts.Unit), Is.EqualTo(ApproximateMatcher.EditDistance(a, b)));
                Assert.That(ApproximateMatcher.EditDistance(a, b, 2, 2, 2), Is.EqualTo(2 * ApproximateMatcher.EditDistance(a, b)));
                Assert.That(ApproximateMatcher.GetEditAlignment(a, b, EditCosts.Unit), Is.EqualTo(ApproximateMatcher.GetEditAlignment(a, b)));
                Assert.That(ApproximateMatcher.GetEditAlignmentLinearSpace(a, b, EditCosts.Unit),
                    Is.EqualTo(ApproximateMatcher.GetEditAlignmentLinearSpace(a, b)));
            });
        }
    }

    [Test]
    public void EditDistance_Weighted_Random_AlignmentsOptimalAndValid()
    {
        var rng = new Random(99);
        for (int iter = 0; iter < 800; iter++)
        {
            string a = RandomString(rng, "ACGT", rng.Next(0, 30));
            string b = rng.Next(2) == 0 ? Mutate(rng, a, "ACGT") : RandomString(rng, "ACGT", rng.Next(0, 30));
            var costs = new EditCosts(rng.Next(0, 7), rng.Next(0, 7), rng.Next(0, 7));
            int d = ApproximateMatcher.EditDistance(a, b, costs);
            var full = ApproximateMatcher.GetEditAlignment(a, b, costs);
            var lin = ApproximateMatcher.GetEditAlignmentLinearSpace(a, b, costs);
            Assert.Multiple(() =>
            {
                Assert.That(d, Is.EqualTo(ReferenceWeighted(a, b, costs)));
                Assert.That(full.Distance, Is.EqualTo(d));
                Assert.That(lin.Distance, Is.EqualTo(d));
                Assert.That(WeightedReplayCost(full.Operations, a, b, costs), Is.EqualTo(d));
                Assert.That(WeightedReplayCost(lin.Operations, a, b, costs), Is.EqualTo(d));
                // Swapping the strings swaps the roles of insertion and deletion.
                Assert.That(ApproximateMatcher.EditDistance(b, a, new EditCosts(costs.Deletion, costs.Insertion, costs.Substitution)), Is.EqualTo(d));
            });
        }
    }

    [Test]
    public void EditDistance_Weighted_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.EditDistance("A", "C", -1, 1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.EditDistance("A", "C", 1, -1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.EditDistance("A", "C", 1, 1, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.GetEditAlignment("A", "C", new EditCosts(1, 1, -2)));
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.GetEditAlignmentLinearSpace("A", "C", new EditCosts(-1, 1, 1)));
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.FindEditEndPositions("ACGT", "AC", 1, new EditCosts(1, -1, 1)).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() => ApproximateMatcher.FindEditEndPositions("ACGT", "AC", -1, EditCosts.Unit).ToList());
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.EditDistance(null!, "C", EditCosts.Unit));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetEditAlignment("A", null!, EditCosts.Unit));
            Assert.Throws<OverflowException>(() => ApproximateMatcher.EditDistance(new string('A', 3), "", 1, int.MaxValue, 1));
            Assert.That(ApproximateMatcher.EditDistance("ACGT", "TTTTTT", 0, 0, 0), Is.Zero);
        });
    }

    #endregion

    #region FindEditEndPositions — weighted Sellers

    [Test]
    public void FindEditEndPositions_Weighted_BruteForceValues()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.FindEditEndPositions("TTACGTAAGGCTACG", "ACGT", 2, EditCosts.Unit).ToList(),
                Is.EqualTo(new[] { (3, 2), (4, 1), (5, 0), (6, 1), (7, 2), (8, 2), (9, 2), (10, 2), (11, 2), (13, 2), (14, 1) }));
            Assert.That(ApproximateMatcher.FindEditEndPositions("TTACGTAAGGCTACG", "ACGT", 3, new EditCosts(2, 1, 3)).ToList(),
                Is.EqualTo(new[] { (0, 3), (1, 3), (2, 3), (3, 2), (4, 1), (5, 0), (6, 2), (7, 3), (8, 2), (9, 3), (10, 3), (11, 2), (12, 3), (13, 2), (14, 1) }));
            Assert.That(ApproximateMatcher.FindEditEndPositions("ttacgtaaggctacg", "acgt", 2, new EditCosts(1, 3, 1)).ToList(),
                Is.EqualTo(new[] { (5, 0), (6, 1), (7, 2), (9, 2), (10, 2), (11, 2) }));
        });
    }

    [Test]
    public void FindEditEndPositions_UnitCosts_EqualsMyersSellers()
    {
        var rng = new Random(31337);
        for (int iter = 0; iter < 300; iter++)
        {
            string t = RandomString(rng, "ACGT", rng.Next(1, 60));
            string p = RandomString(rng, "ACGT", rng.Next(1, 10));
            int k = rng.Next(0, 6);
            Assert.That(ApproximateMatcher.FindEditEndPositions(t, p, k, EditCosts.Unit).ToList(),
                Is.EqualTo(ApproximateMatcher.FindEditEndPositions(t, p, k).ToList()));
        }
    }

    [Test]
    public void FindEditEndPositions_Weighted_Random_EqualsBruteForce()
    {
        var rng = new Random(5);
        for (int iter = 0; iter < 200; iter++)
        {
            string t = RandomString(rng, "ACGT", rng.Next(1, 25));
            string p = RandomString(rng, "ACGT", rng.Next(1, 7));
            var costs = new EditCosts(rng.Next(0, 5), rng.Next(0, 5), rng.Next(0, 5));
            int k = rng.Next(0, 10);
            var expected = new List<(int, int)>();
            for (int j = 0; j < t.Length; j++)
            {
                int best = int.MaxValue;
                for (int i = 0; i <= j + 1; i++)
                    best = Math.Min(best, ReferenceWeighted(p, t.Substring(i, j + 1 - i), costs));
                if (best <= k)
                    expected.Add((j, best));
            }
            Assert.That(ApproximateMatcher.FindEditEndPositions(t, p, k, costs).ToList(), Is.EqualTo(expected));
        }
    }

    #endregion

    #region Helpers

    /// <summary>Plain full-matrix weighted Wagner–Fischer (independent of the library kernels).</summary>
    private static int ReferenceWeighted(string s1, string s2, EditCosts c)
    {
        var d = new int[s1.Length + 1, s2.Length + 1];
        for (int i = 0; i <= s1.Length; i++) d[i, 0] = i * c.Deletion;
        for (int j = 0; j <= s2.Length; j++) d[0, j] = j * c.Insertion;
        for (int i = 1; i <= s1.Length; i++)
            for (int j = 1; j <= s2.Length; j++)
                d[i, j] = Math.Min(Math.Min(d[i - 1, j] + c.Deletion, d[i, j - 1] + c.Insertion),
                    d[i - 1, j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : c.Substitution));
        return d[s1.Length, s2.Length];
    }

    /// <summary>Replays an operation string; fails on an invalid script; returns the weighted cost.</summary>
    private static int WeightedReplayCost(string ops, string query, string target, EditCosts c)
    {
        int q = 0, t = 0, cost = 0;
        foreach (char op in ops)
        {
            switch (op)
            {
                case '=': Assert.That(query[q], Is.EqualTo(target[t])); q++; t++; break;
                case 'X': Assert.That(query[q], Is.Not.EqualTo(target[t])); q++; t++; cost += c.Substitution; break;
                case 'I': q++; cost += c.Deletion; break;
                case 'D': t++; cost += c.Insertion; break;
                default: Assert.Fail($"unknown op {op}"); break;
            }
        }
        Assert.That((q, t), Is.EqualTo((query.Length, target.Length)));
        return cost;
    }

    private static IEnumerable<string> Enumerate(string alphabet, int len)
    {
        if (len == 0) { yield return ""; yield break; }
        foreach (var prefix in Enumerate(alphabet, len - 1))
            foreach (char c in alphabet)
                yield return prefix + c;
    }

    private static string RandomString(Random rng, string alphabet, int len)
    {
        var chars = new char[len];
        for (int i = 0; i < len; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static string Mutate(Random rng, string s, string alphabet)
    {
        var list = s.ToList();
        int edits = rng.Next(0, Math.Max(1, list.Count / 4) + 1);
        for (int e = 0; e < edits; e++)
        {
            int p = rng.Next(0, list.Count + 1);
            int op = rng.Next(4);
            if (op == 0) list.Insert(p, alphabet[rng.Next(alphabet.Length)]);
            else if (op == 1 && p < list.Count) list.RemoveAt(p);
            else if (op == 2 && p + 1 < list.Count) (list[p], list[p + 1]) = (list[p + 1], list[p]);
            else if (p < list.Count) list[p] = alphabet[rng.Next(alphabet.Length)];
        }
        return new string(list.ToArray());
    }

    #endregion
}
