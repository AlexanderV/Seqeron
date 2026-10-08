namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the review 2026-09 B05 audit additions on PAT-APPROX-002
/// (docs/Validation/review-2026-09/B05.md F27 weighted edit costs / linear-space DL, F34 weighted Damerau
/// distances + transposition-aware traceback).
///
/// Properties (independent in-test oracles: weighted Wagner–Fischer 1974 and the Lowrance–Wagner 1975
/// recurrence with W_I, W_D, W_C, W_S; Sellers 1980 free-start search by brute force over substrings):
///   W1 (oracle)       EditDistance(a, b, costs) == weighted Wagner–Fischer; unit costs == EditDistance(a, b);
///                      uniform c == c·EditDistance(a, b).
///   W2 (swap)         EditDistance(a, b, (I, D, S)) == EditDistance(b, a, (D, I, S)) (an insertion of a → b
///                      is a deletion of b → a); likewise for the weighted OSA / DL distances.
///   W3 (scaling)      every cost × c ⇒ distance × c (Levenshtein, OSA, DL).
///   W4 (replay)       GetEditAlignment / GetEditAlignmentLinearSpace with costs: the operation string replays
///                      query → target and its cost Σ X·S + I·D + D·I == Distance == EditDistance(costs).
///   W5 (Sellers)      weighted FindEditEndPositions == brute-force min over substrings ending at j.
///   W6 (Damerau)      weighted OSA / DL == independent oracles; unit DamerauCosts == the unit engines
///                      (DL = Zhao–Sahni linear-space engine, F27); DL ≤ OSA ≤ weighted Levenshtein.
///   W7 (DL replay)    GetDamerauLevenshteinAlignment / GetOptimalStringAlignment: the operations tile s1 and
///                      s2 contiguously, replay s1 → s2 (transposition block a_k … a_i → b_l … b_j with
///                      a_k = b_j, a_i = b_l), Σ Cost == Distance == the distance; OSA blocks are adjacent
///                      swaps (SourceLength = TargetLength = 2).
///   W8 (contract)     DL with 2·W_S &lt; W_I + W_D throws ArgumentException (Lowrance–Wagner condition);
///                      OSA accepts any non-negative costs.
/// Seeded random inputs (fixed seeds), bounded sizes, deterministic.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class WeightedEditDamerauProperties
{
    private static readonly string[] Alphabets = { "AB", "ACGT", "ACGTN", "ABC" };

    private static string RandomString(Random rng, int length, string alphabet) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static EditCosts RandomCosts(Random rng, int max = 6) =>
        new(rng.Next(0, max), rng.Next(0, max), rng.Next(0, max));

    /// <summary>Weighted Wagner–Fischer (1974): D[i,0] = i·W_D, D[0,j] = j·W_I.</summary>
    private static long WeightedLevenshtein(string a, string b, int ins, int del, int sub)
    {
        var d = new long[a.Length + 1, b.Length + 1];
        for (int i = 0; i <= a.Length; i++) d[i, 0] = (long)i * del;
        for (int j = 0; j <= b.Length; j++) d[0, j] = (long)j * ins;
        for (int i = 1; i <= a.Length; i++)
            for (int j = 1; j <= b.Length; j++)
                d[i, j] = Math.Min(d[i - 1, j - 1] + (a[i - 1] == b[j - 1] ? 0 : sub),
                    Math.Min(d[i - 1, j] + del, d[i, j - 1] + ins));
        return d[a.Length, b.Length];
    }

    /// <summary>Optimal string alignment (restricted edit): Levenshtein + adjacent swap a_{i-1}a_i = b_j b_{j-1}.</summary>
    private static long WeightedOsa(string a, string b, DamerauCosts c)
    {
        var d = new long[a.Length + 1, b.Length + 1];
        for (int i = 0; i <= a.Length; i++) d[i, 0] = (long)i * c.Deletion;
        for (int j = 0; j <= b.Length; j++) d[0, j] = (long)j * c.Insertion;
        for (int i = 1; i <= a.Length; i++)
            for (int j = 1; j <= b.Length; j++)
            {
                long v = Math.Min(d[i - 1, j - 1] + (a[i - 1] == b[j - 1] ? 0 : c.Substitution),
                    Math.Min(d[i - 1, j] + c.Deletion, d[i, j - 1] + c.Insertion));
                if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1])
                    v = Math.Min(v, d[i - 2, j - 2] + c.Transposition);
                d[i, j] = v;
            }
        return d[a.Length, b.Length];
    }

    /// <summary>
    /// Lowrance &amp; Wagner (1975) extended edit distance, direct O(m²n²) form of their recurrence:
    /// H[i,j] = min(…, H[k−1, l−1] + (i−k−1)·W_D + W_S + (j−l−1)·W_I) over every k &lt; i, l &lt; j with
    /// a_k = b_j and a_i = b_l (valid when 2·W_S ≥ W_I + W_D).
    /// </summary>
    private static long WeightedDl(string a, string b, DamerauCosts c)
    {
        var h = new long[a.Length + 1, b.Length + 1];
        for (int i = 0; i <= a.Length; i++) h[i, 0] = (long)i * c.Deletion;
        for (int j = 0; j <= b.Length; j++) h[0, j] = (long)j * c.Insertion;
        for (int i = 1; i <= a.Length; i++)
            for (int j = 1; j <= b.Length; j++)
            {
                long v = Math.Min(h[i - 1, j - 1] + (a[i - 1] == b[j - 1] ? 0 : c.Substitution),
                    Math.Min(h[i - 1, j] + c.Deletion, h[i, j - 1] + c.Insertion));
                for (int k = 1; k < i; k++)
                {
                    if (a[k - 1] != b[j - 1]) continue;
                    for (int l = 1; l < j; l++)
                        if (a[i - 1] == b[l - 1])
                            v = Math.Min(v, h[k - 1, l - 1] + (long)(i - k - 1) * c.Deletion + c.Transposition + (long)(j - l - 1) * c.Insertion);
                }
                h[i, j] = v;
            }
        return h[a.Length, b.Length];
    }

    [Test]
    public void W1_WeightedDistance_EqualsOracle_UnitAndUniformCosts()
    {
        var rng = new Random(270_001);
        for (int trial = 0; trial < 1500; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            int maxLen = trial % 10 == 0 ? 90 : 18;
            string a = RandomString(rng, rng.Next(0, maxLen), alphabet);
            string b = RandomString(rng, rng.Next(0, maxLen), alphabet);
            var costs = RandomCosts(rng);
            Assert.That(ApproximateMatcher.EditDistance(a, b, costs),
                Is.EqualTo(WeightedLevenshtein(a, b, costs.Insertion, costs.Deletion, costs.Substitution)), $"{a}|{b}|{costs}");
            Assert.That(ApproximateMatcher.EditDistance(a, b, costs.Insertion, costs.Deletion, costs.Substitution),
                Is.EqualTo(ApproximateMatcher.EditDistance(a, b, costs)));
            Assert.That(ApproximateMatcher.EditDistance(a, b, EditCosts.Unit), Is.EqualTo(ApproximateMatcher.EditDistance(a, b)));
            int c = rng.Next(0, 9);
            Assert.That(ApproximateMatcher.EditDistance(a, b, new EditCosts(c, c, c)), Is.EqualTo(c * ApproximateMatcher.EditDistance(a, b)));
        }
    }

    [Test]
    public void W2_W3_SwapAndScaling_AllWeightedDistances()
    {
        var rng = new Random(270_002);
        for (int trial = 0; trial < 1000; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            string a = RandomString(rng, rng.Next(0, 16), alphabet);
            string b = RandomString(rng, rng.Next(0, 16), alphabet);
            var e = RandomCosts(rng);
            var swapped = new EditCosts(e.Deletion, e.Insertion, e.Substitution);
            int t = rng.Next(0, 8);
            var dc = new DamerauCosts(e, t);
            var dcSwapped = new DamerauCosts(swapped, t);
            int scale = rng.Next(0, 5);
            var eScaled = new EditCosts(e.Insertion * scale, e.Deletion * scale, e.Substitution * scale);
            var dcScaled = new DamerauCosts(eScaled, t * scale);
            string msg = $"{a}|{b}|{dc}";

            Assert.That(ApproximateMatcher.EditDistance(b, a, swapped), Is.EqualTo(ApproximateMatcher.EditDistance(a, b, e)), msg);
            Assert.That(ApproximateMatcher.EditDistance(a, b, eScaled), Is.EqualTo(scale * ApproximateMatcher.EditDistance(a, b, e)), msg);
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(b, a, dcSwapped),
                Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dc)), msg);
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dcScaled),
                Is.EqualTo(scale * ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dc)), msg);
            if (2 * t >= e.Insertion + e.Deletion)
            {
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(b, a, dcSwapped),
                    Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc)), msg);
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, dcScaled),
                    Is.EqualTo(scale * ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc)), msg);
            }
        }
    }

    private static int ReplayCost(EditAlignment al, string query, string target, EditCosts c)
    {
        int q = 0, t = 0, cost = 0;
        foreach (char op in al.Operations)
        {
            switch (op)
            {
                case '=': Assert.That(query[q], Is.EqualTo(target[t])); q++; t++; break;
                case 'X': Assert.That(query[q], Is.Not.EqualTo(target[t])); q++; t++; cost += c.Substitution; break;
                case 'I': q++; cost += c.Deletion; break;   // query character with no target counterpart
                case 'D': t++; cost += c.Insertion; break;  // target character with no query counterpart
                default: Assert.Fail($"unexpected op {op}"); break;
            }
        }
        Assert.That(q, Is.EqualTo(query.Length));
        Assert.That(t, Is.EqualTo(target.Length));
        Assert.That(al.AlignedQuery.Replace("-", ""), Is.EqualTo(query));
        Assert.That(al.AlignedTarget.Replace("-", ""), Is.EqualTo(target));
        return cost;
    }

    [Test]
    public void W4_WeightedAlignments_ReplayWithCostEqualToDistance()
    {
        var rng = new Random(270_004);
        for (int trial = 0; trial < 800; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            int maxLen = trial % 10 == 0 ? 70 : 20;
            string q = RandomString(rng, rng.Next(0, maxLen), alphabet);
            string t = RandomString(rng, rng.Next(0, maxLen), alphabet);
            var c = RandomCosts(rng);
            int d = ApproximateMatcher.EditDistance(q, t, c);
            var full = ApproximateMatcher.GetEditAlignment(q, t, c);
            var linear = ApproximateMatcher.GetEditAlignmentLinearSpace(q, t, c);
            Assert.That(full.Distance, Is.EqualTo(d), $"{q}|{t}|{c}");
            Assert.That(linear.Distance, Is.EqualTo(d), $"{q}|{t}|{c}");
            Assert.That(ReplayCost(full, q, t, c), Is.EqualTo(d), $"full {q}|{t}|{c}");
            Assert.That(ReplayCost(linear, q, t, c), Is.EqualTo(d), $"linear {q}|{t}|{c}");
            Assert.That(ApproximateMatcher.GetEditAlignment(q, t, EditCosts.Unit), Is.EqualTo(ApproximateMatcher.GetEditAlignment(q, t)));
            Assert.That(ApproximateMatcher.GetEditAlignmentLinearSpace(q, t, EditCosts.Unit),
                Is.EqualTo(ApproximateMatcher.GetEditAlignmentLinearSpace(q, t)));
        }
    }

    [Test]
    public void W5_WeightedSellers_EqualsBruteForceOverSubstrings()
    {
        var rng = new Random(270_005);
        for (int trial = 0; trial < 300; trial++)
        {
            string alphabet = trial % 3 == 0 ? "AC" : "ACGT";
            string text = RandomString(rng, rng.Next(1, 22), alphabet);
            string pattern = RandomString(rng, rng.Next(1, 7), alphabet);
            var c = RandomCosts(rng, 5);
            int maxCost = rng.Next(0, 12);
            var expected = new List<(int, int)>();
            for (int j = 0; j < text.Length; j++)
            {
                long best = long.MaxValue;
                for (int i = 0; i <= j + 1; i++)
                    best = Math.Min(best, WeightedLevenshtein(pattern, text.Substring(i, j - i + 1), c.Insertion, c.Deletion, c.Substitution));
                if (best <= maxCost)
                    expected.Add((j, (int)best));
            }
            Assert.That(ApproximateMatcher.FindEditEndPositions(text, pattern, maxCost, c).ToList(), Is.EqualTo(expected),
                $"{text}|{pattern}|{maxCost}|{c}");
            if (c == EditCosts.Unit)
                Assert.That(ApproximateMatcher.FindEditEndPositions(text, pattern, maxCost, c).ToList(),
                    Is.EqualTo(ApproximateMatcher.FindEditEndPositions(text, pattern, maxCost).ToList()));
        }
    }

    [Test]
    public void W6_WeightedDamerau_EqualsOracles_UnitEqualsEngines_AndOrdering()
    {
        var rng = new Random(270_006);
        for (int trial = 0; trial < 1200; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            int maxLen = trial % 12 == 0 ? 30 : 9;
            string a = RandomString(rng, rng.Next(0, maxLen), alphabet);
            string b = RandomString(rng, rng.Next(0, maxLen), alphabet);
            var e = RandomCosts(rng);
            int t = rng.Next(0, 9);
            var dc = new DamerauCosts(e, t);
            string msg = $"{a}|{b}|{dc}";

            int osa = ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dc);
            Assert.That(osa, Is.EqualTo(WeightedOsa(a, b, dc)), msg);
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, e, t), Is.EqualTo(osa), msg);
            int lev = ApproximateMatcher.EditDistance(a, b, e);
            Assert.That(osa, Is.LessThanOrEqualTo(lev), msg);

            if (2 * t >= e.Insertion + e.Deletion)
            {
                int dl = ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc);
                if (maxLen <= 9)
                    Assert.That(dl, Is.EqualTo(WeightedDl(a, b, dc)), msg);
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, e, t), Is.EqualTo(dl), msg);
                Assert.That(dl, Is.LessThanOrEqualTo(osa), msg);
            }

            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, DamerauCosts.Unit),
                Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(a, b)), msg);
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b, DamerauCosts.Unit),
                Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistance(a, b)), msg);
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b),
                Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistanceFullMatrix(a, b)), msg);
            if (maxLen <= 9)
                Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b), Is.EqualTo(WeightedDl(a, b, DamerauCosts.Unit)), msg);
        }
    }

    private static int ReplayDamerau(DamerauAlignment al, string s1, string s2, DamerauCosts c, bool osa)
    {
        int src = 0, tgt = 0, cost = 0, transpositions = 0;
        foreach (var op in al.Operations)
        {
            Assert.That(op.SourcePosition, Is.EqualTo(src));
            Assert.That(op.TargetPosition, Is.EqualTo(tgt));
            switch (op.Kind)
            {
                case DamerauOperationKind.Match:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((1, 1)));
                    Assert.That(s1[src], Is.EqualTo(s2[tgt]));
                    Assert.That(op.Cost, Is.Zero);
                    break;
                case DamerauOperationKind.Substitution:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((1, 1)));
                    Assert.That(s1[src], Is.Not.EqualTo(s2[tgt]));
                    Assert.That(op.Cost, Is.EqualTo(c.Substitution));
                    break;
                case DamerauOperationKind.Insertion:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((0, 1)));
                    Assert.That(op.Cost, Is.EqualTo(c.Insertion));
                    break;
                case DamerauOperationKind.Deletion:
                    Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((1, 0)));
                    Assert.That(op.Cost, Is.EqualTo(c.Deletion));
                    break;
                case DamerauOperationKind.Transposition:
                    transpositions++;
                    Assert.That(op.SourceLength, Is.GreaterThanOrEqualTo(2));
                    Assert.That(op.TargetLength, Is.GreaterThanOrEqualTo(2));
                    if (osa)
                        Assert.That((op.SourceLength, op.TargetLength), Is.EqualTo((2, 2)));
                    int sEnd = src + op.SourceLength - 1, tEnd = tgt + op.TargetLength - 1;
                    Assert.That(s1[src], Is.EqualTo(s2[tEnd]), "a_k = b_j");
                    Assert.That(s1[sEnd], Is.EqualTo(s2[tgt]), "a_i = b_l");
                    Assert.That(op.Cost, Is.EqualTo((op.SourceLength - 2) * c.Deletion + c.Transposition + (op.TargetLength - 2) * c.Insertion));
                    break;
            }
            src += op.SourceLength;
            tgt += op.TargetLength;
            cost += op.Cost;
        }
        Assert.That(src, Is.EqualTo(s1.Length));
        Assert.That(tgt, Is.EqualTo(s2.Length));
        Assert.That(al.TranspositionCount, Is.EqualTo(transpositions));
        return cost;
    }

    [Test]
    public void W7_DamerauAlignments_ReplayWithCostEqualToDistance()
    {
        var rng = new Random(270_007);
        for (int trial = 0; trial < 1000; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            int maxLen = trial % 10 == 0 ? 40 : 12;
            string a = RandomString(rng, rng.Next(0, maxLen), alphabet);
            string b = RandomString(rng, rng.Next(0, maxLen), alphabet);
            var e = RandomCosts(rng);
            int t = rng.Next(0, 9);
            var dc = new DamerauCosts(e, t);
            string msg = $"{a}|{b}|{dc}";

            var osa = ApproximateMatcher.GetOptimalStringAlignment(a, b, dc);
            Assert.That(osa.Distance, Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dc)), msg);
            Assert.That(ReplayDamerau(osa, a, b, dc, osa: true), Is.EqualTo(osa.Distance), msg);

            if (2 * t >= e.Insertion + e.Deletion)
            {
                var dl = ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b, dc);
                Assert.That(dl.Distance, Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc)), msg);
                Assert.That(ReplayDamerau(dl, a, b, dc, osa: false), Is.EqualTo(dl.Distance), msg);
            }

            var unitDl = ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b);
            Assert.That(unitDl, Is.EqualTo(ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b, DamerauCosts.Unit)), msg);
            Assert.That(unitDl.Distance, Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistance(a, b)), msg);
            Assert.That(ApproximateMatcher.GetOptimalStringAlignment(a, b).Distance,
                Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(a, b)), msg);
        }
    }

    [Test]
    public void W8_DlExactnessCondition_ThrowsOnlyWhenViolated()
    {
        var rng = new Random(270_008);
        for (int trial = 0; trial < 400; trial++)
        {
            var e = RandomCosts(rng, 8);
            int t = rng.Next(0, 8);
            var dc = new DamerauCosts(e, t);
            string a = RandomString(rng, rng.Next(0, 8), "ACG");
            string b = RandomString(rng, rng.Next(0, 8), "ACG");
            Assert.That(() => ApproximateMatcher.OptimalStringAlignmentDistance(a, b, dc), NUnit.Framework.Throws.Nothing);
            if (2 * t < e.Insertion + e.Deletion)
            {
                Assert.That(() => ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc), NUnit.Framework.Throws.ArgumentException, dc.ToString());
                Assert.That(() => ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b, dc), NUnit.Framework.Throws.ArgumentException, dc.ToString());
            }
            else
                Assert.That(() => ApproximateMatcher.DamerauLevenshteinDistance(a, b, dc), NUnit.Framework.Throws.Nothing, dc.ToString());
        }
    }
}
