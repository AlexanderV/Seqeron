namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property tests for the review 2026-09 B05 follow-up on PAT-APPROX-002
/// (docs/Validation/review-2026-09/B05.md): Myers bit-parallel engine, traceback/CIGAR,
/// Damerau variants.
///
/// Properties:
///   P1 (Myers == DP)      EditDistance (Myers/Hyyrö, multi-word) == Wagner–Fischer reference
///                          EditDistanceDp; FindEditEndPositions (Myers) == FindEditEndPositionsDp.
///   P2 (CIGAR replay)     GetEditAlignment's CIGAR replays query → target ('=' equal, 'X' different,
///                          consumes both strings) with unit cost == EditDistance; aligned strings minus
///                          gaps are the inputs; STANDARD CIGAR == EXTENDED with '='/'X' → 'M'.
///   P3 (hit alignment)    every FindWithEdits hit's CIGAR replays pattern → window with cost == Distance;
///                          MismatchPositions == pattern indices of 'X'; Substitution ⇔ (|window| = m and
///                          Distance = Hamming), and then MismatchPositions == Hamming mismatch indices.
///   P4 (Damerau order)    DL ≤ OSA ≤ Levenshtein; all ≥ |len difference|; 0 iff equal; symmetric.
///   P5 (DL metric)        DL satisfies the triangle inequality (Lowrance–Wagner distance is a metric).
/// Seeded random inputs (fixed seeds), deterministic.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class EditAlignmentProperties
{
    private static string RandomString(Random rng, int length, string alphabet) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static readonly string[] Alphabets = { "AC", "ACGT", "ACGTN", "acgtACGT", "ACαβ" };

    [Test]
    public void P1_Myers_EqualsDp_GlobalAndSellers()
    {
        var rng = new Random(9301);
        for (int trial = 0; trial < 1500; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            int maxLen = trial % 5 == 0 ? 260 : 40;
            string a = RandomString(rng, rng.Next(0, maxLen), alphabet);
            string b = RandomString(rng, rng.Next(0, maxLen), alphabet);
            Assert.That(ApproximateMatcher.EditDistance(a, b), Is.EqualTo(ApproximateMatcher.EditDistanceDp(a, b)),
                $"EditDistance a={a} b={b}");

            if (a.Length > 0 && b.Length > 0)
            {
                int k = rng.Next(0, a.Length + 2);
                Assert.That(ApproximateMatcher.FindEditEndPositions(b, a, k).ToList(),
                    Is.EqualTo(ApproximateMatcher.FindEditEndPositionsDp(b, a, k).ToList()),
                    $"Sellers text={b} pattern={a} k={k}");
            }
        }
    }

    [Test]
    public void P2_GetEditAlignment_CigarReplaysWithCostEqualToDistance()
    {
        var rng = new Random(9302);
        for (int trial = 0; trial < 1000; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            string q = RandomString(rng, rng.Next(0, 30), alphabet);
            string t = RandomString(rng, rng.Next(0, 30), alphabet);
            var a = ApproximateMatcher.GetEditAlignment(q, t);

            Assert.That(a.Distance, Is.EqualTo(ApproximateMatcher.EditDistance(q, t)));
            Assert.That(Unit.Alignment.ApproximateMatcher_EditAlignment_Tests.ReplayCost(a.Cigar, q, t), Is.EqualTo(a.Distance));
            Assert.That(a.AlignedQuery.Replace("-", ""), Is.EqualTo(q.Replace("-", "")));
            Assert.That(a.AlignedTarget.Replace("-", ""), Is.EqualTo(t.Replace("-", "")));
            Assert.That(a.AlignedQuery.Length, Is.EqualTo(a.Operations.Length));
            Assert.That(a.StandardCigar, Is.EqualTo(RunLength(a.Operations.Replace('=', 'M').Replace('X', 'M'))));
        }
    }

    [Test]
    public void P3_FindWithEdits_HitAlignmentsReplayToWindow()
    {
        var rng = new Random(9303);
        int hits = 0;
        for (int trial = 0; trial < 400; trial++)
        {
            string text = RandomString(rng, rng.Next(1, 50), "ACGT");
            string pattern = RandomString(rng, rng.Next(1, 9), "ACGT");
            int k = rng.Next(0, 4);

            foreach (var h in ApproximateMatcher.FindWithEdits(text, pattern, k))
            {
                hits++;
                var al = h.Alignment!;
                Assert.That(Unit.Alignment.ApproximateMatcher_EditAlignment_Tests.ReplayCost(al.Cigar, pattern, h.MatchedSequence),
                    Is.EqualTo(h.Distance));

                var xPositions = new List<int>();
                int qi = 0;
                foreach (char op in al.Operations)
                {
                    if (op == 'X') xPositions.Add(qi);
                    if (op != 'D') qi++;
                }
                Assert.That(h.MismatchPositions, Is.EqualTo(xPositions));

                bool hammingEqual = h.MatchedSequence.Length == pattern.Length
                    && h.Distance == ApproximateMatcher.HammingDistance(pattern, h.MatchedSequence);
                Assert.That(h.MismatchType == MismatchType.Substitution, Is.EqualTo(hammingEqual),
                    $"text={text} pattern={pattern} window={h.MatchedSequence}");
                if (hammingEqual)
                {
                    var hamming = Enumerable.Range(0, pattern.Length)
                        .Where(j => pattern[j] != h.MatchedSequence[j]).ToList();
                    Assert.That(h.MismatchPositions, Is.EqualTo(hamming));
                }
            }
        }
        Assert.That(hits, Is.GreaterThan(1000), "the property must be exercised on many hits");
    }

    [Test]
    public void P4_DamerauOrdering_DlLeOsaLeLevenshtein()
    {
        var rng = new Random(9304);
        for (int trial = 0; trial < 3000; trial++)
        {
            string alphabet = Alphabets[rng.Next(Alphabets.Length)];
            string a = RandomString(rng, rng.Next(0, 14), alphabet);
            string b = RandomString(rng, rng.Next(0, 14), alphabet);
            int lev = ApproximateMatcher.EditDistance(a, b);
            int osa = ApproximateMatcher.OptimalStringAlignmentDistance(a, b);
            int dl = ApproximateMatcher.DamerauLevenshteinDistance(a, b);

            Assert.That(dl, Is.LessThanOrEqualTo(osa), $"DL ≤ OSA a={a} b={b}");
            Assert.That(osa, Is.LessThanOrEqualTo(lev), $"OSA ≤ Lev a={a} b={b}");
            Assert.That(dl, Is.GreaterThanOrEqualTo(Math.Abs(a.Length - b.Length)));
            Assert.That(dl == 0, Is.EqualTo(a == b));
            Assert.That(osa, Is.EqualTo(ApproximateMatcher.OptimalStringAlignmentDistance(b, a)));
            Assert.That(dl, Is.EqualTo(ApproximateMatcher.DamerauLevenshteinDistance(b, a)));
        }
    }

    [Test]
    public void P5_DamerauLevenshtein_TriangleInequality()
    {
        var rng = new Random(9305);
        for (int trial = 0; trial < 2000; trial++)
        {
            string a = RandomString(rng, rng.Next(0, 9), "ABC");
            string b = RandomString(rng, rng.Next(0, 9), "ABC");
            string c = RandomString(rng, rng.Next(0, 9), "ABC");
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, c),
                Is.LessThanOrEqualTo(ApproximateMatcher.DamerauLevenshteinDistance(a, b)
                                     + ApproximateMatcher.DamerauLevenshteinDistance(b, c)),
                $"a={a} b={b} c={c}");
        }
        // OSA is not a metric: the classic counterexample (Boytsov 2011).
        Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance("CA", "ABC"),
            Is.GreaterThan(ApproximateMatcher.OptimalStringAlignmentDistance("CA", "AC")
                           + ApproximateMatcher.OptimalStringAlignmentDistance("AC", "ABC")));
    }

    private static string RunLength(string ops)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < ops.Length;)
        {
            int run = 1;
            while (i + run < ops.Length && ops[i + run] == ops[i]) run++;
            sb.Append(run).Append(ops[i]);
            i += run;
        }
        return sb.ToString();
    }
}
