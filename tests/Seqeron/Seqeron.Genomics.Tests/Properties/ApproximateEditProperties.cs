namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property tests for the Sellers (1980) k-differences end-position search
/// <see cref="ApproximateMatcher.FindEditEndPositions(string, string, int)"/> added by
/// review 2026-09 B05 F2 (PAT-APPROX-002).
///
/// Definition (Sellers 1980; Navarro 2001 §5.1): end position j is reported with distance
/// C[m, j] = min over i ≤ j+1 of ed(P, T[i..j]) (free start in the text, the empty substring
/// included) iff C[m, j] ≤ k. The oracle below is an independent brute-force Wagner–Fischer
/// over every substring ending at j — it does not reuse the library's DP.
///
/// Properties:
///   P1 (oracle)     end set and distances == brute-force minimum over all substrings ending at j.
///   P2 (monotone)   ends(k) ⊆ ends(k+1) with identical distances (C[m, j] does not depend on k).
///   P3 (k = 0)      ends(0) == { i + m − 1 : i an exact occurrence }, all distances 0.
///   P4 (case)       lower-casing text and pattern does not change the output (inputs are upper-cased).
///
/// Inputs are seeded-random (fixed seeds), bounded sizes, so the fixture is fast and deterministic.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class ApproximateEditProperties
{
    private const int Cases = 250;

    private static string RandomString(Random rng, int length, string alphabet)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++)
            c[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(c);
    }

    /// <summary>Independent unit-cost Levenshtein distance (Wagner–Fischer 1974).</summary>
    private static int Levenshtein(string a, string b)
    {
        var prev = new int[b.Length + 1];
        var curr = new int[b.Length + 1];
        for (int j = 0; j <= b.Length; j++) prev[j] = j;
        for (int i = 1; i <= a.Length; i++)
        {
            curr[0] = i;
            for (int j = 1; j <= b.Length; j++)
            {
                int sub = prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1);
                curr[j] = Math.Min(sub, Math.Min(prev[j] + 1, curr[j - 1] + 1));
            }
            (prev, curr) = (curr, prev);
        }
        return prev[b.Length];
    }

    private static List<(int EndPosition, int Distance)> BruteForceEnds(string text, string pattern, int k)
    {
        var result = new List<(int, int)>();
        for (int j = 0; j < text.Length; j++)
        {
            int best = int.MaxValue;
            for (int i = 0; i <= j + 1; i++) // i = j+1 → empty substring
                best = Math.Min(best, Levenshtein(pattern, text.Substring(i, j - i + 1)));
            if (best <= k)
                result.Add((j, best));
        }
        return result;
    }

    private static IEnumerable<(string Text, string Pattern, int K, int Seed)> RandomCases()
    {
        for (int seed = 0; seed < Cases; seed++)
        {
            var rng = new Random(20260930 + seed);
            string alphabet = seed % 3 == 0 ? "AC" : "ACGT";
            string text = RandomString(rng, rng.Next(1, 26), alphabet);
            string pattern = rng.Next(3) == 0 && text.Length > 2
                ? string.Concat(text.AsSpan(rng.Next(text.Length - 1), 2), RandomString(rng, rng.Next(0, 4), alphabet))
                : RandomString(rng, rng.Next(1, 7), alphabet);
            yield return (text, pattern, rng.Next(0, 4), seed);
        }
    }

    [Test]
    public void FindEditEndPositions_EqualsBruteForceMinimumOverSubstringsEndingAtJ()
    {
        foreach (var (text, pattern, k, seed) in RandomCases())
        {
            var actual = ApproximateMatcher.FindEditEndPositions(text, pattern, k).ToList();
            var expected = BruteForceEnds(text, pattern, k);
            Assert.That(actual, Is.EqualTo(expected),
                $"seed={seed} text={text} pattern={pattern} k={k}");
        }
    }

    [Test]
    public void FindEditEndPositions_IsMonotoneInK_DistancesIndependentOfK()
    {
        foreach (var (text, pattern, _, seed) in RandomCases())
        {
            var previous = new Dictionary<int, int>();
            for (int k = 0; k <= pattern.Length + 1; k++)
            {
                var current = ApproximateMatcher.FindEditEndPositions(text, pattern, k)
                    .ToDictionary(e => e.EndPosition, e => e.Distance);
                foreach (var (end, d) in previous)
                {
                    Assert.That(current.TryGetValue(end, out int dk), Is.True,
                        $"seed={seed} end {end} reported for k={k - 1} but not for k={k}");
                    Assert.That(dk, Is.EqualTo(d), $"seed={seed} end {end}: distance changed with k");
                }
                Assert.That(current.Values.All(d => d <= k), Is.True, $"seed={seed} k={k}: distance > k");
                previous = current;
            }
            // k ≥ m: the empty substring is within m of the pattern, so every end position qualifies.
            Assert.That(previous.Keys.OrderBy(x => x), Is.EqualTo(Enumerable.Range(0, text.Length)),
                $"seed={seed}: k ≥ m must report every end position");
        }
    }

    [Test]
    public void FindEditEndPositions_KZero_EqualsExactOccurrenceEnds()
    {
        foreach (var (text, pattern, _, seed) in RandomCases())
        {
            var expected = new List<(int, int)>();
            for (int i = 0; i + pattern.Length <= text.Length; i++)
                if (string.CompareOrdinal(text, i, pattern, 0, pattern.Length) == 0)
                    expected.Add((i + pattern.Length - 1, 0));

            Assert.That(ApproximateMatcher.FindEditEndPositions(text, pattern, 0).ToList(),
                Is.EqualTo(expected), $"seed={seed} text={text} pattern={pattern}");
        }
    }

    [Test]
    public void FindEditEndPositions_IsCaseInsensitive()
    {
        foreach (var (text, pattern, k, seed) in RandomCases())
        {
            Assert.That(
                ApproximateMatcher.FindEditEndPositions(text.ToLowerInvariant(), pattern.ToLowerInvariant(), k).ToList(),
                Is.EqualTo(ApproximateMatcher.FindEditEndPositions(text, pattern, k).ToList()),
                $"seed={seed}");
        }
    }
}
