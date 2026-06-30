// 08_DIFFERENTIAL_TESTING rows 8-12 (Matching). Production method vs an INDEPENDENT oracle built from
// a different algorithm: String.IndexOf scan, brute char-by-char Hamming, recursive-memoized Levenshtein,
// a dictionary IUPAC membership table, and a dictionary-indexed PWM scan.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Matching_Differential_Tests
{
    private const double Tol = 1e-12;

    // ---- Row 8: PAT-EXACT-001 — suffix-tree FindExactMotif vs String.IndexOf scan ----

    private static List<int> IndexOfScanOracle(string seq, string motif)
    {
        var positions = new List<int>();
        for (int start = seq.IndexOf(motif, StringComparison.Ordinal);
             start >= 0;
             start = seq.IndexOf(motif, start + 1, StringComparison.Ordinal))
        {
            positions.Add(start);
        }
        return positions;
    }

    [Test]
    [Category("PAT-EXACT-001")]
    [TestCase("ACGTACGTACGT", "ACGT")]
    [TestCase("AAAA", "AA")]            // overlapping occurrences
    [TestCase("ACGTACGT", "TT")]        // no match
    [TestCase("GATTACA", "A")]          // single char, multiple positions
    public void FindExactMotif_MatchesIndexOfScan(string seq, string motif)
    {
        var actual = MotifFinder.FindExactMotif(new DnaSequence(seq), motif).ToList();
        Assert.That(actual, Is.EqualTo(IndexOfScanOracle(seq.ToUpperInvariant(), motif.ToUpperInvariant())));
    }

    // ---- Row 9: PAT-APPROX-001 — optimized Hamming vs brute char-by-char (case-insensitive) ----

    private static int HammingOracle(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++)
            if (char.ToUpperInvariant(a[i]) != char.ToUpperInvariant(b[i])) d++;
        return d;
    }

    [Test]
    [Category("PAT-APPROX-001")]
    [TestCase("ACGT", "ACGT")]   // 0
    [TestCase("ACGT", "AGGT")]   // 1
    [TestCase("ACGT", "TGCA")]   // 4
    [TestCase("acgt", "ACGT")]   // case-insensitive -> 0
    [TestCase("AAAA", "AAAT")]   // 1
    public void Hamming_MatchesBruteOracle(string a, string b)
    {
        Assert.That(ApproximateMatcher.HammingDistance(a, b), Is.EqualTo(HammingOracle(a, b)));
    }

    // ---- Row 10: PAT-APPROX-002 — two-row DP edit distance vs recursive-memoized Levenshtein ----

    private static int EditDistanceRecursiveOracle(string a, string b)
    {
        var memo = new Dictionary<(int, int), int>();
        int Rec(int i, int j)
        {
            if (i == 0) return j;
            if (j == 0) return i;
            if (memo.TryGetValue((i, j), out int cached)) return cached;
            int cost = a[i - 1] == b[j - 1] ? 0 : 1;
            int r = Math.Min(Math.Min(Rec(i - 1, j) + 1, Rec(i, j - 1) + 1), Rec(i - 1, j - 1) + cost);
            memo[(i, j)] = r;
            return r;
        }
        return Rec(a.Length, b.Length);
    }

    [Test]
    [Category("PAT-APPROX-002")]
    [TestCase("kitten", "sitting")]   // classic 3
    [TestCase("ACGT", "ACGT")]        // 0
    [TestCase("ACGT", "")]            // 4 (all deletions)
    [TestCase("", "ACGT")]            // 4 (all insertions)
    [TestCase("AAGCTA", "AAGTCA")]    // mixed indel/sub
    [TestCase("GUMBO", "GAMBOL")]     // 2
    public void EditDistance_MatchesRecursiveMemoOracle(string a, string b)
    {
        Assert.That(ApproximateMatcher.EditDistance(a, b), Is.EqualTo(EditDistanceRecursiveOracle(a, b)));
    }

    // ---- Row 11: PAT-IUPAC-001 — switch MatchesIupac vs dictionary membership table ----

    private static readonly Dictionary<char, string> IupacSets = new()
    {
        ['A'] = "A", ['C'] = "C", ['G'] = "G", ['T'] = "T", ['N'] = "ACGT",
        ['R'] = "AG", ['Y'] = "CT", ['S'] = "GC", ['W'] = "AT", ['K'] = "GT", ['M'] = "AC",
        ['B'] = "CGT", ['D'] = "AGT", ['H'] = "ACT", ['V'] = "ACG",
    };

    [Test]
    [Category("PAT-IUPAC-001")]
    public void MatchesIupac_MatchesDictionaryTable_FullMatrix()
    {
        foreach (char code in IupacSets.Keys)
        foreach (char nuc in "ACGT")
        {
            bool expected = IupacSets[code].Contains(nuc);
            Assert.That(IupacHelper.MatchesIupac(nuc, code), Is.EqualTo(expected), $"code {code} vs {nuc}");
        }
    }

    [Test]
    [Category("PAT-IUPAC-001")]
    public void MatchesIupac_InvalidCode_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => IupacHelper.MatchesIupac('A', 'Z'));
    }

    // ---- Row 12: PAT-PWM-001 — array-indexed PWM scan vs dictionary-indexed scan ----

    // Independent oracle: re-scan the sequence reading the SAME PWM via a {base->row} dictionary instead
    // of the production switch, asserting identical match positions AND scores (window/threshold/validity
    // logic is the differential).
    private static readonly Dictionary<char, int> BaseRow = new() { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };

    private static List<(int pos, double score)> PwmScanOracle(string seq, PositionWeightMatrix pwm, double threshold)
    {
        var hits = new List<(int, double)>();
        int len = pwm.Length;
        for (int i = 0; i + len <= seq.Length; i++)
        {
            double score = 0;
            bool valid = true;
            for (int j = 0; j < len; j++)
            {
                if (!BaseRow.TryGetValue(seq[i + j], out int row)) { valid = false; break; }
                score += pwm.Matrix[row, j];
            }
            if (valid && score >= threshold) hits.Add((i, score));
        }
        return hits;
    }

    [Test]
    [Category("PAT-PWM-001")]
    public void ScanWithPwm_MatchesDictionaryIndexedScan()
    {
        var pwm = MotifFinder.CreatePwm(new[] { "ACGT", "ACGA", "ACGT", "TCGT" });
        var seq = new DnaSequence("ACGTACGATCGTTTTT");

        foreach (double threshold in new[] { 0.0, -100.0, 2.0 })
        {
            var actual = MotifFinder.ScanWithPwm(seq, pwm, threshold)
                .Select(m => (m.Position, m.Score)).ToList();
            var expected = PwmScanOracle(seq.Sequence, pwm, threshold);

            Assert.That(actual.Count, Is.EqualTo(expected.Count), $"count @ threshold {threshold}");
            for (int k = 0; k < actual.Count; k++)
            {
                Assert.That(actual[k].Position, Is.EqualTo(expected[k].Item1), $"pos[{k}] @ {threshold}");
                Assert.That(actual[k].Score, Is.EqualTo(expected[k].Item2).Within(Tol), $"score[{k}] @ {threshold}");
            }
        }
    }
}
