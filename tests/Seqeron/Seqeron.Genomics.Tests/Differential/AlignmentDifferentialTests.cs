// 08_DIFFERENTIAL_TESTING rows 35-37 (pairwise alignment). The Needleman-Wunsch / Smith-Waterman /
// fitting DP scores are checked against a BRUTE-FORCE oracle that enumerates every possible alignment
// of short sequences (recursively trying a match/mismatch column or a gap in either sequence) and takes
// the maximum score — an O(exponential) method structurally independent of the production DP.
// SimpleDna scoring: match +1, mismatch -1, linear gap -1.

using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class AlignmentDifferentialTests
{
    private const int Match = 1, Mismatch = -1, Gap = -1;

    // Brute-force global (Needleman-Wunsch) score: max over all alignments of a and b.
    // No memoization (kept thread-safe for parallel test runs); sequences are short enough.
    private static int BruteGlobal(string a, string b)
    {
        if (a.Length == 0) return b.Length * Gap;
        if (b.Length == 0) return a.Length * Gap;
        int diag = BruteGlobal(a[1..], b[1..]) + (a[0] == b[0] ? Match : Mismatch);
        int gapA = BruteGlobal(a[1..], b) + Gap;   // a[0] aligned to a gap
        int gapB = BruteGlobal(a, b[1..]) + Gap;   // b[0] aligned to a gap
        return System.Math.Max(diag, System.Math.Max(gapA, gapB));
    }

    // Smith-Waterman: best global score over any pair of substrings, floored at 0.
    private static int BruteLocal(string a, string b)
    {
        int best = 0;
        for (int i1 = 0; i1 <= a.Length; i1++)
            for (int i2 = i1; i2 <= a.Length; i2++)
                for (int j1 = 0; j1 <= b.Length; j1++)
                    for (int j2 = j1; j2 <= b.Length; j2++)
                        best = System.Math.Max(best, BruteGlobal(a[i1..i2], b[j1..j2]));
        return best;
    }

    // Fitting (semi-global, free end gaps on seq2): align all of seq1 to the best substring of seq2.
    private static int BruteSemi(string seq1, string seq2)
    {
        int best = seq1.Length * Gap; // seq1 vs empty substring
        for (int j1 = 0; j1 <= seq2.Length; j1++)
            for (int j2 = j1; j2 <= seq2.Length; j2++)
                best = System.Math.Max(best, BruteGlobal(seq1, seq2[j1..j2]));
        return best;
    }

    // ---- Row 35: ALIGN-GLOBAL-001 — GlobalAlign (NW) vs brute enumeration ----

    [Test]
    [Category("ALIGN-GLOBAL-001")]
    [TestCase("ACGT", "ACGT")]
    [TestCase("ACGT", "AGT")]
    [TestCase("AAAC", "AGC")]
    [TestCase("GATTAC", "GATCAC")]
    [TestCase("ACGTA", "TGCAT")]
    [TestCase("AAA", "AAAAA")]
    public void GlobalAlign_ScoreMatchesBruteForce(string s1, string s2)
    {
        var r = SequenceAligner.GlobalAlign(new DnaSequence(s1), new DnaSequence(s2));
        Assert.That(r.Score, Is.EqualTo(BruteGlobal(s1, s2)));
    }

    // ---- Row 36: ALIGN-LOCAL-001 — LocalAlign (SW) vs brute over substrings ----

    [Test]
    [Category("ALIGN-LOCAL-001")]
    [TestCase("ACGT", "TACGTA")]
    [TestCase("AAGGTT", "GGTT")]
    [TestCase("ACGTA", "TGCAT")]
    [TestCase("GATTAC", "TTACG")]
    [TestCase("CCCCC", "AAAAA")]   // no positive local alignment -> 0
    public void LocalAlign_ScoreMatchesBruteForce(string s1, string s2)
    {
        var r = SequenceAligner.LocalAlign(new DnaSequence(s1), new DnaSequence(s2));
        Assert.That(r.Score, Is.EqualTo(BruteLocal(s1, s2)));
    }

    // ---- Row 37: ALIGN-SEMI-001 — SemiGlobalAlign (fitting) vs brute over seq2 substrings ----

    [Test]
    [Category("ALIGN-SEMI-001")]
    [TestCase("ACG", "TTACGTT")]
    [TestCase("GATT", "AAGATTCC")]
    [TestCase("ACGT", "ACGT")]
    [TestCase("AAA", "GGGAAAGGG")]
    [TestCase("ACGT", "TTTT")]
    public void SemiGlobalAlign_ScoreMatchesBruteForce(string s1, string s2)
    {
        var r = SequenceAligner.SemiGlobalAlign(new DnaSequence(s1), new DnaSequence(s2));
        Assert.That(r.Score, Is.EqualTo(BruteSemi(s1, s2)));
    }
}
