// 08_DIFFERENTIAL_TESTING rows 13-17 (Repeats). Production finders vs INDEPENDENT oracles:
// regex homopolymer detection, the second in-project tandem implementation (DUAL), a brute
// reverse-complement search, a brute substring-equality scan, and an independent reverse-complement
// palindrome test. The reverse complement used by the oracles is computed from a literal IUPAC table,
// not from DnaSequence.GetReverseComplementString.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class RepeatsDifferentialTests
{
    private static readonly Dictionary<char, char> Comp = new()
    {
        ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G',
    };

    private static string RevComp(string s)
    {
        var arr = s.Select(c => Comp[c]).ToArray();
        Array.Reverse(arr);
        return new string(arr);
    }

    // ---- Row 13: REP-STR-001 — FindMicrosatellites vs regex homopolymer/run oracle ----

    [Test]
    [Category("REP-STR-001")]
    public void Microsatellites_Mononucleotide_MatchesRegexRunOracle()
    {
        const string seq = "GGGAAAAAGGGCT"; // GGG, AAAAA, GGG runs (>=3)
        var actual = RepeatFinder.FindMicrosatellites(seq, minUnitLength: 1, maxUnitLength: 1, minRepeats: 3)
            .Select(m => (m.Position, m.RepeatUnit, m.RepeatCount)).ToList();

        // Independent regex oracle: any base repeated >= 3 times.
        var expected = Regex.Matches(seq, @"([ACGT])\1{2,}")
            .Select(m => (m.Index, m.Value.Substring(0, 1), m.Value.Length)).ToList();

        Assert.That(actual, Is.EqualTo(expected));
    }

    [Test]
    [Category("REP-STR-001")]
    public void Microsatellites_Dinucleotide_MatchesRegexRunOracle()
    {
        const string seq = "GCATATATGC"; // (AT)x3 at position 2
        var actual = RepeatFinder.FindMicrosatellites(seq, minUnitLength: 2, maxUnitLength: 2, minRepeats: 3)
            .Select(m => (m.Position, m.RepeatUnit, m.RepeatCount)).ToList();

        var expected = Regex.Matches(seq, @"(?:AT){3,}")
            .Select(m => (m.Index, "AT", m.Value.Length / 2)).ToList();

        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 14: REP-TANDEM-001 — GenomicAnalyzer.FindTandemRepeats vs RepeatFinder (DUAL) ----

    // On an input with a single unambiguous tandem the two independent in-project implementations must
    // agree on (unit, start, count).
    [Test]
    [Category("REP-TANDEM-001")]
    [TestCase("GCATATATGC", "AT", 2, 3)]
    [TestCase("TTCAGCAGCAGTT", "CAG", 2, 3)]
    public void TandemRepeats_GenomicAnalyzer_AgreesWithRepeatFinder(string seq, string unit, int start, int count)
    {
        var ga = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), minUnitLength: unit.Length, minRepetitions: count)
            .Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
        var rf = RepeatFinder.FindMicrosatellites(seq, minUnitLength: unit.Length, maxUnitLength: unit.Length, minRepeats: count)
            .Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();

        var expected = new List<(string, int, int)> { (unit, start, count) };
        Assert.That(ga, Is.EqualTo(expected), "GenomicAnalyzer");
        Assert.That(rf, Is.EqualTo(expected), "RepeatFinder");
    }

    // ---- Row 15: REP-INV-001 — FindInvertedRepeats vs brute reverse-complement search ----

    private static List<(int i, int j, int arm)> InvertedOracle(string seq, int minArm, int maxLoop, int minLoop)
    {
        var results = new List<(int, int, int)>();
        for (int i = 0; i <= seq.Length - 2 * minArm - minLoop; i++)
        for (int arm = minArm; i + arm <= seq.Length; arm++)
        {
            string leftRc = RevComp(seq.Substring(i, arm));
            int minJ = i + arm + minLoop;
            int maxJ = Math.Min(i + arm + maxLoop, seq.Length - arm);
            for (int j = minJ; j <= maxJ; j++)
                if (seq.Substring(j, arm) == leftRc)
                    results.Add((i, j, arm));
        }
        return results;
    }

    [Test]
    [Category("REP-INV-001")]
    [TestCase("AACCGAGGGTT")]              // arm AACC / loop GAG / arm GGTT (=RC of AACC)
    [TestCase("ACGTACGTAAAACGTACGT")]
    public void InvertedRepeats_MatchesBruteRevCompSearch(string seq)
    {
        var actual = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, maxLoopLength: 50, minLoopLength: 3)
            .Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength)).ToList();
        Assert.That(actual, Is.EqualTo(InvertedOracle(seq.ToUpperInvariant(), 4, 50, 3)));

        // Each result also satisfies the defining hairpin property, checked with the independent RC.
        foreach (var r in RepeatFinder.FindInvertedRepeats(seq, 4, 50, 3))
        {
            Assert.That(r.RightArm, Is.EqualTo(RevComp(r.LeftArm)));
            Assert.That(r.LoopLength, Is.EqualTo(r.RightArmStart - (r.LeftArmStart + r.ArmLength)));
            Assert.That(r.CanFormHairpin, Is.EqualTo(r.LoopLength >= 3));
        }
    }

    // ---- Row 16: REP-DIRECT-001 — suffix-tree FindDirectRepeats vs brute substring scan ----

    private static List<(int i, int j, int len)> DirectOracle(string seq, int minLen, int maxLen, int minSpacing)
    {
        var results = new List<(int, int, int)>();
        for (int len = minLen; len <= maxLen; len++)
        for (int i = 0; i <= seq.Length - len * 2 - minSpacing; i++)
        {
            string repeat = seq.Substring(i, len);
            for (int j = i + len + minSpacing; j + len <= seq.Length; j++)
                if (seq.Substring(j, len) == repeat)
                    results.Add((i, j, len));
        }
        return results;
    }

    [Test]
    [Category("REP-DIRECT-001")]
    [TestCase("ACGTACGTTTT")]
    [TestCase("AAGGAAGGCCAAGG")]
    public void DirectRepeats_MatchesBruteSubstringScan(string seq)
    {
        var actual = RepeatFinder.FindDirectRepeats(seq, minLength: 3, maxLength: 5, minSpacing: 1)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
        Assert.That(actual, Is.EqualTo(DirectOracle(seq.ToUpperInvariant(), 3, 5, 1)));
    }

    // ---- Row 17: REP-PALIN-001 — FindPalindromes vs independent revcomp-equality oracle ----

    private static List<(int pos, string seq, int len)> PalindromeOracle(string seq, int minLen, int maxLen)
    {
        var results = new List<(int, string, int)>();
        for (int len = minLen; len <= maxLen; len += 2)
        for (int i = 0; i + len <= seq.Length; i++)
        {
            string cand = seq.Substring(i, len);
            if (cand == RevComp(cand))
                results.Add((i, cand, len));
        }
        return results;
    }

    [Test]
    [Category("REP-PALIN-001")]
    [TestCase("GAATTC")]            // EcoRI site
    [TestCase("GGGAATTCCCGCGC")]
    [TestCase("ACGTACGT")]
    public void Palindromes_MatchesIndependentRevCompOracle(string seq)
    {
        var actual = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12)
            .Select(p => (p.Position, p.Sequence, p.Length)).ToList();
        Assert.That(actual, Is.EqualTo(PalindromeOracle(seq.ToUpperInvariant(), 4, 12)));
    }
}
