// 08_DIFFERENTIAL_TESTING rows 96, 149, 232, 234. Independent oracles: a stack dot-bracket parse, a
// naive per-window GC profile, a naive per-window Shannon entropy profile, and a manual SBS-96
// trinucleotide context classification + channel enumeration.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class ProfilesSbsDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 149: RNA-DOTBRACKET-001 — ParseDotBracket vs stack parse ----

    [Test]
    [Category("RNA-DOTBRACKET-001")]
    [TestCase("(())")]
    [TestCase("(.)")]
    [TestCase("((..))")]
    [TestCase("(()())")]
    [TestCase("....")]
    public void ParseDotBracket_MatchesStackParse(string db)
    {
        var stack = new Stack<int>();
        var expected = new List<(int, int)>();
        for (int i = 0; i < db.Length; i++)
        {
            if (db[i] == '(') stack.Push(i);
            else if (db[i] == ')' && stack.Count > 0) expected.Add((stack.Pop(), i));
        }
        Assert.That(RnaSecondaryStructure.ParseDotBracket(db).ToList(), Is.EqualTo(expected));
    }

    // ---- Row 234: SEQ-GC-PROFILE-001 — GC profile vs naive per-window count ----

    [Test]
    [Category("SEQ-GC-PROFILE-001")]
    [TestCase("ACGTACGTAC", 4, 1)]
    [TestCase("GGGGCCCCAAAA", 3, 2)]
    public void GcContentProfile_MatchesNaivePerWindow(string seq, int win, int step)
    {
        var s = seq.ToUpperInvariant();
        var expected = new List<double>();
        for (int i = 0; i + win <= s.Length; i += step)
        {
            var w = s.Substring(i, win);
            int gc = w.Count(c => c == 'G' || c == 'C');
            int total = w.Count(c => "ATGCU".Contains(c));
            expected.Add(total > 0 ? (double)gc / total * 100 : 0);
        }
        Assert.That(SequenceStatistics.CalculateGcContentProfile(seq, win, step).ToList(),
            Is.EqualTo(expected).AsCollection.Within(Tol));
    }

    // ---- Row 232: SEQ-ENTROPY-PROFILE-001 — entropy profile vs naive per-window Shannon ----

    private static double ShannonOracle(string s)
    {
        var letters = s.ToUpperInvariant().Where(char.IsLetter).ToList();
        if (letters.Count == 0) return 0;
        return -letters.GroupBy(c => c).Select(g => (double)g.Count() / letters.Count).Sum(p => p * Math.Log2(p));
    }

    [Test]
    [Category("SEQ-ENTROPY-PROFILE-001")]
    [TestCase("ACGTACGTAC", 4, 1)]
    [TestCase("AAAACCCCGG", 5, 2)]
    public void EntropyProfile_MatchesNaivePerWindow(string seq, int win, int step)
    {
        var s = seq;
        var expected = new List<double>();
        for (int i = 0; i + win <= s.Length; i += step)
            expected.Add(ShannonOracle(s.Substring(i, win)));
        Assert.That(SequenceStatistics.CalculateEntropyProfile(seq, win, step).ToList(),
            Is.EqualTo(expected).AsCollection.Within(Tol));
    }

    // ---- Row 96: ONCO-SIG-001 — SBS-96 trinucleotide context classification ----

    private static char Comp(char c) => c switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => c };

    private static string SbsOracle(char five, char refB, char alt, char three)
    {
        if (refB is 'A' or 'G')
        {
            char f = Comp(three), t = Comp(five);
            return $"{f}[{Comp(refB)}>{Comp(alt)}]{t}";
        }
        return $"{five}[{refB}>{alt}]{three}";
    }

    [Test]
    [Category("ONCO-SIG-001")]
    [TestCase('A', 'C', 'A', 'G')]   // pyrimidine ref kept
    [TestCase('A', 'G', 'A', 'T')]   // purine ref folded -> A[C>T]T
    [TestCase('T', 'T', 'C', 'A')]
    [TestCase('C', 'A', 'G', 'T')]
    public void ClassifySbsContext_MatchesManualFold(char five, char refB, char alt, char three)
    {
        Assert.That(OncologyAnalyzer.ClassifySbsContext(five, refB, alt, three),
            Is.EqualTo(SbsOracle(five, refB, alt, three)));
    }

    [Test]
    [Category("ONCO-SIG-001")]
    public void EnumerateSbs96Channels_Matches6x4x4()
    {
        var subs = new[] { ("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G") };
        var expected = new List<string>();
        foreach (var (r, a) in subs)
            foreach (char five in "ACGT")
                foreach (char three in "ACGT")
                    expected.Add($"{five}[{r}>{a}]{three}");
        Assert.That(OncologyAnalyzer.EnumerateSbs96Channels(), Is.EqualTo(expected));
    }
}
