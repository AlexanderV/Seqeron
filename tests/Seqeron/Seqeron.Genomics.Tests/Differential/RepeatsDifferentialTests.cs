// 08_DIFFERENTIAL_TESTING rows 13-17 (Repeats). Production finders vs INDEPENDENT oracles:
// regex homopolymer detection, the second in-project tandem implementation (DUAL), a brute
// reverse-complement search, a brute substring-equality scan, and an independent reverse-complement
// palindrome test. The reverse complement used by the oracles is computed from a literal IUPAC table,
// not from DnaSequence.GetReverseComplementString.

using System.Text.RegularExpressions;

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

    // Brute-force oracle for the documented semantics (EMBOSS palindrome, -nummismatches 0 -overlap Y):
    // enumerate EVERY exact stem (i, j, arm) with RightArm = RC(LeftArm) and the loop in bounds, then keep only
    // the stems that are not contained, in both arms, in another stem (palindrome.c palindrome_AInB).
    // The same definition, as an independent Python script, matched the EMBOSS 6.6.0 palindrome binary on
    // 1000 random ACGT sequences (minLoop = 0) — docs/Evidence/REP-INV-001-Evidence.md.
    private static List<(int i, int j, int arm)> InvertedOracle(string seq, int minArm, int maxLoop, int minLoop)
    {
        var all = new List<(int i, int j, int arm)>();
        for (int i = 0; i <= seq.Length - 2 * minArm - minLoop; i++)
        for (int arm = minArm; i + 2 * arm + minLoop <= seq.Length; arm++)
        {
            string leftRc = RevComp(seq.Substring(i, arm));
            int minJ = i + arm + minLoop;
            int maxJ = Math.Min(i + arm + maxLoop, seq.Length - arm);
            for (int j = minJ; j <= maxJ; j++)
                if (seq.Substring(j, arm) == leftRc)
                    all.Add((i, j, arm));
        }

        return all.Where(a => !all.Any(b => b != a
                && b.i <= a.i && a.i + a.arm <= b.i + b.arm
                && b.j <= a.j && a.j + a.arm <= b.j + b.arm))
            .OrderBy(a => a.i).ThenBy(a => a.j).ToList();
    }

    [Test]
    [Category("REP-INV-001")]
    [TestCase("AACCGAGGGTT")]              // arm AACC / loop GAG / arm GGTT (=RC of AACC)
    [TestCase("ACGTACGTAAAACGTACGT")]
    [TestCase("GGGGGGAAACCCCCC")]          // slipped re-pairings inside the 6-bp stem are dropped
    [TestCase("GAATTCAAAAGAATTCTTTTGAATTC")]
    [TestCase("ATATATATATGCATATATATAT")]
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

    // ---- Row 16: REP-DIRECT-001 — FindDirectRepeats vs brute-force maximal-pair oracle ----
    // Oracle = the definition (Gusfield 1997 §7.12; MUMmer repeat-match -f exhaustive mode -E):
    // i < j, left-maximal (i == 0 or S[i-1] != S[j-1]), L = full common-prefix length over A/C/G/T,
    // minLen <= L <= maxLen, Spacing = j - i - L >= minSpacing. Sorted by (i, j).

    private static List<(int i, int j, int len)> DirectOracle(string seq, int minLen, int maxLen, int minSpacing)
    {
        static bool Acgt(char c) => c is 'A' or 'C' or 'G' or 'T';
        var results = new List<(int, int, int)>();
        for (int i = 0; i < seq.Length; i++)
        for (int j = i + 1; j < seq.Length; j++)
        {
            if (i > 0 && seq[i - 1] == seq[j - 1] && Acgt(seq[i - 1])) continue;
            int len = 0;
            while (j + len < seq.Length && seq[i + len] == seq[j + len] && Acgt(seq[i + len])) len++;
            if (len >= minLen && len <= maxLen && j - i - len >= minSpacing)
                results.Add((i, j, len));
        }
        return results;
    }

    [Test]
    [Category("REP-DIRECT-001")]
    [TestCase("ACGTACGTTTT", 1)]
    [TestCase("AAGGAAGGCCAAGG", 1)]
    [TestCase("ACGTACGTTTTTTTTTACGTACGT", 1)]
    [TestCase("acgtNacgtRRacgtNacgt", 0)]
    [TestCase("ACGTACGTACGTAAACGTACG", -100)]
    public void DirectRepeats_MatchesBruteForceMaximalPairOracle(string seq, int minSpacing)
    {
        var actual = RepeatFinder.FindDirectRepeats(seq, minLength: 3, maxLength: 5, minSpacing: minSpacing)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
        Assert.That(actual, Is.EqualTo(DirectOracle(seq.ToUpperInvariant(), 3, 5, minSpacing)));
    }

    [Test]
    [Category("REP-DIRECT-001")]
    public void DirectRepeats_RandomSequences_MatchBruteForceMaximalPairOracle()
    {
        var rng = new Random(20260930);
        for (int t = 0; t < 300; t++)
        {
            string alphabet = t % 3 == 0 ? "ACGTN" : t % 3 == 1 ? "AC" : "ACGT";
            var seq = new string(Enumerable.Range(0, rng.Next(0, 80)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            int minLen = rng.Next(2, 7), maxLen = minLen + rng.Next(0, 10), minSpacing = rng.Next(-5, 4);
            var actual = RepeatFinder.FindDirectRepeats(seq, minLen, maxLen, minSpacing)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
            Assert.That(actual, Is.EqualTo(DirectOracle(seq, minLen, maxLen, minSpacing)), $"{seq} {minLen} {maxLen} {minSpacing}");
        }
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
