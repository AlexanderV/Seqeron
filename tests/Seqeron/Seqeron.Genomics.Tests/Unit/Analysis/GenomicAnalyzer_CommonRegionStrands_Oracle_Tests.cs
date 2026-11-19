// GENOMIC-COMMON-001 (B09 A2-30 / A2-35) — common regions by mode and strand (MUMmer 3.23 mummer -maxmatch / -mumreference / -mum -b), the
// longest common substring with k mismatches, the longest common substring of N sequences, and the string overloads.
// Oracles (scripts/analysis/generate_common_region_strand_oracle.py):
//   mum   : 4716 (ref, query, L, mode) rows; mummer 3.23 == independent Python brute force on every ACGT row (incl. the -c position convention); rows with N: brute force
//   kmism : 1680 (s1, s2, k) rows, brute force over every (i, j, length)
//   multi : 694 (sequences, support) rows, brute force (Gusfield sec. 7.6 k-common substring)
// Evidence: docs/Evidence/GENOMIC-COMMON-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_CommonRegionStrands_Oracle_Tests
{
    private static List<string[]> Read(string file)
    {
        string name = $"Seqeron.Genomics.Tests.TestData.Analysis.{file}";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) rows.Add(line.Split('\t'));
        return rows;
    }

    private static bool IsAcgt(string s) => s.All(c => "ACGT".Contains(c));

    private static string Fmt(IEnumerable<CommonRegion> regions)
    {
        string s = string.Join(";", regions.OrderBy(r => r.PositionInFirst).ThenBy(r => r.PositionInSecond).ThenBy(r => r.Length)
            .Select(r => $"{r.PositionInFirst}:{r.PositionInSecond}:{r.Length}"));
        return s.Length == 0 ? "-" : s;
    }

    private static string Sorted(string s)
        => s == "-" ? s : string.Join(";", s.Split(';').Select(h => h.Split(':').Select(int.Parse).ToArray()).OrderBy(p => p[0]).ThenBy(p => p[1]).ThenBy(p => p[2])
            .Select(p => string.Join(":", p)));

    private static CommonRegionMode Mode(string m) => m switch
    {
        "maxmatch" => CommonRegionMode.MaximalExactMatches,
        "mumreference" => CommonRegionMode.MaximalUniqueReference,
        "mum" => CommonRegionMode.MaximalUnique,
        _ => throw new ArgumentException(m),
    };

    [Test]
    public void FindCommonRegions_ByMode_BothStrands_MatchesMummerAndBruteForce()
    {
        var rows = Read("common_region_mum_oracle.tsv");
        Assert.That(rows, Has.Count.EqualTo(4716));
        int nonEmpty = 0, withN = 0, reverseRows = 0;
        foreach (var r in rows)
        {
            string reference = r[0], query = r[1];
            int min = int.Parse(r[2]);
            var regions = GenomicAnalyzer.FindCommonRegions(reference, query, min, Mode(r[3]), StrandSelection.Both);
            string forward = Fmt(regions.Where(x => !x.IsReverseComplement)), reverse = Fmt(regions.Where(x => x.IsReverseComplement));
            Assert.That(forward, Is.EqualTo(Sorted(r[4])), $"{reference} {query} L={min} {r[3]} forward");
            Assert.That(reverse, Is.EqualTo(Sorted(r[5])), $"{reference} {query} L={min} {r[3]} reverse");
            if (forward != "-" || reverse != "-") nonEmpty++;
            if (reverse != "-") reverseRows++;
            if (!IsAcgt(reference + query)) withN++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(4000));
        Assert.That(reverseRows, Is.GreaterThan(3000));
        Assert.That(withN, Is.GreaterThan(500));
    }

    [Test]
    public void FindCommonRegions_StrandSelection_AndDnaSequenceOverload()
    {
        foreach (var r in Read("common_region_mum_oracle.tsv").Where(r => IsAcgt(r[0] + r[1])).Take(1500))
        {
            int min = int.Parse(r[2]);
            var mode = Mode(r[3]);
            var a = new DnaSequence(r[0]);
            var b = new DnaSequence(r[1]);
            Assert.That(Fmt(GenomicAnalyzer.FindCommonRegions(a, b, min, mode)), Is.EqualTo(Sorted(r[4])), $"{r[0]} {r[1]} forward-only");
            var onlyReverse = GenomicAnalyzer.FindCommonRegions(a, b, min, mode, StrandSelection.ReverseComplement);
            Assert.That(onlyReverse.All(x => x.IsReverseComplement), Is.True);
            Assert.That(Fmt(onlyReverse), Is.EqualTo(Sorted(r[5])), $"{r[0]} {r[1]} reverse-only");
            foreach (var x in onlyReverse)
            {
                Assert.That(r[0].Substring(x.PositionInFirst, x.Length), Is.EqualTo(x.Sequence));
                Assert.That(DnaSequence.GetReverseComplementString(r[1].Substring(x.PositionInSecond, x.Length)), Is.EqualTo(x.Sequence));
            }
        }
    }

    [Test]
    public void RightMaximalMode_ForwardEqualsTheExistingMethod_AndReverseIsTheSameEngineOnTheReverseComplement()
    {
        var rnd = new Random(31);
        for (int iter = 0; iter < 300; iter++)
        {
            string x = new string(Enumerable.Range(0, rnd.Next(1, 40)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string y = new string(Enumerable.Range(0, rnd.Next(1, 40)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var a = new DnaSequence(x);
            var b = new DnaSequence(y);
            int min = rnd.Next(1, 5);
            var both = GenomicAnalyzer.FindCommonRegions(a, b, min, CommonRegionMode.RightMaximal, StrandSelection.Both);
            Assert.That(both.Where(r => !r.IsReverseComplement).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)),
                Is.EqualTo(GenomicAnalyzer.FindCommonRegions(a, b, min).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond))));
            var rc = new DnaSequence(DnaSequence.GetReverseComplementString(y));
            var expectedReverse = GenomicAnalyzer.FindCommonRegions(a, rc, min)
                .Select(r => (r.Sequence, r.PositionInFirst, y.Length - r.PositionInSecond - r.Length));
            Assert.That(both.Where(r => r.IsReverseComplement).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)), Is.EqualTo(expectedReverse));
        }
    }

    [Test]
    public void FindLongestCommonRegion_Strands_AgreeWithTheMaxmatchOracle()
    {
        int checkedRows = 0;
        foreach (var r in Read("common_region_mum_oracle.tsv").Where(r => r[3] == "maxmatch" && r[2] == "2"))
        {
            int bestF = LongestLength(r[4]), bestR = LongestLength(r[5]);
            var f = GenomicAnalyzer.FindLongestCommonRegion(r[0], r[1]);
            var rv = GenomicAnalyzer.FindLongestCommonRegion(r[0], r[1], StrandSelection.ReverseComplement);
            var both = GenomicAnalyzer.FindLongestCommonRegion(r[0], r[1], StrandSelection.Both);
            // the oracle rows list matches of length >= 2 only: a longest region of length 1 shows up there as nothing
            if (bestF >= 2) Assert.That(f.Length, Is.EqualTo(bestF), $"{r[0]} {r[1]} forward");
            else Assert.That(f.Length, Is.LessThanOrEqualTo(1), $"{r[0]} {r[1]} forward");
            if (bestR >= 2) Assert.That(rv.Length, Is.EqualTo(bestR), $"{r[0]} {r[1]} reverse");
            else Assert.That(rv.Length, Is.LessThanOrEqualTo(1), $"{r[0]} {r[1]} reverse");
            if (Math.Max(bestF, bestR) >= 2)
            {
                Assert.That(both.Length, Is.EqualTo(Math.Max(bestF, bestR)));
                Assert.That(both.IsReverseComplement, Is.EqualTo(bestR > bestF));
            }
            if (rv.Length > 0)
            {
                Assert.That(r[0].Substring(rv.PositionInFirst, rv.Length), Is.EqualTo(rv.Sequence));
                Assert.That(DnaSequence.GetReverseComplementString(r[1].Substring(rv.PositionInSecond, rv.Length)), Is.EqualTo(rv.Sequence));
            }

            checkedRows++;
        }

        Assert.That(checkedRows, Is.GreaterThan(400));

        static int LongestLength(string cell) => cell == "-" ? 0 : cell.Split(';').Max(h => int.Parse(h.Split(':')[2]));
    }

    [Test]
    public void LongestCommonRegionWithMismatches_MatchesBruteForce_AllRows()
    {
        var rows = Read("common_region_kmismatch_oracle.tsv");
        Assert.That(rows, Has.Count.EqualTo(1680));
        int nonEmpty = 0, positiveK = 0;
        foreach (var r in rows)
        {
            string s1 = r[0] == "-" ? "" : r[0], s2 = r[1] == "-" ? "" : r[1];
            int k = int.Parse(r[2]);
            var actual = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(s1, s2, k);
            string text = actual.IsEmpty ? "0" : $"{actual.Length}:{actual.PositionInFirst}:{actual.PositionInSecond}:{actual.Mismatches}";
            Assert.That(text, Is.EqualTo(r[3]), $"{r[0]} {r[1]} k={k}");
            if (!actual.IsEmpty)
            {
                nonEmpty++;
                Assert.That(actual.Sequence1, Is.EqualTo(s1.ToUpperInvariant().Substring(actual.PositionInFirst, actual.Length)));
                Assert.That(actual.Sequence2, Is.EqualTo(s2.ToUpperInvariant().Substring(actual.PositionInSecond, actual.Length)));
                Assert.That(Enumerable.Range(0, actual.Length).Count(i => actual.Sequence1[i] != actual.Sequence2[i] || !"ACGT".Contains(actual.Sequence1[i])), Is.EqualTo(actual.Mismatches));
            }

            if (k > 0) positiveK++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(1300));
        Assert.That(positiveK, Is.GreaterThan(1200));
    }

    [Test]
    public void LongestCommonRegionWithMismatches_ZeroMismatches_IsTheExactMethod()
    {
        var rnd = new Random(9);
        for (int iter = 0; iter < 500; iter++)
        {
            string x = new string(Enumerable.Range(0, rnd.Next(1, 30)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string y = new string(Enumerable.Range(0, rnd.Next(1, 30)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var exact = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(x), new DnaSequence(y));
            var approx = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(new DnaSequence(x), new DnaSequence(y), 0);
            Assert.That((approx.Sequence1, approx.PositionInFirst, approx.PositionInSecond), Is.EqualTo((exact.Sequence, exact.PositionInFirst, exact.PositionInSecond)), $"{x} {y}");
        }
    }

    [Test]
    public void LongestCommonRegionWithMismatches_ReverseStrand()
    {
        // ACGGTTCA vs the reverse complement of its 1-substituted copy
        string s1 = "TTTACGGTTCATTT";
        string s2 = "GGGTGAACCGTAGG"; // revcomp = CCTACGGTTCACCC
        var forward = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(s1, s2, 0, StrandSelection.Forward);
        var reverse = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(s1, s2, 0, StrandSelection.ReverseComplement);
        var both = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(s1, s2, 0, StrandSelection.Both);
        Assert.That(reverse.Length, Is.GreaterThan(forward.Length));
        Assert.That(reverse.IsReverseComplement, Is.True);
        Assert.That(both, Is.EqualTo(reverse));
        Assert.That(DnaSequence.GetReverseComplementString(s2.Substring(reverse.PositionInSecond, reverse.Length)), Is.EqualTo(reverse.Sequence2));
    }

    [Test]
    public void FindLongestCommonRegion_SeveralSequences_MatchesBruteForce_AllRows()
    {
        var rows = Read("common_region_multi_oracle.tsv");
        Assert.That(rows, Has.Count.EqualTo(694));
        int nonEmpty = 0;
        foreach (var r in rows)
        {
            var seqs = r[0].Split(',');
            int q = int.Parse(r[1]);
            var actual = GenomicAnalyzer.FindLongestCommonRegion(seqs, q);
            string text = actual.Count == 0 ? "-" : string.Join(";", actual.Select(a => $"{a.Sequence}@{string.Join(",", a.Positions)}"));
            Assert.That(text, Is.EqualTo(r[2]), $"{r[0]} q={q}");
            if (actual.Count > 0) nonEmpty++;
            if (seqs.All(IsAcgt))
                Assert.That(GenomicAnalyzer.FindLongestCommonRegion(seqs.Select(s => new DnaSequence(s)).ToList(), q).Select(a => a.Sequence),
                    Is.EqualTo(actual.Select(a => a.Sequence)));
        }

        Assert.That(nonEmpty, Is.GreaterThan(600));
    }

    [Test]
    public void FindLongestCommonRegion_SeveralSequences_Validation()
    {
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindLongestCommonRegion(new[] { "ACGT" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindLongestCommonRegion(new[] { "ACGT", "ACGT" }, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindLongestCommonRegion(new[] { "ACGT", "ACGT" }, 3));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestCommonRegion((IReadOnlyList<string>)null!));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestCommonRegion(new string[] { "ACGT", null! }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindLongestCommonRegionWithMismatches("ACGT", "ACGT", -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindCommonRegions("ACGT", "ACGT", 2, (CommonRegionMode)99));
        Assert.That(GenomicAnalyzer.FindCommonRegions(string.Empty, "ACGT", 2), Is.Empty);
        Assert.That(GenomicAnalyzer.FindLongestCommonRegion("NNNN", "NNNN").IsEmpty, Is.True);
    }

    [Test]
    public void StringOverloads_WithN_ReportAnOrdinaryRegionNextToAnNRunUnchanged()
    {
        // MUM / right-maximal regions never contain N and are the same as without the N run
        var withN = GenomicAnalyzer.FindCommonRegions("ACGTACGTNNNNGGATCC", "TTACGTACGTAACCGGATCCTT", 5).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
        Assert.That(withN, Does.Contain(("ACGTACGT", 0, 2)).And.Contain(("GGATCC", 12, 14)));
        Assert.That(withN.Select(r => r.Sequence), Has.None.Contain("N"));
        var maximal = GenomicAnalyzer.FindCommonRegions("ACGTACGTNNNNGGATCC", "TTACGTACGTAACCGGATCCTT", 5, CommonRegionMode.MaximalExactMatches).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
        Assert.That(maximal, Does.Contain(("ACGTACGT", 0, 2)).And.Contain(("GGATCC", 12, 14)));
        Assert.That(maximal.Select(r => r.Sequence), Has.None.Contain("N"));
        var exact = GenomicAnalyzer.FindLongestCommonRegion("ACGTNNNNACGTACGT", "ACGTACGTNACGT");
        Assert.That((exact.Sequence, exact.PositionInFirst, exact.PositionInSecond), Is.EqualTo(("ACGTACGT", 8, 0)));
    }
}
