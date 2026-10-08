// PAT-APPROX-002 — Edit distance: traceback (CIGAR), Myers bit-parallel engine, Damerau variants.
// TestSpec: tests/TestSpecs/PAT-APPROX-002.md
// Review 2026-09 B05 follow-up (docs/Validation/review-2026-09/B05.md).
// References (values computed 2026-09-30):
//   edlib 1.3 (edlib.align(q, t, mode='NW'|'HW', task='path'|'locations'); Šošić & Šikić 2017),
//   rapidfuzz 3.14 (Levenshtein, OSA, DamerauLevenshtein), jellyfish (damerau_levenshtein_distance).
//   Myers (1999) J. ACM 46(3):395; Hyyrö (2003) global form; Lowrance & Wagner (1975) J. ACM 22(2):177.
namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class ApproximateMatcher_EditAlignment_Tests
{
    #region GetEditAlignment — traceback / CIGAR

    // Unique optimal paths: our CIGAR is identical to edlib NW task='path'.
    [TestCase("ACGTACGT", "ACGACGT", 1, "3=1I4=", "3M1I4M", TestName = "GetEditAlignment_ACGTACGT_ACGACGT_EqualsEdlib")]
    [TestCase("survey", "surgery", 2, "3=1X1=1D1=", "5M1D1M", TestName = "GetEditAlignment_SurveySurgery_EqualsEdlib")]
    [TestCase("kitten", "sitting", 3, "1X3=1X1=1D", "6M1D", TestName = "GetEditAlignment_KittenSitting_EqualsEdlib")]
    [TestCase("ACGTTGCA", "ACGTAGCA", 1, "4=1X3=", "8M", TestName = "GetEditAlignment_SingleSubstitution_EqualsEdlib")]
    public void GetEditAlignment_UniqueOptimalPath_EqualsEdlibCigar(
        string query, string target, int distance, string cigar, string standardCigar)
    {
        var a = ApproximateMatcher.GetEditAlignment(query, target);

        Assert.Multiple(() =>
        {
            Assert.That(a.Distance, Is.EqualTo(distance));
            Assert.That(a.Cigar, Is.EqualTo(cigar));
            Assert.That(a.StandardCigar, Is.EqualTo(standardCigar));
            Assert.That(a.Distance, Is.EqualTo(ApproximateMatcher.EditDistance(query, target)));
        });
    }

    // Co-optimal paths: edlib (I, then D, then diagonal on the way back) returns a different path
    // with the same cost; the library's diagonal-first tie-break is deterministic.
    [TestCase("ACGT", "ACGGT", 1, "2=1D2=", "3=1D1=", TestName = "GetEditAlignment_ACGT_ACGGT_DiagonalFirstTieBreak")]
    [TestCase("GATTACA", "GCATGCT", 4, "1=2X1=1X1=1X", "1=1D2=1X1I1=1X", TestName = "GetEditAlignment_GATTACA_GCATGCT_DiagonalFirstTieBreak")]
    [TestCase("AAAA", "AAAAAA", 2, "2D4=", "4=2D", TestName = "GetEditAlignment_Homopolymer_DiagonalFirstTieBreak")]
    public void GetEditAlignment_CoOptimalPaths_UsesDiagonalFirstTieBreak(
        string query, string target, int distance, string expectedCigar, string edlibCigar)
    {
        var a = ApproximateMatcher.GetEditAlignment(query, target);

        Assert.Multiple(() =>
        {
            Assert.That(a.Distance, Is.EqualTo(distance), "edlib editDistance");
            Assert.That(a.Cigar, Is.EqualTo(expectedCigar));
            Assert.That(a.Cigar, Is.Not.EqualTo(edlibCigar), "edlib chose another co-optimal path");
            Assert.That(ReplayCost(edlibCigar, query, target), Is.EqualTo(distance), "edlib path has the same cost");
            Assert.That(ReplayCost(a.Cigar, query, target), Is.EqualTo(distance));
        });
    }

    [Test]
    public void GetEditAlignment_AlignedStringsAndSubstitutions()
    {
        var a = ApproximateMatcher.GetEditAlignment("survey", "surgery");

        Assert.Multiple(() =>
        {
            Assert.That(a.Operations, Is.EqualTo("===X=D="));
            Assert.That(a.AlignedQuery, Is.EqualTo("surve-y"));
            Assert.That(a.AlignedTarget, Is.EqualTo("surgery"));
            Assert.That(a.SubstitutionPositions, Is.EqualTo(new[] { 3 }));
            Assert.That(a.HasIndels, Is.True);
        });
    }

    [Test]
    public void GetEditAlignment_EmptyInputs()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.GetEditAlignment("", "ACG").Cigar, Is.EqualTo("3D"));
            Assert.That(ApproximateMatcher.GetEditAlignment("ACG", "").Cigar, Is.EqualTo("3I"));
            var empty = ApproximateMatcher.GetEditAlignment("", "");
            Assert.That(empty.Cigar, Is.EqualTo(""));
            Assert.That(empty.Distance, Is.EqualTo(0));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetEditAlignment(null!, "A"));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetEditAlignment("A", null!));
        });
    }

    #endregion

    #region FindWithEdits — per-hit alignment and MismatchPositions

    [Test]
    [Description("Navarro 2001 survey/surgery, k=2: per-window CIGARs equal edlib NW on (SURVEY, window)")]
    public void FindWithEdits_SurveySurgery_AlignmentsEqualEdlib()
    {
        var hits = ApproximateMatcher.FindWithEdits("surgery", "survey", 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(hits.Select(h => h.Alignment!.Cigar),
                Is.EqualTo(new[] { "3=1X1=1I", "3=1X1=1X", "3=1X1=1D1=" }));
            Assert.That(hits.Select(h => h.MismatchType),
                Is.EqualTo(new[] { MismatchType.Edit, MismatchType.Substitution, MismatchType.Edit }));
            Assert.That(hits[0].MismatchPositions, Is.EqualTo(new[] { 3 }));
            Assert.That(hits[1].MismatchPositions, Is.EqualTo(new[] { 3, 5 }), "Hamming mismatch indices SURVEY/SURGER");
            Assert.That(hits[2].MismatchPositions, Is.EqualTo(new[] { 3 }));
        });
    }

    [Test]
    [Description("TTAC in GATTACAGATTTACA, k=1 (windows reported with their optimal CIGAR)")]
    public void FindWithEdits_TtacInGattaca_Cigars()
    {
        var hits = ApproximateMatcher.FindWithEdits("GATTACAGATTTACA", "TTAC", 1)
            .Select(h => (h.Position, h.MatchedSequence.Length, h.Distance, h.Alignment!.Cigar)).ToList();

        Assert.That(hits, Is.EqualTo(new[]
        {
            (1, 5, 1, "1D4="), (2, 3, 1, "3=1I"), (2, 4, 0, "4="), (2, 5, 1, "4=1D"), (3, 3, 1, "1I3="),
            (9, 5, 1, "1D4="), (10, 3, 1, "3=1I"), (10, 4, 0, "4="), (10, 5, 1, "4=1D"), (11, 3, 1, "1I3="),
        }));
    }

    [Test]
    public void FindWithEdits_HugeMaxEdits_DoesNotOverflow()
    {
        var hits = ApproximateMatcher.FindWithEdits("ACG", "AC", int.MaxValue).ToList();

        // Every non-empty window: A, AC, ACG, C, CG, G.
        Assert.That(hits.Select(h => h.MatchedSequence), Is.EqualTo(new[] { "A", "AC", "ACG", "C", "CG", "G" }));
        Assert.That(hits.Select(h => h.Distance), Is.EqualTo(new[] { 1, 0, 1, 1, 2, 2 }));
    }

    [Test]
    public void FindWithMismatches_HammingResults_HaveNoAlignment()
    {
        Assert.That(ApproximateMatcher.FindWithMismatches("ACGT", "ACG", 1).All(r => r.Alignment is null), Is.True);
    }

    #endregion

    #region Myers bit-parallel engine

    private static string Repeat(string unit, int times) => string.Concat(Enumerable.Repeat(unit, times));

    [Test]
    [Description("Multi-word (m > 64) global distances equal rapidfuzz/edlib")]
    public void EditDistance_MultiBlock_EqualsRapidfuzz()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.EditDistance(Repeat("ACGT", 40), Repeat("ACGA", 40)), Is.EqualTo(40));
            Assert.That(ApproximateMatcher.EditDistance(Repeat("ACGT", 40), Repeat("ACGT", 50)), Is.EqualTo(40));
            Assert.That(ApproximateMatcher.EditDistance(Repeat("ACGT", 50), Repeat("ACGT", 33)[..130]), Is.EqualTo(70));
            Assert.That(ApproximateMatcher.EditDistance(new string('A', 65), new string('A', 64) + "C"), Is.EqualTo(1));
            Assert.That(ApproximateMatcher.EditDistance(new string('A', 200), new string('C', 130)), Is.EqualTo(200));
            Assert.That(ApproximateMatcher.EditDistance(Repeat("GATTACA", 20), Repeat("GATACA", 20)), Is.EqualTo(20));
        });
    }

    [Test]
    [Description("Exhaustive: Myers == Wagner–Fischer DP for every pair of strings over {A,C} of length ≤ 6")]
    public void EditDistance_Myers_EqualsDp_Exhaustive()
    {
        var all = new List<string> { "" };
        for (int len = 1; len <= 6; len++)
            all.AddRange(Enumerable.Range(0, 1 << len)
                .Select(bits => new string(Enumerable.Range(0, len).Select(i => ((bits >> i) & 1) == 0 ? 'A' : 'C').ToArray())));

        int mismatches = 0;
        foreach (var a in all)
            foreach (var b in all)
                if (ApproximateMatcher.EditDistance(a, b) != ApproximateMatcher.EditDistanceDp(a, b))
                    mismatches++;

        Assert.That(mismatches, Is.EqualTo(0));
        Assert.That(all.Count, Is.EqualTo(127));
    }

    [Test]
    public void EditDistance_NonAsciiSymbols_EqualDp()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.EditDistance("αβγ", "αγ"), Is.EqualTo(1));
            Assert.That(ApproximateMatcher.EditDistance("ACαβ", "βαCA"), Is.EqualTo(ApproximateMatcher.EditDistanceDp("ACαβ", "βαCA")));
        });
    }

    [Test]
    [Description("Sellers end positions (Myers engine) == DP reference, incl. m > 64")]
    public void FindEditEndPositions_Myers_EqualsDp()
    {
        var rng = new Random(20260930);
        for (int trial = 0; trial < 200; trial++)
        {
            int m = trial % 4 == 0 ? rng.Next(60, 140) : rng.Next(1, 20);
            string pattern = new(Enumerable.Range(0, m).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
            string text = new(Enumerable.Range(0, rng.Next(1, 300)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
            int k = rng.Next(0, m + 1);

            Assert.That(ApproximateMatcher.FindEditEndPositions(text, pattern, k).ToList(),
                Is.EqualTo(ApproximateMatcher.FindEditEndPositionsDp(text, pattern, k).ToList()),
                $"trial {trial}");
        }
    }

    #endregion

    #region Damerau variants (OSA, true Damerau–Levenshtein)

    // rapidfuzz.distance.OSA / DamerauLevenshtein / Levenshtein; jellyfish agrees on DL.
    [TestCase("CA", "ABC", 3, 2, 3)]
    [TestCase("ca", "abc", 3, 2, 3)]
    [TestCase("ab", "ba", 1, 1, 2)]
    [TestCase("a cat", "an act", 2, 2, 3)]
    [TestCase("abcdef", "abcfad", 3, 3, 3)]
    [TestCase("ACGT", "CAGT", 1, 1, 2)]
    [TestCase("ACGT", "TGCA", 3, 3, 4)]
    [TestCase("AGCT", "GACT", 1, 1, 2)]
    [TestCase("kitten", "sitting", 3, 3, 3)]
    [TestCase("", "abc", 3, 3, 3)]
    [TestCase("abc", "", 3, 3, 3)]
    public void DamerauVariants_EqualRapidfuzz(string a, string b, int osa, int dl, int lev)
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b), Is.EqualTo(osa), "OSA");
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b), Is.EqualTo(dl), "DL");
            Assert.That(ApproximateMatcher.EditDistance(a, b), Is.EqualTo(lev), "Levenshtein");
        });
    }

    [Test]
    public void DamerauVariants_NullInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.OptimalStringAlignmentDistance(null!, "A"));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.OptimalStringAlignmentDistance("A", null!));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.DamerauLevenshteinDistance(null!, "A"));
            Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.DamerauLevenshteinDistance("A", null!));
        });
    }

    #endregion

    /// <summary>Replays an extended CIGAR over (query, target); returns its unit cost.</summary>
    internal static int ReplayCost(string cigar, string query, string target)
    {
        int q = 0, t = 0, cost = 0;
        foreach (System.Text.RegularExpressions.Match op in System.Text.RegularExpressions.Regex.Matches(cigar, @"(\d+)([=XID])"))
        {
            int n = int.Parse(op.Groups[1].Value);
            for (int i = 0; i < n; i++)
            {
                switch (op.Groups[2].Value[0])
                {
                    case '=': Assert.That(query[q], Is.EqualTo(target[t])); q++; t++; break;
                    case 'X': Assert.That(query[q], Is.Not.EqualTo(target[t])); q++; t++; cost++; break;
                    case 'I': q++; cost++; break;
                    default: t++; cost++; break;
                }
            }
        }
        Assert.That((q, t), Is.EqualTo((query.Length, target.Length)), "CIGAR must consume both strings");
        return cost;
    }
}
