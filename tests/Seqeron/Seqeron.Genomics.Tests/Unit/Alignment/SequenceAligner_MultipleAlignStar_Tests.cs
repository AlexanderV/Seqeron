namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// Center-star MSA (Gusfield 1993, Bull Math Biol 55:141-154) — exact pairwise stage.
/// Test Unit: ALIGN-MULTI-001 (B09 review F4/F5).
///
/// Oracles (independent of the implementation):
/// - Biopython 1.88 <c>PairwiseAligner(mode='global', match=1, mismatch=-1, open_gap_score=-1,
///   extend_gap_score=-1)</c> (= the linear model with GapExtend = -1) for pair optima;
/// - an independent Python implementation of the center-star method (own Needleman-Wunsch with
///   diagonal &gt; up &gt; left traceback, Gusfield center = argmax of total pair score, gap-slot merge)
///   for the exact rows of the 3/4/5-sequence cases; all SP scores recomputed in Python with the
///   column definition (match +1, mismatch -1, residue-gap -1, gap-gap 0).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignStar_Tests
{
    private static List<DnaSequence> Dna(params string[] s) => s.Select(x => new DnaSequence(x)).ToList();

    private static void AssertValidMsa(MultipleAlignmentResult r, IReadOnlyList<string> input)
    {
        Assert.That(r.AlignedSequences, Has.Length.EqualTo(input.Count));
        int len = r.AlignedSequences[0].Length;
        for (int i = 0; i < input.Count; i++)
        {
            Assert.That(r.AlignedSequences[i], Has.Length.EqualTo(len));
            Assert.That(r.AlignedSequences[i].Replace("-", ""), Is.EqualTo(input[i]));
        }

        for (int c = 0; c < len; c++)
            Assert.That(r.AlignedSequences.Any(row => row[c] != '-'), Is.True, $"column {c} is all gaps");
    }

    private static int Sp(string[] rows, ScoringMatrix s)
    {
        int total = 0;
        for (int c = 0; c < rows[0].Length; c++)
            for (int i = 0; i < rows.Length; i++)
                for (int j = i + 1; j < rows.Length; j++)
                {
                    char a = rows[i][c], b = rows[j][c];
                    if (a == '-' && b == '-') continue;
                    total += a == '-' || b == '-' ? s.GapExtend : a == b ? s.Match : s.Mismatch;
                }
        return total;
    }

    // Biopython optimum (linear gap -1, match 1, mismatch -1) for pairs on which the previous
    // suffix-tree-anchored pairwise stage was suboptimal (its MSA scored 9, 22 and -35).
    [TestCase("GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT", "AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT", 24)]
    [TestCase("TAGACGCATCAAACAGCCTAGGCCGTACTGCGATTCAGCTGTCAACAGGTAAAACAGA", "TAGACGCCTAGTCAACAGCCTACGGCCCTATCAGCTGTATACATGTTAGTCCACCGTAACAGT", 27)]
    [TestCase("CACTTGGGGACAGTTAGCTCCGCCCCCGTTAGTTTGATGCGGGTCC", "ACGGCGCGTAGGAAGTATAGTAATTCTTAGATGGGGACACGATACTA", -3)]
    public void MultipleAlign_TwoSequences_IsOptimalPairwiseAlignment(string a, string b, int biopythonOptimum)
    {
        var r = SequenceAligner.MultipleAlign(Dna(a, b));

        AssertValidMsa(r, new[] { a, b });
        Assert.That(r.TotalScore, Is.EqualTo(biopythonOptimum));
        Assert.That(SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b)).Score, Is.EqualTo(biopythonOptimum));
    }

    [Test]
    public void MultipleAlign_ThreeSequences_MatchesIndependentCenterStar()
    {
        // Python reference: center = index 0 (total pair score highest), rows and SP = -3.
        var r = SequenceAligner.MultipleAlign(Dna("ACACG", "AGAGC", "AATGGC"));

        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "-A-CACG-", "-A-GA-GC", "AATG--GC" }));
        Assert.That(r.TotalScore, Is.EqualTo(-3));
    }

    [Test]
    public void MultipleAlign_FourSequences_ParallelPath_MatchesIndependentCenterStar()
    {
        var r = SequenceAligner.MultipleAlign(Dna("TTTGA", "AGGTCA", "AGTGCGA", "ACTGGCA"));

        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "---T-TTGA", "AGGT--C-A", "A-GT-GCGA", "A-CTGGC-A" }));
        Assert.That(r.TotalScore, Is.EqualTo(-5));
    }

    [Test]
    public void MultipleAlign_FiveSequences_MatchesIndependentCenterStar()
    {
        var r = SequenceAligner.MultipleAlign(Dna("CAT", "TAT", "CGAT", "CTAT", "CTA"));

        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "C-AT", "-TAT", "CGAT", "CTAT", "CTA-" }));
        Assert.That(r.TotalScore, Is.EqualTo(10));
    }

    [Test]
    public void MultipleAlign_CenterIsTheSequenceWithBestTotalPairScore_NotTheFirst()
    {
        // Pair scores (SimpleDna, Biopython/python NW): S0-S1 = 0, S0-S2 = -2, S1-S2 = 4 -> totals
        // S0 = -2, S1 = 4, S2 = 2; the center is S1. All three are ungapped equal-length rows.
        var input = new[] { "GGGGAC", "ACGTAC", "ACGTAA" };
        var r = SequenceAligner.MultipleAlign(Dna(input));

        AssertValidMsa(r, input);
        Assert.That(r.AlignedSequences, Is.EqualTo(input));
        Assert.That(r.TotalScore, Is.EqualTo(0 + (-2) + 4));
    }

    [Test]
    public void MultipleAlign_IsDeterministic_AndTotalScoreIsColumnSumOfPairs()
    {
        var input = Dna("ACGTTGCA", "ACGTGCA", "ACTTGCA", "AGTTGGCA", "ACGTTGCC");
        var r1 = SequenceAligner.MultipleAlign(input);
        var r2 = SequenceAligner.MultipleAlign(input);

        Assert.That(r1.AlignedSequences, Is.EqualTo(r2.AlignedSequences));
        Assert.That(r1.TotalScore, Is.EqualTo(Sp(r1.AlignedSequences, SequenceAligner.SimpleDna)));
    }

    [Test]
    public void MultipleAlign_RandomSets_AreValidMsas_AndTwoSequenceSetsAreOptimal()
    {
        var rng = new Random(20261104);
        for (int trial = 0; trial < 150; trial++)
        {
            int k = rng.Next(2, 7);
            var input = Enumerable.Range(0, k).Select(_ => RandomDna(rng, rng.Next(1, 30))).ToArray();
            var r = SequenceAligner.MultipleAlign(Dna(input));

            AssertValidMsa(r, input);
            Assert.That(r.TotalScore, Is.EqualTo(Sp(r.AlignedSequences, SequenceAligner.SimpleDna)));

            if (k == 2)
                Assert.That(r.TotalScore,
                    Is.EqualTo(SequenceAligner.GlobalAlign(new DnaSequence(input[0]), new DnaSequence(input[1])).Score));
        }
    }

    [Test]
    public void MultipleAlign_LargePair_UsesLinearSpaceRecursion_AndStaysOptimal()
    {
        // 2400 x 2400 = 5.8 M cells > the 4 M full-matrix limit -> Hirschberg path; the MSA score
        // must still equal the canonical Needleman-Wunsch optimum.
        var rng = new Random(5150);
        string a = RandomDna(rng, 2400);
        string b = Mutate(rng, a, 0.1);

        var r = SequenceAligner.MultipleAlign(Dna(a, b));

        AssertValidMsa(r, new[] { a, b });
        Assert.That(r.TotalScore, Is.EqualTo(SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b)).Score));
    }

    [Test]
    public void MultipleAlign_EmptyMemberInsideSet_IsAllGapRow()
    {
        var r = SequenceAligner.MultipleAlign(Dna("ACGT", "", "ACGA"));

        AssertValidMsa(r, new[] { "ACGT", "", "ACGA" });
        Assert.That(r.AlignedSequences[1], Is.EqualTo("----"));
    }

    // ----- MultipleAlignAnchored (approximate, suffix-tree anchored variant) ------------------------

    [TestCase("GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT", "AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT", 24)]
    [TestCase("TAGACGCATCAAACAGCCTAGGCCGTACTGCGATTCAGCTGTCAACAGGTAAAACAGA", "TAGACGCCTAGTCAACAGCCTACGGCCCTATCAGCTGTATACATGTTAGTCCACCGTAACAGT", 27)]
    public void MultipleAlignAnchored_TwoSequences_IsValid_AndNeverBeatsTheOptimum(string a, string b, int optimum)
    {
        var r = SequenceAligner.MultipleAlignAnchored(Dna(a, b));

        AssertValidMsa(r, new[] { a, b });
        Assert.That(r.TotalScore, Is.EqualTo(Sp(r.AlignedSequences, SequenceAligner.SimpleDna)));
        Assert.That(r.TotalScore, Is.LessThanOrEqualTo(optimum));
    }

    [Test]
    public void MultipleAlignAnchored_IdenticalSequences_AreGapFree()
    {
        var r = SequenceAligner.MultipleAlignAnchored(Dna("ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC"));

        Assert.That(r.AlignedSequences, Is.All.EqualTo("ACGTACGTAC"));
        Assert.That(r.TotalScore, Is.EqualTo(30));
    }

    [Test]
    public void MultipleAlignAnchored_RelatedSets_AreValidMsas()
    {
        var rng = new Random(77);
        for (int trial = 0; trial < 30; trial++)
        {
            string root = RandomDna(rng, rng.Next(40, 160));
            int k = rng.Next(2, 6);
            var input = Enumerable.Range(0, k).Select(_ => Mutate(rng, root, 0.08)).Select(s => s.Length == 0 ? "A" : s).ToArray();
            var r = SequenceAligner.MultipleAlignAnchored(Dna(input));

            AssertValidMsa(r, input);
            Assert.That(r.TotalScore, Is.EqualTo(Sp(r.AlignedSequences, SequenceAligner.SimpleDna)));
        }
    }

    [Test]
    public void MultipleAlignAnchored_NullAndTrivialInputs()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignAnchored(null!));
        Assert.That(SequenceAligner.MultipleAlignAnchored(Array.Empty<DnaSequence>()), Is.EqualTo(MultipleAlignmentResult.Empty));
        Assert.That(SequenceAligner.MultipleAlignAnchored(Dna("ACG")).AlignedSequences, Is.EqualTo(new[] { "ACG" }));
    }

    private static string RandomDna(Random rng, int n)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[rng.Next(4)];
        return new string(c);
    }

    private static string Mutate(Random rng, string s, double p)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char ch in s)
        {
            double x = rng.NextDouble();
            if (x < p) sb.Append("ACGT"[rng.Next(4)]);
            else if (x < p * 1.5) { }
            else if (x < p * 2) { sb.Append(ch); sb.Append("ACGT"[rng.Next(4)]); }
            else sb.Append(ch);
        }
        return sb.ToString();
    }
}
