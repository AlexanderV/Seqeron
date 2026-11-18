using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class MultipleAlignTests
{
    [Test]
    public void MultipleAlign_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACGT" }));
        // Invalid DNA in any sequence is rejected.
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACXT" }));
        // Empty set yields an empty MSA (no throw).
        var empty = AlignmentTools.MultipleAlign(Array.Empty<string>());
        Assert.That(empty.AlignedSequences, Is.Empty);
    }

    [Test]
    public void MultipleAlign_Binding_InvokesSuccessfully()
    {
        // Three identical 8-mers align without gaps; consensus equals the input.
        // Sum-of-pairs score = C(3,2)=3 pairs * 8 matched columns * (+1) = 24.
        var r = AlignmentTools.MultipleAlign(new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT" });
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Has.Length.EqualTo(3));
            Assert.That(r.AlignedSequences, Is.All.EqualTo("ACGTACGT"));
            Assert.That(r.Consensus, Is.EqualTo("ACGTACGT"));
            Assert.That(r.TotalScore, Is.EqualTo(24));
        });
    }

    // B09 F4: multiple_align is the exact center-star: two sequences give an optimal pairwise alignment.
    // Biopython 1.88 optimum (match 1, mismatch -1, linear gap -1) of this pair is 24; the previous
    // suffix-tree-anchored pairwise stage returned 9.
    [Test]
    public void MultipleAlign_TwoSequences_IsOptimalPairwiseAlignment()
    {
        var r = AlignmentTools.MultipleAlign(new[]
        {
            "GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT",
            "AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT"
        });

        Assert.That(r.TotalScore, Is.EqualTo(24));
    }

    [Test]
    public void MultipleAlignVariants_Schema_ValidatesCorrectly()
    {
        var bad = new[] { "ACGT", "ACXT" };
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignAnchored(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MultipleAlignIterative(new[] { "ACGT", "ACGT" }, maxIterations: 0));

        Assert.That(AlignmentTools.MultipleAlignAnchored(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignProgressive(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignIterative(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignConsistency(Array.Empty<string>()).AlignedSequences, Is.Empty);
    }

    // Hand-derived (progressive unit tests): ACGT/ACGT/AGT -> ACGT, ACGT, A-GT; SP = 3 - 1 + 3 + 3 = 8.
    [Test]
    public void MultipleAlignProgressive_Binding_MatchesHandDerivedAlignment()
    {
        var r = AlignmentTools.MultipleAlignProgressive(new[] { "ACGT", "ACGT", "AGT" });

        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "ACGT", "ACGT", "A-GT" }));
            Assert.That(r.TotalScore, Is.EqualTo(8));
        });
    }

    [Test]
    public void MultipleAlignIterative_Binding_NeverBelowProgressive()
    {
        var input = new[] { "ACGTTGCA", "ACGTGCA", "ACTTGCA", "AGTTGGCA" };
        var progressive = AlignmentTools.MultipleAlignProgressive(input);
        var iterative = AlignmentTools.MultipleAlignIterative(input);

        Assert.That(iterative.TotalScore, Is.GreaterThanOrEqualTo(progressive.TotalScore));
    }

    // Independent Python T-Coffee reference (B09 F7): CTT/ACT/CAA -> -CTT, AC-T, -CAA.
    [Test]
    public void MultipleAlignConsistency_Binding_MatchesIndependentReference()
    {
        var r = AlignmentTools.MultipleAlignConsistency(new[] { "CTT", "ACT", "CAA" });

        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "-CTT", "AC-T", "-CAA" }));
    }

    [Test]
    public void MultipleAlignAnchored_Binding_IdenticalSequencesAreGapFree()
    {
        var r = AlignmentTools.MultipleAlignAnchored(new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT" });

        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.All.EqualTo("ACGTACGT"));
            Assert.That(r.TotalScore, Is.EqualTo(24));
        });
    }

    // B09 F51 (A1-15): guide-tree distance source. The set below is discriminating: ClustalW 2.1's pairwise stage and the
    // repository's linear NW alignment give different distances, hence different guide trees and different rows
    // (rows of the independent Python implementation in scripts/alignment/generate_clustalw_pairwise_oracle.py).
    private static readonly string[] DistanceSet = ["CCAACCAACCAA", "ACACTCCGCCCGCACTCAGAACTTT", "TAA"];

    [Test]
    public void MultipleAlignProgressive_DistanceClustalW_UsesClustalWPairwiseStage()
    {
        var clustalw = AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "clustalw");
        var linear = AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3);
        Assert.Multiple(() =>
        {
            Assert.That(clustalw.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "----------------T--AA----" }));
            Assert.That(linear.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "---------------T---AA----" }));
            Assert.That(AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "LINEAR").AlignedSequences, Is.EqualTo(linear.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignIterativeAndConsistency_Distance_ParsedAndValidated()
    {
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "clustalw"));
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "clustalw"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "upgma"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "bogus"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "bogus"));
    }

    // B09 F55 (A1-15b): distance = clustalw-quick (clustalw -QUICKTREE) delegates to MsaDistanceMode.ClustalWQuick with the
    // ktuple / window / topDiagonals / pairGap parameters; the guide-tree stage validates them.
    [Test]
    public void MultipleAlignProgressive_DistanceClustalWQuick_DelegatesAndValidatesParameters()
    {
        string[] set = ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT", "ACGTAGGATTGACCGTAGCCATCGATGGA"];
        var dna = set.Select(x => new global::Seqeron.Genomics.Core.DnaSequence(x)).ToList();
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(5, -4, 0, -3);
        var expected = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignProgressive(dna, scoring,
            new global::Seqeron.Genomics.Infrastructure.MsaOptions
            {
                Distance = global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.ClustalWQuick,
                QuickPairwise = new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions { KTuple = 3, Window = 5, TopDiagonals = 6, PairGap = 9 },
                GuideTree = global::Seqeron.Genomics.Infrastructure.MsaGuideTree.NeighborJoining,
            });
        var r = AlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "clustalw-quick", "linear", "nj", 3, 5, 6, 9);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(expected.AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "ClustalW_Quick").AlignedSequences, Has.Length.EqualTo(4));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "clustalw-quick", "linear", "upgma", 9));
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignIterative(set, 5, -4, 0, -3, 16, "clustalw-quick"));
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignConsistency(set, 5, -4, 0, -3, "clustalw-quick"));
        });
    }

    // B09 F56 (A1-10b): multiple_align_clustalw delegates to SequenceAligner.MultipleAlignClustalW (ClustalW 2.1 multiple alignment).
    // Rows of the set below = the installed clustalw 2.1 (-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT).
    private static readonly string[] ClustalWSet =
        ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT", "ACGTAGGATTGACCGTAGCCATCGATGGA", "GGGTTTACCAGTAGGCATCGATGGAAACCT"];

    [Test]
    public void MultipleAlignClustalW_DelegatesToTheClustalWPort_AndParsesEveryOption()
    {
        var dna = ClustalWSet.Select(x => new global::Seqeron.Genomics.Core.DnaSequence(x)).ToList();
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(1, -1, -2, -1);
        var expected = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(dna, null, scoring);
        var r = AlignmentTools.MultipleAlignClustalW(ClustalWSet);
        var opt = new global::Seqeron.Genomics.Infrastructure.ClustalWOptions
        {
            GapOpen = 10, GapExtend = 3, PairwiseMode = global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMode.Quick, PairwiseGapOpen = 12, PairwiseGapExtend = 5,
            TransitionWeight = 1f, DivergenceCutoff = 50, NoWeights = true, NegativeMatrix = true,
            Quick = new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions { KTuple = 3, Window = 6, TopDiagonals = 7, PairGap = 10 },
        };
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(dna, opt, new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(2, -3, -5, -1));
        var viaTool = AlignmentTools.MultipleAlignClustalW(ClustalWSet, 10, 3, "ClustalW-Quick", 12, 5, 1f, 50, true, true, 3, 6, 7, 10, 2, -3, -5, -1);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(expected.AlignedSequences));
            Assert.That(r.TotalScore, Is.EqualTo(expected.TotalScore));
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(ClustalWSet));
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(AlignmentTools.MultipleAlignClustalW(ClustalWSet, distance: "full").AlignedSequences, Is.EqualTo(r.AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignClustalW(Array.Empty<string>()).AlignedSequences, Is.Empty);
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWSet, distance: "linear"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(new[] { "ACGT", "ACGX" }));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWSet, gapOpen: -3));
        });
    }

    // B09 F52 (A1-10): gapModel = affine honours gapOpen. Identical pair 33 matches = 66; the two pairs with the 23-nt sequence have
    // 23 matches (46) and ONE 10-nt gap run: affine -8 + 10 * -1 = -18, linear -10.
    private static readonly string[] IndelSet =
        ["ACGTACGTACGTAAAAAAAAAAGGCATCGATCG", "ACGTACGTACGTGGCATCGATCG", "ACGTACGTACGTAAAAAAAAAAGGCATCGATCG"];

    [Test]
    public void MultipleAlign_GapModelAffine_ChargesOneOpeningPerGapRun()
    {
        var linear = AlignmentTools.MultipleAlign(IndelSet, 2, -3, -8, -1);
        var affine = AlignmentTools.MultipleAlign(IndelSet, 2, -3, -8, -1, "affine");
        Assert.Multiple(() =>
        {
            Assert.That(linear.TotalScore, Is.EqualTo(66 + 2 * (46 - 10)));
            Assert.That(affine.TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(AlignmentTools.MultipleAlign(IndelSet, 2, -3, -8, -1, "AFFINE").TotalScore, Is.EqualTo(affine.TotalScore));
        });
    }

    [Test]
    public void MultipleAlignTools_GapModel_AcceptedAndValidated()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MultipleAlignAnchored(IndelSet, 2, -3, -8, -1, "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(AlignmentTools.MultipleAlignProgressive(IndelSet, 2, -3, -8, -1, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(AlignmentTools.MultipleAlignIterative(IndelSet, 2, -3, -8, -1, 16, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(AlignmentTools.MultipleAlignConsistency(IndelSet, 2, -3, -8, -1, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(IndelSet, 2, -3, -8, -1, "quadratic"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignAnchored(IndelSet, 2, -3, -8, -1, "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(IndelSet, 2, -3, -8, -1, "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(IndelSet, 2, -3, -8, -1, 16, "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(IndelSet, 2, -3, -8, -1, "linear", "bogus"));
        });
    }

    // B09 F53 (A1-12): guideTree = upgma | nj. The set is discriminating: UPGMA and ClustalW-style NJ guide trees give different
    // progressive rows (library call SequenceAligner.MultipleAlignProgressive(..., new MsaOptions { GuideTree = ... }), defaults 1/-1/-2/-1).
    // F56 (A1-10b): the NJ tree now uses ClustalW's own NJ arithmetic (exact ties like ClustalW); the former set became a tie where
    // both trees agree, so the set was replaced by one found by a random search (NJ and UPGMA topologies differ).
    private static readonly string[] GuideTreeSet =
        ["AAAAGTCAGGCAT", "AACAGTTAGGCAT", "AAAGTTAGGCCA", "AAATGTTAGGCAT", "AAAGTTAGGCAA"];

    [Test]
    public void MultipleAlignProgressive_GuideTreeNj_SelectsNeighborJoiningTree()
    {
        var upgma = AlignmentTools.MultipleAlignProgressive(GuideTreeSet);
        var nj = AlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "nj");
        Assert.Multiple(() =>
        {
            Assert.That(upgma.AlignedSequences, Is.EqualTo(new[] { "AA-AAGTCAGG-CAT", "AACA-GTTAGG-CAT", "AA-A-GTTAGGCCA-", "AA-ATGTTAGG-CAT", "AA-A-GTTAGG-CAA" }));
            Assert.That(nj.AlignedSequences, Is.EqualTo(new[] { "AA-AAGTCAGGCAT", "AACA-GTTAGGCAT", "AA-A-GTTAGGCCA", "AA-ATGTTAGGCAT", "AA-A-GTTAGGCAA" }));
            Assert.That(upgma.TotalScore, Is.EqualTo(83));
            Assert.That(nj.TotalScore, Is.EqualTo(81));
            Assert.That(AlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "UPGMA").AlignedSequences, Is.EqualTo(upgma.AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "neighbor_joining").AlignedSequences, Is.EqualTo(nj.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignIterativeAndConsistency_GuideTree_ParsedAndValidated()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignIterative(GuideTreeSet, 1, -1, -2, -1, 16, "clustalw", "linear", "nj"));
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignConsistency(GuideTreeSet, 1, -1, -2, -1, "clustalw", "linear", "nj"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(GuideTreeSet, 1, -1, -2, -1, 16, "linear", "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "bogus"));
        });
    }

    // B09 F58 (A1-13): consistencyEngine = 'tcoffee' routes the progressive stage through the T-Coffee 13.41 port; default unchanged.
    [Test]
    public void MultipleAlignConsistency_ConsistencyEngine_TCoffeeIsValid_DefaultUnchanged_BogusRejected()
    {
        var raw = new[] { "ACGTTACGT", "ACGACGT", "ACGTACGA", "TTACGTAG" };
        var tcoffee = AlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "tcoffee");
        var paper = AlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "paper");

        Assert.Multiple(() =>
        {
            Assert.That(tcoffee.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(raw));
            Assert.That(paper.AlignedSequences, Is.EqualTo(AlignmentTools.MultipleAlignConsistency(raw).AlignedSequences));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "bogus"));
        });
    }

    // B09 F61 (A1-14): type = 'protein' / matrix on every multiple_align* tool; the tools delegate to the SequenceAligner matrix overloads.
    private static readonly string[] Proteins =
    [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
    ];

    [Test]
    public void MultipleAlign_Protein_DelegatesToTheMatrixOverloads_AllMethods()
    {
        var blosum = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var affine = new global::Seqeron.Genomics.Infrastructure.MsaOptions { GapModel = global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine };
        var linear = new global::Seqeron.Genomics.Infrastructure.MsaOptions();
        Assert.Multiple(() =>
        {
            var star = AlignmentTools.MultipleAlign(Proteins, 1, -1, -9.5, -0.5, "affine", "protein");
            var starExpected = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlign(Proteins, blosum, affine);
            Assert.That(star.AlignedSequences, Is.EqualTo(starExpected.AlignedSequences));
            Assert.That(star.ExactTotalScore, Is.EqualTo(starExpected.ExactTotalScore).Within(1e-9));
            Assert.That(star.TotalScore, Is.EqualTo(starExpected.TotalScore));
            Assert.That(star.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(Proteins));

            Assert.That(AlignmentTools.MultipleAlignAnchored(Proteins, 1, -1, -9.5, -0.5, "affine", "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignAnchored(Proteins, blosum, affine).AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "linear", "affine", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignProgressive(Proteins, blosum, affine).AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignIterative(Proteins, 1, -1, -9.5, -0.5, 16, "linear", "linear", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignIterative(Proteins, blosum, linear).AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignConsistency(Proteins, 1, -1, -9.5, -0.5, "linear", "affine", "upgma", null, null, null, null, "paper", "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignConsistency(Proteins, blosum, affine).AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignExact_Protein_MatchesTheLibrary_AndReportsTheRealScore()
    {
        string[] small = ["MKVLAAGIVG", "MKVLSAGIVG", "MKVAAGLVG"];
        var r = AlignmentTools.MultipleAlignExact(small, 1, -1, -9.5, -0.5, "affine", true, 50_000_000, "protein", "BLOSUM62");
        var expected = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignExactWithStatistics(small,
            new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62, -9.5, -0.5),
            new global::Seqeron.Genomics.Infrastructure.MsaOptions { GapModel = global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine });
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(expected.Alignment.AlignedSequences));
            Assert.That(r.ExactTotalScore, Is.EqualTo(expected.Alignment.ExactTotalScore).Within(1e-9));
            Assert.That(r.ScoreScale, Is.EqualTo(2), "gap scores in halves");
            Assert.That(r.LowerBound, Is.EqualTo(expected.LowerBound));
            Assert.That(r.LatticeCells, Is.EqualTo(expected.LatticeCells));
        });
        // DNA stays integer: ExactTotalScore = TotalScore, scale 1.
        var dna = AlignmentTools.MultipleAlignExact(new[] { "ACGTACGT", "ACGTACGT", "ACGAACGT" });
        Assert.That(dna.ExactTotalScore, Is.EqualTo(dna.TotalScore));
        Assert.That(dna.ScoreScale, Is.EqualTo(1));
    }

    [Test]
    public void MultipleAlign_TypeAndMatrix_AreValidated()
    {
        Assert.Multiple(() =>
        {
            // protein defaults to BLOSUM62; PAM250 / lower-case type accepted
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlign(Proteins, 1, -1, -10, -1, "affine", "PROTEIN", "pam250"));
            Assert.That(AlignmentTools.MultipleAlign(Proteins, 1, -1, -10, -1, "affine", "protein").AlignedSequences,
                Is.EqualTo(AlignmentTools.MultipleAlign(Proteins, 1, -1, -10, -1, "affine", "protein", "BLOSUM62").AlignedSequences));
            // DNA with EDNAFULL accepts IUPAC codes
            var iupac = AlignmentTools.MultipleAlign(new[] { "ACGTNACGT", "ACGTAACGT", "ACRTAACGT" }, 1, -1, -10, -0.5, "affine", "dna", "EDNAFULL");
            Assert.That(iupac.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "ACGTNACGT", "ACGTAACGT", "ACRTAACGT" }));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(Proteins, 1, -1, -2, -1, "linear", "rna"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACGT" }, 1, -1, -2, -1, "linear", "dna", "BLOSUM62"), "protein matrix with DNA type");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(Proteins, 1, -1, -2, -1, "linear", "protein", "EDNAFULL"), "nucleotide matrix with protein type");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(Proteins, 1, -1, -2, -1, "linear", "protein", "BLOSUM99"), "unknown matrix");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(new[] { "MKVL", "MKXL1" }, 1, -1, -2, -1, "linear", "protein"), "residue outside the matrix alphabet");
            // fractional gap scores need a matrix
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACGT" }, 1, -1, -2, -0.5));
            // empty set, protein: empty MSA
            Assert.That(AlignmentTools.MultipleAlign(Array.Empty<string>(), 1, -1, -2, -1, "linear", "protein").AlignedSequences, Is.Empty);
        });
    }
}
