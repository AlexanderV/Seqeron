using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10c part a2, B09 F70): the remaining options of MUSCLE 3.8.1551 - iterations 3 and above (restricted-partition refinement, <c>refinehoriz.cpp</c> /
/// <c>refinevert.cpp</c> / anchors), the objective scores, the sequence weightings (Henikoff, Henikoff-PB, GSC, three-way), the distances of the first tree (kmer20_3, kmer20_4, kbit20_3),
/// neighbour-joining trees and their rooting, <c>-usetree</c>, <c>-matrix</c>, the program's tree / distance / weight files and the MSA command lines (<c>-profile</c>, <c>-refine</c>,
/// <c>-spscore</c>, <c>-ppscore</c>, <c>-maketree</c>, <c>-refinew</c>, <c>-sw</c>). They are implemented from the source of the program; the verification against the binary is the second
/// Docker round (docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2, AWAITING the user's run), so what is locked here are the published definitions and hand traces of the source on
/// small examples (each test names its source function), and structural invariants of the results.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MuscleOptions_Tests
{
    private static readonly string[] Family =
    [
        "MKVLAAGIVGLLLAQPWTTQAEAKP", "MKVLAGIVGLLAQPWTQAEAKPLL", "MRVLAAGIIGLLLQPWTTQAEAK", "KVLAAGIVALLLAQPWTTAEAKPS", "MKVAAGIVGLLLAQPWTTQEAKP",
        "MKVLAAGIVGLLAQPWTTQAEKPS", "MRVLAGIIGLLLQPWTTQAEAKPL",
    ];

    private static readonly string[] Dna =
    [
        "ACGTACGTAGCTAGCTAGGATCCGATCGATCGGATCGATTAGC", "ACGTACGAGCTAGCTAGGATCCGATCGTCGGATCGATTAGC", "ACGTCGTAGCTAGCTAGGATCCGATCGATCGGATCTATTAGC", "ACGTACGTAGCTAGGATCCGATCGATCGGATCGATTAGCA",
        "ACGAACGTAGCTAGCTAGGATCCGATCGATGGATCGATTAGC",
    ];

    private const string Blosum62Text = """
        # BLOSUM62 (columns in the order of the file, '*' last)
           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
        R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
        N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
        D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
        C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
        Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
        E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
        G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
        H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
        I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
        L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
        K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
        M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
        F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
        P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
        S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
        T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
        W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
        Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
        V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
        B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
        Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
        X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
        * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
        """;

    private static SequenceAligner.MuscleEngine Engine(MuscleOptions? o = null, params string[] seqs) =>
        new(o ?? new MuscleOptions { SequenceType = MuscleSequenceType.Protein }, seqs.Length == 0 ? ["ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK"] : seqs);

    private static SequenceAligner.MuscleMsa Msa(params string[] rows) => new(Enumerable.Range(0, rows.Length).ToArray(), rows.Select(r => r.ToCharArray()).ToArray());

    private static string[] Align(string[] seqs, MuscleOptions o) => SequenceAligner.MultipleAlignMuscle(seqs, o).AlignedSequences.ToArray();

    private static void AssertValidAlignment(string[] seqs, string[] rows, string context)
    {
        Assert.That(rows, Has.Length.EqualTo(seqs.Length), context);
        Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1), context + ": equal row lengths");
        for (int i = 0; i < seqs.Length; i++)
            Assert.That(rows[i].Replace("-", "").ToUpperInvariant(), Is.EqualTo(seqs[i].ToUpperInvariant()), $"{context}: row {i}");
        for (int c = 0; c < rows[0].Length; c++)
            Assert.That(rows.Any(r => r[c] != '-'), Is.True, $"{context}: all-gap column {c}");
    }

    // ---- sequence weights (msa2.cpp SetMSAWeightsMuscle, henikoffweight*.cpp, setgscweights.cpp, threewaywt.cpp) --------------------------------------------------------------------

    /// <summary>Henikoff and Henikoff 1994: a column with r letters gives each sequence 1 / (r x copies of its letter); the weights are normalised to 1. Hand computation for AC / AD / ED.</summary>
    [Test]
    public void HenikoffWeights_AreThePositionBasedWeightsOfHenikoff1994()
    {
        var e = Engine();
        e.UseWeighting(MuscleSeqWeight.Henikoff);
        var m = Msa("AC", "AD", "ED");
        e.SetMsaWeights(m);
        // column 1: A x2, E x1 (r = 2): 1/4, 1/4, 1/2; column 2: C x1, D x2: 1/2, 1/4, 1/4 -> 0.75, 0.5, 0.75 (sum 2)
        Assert.That(m.Weights, Is.EqualTo(new[] { 0.375f, 0.25f, 0.375f }).Within(1e-6f));
        var conserved = Msa("WAC", "WAD", "WED");
        e.SetMsaWeights(conserved);   // the conserved column gives 1/3 to everybody
        Assert.That(conserved.Weights[0], Is.EqualTo((1f / 3 + 0.75f) / 3f).Within(1e-6f));
    }

    /// <summary>The PSI-BLAST variant skips perfectly conserved columns, counts the gap as a letter and gives an all-gap sequence weight 0.</summary>
    [Test]
    public void HenikoffPbWeights_SkipConservedColumns_AndZeroTheGapSequences()
    {
        var e = Engine();
        e.UseWeighting(MuscleSeqWeight.HenikoffPb);
        var m = Msa("WAC", "WAD", "WED");
        e.SetMsaWeights(m);
        Assert.That(m.Weights, Is.EqualTo(new[] { 0.375f, 0.25f, 0.375f }).Within(1e-6f));
        var withGapSeq = Msa("AC", "AD", "ED", "--");
        e.SetMsaWeights(withGapSeq);
        Assert.That(withGapSeq.Weights[3], Is.EqualTo(0f));
        var identical = Msa("AC", "AC", "AC");
        e.SetMsaWeights(identical);   // every column is conserved: all weights 1, normalised
        Assert.That(identical.Weights, Is.EqualTo(new[] { 1f / 3, 1f / 3, 1f / 3 }).Within(1e-6f));
    }

    [Test]
    public void NoneWeights_AreUniform_AndTwoSequencesAreHalfHalf()
    {
        var e = Engine();
        e.UseWeighting(MuscleSeqWeight.None);
        var m = Msa("AC", "AD", "ED", "EE");
        e.SetMsaWeights(m);
        Assert.That(m.Weights, Is.EqualTo(new[] { 0.25f, 0.25f, 0.25f, 0.25f }));
        e.UseWeighting(MuscleSeqWeight.Henikoff);
        var two = Msa("AC", "AD");
        e.SetMsaWeights(two);
        Assert.That(two.Weights, Is.EqualTo(new[] { 0.5f, 0.5f }));
    }

    /// <summary>GSC (Gerstein, Sonnhammer, Chothia; HMMer): two identical sequences and an outlier - the cluster tree joins the identical pair at distance 0, so the outlier takes half the weight.</summary>
    [Test]
    public void GscWeights_SplitTheRootWeightByTheClusterTree()
    {
        var e = Engine();
        e.UseWeighting(MuscleSeqWeight.Gsc);
        var m = Msa("ACDEFGHIKL", "ACDEFGHIKL", "WWWWWWWWWW");
        e.SetMsaWeights(m);
        Assert.That(m.Weights, Is.EqualTo(new[] { 0.25f, 0.25f, 0.5f }).Within(1e-6f));
        var four = Msa("ACDEFGHIKL", "ACDEFGHIKM", "ACDEFGHWWW", "WWWWWWWWWW");
        e.SetMsaWeights(four);
        Assert.That(four.Weights.Sum(), Is.EqualTo(1f).Within(1e-5f));
        Assert.That(four.Weights[3], Is.GreaterThan(four.Weights[0]));
    }

    /// <summary>
    /// threewaywt.cpp: the weight of a sequence is the length of its path from the split edge (the program multiplies the lengths along the way). Tree ((A:1,B:1):1,C:2), split between the root and
    /// the (A,B) node: the root is replaced by its other child C, edge C-(AB) = 2 + 1 = 3, so A = B = 1, C = 3, normalised to 0.2 / 0.2 / 0.6.
    /// </summary>
    [Test]
    public void ThreeWayWeights_FollowTheEdgeLengthsFromTheSplitEdge()
    {
        var e = Engine();
        var tree = new SequenceAligner.MuscleTree
        {
            LeafCount = 3,
            Root = 4,
            Left = [-1, -1, -1, 0, 3],
            Right = [-1, -1, -1, 1, 2],
            Parent = [3, 3, 4, 4, -1],
            EdgeLength = [1.0, 1.0, 2.0, 1.0, 0.0],
            LeafId = [0, 1, 2, -1, -1],
        };
        e.UseWeighting(MuscleSeqWeight.ThreeWay);
        e.SetMuscleTree(tree);
        e.UseSplit(4, 3);
        var m = Msa("AC", "AD", "ED");
        e.SetMsaWeights(m);
        Assert.That(m.Weights, Is.EqualTo(new[] { 0.2f, 0.2f, 0.6f }).Within(1e-6f));
        e.UseSplit(-1, -1);   // no split edge: the program falls back to the Henikoff-PB weights
        var fallback = Msa("WAC", "WAD", "WED");
        e.SetMsaWeights(fallback);
        Assert.That(fallback.Weights, Is.EqualTo(new[] { 0.375f, 0.25f, 0.375f }).Within(1e-6f));
    }

    // ---- objective scores (objscore2.cpp, spfast.cpp, diffobjscore.cpp, anchors.cpp) -----------------------------------------------------------------------------------------------

    private static SequenceAligner.MuscleEngine MatrixEngine(MuscleOptions? extra = null) =>
        Engine((extra ?? new MuscleOptions()) with { SequenceType = MuscleSequenceType.Protein, UserMatrix = Blosum62Text, GapOpen = -10, GapExtend = 0 });

    /// <summary>The sum-of-pairs objective of two sequences: w1 w2 (letter scores + gap scores), here 0.5 x 0.5 x (BLOSUM62 A-A 4 + C-C 9) = 3.25.</summary>
    [Test]
    public void ObjScoreSp_IsTheWeightedSumOfPairScores()
    {
        var e = MatrixEngine();
        e.UseWeighting(MuscleSeqWeight.None);
        var gapFree = Msa("AC", "AC");
        e.SetMsaWeights(gapFree);
        Assert.That(e.ObjScoreSp(gapFree), Is.EqualTo(3.25f).Within(1e-6f));
        // an internal gap costs the full gap open (-10): 0.25 x (4 + 9 - 10) = 0.75
        var internalGap = Msa("A-C", "AGC");
        e.SetMsaWeights(internalGap);
        Assert.That(e.ObjScoreSp(internalGap), Is.EqualTo(0.75f).Within(1e-6f));
    }

    /// <summary>-termgaps: half (default) charges half the gap open at a terminal gap, full charges nothing extra for the terminal (the gap open of the program applies to internal gaps only), ext the extension score.</summary>
    [Test]
    public void ObjScoreSp_TerminalGaps_FollowTheTermGapsOption()
    {
        float Score(MuscleTermGaps t, string a, string b)
        {
            var e = MatrixEngine(new MuscleOptions { TermGaps = t });
            e.UseWeighting(MuscleSeqWeight.None);
            var m = Msa(a, b);
            e.SetMsaWeights(m);
            return e.ObjScoreSp(m);
        }

        Assert.Multiple(() =>
        {
            Assert.That(Score(MuscleTermGaps.Half, "-AC", "GAC"), Is.EqualTo(0.25f * (13 - 5)).Within(1e-6f));
            Assert.That(Score(MuscleTermGaps.Full, "-AC", "GAC"), Is.EqualTo(0.25f * 13).Within(1e-6f));
            Assert.That(Score(MuscleTermGaps.Ext, "-AC", "GAC"), Is.EqualTo(0.25f * 13).Within(1e-6f));
            Assert.That(Score(MuscleTermGaps.Half, "AC-", "ACG"), Is.EqualTo(0.25f * (13 - 5)).Within(1e-6f));
        });
    }

    /// <summary>XP (cross pairs of two groups) equals SP minus the within-group pairs; for gap-free alignments SPF (the fast column formula) equals SP.</summary>
    [Test]
    public void ObjScore_XpSpfAndSp_AreConsistent()
    {
        var e = MatrixEngine();
        e.UseWeighting(MuscleSeqWeight.None);
        var all = Msa("ACDEFG", "ACDEFG", "ACDWFG", "ACKEFG");
        e.SetMsaWeights(all);
        float sp = e.ObjScoreSp(all);
        Assert.That(e.ObjScoreSpDimer(all), Is.EqualTo(sp).Within(1e-4f), "SPF = SP without gaps");
        var g1 = Msa("ACDEFG", "ACDEFG");
        var g2 = Msa("ACDWFG", "ACKEFG");
        e.SetMsaWeights(g1);
        e.SetMsaWeights(g2);
        // the groups are weighted separately (1/2 each), the all-sequence weights are 1/4: the cross pairs of both are w1 x w2 apart by the factor 4 / 4 = 1 (weights 1/2 x 1/2 vs 1/4 x 1/4 -> x4)
        float xp = e.ObjScoreXp(g1, g2);
        float crossSp = 0;
        for (int i = 0; i < 2; i++)
            for (int j = 2; j < 4; j++)
                crossSp += SequenceAlignerTestHelpers.PairScore(all.Rows[i], all.Rows[j]);
        Assert.That(xp, Is.EqualTo(crossSp / 4f).Within(1e-4f));
    }

    /// <summary>The anchor columns (anchors.cpp): the columns scoring above the minimum, smoothed over -smoothwindow, merged within -anchorspacing (the closest to the window start, the later one of a pair).</summary>
    [Test]
    public void FindAnchorCols_MergesTheBestColumnsLikeMergeBestCols()
    {
        var e = MatrixEngine(new MuscleOptions { AnchorSpacing = 8, MinBestColScore = 3, MinSmoothScore = 3, SmoothScoreCeil = 999 });
        e.UseWeighting(MuscleSeqWeight.None);
        var w = new string('W', 40);
        var m = Msa(w, w, w);
        e.SetMsaWeights(m);
        // every column scores 3 x (1/3)^2 x 11 = 3.67 >= 3; the smoothed score is 0 in the first and last 3 columns -> candidates 3..36
        Assert.That(e.FindAnchorCols(m), Is.EqualTo(new List<int> { 4, 12, 20, 28, 36 }));
        var shortAlignment = Msa(new string('W', 15), new string('W', 15), new string('W', 15));
        e.SetMsaWeights(shortAlignment);
        Assert.That(e.FindAnchorCols(shortAlignment), Is.Empty, "alignments shorter than 16 columns have no anchors");
    }

    // ---- distances of the first tree (fastdistkmer.cpp, fastdistjones.cpp, fastdistkbit.cpp) --------------------------------------------------------------------------------------

    [Test]
    public void DistKmer20_4_IsOneMinusTheSharedFractionOfFourMers()
    {
        var d = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance1 = MuscleDistance.Kmer20_4 }, "ACDEFGH", "ACDEFGK", "WWWWWWW").UnalignedDistances();
        // ACDE CDEF DEFG are shared: F = 3 / (7 - 4 + 1) = 0.75, distance 0.25; no shared 4-mer: F = 0 -> 0.01 -> 0.99
        Assert.That(d[1, 0], Is.EqualTo(0.25f).Within(1e-6f));
        Assert.That(d[2, 0], Is.EqualTo(0.99f).Within(1e-6f));
    }

    [Test]
    public void DistKmer20_3_IsTheSharedTripleCountPerMinimumLength()
    {
        var d = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance1 = MuscleDistance.Kmer20_3 }, "ACDEFGH", "ACDEFGK", "WWWWWWW").UnalignedDistances();
        // ACD CDE DEF EFG are shared: 4 / (7 - 2) = 0.8 (the program stores this similarity as the "distance"); nothing shared -> 1
        Assert.That(d[1, 0], Is.EqualTo(0.8f).Within(1e-6f));
        Assert.That(d[2, 0], Is.EqualTo(1f));
    }

    /// <summary>The k-bit distance counts the bits set in either vector (b &amp; 0x101 of the program) per minimum length; identical sequences give the number of their 3-mers (every second residue is skipped).</summary>
    [Test]
    public void DistKbit20_3_CountsTheUnionOfTheBitVectors()
    {
        var d = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance1 = MuscleDistance.Kbit20_3 }, "ACDEFGHIK", "ACDEFGHIK", "WWWWWWWWW").UnalignedDistances();
        // ACDEFGHIK: bits of the 3-mers ACD (22), then F-H step ... : 4 set bits (traced from the loop of SetKmerBitVector); identical sequences: 4 / 9; WWWWWWWWW has 1 3-mer: union 5 / 9
        Assert.That(d[1, 0], Is.EqualTo(4f / 9f).Within(1e-6f));
        Assert.That(d[2, 0], Is.EqualTo(5f / 9f).Within(1e-6f));
    }

    [Test]
    public void ScoreDist_OfIdenticalSequencesIsZero_AndGrowsWithDivergence()
    {
        var e = Engine(null, "ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY");
        var a = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();
        var b = "ACDEFGHIKLMNPQRSTVWF".ToCharArray();
        var c = "WWDEFGHIKLMNPQRSTVAA".ToCharArray();
        Assert.Multiple(() =>
        {
            Assert.That(e.ScoreDist(a, a), Is.EqualTo(0).Within(1e-12));
            Assert.That(e.ScoreDist(a, b), Is.GreaterThan(0));
            Assert.That(e.ScoreDist(a, c), Is.GreaterThan(e.ScoreDist(a, b)));
        });
    }

    // ---- neighbour joining and the rooting of the tree (clust.cpp, phy3.cpp, phy.cpp) ------------------------------------------------------------------------------------------------

    /// <summary>Additive distances (the path lengths of a known tree): neighbour joining reproduces the tree, i.e. the leaf-to-leaf path lengths of its tree equal the input distances.</summary>
    private static readonly float[,] Additive =
    {
        { 0, 5, 9, 9, 8 },
        { 5, 0, 10, 10, 9 },
        { 9, 10, 0, 8, 7 },
        { 9, 10, 8, 0, 3 },
        { 8, 9, 7, 3, 0 },
    };

    private static double PathLength(SequenceAligner.MuscleTree t, int a, int b)
    {
        var up = new Dictionary<int, double>();
        double d = 0;
        for (int n = a; n != -1; n = t.Parent[n])
        {
            up[n] = d;
            d += t.EdgeLength[n];
        }

        double db = 0;
        for (int n = b; n != -1; n = t.Parent[n])
        {
            if (up.TryGetValue(n, out double da))
                return da + db;
            db += t.EdgeLength[n];
        }

        throw new InvalidOperationException("disconnected");
    }

    [Test]
    public void NeighborJoining_ReproducesAnAdditiveMatrix()
    {
        var t = SequenceAligner.MuscleEngine.NeighborJoiningTree(Additive, [0, 1, 2, 3, 4]);
        Assert.Multiple(() =>
        {
            Assert.That(t.LeafCount, Is.EqualTo(5));
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < i; j++)
                    Assert.That(PathLength(t, i, j), Is.EqualTo(Additive[i, j]).Within(1e-4), $"{i},{j}");
        });
    }

    /// <summary>Re-rooting never changes the leaf-to-leaf distances; midlongestspan puts the root in the middle of the longest leaf-leaf path; minavgleafdist balances the average distance to the leaves of the two sides.</summary>
    [Test]
    public void FixRoot_PreservesPathLengths_AndRootsAsDocumented()
    {
        var njTree = SequenceAligner.MuscleEngine.NeighborJoiningTree(Additive, [0, 1, 2, 3, 4]);
        var mid = SequenceAligner.MuscleEngine.FixRoot(njTree, MuscleRoot.MidLongestSpan);
        var avg = SequenceAligner.MuscleEngine.FixRoot(njTree, MuscleRoot.MinAvgLeafDist);
        double Depth(SequenceAligner.MuscleTree t, int leaf)
        {
            double d = 0;
            for (int n = leaf; t.Parent[n] != -1; n = t.Parent[n])
                d += t.EdgeLength[n];
            return d;
        }

        void AvgBelow(SequenceAligner.MuscleTree t, int node, double up, List<double> acc)
        {
            if (t.IsLeaf(node))
            {
                acc.Add(up);
                return;
            }

            AvgBelow(t, t.Left[node], up + t.EdgeLength[t.Left[node]], acc);
            AvgBelow(t, t.Right[node], up + t.EdgeLength[t.Right[node]], acc);
        }

        Assert.Multiple(() =>
        {
            foreach (var t in new[] { mid, avg })
            {
                Assert.That(t.NodeCount, Is.EqualTo(njTree.NodeCount));
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < i; j++)
                        Assert.That(PathLength(t, i, j), Is.EqualTo(Additive[i, j]).Within(1e-4), $"{i},{j}");
            }

            // the longest span is 1-2 (10): the root is at 5 from both leaves
            Assert.That(Depth(mid, 1), Is.EqualTo(5.0).Within(1e-4));
            Assert.That(Depth(mid, 2), Is.EqualTo(5.0).Within(1e-4));
            var left = new List<double>();
            var right = new List<double>();
            AvgBelow(avg, avg.Left[avg.Root], avg.EdgeLength[avg.Left[avg.Root]], left);
            AvgBelow(avg, avg.Right[avg.Root], avg.EdgeLength[avg.Right[avg.Root]], right);
            // phy3.cpp RootByMinAvgLeafDist places the root at x = (ADR2 - ADL1 + L) / 2 from node 1 (ADL1 / ADR2 are the average leaf distances of the two sides measured from their own end of the
            // edge, the program's naming notwithstanding): with e1 / e2 the root edges, 2 e1 - avg(side 1) = 2 e2 - avg(side 2)
            double e1 = avg.EdgeLength[avg.Left[avg.Root]], e2 = avg.EdgeLength[avg.Right[avg.Root]];
            Assert.That(2 * e1 - left.Average(), Is.EqualTo(2 * e2 - right.Average()).Within(1e-4));
            Assert.That(SequenceAligner.MuscleEngine.FixRoot(njTree, MuscleRoot.Pseudo), Is.SameAs(njTree));
        });
    }

    /// <summary>difftrees.cpp: the example of the source (trees ((A,B)k,C)j,D)i and (((A,B)p,D)n,C)m): the diffs are the subtree k and the leaves C and D.</summary>
    [Test]
    public void DiffTrees_FindsTheMarriedSubtreesOfTheSourceExample()
    {
        SequenceAligner.MuscleTree Make(int[] left, int[] right, int root)
        {
            var t = new SequenceAligner.MuscleTree { LeafCount = 4, Root = root, Left = left, Right = right, Parent = new int[7], EdgeLength = new double[7], LeafId = [0, 1, 2, 3, -1, -1, -1] };
            Array.Fill(t.Parent, -1);
            for (int n = 4; n < 7; n++)
            {
                t.Parent[left[n]] = n;
                t.Parent[right[n]] = n;
            }

            return t;
        }

        var t1 = Make([-1, -1, -1, -1, 0, 4, 5], [-1, -1, -1, -1, 1, 2, 3], 6);   // k = (A,B), j = (k,C), i = (j,D)
        var t2 = Make([-1, -1, -1, -1, 0, 4, 5], [-1, -1, -1, -1, 1, 3, 2], 6);   // p = (A,B), n = (p,D), m = (n,C)
        var (diffs, idToLeaf) = SequenceAligner.MuscleEngine.DiffTrees(t1, t2);
        Assert.Multiple(() =>
        {
            Assert.That(diffs, Is.EqualTo(new List<(int, int)> { (1, 2), (3, 4), (-1, -1), (-1, -1), (-1, -1) }));
            Assert.That(idToLeaf, Is.EqualTo(new[] { 3, 3, 4, 2 }), "A, B in the k subtree, C and D alone");
        });
        var same = SequenceAligner.MuscleEngine.DiffTrees(t1, t1);
        Assert.That(same.Diffs, Has.Count.EqualTo(1), "identical trees: the root is the only diff");
    }

    // ---- files of the program: trees, distance matrices, weights, matrix, user tree ------------------------------------------------------------------------------------------------

    [TestCase(0.1, "0.1")]
    [TestCase(1e-5, "1e-05")]
    [TestCase(0.000123456789, "0.000123457")]
    [TestCase(123456789.0, "1.23457e+08")]
    [TestCase(100000.0, "100000")]
    [TestCase(1000000.0, "1e+06")]
    [TestCase(-2.5, "-2.5")]
    [TestCase(0.0, "0")]
    public void FormatG_FollowsPrintfPercentG(double v, string expected) =>
        Assert.That(SequenceAligner.MuscleEngine.FormatG(v, 6), Is.EqualTo(expected));

    /// <summary>The three-sequence tree file (rooted Newick: "(", a newline after every token group, ":%g" edge lengths, ";").</summary>
    [Test]
    public void MuscleRun_WritesTheTreeAndDistanceFilesInTheProgramFormat()
    {
        var r = SequenceAligner.MuscleRun(["ACDEFGHIKLMNPQ", "ACDEFGHIKLMNPR", "ACDEFWHIKLMNPQ"], new MuscleOptions { SequenceType = MuscleSequenceType.Protein, SequenceNames = ["x", "y", "z"] });
        Assert.Multiple(() =>
        {
            Assert.That(r.Tree1, Does.StartWith("(\n").And.EndWith(")\n;\n"));
            Assert.That(r.Tree1!.Split('\n').Count(l => l.StartsWith('x') || l.StartsWith('y') || l.StartsWith('z')), Is.EqualTo(3));
            Assert.That(r.DistanceMatrix1, Does.StartWith("3\n"));
            Assert.That(r.DistanceMatrix1!.Split('\n')[1], Is.EqualTo("         x  "));
            Assert.That(r.Tree2, Is.Not.Null, "iteration 2 ran");
            Assert.That(r.DistanceMatrix2, Is.Not.Null);
            Assert.That(r.Names, Is.EqualTo(new[] { "x", "y", "z" }));
        });
    }

    [Test]
    public void MuscleRun_ComputeWeights_StopsAfterTheProgressiveAlignment()
    {
        var r = SequenceAligner.MuscleRun(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, ComputeWeights = true });
        var lines = r.WeightsFile!.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.Multiple(() =>
        {
            Assert.That(r.Alignment.AlignedSequences, Is.Empty);
            Assert.That(lines, Has.Length.EqualTo(Family.Length));
            Assert.That(lines.Sum(l => double.Parse(l.Split('\t')[1], System.Globalization.CultureInfo.InvariantCulture)), Is.EqualTo(1.0).Within(0.01));
        });
        var c = SequenceAligner.MuscleRun(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, ClusterOnly = true });
        Assert.That(c.Tree1, Is.Not.Null);
        Assert.That(c.Alignment.AlignedSequences, Is.Empty);
    }

    [Test]
    public void ParseUserTree_NumbersTheNodesLikeTheProgramAndChecksTheLabels()
    {
        var e = Engine(new MuscleOptions { SequenceType = MuscleSequenceType.Protein, SequenceNames = ["s0", "s1", "s2"] }, "ACDEFGH", "ACDEFGK", "WWWWWWW");
        var t = e.ParseUserTree("((s0:0.1,s1:0.2):0.3,s2:0.4);");
        Assert.Multiple(() =>
        {
            // creation order of the program: root 0, children 1 and 2, then the children 3, 4 of node 1: leaves 2 (s2), 3 (s0), 4 (s1) -> leaf indexes 0, 1, 2; internals 0 (root), 1 -> 3, 4
            Assert.That(t.LeafId.Take(3), Is.EqualTo(new[] { 2, 0, 1 }));
            Assert.That(t.Root, Is.EqualTo(3));
            Assert.That(t.EdgeLength[0], Is.EqualTo(0.4));
            Assert.That(t.EdgeLength[4], Is.EqualTo(0.3));
            Assert.Throws<ArgumentException>(() => e.ParseUserTree("(s0:0.1,s1:0.2,s2:0.4);"), "trifurcation = unrooted");
            Assert.Throws<ArgumentException>(() => e.ParseUserTree("((s0,s1),zz);"), "unknown label");
            Assert.Throws<ArgumentException>(() => e.ParseUserTree("(s0,s1);"), "tree does not match the sequences");
        });
    }

    [Test]
    public void MuscleRun_UserTree_SkipsTheTreeRefinement()
    {
        var first = SequenceAligner.MuscleRun(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 1 });
        var r = SequenceAligner.MuscleRun(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 3, UserTree = first.Tree1 });
        AssertValidAlignment(Family, r.Alignment.AlignedSequences.ToArray(), "usetree");
        Assert.That(r.Tree2, Is.Null, "iteration 2 does not run with -usetree");
        Assert.Throws<InvalidOperationException>(() => SequenceAligner.MuscleRun(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, UserTree = "((((((s0,s1),s2),s3),s4),s5),s6);" }), "no edge lengths: the ClustalW weights fail like the program");
    }

    /// <summary>readmx.cpp: the heading line starts with a space, a trailing '*' column is dropped; the program reads HeadingCount - 1 columns, so a file without '*' loses its last column.</summary>
    [Test]
    public void ReadMatrixFile_ReadsTheNcbiFormat_AndLosesTheLastColumnWithoutStar()
    {
        var e = Engine();
        var mx = e.ReadMatrixFile(Blosum62Text);
        // letters in the order of MUSCLE's alphabet A C D E F G H I K L M N P Q R S T V W Y: A-A = 4, W-W = 11, C-C = 9, D-E = 2, W-F = 1
        Assert.Multiple(() =>
        {
            Assert.That(mx[0, 0], Is.EqualTo(4f));
            Assert.That(mx[18, 18], Is.EqualTo(11f));
            Assert.That(mx[1, 1], Is.EqualTo(9f));
            Assert.That(mx[2, 3], Is.EqualTo(2f));
            Assert.That(mx[18, 4], Is.EqualTo(1f));
            Assert.That(mx[19, 19], Is.EqualTo(7f), "Y-Y");
        });
        // a file whose heading has only the 20 residue columns: the program reads HeadingCount - 1 = 19 columns, the last one (V) is lost
        var lines = Blosum62Text.Split('\n');
        var twenty = new List<string> { lines[0], "   " + string.Join("  ", lines[1].Split(' ', StringSplitOptions.RemoveEmptyEntries).Take(20)) };
        twenty.AddRange(lines.Skip(2).Take(20).Select(l => string.Join(" ", l.Split(' ', StringSplitOptions.RemoveEmptyEntries).Take(21))));
        var mx2 = e.ReadMatrixFile(string.Join('\n', twenty));
        Assert.Multiple(() =>
        {
            Assert.That(mx2[0, 0], Is.EqualTo(4f));
            Assert.That(mx2[18, 18], Is.EqualTo(11f));
            Assert.That(mx2[17, 17], Is.EqualTo(0f), "V-V (4 in the file) is in the lost last column");
            Assert.That(mx2[17, 18], Is.EqualTo(-3f), "V-W is read");
        });
        Assert.Throws<ArgumentException>(() => e.ReadMatrixFile("# no heading\nA 1 2\n"));
    }

    [Test]
    public void MuscleRun_UserMatrix_ReplacesTheScoreMatrixOfTheSpScore()
    {
        var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 4 };
        var withMatrix = SequenceAligner.MuscleRun(Family, o with { UserMatrix = Blosum62Text });
        var sp = SequenceAligner.MuscleRun(Family, o with { ScoreModel = MuscleScoreModel.SumOfPairs });
        AssertValidAlignment(Family, withMatrix.Alignment.AlignedSequences.ToArray(), "matrix");
        AssertValidAlignment(Family, sp.Alignment.AlignedSequences.ToArray(), "sp");
        Assert.Throws<ArgumentException>(() => SequenceAligner.MuscleRun(Dna, new MuscleOptions { SequenceType = MuscleSequenceType.Dna, UserMatrix = Blosum62Text }));
    }

    // ---- options that the 3.8.1551 binary ignores (compiled out) -------------------------------------------------------------------------------------------------------------------

    [Test]
    public void Diags_Hydro_DoubleAffine_AndResourceOptions_HaveNoEffect()
    {
        var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 4 };
        var baseline = Align(Family, o);
        Assert.Multiple(() =>
        {
            Assert.That(Align(Family, o with { Diags1 = true, Diags2 = true, DiagLength = 12, DiagBreak = 3, DiagMargin = 9 }), Is.EqualTo(baseline), "-diags: the diagonal DP is compiled out (VER_3_52 = 0)");
            Assert.That(Align(Family, o with { Hydro = 9, HydroFactor = 3.0, TomHydro = true }), Is.EqualTo(baseline), "-hydro: #if HYDRO");
            Assert.That(Align(Family, o with { GapOpen2 = -5, GapExtend2 = -1, GapAmbig = 0.5 }), Is.EqualTo(baseline), "double affine gaps: DOUBLE_AFFINE = 0");
            Assert.That(Align(Family, o with { MaxMb = 1, MaxHours = 0.001 }), Is.EqualTo(baseline));
        });
    }

    // ---- the refinement iterations -------------------------------------------------------------------------------------------------------------------------------------------------

    [Test]
    public void Iterations3AndAbove_GiveValidAlignments_ForEveryObjectiveWeightAndAnchorSetting()
    {
        var failures = new List<string>();
        foreach (var objective in new[] { MuscleObjectiveScore.Spm, MuscleObjectiveScore.Sp, MuscleObjectiveScore.Dp, MuscleObjectiveScore.Xp, MuscleObjectiveScore.Spf, MuscleObjectiveScore.Ps })
        foreach (var weight in new[] { MuscleSeqWeight.ClustalW, MuscleSeqWeight.None, MuscleSeqWeight.Henikoff, MuscleSeqWeight.HenikoffPb, MuscleSeqWeight.Gsc, MuscleSeqWeight.ThreeWay })
        foreach (bool anchors in new[] { true, false })
        {
            var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 5, ObjectiveScore = objective, Weight2 = weight, Anchors = anchors };
            try
            {
                AssertValidAlignment(Family, Align(Family, o), $"{objective}/{weight}/{anchors}");
            }
            catch (Exception ex) when (ex is AssertionException or InvalidOperationException)
            {
                failures.Add($"{objective}/{weight}/{anchors}: {ex.Message.Split('\n')[0]}");
            }
        }

        Assert.That(failures, Is.Empty);
    }

    [Test]
    public void Iterations3AndAbove_ChangeSomeAlignmentsOfIteration2_AndNucleotidesAreSupported()
    {
        int changed = 0;
        var rnd = new Random(5);
        for (int trial = 0; trial < 12; trial++)
        {
            var seqs = Enumerable.Range(0, 5).Select(_ => new string(Enumerable.Range(0, 30 + rnd.Next(20)).Select(_ => "ACDEFGHIKLMNPQRSTVWY"[rnd.Next(20)]).ToArray())).ToArray();
            var baseSeq = seqs[0];
            seqs = seqs.Select((s, i) => i == 0 ? s : Mutate(baseSeq, rnd)).ToArray();
            var r2 = Align(seqs, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 2 });
            var r8 = Align(seqs, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 8 });
            AssertValidAlignment(seqs, r8, $"trial {trial}");
            if (!r2.SequenceEqual(r8))
                changed++;
        }

        Assert.That(changed, Is.GreaterThan(0));
        AssertValidAlignment(Dna, Align(Dna, new MuscleOptions { SequenceType = MuscleSequenceType.Dna, MaxIterations = 6 }), "dna");
    }

    private static string Mutate(string s, Random rnd)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            double r = rnd.NextDouble();
            if (r < 0.1)
                sb.Append("ACDEFGHIKLMNPQRSTVWY"[rnd.Next(20)]);
            else if (r is >= 0.15 and < 0.2)
                sb.Append(c).Append("ACDEFGHIKLMNPQRSTVWY"[rnd.Next(20)]);
            else if (r >= 0.2)
                sb.Append(c);
        }

        return sb.ToString();
    }

    // ---- trees, distances and rooting reachable from the options --------------------------------------------------------------------------------------------------------------------

    [Test]
    public void ClusteringDistanceAndRootOptions_GiveValidAlignments()
    {
        var failures = new List<string>();
        var combos = new List<MuscleOptions>();
        var baseO = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 4 };
        foreach (var d in new[] { MuscleDistance.Kmer6_6, MuscleDistance.Kmer20_3, MuscleDistance.Kmer20_4, MuscleDistance.Kbit20_3, MuscleDistance.PwKimura, MuscleDistance.PwScoreDist })
            combos.Add(baseO with { Distance1 = d });
        foreach (var d in new[] { MuscleDistance.PctIdKimura, MuscleDistance.PctIdLog, MuscleDistance.Edit, MuscleDistance.ScoreDist })
            combos.Add(baseO with { Distance2 = d });
        foreach (var c in new[] { MuscleClusterMethod.NeighborJoining, MuscleClusterMethod.Upgma })
        foreach (var root in new[] { MuscleRoot.Pseudo, MuscleRoot.MidLongestSpan, MuscleRoot.MinAvgLeafDist })
            combos.Add(baseO with { Cluster1 = c, Root1 = root, Cluster2 = c, Root2 = root, Distance2 = c == MuscleClusterMethod.NeighborJoining ? MuscleDistance.PctIdKimura : MuscleDistance.Default });
        combos.Add(baseO with { Weight1 = MuscleSeqWeight.Henikoff, MaxTrees = 3 });
        combos.Add(baseO with { Weight1 = MuscleSeqWeight.None, Weight2 = MuscleSeqWeight.ThreeWay, MaxTrees = 2 });
        combos.Add(baseO with { Brenner = true });
        combos.Add(baseO with { TermGaps = MuscleTermGaps.Full });
        combos.Add(baseO with { TermGaps = MuscleTermGaps.Ext });
        foreach (var o in combos)
        {
            try
            {
                AssertValidAlignment(Family, Align(Family, o), "combo");
            }
            catch (Exception ex) when (ex is AssertionException or InvalidOperationException or ArgumentException)
            {
                failures.Add($"{o.Distance1}/{o.Distance2}/{o.Cluster1}/{o.Root1}/{o.Weight1}: {ex.Message.Split('\n')[0]}");
            }
        }

        Assert.That(failures, Is.Empty);
    }

    [Test]
    public void DistanceOfTheWrongKind_IsRejectedLikeTheProgram()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance1 = MuscleDistance.PctIdLog }));
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Distance2 = MuscleDistance.Kmer6_6 }));
            Assert.Throws<ArgumentException>(() => Align(Family, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, Cluster2 = MuscleClusterMethod.NeighborJoining, Distance2 = MuscleDistance.Edit }));
        });
    }

    [Test]
    public void Brenner_AlignsEverySequenceToTheRootProfile_InInputOrder()
    {
        var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 1, Brenner = true };
        var r = SequenceAligner.MuscleRun(Family, o);
        AssertValidAlignment(Family, r.Alignment.AlignedSequences.ToArray(), "brenner");
        Assert.That(r.OutputOrder, Is.EqualTo(Enumerable.Range(0, Family.Length)), "the root MSA of -brenner is in input order");
    }

    // ---- MSA command lines -----------------------------------------------------------------------------------------------------------------------------------------------------------------

    [Test]
    public void ProfileAlign_AlignsTwoProfiles_KeepingBothAlignmentsIntact()
    {
        var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein };
        var r = SequenceAligner.MuscleProfileAlign(["ACDEFGHIK", "ACDEFGHIK"], ["ACDEGHIK", "ACDEGHIK"], o);
        Assert.That(r, Is.EqualTo(new[] { "ACDEFGHIK", "ACDEFGHIK", "ACDE-GHIK", "ACDE-GHIK" }));
        var gapped = SequenceAligner.MuscleProfileAlign(["AC-DEFGHIK", "ACGDEFGHIK"], ["ACDEFGHIK"], o);
        Assert.That(gapped.Take(2).Select(x => x.Replace("-", "")), Is.EqualTo(new[] { "ACDEFGHIK", "ACGDEFGHIK" }));
        Assert.That(gapped.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
        Assert.Throws<ArgumentException>(() => SequenceAligner.MuscleProfileAlign(["AC", "ACD"], ["AC"], o));
    }

    [Test]
    public void Refine_RealignsAMisalignedProfile_AndKeepsTheResidues()
    {
        var bad = new[] { "ACDEFGHIKLMNPQRSTVWY--", "ACDEFGHIKLMNPQRSTVWY--", "--ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY--", "-ACDEFGHIKLMNPQRSTVWY-" };
        var r = SequenceAligner.MuscleRefine(bad, new MuscleOptions { SequenceType = MuscleSequenceType.Protein, MaxIterations = 8 });
        Assert.Multiple(() =>
        {
            for (int i = 0; i < bad.Length; i++)
                Assert.That(r.First(x => x.Replace("-", "") == bad[i].Replace("-", "")), Is.Not.Null);
            Assert.That(r.Select(x => x.Length).Distinct().Count(), Is.EqualTo(1));
            Assert.That(r[0].Length, Is.LessThan(bad[0].Length), "the shifted rows are brought into register: fewer columns");
        });
    }

    [Test]
    public void ScoreAlignment_IsTheSumOfPairsObjective_AndGuideTreeIsNewick()
    {
        var o = new MuscleOptions { SequenceType = MuscleSequenceType.Protein, UserMatrix = Blosum62Text, GapOpen = -10, GapExtend = 0, Weight1 = MuscleSeqWeight.None };
        Assert.That(SequenceAligner.MuscleScoreAlignment(["AC", "AC"], o), Is.EqualTo(3.25).Within(1e-6));
        Assert.That(SequenceAligner.MuscleScoreAlignment(["A-C", "AGC"], o), Is.EqualTo(0.75).Within(1e-6));
        var tree = SequenceAligner.MuscleGuideTree(["ACDEFG", "ACDEFG", "ACDEWW", "WWDEFG"], new MuscleOptions { SequenceType = MuscleSequenceType.Protein });
        Assert.That(tree, Does.StartWith("(\n").And.EndWith(")\n;\n"));
    }

    [Test]
    public void NucleotideOnlyModes_AreRejectedForProteins_AndWorkForNucleotides()
    {
        var p = new MuscleOptions { SequenceType = MuscleSequenceType.Protein };
        var d = new MuscleOptions { SequenceType = MuscleSequenceType.Dna };
        string[] a = ["ACGTACGTAGCTAGCT", "ACGTACGAAGCTAGCT"], b = ["ACGTACGTAGCTAGCT", "ACGAACGTAGCTAGCT"];
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MuscleProfileScore(["ACDE"], ["ACDE"], p));
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MuscleLocalAlign(["ACDE"], ["ACDE"], p));
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MuscleRefineWindows(["ACDE"], p));
            Assert.That(SequenceAligner.MuscleProfileScore(a, b, d), Is.Not.EqualTo(0));
            Assert.That(SequenceAligner.MuscleProfileScore(a, a, d), Is.GreaterThan(SequenceAligner.MuscleProfileScore(a, ["TTTTTTTTTTTTTTTT", "TTTTTTTTTTTTTTTT"], d)), "a profile scores better against itself than against unrelated sequences");
        });
    }

    /// <summary>-sw: the local alignment of GGGGTACGTACGTAAAA and ACGTACGTACGTTTGACCA contains the common core GTACGTACGT; the residues outside the local region are lower case, padded with '.'.</summary>
    [Test]
    public void LocalAlign_AlignsTheCommonCore_AndLowerCasesTheTerminals()
    {
        var r = SequenceAligner.MuscleLocalAlign(["ACGTACGTACGTTTGACCA"], ["GGGGTACGTACGTAAAA"], new MuscleOptions { SequenceType = MuscleSequenceType.Dna });
        Assert.Multiple(() =>
        {
            Assert.That(r, Has.Length.EqualTo(2));
            Assert.That(r[0].Length, Is.EqualTo(r[1].Length));
            Assert.That(r[0], Does.Contain("GTACGTACGT"));
            Assert.That(r[1], Does.Contain("GTACGTACGT"));
            Assert.That(r[0].Replace(".", "").ToUpperInvariant().Replace("-", ""), Is.EqualTo("ACGTACGTACGTTTGACCA"));
            Assert.That(r[1].Replace(".", "").ToUpperInvariant().Replace("-", ""), Is.EqualTo("GGGGTACGTACGTAAAA"));
            Assert.That(r[0].Any(char.IsLower), Is.True);
        });
    }

    [Test]
    public void RefineWindows_RealignsEveryWindow_AndConcatenatesThem()
    {
        var aligned = SequenceAligner.MuscleRun(Dna, new MuscleOptions { SequenceType = MuscleSequenceType.Dna, MaxIterations = 2 }).Alignment.AlignedSequences.ToArray();
        var r = SequenceAligner.MuscleRefineWindows(aligned, new MuscleOptions { SequenceType = MuscleSequenceType.Dna, MaxIterations = 3, RefineWindow = 20 });
        AssertValidAlignment(Dna, r, "refinew");
    }
}

/// <summary>Small independent helpers of the option tests.</summary>
internal static class SequenceAlignerTestHelpers
{
    private const string Order = "ACDEFGHIKLMNPQRSTVWY";

    // BLOSUM62 pair score of two gap-free rows (the matrix of the tests; letters in MUSCLE's order)
    private static readonly int[,] B62 =
    {
        { 4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2 },
        { 0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 },
        { -2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3 },
        { -1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2 },
        { -2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3 },
        { 0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3 },
        { -2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2 },
        { -1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1 },
        { -1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2 },
        { -1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1 },
        { -1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1 },
        { -2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2 },
        { -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3 },
        { -1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1 },
        { -1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2 },
        { 1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2 },
        { 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2 },
        { 0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1 },
        { -3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2 },
        { -2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7 },
    };

    public static int PairScore(char[] a, char[] b)
    {
        int s = 0;
        for (int i = 0; i < a.Length; i++)
            s += B62[Order.IndexOf(a[i]), Order.IndexOf(b[i])];
        return s;
    }
}
