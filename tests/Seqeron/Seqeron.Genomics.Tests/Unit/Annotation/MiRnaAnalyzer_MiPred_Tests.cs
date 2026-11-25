using System.Globalization;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;
using static Seqeron.Genomics.Tests.Unit.Annotation.MiPredReferenceData;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// MIRNA-CLASSIFY-001: published MiPred pre-miRNA classifier (Jiang et al. 2007, Nucleic Acids Res 35:W339) with the authors'
/// 1000-tree randomForest model and the feature definitions of the package script microRNAcheck_parallel.pl.
/// Oracles (see MiPredReferenceData / MiPredOracleData): (1) forests trained locally with the installed R randomForest
/// package, exported with getTree and compared with R predict(type = "vote") on 560 rows incl. rows at the split points
/// and exact 30:30 ties; (2) an independent Python port of the Perl features (ViennaRNA 2.8 folds, numpy shuffles) and an
/// independent numpy traversal of the forest arrays parsed from model.RData; (3) the 331 training rows stored in model.RData.
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_MiPred_Tests
{
    private static Hairpin H(string name) => Hairpins.Single(h => h.Name == name);

    private static double[] Features(Hairpin h) =>
        ComputeMiPredFeatures(h.Sequence, h.ViennaStructure, h.ViennaMfe, h.PythonPValue ?? 1.0);

    // MP1 — hand-computed triplet frequencies of the package: GAUCC / "((.((": k=1 'A' "((." -> A,"((." (index 1);
    // k=2 'U' "(.(" -> U,"(.(" (8+2=10); k=3 'C' ".((" -> C,".((" (24+4=28); counts / LENGTH (5), positions 1 and n are not counted.
    [Test]
    public void MP1_TripletFrequencies_HandComputed()
    {
        var f = ComputeMiPredTripletFrequencies("GAUCC", "((.((");
        Assert.That(f, Has.Length.EqualTo(32));
        for (int i = 0; i < 32; i++)
            Assert.That(f[i], Is.EqualTo(i is 1 or 10 or 28 ? 0.2 : 0.0), $"V{i + 1}");
    }

    // MP1b — closing brackets are never part of a counted triplet (the package keeps "(((" .. "..." only) and DNA is read as RNA.
    [Test]
    public void MP1b_ClosingBracketTripletsAreNotCounted_DnaEqualsRna()
    {
        Assert.That(ComputeMiPredTripletFrequencies("GAUCC", "(.)).").All(v => v == 0), Is.True);
        Assert.That(ComputeMiPredTripletFrequencies("gatcc", "((.(("), Is.EqualTo(ComputeMiPredTripletFrequencies("GAUCC", "((.((")));
    }

    // MP2 — the 32 frequencies of all reference hairpins = the Python port of the Perl counting (integer counts / length), exactly.
    [TestCaseSource(nameof(AllHairpinNames))]
    public void MP2_TripletFrequencies_EqualPerlPort(string name)
    {
        var h = H(name);
        var f = ComputeMiPredTripletFrequencies(h.Sequence, h.ViennaStructure);
        for (int i = 0; i < 32; i++)
            Assert.That(f[i], Is.EqualTo(h.TripletCounts[i] / (double)h.Sequence.Length), $"{name} V{i + 1}");
    }

    private static IEnumerable<string> AllHairpinNames => Hairpins.Select(h => h.Name);

    // MP3 — pre-judgement of the package (MFE <= -15, >= 18 pairs, no multiloop, 51..137 nt) on the reference structures.
    [TestCaseSource(nameof(AllHairpinNames))]
    public void MP3_Prejudgement_EqualsPerlPort(string name)
    {
        var h = H(name);
        var c = AssessMiPredHairpin(h.ViennaStructure, h.ViennaMfe);
        Assert.That((c.BasePairs, c.HasMultiloop, c.IsHairpinLike), Is.EqualTo((h.Pairs, h.Multiloop, h.IsHairpinLike)), name);
    }

    // MP3b — each boundary of the package rule (>= 18 pairs, length 51..137, MFE <= -15, no second loop / multiloop).
    [Test]
    public void MP3b_PrejudgementBoundaries()
    {
        static string Hp(int pairs, int loop) => new string('(', pairs) + new string('.', loop) + new string(')', pairs);
        Assert.That(AssessMiPredHairpin(Hp(18, 15), -15.0).IsHairpinLike, Is.True, "51 nt, 18 pairs, -15.00");
        Assert.That(AssessMiPredHairpin(Hp(18, 15), -14.99).IsHairpinLike, Is.False, "MFE above -15");
        Assert.That(AssessMiPredHairpin(Hp(18, 14), -20).IsHairpinLike, Is.False, "50 nt");
        Assert.That(AssessMiPredHairpin(Hp(17, 17), -20).IsHairpinLike, Is.False, "17 pairs");
        Assert.That(AssessMiPredHairpin(Hp(18, 101), -20).IsHairpinLike, Is.True, "137 nt");
        Assert.That(AssessMiPredHairpin(Hp(18, 102), -20).IsHairpinLike, Is.False, "138 nt");
        string two = Hp(10, 4) + Hp(10, 4) + new string('.', 10);
        var c = AssessMiPredHairpin(two, -30);
        Assert.That((c.BasePairs, c.HasMultiloop, c.IsHairpinLike), Is.EqualTo((0, true, false)), "second hairpin: pair count undefined (0), multiloop");
        var m = AssessMiPredHairpin("(((((((((((((((((((....(((.....)))....)))))))))))))))))))....(((....)))", -40);
        Assert.That(m.HasMultiloop, Is.True, "'(' after the first ')'");
        Assert.That(AssessMiPredHairpin(new string('.', 60), -20).IsHairpinLike, Is.False, "no pairs");
    }

    // MP4 — forest votes on the Python-derived feature vectors (Vienna fold, numpy shuffle p-value) = independent numpy traversal.
    [TestCaseSource(nameof(ClassifiedHairpinNames))]
    public void MP4_ForestVotes_EqualIndependentTraversal(string name)
    {
        var h = H(name);
        var (votes, trees) = MiPredVotes(Features(h));
        Assert.That(trees, Is.EqualTo(1000));
        Assert.That(votes, Is.EqualTo(h.PythonVotes!.Value), name);
    }

    private static IEnumerable<string> ClassifiedHairpinNames => Hairpins.Where(h => h.PythonVotes is not null).Select(h => h.Name);

    // MP5 — the 60 stored training rows of model.RData: C# votes = independent traversal; in-sample class agreement.
    [Test]
    public void MP5_StoredTrainingRows_VotesAndClasses()
    {
        Assert.That(TrainRows, Has.Length.EqualTo(60));
        foreach (var r in TrainRows)
        {
            var (votes, _) = MiPredVotes(r.Features);
            Assert.That(votes, Is.EqualTo(r.Votes));
            Assert.That(votes > 500, Is.EqualTo(r.Class == "P"), r.Class);
        }
    }

    // MP6 — structure of the bundled forest: ntree 1000, nrnodes 103 (ndbigtree 9..103), 56,324 nodes = 28,662 terminal + 27,662 interior.
    [Test]
    public void MP6_BundledForest_StructuralInvariants()
    {
        var forest = BundledMiPredForest;
        Assert.That(forest.TreeCount, Is.EqualTo(1000));
        Assert.That(forest.NodeCount, Is.EqualTo(56324));
        Assert.That(forest.TerminalCount, Is.EqualTo(28662));
        Assert.That(forest.TerminalCount, Is.EqualTo(forest.NodeCount - forest.TerminalCount + forest.TreeCount), "terminals = interior + trees");
        for (int t = 0; t < forest.TreeCount; t++)
        {
            Assert.That(forest.TreeNodeCount(t), Is.InRange(9, 103), $"tree {t}");
            Assert.That(forest.TreeNodeCount(t) % 2, Is.EqualTo(1), $"tree {t}: full binary tree has an odd node count");
        }
    }

    // MP6b — corrupt tables are rejected (cycle, out-of-range daughter, trailing bytes, bad code).
    [Test]
    public void MP6b_Decode_RejectsCorruptTables()
    {
        byte[] counts = { 3, 0 };
        byte[] ok = Node(1, 0.5, 1, 2).Concat(new byte[] { 101, 102 }).ToArray();
        Assert.That(MiPredForest.Decode(ok, counts, 34).CountRealVotes(new double[34]), Is.EqualTo(0));
        Assert.Throws<InvalidDataException>(() => MiPredForest.Decode(Node(1, 0.5, 0, 2).Concat(new byte[] { 101, 102 }).ToArray(), counts, 34), "daughter not after parent");
        Assert.Throws<InvalidDataException>(() => MiPredForest.Decode(Node(1, 0.5, 1, 3).Concat(new byte[] { 101, 102 }).ToArray(), counts, 34), "daughter out of range");
        Assert.Throws<InvalidDataException>(() => MiPredForest.Decode(ok.Concat(new byte[] { 1 }).ToArray(), counts, 34), "trailing byte");
        Assert.Throws<InvalidDataException>(() => MiPredForest.Decode(Node(35, 0.5, 1, 2).Concat(new byte[] { 101, 102 }).ToArray(), counts, 34), "variable 35");
        Assert.Throws<InvalidDataException>(() => MiPredForest.Decode(Node(1, 0.5, 1, 1).Concat(new byte[] { 101, 102 }).ToArray(), counts, 34), "both daughters equal");

        static byte[] Node(byte v, double split, ushort l, ushort r)
        {
            var b = new List<byte> { v };
            b.AddRange(BitConverter.GetBytes(split));
            b.AddRange(BitConverter.GetBytes(l));
            b.AddRange(BitConverter.GetBytes(r));
            return b.ToArray();
        }
    }

    // MP7 — traversal semantics against R itself: a 60-tree forest trained with R randomForest; predict(type = "vote") on 560 rows
    // (300 fresh, 100 training, 160 with features equal to split points = "<=" rule; 4 rows are exact 30:30 ties).
    [Test]
    public void MP7_OracleForest_VotesEqualRandomForestPredict()
    {
        var forest = MiPredForest.Decode(
            Convert.FromBase64String(string.Concat(MiPredOracleData.ForestBase64)),
            Convert.FromBase64String(string.Concat(MiPredOracleData.NodeCountsBase64)), MiPredFeatureCount);
        Assert.That(forest.TreeCount, Is.EqualTo(MiPredOracleData.TreeCount));
        var lines = MiPredOracleData.RowsCsv.Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
        Assert.That(lines, Has.Length.EqualTo(560));
        int ties = 0, boundaryHits = 0;
        for (int i = 0; i < lines.Length; i++)
        {
            double[] x = lines[i].Split(',').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
            Assert.That(x, Has.Length.EqualTo(MiPredFeatureCount));
            int votes = forest.CountRealVotes(x);
            Assert.That(votes, Is.EqualTo(MiPredOracleData.RVotesP[i]), $"row {i}");
            if (votes * 2 == MiPredOracleData.TreeCount) ties++;
            if (i >= 400 && forest.CountSplitPointHits(x) > 0) boundaryHits++;
        }
        Assert.That(ties, Is.GreaterThanOrEqualTo(4), "exact 30:30 ties are in the oracle set");
        Assert.That(boundaryHits, Is.GreaterThanOrEqualTo(100), "most boundary rows (131 of 160 in the stored set) have a feature exactly at a split point on a tree path: the \"<=\" rule is exercised");
    }

    // MP8 — end-to-end (canonical fold + 999 seeded shuffles) on the five miRBase precursors: hairpin-like, p <= 0.01, voted real.
    // Locked votes of the seeded run; |Python - C#| <= 5 trees (folds differ slightly: ViennaRNA 1999/2004 vs canonical Turner 2004).
    [TestCase("hsa-mir-21", 781)]
    [TestCase("hsa-let-7a-1", 823)]
    [TestCase("hsa-mir-16-1", 836)]
    [TestCase("hsa-mir-155", 802)]
    [TestCase("hsa-mir-30a", 784)]
    public void MP8_RealPreMiRnas_ClassifiedReal(string name, int lockedVotes)
    {
        var h = H(name);
        var r = ClassifyPreMiRnaMiPred(h.Sequence)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r.Check.IsHairpinLike, Is.True);
            Assert.That(r.PValue!.Value, Is.LessThanOrEqualTo(0.01), "p-value of a real pre-miRNA");
            Assert.That(r.Features, Has.Length.EqualTo(34));
            Assert.That(r.Features![32], Is.EqualTo(r.FreeEnergy));
            Assert.That(r.Features[33], Is.EqualTo(r.PValue));
            Assert.That(r.RealVotes, Is.EqualTo(lockedVotes));
            Assert.That(Math.Abs(r.RealVotes - h.PythonVotes!.Value), Is.LessThanOrEqualTo(5));
            Assert.That(r.RealProbability, Is.EqualTo(r.RealVotes / 1000.0));
            Assert.That(r.TreeCount, Is.EqualTo(1000));
            Assert.That(r.IsPreMiRna, Is.True);
        });
    }

    // MP9 — pseudo hairpins (shuffles of let-7a-1 passing the pre-judgement): voted pseudo, p-value consistent with the numpy/Perl-shuffle p
    // (two independent shuffle algorithms and RNGs: |Δp| within 4 sampling sigma of the difference). s7/s9 are shuffles whose canonical fold
    // agrees with ViennaRNA (same structure, ΔG° within 0.5 kcal/mol); neg-let7a-s3 is excluded here because the canonical engine scores
    // its 1-nt bulge (30,53)->(32,52) 1.26 kcal/mol lower than ViennaRNA (cross-batch note, RnaSecondaryStructure), which moves its p-value.
    [TestCase("neg-let7a-s7")]
    [TestCase("neg-let7a-s9")]
    public void MP9_PseudoHairpins_ClassifiedPseudo_PValueConsistent(string name)
    {
        var h = H(name);
        var r = ClassifyPreMiRnaMiPred(h.Sequence)!.Value;
        Assert.That(r.Check.IsHairpinLike, Is.True);
        double p = r.PValue!.Value, q = h.PythonPValue!.Value;
        double sigma = Math.Sqrt(((p * (1 - p)) + (q * (1 - q))) / 999.0) + 0.01; // +0.01: fold-energy differences between the engines
        Assert.That(Math.Abs(p - q), Is.LessThanOrEqualTo(4 * sigma), $"C# p={p}, Python p={q}");
        Assert.That(r.IsPreMiRna, Is.False);
        Assert.That(r.RealVotes, Is.LessThan(500));
        Assert.That(Math.Abs(r.RealVotes - h.PythonVotes!.Value), Is.LessThanOrEqualTo(60));
    }

    // MP10 — not hairpin-like: the package does not classify (no shuffling, no features).
    [TestCase("pri-mir21")]
    [TestCase("pri-let7a1")]
    public void MP10_NotHairpinLike_IsNotClassified(string name)
    {
        var r = ClassifyPreMiRnaMiPred(H(name).Sequence)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r.Check.IsHairpinLike, Is.False);
            Assert.That(r.PValue, Is.Null);
            Assert.That(r.Features, Is.Null);
            Assert.That(r.RealProbability, Is.Null);
            Assert.That(r.RealVotes, Is.EqualTo(0));
            Assert.That(r.IsPreMiRna, Is.False);
        });
    }

    // MP11 — a long pri-miRNA with flanks that still passes the pre-judgement (hsa-mir-155, 125 nt): classified real like the precursor.
    [Test]
    public void MP11_PriMiR155_PassesPrejudgementAndVotesReal()
    {
        var r = ClassifyPreMiRnaMiPred(H("pri-mir155").Sequence)!.Value;
        Assert.That(r.Check.IsHairpinLike, Is.True);
        Assert.That(r.IsPreMiRna, Is.True);
        Assert.That(Math.Abs(r.RealVotes - H("pri-mir155").PythonVotes!.Value), Is.LessThanOrEqualTo(15));
    }

    // MP12 — same seed => identical result; the default RNG is Random(MiPredDefaultSeed); other seeds change the p-value feature only within sampling error.
    [Test]
    public void MP12_Determinism_And_Seeding()
    {
        string seq = H("hsa-mir-155").Sequence;
        var a = ClassifyPreMiRnaMiPred(seq, new Random(MiPredDefaultSeed), 99)!.Value;
        var b = ClassifyPreMiRnaMiPred(seq, null, 99)!.Value;
        var c = ClassifyPreMiRnaMiPred(seq, new Random(12345), 99)!.Value;
        Assert.That(b.PValue, Is.EqualTo(a.PValue));
        Assert.That(b.RealVotes, Is.EqualTo(a.RealVotes));
        Assert.That(b.Features, Is.EqualTo(a.Features));
        Assert.That(c.PValue, Is.LessThanOrEqualTo(0.05));
        Assert.That(a.PValue, Is.EqualTo(0.01).Within(1e-12), "(1 + 0)/(99 + 1) for a precursor whose shuffles never fold as low");
    }

    // MP13 — clean-up of the package (FASTA): other characters dropped, DNA read as RNA, case-insensitive; nothing left => null.
    [Test]
    public void MP13_QueryCleanUp()
    {
        string seq = H("hsa-mir-155").Sequence;
        string messy = seq.ToLowerInvariant().Replace('u', 't').Insert(10, "- 12\n").Insert(3, "NN");
        var a = ClassifyPreMiRnaMiPred(messy, new Random(3), 49)!.Value;
        var b = ClassifyPreMiRnaMiPred(seq, new Random(3), 49)!.Value;
        Assert.That(a.Sequence, Is.EqualTo(seq));
        Assert.That(a.DotBracket, Is.EqualTo(b.DotBracket));
        Assert.That(a.RealVotes, Is.EqualTo(b.RealVotes));
        Assert.That(ClassifyPreMiRnaMiPred(""), Is.Null);
        Assert.That(ClassifyPreMiRnaMiPred("NNNN--"), Is.Null);
    }

    // MP14 — argument contract.
    [Test]
    public void MP14_ArgumentValidation()
    {
        Assert.Throws<ArgumentNullException>(() => ClassifyPreMiRnaMiPred(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ClassifyPreMiRnaMiPred("GGGGAAAACCCC", null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => ClassifyPreMiRnaMiPred("GGGGAAAACCCC", null, 10, -1));
        Assert.Throws<ArgumentNullException>(() => AssessMiPredHairpin(null!, -20));
        Assert.Throws<ArgumentNullException>(() => MiPredVotes(null!));
        Assert.Throws<ArgumentException>(() => MiPredVotes(new double[33]));
        Assert.Throws<ArgumentException>(() => ComputeMiPredTripletFrequencies("GGGAAACCC", "(((...))"));
        Assert.Throws<ArgumentException>(() => ComputeMiPredTripletFrequencies("GGGAAACCC", "(((...)]]"));
        Assert.Throws<ArgumentException>(() => ComputeMiPredTripletFrequencies("GGGNAACCC", "(((...)))"));
        Assert.That(ComputeMiPredFeatures("GGGAAACCC", "(((...)))", -3.1, 0.5), Has.Length.EqualTo(34));
    }

    // MP15 — the in-repo logistic model and Triplet-SVM stay untouched by the opt-in MiPred (same call, same values as before).
    [Test]
    public void MP15_OtherClassifiersUnchanged()
    {
        var f = ExtractPreMiRnaFeatures(H("hsa-mir-21").Sequence)!.Value;
        Assert.That(f.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
        Assert.That(f.Mfei, Is.EqualTo(1.0037142857142858).Within(1e-12));
    }
}
