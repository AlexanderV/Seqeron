using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class MultipleAlignTests
{

    // B09 F51 (A1-15): guide-tree distance source. The set below is discriminating: ClustalW 2.1's pairwise stage and the
    // repository's linear NW alignment give different distances, hence different guide trees and different rows
    // (rows of the independent Python implementation in scripts/alignment/generate_clustalw_pairwise_oracle.py).
    private static readonly string[] DistanceSet = ["CCAACCAACCAA", "ACACTCCGCCCGCACTCAGAACTTT", "TAA"];

    // B09 F56 (A1-10b): multiple_align_clustalw delegates to ClustalWAligner.MultipleAlignClustalW (ClustalW 2.1 multiple alignment).
    // Rows of the set below = the installed clustalw 2.1 (-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT).
    private static readonly string[] ClustalWSet =
        ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT", "ACGTAGGATTGACCGTAGCCATCGATGGA", "GGGTTTACCAGTAGGCATCGATGGAAACCT"];

    // B09 F63 (A1-10c b): multiple_align_clustalw type = protein delegates to the PROTEIN overload of ClustalWAligner.MultipleAlignClustalW.
    // Rows of the default run = the installed clustalw 2.1 (-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT), set 20 of clustalw_protein_msa_oracle.tsv.
    private static readonly string[] ClustalWProteinSet =
        ["LFFFFFFFMFFFFFDKFFFFFFFFFFFFFFFF", "FLAFLVKTYAHTGLKDKEGF", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", "FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF"];

    private static readonly string[] ClustalWProteinRows =
        ["LFFF---FFFFMFFFFFDKFFFFFFFFFFFFFFFF", "------------FLAFLVKTYAHTGLKDKEGF---", "-FFF---FFFFFFFFFFFFFFFFFFFFFFFFFFFF", "-FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF-"];

    // B09 F67 (A1-10d part 2): multiple_align_clustalw useEndGaps / endGapPenalties. Rows = the patched clustalw 2.1 of the Docker package (clustalw_end_gaps_oracle.tsv, set D84 DNA):
    // control, useendgaps, endgappen, both; set D17 aborts with endGapPenalties (no output file) -> an explicit error.
    private static readonly string[] EndGapSet = ["AAGAGCCCT", "CGATAAAGCT", "CGATAAGTT", "CGATAAGTT", "CGCATAAGTT", "CGTATATTC", "CGAAATCT"];

    // B09 F68 (A1-10d part 3): multiple_align_clustalw dnaMatrixFile / pwDnaMatrixFile / dnaMatrix / pwDnaMatrix / userDnaMatrixUse. Rows = the installed clustalw 2.1 driven through its
    // interactive menu (the file is used; set 244 of clustalw_dna_matrix_oracle.tsv, pwgapopen 10); the command line reads but ignores the file (userDnaMatrixUse = commandline).
    private static readonly string[] DnaMatrixSet = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "GAACTCGGAATGGGAGTGG", "GGACTCGGGATGGGAGTGG"];
    private static readonly string[] DnaMatrixRows = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "-GAACTCGGAATGGGAGTGG", "-GGACTCGGGATGGGAGTGG"];
    private const string DnaMatrixFile = "# DNA matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n A T G C S W R Y K M B V H D N U *\nA 1 -9 1 -9 -9 -1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -4\nT -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\nG 1 -9 1 -9 -1 -9 -1 -9 -1 -9 -1 -1 -9 -1 -1 -9 -4\nC -9 1 -9 1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -1 1 -4\nS -9 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nW -1 -1 -9 -9 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nR -1 -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -9 -4\nY -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nK -9 -1 -1 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -4\nM -1 -9 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -9 -4\nB -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nV -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nH -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nD -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nN -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nU -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\n* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4\n";

    // B09 F64 (A1-10d file part): multiple_align_clustalw matrixFile / matrixFiles / pwMatrixFile = clustalw -MATRIX=file (here a CLUSTAL_SERIES file with an unlabelled
    // tab-separated integer matrix and a labelled floating-point matrix) / -PWMATRIX=file. Rows = the installed clustalw 2.1, set 11 of clustalw_user_matrix_oracle.tsv.
    private static readonly string[] UserMatrixSet = ["ALAAALLLAALAAA", "RPKRQESPKNPGRQEPKEPENSNRRQKK", "ALAAAECLAALAAA", "AGAAALLSEPLAALAAA", "ALAAAPLLAALAAA", "ALGGALLNAALAAA"];

    private static readonly string[] UserMatrixRows = ["------------------------------ALAAALLLAALAAA", "-RPKRQESPKNPGRQEPKEPENSNRRQKK---------------", "ALAAAEC--LAALAAA----------------------------", "----AGAAALLSEPLAALAAA-----------------------", "----------------ALAAAPLLAALAAA--------------", "---ALGG--ALLNAALAAA-------------------------"];

    private const string UserSeries = "# series of matrices\nCLUSTAL_SERIES\nMATRIX 0 23 m0.mat\nMATRIX 23 100 m1.mat\n# comment between entries\n";
    private const string UserMatrixM0 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n\tF\tA\tQ\tK\tC\tI\tY\tW\tH\tD\tT\tE\tL\tN\t-\n1\t-1\t3\t-6\t5\t7\t-5\t-5\t-3\t-8\t3\t-8\t5\t3\t-4\n0\t0\t-1\t-9\t-3\t1\t0\t1\t-5\t8\t0\t-1\t-3\t-1\t-4\n0\t0\t0\t3\t5\t-9\t-2\t-8\t9\t3\t-9\t-7\t-3\t-9\t-4\n-1\t0\t0\t0\t5\t-7\t-9\t2\t1\t2\t5\t0\t-9\t1\t-4\n0\t0\t-1\t-1\t1\t-5\t7\t-3\t0\t-4\t8\t6\t-4\t4\t-4\n0\t0\t0\t0\t0\t0\t-1\t0\t-3\t3\t-2\t-5\t-1\t0\t-4\n1\t0\t0\t0\t0\t0\t1\t-5\t-1\t8\t-3\t-3\t-1\t9\t-4\n0\t-1\t0\t0\t-1\t-1\t0\t2\t6\t9\t-7\t5\t7\t4\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t1\t-8\t-4\t-3\t-6\t-9\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t6\t-9\t6\t-4\t-4\n0\t0\t0\t0\t0\t0\t0\t-1\t0\t0\t0\t-7\t-2\t-5\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t0\t0\t1\t-6\t-4\n0\t0\t0\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t1\t5\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t-4\n-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\n";
    private const string UserMatrixM1 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\nR M E G T Q W N *\nR 2 -0.5 -1 2.5 1 4 2 -2.5 -2\nM 0.5 3.5 2.5 -2.5 0.5 -3 -1 -4.5 -2\nE -2 -2 2.5 -3 -4.5 -2 -2.5 -2.5 -2\nG -1 -1.5 0.5 4 1 -3.5 -2.5 4.5 -2\nT -1 -0.5 0 0.5 4 -4 1 1.5 -2\nQ 0.5 -1.5 1 -2.5 0.5 2 3.5 4.5 -2\nW 1.5 0 -1.5 0 0 1.5 4 -0.5 -2\nN 1 1.5 0 1.5 -0.5 1 -2.5 3 -2\n* -2 -2 -2 -2 -2 -2 -2 -2 -2\n";

    // B09 F65 (A1-10c a): multiple_align_muscle delegates to SequenceAligner.MultipleAlignMuscle (MUSCLE 3.8.1551 progressive alignment, maxiters 1 / 2).
    private static readonly string[] MuscleSet =
        ["MKVLAAGIVGLLLAQPWTTQAEAKP", "MKVLAGIVGLLAQPWTQAEAKPLL", "MRVLAAGIIGLLLQPWTTQAEAK", "KVLAAGIVALLLAQPWTTAEAKPS", "MKVAAGIVGLLLAQPWTTQEAKP"];

    [Test]
    public void MultipleAlignMuscle_DelegatesToThePort_ParsesEveryOption_AndValidates()
    {
        var viaTool = AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein");
        var opt = new global::Seqeron.Genomics.Infrastructure.MuscleOptions
        {
            SequenceType = global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Protein,
            MaxIterations = 1,
            ScoreModel = global::Seqeron.Genomics.Infrastructure.MuscleScoreModel.SumOfPairs,
            GapOpen = -1200, GapExtend = -5, Center = 1,
            Cluster1 = global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.Upgma,
            Cluster2 = global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.UpgmaMax,
            Sueff = 0.3f, MaxTrees = 0, TreeOrder = true,
            Distance1 = global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer6_6,
        };
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignMuscle(MuscleSet, opt);
        var viaToolAll = AlignmentTools.MultipleAlignMuscle(MuscleSet, "PROTEIN", 1, "sp", -1200, -5, 1, "kmer6_6", "upgma", "upgmamax", 0.3f, 0, true);
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(MuscleSet));
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignMuscle(MuscleSet, new global::Seqeron.Genomics.Infrastructure.MuscleOptions { SequenceType = global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Protein }).AlignedSequences));
            Assert.That(viaToolAll.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(viaToolAll.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(AlignmentTools.MultipleAlignMuscle(new[] { "ACGTACGTTAGCTAGC", "ACGTACGTAGCTAGC", "ACGTTCGTTAGCTAGC" }).AlignedSequences.Select(r => r.Replace("-", "")),
                Is.EqualTo(new[] { "ACGTACGTTAGCTAGC", "ACGTACGTAGCTAGC", "ACGTTCGTTAGCTAGC" }), "type auto: nucleotides");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(Array.Empty<string>()), "A3-4: the program refuses an empty input");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, maxIters: 0), "at least one iteration");
            Assert.That(AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", 4).AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(MuscleSet), "iterations >= 3");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, "rna2"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, scoreModel: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, distance1: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, cluster1: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, cluster2: "nj"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(new[] { "ACD", "" }));
        });
    }

    // B09 F76 (A1-10c a2): -stable is refused by the program ("not supported in this version of muscle"); the tool rejects it as well.
    [Test]
    public void MultipleAlignMuscle_Stable_IsRejectedLikeTheProgram()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", stable: true));
        Assert.That(AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", stable: false).AlignedSequences, Has.Length.EqualTo(MuscleSet.Length));
    }

    // B09 F70 (A1-10c a2): the options of MUSCLE 3.8 beyond the first Docker round are parsed and delegated to SequenceAligner.MuscleRun / the MSA command lines.
    [Test]
    public void MultipleAlignMuscle_NewOptions_AreParsedAndDelegated()
    {
        var viaTool = AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", 5, distance1: "kmer20_4", cluster1: "neighborjoining", root1: "midlongestspan", distance2: "pctidlog", weight1: "henikoff", weight2: "gsc",
            objScore: "sp", termGaps: "ext", anchors: false, brenner: true, diags1: true, hydro: 9, names: ["a", "b", "c", "d", "e"], cluster2: "upgma");
        var opt = new global::Seqeron.Genomics.Infrastructure.MuscleOptions
        {
            SequenceType = global::Seqeron.Genomics.Infrastructure.MuscleSequenceType.Protein,
            MaxIterations = 5,
            Distance1 = global::Seqeron.Genomics.Infrastructure.MuscleDistance.Kmer20_4,
            Cluster1 = global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.NeighborJoining,
            Root1 = global::Seqeron.Genomics.Infrastructure.MuscleRoot.MidLongestSpan,
            Distance2 = global::Seqeron.Genomics.Infrastructure.MuscleDistance.PctIdLog,
            Cluster2 = global::Seqeron.Genomics.Infrastructure.MuscleClusterMethod.Upgma,
            Weight1 = global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.Henikoff,
            Weight2 = global::Seqeron.Genomics.Infrastructure.MuscleSeqWeight.Gsc,
            ObjectiveScore = global::Seqeron.Genomics.Infrastructure.MuscleObjectiveScore.Sp,
            TermGaps = global::Seqeron.Genomics.Infrastructure.MuscleTermGaps.Ext,
            Anchors = false,
            Brenner = true,
            Diags1 = true,
            Hydro = 9,
            SequenceNames = ["a", "b", "c", "d", "e"],
        };
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MuscleRun(MuscleSet, opt);
        var trees = AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", 3, names: ["a", "b", "c", "d", "e"]);
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(viaApi.Alignment.AlignedSequences));
            Assert.That(viaTool.Tree1, Is.EqualTo(viaApi.Tree1));
            Assert.That(viaTool.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(MuscleSet));
            Assert.That(trees.Tree1, Does.Contain("a:").And.StartWith("(\n"));
            Assert.That(trees.Tree2, Is.Not.Null);
            Assert.That(trees.DistanceMatrix1, Does.StartWith("5\n"));
            Assert.That(AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", computeWeights: true).WeightsFile, Is.Not.Null);
            Assert.That(AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein", clusterOnly: true).AlignedSequences, Is.Empty);
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, weight1: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, objScore: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, root1: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMuscle(MuscleSet, termGaps: "bogus"));
        });
    }

    [Test]
    public void MuscleMsaCommand_RunsEveryMode_AndRejectsProteinsWhereTheProgramCrashes()
    {
        var aligned = AlignmentTools.MultipleAlignMuscle(MuscleSet, "protein").AlignedSequences;
        string[] dna = ["ACGTACGTAGCTAGCT", "ACGTACGAAGCTAGCT"], dna2 = ["ACGTACGTAGCTAGCT", "ACGAACGTAGCTAGCT"];
        Assert.Multiple(() =>
        {
            var profile = AlignmentTools.MuscleMsaCommand("profile", aligned.Take(3).ToArray(), aligned.Skip(3).ToArray(), "protein");
            Assert.That(profile.AlignedSequences, Has.Length.EqualTo(5));
            Assert.That(AlignmentTools.MuscleMsaCommand("refine", aligned, type: "protein").AlignedSequences, Has.Length.EqualTo(5));
            Assert.That(AlignmentTools.MuscleMsaCommand("spscore", aligned, type: "protein").Score, Is.Not.Null);
            Assert.That(AlignmentTools.MuscleMsaCommand("maketree", aligned, type: "protein").Tree, Does.StartWith("(\n"));
            Assert.That(AlignmentTools.MuscleMsaCommand("ppscore", dna, dna2, "dna").Score, Is.Not.Null);
            Assert.That(AlignmentTools.MuscleMsaCommand("sw", ["ACGTACGTACGTTTGACCA"], ["GGGGTACGTACGTAAAA"], "dna").AlignedSequences, Has.Length.EqualTo(2));
            Assert.That(AlignmentTools.MuscleMsaCommand("refinew", dna, type: "dna", refineWindow: 8).AlignedSequences, Has.Length.EqualTo(2));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("ppscore", aligned, aligned, "protein"), "the program crashes: no protein score matrix in this mode");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("profdb", aligned, aligned, "protein"), "the program crashes in -profdb: the alphabet is never initialised");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("bogus", aligned));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("profile", aligned));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("refine", Array.Empty<string>()));
        });
    }

    // B09 F66 (A1-13b): multiple_align_tcoffee delegates to TCoffeeAligner.MultipleAlignTCoffeeProgram (the T-Coffee 13.41 program: its own libraries, guide tree, progressive stage).
    private static readonly string[] TCoffeeProgramSet = ["ACGTTGCATGCCGATAGC", "ACGTGCATGCGATAGGC", "ACTTGCATCCGATTGC", "CGTTGCATGGCCGATAGC"];

    private const string BlastTwoSequenceXml = """
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.12.0+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>s.fa</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>q000</BlastOutput_query-def>
  <BlastOutput_query-len>50</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>L;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
<BlastOutput_iterations>
<Iteration>
  <Iteration_iter-num>1</Iteration_iter-num>
  <Iteration_query-ID>Query_1</Iteration_query-ID>
  <Iteration_query-def>q000</Iteration_query-def>
  <Iteration_query-len>50</Iteration_query-len>
<Iteration_hits>
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
  <Hit_def>q000</Hit_def>
  <Hit_accession>0</Hit_accession>
  <Hit_len>50</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>109.768</Hsp_bit-score>
      <Hsp_score>273</Hsp_score>
      <Hsp_evalue>1.76494e-39</Hsp_evalue>
      <Hsp_query-from>1</Hsp_query-from>
      <Hsp_query-to>50</Hsp_query-to>
      <Hsp_hit-from>1</Hsp_hit-from>
      <Hsp_hit-to>50</Hsp_hit-to>
      <Hsp_query-frame>0</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>50</Hsp_identity>
      <Hsp_positive>50</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>50</Hsp_align-len>
      <Hsp_qseq>QMPLHKEGVTYTDFMQAPDRKNYNQRFHTNPEHIGKAAIWDSMEDNTQLK</Hsp_qseq>
      <Hsp_hseq>QMPLHKEGVTYTDFMQAPDRKNYNQRFHTNPEHIGKAAIWDSMEDNTQLK</Hsp_hseq>
      <Hsp_midline>QMPLHKEGVTYTDFMQAPDRKNYNQRFHTNPEHIGKAAIWDSMEDNTQLK</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
</Iteration_hits>
  <Iteration_stat>
    <Statistics>
      <Statistics_db-num>2</Statistics_db-num>
      <Statistics_db-len>84</Statistics_db-len>
      <Statistics_hsp-len>4</Statistics_hsp-len>
      <Statistics_eff-space>3496</Statistics_eff-space>
      <Statistics_kappa>0.041</Statistics_kappa>
      <Statistics_lambda>0.267</Statistics_lambda>
      <Statistics_entropy>0.14</Statistics_entropy>
    </Statistics>
  </Iteration_stat>
</Iteration>
<Iteration>
  <Iteration_iter-num>2</Iteration_iter-num>
  <Iteration_query-ID>Query_2</Iteration_query-ID>
  <Iteration_query-def>q001</Iteration_query-def>
  <Iteration_query-len>34</Iteration_query-len>
<Iteration_hits>
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gnl|BL_ORD_ID|1</Hit_id>
  <Hit_def>q001</Hit_def>
  <Hit_accession>1</Hit_accession>
  <Hit_len>34</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>73.9442</Hsp_bit-score>
      <Hsp_score>180</Hsp_score>
      <Hsp_evalue>8.25478e-26</Hsp_evalue>
      <Hsp_query-from>1</Hsp_query-from>
      <Hsp_query-to>34</Hsp_query-to>
      <Hsp_hit-from>1</Hsp_hit-from>
      <Hsp_hit-to>34</Hsp_hit-to>
      <Hsp_query-frame>0</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>34</Hsp_identity>
      <Hsp_positive>34</Hsp_positive>
      <Hsp_gaps>0</Hsp_gaps>
      <Hsp_align-len>34</Hsp_align-len>
      <Hsp_qseq>WEFDMPGSGRKNLESPVYYDAKRIGVAFVGHFIY</Hsp_qseq>
      <Hsp_hseq>WEFDMPGSGRKNLESPVYYDAKRIGVAFVGHFIY</Hsp_hseq>
      <Hsp_midline>WEFDMPGSGRKNLESPVYYDAKRIGVAFVGHFIY</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
</Iteration_hits>
  <Iteration_stat>
    <Statistics>
      <Statistics_db-num>2</Statistics_db-num>
      <Statistics_db-len>84</Statistics_db-len>
      <Statistics_hsp-len>2</Statistics_hsp-len>
      <Statistics_eff-space>2560</Statistics_eff-space>
      <Statistics_kappa>0.041</Statistics_kappa>
      <Statistics_lambda>0.267</Statistics_lambda>
      <Statistics_entropy>0.14</Statistics_entropy>
    </Statistics>
  </Iteration_stat>
</Iteration>
</BlastOutput_iterations>
</BlastOutput>

""";

    // B09 F52 (A1-10): gapModel = affine honours gapOpen. Identical pair 33 matches = 66; the two pairs with the 23-nt sequence have
    // 23 matches (46) and ONE 10-nt gap run: affine -8 + 10 * -1 = -18, linear -10.
    private static readonly string[] IndelSet =
        ["ACGTACGTACGTAAAAAAAAAAGGCATCGATCG", "ACGTACGTACGTGGCATCGATCG", "ACGTACGTACGTAAAAAAAAAAGGCATCGATCG"];

    // B09 F53 (A1-12): guideTree = upgma | nj. The set is discriminating: UPGMA and ClustalW-style NJ guide trees give different
    // progressive rows (library call SequenceAligner.MultipleAlignProgressive(..., new MsaOptions { GuideTree = ... }), defaults 1/-1/-2/-1).
    // F56 (A1-10b): the NJ tree now uses ClustalW's own NJ arithmetic (exact ties like ClustalW); the former set became a tie where
    // both trees agree, so the set was replaced by one found by a random search (NJ and UPGMA topologies differ).
    private static readonly string[] GuideTreeSet =
        ["AAAAGTCAGGCAT", "AACAGTTAGGCAT", "AAAGTTAGGCCA", "AAATGTTAGGCAT", "AAAGTTAGGCAA"];

    // B09 F61 (A1-14): type = 'protein' / matrix on every multiple_align* tool; the tools delegate to the SequenceAligner matrix overloads.
    private static readonly string[] Proteins =
    [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
    ];

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
}

