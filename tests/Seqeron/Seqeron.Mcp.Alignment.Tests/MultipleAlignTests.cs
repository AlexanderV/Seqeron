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

    // B09 F63 (A1-10c b): multiple_align_clustalw type = protein delegates to the PROTEIN overload of SequenceAligner.MultipleAlignClustalW.
    // Rows of the default run = the installed clustalw 2.1 (-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT), set 20 of clustalw_protein_msa_oracle.tsv.
    private static readonly string[] ClustalWProteinSet =
        ["LFFFFFFFMFFFFFDKFFFFFFFFFFFFFFFF", "FLAFLVKTYAHTGLKDKEGF", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", "FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF"];

    private static readonly string[] ClustalWProteinRows =
        ["LFFF---FFFFMFFFFFDKFFFFFFFFFFFFFFFF", "------------FLAFLVKTYAHTGLKDKEGF---", "-FFF---FFFFFFFFFFFFFFFFFFFFFFFFFFFF", "-FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF-"];

    [Test]
    public void MultipleAlignClustalW_Protein_DelegatesToTheProteinClustalWPort_AndParsesEveryOption()
    {
        var def = AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein");
        var opt = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinOptions
        {
            Matrix = global::Seqeron.Genomics.Infrastructure.ClustalWProfileMatrix.Pam, GapOpen = 12, GapExtend = 0.4,
            PairwiseMode = global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMode.Quick,
            Pairwise = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions
            {
                Matrix = global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Blosum, GapOpen = 8, GapExtend = 0.3,
                Quick = new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions { KTuple = 2, Window = 6, TopDiagonals = 7, PairGap = 5 },
            },
            DivergenceCutoff = 50, NoWeights = true, NegativeMatrix = true, NoPreferredGapPenalties = true, NoHydrophilicGapPenalties = true,
            VariableRegionGapPenalties = true, GapSeparationDistance = 6, HydrophilicResidues = "ACDEFGHIK",
        };
        var scoring = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Pam250, -5.0, -0.5);
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(ClustalWProteinSet, opt, scoring);
        var viaTool = AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, 12, 0.4, "clustalw-quick", 8, 0.3, 0.5f, 50, true, true, 2, 6, 7, 5, 1, -1, -5.0, -0.5,
            "protein", "pam", "blosum", true, true, true, 6, "ACDEFGHIK", "PAM250");
        Assert.Multiple(() =>
        {
            Assert.That(def.AlignedSequences, Is.EqualTo(ClustalWProteinRows));
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(viaTool.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(ClustalWProteinSet));
            Assert.That(AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", matrix: "gonnet", pwMatrix: "gonnet", distance: "full").AlignedSequences, Is.EqualTo(def.AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignClustalW(Array.Empty<string>(), type: "protein").AlignedSequences, Is.Empty);
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "rna"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", matrix: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", pwMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", scoreMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(new[] { "MKVL", "MKJL" }, type: "protein"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", gapOpen: -1));
            // the DNA defaults are unchanged by the protein parameters
            Assert.That(AlignmentTools.MultipleAlignClustalW(ClustalWSet, matrix: "pam", noPreferredGaps: true).AlignedSequences, Is.EqualTo(AlignmentTools.MultipleAlignClustalW(ClustalWSet).AlignedSequences));
        });
    }

    // B09 F67 (A1-10d part 2): multiple_align_clustalw useEndGaps / endGapPenalties. Rows = the patched clustalw 2.1 of the Docker package (clustalw_end_gaps_oracle.tsv, set D84 DNA):
    // control, useendgaps, endgappen, both; set D17 aborts with endGapPenalties (no output file) -> an explicit error.
    private static readonly string[] EndGapSet = ["AAGAGCCCT", "CGATAAAGCT", "CGATAAGTT", "CGATAAGTT", "CGCATAAGTT", "CGTATATTC", "CGAAATCT"];

    [Test]
    public void MultipleAlignClustalW_EndGapOptions_DelegateAndAbortIsExplicit()
    {
        var control = AlignmentTools.MultipleAlignClustalW(EndGapSet);
        var pen = AlignmentTools.MultipleAlignClustalW(EndGapSet, useEndGaps: false, endGapPenalties: true);
        var both = AlignmentTools.MultipleAlignClustalW(EndGapSet, useEndGaps: true, endGapPenalties: true);
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(EndGapSet,
            new global::Seqeron.Genomics.Infrastructure.ClustalWOptions { EndGapPenalties = true });
        Assert.Multiple(() =>
        {
            Assert.That(control.AlignedSequences, Is.EqualTo(new[] { "---AAGAGCCCT", "CGATAAAGCT--", "CG-ATAAGTT--", "CG-ATAAGTT--", "CGCATAAGTT--", "CGTATATTC---", "--CGAAATCT--" }));
            Assert.That(pen.AlignedSequences, Is.EqualTo(new[] { "AAGA-GCCCT", "CGATAAAGCT", "CG-ATAAGTT", "CG-ATAAGTT", "CGCATAAGTT", "CGTATAT-TC", "CGA--AATCT" }));
            Assert.That(pen.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(both.AlignedSequences, Is.EqualTo(pen.AlignedSequences));
            var ex = Assert.Throws<global::Seqeron.Genomics.Infrastructure.ClustalWAbortedException>(() => AlignmentTools.MultipleAlignClustalW(new[] { "AT", "AAACA", "AA" }, endGapPenalties: true))!;
            Assert.That(ex.Message, Does.Contain("abandons the alignment"));
            Assert.That(AlignmentTools.MultipleAlignClustalW(new[] { "AT", "AAACA", "AA" }, useEndGaps: true).AlignedSequences, Has.Length.EqualTo(3));
        });
    }

    // B09 F68 (A1-10d part 3): multiple_align_clustalw dnaMatrixFile / pwDnaMatrixFile / dnaMatrix / pwDnaMatrix / userDnaMatrixUse. Rows = the installed clustalw 2.1 driven through its
    // interactive menu (the file is used; set 244 of clustalw_dna_matrix_oracle.tsv, pwgapopen 10); the command line reads but ignores the file (userDnaMatrixUse = commandline).
    private static readonly string[] DnaMatrixSet = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "GAACTCGGAATGGGAGTGG", "GGACTCGGGATGGGAGTGG"];
    private static readonly string[] DnaMatrixRows = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "-GAACTCGGAATGGGAGTGG", "-GGACTCGGGATGGGAGTGG"];
    private const string DnaMatrixFile = "# DNA matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n A T G C S W R Y K M B V H D N U *\nA 1 -9 1 -9 -9 -1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -4\nT -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\nG 1 -9 1 -9 -1 -9 -1 -9 -1 -9 -1 -1 -9 -1 -1 -9 -4\nC -9 1 -9 1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -1 1 -4\nS -9 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nW -1 -1 -9 -9 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nR -1 -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -9 -4\nY -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nK -9 -1 -1 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -4\nM -1 -9 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -9 -4\nB -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nV -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nH -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nD -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nN -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nU -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\n* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4\n";

    [Test]
    public void MultipleAlignClustalW_DnaMatrices_DelegateAndFollowClustalWsCommandLineAndMenu()
    {
        var applied = AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, dnaMatrixFile: DnaMatrixFile);
        var cli = AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, dnaMatrixFile: DnaMatrixFile, userDnaMatrixUse: "commandline");
        var plain = AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10);
        var clustalw = AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrix: "CLUSTALW", pwDnaMatrix: "clustalw");
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(DnaMatrixSet, new global::Seqeron.Genomics.Infrastructure.ClustalWOptions
        {
            DnaMatrix = global::Seqeron.Genomics.Infrastructure.ClustalWDnaMatrix.Clustalw,
            PairwiseDnaMatrix = global::Seqeron.Genomics.Infrastructure.ClustalWDnaMatrix.Clustalw,
        });
        Assert.Multiple(() =>
        {
            Assert.That(applied.AlignedSequences, Is.EqualTo(DnaMatrixRows));
            Assert.That(cli.AlignedSequences, Is.EqualTo(plain.AlignedSequences));
            Assert.That(applied.AlignedSequences, Is.Not.EqualTo(plain.AlignedSequences));
            Assert.That(clustalw.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, pwDnaMatrixFile: DnaMatrixFile).AlignedSequences, Has.Length.EqualTo(4));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwDnaMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, userDnaMatrixUse: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrixFile: "CLUSTAL_SERIES\nMATRIX 0 100 x\n"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwDnaMatrixFile: "# nothing\n"));
        });
    }

    // B09 F64 (A1-10d file part): multiple_align_clustalw matrixFile / matrixFiles / pwMatrixFile = clustalw -MATRIX=file (here a CLUSTAL_SERIES file with an unlabelled
    // tab-separated integer matrix and a labelled floating-point matrix) / -PWMATRIX=file. Rows = the installed clustalw 2.1, set 11 of clustalw_user_matrix_oracle.tsv.
    private static readonly string[] UserMatrixSet = ["ALAAALLLAALAAA", "RPKRQESPKNPGRQEPKEPENSNRRQKK", "ALAAAECLAALAAA", "AGAAALLSEPLAALAAA", "ALAAAPLLAALAAA", "ALGGALLNAALAAA"];

    private static readonly string[] UserMatrixRows = ["------------------------------ALAAALLLAALAAA", "-RPKRQESPKNPGRQEPKEPENSNRRQKK---------------", "ALAAAEC--LAALAAA----------------------------", "----AGAAALLSEPLAALAAA-----------------------", "----------------ALAAAPLLAALAAA--------------", "---ALGG--ALLNAALAAA-------------------------"];

    private const string UserSeries = "# series of matrices\nCLUSTAL_SERIES\nMATRIX 0 23 m0.mat\nMATRIX 23 100 m1.mat\n# comment between entries\n";
    private const string UserMatrixM0 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n\tF\tA\tQ\tK\tC\tI\tY\tW\tH\tD\tT\tE\tL\tN\t-\n1\t-1\t3\t-6\t5\t7\t-5\t-5\t-3\t-8\t3\t-8\t5\t3\t-4\n0\t0\t-1\t-9\t-3\t1\t0\t1\t-5\t8\t0\t-1\t-3\t-1\t-4\n0\t0\t0\t3\t5\t-9\t-2\t-8\t9\t3\t-9\t-7\t-3\t-9\t-4\n-1\t0\t0\t0\t5\t-7\t-9\t2\t1\t2\t5\t0\t-9\t1\t-4\n0\t0\t-1\t-1\t1\t-5\t7\t-3\t0\t-4\t8\t6\t-4\t4\t-4\n0\t0\t0\t0\t0\t0\t-1\t0\t-3\t3\t-2\t-5\t-1\t0\t-4\n1\t0\t0\t0\t0\t0\t1\t-5\t-1\t8\t-3\t-3\t-1\t9\t-4\n0\t-1\t0\t0\t-1\t-1\t0\t2\t6\t9\t-7\t5\t7\t4\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t1\t-8\t-4\t-3\t-6\t-9\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t6\t-9\t6\t-4\t-4\n0\t0\t0\t0\t0\t0\t0\t-1\t0\t0\t0\t-7\t-2\t-5\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t0\t0\t1\t-6\t-4\n0\t0\t0\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t1\t5\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t-4\n-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\n";
    private const string UserMatrixM1 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\nR M E G T Q W N *\nR 2 -0.5 -1 2.5 1 4 2 -2.5 -2\nM 0.5 3.5 2.5 -2.5 0.5 -3 -1 -4.5 -2\nE -2 -2 2.5 -3 -4.5 -2 -2.5 -2.5 -2\nG -1 -1.5 0.5 4 1 -3.5 -2.5 4.5 -2\nT -1 -0.5 0 0.5 4 -4 1 1.5 -2\nQ 0.5 -1.5 1 -2.5 0.5 2 3.5 4.5 -2\nW 1.5 0 -1.5 0 0 1.5 4 -0.5 -2\nN 1 1.5 0 1.5 -0.5 1 -2.5 3 -2\n* -2 -2 -2 -2 -2 -2 -2 -2 -2\n";

    [Test]
    public void MultipleAlignClustalW_Protein_UserMatrixFiles_EqualClustalW_AndAreValidated()
    {
        var files = new Dictionary<string, string> { ["m0.mat"] = UserMatrixM0, ["m1.mat"] = UserMatrixM1 };
        var viaTool = AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserSeries, matrixFiles: files);
        var series = global::Seqeron.Genomics.Infrastructure.ClustalWMatrixFile.Parse(UserSeries, n => files[n]);
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(UserMatrixSet,
            new global::Seqeron.Genomics.Infrastructure.ClustalWProteinOptions { UserMatrixSeries = series.Series });
        // a single matrix file (-MATRIX=file) and the pairwise matrix file (-PWMATRIX=file) reach the same code
        var viaSingle = AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserMatrixM1, pwMatrixFile: UserMatrixM0);
        var viaSingleApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(UserMatrixSet,
            new global::Seqeron.Genomics.Infrastructure.ClustalWProteinOptions
            {
                UserMatrix = global::Seqeron.Genomics.Infrastructure.ClustalWUserMatrix.Parse(UserMatrixM1),
                Pairwise = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions { UserMatrix = global::Seqeron.Genomics.Infrastructure.ClustalWUserMatrix.Parse(UserMatrixM0) },
            });
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(UserMatrixRows));
            Assert.That(viaApi.AlignedSequences, Is.EqualTo(UserMatrixRows));
            Assert.That(viaSingle.AlignedSequences, Is.EqualTo(viaSingleApi.AlignedSequences));
            // the user matrices change the result against the built-in series
            Assert.That(AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein").AlignedSequences, Is.Not.EqualTo(UserMatrixRows));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserSeries), "the series names files that are not supplied");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: "# nothing\n"), "no header");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", pwMatrixFile: UserSeries), "a series where one matrix is expected");
        });
    }

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
            Assert.That(AlignmentTools.MultipleAlignMuscle(Array.Empty<string>()).AlignedSequences, Is.Empty);
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
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("bogus", aligned));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("profile", aligned));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MuscleMsaCommand("refine", Array.Empty<string>()));
        });
    }

    // B09 F66 (A1-13b): multiple_align_tcoffee delegates to SequenceAligner.MultipleAlignTCoffeeProgram (the T-Coffee 13.41 program: its own libraries, guide tree, progressive stage).
    private static readonly string[] TCoffeeProgramSet = ["ACGTTGCATGCCGATAGC", "ACGTGCATGCGATAGGC", "ACTTGCATCCGATTGC", "CGTTGCATGGCCGATAGC"];

    [Test]
    public void MultipleAlignTCoffee_DelegatesToTheProgramPort_ParsesEveryOption_AndValidates()
    {
        // rows written by t_coffee -seq s.fa -method slow_pair,lalign_id_pair -output fasta_aln (installed 13.41)
        string[] program = ["ACGTTGCAT-GCCGATAG-C", "ACGT-GCAT-GC-GATAGGC", "AC-TTGCAT--CCGATTG-C", "-CGTTGCATGGCCGATAG-C"];
        var viaTool = AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["slow_pair", "LALIGN_ID_PAIR"]);
        var viaApi = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignTCoffeeProgram(TCoffeeProgramSet, new global::Seqeron.Genomics.Alignment.TCoffeeProgramOptions
        {
            Methods = [global::Seqeron.Genomics.Alignment.TCoffeePairMethod.SlowPair, global::Seqeron.Genomics.Alignment.TCoffeePairMethod.LalignIdPair],
        });
        var plain = AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet);
        var relax0 = AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["proba_pair"], 3, false, "dna");
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(program));
            Assert.That(viaApi.AlignedSequences, Is.EqualTo(program));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(plain.AlignedSequences, Is.EqualTo(AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["proba_pair"]).AlignedSequences), "the default method is proba_pair");
            Assert.That(plain.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(TCoffeeProgramSet));
            Assert.That(relax0.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(TCoffeeProgramSet));
            Assert.That(AlignmentTools.MultipleAlignTCoffee(["MKVLAAGIVG", "MKVLAGIVGL", "MRVLAAGIIG"], null, 10, true, "protein").AlignedSequences.Select(r => r.Replace("-", "")),
                Is.EqualTo(new[] { "MKVLAAGIVG", "MKVLAGIVGL", "MRVLAAGIIG" }));
            Assert.That(AlignmentTools.MultipleAlignTCoffee(Array.Empty<string>()).AlignedSequences, Is.Empty);
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["blast_pair"]));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, type: "rna"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, lalignTopAlignments: 0));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignTCoffee(new[] { "ACGT", "" }));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignTCoffee(new[] { "ACGUACGU", "ACGUUCGU" }), "proba_pair of RNA is not available");
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

    // B09 F62 (A1-15b b): distance = clustalw | clustalw-quick with type = protein uses ClustalW's PROTEIN pairwise stages; pwMatrix / pwGapOpen / pwGapExtend.
    [Test]
    public void MultipleAlign_Protein_ClustalWDistance_DelegatesWithTheProteinPairwiseOptions()
    {
        var blosum = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var affine = global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine;
        var full = new global::Seqeron.Genomics.Infrastructure.MsaOptions
        {
            GapModel = affine,
            Distance = global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.ClustalWAffine,
            ClustalWProtein = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions
            {
                Matrix = global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Pam, GapOpen = 5, GapExtend = 0.5,
            },
        };
        var quick = new global::Seqeron.Genomics.Infrastructure.MsaOptions
        {
            GapModel = affine,
            Distance = global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.ClustalWQuick,
            ClustalWProtein = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions
            {
                Quick = new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions { KTuple = 2, Window = 6, TopDiagonals = 7, PairGap = 5 },
            },
        };
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw", "affine", "upgma", null, null, null, null, "protein", null, "pam", 5, 0.5).AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignProgressive(Proteins, blosum, full).AlignedSequences));
            Assert.That(AlignmentTools.MultipleAlignIterative(Proteins, 1, -1, -9.5, -0.5, 16, "clustalw", "affine", "nj", null, null, null, null, "protein", null, "pam", 5, 0.5).AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignIterative(Proteins, blosum, full with { GuideTree = global::Seqeron.Genomics.Infrastructure.MsaGuideTree.NeighborJoining }).AlignedSequences));
            // quick: ktuple 2 given -> window / topDiagonals / pairGap explicit
            Assert.That(AlignmentTools.MultipleAlignConsistency(Proteins, 1, -1, -9.5, -0.5, "clustalw-quick", "affine", "upgma", 2, 6, 7, 5, "paper", "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignConsistency(Proteins, blosum, quick).AlignedSequences));
            // defaults: ClustalW's protein defaults (GONNET, 10 / 0.1; quick 1 / 5 / 5 / 3)
            Assert.That(AlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw-quick", "affine", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignProgressive(Proteins, blosum,
                    new global::Seqeron.Genomics.Infrastructure.MsaOptions { GapModel = affine, Distance = global::Seqeron.Genomics.Infrastructure.MsaDistanceMode.ClustalWQuick }).AlignedSequences));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw", "affine", "upgma", null, null, null, null, "protein", null, "bogus"));
        });
    }
}

