// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

[TestFixture]
public class MultipleAlignTests
{

    private static readonly string[] UserMatrixRows = ["------------------------------ALAAALLLAALAAA", "-RPKRQESPKNPGRQEPKEPENSNRRQKK---------------", "ALAAAEC--LAALAAA----------------------------", "----AGAAALLSEPLAALAAA-----------------------", "----------------ALAAAPLLAALAAA--------------", "---ALGG--ALLNAALAAA-------------------------"];

    private static readonly string[] ClustalWProteinRows =
        ["LFFF---FFFFMFFFFFDKFFFFFFFFFFFFFFFF", "------------FLAFLVKTYAHTGLKDKEGF---", "-FFF---FFFFFFFFFFFFFFFFFFFFFFFFFFFF", "-FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF-"];

    // B09 F66 (A1-13b): multiple_align_tcoffee delegates to TCoffeeAligner.MultipleAlignTCoffeeProgram (the T-Coffee 13.41 program: its own libraries, guide tree, progressive stage).
    private static readonly string[] TCoffeeProgramSet = ["ACGTTGCATGCCGATAGC", "ACGTGCATGCGATAGGC", "ACTTGCATCCGATTGC", "CGTTGCATGGCCGATAGC"];

    // B09 F64 (A1-10d file part): multiple_align_clustalw matrixFile / matrixFiles / pwMatrixFile = clustalw -MATRIX=file (here a CLUSTAL_SERIES file with an unlabelled
    // tab-separated integer matrix and a labelled floating-point matrix) / -PWMATRIX=file. Rows = the installed clustalw 2.1, set 11 of clustalw_user_matrix_oracle.tsv.
    private static readonly string[] UserMatrixSet = ["ALAAALLLAALAAA", "RPKRQESPKNPGRQEPKEPENSNRRQKK", "ALAAAECLAALAAA", "AGAAALLSEPLAALAAA", "ALAAAPLLAALAAA", "ALGGALLNAALAAA"];

    // B09 F51 (A1-15): guide-tree distance source. The set below is discriminating: ClustalW 2.1's pairwise stage and the
    // repository's linear NW alignment give different distances, hence different guide trees and different rows
    // (rows of the independent Python implementation in scripts/alignment/generate_clustalw_pairwise_oracle.py).
    private static readonly string[] DistanceSet = ["CCAACCAACCAA", "ACACTCCGCCCGCACTCAGAACTTT", "TAA"];

    // B09 F52 (A1-10): gapModel = affine honours gapOpen. Identical pair 33 matches = 66; the two pairs with the 23-nt sequence have
    // 23 matches (46) and ONE 10-nt gap run: affine -8 + 10 * -1 = -18, linear -10.
    private static readonly string[] IndelSet =
        ["ACGTACGTACGTAAAAAAAAAAGGCATCGATCG", "ACGTACGTACGTGGCATCGATCG", "ACGTACGTACGTAAAAAAAAAAGGCATCGATCG"];

    // B09 F53 (A1-12): guideTree = upgma | nj. The set is discriminating: UPGMA and ClustalW-style NJ guide trees give different
    // progressive rows (library call Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(..., new ConsistencyMsaOptions { GuideTree = ... }), defaults 1/-1/-2/-1).
    // F56 (A1-10b): the NJ tree now uses ClustalW's own NJ arithmetic (exact ties like ClustalW); the former set became a tie where
    // both trees agree, so the set was replaced by one found by a random search (NJ and UPGMA topologies differ).
    private static readonly string[] GuideTreeSet =
        ["AAAAGTCAGGCAT", "AACAGTTAGGCAT", "AAAGTTAGGCCA", "AAATGTTAGGCAT", "AAAGTTAGGCAA"];
    private static readonly string[] DnaMatrixRows = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "-GAACTCGGAATGGGAGTGG", "-GGACTCGGGATGGGAGTGG"];

    // B09 F67 (A1-10d part 2): multiple_align_clustalw useEndGaps / endGapPenalties. Rows = the patched clustalw 2.1 of the Docker package (clustalw_end_gaps_oracle.tsv, set D84 DNA):
    // control, useendgaps, endgappen, both; set D17 aborts with endGapPenalties (no output file) -> an explicit error.
    private static readonly string[] EndGapSet = ["AAGAGCCCT", "CGATAAAGCT", "CGATAAGTT", "CGATAAGTT", "CGCATAAGTT", "CGTATATTC", "CGAAATCT"];

    // B09 F68 (A1-10d part 3): multiple_align_clustalw dnaMatrixFile / pwDnaMatrixFile / dnaMatrix / pwDnaMatrix / userDnaMatrixUse. Rows = the installed clustalw 2.1 driven through its
    // interactive menu (the file is used; set 244 of clustalw_dna_matrix_oracle.tsv, pwgapopen 10); the command line reads but ignores the file (userDnaMatrixUse = commandline).
    private static readonly string[] DnaMatrixSet = ["GGACTCGGGAATGGGAGTGG", "CGTCTCGGGAATGGGAGTGG", "GAACTCGGAATGGGAGTGG", "GGACTCGGGATGGGAGTGG"];
    private const string DnaMatrixFile = "# DNA matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n A T G C S W R Y K M B V H D N U *\nA 1 -9 1 -9 -9 -1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -4\nT -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\nG 1 -9 1 -9 -1 -9 -1 -9 -1 -9 -1 -1 -9 -1 -1 -9 -4\nC -9 1 -9 1 -1 -9 -9 -1 -9 -1 -1 -1 -1 -9 -1 1 -4\nS -9 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nW -1 -1 -9 -9 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nR -1 -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -9 -4\nY -9 -1 -9 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nK -9 -1 -1 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -4\nM -1 -9 -9 -1 -1 -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -9 -4\nB -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nV -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -9 -4\nH -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nD -1 -1 -1 -9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nN -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\nU -9 1 -9 1 -9 -1 -9 -1 -1 -9 -1 -9 -1 -1 -1 1 -4\n* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4\n";

    // B09 F63 (A1-10c b): multiple_align_clustalw type = protein delegates to the PROTEIN overload of ClustalWAligner.MultipleAlignClustalW.
    // Rows of the default run = the installed clustalw 2.1 (-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT), set 20 of clustalw_protein_msa_oracle.tsv.
    private static readonly string[] ClustalWProteinSet =
        ["LFFFFFFFMFFFFFDKFFFFFFFFFFFFFFFF", "FLAFLVKTYAHTGLKDKEGF", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", "FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF"];

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

    private const string UserSeries = "# series of matrices\nCLUSTAL_SERIES\nMATRIX 0 23 m0.mat\nMATRIX 23 100 m1.mat\n# comment between entries\n";
    private const string UserMatrixM0 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\n\tF\tA\tQ\tK\tC\tI\tY\tW\tH\tD\tT\tE\tL\tN\t-\n1\t-1\t3\t-6\t5\t7\t-5\t-5\t-3\t-8\t3\t-8\t5\t3\t-4\n0\t0\t-1\t-9\t-3\t1\t0\t1\t-5\t8\t0\t-1\t-3\t-1\t-4\n0\t0\t0\t3\t5\t-9\t-2\t-8\t9\t3\t-9\t-7\t-3\t-9\t-4\n-1\t0\t0\t0\t5\t-7\t-9\t2\t1\t2\t5\t0\t-9\t1\t-4\n0\t0\t-1\t-1\t1\t-5\t7\t-3\t0\t-4\t8\t6\t-4\t4\t-4\n0\t0\t0\t0\t0\t0\t-1\t0\t-3\t3\t-2\t-5\t-1\t0\t-4\n1\t0\t0\t0\t0\t0\t1\t-5\t-1\t8\t-3\t-3\t-1\t9\t-4\n0\t-1\t0\t0\t-1\t-1\t0\t2\t6\t9\t-7\t5\t7\t4\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t1\t-8\t-4\t-3\t-6\t-9\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t6\t-9\t6\t-4\t-4\n0\t0\t0\t0\t0\t0\t0\t-1\t0\t0\t0\t-7\t-2\t-5\t-4\n-1\t0\t0\t0\t-1\t0\t0\t-1\t0\t0\t0\t0\t1\t-6\t-4\n0\t0\t0\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t1\t5\t-4\n0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t-4\n-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\t-4\n";

    // B09 F61 (A1-14): type = 'protein' / matrix on every multiple_align* tool; the tools delegate to the SequenceAligner matrix overloads.
    private static readonly string[] Proteins =
    [
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF",
        "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV",
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
        "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
    ];
    private const string UserMatrixM1 = "# Matrix made for the user-matrix oracle\n# Entries for the matrix at a scale of 1/2 bit\n\nR M E G T Q W N *\nR 2 -0.5 -1 2.5 1 4 2 -2.5 -2\nM 0.5 3.5 2.5 -2.5 0.5 -3 -1 -4.5 -2\nE -2 -2 2.5 -3 -4.5 -2 -2.5 -2.5 -2\nG -1 -1.5 0.5 4 1 -3.5 -2.5 4.5 -2\nT -1 -0.5 0 0.5 4 -4 1 1.5 -2\nQ 0.5 -1.5 1 -2.5 0.5 2 3.5 4.5 -2\nW 1.5 0 -1.5 0 0 1.5 4 -0.5 -2\nN 1 1.5 0 1.5 -0.5 1 -2.5 3 -2\n* -2 -2 -2 -2 -2 -2 -2 -2 -2\n";

    // B09 F56 (A1-10b): multiple_align_clustalw delegates to ClustalWAligner.MultipleAlignClustalW (ClustalW 2.1 multiple alignment).
    // Rows of the set below = the installed clustalw 2.1 (-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT).
    private static readonly string[] ClustalWSet =
        ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT", "ACGTAGGATTGACCGTAGCCATCGATGGA", "GGGTTTACCAGTAGGCATCGATGGAAACCT"];
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
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProgressive(bad));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignIterative(bad));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignConsistency(bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentTools.MultipleAlignIterative(new[] { "ACGT", "ACGT" }, maxIterations: 0));

        Assert.That(AlignmentTools.MultipleAlignAnchored(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(MultipleAlignmentTools.MultipleAlignIterative(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(MultipleAlignmentTools.MultipleAlignConsistency(Array.Empty<string>()).AlignedSequences, Is.Empty);
    }

    // Hand-derived (progressive unit tests): ACGT/ACGT/AGT -> ACGT, ACGT, A-GT; SP = 3 - 1 + 3 + 3 = 8.
    [Test]
    public void MultipleAlignProgressive_Binding_MatchesHandDerivedAlignment()
    {
        var r = MultipleAlignmentTools.MultipleAlignProgressive(new[] { "ACGT", "ACGT", "AGT" });

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
        var progressive = MultipleAlignmentTools.MultipleAlignProgressive(input);
        var iterative = MultipleAlignmentTools.MultipleAlignIterative(input);

        Assert.That(iterative.TotalScore, Is.GreaterThanOrEqualTo(progressive.TotalScore));
    }

    // Independent Python T-Coffee reference (B09 F7): CTT/ACT/CAA -> -CTT, AC-T, -CAA.
    [Test]
    public void MultipleAlignConsistency_Binding_MatchesIndependentReference()
    {
        var r = MultipleAlignmentTools.MultipleAlignConsistency(new[] { "CTT", "ACT", "CAA" });

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

    [Test]
    public void MultipleAlignProgressive_DistanceClustalW_UsesClustalWPairwiseStage()
    {
        var clustalw = MultipleAlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "clustalw");
        var linear = MultipleAlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3);
        Assert.Multiple(() =>
        {
            Assert.That(clustalw.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "----------------T--AA----" }));
            Assert.That(linear.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "---------------T---AA----" }));
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "LINEAR").AlignedSequences, Is.EqualTo(linear.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignIterativeAndConsistency_Distance_ParsedAndValidated()
    {
        Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "clustalw"));
        Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "clustalw"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "upgma"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "bogus"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "bogus"));
    }

    // B09 F55 (A1-15b): distance = clustalw-quick (clustalw -QUICKTREE) delegates to MsaDistanceMode.ClustalWQuick with the
    // ktuple / window / topDiagonals / pairGap parameters; the guide-tree stage validates them.
    [Test]
    public void MultipleAlignProgressive_DistanceClustalWQuick_DelegatesAndValidatesParameters()
    {
        string[] set = ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "TTGACCAGTAGGCATCGATCCCAAGT", "ACGTAGGATTGACCGTAGCCATCGATGGA"];
        var dna = set.Select(x => new global::Seqeron.Genomics.Core.DnaSequence(x)).ToList();
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(5, -4, 0, -3);
        var expected = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(dna, scoring,
            new ConsistencyMsaOptions
            {
                Distance = global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWQuick,
                QuickPairwise = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions { KTuple = 3, Window = 5, TopDiagonals = 6, PairGap = 9 },
                GuideTree = global::Seqeron.Genomics.ClustalW.MsaGuideTree.NeighborJoining,
            });
        var r = MultipleAlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "clustalw-quick", "linear", "nj", 3, 5, 6, 9);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(expected.AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "ClustalW_Quick").AlignedSequences, Has.Length.EqualTo(4));
            Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentTools.MultipleAlignProgressive(set, 5, -4, 0, -3, "clustalw-quick", "linear", "upgma", 9));
            Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignIterative(set, 5, -4, 0, -3, 16, "clustalw-quick"));
            Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignConsistency(set, 5, -4, 0, -3, "clustalw-quick"));
        });
    }

    [Test]
    public void MultipleAlignClustalW_DelegatesToTheClustalWPort_AndParsesEveryOption()
    {
        var dna = ClustalWSet.Select(x => new global::Seqeron.Genomics.Core.DnaSequence(x)).ToList();
        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(1, -1, -2, -1);
        var expected = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(dna, null, scoring);
        var r = MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet);
        var opt = new global::Seqeron.Genomics.ClustalW.ClustalWOptions
        {
            GapOpen = 10, GapExtend = 3, PairwiseMode = global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMode.Quick, PairwiseGapOpen = 12, PairwiseGapExtend = 5,
            TransitionWeight = 1f, DivergenceCutoff = 50, NoWeights = true, NegativeMatrix = true,
            Quick = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions { KTuple = 3, Window = 6, TopDiagonals = 7, PairGap = 10 },
        };
        var viaApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(dna, opt, new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(2, -3, -5, -1));
        var viaTool = MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet, 10, 3, "ClustalW-Quick", 12, 5, 1f, 50, true, true, 3, 6, 7, 10, 2, -3, -5, -1);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(expected.AlignedSequences));
            Assert.That(r.TotalScore, Is.EqualTo(expected.TotalScore));
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(ClustalWSet));
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet, distance: "full").AlignedSequences, Is.EqualTo(r.AlignedSequences));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(Array.Empty<string>()), "A3-4: the program refuses an empty input");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet, distance: "linear"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(new[] { "ACGT", "ACGX" }));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet, gapOpen: -3));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Protein_DelegatesToTheProteinClustalWPort_AndParsesEveryOption()
    {
        var def = MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein");
        var opt = new global::Seqeron.Genomics.ClustalW.ClustalWProteinOptions
        {
            Matrix = global::Seqeron.Genomics.ClustalW.ClustalWProfileMatrix.Pam, GapOpen = 12, GapExtend = 0.4,
            PairwiseMode = global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMode.Quick,
            Pairwise = new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions
            {
                Matrix = global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Blosum, GapOpen = 8, GapExtend = 0.3,
                Quick = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions { KTuple = 2, Window = 6, TopDiagonals = 7, PairGap = 5 },
            },
            DivergenceCutoff = 50, NoWeights = true, NegativeMatrix = true, NoPreferredGapPenalties = true, NoHydrophilicGapPenalties = true,
            VariableRegionGapPenalties = true, GapSeparationDistance = 6, HydrophilicResidues = "ACDEFGHIK",
        };
        var scoring = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Pam250, -5.0, -0.5);
        var viaApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(ClustalWProteinSet, opt, scoring);
        var viaTool = MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, 12, 0.4, "clustalw-quick", 8, 0.3, 0.5f, 50, true, true, 2, 6, 7, 5, 1, -1, -5.0, -0.5,
            "protein", "pam", "blosum", true, true, true, 6, "ACDEFGHIK", "PAM250");
        Assert.Multiple(() =>
        {
            Assert.That(def.AlignedSequences, Is.EqualTo(ClustalWProteinRows));
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(viaTool.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(ClustalWProteinSet));
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", matrix: "gonnet", pwMatrix: "gonnet", distance: "full").AlignedSequences, Is.EqualTo(def.AlignedSequences));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(Array.Empty<string>(), type: "protein"), "A3-4: the program refuses an empty input");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "rna"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", matrix: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", pwMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", scoreMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(new[] { "MKVL", "MKJL" }, type: "protein"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(ClustalWProteinSet, type: "protein", gapOpen: -1));
            // the DNA defaults are unchanged by the protein parameters
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet, matrix: "pam", noPreferredGaps: true).AlignedSequences, Is.EqualTo(MultipleAlignmentTools.MultipleAlignClustalW(ClustalWSet).AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignClustalW_EndGapOptions_DelegateAndAbortIsExplicit()
    {
        var control = MultipleAlignmentTools.MultipleAlignClustalW(EndGapSet);
        var pen = MultipleAlignmentTools.MultipleAlignClustalW(EndGapSet, useEndGaps: false, endGapPenalties: true);
        var both = MultipleAlignmentTools.MultipleAlignClustalW(EndGapSet, useEndGaps: true, endGapPenalties: true);
        var viaApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(EndGapSet,
            new global::Seqeron.Genomics.ClustalW.ClustalWOptions { EndGapPenalties = true });
        Assert.Multiple(() =>
        {
            Assert.That(control.AlignedSequences, Is.EqualTo(new[] { "---AAGAGCCCT", "CGATAAAGCT--", "CG-ATAAGTT--", "CG-ATAAGTT--", "CGCATAAGTT--", "CGTATATTC---", "--CGAAATCT--" }));
            Assert.That(pen.AlignedSequences, Is.EqualTo(new[] { "AAGA-GCCCT", "CGATAAAGCT", "CG-ATAAGTT", "CG-ATAAGTT", "CGCATAAGTT", "CGTATAT-TC", "CGA--AATCT" }));
            Assert.That(pen.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(both.AlignedSequences, Is.EqualTo(pen.AlignedSequences));
            var ex = Assert.Throws<global::Seqeron.Genomics.ClustalW.ClustalWAbortedException>(() => MultipleAlignmentTools.MultipleAlignClustalW(new[] { "AT", "AAACA", "AA" }, endGapPenalties: true))!;
            Assert.That(ex.Message, Does.Contain("abandons the alignment"));
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(new[] { "AT", "AAACA", "AA" }, useEndGaps: true).AlignedSequences, Has.Length.EqualTo(3));
        });
    }

    [Test]
    public void MultipleAlignClustalW_DnaMatrices_DelegateAndFollowClustalWsCommandLineAndMenu()
    {
        var applied = MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, dnaMatrixFile: DnaMatrixFile);
        var cli = MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, dnaMatrixFile: DnaMatrixFile, userDnaMatrixUse: "commandline");
        var plain = MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10);
        var clustalw = MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrix: "CLUSTALW", pwDnaMatrix: "clustalw");
        var viaApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(DnaMatrixSet, new global::Seqeron.Genomics.ClustalW.ClustalWOptions
        {
            DnaMatrix = global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix.Clustalw,
            PairwiseDnaMatrix = global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix.Clustalw,
        });
        Assert.Multiple(() =>
        {
            Assert.That(applied.AlignedSequences, Is.EqualTo(DnaMatrixRows));
            Assert.That(cli.AlignedSequences, Is.EqualTo(plain.AlignedSequences));
            Assert.That(applied.AlignedSequences, Is.Not.EqualTo(plain.AlignedSequences));
            Assert.That(clustalw.AlignedSequences, Is.EqualTo(viaApi.AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwGapOpen: 10, pwDnaMatrixFile: DnaMatrixFile).AlignedSequences, Has.Length.EqualTo(4));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwDnaMatrix: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, userDnaMatrixUse: "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, dnaMatrixFile: "CLUSTAL_SERIES\nMATRIX 0 100 x\n"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(DnaMatrixSet, pwDnaMatrixFile: "# nothing\n"));
        });
    }

    [Test]
    public void MultipleAlignClustalW_Protein_UserMatrixFiles_EqualClustalW_AndAreValidated()
    {
        var files = new Dictionary<string, string> { ["m0.mat"] = UserMatrixM0, ["m1.mat"] = UserMatrixM1 };
        var viaTool = MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserSeries, matrixFiles: files);
        var series = global::Seqeron.Genomics.ClustalW.ClustalWMatrixFile.Parse(UserSeries, n => files[n]);
        var viaApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(UserMatrixSet,
            new global::Seqeron.Genomics.ClustalW.ClustalWProteinOptions { UserMatrixSeries = series.Series });
        // a single matrix file (-MATRIX=file) and the pairwise matrix file (-PWMATRIX=file) reach the same code
        var viaSingle = MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserMatrixM1, pwMatrixFile: UserMatrixM0);
        var viaSingleApi = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(UserMatrixSet,
            new global::Seqeron.Genomics.ClustalW.ClustalWProteinOptions
            {
                UserMatrix = global::Seqeron.Genomics.ClustalW.ClustalWUserMatrix.Parse(UserMatrixM1),
                Pairwise = new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions { UserMatrix = global::Seqeron.Genomics.ClustalW.ClustalWUserMatrix.Parse(UserMatrixM0) },
            });
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(UserMatrixRows));
            Assert.That(viaApi.AlignedSequences, Is.EqualTo(UserMatrixRows));
            Assert.That(viaSingle.AlignedSequences, Is.EqualTo(viaSingleApi.AlignedSequences));
            // the user matrices change the result against the built-in series
            Assert.That(MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein").AlignedSequences, Is.Not.EqualTo(UserMatrixRows));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: UserSeries), "the series names files that are not supplied");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", matrixFile: "# nothing\n"), "no header");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignClustalW(UserMatrixSet, type: "protein", pwMatrixFile: UserSeries), "a series where one matrix is expected");
        });
    }

    // B09 F107 (A2-3 / A2-4): the further pair-wise methods of the program through the tool; rows written by t_coffee -seq s.fa -method M -output fasta_aln (installed 13.41); the refused methods raise TCoffeeAbortedException.
    [Test]
    public void MultipleAlignTCoffee_FurtherPairMethods_EqualTheProgram_AndRefusedMethodsAbort()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["biphasic_pair"]).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "ACG-TGCATG--CGATAGGC", "ACT-TGCATC--CGATT-GC", "-CGTTGCATGGCCGATA-GC" }));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["lalign_rs_s_pair"]).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "ACGT-GCATG-C-GATAGGC", "ACTT-GCATC-C-GATT-GC", "-CGTTGCATGGCCGATA-GC" }));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["promo_pair"]).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATGCCGATAGC-", "-ACGTGCATGCGATAGGC-", "-ACTTGCATCCGATTGC--", "-CGTTGCATGGCCGATAGC" }));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["blastr_pair"]).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCAT-GCCGATAGC-", "-ACGTGCAT-G-CGATAGGC", "-ACTTGCAT-C-CGATTGC-", "-CGTTGCATGGCCGATAGC-" }));
            var aborted = Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["hash_pair"]));
            Assert.That(aborted!.ExitStatus, Is.EqualTo(1));
            Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["clean_slow_pair"]));
        });
    }

    // B09 F114-F117 (WP30): the methods that run an external program and the further dynamic programs through the tool. Rows written by the installed t_coffee 13.41 (muscle3 3.8.1551 as `muscle`, TMalign 20190822,
    // blastall 2.12.0 with a local database of the two sequences); the structure aligner's / BLAST's own output is the input.
    [Test]
    public void MultipleAlignTCoffee_ExternalProgramMethods_AndFurtherDynamicPrograms()
    {
        const string tmalign = "TM-score= 0.59749 (if normalized by length of Chain_2)\n(You should use TM-score normalized by length of the reference protein)\n\n(\":\" denotes aligned residue pairs of d < 5.0 A, \".\" denotes other aligned residues)\nDPY-NIVIHMQLAQ\n::: :: :::::::\nDPYYNI-IHMQLAQ\n\n";
        string[] pair = ["DPYNIVIHMQLAQ", "DPYYNIIHMQLAQ"];
        var viaTool = MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["TMalign_pair"], structureTemplates: pair, structureAlignerOutputs: [new StructureAlignerOutputInput("TMalign_pair", 0, 1, tmalign)]);
        string[] blastSet = ["QMPLHKEGVTYTDFMQAPDRKNYNQRFHTNPEHIGKAAIWDSMEDNTQLK", "WEFDMPGSGRKNLESPVYYDAKRIGVAFVGHFIY"];
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(new[] { "DPY-NIVIHMQLAQ", "DPYYNI-IHMQLAQ" }), "TMalign_pair: the alignment TM-align wrote");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["sap_pair"], structureTemplates: pair, structureAlignerOutputs: [new StructureAlignerOutputInput("TMalign_pair", 0, 1, tmalign)]).AlignedSequences,
                Is.EqualTo(viaTool.AlignedSequences), "sap_pair = TM-align (SAP is not installed)");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["TMalign_pair"]).AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(pair), "no structures: proba_pair");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["muscle_pair"]).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCAT-GCCGAT-AGC", "AC-GTGCATG--CGATAGGC", "AC-TTGCAT--CCGAT-TGC", "CG-TTGCATGGCCGAT-AGC" }), "muscle_pair = MUSCLE 3.8.1551 on every pair");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-dp_mode ktup_pair_wise").AlignedSequences, Is.EqualTo(TCoffeeProgramSet), "ktup_pair_wise leaves the rows as they are");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2").AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "AC-GTGCATG--CGATAGGC", "AC-TTGCATC--CGATT-GC", "CG-TTGCATGGCCGATA-GC" }));
            Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-dp_mode precomputed_pair_wise"), "unequal first rows: the program stops");
            Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["fugue_pair"]));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["mafft_pair"]).AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(TCoffeeProgramSet), "mafft_pair = MAFFT 7.505 (pair stage, A2-19 / F119) on every pair; the 460 oracle sets are compared in the library tests");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-mode expresso"), "downloads structures: network");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["blastp_msa"]), "blastp_msa needs the XML of the BLAST search");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["TMalign_pair"], structureTemplates: ["DPYNIVIHMQLAQ"]), "one structure entry per sequence");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["TMalign_pair"], structureTemplates: pair, structureAlignerOutputs: [new StructureAlignerOutputInput("dali_pair", 0, 1, tmalign)]));
        });
        string[] expectedBlast = ["QMPLHKEGVTYTDFMQAPDRKNYNQRFHTNPEHIGKAAIWDSMEDNTQLK", "----------------WEFDMPGSGRKNLESPVYYDAKRIGVAFVGHFIY"];
        string xml = BlastTwoSequenceXml;
        Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(blastSet, ["blastp_msa"], blastXml: xml).AlignedSequences, Is.EqualTo(expectedBlast), "blastp_msa from the HSPs of the local all-against-all search");
    }

    // B09 F298 (A2-12): sap_pair with the REAL SAP 1.1.3 installed: the output of the tool item 'sap_pair' is SAP's own results file (an oracle set of the user-run Docker package A2-12: rows written by
    // t_coffee 13.41 with SAP 1.1.3 on the PATH of its plug-in directory), read by the program's glue; without such an item the method is the TM-align replacement.
    [Test]
    public void MultipleAlignTCoffee_SapPair_WithTheRealSapResults()
    {
        string[] pair = ["TFSVMQWIMNHWCITMGQPMCD", "IHQGVRCMPWSHNKIPNQHWNSSLDGVDTDSQWGSY"];
        const string sap = " file1 = 571_1.4135.sap_tmp\n prot1->Compound =  CHAIN: A;\n file2 = 571_2.2450.sap_tmp\n prot2->Compound =  CHAIN: A;\nARNDCQEGHILKMFPSTWYVBZX\nmatrix constant = 8\n\n score = 0.132805\n\nPercent sel on aln =   41.18\nPercent aln in sel =   31.82\n"
            + " *K   14   0.0   1 T**\n *I   15   0.0   2 F  \n  P   16   0.0   3 S* \n  N   17   0.0   4 V* \n  Q   18   0.1   5 M* \n *H   19   0.1   6 Q* \n**W   20   0.0   7 W* \n *N   21   0.0   8 I* \n *S   22   0.1   9 M* \n *S   23   0.1  10 N  \n  L   24   0.1  11 H* \n  D   25   0.0  12 W* \n"
            + " *G   26   0.0  13 C* \n *V   27   0.0  14 I  \n**D   28   0.0  15 T  \n  T   29   0.0  16 M  \n *D   30   0.0  17 G* \n  S   31   0.1  18 Q  \n *Q   32   0.2  19 P* \n *W   33   0.0  20 M* \n  G   34   0.0  21 C  \n *S   35   0.0  22 D  \n"
            + "571_2.2450.sap_tmp Weighted RMSd =  10.788 over 22 atoms 15.028649 1.041443   4.5 4.087417\n";
        var result = MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["sap_pair"], structureTemplates: pair, structureAlignerOutputs: [new StructureAlignerOutputInput("sap_pair", 0, 1, sap)]);
        Assert.That(result.AlignedSequences, Is.EqualTo(new[] { "-------------TFSVMQWIMNHWCITMGQPMCD-", "IHQGVRCMPWSHNKIPNQHWNSSLDGVDTDSQWGSY" }), "rows written by t_coffee 13.41 with SAP 1.1.3");
        Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(
            () => MultipleAlignmentTools.MultipleAlignTCoffee(pair, ["sap_pair"], structureTemplates: pair, structureAlignerOutputs: [new StructureAlignerOutputInput("sap_pair", 1, 0, sap)]),
            "SAP installed, no results file for the pair of two usable structures: the program stopped in every such oracle set");
    }

    [Test]
    public void MultipleAlignTCoffee_DelegatesToTheProgramPort_ParsesEveryOption_AndValidates()
    {
        // rows written by t_coffee -seq s.fa -method slow_pair,lalign_id_pair -output fasta_aln (installed 13.41)
        string[] program = ["ACGTTGCAT-GCCGATAG-C", "ACGT-GCAT-GC-GATAGGC", "AC-TTGCAT--CCGATTG-C", "-CGTTGCATGGCCGATAG-C"];
        var viaTool = MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["slow_pair", "LALIGN_ID_PAIR"]);
        var viaApi = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram(TCoffeeProgramSet, new global::Seqeron.Genomics.TCoffee.TCoffeeProgramOptions
        {
            Methods = [global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.SlowPair, global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.LalignIdPair],
        });
        var plain = MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet);
        var relax0 = MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["proba_pair"], 3, false, "dna");
        Assert.Multiple(() =>
        {
            Assert.That(viaTool.AlignedSequences, Is.EqualTo(program));
            Assert.That(viaApi.AlignedSequences, Is.EqualTo(program));
            Assert.That(viaTool.TotalScore, Is.EqualTo(viaApi.TotalScore));
            Assert.That(plain.AlignedSequences, Is.EqualTo(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["proba_pair"]).AlignedSequences), "the default method is proba_pair");
            Assert.That(plain.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(TCoffeeProgramSet));
            Assert.That(relax0.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(TCoffeeProgramSet));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(["MKVLAAGIVG", "MKVLAGIVGL", "MRVLAAGIIG"], null, 10, true, "protein").AlignedSequences.Select(r => r.Replace("-", "")),
                Is.EqualTo(new[] { "MKVLAAGIVG", "MKVLAGIVGL", "MRVLAAGIIG" }));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(Array.Empty<string>()), "A3-4: the program refuses an empty input");
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, ["blast_pair"]));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, type: "rna"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, lalignTopAlignments: 0));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(new[] { "ACGT", "" }));
            // A1-13c: RNA (a U) uses the RNA pair-HMM; rows written by t_coffee 13.41 for this set
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(new[] { "ACGUACGUUG", "ACGUUCGUUG", "ACGAACGUUUG" }).AlignedSequences, Is.EqualTo(new[] { "ACGUACGU-UG", "ACGUUCGU-UG", "ACGAACGUUUG" }));
        });
    }

    // B09 F71-F75 (A1-13c / A1-13d): the 'options' flag string and 'names' of multiple_align_tcoffee (rows written by the installed t_coffee 13.41 for the same flags).
    [Test]
    public void MultipleAlignTCoffee_Options_AndNames_EqualTheProgram_AndAreValidated()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee").AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "AC-GTGCATG--CGATAGGC", "AC-TTGCATC--CGATT-GC", "CG-TTGCATGGCCGATA-GC" }));
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-relax_lib 2 -gapext -2 -nomatch 5").AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "AC-GTGCATG--CGATAGGC", "AC-TTGCATC--CGATT-GC", "CG-TTGCATGGCCGATA-GC" }));
            string[] names = ["zeta", "alpha", "mid", "beta"];
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-outorder aligned", names: names).AlignedSequences,
                Is.EqualTo(new[] { "AC-GTGCATG--CGATAGGC", "CG-TTGCATGGCCGATA-GC", "AC-TTGCATC--CGATT-GC", "ACGTTGCATG-CCGATA-GC" }), "the program's internal (sorted-name) order");
            Assert.That(MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, names: names).AlignedSequences,
                Is.EqualTo(new[] { "ACGTTGCATG-CCGATA-GC", "AC-GTGCATG--CGATAGGC", "AC-TTGCATC--CGATT-GC", "CG-TTGCATGGCCGATA-GC" }));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-no_such_flag 1"));
            // -tree_mode kmeans is accepted by the program, which then ends WITHOUT writing an alignment (exit status 0): TCoffeeAbortedException (A2-6)
            Assert.That(Assert.Throws<global::Seqeron.Genomics.TCoffee.TCoffeeAbortedException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, options: "-tree_mode kmeans"))!.ExitStatus, Is.Zero);
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(TCoffeeProgramSet, names: ["a", "b"]));
        });
    }

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
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(IndelSet, 2, -3, -8, -1, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(MultipleAlignmentTools.MultipleAlignIterative(IndelSet, 2, -3, -8, -1, 16, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.That(MultipleAlignmentTools.MultipleAlignConsistency(IndelSet, 2, -3, -8, -1, "linear", "affine").TotalScore, Is.EqualTo(66 + 2 * (46 - 18)));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(IndelSet, 2, -3, -8, -1, "quadratic"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignAnchored(IndelSet, 2, -3, -8, -1, "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProgressive(IndelSet, 2, -3, -8, -1, "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignIterative(IndelSet, 2, -3, -8, -1, 16, "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignConsistency(IndelSet, 2, -3, -8, -1, "linear", "bogus"));
        });
    }

    [Test]
    public void MultipleAlignProgressive_GuideTreeNj_SelectsNeighborJoiningTree()
    {
        var upgma = MultipleAlignmentTools.MultipleAlignProgressive(GuideTreeSet);
        var nj = MultipleAlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "nj");
        Assert.Multiple(() =>
        {
            Assert.That(upgma.AlignedSequences, Is.EqualTo(new[] { "AA-AAGTCAGG-CAT", "AACA-GTTAGG-CAT", "AA-A-GTTAGGCCA-", "AA-ATGTTAGG-CAT", "AA-A-GTTAGG-CAA" }));
            Assert.That(nj.AlignedSequences, Is.EqualTo(new[] { "AA-AAGTCAGGCAT", "AACA-GTTAGGCAT", "AA-A-GTTAGGCCA", "AA-ATGTTAGGCAT", "AA-A-GTTAGGCAA" }));
            Assert.That(upgma.TotalScore, Is.EqualTo(83));
            Assert.That(nj.TotalScore, Is.EqualTo(81));
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "UPGMA").AlignedSequences, Is.EqualTo(upgma.AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "neighbor_joining").AlignedSequences, Is.EqualTo(nj.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignIterativeAndConsistency_GuideTree_ParsedAndValidated()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignIterative(GuideTreeSet, 1, -1, -2, -1, 16, "clustalw", "linear", "nj"));
            Assert.DoesNotThrow(() => MultipleAlignmentTools.MultipleAlignConsistency(GuideTreeSet, 1, -1, -2, -1, "clustalw", "linear", "nj"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProgressive(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignIterative(GuideTreeSet, 1, -1, -2, -1, 16, "linear", "linear", "bogus"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignConsistency(GuideTreeSet, 1, -1, -2, -1, "linear", "linear", "bogus"));
        });
    }

    // B09 F58 (A1-13): consistencyEngine = 'tcoffee' routes the progressive stage through the T-Coffee 13.41 port; default unchanged.
    [Test]
    public void MultipleAlignConsistency_ConsistencyEngine_TCoffeeIsValid_DefaultUnchanged_BogusRejected()
    {
        var raw = new[] { "ACGTTACGT", "ACGACGT", "ACGTACGA", "TTACGTAG" };
        var tcoffee = MultipleAlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "tcoffee");
        var paper = MultipleAlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "paper");

        Assert.Multiple(() =>
        {
            Assert.That(tcoffee.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(raw));
            Assert.That(paper.AlignedSequences, Is.EqualTo(MultipleAlignmentTools.MultipleAlignConsistency(raw).AlignedSequences));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignConsistency(raw, 1, -1, -2, -1, "linear", "linear", "upgma", null, null, null, null, "bogus"));
        });
    }

    [Test]
    public void MultipleAlign_Protein_DelegatesToTheMatrixOverloads_AllMethods()
    {
        var blosum = new global::Seqeron.Genomics.Infrastructure.MatrixScoring(global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var affine = new ConsistencyMsaOptions { GapModel = global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine };
        var linear = new ConsistencyMsaOptions();
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
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "linear", "affine", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(Proteins, blosum, affine).AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignIterative(Proteins, 1, -1, -9.5, -0.5, 16, "linear", "linear", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(Proteins, blosum, linear).AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignConsistency(Proteins, 1, -1, -9.5, -0.5, "linear", "affine", "upgma", null, null, null, null, "paper", "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Proteins, blosum, affine).AlignedSequences));
        });
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
        var full = new ConsistencyMsaOptions
        {
            GapModel = affine,
            Distance = global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWAffine,
            ClustalWProtein = new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions
            {
                Matrix = global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Pam, GapOpen = 5, GapExtend = 0.5,
            },
        };
        var quick = new ConsistencyMsaOptions
        {
            GapModel = affine,
            Distance = global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWQuick,
            ClustalWProtein = new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions
            {
                Quick = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions { KTuple = 2, Window = 6, TopDiagonals = 7, PairGap = 5 },
            },
        };
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw", "affine", "upgma", null, null, null, null, "protein", null, "pam", 5, 0.5).AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(Proteins, blosum, full).AlignedSequences));
            Assert.That(MultipleAlignmentTools.MultipleAlignIterative(Proteins, 1, -1, -9.5, -0.5, 16, "clustalw", "affine", "nj", null, null, null, null, "protein", null, "pam", 5, 0.5).AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(Proteins, blosum, full with { GuideTree = global::Seqeron.Genomics.ClustalW.MsaGuideTree.NeighborJoining }).AlignedSequences));
            // quick: ktuple 2 given -> window / topDiagonals / pairGap explicit
            Assert.That(MultipleAlignmentTools.MultipleAlignConsistency(Proteins, 1, -1, -9.5, -0.5, "clustalw-quick", "affine", "upgma", 2, 6, 7, 5, "paper", "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Proteins, blosum, quick).AlignedSequences));
            // defaults: ClustalW's protein defaults (GONNET, 10 / 0.1; quick 1 / 5 / 5 / 3)
            Assert.That(MultipleAlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw-quick", "affine", "upgma", null, null, null, null, "protein").AlignedSequences,
                Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(Proteins, blosum,
                    new ConsistencyMsaOptions { GapModel = affine, Distance = global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWQuick }).AlignedSequences));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignProgressive(Proteins, 1, -1, -9.5, -0.5, "clustalw", "affine", "upgma", null, null, null, null, "protein", null, "bogus"));
        });
    }}

