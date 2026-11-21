using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{
    private const string TandemSeq = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACACGTTAACGTTAACGTTAACGTTAACGTTAacgtNN";

    [Test]
    public void Consambig_UnionOfTheColumnCodes_LowerCaseAtGaps()
    {
        // consambig 6.6.0: AC-GTN / AC-GTA / ARGGT- -> AVgGTn; D|N = B, I|L = J for proteins; an all-gap protein column is the NUL character (B09 A3-13d, F147)
        Assert.That(AnalysisTools.Consambig(["AC-GTN", "AC-GTA", "ARGGT-"]).Consensus, Is.EqualTo("AVgGTn"));
        Assert.That(AnalysisTools.Consambig(["DI", "NL"], "protein").Consensus, Is.EqualTo("BJ"));
        Assert.That(AnalysisTools.Consambig(["MK-LA", "MK-LD"]).Consensus, Is.EqualTo("MK\0LX"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Consambig(["ACGT"]));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Consambig(["ACGT", "ACGT"], "rna"));
    }

    [Test]
    public void ShowAlign_RendersTheAlignmentAsShowalign()
    {
        // showalign 6.6.0 on s1 acgtac-gtaGT / s2 acgtaccgtaGA / s3 acctaccgtTGA (default options) and with -show I -refseq 2 (B09 A3-13c, F148)
        var r = AnalysisTools.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"]);
        Assert.That(r.Text, Is.EqualTo("                  10  \n          ----:----|--\ns1        ......-....T\ns2        ............\ns3        ..C......T..\nConsensus acgtaccgtaGA\n\n"));
        var identities = AnalysisTools.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"], refseq: "2", show: "I");
        Assert.That(identities.Text, Does.Contain("s1        acgtac-gtaG.\ns3        ac.taccgt.GA\nConsensus acgtaccgtaGA\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowAlign(["a", "b"], ["AC", "AC"], show: "Z"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowAlign(["a", "b"], ["AC", "AC"], uppercase: "5"));
    }

    [Test]
    public void ShowAlign_Html_WrapsAndColoursAsShowalign()
    {
        // showalign 6.6.0 -html -highlight "3-8 blue 10-12 orange" -width 30 (A3-13c2): the ranges are coloured on the sequence lines and the consensus, not on the names or ticks
        var r = AnalysisTools.ShowAlign(["a", "b"], ["ACGTACGTAAACCCGGG", "ACGTACGAAAACCCGGT"], html: true, highlight: "3-8 blue 10-12 orange", width: 30);
        Assert.That(r.Text, Does.StartWith("<pre>\n").And.EndWith("</pre>\n"));
        Assert.That(r.Text, Does.Contain("a         ..<font color=blue>......</font>.<font color=orange>...</font>").And.Contain("Consensus AC<font color=blue>GTACGt</font>A<font color=orange>AAC</font>CCGGg"));
        Assert.That(AnalysisTools.ShowAlign(["a", "b"], ["AC", "AC"], highlight: "1-2 red").Text, Does.Not.Contain("<font"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowAlign(["a", "b"], ["AC", "AC"], html: true, highlight: "0-9 red"));
    }

    [Test]
    public void Diffseq_SnpsAndTheSevenTagReport_AsTheProgram()
    {
        // diffseq 6.6.0 on ACGTACGTACGTAAGGCCTTAACCGGTT / ACGTACGTACGTAAGGCATTAACCGGTT -wordsize 8 (A5-8c3, F196): one SNP, C -> A, a transversion
        var r = AnalysisTools.Diffseq("ACGTACGTACGTAAGGCCTTAACCGGTT", "ACGTACGTACGTAAGGCATTAACCGGTT", 8, false, "a", "b");
        Assert.That(r.Differences, Has.Length.EqualTo(1));
        Assert.That(r.Differences[0], Is.EqualTo(new DiffseqDifferenceItem(18, 18, 1, 18, 18, 1)));
        Assert.That((r.Snps, r.Transitions, r.Transversions, r.OverlapStart1, r.OverlapEnd1), Is.EqualTo((1, 0, 1, 1, 28)));
        Assert.That(r.Report, Does.Contain("\n\na 18-18 Length: 1\nSequence: C\nSequence: A\nb 18-18 Length: 1\n"));
        Assert.That(r.Report, Does.Contain("# SNP_count: 1\n# Transitions: 0\n# Transversions: 1\n"));
        Assert.That(r.FeatureFile1, Does.Contain("a\tdiffseq\tsequence_conflict\t18\t18\t1\t+\t.\tID=a.1;note=SNP in b;replace=A"));
        Assert.That(AnalysisTools.Diffseq("ACGTACGTACGTAAGGCCTTAACCGGTT", "ACGTACGTACGTAAGGCATTAACCGGTT", 8, false, "a", "b", "-rformat excel").Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\tstart\tend\tlength\tname\tsequence\tfirst_feature\tsecond_feature\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Diffseq("ACGT", "ACGT", 8, false, "a", "b", "pir"));
    }

    [Test]
    public void Diffseq_SequenceQualifiers_AsTheProgram()
    {
        // diffseq 6.6.0 qa.fa qb.fa -wordsize 5 -sbegin1 11 -send1 70 -sbegin2 11 -send2 70 -rformat table (A5-8c3b, F216): the SNP of position 31 is at 21 of the cut sequences, the header keeps the original range of the first and 1 - 60 for the second
        const string a = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        const string b = "TTTCCTCATGCAATTCAAAACCATGTCCGTCATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        var r = AnalysisTools.Diffseq(a, b, 5, false, "seqa", "seqb", "-sbegin1 11 -send1 70 -sbegin2 11 -send2 70 -rformat table");
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(21, 21, 1, 21, 21, 1) }));
        Assert.That(r.Report, Does.Contain("# Sequence: seqa     from: 11   to: 70\n"));
        Assert.That(r.Report, Does.Contain("# Compare: seqb     from: 1   to: 60\n"));
        Assert.That(r.Report, Does.Contain("     21      21       +   0.000     21     21      1 seqb   C        .             .             \n"));
        Assert.That(r.FeatureFile1, Does.Contain("seqa\tdiffseq\tsequence_conflict\t21\t21\t1\t+\t.\tID=seqa.1;note=SNP in seqb;replace=C"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Diffseq(a, b, 5, false, "seqa", "seqb", "-sbegin3 4"));
    }

    [Test]
    public void Diffseq_GenbankAndEmblFeatures_AsTheProgram()
    {
        // diffseq 6.6.0 a.gb b.embl -wordsize 5 (A5-8c3b, F218): the features of both entries that overlap the difference (the join is one feature from its first to its last base, complement is a strand only)
        const string seq = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        string genbank = "LOCUS       SEQA 80 bp    DNA     linear   UNK 01-JAN-2000\nDEFINITION  first.\nACCESSION   SEQA\nVERSION     SEQA.1\nREFERENCE   1  (bases 1 to 80)\nFEATURES             Location/Qualifiers\n"
            + "     source          1..80\n     gene            join(20..40,45..60)\n                     /gene=\"abc\"\n     CDS             complement(25..35)\n                     /note=\"a note\"\n                     /translation=\"MKV\"\n"
            + "ORIGIN\n        1 " + string.Join(' ', Enumerable.Range(0, 8).Select(i => seq.Substring(i * 10, 10).ToLowerInvariant())) + "\n//\n";
        string snp = seq[..30] + "C" + seq[31..];
        string embl = "ID   SEQB; SV 1; linear; genomic DNA; STD; UNC; 80 BP.\nXX\nDE   second\nXX\nFT   source          1..80\nFT   misc_feature    28..33\nFT                   /gene=\"x1\"\nXX\nSQ   Sequence 80 BP;\n     "
            + string.Join(' ', Enumerable.Range(0, 8).Select(i => snp.Substring(i * 10, 10).ToLowerInvariant())) + "   80\n//\n";
        var r = AnalysisTools.Diffseq(genbank, embl, 5, false);
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(31, 31, 1, 31, 31, 1) }));
        Assert.That(r.Report, Does.Contain("Feature: gene 20-60 gene='abc'\nFeature: CDS 25-35 note='a note'\n"));
        Assert.That(r.Report, Does.Contain("Feature: misc_feature 28-33 gene='x1'\n"));
        Assert.That(r.Report, Does.Contain("SEQA 31-31 Length: 1\n"));
        // a special qualifier whose check is not ported is refused
        string special = genbank.Replace("/translation=\"MKV\"", "/collection_date=\"2020\"", StringComparison.Ordinal);
        Assert.Throws<ArgumentException>(() => AnalysisTools.Diffseq(special, embl, 5, false));
    }

    [Test]
    public void Wordmatch_RabinKarpMatches_AsTheProgram()
    {
        // wordmatch 6.6.0 on p = ACGTACGTAC against g = TTACGTACGTACTT, -wordsize 6 (A2-28d, F198): the two words TACGTA / CGTACG are extended to the right, the match of the second starts inside the first
        var r = AnalysisTools.Wordmatch(">p\nACGTACGTAC\n", ">g\nTTACGTACGTACTT\n", 6);
        Assert.That(r.Matches, Is.EqualTo(new[] { new WordmatchHitItem("p", 3, "g", 1, 7), new WordmatchHitItem("p", 1, "g", 3, 9) }));
        Assert.That((r.PatternCount, r.AllMatches, r.SumLength), Is.EqualTo((4, 2L, 16L)));
        Assert.That(r.Alignments, Does.Contain("     7 p               +        4..10       g               +        2..8\n     9 p               +        2..10       g               +        4..12\n"));
        Assert.That(r.Log, Does.Contain("Average match length: 8.00").And.Contain("TACGTA :            1             1              7.00"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Wordmatch("", ">g\nAC\n"));
    }

    [Test]
    public void Fuzzpro_WithMismatches_AsTheProgram()
    {
        // fuzzpro 6.6.0 -pmismatch 1 'H-H-K' on MKHHKLLHHKAAHHKMMHHRPPQQ (A5-8c2, F195): the last window HHR has one mismatch (score 2)
        var r = AnalysisTools.Fuzzpro("MKHHKLLHHKAAHHKMMHHRPPQQ", "H-H-K", mismatch: 1);
        Assert.That(r.Hits.Select(h => (h.Position, h.Length, h.Mismatches)), Is.EqualTo(new[] { (2, 3, 0), (7, 3, 0), (12, 3, 0), (17, 3, 1) }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzzpro("MKHHK", "[HK]-[HK]", mismatch: 2));      // type 6: the mismatches must be below the real pattern length
    }

    [Test]
    public void Megamerger_MergesByWordMatchAnchors_AsTheProgram()
    {
        // megamerger a=ACGTACGTAC b=ACGTACGTACGGGGGG -wordsize 4 (B09 A3-13b, F149)
        var r = AnalysisTools.Megamerger("ACGTACGTAC", "ACGTACGTACGGGGGG", 4, name1: "a", name2: "b");
        Assert.That(r.MergedSequence, Is.EqualTo("acgtacgtacgggggg"));
        Assert.That(r.Report, Does.Contain("Matching region a 1-10 : b 1-10\nLength of match: 10\n").And.Contain("Using b 11-16 as the final sequence"));
        Assert.That(r.Matches.Single(), Is.EqualTo(new WordMatchItem(0, 1, 0, 0, 10)));
        Assert.That(AnalysisTools.Megamerger("ACGTACGTACGGGGGG", "ACGTACGTAC", 12).MergedSequence, Is.Empty);
        Assert.Throws<ArgumentException>(() => AnalysisTools.Megamerger("ACGTACGTAC", "ACGTACGTACGGGGGG", 20));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Megamerger("", "ACGT"));
    }

    [Test]
    public void NucmerClusters_ClusterTheMummerMatchesLikeMgaps()
    {
        // mummer -mum -b -l 15 | mgaps -l 65 -s 90 -d 5 -f 0.12: "6 1 65 none - -" and "75 70 66 none 4 4" (B09 A3-9a, F150)
        const string reference = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA";
        const string query = "AGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA";
        var r = AnalysisTools.NucmerClusters(reference, query, 15);
        var cluster = r.Clusters.Single();
        Assert.That(cluster.Strand, Is.EqualTo("forward"));
        Assert.That(cluster.Matches, Is.EqualTo(new[] { new NucmerClusterMatchItem(6, 1, 65, 0), new NucmerClusterMatchItem(75, 70, 66, 0) }));
        Assert.That(AnalysisTools.NucmerClusters(reference, query, 15, minClusterLength: 200).Clusters, Is.Empty);
        Assert.Throws<ArgumentException>(() => AnalysisTools.NucmerClusters(reference, query, mode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.NucmerClusters("", query));
    }

    [Test]
    public void Nucmer_WritesTheDeltaFileOfTheProgram()
    {
        // nucmer -l 15 ref.fa qry.fa (MUMmer 3.23): ">ref qry 140 135", "6 140 1 135 3 3 0" (three substitutions, no indel), "0" (B09 A3-9b, F199)
        const string reference = ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";
        const string query = ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";
        var r = AnalysisTools.Nucmer(reference, query, minMatch: 15, referenceName: "/tmp/b09nuc/ref.fa", queryName: "/tmp/b09nuc/qry.fa");
        Assert.That(r.Delta, Is.EqualTo("/tmp/b09nuc/ref.fa /tmp/b09nuc/qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n"));
        Assert.That(r.Cluster, Is.Null);
        Assert.That(AnalysisTools.Nucmer(reference, query, minMatch: 15, delta: false).Cluster, Does.StartWith("reference.fa query.fa\nNUCMER\n>ref qry 140 135\n"));
        Assert.That(AnalysisTools.Nucmer(reference, query, minMatch: 15, breakLength: 20000).Messages, Does.Contain("WARNING: Invalid break length 20000, ignoring"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Nucmer(reference, query, matchMode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Nucmer("", query));
        Assert.Throws<FormatException>(() => AnalysisTools.Nucmer(reference, ">q\n"));
    }

    [Test]
    public void MummerDeltaPrograms_WriteTheOutputOfTheInstalledPrograms()
    {
        // MUMmer 3.23 on the three-substitution example (delta-filter -1, show-coords, show-snps -x 2, show-aligns -w 50) and on two contigs (show-diff -r, show-tiling, delta-filter -q) (B09 A3-9b, F200 - F202)
        const string delta1 = "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n";
        const string ref1 = ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";
        const string qry1 = ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";
        const string delta2 = "ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n";
        Assert.That(AnalysisTools.DeltaFilter(delta1, oneToOne: true).Text, Is.EqualTo("ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n"));
        Assert.That(AnalysisTools.DeltaFilter(delta2, queryLis: true).Text, Is.EqualTo("ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n"));
        Assert.That(AnalysisTools.ShowCoords(delta1).Text, Is.EqualTo("ref.fa qry.fa\nNUCMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]\n=====================================================================================\n       6      140  |        1      135  |      135      135  |    97.78  | ref\tqry\n"));
        Assert.That(AnalysisTools.ShowCoords(delta1, header: false, coverage: true, lengths: true, tabular: true, sortByReference: true).Text, Is.EqualTo("6\t140\t1\t135\t135\t135\t97.78\t140\t135\t96.43\t100.00\tref\tqry\n"));
        Assert.That(AnalysisTools.ShowSnps(delta1, ref1, qry1, context: 2).Text, Is.EqualTo("ref.fa qry.fa\nNUCMER\n\n    [P1]  [SUB]  [P2]      |   [BUFF]   [DIST]  |  [R]  [Q]  |  [CTX R]  [CTX Q]  | [FRM]  [TAGS]\n=====================================================================================================\n      71   A G   66        |        1       66  |    0    0  |    TTAAG    TTGGG  |  1  1  ref\tqry\n      72   A G   67        |        1       67  |    0    0  |    TAAGT    TGGGG  |  1  1  ref\tqry\n      74   T G   69        |        2       67  |    0    0  |    AGTAA    GGGAA  |  1  1  ref\tqry\n"));
        Assert.That(AnalysisTools.ShowAligns(delta1, ref1, qry1, "ref", "qry", screenWidth: 50).Text, Is.EqualTo("ref.fa qry.fa\n\n==================================================\n-- Alignments between ref and qry\n\n-- BEGIN alignment [ +1 6 - 140 | +1 1 - 135 ]\n\n\n6          agacaattacataacatacacgtcagcacgaaacttgtt\n1          agacaattacataacatacacgtcagcacgaaacttgtt\n                                                  \n\n45         ggcccagtgtgaatcgcttaagggttaagtaagtgtgat\n40         ggcccagtgtgaatcgcttaagggttggggaagtgtgat\n                                     ^^ ^         \n\n84         gcatacgcctttacttgctgtgtccaccccatcggactg\n79         gcatacgcctttacttgctgtgtccaccccatcggactg\n                                                  \n\n123        gcatttttattacactca\n118        gcatttttattacactca\n                             \n\n\n--   END alignment [ +1 6 - 140 | +1 1 - 135 ]\n\n==================================================\n"));
        Assert.That(AnalysisTools.ShowDiff(delta2, referenceDiff: true).Text, Is.EqualTo("ref2.fa qry2.fa\nNUCMER\n\n[SEQ]\t[TYPE]\t[S1]\t[E1]\t[LEN 1]\nr\tSEQ\t321\t280\t-40\tq1\tq2\n"));
        Assert.That(AnalysisTools.ShowTiling(delta2).Text, Is.EqualTo(">r 600 bases\n1\t320\t-40\t320\t100.00\t100.00\t+\tq1\n281\t600\t0\t320\t100.00\t100.00\t-\tq2\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowDiff(delta2, referenceDiff: true, queryDiff: true));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.DeltaFilter(delta1, minIdentity: 101));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowCoords(""));
    }

    // perl's hash order shuffles the per-substitution lines of the report from run to run: the blocks are compared as sets of lines
    private static string CanonicalReport(string report)
        => string.Join("\n\n", report.Split("\n\n").Select(block => string.Join("\n", block.Split('\n').Order(StringComparer.Ordinal))));

    [Test]
    public void Dnadiff_ReportsTheAlignmentsOfTwoSequenceSets()
    {
        // dnadiff -d on two contigs of one reference: report of the installed program (nothing unaligned: no .unref / .unqry)
        var r = AnalysisTools.Dnadiff(">r\nTTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAACGTTACTAAGGGGTATAATCTTCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTGCGC\n", ">q1\nTTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGA\n>q2\nGCGCACAAGACTGCTGTAATTGGATTGCAAGTCTACTAAGTGTTCCCACCCACAAATAGAAGATTATACCCCTTAGTAACGTTTATGAATGAGCCCGTCTTACCTGGAGACCGAGTCCCGATCTTGGACCAAGGTTCGACAAACTTAGTTCGGGTGGGATTTATGCTCCGTTAGGAATACTGGGTGGGCCGCGGGGGGAATAGCTATAGCCCTGGGTTGCTAGTCAACGAAGCAAATGTATATTATTAGACGCGACCGCAGTAAGAGCGTATGTACACATTCTCCCTAGGTTGCACATGAAGAATGAGCCGGTATGCAGG\n", "ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n");
        Assert.That(CanonicalReport(r.Report), Is.EqualTo(CanonicalReport("ref2.fa qry2.fa\nNUCMER\n\n                               [REF]                [QRY]\n[Sequences]\nTotalSeqs                          1                    2\nAlignedSeqs               1(100.00%)           2(100.00%)\nUnalignedSeqs               0(0.00%)             0(0.00%)\n\n[Bases]\nTotalBases                       600                  640\nAlignedBases            600(100.00%)         640(100.00%)\nUnalignedBases              0(0.00%)             0(0.00%)\n\n[Alignments]\n1-to-1                             2                    2\nTotalLength                      640                  640\nAvgLength                     320.00               320.00\nAvgIdentity                   100.00               100.00\n\nM-to-M                             2                    2\nTotalLength                      640                  640\nAvgLength                     320.00               320.00\nAvgIdentity                   100.00               100.00\n\n[Feature Estimates]\nBreakpoints                        2                    0\nRelocations                        0                    0\nTranslocations                     1                    0\nInversions                         0                    0\n\nInsertions                         0                    0\nInsertionSum                       0                    0\nInsertionAvg                    0.00                 0.00\n\nTandemIns                          0                    0\nTandemInsSum                       0                    0\nTandemInsAvg                    0.00                 0.00\n\n[SNPs]\nTotalSNPs                          0                    0\nAC                          0(0.00%)             0(0.00%)\nAT                          0(0.00%)             0(0.00%)\nAG                          0(0.00%)             0(0.00%)\nTA                          0(0.00%)             0(0.00%)\nTG                          0(0.00%)             0(0.00%)\nTC                          0(0.00%)             0(0.00%)\nGT                          0(0.00%)             0(0.00%)\nGC                          0(0.00%)             0(0.00%)\nGA                          0(0.00%)             0(0.00%)\nCT                          0(0.00%)             0(0.00%)\nCG                          0(0.00%)             0(0.00%)\nCA                          0(0.00%)             0(0.00%)\n\nTotalGSNPs                         0                    0\nTA                          0(0.00%)             0(0.00%)\nTG                          0(0.00%)             0(0.00%)\nTC                          0(0.00%)             0(0.00%)\nGT                          0(0.00%)             0(0.00%)\nGC                          0(0.00%)             0(0.00%)\nGA                          0(0.00%)             0(0.00%)\nCG                          0(0.00%)             0(0.00%)\nCT                          0(0.00%)             0(0.00%)\nCA                          0(0.00%)             0(0.00%)\nAG                          0(0.00%)             0(0.00%)\nAT                          0(0.00%)             0(0.00%)\nAC                          0(0.00%)             0(0.00%)\n\nTotalIndels                        0                    0\nA.                          0(0.00%)             0(0.00%)\nT.                          0(0.00%)             0(0.00%)\nG.                          0(0.00%)             0(0.00%)\nC.                          0(0.00%)             0(0.00%)\n.C                          0(0.00%)             0(0.00%)\n.T                          0(0.00%)             0(0.00%)\n.G                          0(0.00%)             0(0.00%)\n.A                          0(0.00%)             0(0.00%)\n\nTotalGIndels                       0                    0\nT.                          0(0.00%)             0(0.00%)\nG.                          0(0.00%)             0(0.00%)\nC.                          0(0.00%)             0(0.00%)\nA.                          0(0.00%)             0(0.00%)\n.A                          0(0.00%)             0(0.00%)\n.C                          0(0.00%)             0(0.00%)\n.T                          0(0.00%)             0(0.00%)\n.G                          0(0.00%)             0(0.00%)\n")));
        Assert.That(r.UnalignedReference, Is.Null);
        Assert.That(r.UnalignedQuery, Is.Null);
        Assert.Throws<ArgumentException>(() => AnalysisTools.Dnadiff("", ">q\nACGT\n"));
    }

    [Test]
    public void Est2genome_AlignsASplicedEstLikeTheProgram()
    {
        // est2genome -estsequence e.fa -genomesequence g.fa (EMBOSS 6.6.0): three exons, two GT/AG introns (B09 A2-28d, F203)
        string est = ">est1 an est\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCT\n";
        string genome = ">genomic some description\nACGTGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTTTTT\n";
        var r = AnalysisTools.Est2genome(est, genome);
        Assert.That(r.Text, Is.EqualTo("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\nExon        60 100.0     5    64 genomic          1    60 est1          an est\n+Intron    -20   0.0    65   138 genomic     \nExon        80 100.0   139   218 genomic         61   140 est1          an est\n+Intron    -20   0.0   219   262 genomic     \nExon        50 100.0   263   312 genomic        141   190 est1          an est\n\nSpan       150 100.0     5   312 genomic          1   190 est1          an est\n\nSegment     60 100.0     5    64 genomic          1    60 est1          an est\nSegment     80 100.0   139   218 genomic         61   140 est1          an est\nSegment     50 100.0   263   312 genomic        141   190 est1          an est\n"));
        Assert.That(r.Messages, Is.Empty);
        Assert.That(AnalysisTools.Est2genome(est, genome, minScore: 200).Text, Is.EqualTo("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\n"));
        Assert.Throws<InvalidOperationException>(() => AnalysisTools.Est2genome(est, genome, mode: "forward", best: false));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Est2genome(est, genome, mode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Est2genome("", genome));
    }

    [Test]
    public void FindTandemRepeatsEmboss_Etandem_ReportsScoreIdentityAndConsensus()
    {
        // etandem -minrepeat 2 -maxrepeat 8 -threshold 10: "33 61 96 3 12 100.0 cag" and "26 137 172 6 6 94.4 acgtta"
        var r = AnalysisTools.FindTandemRepeatsEmboss(TandemSeq, "etandem", 2, 8, 10);
        Assert.That(r.Program, Is.EqualTo("etandem"));
        Assert.That(r.Repeats.Select(h => (h.Position + 1, h.Position + h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus)),
            Is.EqualTo(new[] { (61, 96, 33, 3, 12, (double?)100.0, "cag"), (137, 172, 26, 6, 6, (double?)94.4, "acgtta") }));
    }

    [Test]
    public void FindTandemRepeatsEmboss_Equicktandem_ReportsDiagonalScores()
    {
        // equicktandem -maxrepeat 6 -threshold 10 on (ACG)x8: "21 1 24 3 8"
        var r = AnalysisTools.FindTandemRepeatsEmboss("ACGACGACGACGACGACGACGACG", "EquickTandem", maxRepeat: 6, threshold: 10);
        var h = r.Repeats.Single();
        Assert.That((h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus), Is.EqualTo((0, 24, 21, 3, 8, (double?)null, (string?)null)));
    }

    [Test]
    public void FindTandemRepeatsEmboss_RejectsBadArguments()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "tandem"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeatsEmboss("", "etandem"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "etandem", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "equicktandem", threshold: 0));
    }

    [Test]
    public void SeqMatchAll_ReportsMaximalExactMatchesOfEveryPair()
    {
        // seqmatchall -wordsize 6 on these three sequences: "10 a 1..10 b 3..12", "6 a 1..6 b 7..12", "6 a 5..10 b 3..8" (1-based); the third is shorter than the word size
        var r = AnalysisTools.SeqMatchAll(["ACGTACGTAC", "GGACGTACGTACTT", "TTTT"], 6);
        Assert.That(r.Matches.Select(m => (m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray(),
            Is.EqualTo(new[] { (0, 1, 0, 2, 10), (0, 1, 0, 6, 6), (0, 1, 4, 2, 6) }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SeqMatchAll([], 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SeqMatchAll(["ACGT", ""], 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.SeqMatchAll(["ACGT", "ACGT"], 0));
    }

    [Test]
    public void Fuzzpro_ReportsEveryWindowOfTheProteinPattern()
    {
        // fuzzpro -pattern "[AG]-x(2)-[IV]" on MKVLAAGIVGLLLAQPSTAEEKV prints the 1-based windows 5-8 and 6-9
        var r = AnalysisTools.Fuzzpro("MKVLAAGIVGLLLAQPSTAEEKV", "[AG]-x(2)-[IV]");
        Assert.That(r.Pattern, Is.EqualTo("[AG]x(2)[IV]"));
        Assert.That(r.Hits.Select(h => (h.Position, h.Length)).ToArray(), Is.EqualTo(new[] { (4, 4), (5, 4) }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzzpro("MKV", "[BA]"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzzpro("", "A"));
    }

    [Test]
    public void Fuzztran_SearchesTheTranslatedFrames()
    {
        // fuzztran -pattern KV on ATGAAAGTTCTGGCT (frame 1 = MKVLA) prints the feature 4 9 with *start 2 *end 3 *translation KV
        var r = AnalysisTools.Fuzztran("ATGAAAGTTCTGGCT", "KV");
        var h = r.Hits.Single();
        Assert.That((h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd), Is.EqualTo((1, 2, 3, "KV", 3, 6, (int?)4, (int?)9)));
        // reverse frames carry the forward-strand location and the numbers the program prints for them (S - |frame| - 3 start + 4 / + 2, A5-7)
        var rev = AnalysisTools.Fuzztran("GGATTA", "S", "-1").Hits.Single();
        Assert.That((rev.Frame, rev.Translation, rev.Position, rev.Length, rev.EmbossStart, rev.EmbossEnd), Is.EqualTo((-1, "S", 0, 3, (int?)(-2), (int?)(-4))));        // protein position 2: 1 - 1 - 3 * 2 + 4 and 1 - 1 - 3 * 2 + 2
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzztran("ACGT", "A", "7"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzztran("ACGT", "[X]"));
    }

    /// <summary>A4-3 (B09 F153): the EMBOSS tools read the sequence like the programs (gap characters and other non-letters dropped first; positions count the cleaned sequence).</summary>
    [Test]
    public void EmbossTools_DropGapCharactersFromThePlainSequence_LikeThePrograms()
    {
        // fuzznuc -pattern ACGT on `AC-GT.AACGT` (cleaned ACGTAACGT): hits 1-4 and 6-9
        var hits = AnalysisTools.FindPatternHits("AC-GT.AACGT", "ACGT").Hits.Select(h => (h.Position, h.Length)).ToArray();
        Assert.That(hits, Is.EqualTo(new[] { (0, 4), (5, 4) }));
        // fuzztran -pattern KV on `ATG-AAA GTT.CTGGCT` is the clean-text hit
        var h = AnalysisTools.Fuzztran("ATG-AAA GTT.CTGGCT", "KV").Hits.Single();
        Assert.That((h.EmbossStart, h.EmbossEnd), Is.EqualTo(((int?)4, (int?)9)));
        Assert.That(AnalysisTools.FindTandemRepeatsEmboss("ACACAC-ACACAC.ACAC", "etandem", 2, 4, 5).Repeats.Select(r => (r.Position, r.Length)).ToArray(),
            Is.EqualTo(AnalysisTools.FindTandemRepeatsEmboss("ACACACACACACACAC", "etandem", 2, 4, 5).Repeats.Select(r => (r.Position, r.Length)).ToArray()));
        // consambig / showalign: `? # +` are the unknown residue, digits and blanks vanish
        Assert.That(AnalysisTools.Consambig(["ACGT#CGT", "ACGT?CGT"]).Consensus, Is.EqualTo(AnalysisTools.Consambig(["ACGTNCGT", "ACGTNCGT"]).Consensus));
        Assert.That(AnalysisTools.Consambig(["AC1GT CGT", "ACG_TACGT"]).Consensus, Is.EqualTo(AnalysisTools.Consambig(["ACGTCGT", "ACGTACGT"]).Consensus));
        Assert.That(AnalysisTools.SeqMatchAll(["ACGT-ACGTA", "ACGT.ACGTA"], 4).Matches, Is.EqualTo(AnalysisTools.SeqMatchAll(["ACGTACGTA", "ACGTACGTA"], 4).Matches));
    }

    /// <summary>A4-1 (B09 F157): fuzztran, showalign and consambig took no token; they now hand the request token to the ports.</summary>
    [Test]
    public void A4_1_FuzztranShowalignConsambig_StopWithACancelledToken()
    {
        var cancelled = new CancellationToken(true);
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => AnalysisTools.Fuzztran("ATGAAAGTTCTGGCT", "KV", cancellationToken: cancelled), "fuzztran");
            Assert.Catch<OperationCanceledException>(() => AnalysisTools.ShowAlign(["a", "b"], ["ACGTACGT", "ACGTACGA"], cancellationToken: cancelled), "showalign");
            Assert.Catch<OperationCanceledException>(() => AnalysisTools.Consambig(["ACGTACGT", "ACGTACGA"], cancellationToken: cancelled), "consambig");
        });
        Assert.That(AnalysisTools.Consambig(["ACGTACGT", "ACGTACGA"], cancellationToken: CancellationToken.None).Consensus, Is.EqualTo("ACGTACGW"));
    }
}
