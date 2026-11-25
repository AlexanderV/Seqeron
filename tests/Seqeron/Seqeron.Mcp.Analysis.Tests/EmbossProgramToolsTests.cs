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

    private const string EmblWithFeatures = "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             join(3..12,20..29)\nFT                   /translation=\"MKV\"\nFT                   /gene=\"g1\"\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n";

    // the index files dbxfasta 6.6.0 made for a FASTA file of two entries (the field id, compressed pages), gzip + base64 of their bytes
    private static readonly Dictionary<string, string> DbxIndexFiles = new()
    {
        ["idx/tdb.ent"] = "H4sIAAAAAAAC/1NW8CvNTUotUshPU0jLzEkttlIw5FJWCErNSU0sTrVSMNAzAHJdEkuAbCAwMNAHI67gzLz0nFSwlrzE3FSFlMSSxCSgDq5EvbRELgC3kTO+VgAAAA==",
        ["idx/tdb.pxid"] = "H4sIAAAAAAAC/z2PywrCMBBF9/MV/YRkbNR9QRAFRd24LOlUB9IHSavo12vbZLLLublnJrdPT1k6+4ragWsmD0XX9J5CmPmdApzLB4X0cAVXsr0QBSdfkU/pRsOOnROtWcORXiRALTL+LpNR5VsoSvsUhCoaMTYQZ2O6ZrmWBeaK0TgBK5LJUHRjO6TGXzA6ZxNCOLwdNxxzbeBCdc1OPvQDgBCNYRgBAAA=",
        ["idx/tdb.xid"] = "H4sIAAAAAAAC/2NkQABGIGZhIAxAaopTCw0ZfKECsVCaA00hiM+ExBeFYrBmdAAUNIJzFBlGwSgYBaNgFIyCUTA0AQAD30AfXQgAAA==",
    };

    private static Dictionary<string, string> DbxFiles()
    {
        var files = new Dictionary<string, string>();
        foreach (var (name, text) in DbxIndexFiles)
        {
            using var gz = new System.IO.Compression.GZipStream(new MemoryStream(Convert.FromBase64String(text)), System.IO.Compression.CompressionMode.Decompress);
            using var ms = new MemoryStream();
            gz.CopyTo(ms);
            files[name] = System.Text.Encoding.Latin1.GetString(ms.ToArray());
        }

        files["db/a.fa"] = ">seq1 first one\nACGTACGTACGTAAAA\n>seq2 second one\nGGGGCCCCAAAATTTT\n";
        return files;
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
}
