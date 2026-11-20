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
        Assert.Throws<ArgumentException>(() => AnalysisTools.ShowAlign(["a", "b"], ["AC", "AC"], uppercase: "x-y"));
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
        // reverse frames carry the forward-strand location and no program coordinates
        var rev = AnalysisTools.Fuzztran("GGATTA", "S", "-1").Hits.Single();
        Assert.That((rev.Frame, rev.Translation, rev.Position, rev.Length, rev.EmbossStart), Is.EqualTo((-1, "S", 0, 3, (int?)null)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzztran("ACGT", "A", "7"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Fuzztran("ACGT", "[X]"));
    }
}
