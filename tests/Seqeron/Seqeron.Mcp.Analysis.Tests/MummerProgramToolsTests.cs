using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>MCP wrappers of the other programs of the MUMmer 3.23 package (B09 A3-9d / A3-9e / A3-9f): the output of the installed program on small examples.</summary>
[TestFixture]
public class MummerProgramToolsTests
{
    [Test]
    public void Mummer_WritesTheMatchesOfTheProgram()
    {
        // mummer -mum -l 8 -b -s r.fa q.fa (MUMmer 3.23)
        const string reference = ">r\nacgtacgtacgtaaccggtt\n";
        const string query = ">q\nacgtacgtaacctacgtacgtacgtaaccggtta\n";
        var r = AnalysisTools.Mummer(reference, [query], "mum", minLength: 8, both: true, showString: true);
        Assert.That(r.Text, Is.EqualTo("> q\n       1        14        20\nacgtacgtacgtaaccggtt\n> q Reverse\n       1        10        13\nacgtacgtacgta\n      12         1         9\ntaaccggtt\n"));
        Assert.That(r.Messages, Does.StartWith("# reading input file \"reference.fa\" of length 20\n"));
        // mummer -maxmatch -l 4 (every maximal match, the suffix-tree order: ref 9, 5, 1 for the first query position)
        Assert.That(AnalysisTools.Mummer(reference, [">q\nacgtacgtaacc\ntacgtacgtacgtaaccggtta\n"], "maxmatch", minLength: 4).Text,
            Is.EqualTo("> q\n       9         1         5\n       5         1        12\n       1         1         9\n       1         5         5\n       8        13         6\n       4        13        10\n       1        14        20\n       1        18         9\n       1        22         5\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, [query], both: true, reverseOnly: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, [query], matchMode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer("", [query]));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, []));
        Assert.Throws<FormatException>(() => AnalysisTools.Mummer(reference, [">q\n"]));
    }

    // promer -l 5 -c 10 -o ref.fa qry.fa (MUMmer 3.23): one gene with two amino-acid substitutions and an insertion of two residues, on the reverse strand too (B09 A3-9d, F205)
    private const string PromerRef = ">ref gene\nACGTTGCATATAAAAACACTGCACGGATTTGCGGAGAGAATTCGATTCCCGTCGAATGGATCGCACACCAGAAGGGTCCAGGAGATTTCTATTATCGATTCTTTGCAGCACATCACGGAGGATTGATGCATGTACACGGCCATCCTCTAGCGAATCAAGGATTCAAGGACATGTTATATTGGGCATATATGGATCTGAACTTATCCATGGGATCCTCTGGATC\n";
    private const string PromerQry = ">qry\nTTTACAAGAATACGGCCAGAATCTGCGGAGAGTGGTCAATCCCGGTCGAATGGATAGCCCACCAGAAAGGTCCGGGAGACTTTTACTATCGATTCTTCGCGGGACACCACGGCGGTCTCATGCACGTACACGGACACCCACTTGCAAACCAGAAAAAGGGATTCAAAGATATGCTTTACTGGGCTTACATGGACCTTAACCTGTCTATGGGCAGTTCGACGTA\n";

    [Test]
    public void Promer_WritesTheDeltaFileOfTheProgram()
    {
        var r = AnalysisTools.Promer(PromerRef, PromerQry, minMatch: 5, minCluster: 10, coords: true, referenceName: "ref.fa", queryName: "qry.fa");
        Assert.That(r.Delta, Is.EqualTo("ref.fa qry.fa\nPROMER\n>ref qry 223 223\n9 218 3 218 4 4 0\n-51\n-1\n0\n217 11 217 5 22 15 8\n-21\n-1\n0\n"));
        Assert.That(r.Coords, Is.EqualTo("ref.fa qry.fa\nPROMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  [% SIM]  [% STP]  | [FRM]  [TAGS]\n"
            + "==============================================================================================================\n"
            + "       9      218  |        3      218  |      210      216  |    94.44    94.44     0.00  |  3  3  ref\tqry\n     217       11  |      217        5  |      207      213  |    69.01    78.87     5.63  | -1 -1  ref\tqry\n"));
        Assert.That(r.Cluster, Is.Null);
        Assert.That(r.Messages, Does.StartWith("1: PREPARING DATA\n2,3: RUNNING mummer AND CREATING CLUSTERS\n# reading input file \"out.aaref\" of length 448\n"));
        Assert.That(r.Messages, Does.EndWith("4: FINISHING DATA\n5: GENERATING COORDS FILE\n"));
        Assert.That(r.AminoReference, Does.StartWith(">allcontigs ref.fa\n"));
        Assert.That(AnalysisTools.Promer(PromerRef, PromerQry, minMatch: 5, minCluster: 10, delta: false, referenceName: "ref.fa", queryName: "qry.fa").Cluster,
            Does.StartWith("ref.fa qry.fa\nPROMER\n>ref qry 223 223\n"));
        Assert.That(AnalysisTools.Promer(PromerRef, PromerQry, breakLength: 20000).Messages, Does.Contain("WARNING: Invalid break length 20000, ignoring"));
        Assert.Throws<FormatException>(() => AnalysisTools.Promer(PromerRef, PromerQry, minMatch: 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Promer(PromerRef, PromerQry, matchMode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Promer("", PromerQry));
        Assert.Throws<FormatException>(() => AnalysisTools.Promer("acgt\n", PromerQry));
    }

    [Test]
    public void ShowAlignsAndShowSnps_ReadPromerDeltaFiles()
    {
        const string delta = "ref.fa qry.fa\nPROMER\n>ref qry 223 223\n9 218 3 218 4 4 0\n-51\n-1\n0\n217 11 217 5 22 15 8\n-21\n-1\n0\n";
        var aligns = AnalysisTools.ShowAligns(delta, PromerRef, PromerQry, "ref", "qry", screenWidth: 60);
        Assert.That(aligns.Text, Does.StartWith("ref.fa qry.fa\n\n============================================================\n-- Alignments between ref and qry\n\n-- BEGIN alignment [ +3 9 - 218 | +3 3 - 218 ]\n\n\n"
            + "9          YKNTARICGENSIPVEWIAHQKGPGDFYYRFFAAHHGGLMHVHGHPLAN\n           YKNTARICGE SIPVEWIAHQKGPGDFYYRFFA HHGGLMHVHGHPLAN\n3          YKNTARICGEWSIPVEWIAHQKGPGDFYYRFFAGHHGGLMHVHGHPLAN\n\n"
            + "156        Q..GFKDMLYWAYMDLNLSMGSS\n           Q  GFKDMLYWAYMDLNLSMGSS\n150        QKKGFKDMLYWAYMDLNLSMGSS\n\n\n--   END alignment [ +3 9 - 218 | +3 3 - 218 ]\n"));
        Assert.That(aligns.Text, Does.Contain("217        EDPMDKFRSIYAQYNMSLNP..*FARGWPCTCINPP*CAAKNR**KSPG\n           E P+D+ RS+ AQ ++SLNP   FA G PCTC+ PP C AKNR**KSPG\n217        ELPIDRLRSM*AQ*SISLNPFFWFASGCPCTCMRPPWCPAKNR**KSPG\n"));
        Assert.That(AnalysisTools.ShowAligns(delta, PromerRef, PromerQry, "ref", "qry", matrixType: 7).Messages, Is.EqualTo("WARNING: invalid matrix type 7, using default\n"));
        var snps = AnalysisTools.ShowSnps(delta, PromerRef, PromerQry, context: 2);
        Assert.That(snps.Text, Does.StartWith("ref.fa qry.fa\nPROMER\n\n    [P1]  [SUB]  [P2]      |   [BUFF]   [DIST]  |  [R]  [Q]  |  [CTX R]  [CTX Q]  | [FRM]  [TAGS]\n"
            + "=====================================================================================================\n"
            + "      25   R L   19        |       15       19  |    1    1  |    QIRAV    QILAV  | -1 -1  ref\tqry\n      39   N W   33        |       31       33  |    1    1  |    GENSI    GEWSI  |  3  3  ref\tqry\n"));
    }

    // gaps g1.fa -r, delta2blocks, mummer-annotate and run-mummer1 of MUMmer 3.23 on the three-substitution example (B09 A3-9e, F206)
    private const string Ref1 = ">ref\nGCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";
    private const string Qry1 = ">qry\nAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTGGGGAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCA\n";

    [Test]
    public void MummerGaps_PrintsTheConsistentAndTheOtherMatches()
    {
        Assert.That(AnalysisTools.MummerGaps("       3         9        20\n      40        50        15\n      70        80        30\n      30        20         5\n", "g1.fa", true).Text,
            Is.EqualTo("> g1.fa reverse Consistent matches\n       3        9     20    none      -      -\n      40       50     15    none     17     21\n      70       80     30    none     15     15\n> g1.fa reverse Other matches\n      30       20      5    none      -      -\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MummerGaps(" "));
        Assert.Throws<FormatException>(() => AnalysisTools.MummerGaps("> q\n1 2 3\n"));
    }

    [Test]
    public void Delta2Blocks_PrintsTheCigarOfEveryAlignment()
    {
        const string delta = "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n";
        Assert.That(AnalysisTools.Delta2Blocks(delta, Ref1, Qry1, "ref", "qry").Text, Is.EqualTo("ref qry 6 140 1 135 65M2X1M1X66M\n"));
        Assert.Throws<NotSupportedException>(() => AnalysisTools.Delta2Blocks("ref.fa qry.fa\nPROMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", Ref1, Qry1, "ref", "qry"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Delta2Blocks("", Ref1, Qry1, "ref", "qry"));
    }

    [Test]
    public void MummerAnnotate_AndRunMummer1_AlignTheGaps()
    {
        const string gaps = "> ref.fa  Consistent matches\n       6        1     65    none      -      -\n      75       70     66    none      4      4\n> ref.fa  Other matches\n";
        var r = AnalysisTools.MummerAnnotate(gaps, Qry1, new Dictionary<string, string> { ["ref.fa"] = Ref1 });
        Assert.That(r.Text, Is.EqualTo(gaps.Replace("      75       70     66    none      4      4\n", "      75       70     66    none      4      4\n    Errors = 3\nT:  cttaagggttaagtaagtgtgatg\nS:  cttaagggttggggaagtgtgatg\n              ^^ ^          \n", StringComparison.Ordinal)));
        Assert.That(r.WithErrors, Is.EqualTo("> ref.fa  Consistent matches\n       6        1     65    none      -      -\n      75       70     66    none      4      4       3\n> ref.fa  Other matches\n"));
        Assert.Throws<NotSupportedException>(() => AnalysisTools.MummerAnnotate("> ref.fa x\n   9   9  30  none  2  2\n", Qry1, new Dictionary<string, string> { ["ref.fa"] = Ref1 }));
        Assert.Throws<FormatException>(() => AnalysisTools.RunMummer1(Ref1, Qry1, reverse: false));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RunMummer1("", Qry1));
    }

    [Test]
    public void Delta2Maf_WritesTheMafOfTheAlignments()
    {
        const string delta = "ref.fa qry.fa\nNUCMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n";
        string text = AnalysisTools.Delta2Maf(delta, Ref1, Qry1).Text;
        Assert.That(text, Is.EqualTo("##maf version=1 scoring=single_cov2\na score=134\ns ref 5 135 + 140 " + Ref1.Split('\n')[1].ToLowerInvariant()[5..140]
            + "\ns qry 0 135 + 135 " + Qry1.Split('\n')[1].ToLowerInvariant() + "\n\n##eof maf"));
        Assert.Throws<FormatException>(() => AnalysisTools.Delta2Maf(delta, Ref1, Qry1, "nonexistent"));
        Assert.Throws<NotSupportedException>(() => AnalysisTools.Delta2Maf("ref.fa qry.fa\nPROMER\n>ref qry 140 135\n6 140 1 135 3 3 0\n0\n", Ref1, Qry1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Delta2Maf("", Ref1, Qry1));
    }

    [Test]
    public void Nucmer2Xfig_PlotsACoordsFile()
    {
        // nucmer2xfig c.coords (MUMmer 3.23, B09 A3-9e F208)
        string text = AnalysisTools.Nucmer2Xfig("  2551   2577 | 240 266 | 27 27 | 96.30 | 30000 1424 | 0.09 1.90 | 2R ctg1\n  2600   9000 | 300 6700 | 6401 6401 | 99.00 | 30000 1424 | 21.33 99.00 | 2R ctg1\n").Text;
        Assert.That(text, Does.StartWith("#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nSingle\n-2\n1200 2\n2 1 0 2 2 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 12 926 13 926\n2 1 0 2 2 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 13 980 45 980\n"
            + "2 1 0 2 1 7 100 0 -1 0.000 0 0 -1 0 0 2\n\t 13 1000 45 1000\n4 0 0 100 0 0 12 0.0000 4 135 405 150 0 2R\\001\n"));
        Assert.That(text, Does.EndWith("4 0 0 100 0 0 5 0.0000 4 135 405 -150 1010 100\\001\n2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 0 0 2\n-50 1000 0 1000\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Nucmer2Xfig(null!));
    }

    [Test]
    public void Nucmer_Prefix_NamesTheFilesTheScriptLeavesBehind()
    {
        // nucmer -p pre -l 15 -o ref.fa qry.fa leaves pre.coords and pre.delta (the temporary pre.ntref and pre.mgaps are deleted) and writes the progress lines to standard error (B09 A3-9f, F209)
        var r = AnalysisTools.Nucmer(Ref1, Qry1, minMatch: 15, coords: true, prefix: "pre", referenceName: "/w/ref.fa", queryName: "/w/qry.fa");
        Assert.That(r.Files!.Keys, Is.EquivalentTo(new[] { "pre.delta", "pre.coords" }));
        Assert.That(r.Files["pre.delta"], Is.EqualTo(r.Delta));
        Assert.That(r.Log, Is.EqualTo("1: PREPARING DATA\n2,3: RUNNING mummer AND CREATING CLUSTERS\n# reading input file \"pre.ntref\" of length 141\n# construct suffix tree for sequence of length 141\n"
            + "# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n# reading input file \"/w/qry.fa\" of length 135\n# matching query-file \"/w/qry.fa\"\n# against subject-file \"pre.ntref\"\n4: FINISHING DATA\n5: GENERATING COORDS FILE\n"));
        Assert.That(AnalysisTools.Nucmer(Ref1, Qry1, minMatch: 15, delta: false, prefix: "x").Files!.Keys, Is.EquivalentTo(new[] { "x.cluster" }));
    }

    // combineMUMs -x ref.fa qry.fa pf.gaps and run-mummer3 ref.fa qry.fa pf (MUMmer 3.23): one query with a substitution and two single-base indels (B09 A3-9e2, F210)
    private const string Cm3Ref = ">ref demo\nggcccagtccagatcctcggaagtcccaactaacgaataagtagatccttctaaatagta\ngtatacgaatgtaatcaacatattaccgcacccaactttgattgggtcataaacaaacat\ncatgtttctggaggcggtttgtgagaccgacaaacgcgatagagtaatcttataaatcag\ngcaaatactcggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttc\ntttgaggttaatgggaggca\n";
    private const string Cm3Qry = ">qry demo\nagatcctcggaagtcccaactaacgaataagtagatccttctaaatagtagtatacgaat\nataatcaacatattaccgcacccaactttgattgggtcataaacaaacatcatgtttctg\ngaggcggttttgagaccgacaaacgcgatagagtaatcttataaatcaggcaaatactcg\nggtcctttctacgggtcgcggccgcgtggtcataggcagttttatgtttctttgaggtta\n";

    [Test]
    public void CombineMums_AlignsTheGapsOfTheProgram()
    {
        var r = AnalysisTools.CombineMums(Cm3Ref, Cm3Qry, "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n", noCoverFiles: true, referenceName: "ref.fa", queryName: "qry.fa");
        Assert.That(r.Text, Is.EqualTo("> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     Errors = 1\nA: gtatacgaatgtaatcaacat\nB: gtatacgaatataatcaacat\n   ==========^==========\n     142      131     51    none      1      0\n     Errors = 1\nA: gaggcggtttgtgagaccgac\nB: gaggcggttt-tgagaccgac\n   ==========^==========\n     193      183     58      -2      0      1\n     Errors = 1\nA: aaatactcgg-tcctttctac\nB: aaatactcgggtcctttctac\n   ==========^==========\nRegion:         11 .. 250                1 .. 240              3 / 240         1.25%\n> qry Reverse\n"));
        Assert.That(r.ErrorsFile, Is.EqualTo("> qry\n      11        1     60    none      -      -       -\n      72       62     69    none      1      1       1\n     142      131     51    none      1      0       1\n     193      183     58      -2      0      1       1\n> qry Reverse\n"));
        Assert.That(r.Messages, Is.EqualTo("           Ref len = 260\n            acgt's = 260\n        Non acgt's = 0\n Number of matches = 1\nSum of match bases = 240\n   Avg match bases = 240\n"));
        Assert.That(r.Files.Keys, Is.EquivalentTo(new[] { "witherrors.gaps" }));
        // without -x the cover files are named after the files: <ref base name>Ref.cover, <query base name>Query.cover
        Assert.That(AnalysisTools.CombineMums(Cm3Ref, Cm3Qry, "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n", referenceName: "ref.fa", queryName: "qry.fa", errorsFileName: "my.errors").Files.Keys, Is.EquivalentTo(new[] { "my.errors", "refRef.cover", "qryQuery.cover" }));
        Assert.That(AnalysisTools.CombineMums(Cm3Ref, Cm3Qry, "> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n", onlyDifferences: true, noCoverFiles: true).Text, Does.Not.Contain("Region:"));
        Assert.Throws<FormatException>(() => AnalysisTools.CombineMums(Cm3Ref, Cm3Qry, "> a\n> a Reverse\n> b\n", noCoverFiles: true));    // more forward blocks than query records
        Assert.Throws<FormatException>(() => AnalysisTools.CombineMums(Cm3Ref, Cm3Qry, "> a\nnot numbers\n", noCoverFiles: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CombineMums("", Cm3Qry, ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CombineMums(Cm3Ref, "", ""));
    }

    [Test]
    public void RunMummer3_RunsTheThreeProgramsOfTheScript()
    {
        var r = AnalysisTools.RunMummer3(Cm3Ref, Cm3Qry, "ref.fa", "qry.fa", "pf");
        Assert.That(r.Out, Is.EqualTo("> qry\n      11         1        60\n      72        62        69\n     142       131        51\n     191       181        60\n> qry Reverse\n"));
        Assert.That(r.Gaps, Is.EqualTo("> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     142      131     51    none      1      0\n     193      183     58      -2      0      1\n> qry Reverse\n"));
        Assert.That(r.Align, Is.EqualTo("> qry\n      11        1     60    none      -      -\n      72       62     69    none      1      1\n     Errors = 1\nA: gtatacgaatgtaatcaacat\nB: gtatacgaatataatcaacat\n   ==========^==========\n     142      131     51    none      1      0\n     Errors = 1\nA: gaggcggtttgtgagaccgac\nB: gaggcggttt-tgagaccgac\n   ==========^==========\n     193      183     58      -2      0      1\n     Errors = 1\nA: aaatactcgg-tcctttctac\nB: aaatactcgggtcctttctac\n   ==========^==========\nRegion:         11 .. 250                1 .. 240              3 / 240         1.25%\n> qry Reverse\n"));
        Assert.That(r.ErrorsGaps, Is.EqualTo("> qry\n      11        1     60    none      -      -       -\n      72       62     69    none      1      1       1\n     142      131     51    none      1      0       1\n     193      183     58      -2      0      1       1\n> qry Reverse\n"));
        Assert.That(r.Stdout, Is.EqualTo("Find MUMs\nDetermine gaps\nAlign gaps\n"));
        Assert.That(r.Messages, Is.EqualTo("# reading input file \"ref.fa\" of length 260\n# construct suffix tree for sequence of length 260\n# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n# reading input file \"qry.fa\" of length 240\n# matching query-file \"qry.fa\"\n# against subject-file \"ref.fa\"\n           Ref len = 260\n            acgt's = 260\n        Non acgt's = 0\n Number of matches = 1\nSum of match bases = 240\n   Avg match bases = 240\n"));
        Assert.Throws<FormatException>(() => AnalysisTools.RunMummer3(Cm3Ref, ">q\n"));      // mummer: empty sequence in multiple fasta file; sh -e stops the script
        Assert.Throws<ArgumentException>(() => AnalysisTools.RunMummer3("", Cm3Qry));
    }

    // repeat-match -n 12 t.fa and exact-tandems t.fa 6 (MUMmer 3.23) (B09 A3-9e3, F211)
    private const string TandemFasta = ">tandem demo\nacactgctcactccaaccccggcccctgagtccgaggagaggatcacagtctggatcaca\ngtctggatcacagtctggatcacagtctggatcacagtctgggtgcttcagagtatgtat\naccactgggtggatcaggatcaggatcaaggatacggcggagggcacg\n";

    [Test]
    public void RepeatMatch_ListsTheRepeatsInTheOrderOfTheProgram()
    {
        var r = AnalysisTools.RepeatMatch(TandemFasta, minLength: 12);
        Assert.That(r.Text, Is.EqualTo("Long Exact Matches:\n   Start1     Start2    Length\n       41         53        50\n       41         65        38\n       41         77        26\n       41         89        14\n      131        137        12\n"));
        Assert.That(r.Messages, Is.EqualTo("Genome Length = 168   Used 248 internal nodes\n"));
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(AnalysisTools.RepeatMatch(TandemFasta, minLength: 12, tandemOnly: true).Text, Does.StartWith("Long Exact Matches:\n   Start1     Start2    Length\n       41         53        50\n"));
        Assert.That(AnalysisTools.RepeatMatch(">pal\nacgtagcta\n", minLength: 3).Text, Is.EqualTo("Long Exact Matches:\n   Start1     Start2    Length\n        4          9r        6\n        1          4r        4\n"));
        Assert.That(AnalysisTools.RepeatMatch(">p\nacgtacgt\n", minLength: 3).ExitCode, Is.EqualTo(255));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RepeatMatch("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RepeatMatch(TandemFasta, minLength: 0));
        Assert.Throws<FormatException>(() => AnalysisTools.RepeatMatch("acgtacgtacgt\n"));      // no header line: the program works on an uninitialised buffer
    }

    [Test]
    public void ExactTandems_WritesTheTableOfTheIntendedPipeline_AndTheInstalledScriptFails()
    {
        var r = AnalysisTools.ExactTandems(TandemFasta, 6);
        Assert.That(r.Text, Is.EqualTo("   Start   Extent  UnitLen     Copies\n      41       62       12        5.2\n     131       18        6        3.0\n"));
        Assert.That(r.Messages, Is.EqualTo("Genome Length = 168   Used 129 internal nodes\n"));
        Assert.That(r.ExitCode, Is.Zero);
        var installed = AnalysisTools.ExactTandems(TandemFasta, 6, installedScript: true);
        Assert.That(installed.Text, Is.EqualTo("Finding matches\nTandem repeats\n"));
        Assert.That(installed.ExitCode, Is.EqualTo(2));
        Assert.That(installed.Messages, Does.EndWith("awk: fatal: cannot open source file `$' for reading: No such file or directory\n"));
        Assert.That(installed.FilesLeft, Is.EqualTo(new[] { "$$.tmp.matches" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ExactTandems("", 6));
    }

    // mapview -n 1 -p demo x.coords (MUMmer 3.23) (B09 A3-9e4, F212)
    private const string DemoCoords = "ref.fa qry.fa\nNUCMER\n\n    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  | [TAGS]\n==========================================================================================================\n      11      250  |        1      240  |      240      240  |    98.76  |      260      240  | ref\tqry\n";

    [Test]
    public void Mapview_DrawsTheFigFilesOfTheScript()
    {
        var r = AnalysisTools.Mapview(DemoCoords, numberOfFiles: 1, prefix: "demo");
        Assert.That(r.Figs.Keys, Is.EquivalentTo(new[] { "demo_0.fig" }));
        Assert.That(r.Figs["demo_0.fig"], Is.EqualTo("#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nMultiple\n-2\n1200 2\n2 1 0 10 11 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 0 3050 17 3050\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 5250 17 5250\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 4250 17 4250\n2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n\t0 3250 17 3250\n2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n\t0 5250 0 5300\n4 0 0 100 0 0 8 0.0000 4 135 405 \t 0 5400 0k\\001\n4 0 0 50 0 0 14 0.0000 4 135 450 0 3000 0\\001\n4 0 0 50 0 0 14 0.0000 4 135 810 -233 3000 259\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 5270 50%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 5270 50%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 4270 75%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 4270 75%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t-200 3270 100%\\001\n4 0 0 100 0 0 8 0.0000 4 135 405 \t37 3270 100%\\001\n2 1 0 2 4 0 40 0 -1 0.000 0 0 -1 0 0 2\n\t0 3330 17 3330\n2 1 0 5 20 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 0 3150 17 3150\n4 0 0 100 0 0 8 0.0000 4 135 405 \t100 6100  Legend \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6300 99 6300\n2 1 0 5 27 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6300 200 6300\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6300 230 6300\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6330  5' utr \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6450 99 6450\n2 1 0 5 2 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6450 200 6450\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6450 230 6450\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6480  cds \\001\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 70 6600 99 6600\n2 1 0 5 1 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 100 6600 200 6600\n2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t 200 6600 230 6600\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6630  3' utr \\001\n2 1 0 2 4 0 50 0 -1 0.000 0 0 -1 0 0 2\n\t100 6750 200 6750\n4 0 0 100 0 0 8 0.0000 4 135 405 \t300 6780  match found by NUCMER \\001\n"));
        Assert.That(r.Stdout, Is.Empty);
        Assert.That(r.ExitCode, Is.Zero);
        var pdf = AnalysisTools.Mapview(DemoCoords, numberOfFiles: 1, prefix: "demo", format: "pdf");
        Assert.That(pdf.Commands, Is.EqualTo(new[] { "fig2dev -L pdf -x 100 -M tmpdemo_0.fig demo_0.pdf" }));
        Assert.That(pdf.Figs.Keys, Is.EquivalentTo(new[] { "tmpdemo_0.fig" }));
        Assert.That(AnalysisTools.Mapview(DemoCoords, "other\tx\tsingle-exon\t5\t50\t.\t+\t.\tg1\n").ExitCode, Is.EqualTo(1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mapview("", numberOfFiles: 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mapview(DemoCoords, format: "png"));
    }
}
