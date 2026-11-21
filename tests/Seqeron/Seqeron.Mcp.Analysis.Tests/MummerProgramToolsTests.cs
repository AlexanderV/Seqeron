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
}
