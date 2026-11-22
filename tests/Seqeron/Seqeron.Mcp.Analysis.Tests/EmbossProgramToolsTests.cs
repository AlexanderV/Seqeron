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
        // the special qualifiers (A5-8c3c, F219): a valid value is kept, a bad one is a note, a correctable one is corrected (the program also warns on stderr)
        string special = genbank.Replace("/translation=\"MKV\"", "/collection_date=\"2020\"\n                     /lat_lon=\"12.5 N 33.1\"\n                     /replace=\"a c g\"", StringComparison.Ordinal);
        var s2 = AnalysisTools.Diffseq(special, embl, 5, false);
        Assert.That(s2.Report, Does.Contain("Feature: CDS 25-35 note='a note' collection_date='2020' note='*lat_lon: 12.5 N 33.1' replace='acg'\n"));
    }

    [Test]
    public void Featcopy_ReadsAnyFormatAndWritesTheOthers_AsTheProgram()
    {
        // featcopy 6.6.0 (A5-8c3e, F241), dates as of the run: a GFF3 file with Name / Parent / product / percent-encoded values written as EMBL, GenBank, GFF2 (the sub-feature is written as its own feature)
        const string gff = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc;product=foo\nseq1\tsrc\tmRNA\t20\t40\t5.5\t-\t.\tID=m1;Parent=g1;note=hello%20world\n";
        var embl = AnalysisTools.Featcopy(gff, null, "embl");
        Assert.That(embl.Text, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   gene            10..50\nFT                   /note=\"*Name abc\"\nFT                   /product=\"foo\"\nFT   mRNA            complement(20..40)\nFT                   /note=\"hello world\"\n"));
        Assert.That(embl.Tables, Is.EqualTo(1));
        Assert.That(AnalysisTools.Featcopy(gff, null, "gb").Text, Does.StartWith("FEATURES             Location/Qualifiers\n     gene            10..50\n                     /note=\"*Name abc\"\n"));
        var gff2 = AnalysisTools.Featcopy(gff, null, "gff2");
        Assert.That(gff2.Text, Does.EndWith("##Type unknown <> seq1\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tSequence \"seq1.2\" ; product \"foo\"\n"));
        // the format:: prefix forces a reader: SwissProt feature lines, PIR feature lines (a text with no feature of the format is refused as the program refuses it)
        const string swiss = "ID   KIN1_HUMAN   Reviewed;   120 AA.\nAC   P12345;\nFT   DOMAIN         5     50       Kinase domain.\nFT   ACT_SITE     30     30       Proton acceptor (By similarity).\n//\n";
        var sw = AnalysisTools.Featcopy(swiss, "swiss::features.dat", "gff3");
        Assert.That(sw.Text, Does.EndWith("#!Type Protein\n#!Source-version EMBOSS 6.6.0.0\nKIN1_HUMAN\tSWISSPROT\tpolypeptide_domain\t5\t50\t.\t.\t.\tID=KIN1_HUMAN.1;note=Kinase domain\nKIN1_HUMAN\tSWISSPROT\tactive_site_residue\t30\t30\t.\t.\t.\tID=KIN1_HUMAN.2;note=Proton acceptor;comment=By similarity\n"));
        Assert.That(AnalysisTools.Featcopy(swiss, "swiss::features.dat", "gff2").Text, Does.Contain("Sequence \"KIN1_HUMAN.2\" ; note \"Proton acceptor\" ; comment \"By similarity\"\n"));
        Assert.That(() => AnalysisTools.Featcopy("F;1-10/Domain: x\n", "gff2::features.dat", "gff3"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        // an unknown format name is an error message of the library, the formats are tried as for an unknown input
        var bogus = AnalysisTools.Featcopy(gff, "bogus::features.dat", "gff3");
        Assert.That(bogus.Messages, Does.Contain("Error: Unknown input format 'bogus'"));
        Assert.That(bogus.Tables, Is.EqualTo(1));
        // the debug trace and a range (-fbegin / -fend cut the table: positions shift, a feature cut at the end is flagged)
        var cut = AnalysisTools.Featcopy(gff, null, "debug", 10, 30);
        Assert.That(cut.Text, Does.StartWith("  Format: 0 (gff3)  Start: 10\n  End: 30\n  Len: 21\n  Offset: 0\n  Groups: 4\n"));
        Assert.That(cut.Text, Does.Contain("    Start: 1\n    End: 21\n").And.Contain("      END_AFTER_SEQ\n"));
        // a format the program refuses for the kind of sequence: an error message of the library and nothing written (embl for a protein table)
        var refused = AnalysisTools.Featcopy(swiss, "swiss::features.dat", "embl");
        Assert.That(refused.Messages, Does.Contain("Error: Feature format 'embl' not supported for protein sequence features"));
        Assert.That(refused.Tables, Is.EqualTo(0));
        // the other writers of the feature library (A5-8c3f, F244): bed with its blocks, Swiss-Prot FT lines, PIR, DAS GFF, draw
        const string gb = "LOCUS       SEQA 300 bp    DNA     linear   UNK 01-JAN-2000\nFEATURES             Location/Qualifiers\n     gene            complement(10..50)\n                     /gene=\"abc\"\n     CDS             join(60..80,90..120)\n                     /product=\"foo\"\nORIGIN\n        1 acgt\n//\n";
        Assert.That(AnalysisTools.Featcopy(gb, null, "bed").Text, Is.EqualTo("track name= description=\"EMBOSS:6.6.0.0\" useScore=0\n\t9\t49\tgene\t0\t-\t9\t49\t0\t1\t0\t41\n\t59\t119\tCDS\t0\t+\t59\t119\t0\t2\t21,31\t0,30\n"));
        Assert.That(AnalysisTools.Featcopy(gb, null, "draw").Text, Does.Contain("Range 10 50 | | 9 H\ngene\nendlabel\n").And.EndWith("endgroup\n"));
        Assert.That(AnalysisTools.Featcopy(swiss, "swiss::features.dat", "swiss").Text, Is.EqualTo("FT   DOMAIN        5     50       Kinase domain.\nFT   ACT_SITE     30     30       Proton acceptor (By similarity).\n"));
        Assert.That(AnalysisTools.Featcopy(swiss, "swiss::features.dat", "pir").Text, Is.EqualTo("F;5-50/Domain: Kinase domain\nF;30/Active site: Proton acceptor #By similarity\n"));
        Assert.That(AnalysisTools.Featcopy(gb, null, "dasgff").Text, Does.Contain("<TYPE id=\"SO:0000704\" category=\"other\" reference=\"no\" subparts=\"no\">gene</TYPE>"));
        Assert.That(() => AnalysisTools.Featcopy(" ", null, "gff3"), Throws.ArgumentException);
    }

    [Test]
    public void Featcopy_ListFilesAsFeatureInput_AsTheProgram()
    {
        // featcopy 6.6.0 -features @list (A5-8c3f, F247): the words of a list file are queries, a '#' starts a comment, a missing file is an error message and the list goes on
        const string one = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1\n";
        const string two = "##gff-version 3\n##sequence-region seq2 1 100\nseq2\tsrc\texon\t20\t30\t.\t+\t.\tID=e1\n";
        var files = new Dictionary<string, string> { ["a.gff"] = one, ["b.gff"] = two, ["l.txt"] = "# the files\na.gff\nmissing.gff\nb.gff # last\n" };
        var bed = AnalysisTools.Featcopy(one, "@l.txt", "bed", 1, null, default, files);
        Assert.That(bed.Text, Is.EqualTo("track name=seq1 description=\"EMBOSS:6.6.0.0\" useScore=0\nseq1\t9\t49\tgene\t0\t+\t9\t49\t0\t1\t0\t41\ntrack name=seq2 description=\"EMBOSS:6.6.0.0\" useScore=0\nseq2\t19\t29\texon\t0\t+\t19\t29\t0\t1\t0\t11\n"));
        var all = AnalysisTools.Featcopy(one, "@l.txt", "gff3", 1, null, default, files);
        Assert.That(all.Tables, Is.EqualTo(2));
        Assert.That(all.Messages, Does.Contain("Error: Failed to open filename 'missing.gff'").And.Contain("Error: Failed to read data data 'missing.gff'"));
        Assert.That(all.Text, Does.Contain("seq1\tsrc\tgene\t10\t50").And.Contain("seq2\tsrc\texon\t20\t30"));
        Assert.That(() => AnalysisTools.Featcopy(one, "@nolist.txt", "gff3"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        Assert.That(() => AnalysisTools.Featcopy(one, "asis::ACGT", "gff3", 1, null, default, files), Throws.ArgumentException.With.Message.Contains("asis:: not supported by datatype"));
    }

    [Test]
    public void Featmerge_CombinesTwoFeatureTables_AsTheProgram()
    {
        // featmerge 6.6.0 (A5-8c3f, F245): the features of the second table are copied into the first, which is written in the output format
        const string a = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc\n";
        const string b = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\texon\t20\t30\t.\t+\t.\tID=e1\n";
        var merged = AnalysisTools.Featmerge(a, b);
        Assert.That(merged.Tables, Is.EqualTo(1));
        Assert.That(merged.Text, Does.Contain("seq1\tsrc\tgene\t10\t50").And.Contain("seq1\tsrc\texon\t20\t30"));
        Assert.That(AnalysisTools.Featmerge(a, b, null, null, "bed").Text, Does.StartWith("track name=seq1 description=\"EMBOSS:6.6.0.0\" useScore=0\n"));
        Assert.That(() => AnalysisTools.Featmerge(a, b, "missing.dat"), Throws.ArgumentException.With.Message.Contains("Bad value for '-afeatures'"));
        Assert.That(() => AnalysisTools.Featmerge(" ", b), Throws.ArgumentException);
    }

    [Test]
    public void Feattext_ReturnsTheStoredEntryText_AsTheProgram()
    {
        // feattext 6.6.0 (A5-8c3f, F246): the EMBL reader keeps the lines of the entry, the PIR reader none
        const string embl = "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   gene            1..30\nXX\n//\n";
        var r = AnalysisTools.Feattext(embl);
        Assert.That(r.Tables, Is.EqualTo(1));
        Assert.That(r.Text, Is.EqualTo(embl));
        Assert.That(AnalysisTools.Feattext("F;1-10/Domain: x\n", "pir::features.dat").Text, Is.EqualTo(""));
        Assert.That(() => AnalysisTools.Feattext("nothing here\n"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        Assert.That(() => AnalysisTools.Feattext(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqret_ReadsTheSequenceFormats_AsTheProgram()
    {
        // seqret 6.6.0 (A5-8c3g, F248): the DE description of a Swiss-Prot entry, its features, the type letters of a FASTA file, an asis:: sequence
        const string sw = "ID   ABC_HUMAN   Reviewed;   10 AA.\nAC   P12345;\nDE   RecName: Full=Kinase alpha; Short=KA;\nDE   AltName: Full=Other name;\nFT   DOMAIN       1      5       Kinase domain.\nSQ   SEQUENCE   10 AA;  1000 MW;  0 CRC64;\n     MKVLAAGIVG\n//\n";
        var r = AnalysisTools.Seqret(sw, "swiss", feature: true);
        Assert.That(r.SequenceText, Is.EqualTo(">ABC_HUMAN P12345 Kinase alpha (KA) (Other name)\nMKVLAAGIVG\n"));
        Assert.That(r.Sequences, Is.EqualTo(1));
        Assert.That(r.FeatureText, Does.Contain("ABC_HUMAN\tSWISSPROT\tpolypeptide_domain\t1\t5\t.\t.\t.\tID=ABC_HUMAN.1;note=Kinase domain"));
        var fa = AnalysisTools.Seqret(">s1 first\nACGTXX\n>s2\nMKVL*\n");
        Assert.That(fa.SequenceText, Is.EqualTo(">s1 first\nACGTNN\n>s2\nMKVL*\n"));
        Assert.That(AnalysisTools.Seqret("asis::ACGT.ACGT").SequenceText, Is.EqualTo(">asis\nACGT-ACGT\n"));
        Assert.That(AnalysisTools.Seqret(">a\nACGT\n>b\nTT\n", firstOnly: true).Sequences, Is.EqualTo(1));
        Assert.That(AnalysisTools.Seqret("CLUSTAL W (1.83) multiple sequence alignment\n\n\nseq1      ACGT\nseq2      ACGA\n\n").SequenceText, Is.EqualTo(">seq1\nACGT\n>seq2\nACGA\n"));
        Assert.That(() => AnalysisTools.Seqret("AS 1 1\n"), Throws.ArgumentException);          // the ACE reader of the library dies by a segmentation fault
        Assert.That(() => AnalysisTools.Seqret(" "), Throws.ArgumentException);
    }

    [Test]
    public void Maskfeat_MasksTheFeatureRegions_AsTheProgram()
    {
        // maskfeat 6.6.0 (A5-8c3g, F249): the default type repeat*, a type list, a mask character, lower case
        const string embl = "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   repeat_region   3..9\nFT   gene            12..18\nXX\nSQ   Sequence 20 BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n     aacgtacgta cgtacgtacg                                                20\n//\n";
        Assert.That(AnalysisTools.Maskfeat(embl, "embl").SequenceText, Is.EqualTo(">SEQA SEQA.1\naaNNNNNNNacgtacgtacg\n"));
        Assert.That(AnalysisTools.Maskfeat(embl, "embl", "repeat* gene", maskChar: "Q").SequenceText, Is.EqualTo(">SEQA SEQA.1\naaQQQQQQQacQQQQQQQcg\n"));
        Assert.That(AnalysisTools.Maskfeat(embl, "embl", "gene", toLower: true, qualifiers: "-supper1").SequenceText, Is.EqualTo(">SEQA SEQA.1\nAACGTACGTACgtacgtaCG\n"));
        Assert.That(() => AnalysisTools.Maskfeat(embl, "embl", maskChar: "AB"), Throws.ArgumentException);
        Assert.That(() => AnalysisTools.Maskfeat(" "), Throws.ArgumentException);
    }

    [Test]
    public void Descseq_ReplacesOrAppendsNameAndDescription_AsTheProgram()
    {
        // descseq 6.6.0 (A5-8c3g, F250)
        const string fa = ">seq1 old words\nACGTACGT\n>seq2\nTTTT\n";
        Assert.That(AnalysisTools.Descseq(fa, name: "New", description: "new text").SequenceText, Is.EqualTo(">New new text\nACGTACGT\n"));
        Assert.That(AnalysisTools.Descseq(fa, name: "_x", description: " more", append: true).SequenceText, Is.EqualTo(">seq1_x old words more\nACGTACGT\n"));
        Assert.That(AnalysisTools.Descseq(fa).SequenceText, Is.EqualTo(">seq1 old words\nACGTACGT\n"));
        Assert.That(() => AnalysisTools.Descseq(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqret_WritesTheOutputFormats_AsTheProgram()
    {
        // seqret 6.6.0 sequence writers (A5-8c3h, F252): the entries are the ones of the installed program
        const string fa = ">s1 desc one\nACGTACGTAC\n>s2\nTTGGCCAATT\n";
        Assert.That(AnalysisTools.Seqret(fa, outputFormat: "pir").SequenceText, Is.EqualTo(">D1;s1\ndesc one, 10 bases\n ACGTACGTAC*\n\n>D1;s2\n, 10 bases\n TTGGCCAATT*\n\n"));
        Assert.That(AnalysisTools.Seqret(">p1 prot\nMKVLAAGIVG\n", outputFormat: "gcg").SequenceText, Is.EqualTo("!!AA_SEQUENCE 1.0\n\nprot\n\np1  Length: 10  Type: P  Check: 4069 ..\n\n  1 MKVLAAGIVG\n\n"));
        var single = AnalysisTools.Seqret(fa, outputFormat: "gcg", outputSingle: true, outputDirectory: "d", outputExtension: "X1");
        Assert.That(single.SingleFiles!.Keys, Is.EquivalentTo(new[] { "d/s1.x1", "d/s2.x1" }));
        Assert.That(() => AnalysisTools.Seqret(fa, outputFormat: "embl"), Throws.ArgumentException);          // a writer this port does not have yet
    }

    [Test]
    public void Seqret_ReadsTheUsaForms_AsTheProgram()
    {
        // seqret 6.6.0 USA forms (A5-8c3h, F255): ranges, strands, list files, queries on files, -sbegin / -send / -sreverse, -sid; the texts are the ones of the installed program
        const string fa = ">a desc\nACGTACGTAAACCCGGGTTT\n>b\nGGGGCCCCAAAATTTT\n";
        Assert.That(AnalysisTools.Seqret(fa, sequenceValue: "sequence.dat[3:10]").SequenceText, Is.EqualTo(">a desc\nGTACGTAA\n>b\nGGCCCCAA\n"));
        Assert.That(AnalysisTools.Seqret(fa, sequenceValue: "sequence.dat[3:10:r]").SequenceText, Is.EqualTo(">a desc\nTTACGTAC\n>b\nTTGGGGCC\n"));
        Assert.That(AnalysisTools.Seqret(fa, qualifiers: "-sbegin1 3 -send1 10 -sreverse1").SequenceText, Is.EqualTo(">a desc\nTTACGTAC\n>b\nTTGGGGCC\n"));
        var files = new Dictionary<string, string> { ["l.lst"] = "sequence.dat[2:5]\nfasta::sequence.dat[1:3:r] # comment\n" };
        Assert.That(AnalysisTools.Seqret(fa, sequenceValue: "@l.lst", files: files).SequenceText, Is.EqualTo(">a desc\nCGTA\n>b\nGGGGCCCCAAAATTTT\n>a desc\nCGT\n>b\nCCC\n"));
        Assert.That(AnalysisTools.Seqret(fa, sequenceValue: "sequence.dat:b").SequenceText, Is.EqualTo(">b\nGGGGCCCCAAAATTTT\n"));
        Assert.That(AnalysisTools.Seqret(fa, qualifiers: "-sid1 zzz").SequenceText, Is.EqualTo(">zzz desc\nACGTACGTAAACCCGGGTTT\n>zzz\nGGGGCCCCAAAATTTT\n"));
        Assert.That(() => AnalysisTools.Seqret(fa, sequenceValue: "missing.dat"), Throws.ArgumentException);
    }

    [Test]
    public void Entret_ReturnsTheTextOfTheEntries_AsTheProgram()
    {
        // entret 6.6.0 (A5-8c3h, F257): the text of every entry as the reader stored it, a query on the file, -firstonly
        const string fa = ">s1 desc one\nACGTACGTAC\n>s2\nTTGGCCAATT\n";
        Assert.That(AnalysisTools.Entret(fa).Text, Is.EqualTo(fa));
        Assert.That(AnalysisTools.Entret(fa, sequenceValue: "sequence.dat:s2").Text, Is.EqualTo(">s2\nTTGGCCAATT\n"));
        Assert.That(AnalysisTools.Entret(fa, firstOnly: true).Entries, Is.EqualTo(1));
        const string embl = "ID   X1; SV 1; linear; mRNA; STD; HUM; 10 BP.\nXX\nSQ   Sequence 10 BP;\n     acgtacgtac                                                    10\n//\n";
        Assert.That(AnalysisTools.Entret(embl).Text, Is.EqualTo(embl));
        const string defs = "DBNAME tdb [ type: N format: fasta method: direct directory: \"d\" filename: \"*.fa\" ]\n";
        var files = new Dictionary<string, string> { ["d/a.fa"] = fa };
        Assert.That(AnalysisTools.Entret("x", sequenceValue: "tdb:s1", files: files, definitions: defs).Text, Is.EqualTo(">s1 desc one\nACGTACGTAC\n"));
        Assert.That(() => AnalysisTools.Entret(" "), Throws.ArgumentException);
    }

    [Test]
    public void Union_JoinsTheSequences_AsTheProgram()
    {
        // union 6.6.0 (A5-8c3h, F253)
        const string fa = ">f1\nACGTACGTAC\n>f2\nCGTACGTACGGG\n>f3\nGGGTTTAAA\n";
        var r = AnalysisTools.Union(fa, findOverlap: true, overlapFile: true);
        Assert.That(r.SequenceText, Is.EqualTo(">f1\nACGTACGTACGGGTTTAAA\n"));
        Assert.That(r.OverlapText, Is.EqualTo("0\n9\n3\n"));
        Assert.That(AnalysisTools.Union(fa).SequenceText, Is.EqualTo(">f1\nACGTACGTACCGTACGTACGGGGGGTTTAAA\n"));
        Assert.That(() => AnalysisTools.Union(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqretsplit_WritesOneFilePerSequence_AsTheProgram()
    {
        // seqretsplit 6.6.0 (A5-8c3h, F254)
        var r = AnalysisTools.Seqretsplit(">a\nACGT\n>A\nGG\n>b\nTT\n");
        Assert.That(r.SingleFiles.Keys, Is.EquivalentTo(new[] { "a.fasta", "a.002.fasta", "b.fasta" }));
        Assert.That(r.SingleFiles["a.002.fasta"], Is.EqualTo(">A.002\nGG\n"));
        Assert.That(r.Messages, Does.Contain("Warning: Duplicate name 'A' changed to 'A.002'"));
        Assert.That(AnalysisTools.Seqretsplit(">a\nACGT\n>b\nTT\n", outputSingle: false).SequenceText, Is.EqualTo(">a\nACGT\n>b\nTT\n"));
        Assert.That(() => AnalysisTools.Seqretsplit(" "), Throws.ArgumentException);
    }

    [Test]
    public void Diffseq_FeaturesFromSeparateFiles_AsTheProgram()
    {
        // diffseq 6.6.0 a.fa b.fa -wordsize 5 -ufo1 f1.gff -ufo2 swiss::f2.dat (A5-8c3e, F241): the features of the sequences come from the files, whatever their format
        const string seq = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        string snp = seq[..30] + "C" + seq[31..];
        const string f1 = "##gff-version 3\nseqa\tsrc\tgene\t20\t40\t.\t+\t.\tID=g1;Name=abc\n";
        const string f2 = "ID   SEQB_HUMAN   Reviewed;   80 AA.\nFT   DOMAIN         25     35       a domain.\n//\n";
        var files = new Dictionary<string, string> { ["f1.gff"] = f1, ["f2.dat"] = f2 };
        var r = AnalysisTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-ufo1 f1.gff -ufo2 swiss::f2.dat", default, files);
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(31, 31, 1, 31, 31, 1) }));
        Assert.That(r.Report, Does.Contain("Feature: gene 20-40\n"));
        Assert.That(r.Report, Does.Contain("Feature: polypeptide_domain 25-35 note='a domain'"));
        // a UFO that names no file gives no features
        var none = AnalysisTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-ufo1 missing.gff", default, files);
        Assert.That(none.Report, Does.Not.Contain("Feature: "));
    }

    [Test]
    public void Diffseq_FeatureOutputFormats_AsTheProgram()
    {
        // diffseq 6.6.0 a.fa b.fa -wordsize 5 -offormat4 embl -offormat5 bed (A5-8c3c, F220): each feature file in its own format, and the files the tool names
        const string seq = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        string snp = seq[..30] + "C" + seq[31..];
        var r = AnalysisTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-offormat4 embl -offormat5 bed");
        Assert.That(r.FeatureFile1, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   conflict        31..31\nFT                   /note=\"SNP in seqb\"\nFT                   /replace=\"C\"\n"));
        Assert.That(r.FeatureFile2, Is.EqualTo("track name=seqb description=\"EMBOSS:6.6.0.0\" useScore=0\nseqb\t30\t30\tsequence_conflict\t0\t+\t30\t30\t0\t1\t0\t1\n"));
        Assert.That(r.FeatureFileName1, Is.EqualTo("seqa.diffgff"));
        // a protein pair in the Swiss-Prot and PIR formats; a format the program refuses for the kind of sequence writes an empty file
        var p = AnalysisTools.Diffseq("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKQRQISFVKTHFSRQ", 3, false, "pa", "pb", "-offormat4 swiss -offormat5 pir");
        Assert.That(p.FeatureFile1, Is.EqualTo("FT   CONFLICT     17     17       S -> T.\n"));
        Assert.That(p.FeatureFile2, Is.EqualTo("F;17/Region: T -> S, sequence_conflict\n"));
        var refused = AnalysisTools.Diffseq("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKQRQISFVKTHFSRQ", 3, false, "pa", "pb", "-offormat4 embl");
        Assert.That(refused.FeatureFile1, Is.Empty);
        Assert.That(refused.FeatureFileName1, Is.EqualTo("pa.diffgff"));
    }

    [Test]
    public void Wordmatch_FeatureOutputFormats_AsTheProgram()
    {
        // wordmatch 6.6.0 -offormat4 embl -offormat5 genbank on the same pair (A5-8c3d, F225): the feature files in the EMBL feature table and the GenBank FEATURES layout
        var r = AnalysisTools.Wordmatch(">p\nACGTACGTAC\n", ">g\nTTACGTACGTACTT\n", 6, featureQualifiers: "-offormat4 embl -offormat5 genbank", featureFile1: "pa.out", featureFile2: "gb.out");
        Assert.Multiple(() =>
        {
            Assert.That(r.FeatureFile1, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   misc_feature    4..10\nFT                   /note=\"g\"\nFT   misc_feature    2..10\nFT                   /note=\"g\"\n"));
            Assert.That(r.FeatureFile2, Is.EqualTo("FEATURES             Location/Qualifiers\n     misc_feature    2..8\n                     /note=\"p\"\n     misc_feature    4..12\n                     /note=\"p\"\n"));
            Assert.That((r.FeatureFileName1, r.FeatureFileName2), Is.EqualTo(("pa.out", "gb.out")));
        });
    }

    [Test]
    public void Wordmatch_SequenceQualifiers_AsTheProgram()
    {
        // wordmatch 6.6.0 -sbegin1 3 -send1 18 -sreverse2 -sbegin2 4 -send2 40 (A5-8c3e, F239): the set is cut to residues 3..18, the scanned sequences are reversed (their -sbegin2 / -send2 are ignored by the program)
        var r = AnalysisTools.Wordmatch(">p1 desc one\nACGTACGTTAGCATCGATCG\n>p2\nTTAGCATCGGGATTACGAT\n", ">g1 genome\nGGACGTACGTTAGCATCGATCGAAATTAGCATCGGGATTACGATTT\n>g2\nCCCTTAGCATCGATCGATCGATGGG\n", 5,
            featureQualifiers: "-sbegin1 3 -send1 18 -sreverse2 -sbegin2 4 -send2 40");
        Assert.Multiple(() =>
        {
            Assert.That(r.AllMatches, Is.EqualTo(7));
            Assert.That(r.Alignments, Does.Contain("     5 p1              +       12..16       g1              -       28..24\n     6 p1              +       11..16       g1              -       32..27\n     6 p1              +        1..6        g1              -       44..39\n"));
            Assert.That(r.FeatureFile1, Does.Contain("##sequence-region p1 1 16\n"));
        });
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
