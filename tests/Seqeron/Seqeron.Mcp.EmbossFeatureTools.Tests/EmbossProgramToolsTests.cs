// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossFeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossFeatureTools.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{

    [Test]
    public void Featcopy_ReadsAnyFormatAndWritesTheOthers_AsTheProgram()
    {
        // featcopy 6.6.0 (A5-8c3e, F241), dates as of the run: a GFF3 file with Name / Parent / product / percent-encoded values written as EMBL, GenBank, GFF2 (the sub-feature is written as its own feature)
        const string gff = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc;product=foo\nseq1\tsrc\tmRNA\t20\t40\t5.5\t-\t.\tID=m1;Parent=g1;note=hello%20world\n";
        var embl = EmbossFeatureToolsTools.Featcopy(gff, null, "embl");
        Assert.That(embl.Text, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   gene            10..50\nFT                   /note=\"*Name abc\"\nFT                   /product=\"foo\"\nFT   mRNA            complement(20..40)\nFT                   /note=\"hello world\"\n"));
        Assert.That(embl.Tables, Is.EqualTo(1));
        Assert.That(EmbossFeatureToolsTools.Featcopy(gff, null, "gb").Text, Does.StartWith("FEATURES             Location/Qualifiers\n     gene            10..50\n                     /note=\"*Name abc\"\n"));
        var gff2 = EmbossFeatureToolsTools.Featcopy(gff, null, "gff2");
        Assert.That(gff2.Text, Does.EndWith("##Type unknown <> seq1\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tSequence \"seq1.2\" ; product \"foo\"\n"));
        // the format:: prefix forces a reader: SwissProt feature lines, PIR feature lines (a text with no feature of the format is refused as the program refuses it)
        const string swiss = "ID   KIN1_HUMAN   Reviewed;   120 AA.\nAC   P12345;\nFT   DOMAIN         5     50       Kinase domain.\nFT   ACT_SITE     30     30       Proton acceptor (By similarity).\n//\n";
        var sw = EmbossFeatureToolsTools.Featcopy(swiss, "swiss::features.dat", "gff3");
        Assert.That(sw.Text, Does.EndWith("#!Type Protein\n#!Source-version EMBOSS 6.6.0.0\nKIN1_HUMAN\tSWISSPROT\tpolypeptide_domain\t5\t50\t.\t.\t.\tID=KIN1_HUMAN.1;note=Kinase domain\nKIN1_HUMAN\tSWISSPROT\tactive_site_residue\t30\t30\t.\t.\t.\tID=KIN1_HUMAN.2;note=Proton acceptor;comment=By similarity\n"));
        Assert.That(EmbossFeatureToolsTools.Featcopy(swiss, "swiss::features.dat", "gff2").Text, Does.Contain("Sequence \"KIN1_HUMAN.2\" ; note \"Proton acceptor\" ; comment \"By similarity\"\n"));
        Assert.That(() => EmbossFeatureToolsTools.Featcopy("F;1-10/Domain: x\n", "gff2::features.dat", "gff3"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        // an unknown format name is an error message of the library, the formats are tried as for an unknown input
        var bogus = EmbossFeatureToolsTools.Featcopy(gff, "bogus::features.dat", "gff3");
        Assert.That(bogus.Messages, Does.Contain("Error: Unknown input format 'bogus'"));
        Assert.That(bogus.Tables, Is.EqualTo(1));
        // the debug trace and a range (-fbegin / -fend cut the table: positions shift, a feature cut at the end is flagged)
        var cut = EmbossFeatureToolsTools.Featcopy(gff, null, "debug", 10, 30);
        Assert.That(cut.Text, Does.StartWith("  Format: 0 (gff3)  Start: 10\n  End: 30\n  Len: 21\n  Offset: 0\n  Groups: 4\n"));
        Assert.That(cut.Text, Does.Contain("    Start: 1\n    End: 21\n").And.Contain("      END_AFTER_SEQ\n"));
        // a format the program refuses for the kind of sequence: an error message of the library and nothing written (embl for a protein table)
        var refused = EmbossFeatureToolsTools.Featcopy(swiss, "swiss::features.dat", "embl");
        Assert.That(refused.Messages, Does.Contain("Error: Feature format 'embl' not supported for protein sequence features"));
        Assert.That(refused.Tables, Is.EqualTo(0));
        // the other writers of the feature library (A5-8c3f, F244): bed with its blocks, Swiss-Prot FT lines, PIR, DAS GFF, draw
        const string gb = "LOCUS       SEQA 300 bp    DNA     linear   UNK 01-JAN-2000\nFEATURES             Location/Qualifiers\n     gene            complement(10..50)\n                     /gene=\"abc\"\n     CDS             join(60..80,90..120)\n                     /product=\"foo\"\nORIGIN\n        1 acgt\n//\n";
        Assert.That(EmbossFeatureToolsTools.Featcopy(gb, null, "bed").Text, Is.EqualTo("track name= description=\"EMBOSS:6.6.0.0\" useScore=0\n\t9\t49\tgene\t0\t-\t9\t49\t0\t1\t0\t41\n\t59\t119\tCDS\t0\t+\t59\t119\t0\t2\t21,31\t0,30\n"));
        Assert.That(EmbossFeatureToolsTools.Featcopy(gb, null, "draw").Text, Does.Contain("Range 10 50 | | 9 H\ngene\nendlabel\n").And.EndWith("endgroup\n"));
        Assert.That(EmbossFeatureToolsTools.Featcopy(swiss, "swiss::features.dat", "swiss").Text, Is.EqualTo("FT   DOMAIN        5     50       Kinase domain.\nFT   ACT_SITE     30     30       Proton acceptor (By similarity).\n"));
        Assert.That(EmbossFeatureToolsTools.Featcopy(swiss, "swiss::features.dat", "pir").Text, Is.EqualTo("F;5-50/Domain: Kinase domain\nF;30/Active site: Proton acceptor #By similarity\n"));
        Assert.That(EmbossFeatureToolsTools.Featcopy(gb, null, "dasgff").Text, Does.Contain("<TYPE id=\"SO:0000704\" category=\"other\" reference=\"no\" subparts=\"no\">gene</TYPE>"));
        Assert.That(() => EmbossFeatureToolsTools.Featcopy(" ", null, "gff3"), Throws.ArgumentException);
    }

    [Test]
    public void Featcopy_ListFilesAsFeatureInput_AsTheProgram()
    {
        // featcopy 6.6.0 -features @list (A5-8c3f, F247): the words of a list file are queries, a '#' starts a comment, a missing file is an error message and the list goes on
        const string one = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1\n";
        const string two = "##gff-version 3\n##sequence-region seq2 1 100\nseq2\tsrc\texon\t20\t30\t.\t+\t.\tID=e1\n";
        var files = new Dictionary<string, string> { ["a.gff"] = one, ["b.gff"] = two, ["l.txt"] = "# the files\na.gff\nmissing.gff\nb.gff # last\n" };
        var bed = EmbossFeatureToolsTools.Featcopy(one, "@l.txt", "bed", 1, null, default, files);
        Assert.That(bed.Text, Is.EqualTo("track name=seq1 description=\"EMBOSS:6.6.0.0\" useScore=0\nseq1\t9\t49\tgene\t0\t+\t9\t49\t0\t1\t0\t41\ntrack name=seq2 description=\"EMBOSS:6.6.0.0\" useScore=0\nseq2\t19\t29\texon\t0\t+\t19\t29\t0\t1\t0\t11\n"));
        var all = EmbossFeatureToolsTools.Featcopy(one, "@l.txt", "gff3", 1, null, default, files);
        Assert.That(all.Tables, Is.EqualTo(2));
        Assert.That(all.Messages, Does.Contain("Error: Failed to open filename 'missing.gff'").And.Contain("Error: Failed to read data data 'missing.gff'"));
        Assert.That(all.Text, Does.Contain("seq1\tsrc\tgene\t10\t50").And.Contain("seq2\tsrc\texon\t20\t30"));
        Assert.That(() => EmbossFeatureToolsTools.Featcopy(one, "@nolist.txt", "gff3"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        Assert.That(() => EmbossFeatureToolsTools.Featcopy(one, "asis::ACGT", "gff3", 1, null, default, files), Throws.ArgumentException.With.Message.Contains("asis:: not supported by datatype"));
    }

    [Test]
    public void Featmerge_CombinesTwoFeatureTables_AsTheProgram()
    {
        // featmerge 6.6.0 (A5-8c3f, F245): the features of the second table are copied into the first, which is written in the output format
        const string a = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\tgene\t10\t50\t.\t+\t.\tID=g1;Name=abc\n";
        const string b = "##gff-version 3\n##sequence-region seq1 1 100\nseq1\tsrc\texon\t20\t30\t.\t+\t.\tID=e1\n";
        var merged = EmbossFeatureToolsTools.Featmerge(a, b);
        Assert.That(merged.Tables, Is.EqualTo(1));
        Assert.That(merged.Text, Does.Contain("seq1\tsrc\tgene\t10\t50").And.Contain("seq1\tsrc\texon\t20\t30"));
        Assert.That(EmbossFeatureToolsTools.Featmerge(a, b, null, null, "bed").Text, Does.StartWith("track name=seq1 description=\"EMBOSS:6.6.0.0\" useScore=0\n"));
        Assert.That(() => EmbossFeatureToolsTools.Featmerge(a, b, "missing.dat"), Throws.ArgumentException.With.Message.Contains("Bad value for '-afeatures'"));
        Assert.That(() => EmbossFeatureToolsTools.Featmerge(" ", b), Throws.ArgumentException);
    }

    [Test]
    public void Feattext_ReturnsTheStoredEntryText_AsTheProgram()
    {
        // feattext 6.6.0 (A5-8c3f, F246): the EMBL reader keeps the lines of the entry, the PIR reader none
        const string embl = "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   gene            1..30\nXX\n//\n";
        var r = EmbossFeatureToolsTools.Feattext(embl);
        Assert.That(r.Tables, Is.EqualTo(1));
        Assert.That(r.Text, Is.EqualTo(embl));
        Assert.That(EmbossFeatureToolsTools.Feattext("F;1-10/Domain: x\n", "pir::features.dat").Text, Is.EqualTo(""));
        Assert.That(() => EmbossFeatureToolsTools.Feattext("nothing here\n"), Throws.ArgumentException.With.Message.Contains("Bad value for '-features'"));
        Assert.That(() => EmbossFeatureToolsTools.Feattext(" "), Throws.ArgumentException);
    }

    [Test]
    public void Maskfeat_MasksTheFeatureRegions_AsTheProgram()
    {
        // maskfeat 6.6.0 (A5-8c3g, F249): the default type repeat*, a type list, a mask character, lower case
        const string embl = "ID   SEQA; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   repeat_region   3..9\nFT   gene            12..18\nXX\nSQ   Sequence 20 BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n     aacgtacgta cgtacgtacg                                                20\n//\n";
        Assert.That(EmbossFeatureToolsTools.Maskfeat(embl, "embl").SequenceText, Is.EqualTo(">SEQA SEQA.1\naaNNNNNNNacgtacgtacg\n"));
        Assert.That(EmbossFeatureToolsTools.Maskfeat(embl, "embl", "repeat* gene", maskChar: "Q").SequenceText, Is.EqualTo(">SEQA SEQA.1\naaQQQQQQQacQQQQQQQcg\n"));
        Assert.That(EmbossFeatureToolsTools.Maskfeat(embl, "embl", "gene", toLower: true, qualifiers: "-supper1").SequenceText, Is.EqualTo(">SEQA SEQA.1\nAACGTACGTACgtacgtaCG\n"));
        Assert.That(() => EmbossFeatureToolsTools.Maskfeat(embl, "embl", maskChar: "AB"), Throws.ArgumentException);
        Assert.That(() => EmbossFeatureToolsTools.Maskfeat(" "), Throws.ArgumentException);
    }

    [Test]
    public void Featreport_NamesTheReportFile_AsTheProgram()
    {
        // featreport 6.6.0 (B09 WP60, F278): -rname3 -rextension3 -rdirectory3 give the report file (a directory is put in front of a name that has none); the default name is the name of the sequence + .featreport
        var r = EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "-rformat excel -rname rep -rextension xl -rdirectory d");
        Assert.That(r.ReportFile, Is.EqualTo("d/rep.xl"));
        Assert.That(r.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\nX1\t3\t29\t0.000\t+\nX1\t1\t40\t0.000\t+\nX1\t45\t55\t0.000\t+\n"));
        Assert.That(EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "-rformat excel -rdirectory d").ReportFile, Is.EqualTo("d/x1.featreport"));
        Assert.That(EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "excel").ReportFile, Is.EqualTo("stdout"));
        // the qualifiers of the features input have no effect on the report (measured)
        Assert.That(EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "excel", qualifiers: "-fformat2 genbank -freverse2 -fcircular2 -fask2").Report, Is.EqualTo(r.Report));
    }

    [Test]
    public void Extractfeat_ExtractsTheFeatures_AsTheProgram()
    {
        // extractfeat 6.6.0 (A5-8c3h, F258): the features that pass the filters as sequences, with the type in the name and the tags in the description
        var r = EmbossFeatureToolsTools.Extractfeat(EmblWithFeatures, type: "mRNA exon", describe: "gene", featInName: true);
        Assert.That(r.SequenceText, Is.EqualTo(">X1_1_40_mRNA [mRNA] test\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n>X1_45_55_exon [exon] test\nacgtacgtacg\n"));
        Assert.That(r.Sequences, Is.EqualTo(2));
        Assert.That(EmbossFeatureToolsTools.Extractfeat(EmblWithFeatures, type: "exon", before: 2, after: 3).SequenceText, Is.EqualTo(">X1_45_55 [exon] test\ngtacgtacgtacgtac\n"));
        Assert.That(EmbossFeatureToolsTools.Extractfeat(EmblWithFeatures, type: "gene").Messages, Is.EqualTo(new[] { "Warning: No sequences written to output file 'stdout'" }));
        // -ossingle2 with -osdirectory2 and -osextension2 (F273): one file per extracted sequence, name.extension in lower case
        var single = EmbossFeatureToolsTools.Extractfeat(EmblWithFeatures, type: "mRNA exon", featInName: true, outputSingle: true, outputDirectory: "d", outputExtension: "X1");
        Assert.That(single.SingleFiles!.Keys, Is.EquivalentTo(new[] { "d/x1_1_40_mrna.x1", "d/x1_45_55_exon.x1" }));
        Assert.That(single.SingleFiles["d/x1_45_55_exon.x1"], Is.EqualTo(">X1_45_55_exon [exon] test\nacgtacgtacg\n"));
        Assert.That(() => EmbossFeatureToolsTools.Extractfeat(" "), Throws.ArgumentException);
    }

    [Test]
    public void Coderet_WritesTheCodingSequences_AsTheProgram()
    {
        // coderet 6.6.0 (A5-8c3h, F259): the CDS, mRNA, translation and the non-coding regions of a sequence with features, the log of the counts
        var r = EmbossFeatureToolsTools.Coderet(EmblWithFeatures);
        Assert.That(r.Log, Is.EqualTo("   CDS  mRNA non-c Trans Total Sequence\n ===== ===== ===== ===== ===== ========\n     1     1     0     1     2 X1\n"));
        Assert.That(r.Files["x1.cds"], Is.EqualTo(">x1_cds_1\ngtacgtacgttacgtacgta\n"));
        Assert.That(r.Files["x1.prot"], Is.EqualTo(">x1_pro_1\nMKV\n"));
        Assert.That(r.Files["x1.noncoding"], Is.EqualTo(">x1_noncoding_41\nacgt\n>x1_noncoding_56\ntacgt\n"));
        Assert.That(EmbossFeatureToolsTools.Coderet(EmblWithFeatures, cds: false, mrna: false, translation: false).Files.Keys, Is.EqualTo(new[] { "x1.noncoding" }));
        // -osname, -osextension, -ossingle and -osdirectory of the outputs (F273)
        Assert.That(EmbossFeatureToolsTools.Coderet(EmblWithFeatures, mrna: false, translation: false, nonCoding: false, outputName: "a b", outputExtension: "ext").Files.Keys, Is.EqualTo(new[] { "a_b.ext" }));
        Assert.That(EmbossFeatureToolsTools.Coderet(EmblWithFeatures, mrna: false, translation: false, nonCoding: false, outputSingle: true, outputDirectory: "d").Files.Keys, Is.EqualTo(new[] { "d/x1_cds_1.cds" }));
        Assert.That(() => EmbossFeatureToolsTools.Coderet(" "), Throws.ArgumentException);
    }

    [Test]
    public void Showfeat_DrawsTheFeatures_AsTheProgram()
    {
        // showfeat 6.6.0 (A5-8c3h, F260): the features as lines of ASCII graphics, the options of the line, the parts of a join on one line
        var r = EmbossFeatureToolsTools.Showfeat(EmblWithFeatures, width: 30, joinFeatures: true, showStrand: true, showTags: true);
        Assert.That(r.Text, Is.EqualTo("X1\ntest\n|============================| 60\n|------------------>           + mRNA \n|---->   |--->                 + CDS  gene=\"g1\"\n                     |---->    + exon \n\n"));
        Assert.That(r.Sequences, Is.EqualTo(1));
        Assert.That(EmbossFeatureToolsTools.Showfeat(EmblWithFeatures, id: false, description: false, scale: false, width: 3, showPosition: true).Text, Is.EqualTo("1-40 mRNA\n3-12 CDS\n20-29 CDS\n45-55 exon\n\n"));
        Assert.That(EmbossFeatureToolsTools.Showfeat(EmblWithFeatures, width: -5).Messages, Is.EqualTo(new[] { "Warning: integer value out of range -5 less than (reset to) 0" }));
        Assert.That(() => EmbossFeatureToolsTools.Showfeat(" "), Throws.ArgumentException);
    }

    [Test]
    public void Featreport_ReportsTheFeatureTable_AsTheProgram()
    {
        // featreport 6.6.0 (A5-8c3h, F262): the feature table of a feature file reported for a sequence (excel: no header, no date)
        var r = EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "excel");
        Assert.That(r.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\nX1\t3\t29\t0.000\t+\nX1\t1\t40\t0.000\t+\nX1\t45\t55\t0.000\t+\n"));
        Assert.That(EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures).Report, Does.StartWith("##gff-version 3\n##sequence-region X1 1 55\n#!Date "));
        Assert.That(EmbossFeatureToolsTools.Featreport(EmblWithFeatures, EmblWithFeatures, reportFormat: "gff").Report, Does.Contain("X1\tEMBL\tbiological_region\t3\t29\t.\t+\t0\tID=X1.1;featflags=type:CDS;translation=MKV;gene=g1\n"));
        Assert.That(() => EmbossFeatureToolsTools.Featreport(EmblWithFeatures, " "), Throws.ArgumentException);
    }

    [Test]
    public void Twofeat_ReportsThePairs_AsTheProgram()
    {
        // twofeat 6.6.0 (A5-8c3h, F263): a mRNA feature A and an exon feature B, reported as one feature with the tags startA ... (excel: no header) or as the two features
        var r = EmbossFeatureToolsTools.Twofeat(EmblWithFeatures, aType: "mRNA", bType: "exon", reportFormat: "excel");
        Assert.That(r.Report, Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tstarta\tenda\tstartb\tendb\nX1\t1\t55\t0.000\t+\t1\t40\t45\t55\n"));
        Assert.That(EmbossFeatureToolsTools.Twofeat(EmblWithFeatures, aType: "mRNA", bType: "exon", twoOut: true, reportFormat: "excel").Report,
            Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tstarta\tenda\tstartb\tendb\nX1\t1\t40\t0.000\t+\t.\t.\t.\t.\nX1\t45\t55\t0.000\t+\t.\t.\t.\t.\n"));
        Assert.That(EmbossFeatureToolsTools.Twofeat(EmblWithFeatures, aType: "mRNA", bType: "exon", overlap: "O", reportFormat: "excel").Messages, Is.EqualTo(new[] { "Warning: No hits in output report file 'stdout'" }));
        Assert.That(() => EmbossFeatureToolsTools.Twofeat(" "), Throws.ArgumentException);
    }

    private const string EmblWithFeatures = "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             join(3..12,20..29)\nFT                   /translation=\"MKV\"\nFT                   /gene=\"g1\"\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n";
}
