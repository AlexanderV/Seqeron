// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{

    [Test]
    public void ShowAlign_RendersTheAlignmentAsShowalign()
    {
        // showalign 6.6.0 on s1 acgtac-gtaGT / s2 acgtaccgtaGA / s3 acctaccgtTGA (default options) and with -show I -refseq 2 (B09 A3-13c, F148)
        var r = EmbossPairwiseAlignmentTools.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"]);
        Assert.That(r.Text, Is.EqualTo("                  10  \n          ----:----|--\ns1        ......-....T\ns2        ............\ns3        ..C......T..\nConsensus acgtaccgtaGA\n\n"));
        var identities = EmbossPairwiseAlignmentTools.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"], refseq: "2", show: "I");
        Assert.That(identities.Text, Does.Contain("s1        acgtac-gtaG.\ns3        ac.taccgt.GA\nConsensus acgtaccgtaGA\n"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["AC", "AC"], show: "Z"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["AC", "AC"], uppercase: "5"));
    }

    [Test]
    public void ShowAlign_Html_WrapsAndColoursAsShowalign()
    {
        // showalign 6.6.0 -html -highlight "3-8 blue 10-12 orange" -width 30 (A3-13c2): the ranges are coloured on the sequence lines and the consensus, not on the names or ticks
        var r = EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["ACGTACGTAAACCCGGG", "ACGTACGAAAACCCGGT"], html: true, highlight: "3-8 blue 10-12 orange", width: 30);
        Assert.That(r.Text, Does.StartWith("<pre>\n").And.EndWith("</pre>\n"));
        Assert.That(r.Text, Does.Contain("a         ..<font color=blue>......</font>.<font color=orange>...</font>").And.Contain("Consensus AC<font color=blue>GTACGt</font>A<font color=orange>AAC</font>CCGGg"));
        Assert.That(EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["AC", "AC"], highlight: "1-2 red").Text, Does.Not.Contain("<font"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["AC", "AC"], html: true, highlight: "0-9 red"));
    }

    [Test]
    public void Diffseq_SnpsAndTheSevenTagReport_AsTheProgram()
    {
        // diffseq 6.6.0 on ACGTACGTACGTAAGGCCTTAACCGGTT / ACGTACGTACGTAAGGCATTAACCGGTT -wordsize 8 (A5-8c3, F196): one SNP, C -> A, a transversion
        var r = EmbossPairwiseAlignmentTools.Diffseq("ACGTACGTACGTAAGGCCTTAACCGGTT", "ACGTACGTACGTAAGGCATTAACCGGTT", 8, false, "a", "b");
        Assert.That(r.Differences, Has.Length.EqualTo(1));
        Assert.That(r.Differences[0], Is.EqualTo(new DiffseqDifferenceItem(18, 18, 1, 18, 18, 1)));
        Assert.That((r.Snps, r.Transitions, r.Transversions, r.OverlapStart1, r.OverlapEnd1), Is.EqualTo((1, 0, 1, 1, 28)));
        Assert.That(r.Report, Does.Contain("\n\na 18-18 Length: 1\nSequence: C\nSequence: A\nb 18-18 Length: 1\n"));
        Assert.That(r.Report, Does.Contain("# SNP_count: 1\n# Transitions: 0\n# Transversions: 1\n"));
        Assert.That(r.FeatureFile1, Does.Contain("a\tdiffseq\tsequence_conflict\t18\t18\t1\t+\t.\tID=a.1;note=SNP in b;replace=A"));
        Assert.That(EmbossPairwiseAlignmentTools.Diffseq("ACGTACGTACGTAAGGCCTTAACCGGTT", "ACGTACGTACGTAAGGCATTAACCGGTT", 8, false, "a", "b", "-rformat excel").Report, Does.StartWith("SeqName\tStart\tEnd\tScore\tStrand\tstart\tend\tlength\tname\tsequence\tfirst_feature\tsecond_feature\n"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Diffseq("ACGT", "ACGT", 8, false, "a", "b", "pir"));
    }

    [Test]
    public void Diffseq_SequenceQualifiers_AsTheProgram()
    {
        // diffseq 6.6.0 qa.fa qb.fa -wordsize 5 -sbegin1 11 -send1 70 -sbegin2 11 -send2 70 -rformat table (A5-8c3b, F216): the SNP of position 31 is at 21 of the cut sequences, the header keeps the original range of the first and 1 - 60 for the second
        const string a = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        const string b = "TTTCCTCATGCAATTCAAAACCATGTCCGTCATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        var r = EmbossPairwiseAlignmentTools.Diffseq(a, b, 5, false, "seqa", "seqb", "-sbegin1 11 -send1 70 -sbegin2 11 -send2 70 -rformat table");
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(21, 21, 1, 21, 21, 1) }));
        Assert.That(r.Report, Does.Contain("# Sequence: seqa     from: 11   to: 70\n"));
        Assert.That(r.Report, Does.Contain("# Compare: seqb     from: 1   to: 60\n"));
        Assert.That(r.Report, Does.Contain("     21      21       +   0.000     21     21      1 seqb   C        .             .             \n"));
        Assert.That(r.FeatureFile1, Does.Contain("seqa\tdiffseq\tsequence_conflict\t21\t21\t1\t+\t.\tID=seqa.1;note=SNP in seqb;replace=C"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Diffseq(a, b, 5, false, "seqa", "seqb", "-sbegin3 4"));
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
        var r = EmbossPairwiseAlignmentTools.Diffseq(genbank, embl, 5, false);
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(31, 31, 1, 31, 31, 1) }));
        Assert.That(r.Report, Does.Contain("Feature: gene 20-60 gene='abc'\nFeature: CDS 25-35 note='a note'\n"));
        Assert.That(r.Report, Does.Contain("Feature: misc_feature 28-33 gene='x1'\n"));
        Assert.That(r.Report, Does.Contain("SEQA 31-31 Length: 1\n"));
        // the special qualifiers (A5-8c3c, F219): a valid value is kept, a bad one is a note, a correctable one is corrected (the program also warns on stderr)
        string special = genbank.Replace("/translation=\"MKV\"", "/collection_date=\"2020\"\n                     /lat_lon=\"12.5 N 33.1\"\n                     /replace=\"a c g\"", StringComparison.Ordinal);
        var s2 = EmbossPairwiseAlignmentTools.Diffseq(special, embl, 5, false);
        Assert.That(s2.Report, Does.Contain("Feature: CDS 25-35 note='a note' collection_date='2020' note='*lat_lon: 12.5 N 33.1' replace='acg'\n"));
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
        var r = EmbossPairwiseAlignmentTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-ufo1 f1.gff -ufo2 swiss::f2.dat", default, files);
        Assert.That(r.Differences, Is.EqualTo(new[] { new DiffseqDifferenceItem(31, 31, 1, 31, 31, 1) }));
        Assert.That(r.Report, Does.Contain("Feature: gene 20-40\n"));
        Assert.That(r.Report, Does.Contain("Feature: polypeptide_domain 25-35 note='a domain'"));
        // a UFO that names no file gives no features
        var none = EmbossPairwiseAlignmentTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-ufo1 missing.gff", default, files);
        Assert.That(none.Report, Does.Not.Contain("Feature: "));
    }

    [Test]
    public void Diffseq_FeatureOutputFormats_AsTheProgram()
    {
        // diffseq 6.6.0 a.fa b.fa -wordsize 5 -offormat4 embl -offormat5 bed (A5-8c3c, F220): each feature file in its own format, and the files the tool names
        const string seq = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCT";
        string snp = seq[..30] + "C" + seq[31..];
        var r = EmbossPairwiseAlignmentTools.Diffseq(seq, snp, 5, false, "seqa", "seqb", "-offormat4 embl -offormat5 bed");
        Assert.That(r.FeatureFile1, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   conflict        31..31\nFT                   /note=\"SNP in seqb\"\nFT                   /replace=\"C\"\n"));
        Assert.That(r.FeatureFile2, Is.EqualTo("track name=seqb description=\"EMBOSS:6.6.0.0\" useScore=0\nseqb\t30\t30\tsequence_conflict\t0\t+\t30\t30\t0\t1\t0\t1\n"));
        Assert.That(r.FeatureFileName1, Is.EqualTo("seqa.diffgff"));
        // a protein pair in the Swiss-Prot and PIR formats; a format the program refuses for the kind of sequence writes an empty file
        var p = EmbossPairwiseAlignmentTools.Diffseq("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKQRQISFVKTHFSRQ", 3, false, "pa", "pb", "-offormat4 swiss -offormat5 pir");
        Assert.That(p.FeatureFile1, Is.EqualTo("FT   CONFLICT     17     17       S -> T.\n"));
        Assert.That(p.FeatureFile2, Is.EqualTo("F;17/Region: T -> S, sequence_conflict\n"));
        var refused = EmbossPairwiseAlignmentTools.Diffseq("MKTAYIAKQRQISFVKSHFSRQ", "MKTAYIAKQRQISFVKTHFSRQ", 3, false, "pa", "pb", "-offormat4 embl");
        Assert.That(refused.FeatureFile1, Is.Empty);
        Assert.That(refused.FeatureFileName1, Is.EqualTo("pa.diffgff"));
    }

    [Test]
    public void Megamerger_MergesByWordMatchAnchors_AsTheProgram()
    {
        // megamerger a=ACGTACGTAC b=ACGTACGTACGGGGGG -wordsize 4 (B09 A3-13b, F149)
        var r = EmbossPairwiseAlignmentTools.Megamerger("ACGTACGTAC", "ACGTACGTACGGGGGG", 4, name1: "a", name2: "b");
        Assert.That(r.MergedSequence, Is.EqualTo("acgtacgtacgggggg"));
        Assert.That(r.Report, Does.Contain("Matching region a 1-10 : b 1-10\nLength of match: 10\n").And.Contain("Using b 11-16 as the final sequence"));
        Assert.That(r.Matches.Single(), Is.EqualTo(new WordMatchItem(0, 1, 0, 0, 10)));
        Assert.That(EmbossPairwiseAlignmentTools.Megamerger("ACGTACGTACGGGGGG", "ACGTACGTAC", 12).MergedSequence, Is.Empty);
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Megamerger("ACGTACGTAC", "ACGTACGTACGGGGGG", 20));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Megamerger("", "ACGT"));
    }

    [Test]
    public void Est2genome_AlignsASplicedEstLikeTheProgram()
    {
        // est2genome -estsequence e.fa -genomesequence g.fa (EMBOSS 6.6.0): three exons, two GT/AG introns (B09 A2-28d, F203)
        string est = ">est1 an est\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCT\n";
        string genome = ">genomic some description\nACGTGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTTTTT\n";
        var r = EmbossPairwiseAlignmentTools.Est2genome(est, genome);
        Assert.That(r.Text, Is.EqualTo("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\nExon        60 100.0     5    64 genomic          1    60 est1          an est\n+Intron    -20   0.0    65   138 genomic     \nExon        80 100.0   139   218 genomic         61   140 est1          an est\n+Intron    -20   0.0   219   262 genomic     \nExon        50 100.0   263   312 genomic        141   190 est1          an est\n\nSpan       150 100.0     5   312 genomic          1   190 est1          an est\n\nSegment     60 100.0     5    64 genomic          1    60 est1          an est\nSegment     80 100.0   139   218 genomic         61   140 est1          an est\nSegment     50 100.0   263   312 genomic        141   190 est1          an est\n"));
        Assert.That(r.Messages, Is.Empty);
        Assert.That(EmbossPairwiseAlignmentTools.Est2genome(est, genome, minScore: 200).Text, Is.EqualTo("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\n"));
        Assert.Throws<InvalidOperationException>(() => EmbossPairwiseAlignmentTools.Est2genome(est, genome, mode: "forward", best: false));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Est2genome(est, genome, mode: "x"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Est2genome("", genome));
    }

    /// <summary>A4-1 (B09 F157): fuzztran, showalign and consambig took no token; they now hand the request token to the ports.</summary>
    [Test]
    public void A4_1_FuzztranShowalignConsambig_StopWithACancelledToken()
    {
        var cancelled = new CancellationToken(true);
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => EmbossMotifsTools.Fuzztran("ATGAAAGTTCTGGCT", "KV", cancellationToken: cancelled), "fuzztran");
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.ShowAlign(["a", "b"], ["ACGTACGT", "ACGTACGA"], cancellationToken: cancelled), "showalign");
            Assert.Catch<OperationCanceledException>(() => EmbossMotifsTools.Consambig(["ACGTACGT", "ACGTACGA"], cancellationToken: cancelled), "consambig");
        });
        Assert.That(EmbossMotifsTools.Consambig(["ACGTACGT", "ACGTACGA"], cancellationToken: CancellationToken.None).Consensus, Is.EqualTo("ACGTACGW"));
    }}
