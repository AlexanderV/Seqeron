// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{

    [Test]
    public void Seqret_ReadsTheSequenceFormats_AsTheProgram()
    {
        // seqret 6.6.0 (A5-8c3g, F248): the DE description of a Swiss-Prot entry, its features, the type letters of a FASTA file, an asis:: sequence
        const string sw = "ID   ABC_HUMAN   Reviewed;   10 AA.\nAC   P12345;\nDE   RecName: Full=Kinase alpha; Short=KA;\nDE   AltName: Full=Other name;\nFT   DOMAIN       1      5       Kinase domain.\nSQ   SEQUENCE   10 AA;  1000 MW;  0 CRC64;\n     MKVLAAGIVG\n//\n";
        var r = EmbossSequenceFormatsTools.Seqret(sw, "swiss", feature: true);
        Assert.That(r.SequenceText, Is.EqualTo(">ABC_HUMAN P12345 Kinase alpha (KA) (Other name)\nMKVLAAGIVG\n"));
        Assert.That(r.Sequences, Is.EqualTo(1));
        Assert.That(r.FeatureText, Does.Contain("ABC_HUMAN\tSWISSPROT\tpolypeptide_domain\t1\t5\t.\t.\t.\tID=ABC_HUMAN.1;note=Kinase domain"));
        var fa = EmbossSequenceFormatsTools.Seqret(">s1 first\nACGTXX\n>s2\nMKVL*\n");
        Assert.That(fa.SequenceText, Is.EqualTo(">s1 first\nACGTNN\n>s2\nMKVL*\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret("asis::ACGT.ACGT").SequenceText, Is.EqualTo(">asis\nACGT-ACGT\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(">a\nACGT\n>b\nTT\n", firstOnly: true).Sequences, Is.EqualTo(1));
        Assert.That(EmbossSequenceFormatsTools.Seqret("CLUSTAL W (1.83) multiple sequence alignment\n\n\nseq1      ACGT\nseq2      ACGA\n\n").SequenceText, Is.EqualTo(">seq1\nACGT\n>seq2\nACGA\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Seqret("AS 1 1\n"), Throws.ArgumentException);          // the ACE reader of the library dies by a segmentation fault
        Assert.That(() => EmbossSequenceFormatsTools.Seqret(" "), Throws.ArgumentException);
    }

    [Test]
    public void Descseq_ReplacesOrAppendsNameAndDescription_AsTheProgram()
    {
        // descseq 6.6.0 (A5-8c3g, F250)
        const string fa = ">seq1 old words\nACGTACGT\n>seq2\nTTTT\n";
        Assert.That(EmbossSequenceFormatsTools.Descseq(fa, name: "New", description: "new text").SequenceText, Is.EqualTo(">New new text\nACGTACGT\n"));
        Assert.That(EmbossSequenceFormatsTools.Descseq(fa, name: "_x", description: " more", append: true).SequenceText, Is.EqualTo(">seq1_x old words more\nACGTACGT\n"));
        Assert.That(EmbossSequenceFormatsTools.Descseq(fa).SequenceText, Is.EqualTo(">seq1 old words\nACGTACGT\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Descseq(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqret_WritesTheOutputFormats_AsTheProgram()
    {
        // seqret 6.6.0 sequence writers (A5-8c3h, F252): the entries are the ones of the installed program
        const string fa = ">s1 desc one\nACGTACGTAC\n>s2\nTTGGCCAATT\n";
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "pir").SequenceText, Is.EqualTo(">D1;s1\ndesc one, 10 bases\n ACGTACGTAC*\n\n>D1;s2\n, 10 bases\n TTGGCCAATT*\n\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(">p1 prot\nMKVLAAGIVG\n", outputFormat: "gcg").SequenceText, Is.EqualTo("!!AA_SEQUENCE 1.0\n\nprot\n\np1  Length: 10  Type: P  Check: 4069 ..\n\n  1 MKVLAAGIVG\n\n"));
        var single = EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "gcg", outputSingle: true, outputDirectory: "d", outputExtension: "X1");
        Assert.That(single.SingleFiles!.Keys, Is.EquivalentTo(new[] { "d/s1.x1", "d/s2.x1" }));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "experiment").SequenceText, Is.EqualTo("ID   s1         standard; DNA; UNC; 10 BP.\nEX   desc one\nSQ   Sequence 10 BP; 3 A; 3 C; 2 G; 2 T; 0 other;\n     ACGTACGTAC                                                               10\n//\nID   s2         standard; DNA; UNC; 10 BP.\nSQ   Sequence 10 BP; 2 A; 2 C; 2 G; 4 T; 0 other;\n     TTGGCCAATT                                                               10\n//\n"));
        // bam: the BGZF file of the installed program (blocks deflated by zlib 1.3), byte for byte
        var bam = System.Text.Encoding.Latin1.GetBytes(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "bam").SequenceText);
        Assert.That(Convert.ToHexString(bam).ToLowerInvariant(), Is.EqualTo("1f8b08040000000000ff06004243020061007372f4651467606070f070e10cf3b332d433e60cf6b72acd2bce2f2a494de102ca301801f17f2860060930b0307041c518a0a0d89041c8030819e180285d460c1d2e4a821d085d00ea1c83b28f0000001f8b08040000000000ff0600424302001b0003000000000000000000"));
    }

    [Test]
    public void Seqret_WritesTheNcbiDasAndAsn1Formats_AsTheProgram()
    {
        // seqret 6.6.0 sequence writers ncbi, gifasta, das, dasdna and asn1 (A5-8c3h, F264): the entries are the ones of the installed program
        const string fa = ">a desc\nACGTACGTAC\n";
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "ncbi").SequenceText, Is.EqualTo(">gnl|unk|a desc\nACGTACGTAC\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "gifasta").SequenceText, Is.EqualTo(">gi|000000|gnl|unk|a desc\nACGTACGTAC\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "dasdna").SequenceText, Is.EqualTo("<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE DASDNA SYSTEM \"http://www.biodas.org/dtd/dasdna.dtd\">\n<DASDNA>\n  <SEQUENCE id=\"a\" start=\"1\" stop=\"10\" version=\"\">\n    <DNA length=\"10\">\n      ACGTACGTAC\n    </DNA>\n  </SEQUENCE>\n</DASDNA>\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "das").SequenceText, Is.EqualTo("<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE DASSEQUENCE SYSTEM \"http://www.biodas.org/dtd/dassequence.dtd\">\n<DASSEQUENCE>\n  <SEQUENCE id=\"a\" start=\"1\" stop=\"10\"\n               moltype=\"DNA\" version=\"0.0\">\n      ACGTACGTAC\n  </SEQUENCE>\n</DASSEQUENCE>\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "asn1").SequenceText, Is.EqualTo("  seq {\n    id { local id 1 },\n    descr { title \"desc\" },\n    inst {\n      repr raw, mol dna, length 10, topology linear,\n      seq-data\n        iupacna \"ACGTACGTAC\"\n      } } ,\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "nib").SequenceText, Is.EqualTo("\u003A\u00ED\u00E9\u006B\u000A\0\0\0\u0021\u0030\u0021\u0030\u0021"));
    }

    [Test]
    public void Seqret_WritesTheEntryQualityAndSetFormats_AsTheProgram()
    {
        // seqret 6.6.0 sequence writers embl, fastq, phylip, clustal and treecon (A5-8c3h, WP59 F265 - F271): the texts are the ones of the installed program
        const string fa = ">a desc\nACGTACGTAC\n>b\nACGTA\n";
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "embl").SequenceText, Is.EqualTo("ID   a; SV 1; linear; unassigned DNA; STD; UNC; 10 BP.\nXX\nDE   desc\nXX\nSQ   Sequence 10 BP; 3 A; 3 C; 2 G; 2 T; 0 other;\n     ACGTACGTAC                                                               10\n//\nID   b; SV 1; linear; unassigned DNA; STD; UNC; 5 BP.\nXX\nSQ   Sequence 5 BP; 2 A; 1 C; 1 G; 1 T; 0 other;\n     ACGTA                                                                     5\n//\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "fastq-sanger").SequenceText, Is.EqualTo("@a desc\nACGTACGTAC\n+\n\"\"\"\"\"\"\"\"\"\"\n@b\nACGTA\n+\n\"\"\"\"\"\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "phylip").SequenceText, Is.EqualTo(" 2 10\na         ACGTACGTAC\nb         ACGTA-----\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "clustal").SequenceText, Is.EqualTo("CLUSTAL W (1.83) multiple sequence alignment\n\n\na               ACGTACGTAC\nb               ACGTA-----\n                          \n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "treecon").SequenceText, Is.EqualTo("10\na\nACGTACGTAC\nb\nACGTA\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, outputFormat: "fasta", outputDbName: "MYDB").SequenceText, Is.EqualTo(">MYDB:a desc\nACGTACGTAC\n>MYDB:b\nACGTA\n"));
    }

    [Test]
    public void Seqret_ReadsTheUsaForms_AsTheProgram()
    {
        // seqret 6.6.0 USA forms (A5-8c3h, F255): ranges, strands, list files, queries on files, -sbegin / -send / -sreverse, -sid; the texts are the ones of the installed program
        const string fa = ">a desc\nACGTACGTAAACCCGGGTTT\n>b\nGGGGCCCCAAAATTTT\n";
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, sequenceValue: "sequence.dat[3:10]").SequenceText, Is.EqualTo(">a desc\nGTACGTAA\n>b\nGGCCCCAA\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, sequenceValue: "sequence.dat[3:10:r]").SequenceText, Is.EqualTo(">a desc\nTTACGTAC\n>b\nTTGGGGCC\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, qualifiers: "-sbegin1 3 -send1 10 -sreverse1").SequenceText, Is.EqualTo(">a desc\nTTACGTAC\n>b\nTTGGGGCC\n"));
        var files = new Dictionary<string, string> { ["l.lst"] = "sequence.dat[2:5]\nfasta::sequence.dat[1:3:r] # comment\n" };
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, sequenceValue: "@l.lst", files: files).SequenceText, Is.EqualTo(">a desc\nCGTA\n>b\nGGGGCCCCAAAATTTT\n>a desc\nCGT\n>b\nCCC\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, sequenceValue: "sequence.dat:b").SequenceText, Is.EqualTo(">b\nGGGGCCCCAAAATTTT\n"));
        Assert.That(EmbossSequenceFormatsTools.Seqret(fa, qualifiers: "-sid1 zzz").SequenceText, Is.EqualTo(">zzz desc\nACGTACGTAAACCCGGGTTT\n>zzz\nGGGGCCCCAAAATTTT\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Seqret(fa, sequenceValue: "missing.dat"), Throws.ArgumentException);
    }

    [Test]
    public void Entret_ReturnsTheTextOfTheEntries_AsTheProgram()
    {
        // entret 6.6.0 (A5-8c3h, F257): the text of every entry as the reader stored it, a query on the file, -firstonly
        const string fa = ">s1 desc one\nACGTACGTAC\n>s2\nTTGGCCAATT\n";
        Assert.That(EmbossSequenceFormatsTools.Entret(fa).Text, Is.EqualTo(fa));
        Assert.That(EmbossSequenceFormatsTools.Entret(fa, sequenceValue: "sequence.dat:s2").Text, Is.EqualTo(">s2\nTTGGCCAATT\n"));
        Assert.That(EmbossSequenceFormatsTools.Entret(fa, firstOnly: true).Entries, Is.EqualTo(1));
        const string embl = "ID   X1; SV 1; linear; mRNA; STD; HUM; 10 BP.\nXX\nSQ   Sequence 10 BP;\n     acgtacgtac                                                    10\n//\n";
        Assert.That(EmbossSequenceFormatsTools.Entret(embl).Text, Is.EqualTo(embl));
        const string defs = "DBNAME tdb [ type: N format: fasta method: direct directory: \"d\" filename: \"*.fa\" ]\n";
        var files = new Dictionary<string, string> { ["d/a.fa"] = fa };
        Assert.That(EmbossSequenceFormatsTools.Entret("x", sequenceValue: "tdb:s1", files: files, definitions: defs).Text, Is.EqualTo(">s1 desc one\nACGTACGTAC\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Entret(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqret_ReadsTheIndexOfTheMethodEmboss_AsTheProgram()
    {
        // B09 WP60 (F275): a database of the method emboss: the B+tree index files of dbxfasta are searched (seqret -sequence tdb:seq2 / tdb:seq* with the installed 6.6.0 gives the same entries)
        const string definitions = "RES tdb [ type: Index idlen: 15 acclen: 15 svlen: 15 keylen: 15 deslen: 15 orglen: 15 ]\nSET PAGESIZE 2048\nDBNAME tdb [ type: N format: fasta method: emboss directory: \"db\" indexdirectory: \"idx\" fields: \"id\" hasaccession: N filename: \"*.fa\" ]\n";
        var one = EmbossSequenceFormatsTools.Seqret("x", sequenceValue: "tdb:seq2", files: DbxFiles(), definitions: definitions);
        Assert.That(one.SequenceText, Is.EqualTo(">seq2 second one\nGGGGCCCCAAAATTTT\n"));
        var all = EmbossSequenceFormatsTools.Seqret("x", sequenceValue: "tdb:seq*", files: DbxFiles(), definitions: definitions);
        Assert.That(all.SequenceText, Is.EqualTo(">seq1 first one\nACGTACGTACGTAAAA\n>seq2 second one\nGGGGCCCCAAAATTTT\n"));
        Assert.That(EmbossSequenceFormatsTools.Entret("x", sequenceValue: "tdb:seq1", files: DbxFiles(), definitions: definitions).Text, Is.EqualTo(">seq1 first one\nACGTACGTACGTAAAA\n"));
    }

    [Test]
    public void Union_WritesTheFileOfOutputSingle_AsTheProgram()
    {
        // union 6.6.0 (B09 WP60, F277): -ossingle2 -osdirectory2 -osextension2: one file, name.extension in lower case, in the directory (measured on the installed program)
        var r = EmbossSequenceFormatsTools.Union(EmblWithFeatures, outputSingle: true, outputDirectory: "d", outputExtension: "txt");
        Assert.That(r.SingleFiles, Is.Not.Null);
        Assert.That(r.SingleFiles!.Keys, Is.EqualTo(new[] { "d/x1.txt" }));
        Assert.That(r.SingleFiles["d/x1.txt"], Is.EqualTo(">X1 X1.1 test\nacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n"));
        Assert.That(r.SequenceText, Is.Empty);
    }

    [Test]
    public void Union_JoinsTheFeatures_AsTheProgram()
    {
        // union 6.6.0 -feature -source (A5-8c3h, F261): the features of the entries moved by the length of the entries before, one source feature per entry (tag origid)
        const string two = "ID   Y1; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\nXX\nAC   Y1;\nXX\nDE   first\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             3..12\nFT                   /gene=\"a\"\nSQ   Sequence 20 BP; 5 A; 5 C; 5 G; 5 T; 0 other;\n     acgtacgtac gtacgtacgt                                                    20\n//\n"
            + "ID   Y2; SV 1; linear; genomic DNA; STD; UNC; 20 BP.\nXX\nAC   Y2;\nXX\nDE   second\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   exon            complement(1..5)\nSQ   Sequence 20 BP; 5 A; 5 C; 5 G; 5 T; 0 other;\n     gtacgtacgt gtacgtacgt                                                    20\n//\n";
        var r = EmbossSequenceFormatsTools.Union(two, feature: true, source: true, featureFormat: "embl");
        Assert.That(r.SequenceText, Is.EqualTo(">Y1 Y1.1 first\nacgtacgtacgtacgtacgtgtacgtacgtgtacgtacgt\n"));
        Assert.That(r.FeatureText, Is.EqualTo("FH   Key             Location/Qualifiers\nFH\nFT   source          1..20\nFT                   /note=\"*origid: Y1\"\nFT   CDS             3..12\nFT                   /gene=\"a\"\nFT   source          21..40\nFT                   /note=\"*origid: Y2\"\nFT   exon            complement(21..25)\n"));
        Assert.That(EmbossSequenceFormatsTools.Union(two).FeatureText, Is.Empty);
    }

    [Test]
    public void Union_JoinsTheSequences_AsTheProgram()
    {
        // union 6.6.0 (A5-8c3h, F253)
        const string fa = ">f1\nACGTACGTAC\n>f2\nCGTACGTACGGG\n>f3\nGGGTTTAAA\n";
        var r = EmbossSequenceFormatsTools.Union(fa, findOverlap: true, overlapFile: true);
        Assert.That(r.SequenceText, Is.EqualTo(">f1\nACGTACGTACGGGTTTAAA\n"));
        Assert.That(r.OverlapText, Is.EqualTo("0\n9\n3\n"));
        Assert.That(EmbossSequenceFormatsTools.Union(fa).SequenceText, Is.EqualTo(">f1\nACGTACGTACCGTACGTACGGGGGGTTTAAA\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Union(" "), Throws.ArgumentException);
    }

    [Test]
    public void Seqretsplit_WritesOneFilePerSequence_AsTheProgram()
    {
        // seqretsplit 6.6.0 (A5-8c3h, F254)
        var r = EmbossSequenceFormatsTools.Seqretsplit(">a\nACGT\n>A\nGG\n>b\nTT\n");
        Assert.That(r.SingleFiles.Keys, Is.EquivalentTo(new[] { "a.fasta", "a.002.fasta", "b.fasta" }));
        Assert.That(r.SingleFiles["a.002.fasta"], Is.EqualTo(">A.002\nGG\n"));
        Assert.That(r.Messages, Does.Contain("Warning: Duplicate name 'A' changed to 'A.002'"));
        Assert.That(EmbossSequenceFormatsTools.Seqretsplit(">a\nACGT\n>b\nTT\n", outputSingle: false).SequenceText, Is.EqualTo(">a\nACGT\n>b\nTT\n"));
        Assert.That(() => EmbossSequenceFormatsTools.Seqretsplit(" "), Throws.ArgumentException);
    }

    private const string EmblWithFeatures = "ID   X1; SV 1; linear; genomic DNA; STD; UNC; 60 BP.\nXX\nAC   X1;\nXX\nDE   test\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             join(3..12,20..29)\nFT                   /translation=\"MKV\"\nFT                   /gene=\"g1\"\nFT   mRNA            1..40\nFT   exon            45..55\nSQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n//\n";

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

    // the index files dbxfasta 6.6.0 made for a FASTA file of two entries (the field id, compressed pages), gzip + base64 of their bytes
    private static readonly Dictionary<string, string> DbxIndexFiles = new()
    {
        ["idx/tdb.ent"] = "H4sIAAAAAAAC/1NW8CvNTUotUshPU0jLzEkttlIw5FJWCErNSU0sTrVSMNAzAHJdEkuAbCAwMNAHI67gzLz0nFSwlrzE3FSFlMSSxCSgDq5EvbRELgC3kTO+VgAAAA==",
        ["idx/tdb.pxid"] = "H4sIAAAAAAAC/z2PywrCMBBF9/MV/YRkbNR9QRAFRd24LOlUB9IHSavo12vbZLLLublnJrdPT1k6+4ragWsmD0XX9J5CmPmdApzLB4X0cAVXsr0QBSdfkU/pRsOOnROtWcORXiRALTL+LpNR5VsoSvsUhCoaMTYQZ2O6ZrmWBeaK0TgBK5LJUHRjO6TGXzA6ZxNCOLwdNxxzbeBCdc1OPvQDgBCNYRgBAAA=",
        ["idx/tdb.xid"] = "H4sIAAAAAAAC/2NkQABGIGZhIAxAaopTCw0ZfKECsVCaA00hiM+ExBeFYrBmdAAUNIJzFBlGwSgYBaNgFIyCUTA0AQAD30AfXQgAAA==",
    };
}
