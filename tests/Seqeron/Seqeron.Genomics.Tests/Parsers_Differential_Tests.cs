// 08_DIFFERENTIAL_TESTING rows 64-68 (FileIO parsers). Each production parser is checked against an
// INDEPENDENT manual line/tab-split parse of the same content (a different parsing strategy), asserting
// the same records.

using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Parsers_Differential_Tests
{
    // ---- Row 64: PARSE-FASTA-001 — FastaParser vs manual '>'-block parse ----

    [Test]
    [Category("PARSE-FASTA-001")]
    public void FastaParser_MatchesManualBlockParse()
    {
        const string content = ">seq1 first sequence\nACGT\nAACC\n>seq2\nGGGGTTTT\n";
        var actual = FastaParser.Parse(content)
            .Select(e => (e.Id, e.Description, Seq: e.Sequence.Sequence)).ToList();

        var expected = new List<(string, string?, string)>();
        foreach (var block in content.Split('>', StringSplitOptions.RemoveEmptyEntries))
        {
            var lines = block.Split('\n');
            var hdr = lines[0].Trim().Split(new[] { ' ', '\t' }, 2);
            string id = hdr[0];
            string? desc = hdr.Length > 1 ? hdr[1] : null;
            string seq = string.Concat(lines.Skip(1).Select(l => new string(l.Where(c => !char.IsWhiteSpace(c)).ToArray()))).ToUpperInvariant();
            if (seq.Length > 0) expected.Add((id, desc, seq));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 65: PARSE-FASTQ-001 — FastqParser vs manual 4-line block parse ----

    [Test]
    [Category("PARSE-FASTQ-001")]
    public void FastqParser_MatchesManual4LineBlocks()
    {
        const string content = "@r1 read one\nACGT\n+\nIIII\n@r2\nGGGG\n+\n!!!!\n";
        var actual = FastqParser.Parse(content)
            .Select(r => (r.Id, r.Description, r.Sequence, r.QualityString)).ToList();

        var lines = content.Split('\n').Where(l => l.Length > 0).ToList();
        var expected = new List<(string, string, string, string)>();
        for (int i = 0; i + 3 < lines.Count + 1 && i + 3 <= lines.Count - 1; i += 4)
        {
            var hdr = lines[i][1..];
            int sp = hdr.IndexOf(' ');
            string id = sp > 0 ? hdr[..sp] : hdr;
            string desc = sp > 0 ? hdr[(sp + 1)..] : "";
            expected.Add((id, desc, lines[i + 1], lines[i + 3]));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 66: PARSE-BED-001 — BedParser vs manual tab-split ----

    [Test]
    [Category("PARSE-BED-001")]
    public void BedParser_MatchesManualTabSplit()
    {
        const string content = "chr1\t100\t200\tfeatA\t500\t+\nchr2\t50\t75\tfeatB\t0\t-\n";
        var actual = BedParser.Parse(content)
            .Select(b => (b.Chrom, b.ChromStart, b.ChromEnd, b.Name, b.Strand)).ToList();

        var expected = new List<(string, int, int, string?, char?)>();
        foreach (var line in content.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            var f = line.Split('\t');
            expected.Add((f[0], int.Parse(f[1]), int.Parse(f[2]),
                f.Length > 3 ? f[3] : null, f.Length > 5 ? f[5][0] : (char?)null));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 67: PARSE-VCF-001 — VcfParser vs manual tab-split ----

    [Test]
    [Category("PARSE-VCF-001")]
    public void VcfParser_MatchesManualTabSplit()
    {
        const string content =
            "##fileformat=VCFv4.2\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\n" +
            "chr2\t200\t.\tC\tT,A\t.\tq10\tNS=3\n";
        var actual = VcfParser.Parse(content)
            .Select(v => (v.Chrom, v.Pos, v.Ref, Alt: string.Join(",", v.Alt))).ToList();

        var expected = new List<(string, int, string, string)>();
        foreach (var line in content.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            if (line.StartsWith("#")) continue;
            var f = line.Split('\t');
            expected.Add((f[0], int.Parse(f[1]), f[3], f[4]));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }

    // ---- Row 68: PARSE-GFF-001 — GffParser vs manual tab-split ----

    [Test]
    [Category("PARSE-GFF-001")]
    public void GffParser_MatchesManualTabSplit()
    {
        const string content =
            "##gff-version 3\n" +
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n" +
            "chr2\tsrc\tCDS\t10\t50\t0.5\t-\t0\tID=c1\n";
        var actual = GffParser.Parse(content)
            .Select(g => (g.Seqid, g.Source, g.Type, g.Start, g.End, g.Strand)).ToList();

        var expected = new List<(string, string, string, int, int, char)>();
        foreach (var line in content.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            if (line.StartsWith("#")) continue;
            var f = line.Split('\t');
            expected.Add((f[0], f[1], f[2], int.Parse(f[3]), int.Parse(f[4]), f[6][0]));
        }
        Assert.That(actual, Is.EqualTo(expected));
    }
}
