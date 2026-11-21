using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP <c>minimap2_map</c> (B09 A3-9c2 / A5-9b / A5-9c / A5-9d, F187 - F189): minimap2 2.26 as a program (<c>-c</c> / <c>-a</c>, paired reads, the options) on the same synthetic genome as <see cref="Minimap2PresetToolsTests"/>; the expected lines are
/// those of the installed program (<c>minimap2 [options] ref.fa query.fa</c>); the library tests compare 7 201 whole outputs.
/// </summary>
[TestFixture]
public class Minimap2MapToolTests
{
    // a 64-bit LCG (x <- x * 6364136223846793005 + 1442695040888963407, base = x >> 62), the same as the Python script that produced the expected lines
    private static string Genome(int n, int seed)
    {
        ulong x = (ulong)seed;
        var chars = new char[n];
        for (int i = 0; i < n; i++)
        {
            x = unchecked(x * 6364136223846793005UL + 1442695040888963407UL);
            chars[i] = "ACGT"[(int)(x >> 62)];
        }

        return new string(chars);
    }

    private static readonly string Target = Genome(6000, 1);
    private static readonly string Reference = ">t1\n" + Target + "\n";

    /// <summary>target[1000, 3500) with every 97th base changed (A becomes C, any other base A).</summary>
    private static readonly string Query = ">query\n" + string.Concat(Target[1000..3500].Select((c, i) => i % 97 == 0 ? (c != 'A' ? 'A' : 'C') : c)) + "\n";

    private static string Revcomp(string s) => string.Concat(s.Reverse().Select(c => c switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', _ => 'A' }));

    private static readonly string Read1 = string.Concat(">pair/1\n", Target.AsSpan(2000, 150), "\n");
    private static readonly string Read2 = ">pair/2\n" + Revcomp(Target[2300..2450]) + "\n";

    private static string[] Columns(string line, params int[] columns) => [.. columns.Select(c => line.Split('\t')[c])];

    [Test]
    public void BaseLevelAlignment_Paf_PrintsTheCigarAndTheAlignmentTags()
    {
        var r = AlignmentTools.Minimap2Map(Reference, [Query], ["-c", "--MD"]);
        Assert.That(r.ExitCode, Is.EqualTo(0));
        string[] f = r.Output.TrimEnd('\n').Split('\t');
        // minimap2 -c --MD ref.fa query.fa: the first base is a mismatch and clipped (the alignment starts at query position 1), then 2499 columns of M with 25 differences
        Assert.That(f[..21], Is.EqualTo(new[]
        {
            "query", "2500", "1", "2500", "+", "t1", "6000", "1001", "3500", "2474", "2499", "60", "NM:i:25", "ms:i:4848", "AS:i:4848", "nn:i:0", "tp:A:P", "cm:i:386", "s1:i:2302", "s2:i:0", "de:f:0.0100",
        }));
        Assert.That(f[21..23], Is.EqualTo(new[] { "rl:i:0", "cg:Z:2499M" }));
        Assert.That(f[23], Does.StartWith("MD:Z:"));
    }

    [Test]
    public void Sam_PrintsTheHeaderWithTheCommandLineAndTheClippedAlignment()
    {
        var r = AlignmentTools.Minimap2Map(Reference, [Query], ["-a"]);
        string[] lines = r.Output.TrimEnd('\n').Split('\n');
        Assert.That(lines[..3], Is.EqualTo(new[] { "@HD\tVN:1.6\tSO:unsorted\tGO:query", "@SQ\tSN:t1\tLN:6000", "@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -a ref.fa query.fa" }));
        Assert.That(Columns(lines[3], 0, 1, 2, 3, 4, 5), Is.EqualTo(new[] { "query", "0", "t1", "1002", "60", "1S2499M" }));
    }

    [Test]
    public void PairedReads_TwoFiles_AreAProperPair()
    {
        var sam = AlignmentTools.Minimap2Map(Reference, [Read1, Read2], ["-x", "sr", "-a"]);
        string[] lines = sam.Output.TrimEnd('\n').Split('\n');
        Assert.That(lines[2], Does.EndWith("CL:minimap2 -x sr -a ref.fa query1.fa query2.fa"));
        // flags 99 / 147 (paired, proper, mate reversed / reversed, first / second), mate position, template length +-450, mapq 60
        Assert.That(Columns(lines[3], 0, 1, 2, 3, 4, 5, 6, 7, 8), Is.EqualTo(new[] { "pair", "99", "t1", "2001", "60", "150M", "=", "2301", "450" }));
        Assert.That(Columns(lines[4], 0, 1, 2, 3, 4, 5, 6, 7, 8), Is.EqualTo(new[] { "pair", "147", "t1", "2301", "60", "150M", "=", "2001", "-450" }));
        var paf = AlignmentTools.Minimap2Map(Reference, [Read1, Read2], ["-x", "sr", "-c"]);
        string[] p = paf.Output.TrimEnd('\n').Split('\n');
        Assert.That(Columns(p[0], 0, 4, 7, 8, 11), Is.EqualTo(new[] { "pair/1", "+", "2000", "2150", "60" }));
        Assert.That(Columns(p[1], 0, 4, 7, 8, 11), Is.EqualTo(new[] { "pair/2", "-", "2300", "2450", "60" }));
    }

    [Test]
    public void InvalidOptions_ReturnTheProgramsExitStatusAndMessage()
    {
        var r = AlignmentTools.Minimap2Map(Reference, [Query], ["-v", "2", "-k", "0"]);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Output, Is.Empty);
        Assert.That(r.StandardError, Is.EqualTo("[ERROR]\u001b[1;31m -k and -w must be positive\u001b[0m\n"));
        Assert.That(AlignmentTools.Minimap2Map(Reference, [Query], ["--splice", "--frag=yes"]).ExitCode, Is.EqualTo(1));
        Assert.That(AlignmentTools.Minimap2Map(Reference, [Query], ["-V"]).Output, Is.EqualTo("2.26-r1175\n"));
    }

    [Test]
    public void OutputFile_IsReturnedInFilesAndStandardOutputStaysEmpty()
    {
        var r = AlignmentTools.Minimap2Map(Reference, [Query], ["-c", "-o", "hits.paf"]);
        Assert.That(r.Output, Is.Empty);
        Assert.That(r.Files!["hits.paf"], Does.StartWith("query\t2500\t1\t2500\t+\tt1\t6000\t1001\t3500"));
    }

    [Test]
    public void IndexInParts_SamHeaderHasNoSequenceLines()
    {
        // -I 3k: the 6 kb target is one sequence (a part ends only between sequences), so it stays one part; two targets split into two parts
        string two = Reference + ">t2\n" + Genome(2000, 7) + "\n";
        var r = AlignmentTools.Minimap2Map(two, [Query], ["-a", "-I", "3k", "-v", "2"]);
        Assert.That(r.Output, Does.Not.Contain("@SQ"));
        Assert.That(r.StandardError, Does.Contain("For a multi-part index, no @SQ lines will be outputted"));
    }

    [Test]
    public void Arguments_AreValidated()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Map("", [Query]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Map(Reference, []));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Map(Reference, [Query, Query, Query, Query]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Map(Reference, [""]));
    }

    [Test]
    public void Ksw2Extend_ExtendS2_AndLocalScore_Variants()
    {
        string e1 = "ACGTTGCAGGCTAAGCTTCGAGGATCCTAG", e2 = "TTGACCGATAGGCTCATCGGATCCAGTCAA";
        var spliced = AlignmentTools.Ksw2Extend(e1 + e2, e1 + "GT" + new string('A', 96) + "AG" + e2, "exts2", gapOpen: 4, gapExtend: 2, gapOpen2: 24, nonCanonical: 9, junctionBonus: 9, spliceStrand: "forward");
        Assert.That(spliced.Cigar, Is.EqualTo("30M100N30M"));
        var local = AlignmentTools.Ksw2Extend("ACGTACGTTT", "ACGTACGTTT", "ll", gapOpen: 4, gapExtend: 2);
        Assert.That(local.Max, Is.EqualTo(20));
        Assert.That(local.MaxT, Is.EqualTo(9));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Ksw2Extend("ACGT", "ACGT", "exts2", spliceStrand: "sideways"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Ksw2Extend("ACGT", "ACGT", "exts2", junctions: "00"));
    }
}
