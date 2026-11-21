using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A2-28d (B09 F203): EMBOSS 6.6.0 <c>est2genome</c> (TestData/Alignment/est2genome_oracle.jsonl.gz, scripts/alignment/generate_est2genome_oracle.py): the -outfile text of spliced EST alignments to genomic sequences with forward / reverse
/// genes, ESTs on either strand, all scoring options, -nousesplice, -mode, -nobest, -space (small values force the linear-space recursion), -shuffle, -seed, -align / -width and the sequence qualifiers, including the runs on which the program crashes.
/// </summary>
[TestFixture]
public class EmbossEst2genome_Oracle_Tests
{
    [Test]
    public void Est2genome_OutfileEqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.est2genome_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, crashes = 0, align = 0, nobest = 0, recursion = 0, reverseGene = 0, forwardEst = 0, reversedEst = 0, noSplice = 0, qualified = 0, shuffled = 0, warnings = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string est = j.GetProperty("est").GetString()!, genome = j.GetProperty("genome").GetString()!, err = j.GetProperty("err").GetString()!;
            string? expected = j.GetProperty("out").ValueKind == JsonValueKind.Null ? null : j.GetProperty("out").GetString();
            var got = Est2genomeOracleRunner.Run(args, est, genome);
            string label = $"row {rows} [{string.Join(' ', args)}]";
            if (expected is null)
            {
                crashes++;
                if (got is not null)
                    failures.Add(label + ": the program fails (rc " + j.GetProperty("rc").GetInt32() + "), the port returns text");
                continue;
            }

            if (got is null || got.Text != expected)
            {
                failures.Add($"{label}: {(got is null ? "<error>" : MummerShowOracleRunner.Diff(got.Text, expected))}");
                continue;
            }

            if (err.Contains("increasing space", StringComparison.Ordinal) != got.Messages.Contains("increasing space", StringComparison.Ordinal))
                failures.Add(label + ": the space warning differs");
            if (args.Contains("-align")) align++;
            if (args.Contains("-nobest")) nobest++;
            if (args.Contains("-nousesplice")) noSplice++;
            if (args.Any(a => a.StartsWith("-s", StringComparison.Ordinal) && a.Length > 5 && a[^1] is '1' or '2')) qualified++;
            if (args.Contains("-shuffle")) shuffled++;
            if (err.Contains("increasing space", StringComparison.Ordinal)) warnings++;
            if (expected.Contains("REVERSED GENE", StringComparison.Ordinal)) reverseGene++;
            if (expected.Contains("between forward est", StringComparison.Ordinal)) forwardEst++;
            if (expected.Contains("reversed est", StringComparison.Ordinal)) reversedEst++;
            int space = args.IndexOf("-space");
            if (space >= 0 && float.Parse(args[space + 1], System.Globalization.CultureInfo.InvariantCulture) <= 0.05f)
                recursion++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(300));
        Assert.That(crashes, Is.GreaterThan(5));
        Assert.That(new[] { align, nobest, recursion, reverseGene, forwardEst, reversedEst, noSplice, qualified, shuffled, warnings }, Has.All.GreaterThan(5));
    }

    [Test]
    public void Est2genome_WorkedExample_PrintsExonsIntronsAndSpan()
    {
        // est2genome -estsequence e.fa -genomesequence g.fa (EMBOSS 6.6.0): three exons, two forward-spliced introns
        const string exon1 = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTT";
        const string exon2 = "CAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGT";
        const string exon3 = "TCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCT";
        string genome = "ACGT" + exon1 + "GT" + new string('A', 70) + "AG" + exon2 + "GT" + new string('C', 40) + "AG" + exon3 + "TTTT";
        var o = EmbossEst2genomePrograms.Run(new EmbossEst2genomeRequest(">est1 an est\n" + exon1 + exon2 + exon3 + "\n", ">genomic some description\n" + genome + "\n"));
        Assert.That(o.Text, Is.EqualTo("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\nExon        60 100.0     5    64 genomic          1    60 est1          an est\n+Intron    -20   0.0    65   138 genomic     \nExon        80 100.0   139   218 genomic         61   140 est1          an est\n+Intron    -20   0.0   219   262 genomic     \nExon        50 100.0   263   312 genomic        141   190 est1          an est\n\nSpan       150 100.0     5   312 genomic          1   190 est1          an est\n\nSegment     60 100.0     5    64 genomic          1    60 est1          an est\nSegment     80 100.0   139   218 genomic         61   140 est1          an est\nSegment     50 100.0   263   312 genomic        141   190 est1          an est\n"));
        Assert.That(o.Messages, Is.Empty);
        Assert.Throws<ArgumentException>(() => EmbossEst2genomePrograms.Run(new EmbossEst2genomeRequest(">e\nACGT\n", ">g\nACGT\n", Mode: "x")));
        Assert.Throws<ArgumentException>(() => EmbossEst2genomePrograms.Run(new EmbossEst2genomeRequest(">e\nACGT\n", ">g\nACGT\n", Align: true, Width: 0)));
        Assert.Throws<InvalidOperationException>(() => EmbossEst2genomePrograms.Run(new EmbossEst2genomeRequest(">e\n" + exon1 + "\n", ">g\n" + genome + "\n", Mode: "forward", Best: false)));
    }
}
