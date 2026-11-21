using System.IO.Compression;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3e (B09 F241): the feature readers of EMBOSS 6.6.0 (<c>ajfeatread.c</c>: embl, genbank, refseq, genpept, refseqp, gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff) read every feature table of a file as the
/// installed <c>featcopy -features FILE -outfeat out -offormat debug</c> does (the debug trace prints type, source, score, strand, positions, flags, group, exon, tags and sub-features of every feature):
/// TestData/Alignment/featread_oracle.jsonl.gz, scripts/alignment/generate_featread_oracle.py; no format (all the formats the library tries for an unknown input), the <c>format::file</c> prefix (a format that does not read the
/// file, an unknown name), <c>-fbegin1</c> / <c>-fend1</c> (the trimmed first table), damaged and mixed files, the group field the GFF2 reader loops on forever.
/// </summary>
[TestFixture]
public class EmbossFeatureReaders_Oracle_Tests
{
    // ajFeattableWrite of the debug format
    private static string WriteDebug(EmbossSeqFeatureTable table)
    {
        if (table.Rev || table.GetBegin() > 1 || table.GetEnd() < table.Length)
        {
            var copy = table.NewFtable();
            copy.Trim();
            return copy.Debug();
        }

        return table.Debug();
    }

    internal static (int Rc, string Output, string Error) RunDebug(string text, string ufo, IReadOnlyList<string> extra)
    {
        int fbegin = 1;
        int? fend = null;
        for (int i = 0; i + 1 < extra.Count; i += 2)
        {
            if (extra[i] == "-fbegin1")
                fbegin = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
            else if (extra[i] == "-fend1")
                fend = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
        }

        var context = new EmbossFeatureReadContext();
        EmbossFeatureInput input;
        try
        {
            input = new EmbossFeatureInput(ufo, name => name == "in.dat" ? text : null, fbegin, fend, context);
        }
        catch (EmbossFeatureDiedException e)
        {
            return (1, "", string.Join("\n", [.. e.Errors, e.Message]) + "\n");
        }

        var sb = new StringBuilder();
        EmbossSeqFeatureTable? table;
        while ((table = input.Next()) is not null)
            sb.Append(WriteDebug(table));
        return (0, sb.ToString(), string.Join("\n", input.Errors) + (input.Errors.Count > 0 ? "\n" : ""));
    }

    [Test]
    public void FeatureReaders_EqualFeatcopyDebug_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.featread_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var perKind = new Dictionary<string, int>();
        var allKind = new Dictionary<string, int>();
        int rows = 0, hang = 0, died = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string kind = j.GetProperty("kind").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            string ufo = j.GetProperty("ufo").GetString()!;
            var extra = j.GetProperty("extra").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            string expected = j.GetProperty("debug").GetString()!;
            string err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} ({kind} {ufo} {string.Join(' ', extra)})";
            allKind[kind] = allKind.GetValueOrDefault(kind) + 1;
            if (j.GetProperty("hang").GetBoolean())
            {
                hang++;
                try
                {
                    RunDebug(text, ufo, extra);
                    failures.Add($"{id}: the program loops forever, the port returned");
                }
                catch (InvalidOperationException)
                {
                }

                continue;
            }

            try
            {
                var (grc, output, gerr) = RunDebug(text, ufo, extra);
                string[] eLines = [.. err.Split('\n').Where(l => l.Length > 0)];
                string[] gLines = [.. gerr.Split('\n').Where(l => l.Length > 0)];
                if (grc != rc)
                    failures.Add($"{id}: exit code {grc} expected {rc}");
                else if (output != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(output, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: errors '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
                else if (rc == 0 && expected.Length > 0)
                    perKind[kind] = perKind.GetValueOrDefault(kind) + 1;
                if (rc != 0)
                    died++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(2600));
        foreach (string kind in new[] { "embl", "genbank", "refseq", "genpept", "refseqp", "gff2", "gff3", "gff3emboss", "swiss", "pir", "dasgff" })
        {
            Assert.That(allKind.GetValueOrDefault(kind), Is.GreaterThanOrEqualTo(200), $"files of {kind}");
            Assert.That(perKind.GetValueOrDefault(kind), Is.GreaterThanOrEqualTo(80), $"read files of {kind}");
        }
        Assert.That(hang, Is.GreaterThan(5), "files on which the program never stops");
        Assert.That(died, Is.GreaterThan(100), "files the program refuses");
    }

    /// <summary>
    /// The features of a sequence read from a SEPARATE file (<c>-ufo</c> = <c>[format::]file</c>, <c>-fformat</c>, <c>-fopenfile</c>, a file named after the sequence: <c>ajFeattableNewReadUfo</c>) of every input format equal the trace of
    /// <c>seqret -feature -ufo1 ... -offormat2 debug</c> (TestData/Alignment/featread_ufo_oracle.jsonl.gz, <c>generate_featread_oracle.py --ufo</c>, 1 200 runs: no format = every format of the table is tried, a format that is none, a missing file).
    /// </summary>
    [Test]
    public void FeatureFileOfASequence_EqualsSeqretDebug_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.featread_ufo_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, hang = 0, read = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string ufo = "", fformat = "", fopen = "";
            for (int i = 0; i + 1 < args.Count; i += 2)
            {
                if (args[i] == "-ufo1") ufo = args[i + 1];
                else if (args[i] == "-fformat1") fformat = args[i + 1];
                else if (args[i] == "-fopenfile1") fopen = args[i + 1];
            }

            string seq = j.GetProperty("seq").GetString()!;
            string expected = j.GetProperty("debug").GetString()!;
            string[] eLines = [.. j.GetProperty("err").GetString()!.Split('\n').Where(l => l.Length > 0 && !l.StartsWith("Warning: No features written", StringComparison.Ordinal))];
            string id = $"row {j.GetProperty("id")} ({j.GetProperty("kind").GetString()} {string.Join(' ', args)})";
            var context = new EmbossFeatureReadContext { Program = "seqret" };
            try
            {
                var table = EmbossFeatureUfo.Read(ufo, fformat, fopen, j.GetProperty("name").GetString()!, n => files.GetValueOrDefault(n), context);
                if (j.GetProperty("hang").GetBoolean())
                {
                    hang++;
                    failures.Add($"{id}: the program never stops, the port returned");
                    continue;
                }

                string got = "";
                if (table is not null)
                {
                    table.Length = seq.Length;
                    got = WriteDebug(table);
                    read++;
                }

                string[] gLines = [.. context.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0)];
                if (got != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(got, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
            }
            catch (InvalidOperationException) when (j.GetProperty("hang").GetBoolean())
            {
                hang++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1200));
        Assert.That(read, Is.GreaterThan(600), "files that gave features");
        Assert.That(hang, Is.GreaterThan(5));
    }

    private static (int Rows, Dictionary<string, int> PerFormat, int Hang) RunWriteOracle(string name)
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var perFormat = new Dictionary<string, int>();
        int rows = 0, hang = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string text = j.GetProperty("text").GetString()!;
            var extra = j.GetProperty("extra").EnumerateArray().Select(x => x.GetString()!).ToList();
            string format = j.GetProperty("format").GetString()!;
            int fbegin = 1;
            int? fend = null;
            for (int i = 0; i + 1 < extra.Count; i += 2)
            {
                if (extra[i] == "-fbegin1") fbegin = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
                else if (extra[i] == "-fend1") fend = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
            }

            string expected = j.GetProperty("out").GetString()!;
            string id = $"row {j.GetProperty("id")} ({j.GetProperty("kind").GetString()} {format} {string.Join(' ', extra)})";
            var request = new EmbossFeatcopyRequest(j.GetProperty("ufo").GetString()!, new Dictionary<string, string> { ["in.dat"] = text }, fbegin, fend, "out.dat", format, null,
                DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime);
            string[] eLines = [.. j.GetProperty("err").GetString()!.Split('\n').Where(l => l.Length > 0)];
            try
            {
                var got = EmbossFeaturePrograms.Featcopy(request);
                if (j.GetProperty("hang").GetBoolean())
                {
                    failures.Add($"{id}: the program never stops, the port returned");
                    continue;
                }

                string[] gLines = [.. got.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 0)
                    failures.Add($"{id}: the program died, the port returned");
                else if (got.Text != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(got.Text, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
                else if (expected.Length > 0)
                    perFormat[format.ToLowerInvariant()] = perFormat.GetValueOrDefault(format.ToLowerInvariant()) + 1;
            }
            catch (EmbossFeatureDiedException e)
            {
                string[] gLines = [.. e.Errors.SelectMany(m => m.Split('\n')).Append(e.Message).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 1)
                    failures.Add($"{id}: the port died, the program did not");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
            }
            catch (InvalidOperationException) when (j.GetProperty("hang").GetBoolean())
            {
                hang++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        return (rows, perFormat, hang);
    }

    /// <summary>
    /// <c>featcopy</c> (<see cref="EmbossFeaturePrograms.Featcopy"/>): every table of a file in any input format written in the output formats gff3 (gff, unknown), gff2, embl (em), genbank (gb, ddbj) of the feature library - sub-features of join / order locations,
    /// flags, remote ids, GFF3 attributes, the type notes of the ontology, the group order, the table range (<c>-fbegin1</c> / <c>-fend1</c>), the formats the program refuses for a protein table - equal the installed
    /// <c>featcopy -features FILE -outfeat out.dat -offormat F</c> byte for byte, with its messages (TestData/Alignment/featwrite_oracle.jsonl.gz, <c>generate_featread_oracle.py --write</c>, dates faked).
    /// </summary>
    [Test]
    public void Featcopy_WritesEveryFormat_EqualsTheProgram_AllRows()
    {
        var (rows, perFormat, hang) = RunWriteOracle("Seqeron.Genomics.Tests.TestData.Alignment.featwrite_oracle.jsonl.gz");
        Assert.That(rows, Is.GreaterThanOrEqualTo(1800));
        foreach (string f in new[] { "gff3", "gff2", "embl", "genbank" })
            Assert.That(perFormat.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(150), $"tables written in {f}");
        Assert.That(perFormat.GetValueOrDefault("gff") + perFormat.GetValueOrDefault("em") + perFormat.GetValueOrDefault("gb") + perFormat.GetValueOrDefault("ddbj"), Is.GreaterThan(100), "the aliases");
        Assert.That(hang, Is.GreaterThan(3));
    }

    /// <summary>
    /// A5-8c3f (F244): the feature writers of bed, refseq, refseqp, pir (nbrf), swiss (sw, swissprot, uniprot), dasgff and draw for a full table equal <c>featcopy -features FILE -outfeat out.dat -offormat F</c>
    /// (TestData/Alignment/featwrite2_oracle.jsonl.gz, <c>generate_featread_oracle.py --write2</c>, 250 runs per format; the runs the program never ends or whose file it fills are rows that must be refused).
    /// </summary>
    [Test]
    public void Featcopy_WritesBedRefseqRefseqpPirSwissDasgffDraw_EqualsTheProgram_AllRows()
    {
        var (rows, perFormat, hang) = RunWriteOracle("Seqeron.Genomics.Tests.TestData.Alignment.featwrite2_oracle.jsonl.gz");
        Assert.That(rows, Is.GreaterThanOrEqualTo(2750));
        foreach (string f in new[] { "bed", "refseq", "refseqp", "pir", "nbrf", "swiss", "sw", "swissprot", "uniprot", "dasgff", "draw" })
            Assert.That(perFormat.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(20), $"tables written in {f}");
        Assert.That(hang, Is.GreaterThan(3));
    }

    /// <summary>
    /// A5-8c3f (F245): <c>featmerge</c> (<see cref="EmbossFeaturePrograms.Featmerge"/>) equals <c>featmerge -afeatures A -bfeatures B -outfeat out.dat -offormat F</c> byte for byte with its messages, over every pair of input formats,
    /// ranges of both inputs and every output format (TestData/Alignment/featmerge_oracle.jsonl.gz, <c>generate_featread_oracle.py --merge</c>).
    /// </summary>
    [Test]
    public void Featmerge_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.featmerge_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, hang = 0, written = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var extra = j.GetProperty("extra").EnumerateArray().Select(x => x.GetString()!).ToList();
            var range = new Dictionary<string, int> { ["-fbegin1"] = 1, ["-fbegin2"] = 1 };
            var ends = new Dictionary<string, int?> { ["-fend1"] = null, ["-fend2"] = null };
            for (int i = 0; i + 1 < extra.Count; i += 2)
            {
                if (range.ContainsKey(extra[i])) range[extra[i]] = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
                else ends[extra[i]] = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
            }

            string format = j.GetProperty("format").GetString()!;
            string expected = j.GetProperty("out").GetString()!;
            string id = $"row {j.GetProperty("id")} ({j.GetProperty("kind").GetString()} {format} {string.Join(' ', extra)})";
            var request = new EmbossFeatmergeRequest(j.GetProperty("ufoa").GetString()!, j.GetProperty("ufob").GetString()!,
                new Dictionary<string, string> { ["a.dat"] = j.GetProperty("a").GetString()!, ["b.dat"] = j.GetProperty("b").GetString()! },
                range["-fbegin1"], ends["-fend1"], range["-fbegin2"], ends["-fend2"], "out.dat", format, null, DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime);
            string[] eLines = [.. j.GetProperty("err").GetString()!.Split('\n').Where(l => l.Length > 0)];
            try
            {
                var got = EmbossFeaturePrograms.Featmerge(request);
                if (j.GetProperty("hang").GetBoolean())
                {
                    failures.Add($"{id}: the program never stops, the port returned");
                    continue;
                }

                string[] gLines = [.. got.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 0)
                    failures.Add($"{id}: the program died, the port returned");
                else if (got.Text != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(got.Text, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
                else if (expected.Length > 0)
                    written++;
            }
            catch (EmbossFeatureDiedException e)
            {
                string[] gLines = [.. e.Errors.SelectMany(m => m.Split('\n')).Append(e.Message).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 1)
                    failures.Add($"{id}: the port died, the program did not");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
            }
            catch (InvalidOperationException) when (j.GetProperty("hang").GetBoolean())
            {
                hang++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(400));
        Assert.That(written, Is.GreaterThanOrEqualTo(200), "tables written");
        Assert.That(hang, Is.GreaterThan(0));
    }

    /// <summary>
    /// A5-8c3f (F246): <c>feattext</c> (<see cref="EmbossFeaturePrograms.Feattext"/>) equals <c>feattext -features FILE -outfile out.txt</c> byte for byte, with its messages: the text of the entries the EMBL / GenBank and the GFF3 readers stored, nothing for the
    /// other formats (TestData/Alignment/feattext_oracle.jsonl.gz, <c>generate_featread_oracle.py --text</c>).
    /// </summary>
    [Test]
    public void Feattext_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.feattext_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, hang = 0, written = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string expected = j.GetProperty("out").GetString()!;
            string id = $"row {j.GetProperty("id")} ({j.GetProperty("kind").GetString()} {j.GetProperty("ufo").GetString()})";
            var request = new EmbossFeattextRequest(j.GetProperty("ufo").GetString()!, new Dictionary<string, string> { ["in.dat"] = j.GetProperty("text").GetString()! });
            string[] eLines = [.. j.GetProperty("err").GetString()!.Split('\n').Where(l => l.Length > 0)];
            try
            {
                var got = EmbossFeaturePrograms.Feattext(request);
                if (j.GetProperty("hang").GetBoolean())
                {
                    failures.Add($"{id}: the program never stops, the port returned");
                    continue;
                }

                string[] gLines = [.. got.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 0)
                    failures.Add($"{id}: the program died, the port returned");
                else if (got.Text != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(got.Text, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
                else if (expected.Length > 0)
                    written++;
            }
            catch (EmbossFeatureDiedException e)
            {
                string[] gLines = [.. e.Errors.SelectMany(m => m.Split('\n')).Append(e.Message).Where(l => l.Length > 0)];
                if (j.GetProperty("rc").GetInt32() != 1)
                    failures.Add($"{id}: the port died, the program did not");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: messages '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
            }
            catch (InvalidOperationException) when (j.GetProperty("hang").GetBoolean())
            {
                hang++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(400));
        Assert.That(written, Is.GreaterThanOrEqualTo(150), "entries written");
        Assert.That(hang, Is.GreaterThan(0));
    }

    private static string FirstDifference(string got, string expected)
    {
        int i = 0;
        while (i < got.Length && i < expected.Length && got[i] == expected[i])
            i++;
        string Cut(string s) => s[Math.Max(0, i - 40)..Math.Min(s.Length, i + 60)].Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: got '{Cut(got)}' expected '{Cut(expected)}'";
    }
}
