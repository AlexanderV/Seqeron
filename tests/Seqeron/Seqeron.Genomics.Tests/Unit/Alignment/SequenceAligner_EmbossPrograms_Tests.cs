using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 / ALIGN-GLOBAL-001 (audit A2-28c, B09 F124): EMBOSS 6.6.0 <c>supermatcher</c> (Rabin-Karp word seeds + banded Smith-Waterman) and <c>needleall</c> (all-vs-all needle) against the installed programs:
/// TestData/Alignment/emboss_multi_oracle.tsv (scripts/alignment/generate_emboss_multi_oracle.py), whole output files (header, every -aformat) and the error files, byte for byte.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-LOCAL-001")]
public class SequenceAligner_EmbossPrograms_Tests
{
    private sealed record Row(string Id, string Program, string[] Args, string FastaA, string FastaB, string Format, DateTime Date, string File, string ErrFile);

    private static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        bool escaped = false;
        foreach (char c in s)
        {
            if (escaped)
            {
                sb.Append(c switch { 'n' => '\n', 't' => '\t', _ => c });
                escaped = false;
            }
            else if (c == '\\')
            {
                escaped = true;
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_multi_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2].Split(' ', StringSplitOptions.RemoveEmptyEntries), Unescape(f[3]), Unescape(f[4]), f[5],
                DateTime.ParseExact(f[6], "yyyy-MM-dd'T'HH:mm:ss", CultureInfo.InvariantCulture), Unescape(f[7]), Unescape(f[8])));
        }

        return list;
    });

    private static List<EmbossSequence> ParseFasta(string fasta)
    {
        var list = new List<EmbossSequence>();
        foreach (string block in fasta.Split('>', StringSplitOptions.RemoveEmptyEntries))
        {
            var lines = block.TrimEnd('\n').Split('\n');
            string h = lines[0];
            int sp = h.IndexOf(' ', StringComparison.Ordinal);
            list.Add(new EmbossSequence(sp < 0 ? h : h[..sp], string.Concat(lines.Skip(1)), sp < 0 ? "" : h[(sp + 1)..]));
        }

        return list;
    }

    private static List<string> Groups(string[] args)
    {
        var groups = new List<string>();
        foreach (string t in args)
        {
            if ((t.StartsWith('-') && !double.TryParse(t, NumberStyles.Float, CultureInfo.InvariantCulture, out _)) || groups.Count == 0)
                groups.Add(t);
            else
                groups[^1] += " " + t;
        }

        return groups;
    }

    private static (string Text, string Err) Produce(Row r)
    {
        string? Opt(string name)
        {
            int i = Array.IndexOf(r.Args, name);
            return i >= 0 ? r.Args[i + 1] : null;
        }

        double D(string name, double dflt) => Opt(name) is { } v ? double.Parse(v, CultureInfo.InvariantCulture) : dflt;
        var a = ParseFasta(r.FastaA);
        var b = ParseFasta(r.FastaB);
        string? matrixName = Opt("-datafile");
        bool nuc = a.All(s => s.Residues.ToUpperInvariant().All(c => "ACGTUBDHKMNRSVWXY".Contains(c)))
            && (r.Program == "needleall" || b.All(s => s.Residues.ToUpperInvariant().All(c => "ACGTUBDHKMNRSVWXY".Contains(c))));
        matrixName ??= nuc ? "EDNAFULL" : "EBLOSUM62";
        var matrix = SubstitutionMatrix.Get(matrixName);
        Assert.That(EmbossAlignmentWriter.TryParseFormat(r.Format, out var format), r.Format);
        EmbossProgramResult result;
        if (r.Program == "supermatcher")
        {
            result = SequenceAligner.Supermatcher(a, b, matrix, D("-gapopen", 10), D("-gapextend", 0.5), (int)D("-wordlen", 6), (int)D("-width", 16), D("-minscore", 0), format);
        }
        else
        {
            result = SequenceAligner.NeedleAll(a, b, matrix, D("-gapopen", 10), D("-gapextend", 0.5), r.Args.Contains("-endweight"), D("-endopen", 10), D("-endextend", 0.5), D("-minscore", 0), !r.Args.Contains("-nobrief"));
        }

        // the formats that do not show sequences are written from plain rows: no matrix, no gap penalties in the header
        bool plain = format is EmbossAlignFormat.Match or EmbossAlignFormat.Score;
        var options = new EmbossAlignOptions
        {
            Program = r.Program,
            CommandLine = Groups(r.Args),
            RunDate = r.Date,
            Global = r.Program == "needleall",
            Matrix = plain ? null : matrix,
            MatrixName = plain ? "" : matrixName,
            GapPenalty = plain ? null : D("-gapopen", 10),
            ExtendPenalty = plain ? null : D("-gapextend", 0.5),
        };
        var writer = new EmbossAlignmentWriter(format, options);
        foreach (var al in result.Alignments)
            writer.Write(al);
        return (writer.Close(), result.ErrorText);
    }

    private static string FirstDifference(string got, string expected)
    {
        var g = got.Split('\n');
        var e = expected.Split('\n');
        for (int i = 0; i < Math.Max(g.Length, e.Length); i++)
        {
            string x = i < g.Length ? g[i] : "<end>";
            string y = i < e.Length ? e[i] : "<end>";
            if (x != y)
                return $"\n    line {i + 1}: got      '{x}'\n    line {i + 1}: expected '{y}'";
        }

        return "";
    }

    private static void Check(string program, int minimum)
    {
        var failures = new List<string>();
        var rows = Oracle.Value.Where(r => r.Program == program).ToList();
        foreach (var row in rows)
        {
            try
            {
                var (text, err) = Produce(row);
                if (text != row.File)
                    failures.Add($"{row.Id} {program} {row.Format} [{string.Join(' ', row.Args)}]{FirstDifference(text, row.File)}");
                else if (err != row.ErrFile)
                    failures.Add($"{row.Id} {program} errfile{FirstDifference(err, row.ErrFile)}");
            }
            catch (Exception e) when (e is not AssertionException)
            {
                failures.Add($"{row.Id} {program} {row.Format}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}: " + string.Join("\n", failures.Take(400)));
    }

    [Test]
    public void Supermatcher_EveryOutputFileAndErrorFile_IsByteIdenticalToEmbossOnAllCases() => Check("supermatcher", 150);

    [Test]
    public void NeedleAll_EveryOutputFileAndErrorFile_IsByteIdenticalToEmbossOnAllCases() => Check("needleall", 100);

    [Test]
    public void Supermatcher_NoSharedWord_WritesTheErrorLineAndNoAlignment()
    {
        var q = new[] { new EmbossSequence("q", "AAAAAAAAAAAAAAAAAAAA") };
        var t = new[] { new EmbossSequence("t", "CCCCCCCCCCCCCCCCCCCC") };
        var r = SequenceAligner.Supermatcher(q, t);
        Assert.That(r.Alignments, Is.Empty);
        Assert.That(r.ErrorText, Is.EqualTo("No wordmatch start points for q vs t. No alignment\n"));
    }
}
