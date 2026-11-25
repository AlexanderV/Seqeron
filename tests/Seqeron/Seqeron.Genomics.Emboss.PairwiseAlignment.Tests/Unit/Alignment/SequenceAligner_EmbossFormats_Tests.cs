// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.FeatureTools;
using Seqeron.Genomics.Emboss.SequenceFormats;
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-27, B09 F120): the EMBOSS 6.6.0 alignment output formats of <see cref="EmbossAlignmentWriter"/> (<c>-aformat3 unknown | fasta | a2m | msf | clustal |
/// mega | meganon | nexus | paup | nexusnon | paupnon | phylip | phylipnon | selex | treecon | markx0 | markx1 | markx2 | markx3 | markx10 | match | multiple | pair | simple | sam | score | srs | srspair | tcoffee</c>)
/// against the installed <c>needle</c> / <c>water</c> / <c>stretcher</c> 6.6.0: TestData/Alignment/emboss_aformat_oracle.tsv (scripts/alignment/generate_emboss_aformat_oracle.py; faked run dates), 110 alignment sets x 29 formats; every file byte for byte.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_EmbossFormats_Tests
{
    private sealed record Row(string Id, string Program, string Format, string[] Args, string FastaA, string FastaB, string Matrix, string Alignments, DateTime Date, string File);

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
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_aformat_oracle.tsv";
        using var stream = TestDataResource.Open(name)
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
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(' ', StringSplitOptions.RemoveEmptyEntries), Unescape(f[4]), Unescape(f[5]), f[6], f[7],
                DateTime.ParseExact(f[8], "yyyy-MM-dd'T'HH:mm:ss", CultureInfo.InvariantCulture), Unescape(f[9])));
        }

        return list;
    });

    // ">name some description" records of a FASTA text
    private static List<(string Name, string Desc)> Headers(string fasta) =>
        fasta.Split('\n').Where(l => l.StartsWith('>')).Select(l =>
        {
            string h = l[1..];
            int sp = h.IndexOf(' ', StringComparison.Ordinal);
            return sp < 0 ? (h, "") : (h[..sp], h[(sp + 1)..]);
        }).ToList();

    // the header prints a qualifier with its value on one line: "#    -gapopen 10"
    private static List<string> Groups(string[] args)
    {
        var groups = new List<string>();
        foreach (string t in args)
        {
            if (t.StartsWith('-') || groups.Count == 0)
                groups.Add(t);
            else
                groups[^1] += " " + t;
        }

        return groups;
    }

    private static string Produce(Row r)
    {
        string Opt(string name, string dflt)
        {
            int i = Array.IndexOf(r.Args, name);
            return i >= 0 ? r.Args[i + 1] : dflt;
        }

        Assert.That(EmbossAlignmentWriter.TryParseFormat(r.Format, out var format), r.Format);
        bool global = r.Program != "water";
        if (r.Args.Contains("-aglobal3"))
            global = true;
        if (r.Args.Contains("-noaglobal3"))
            global = false;
        int width = int.Parse(Opt("-awidth3", "0"), CultureInfo.InvariantCulture);
        var options = new EmbossAlignOptions
        {
            Program = r.Program,
            CommandLine = Groups(r.Args),
            RunDate = r.Date,
            Width = width,
            Global = global,
            ShowAccession = r.Args.Contains("-aaccshow3"),
            ShowDescription = r.Args.Contains("-adesshow3"),
            ShowUsa = r.Args.Contains("-ausashow3"),
            IntegerScores = r.Program == "stretcher",
            StatisticsSetByProgram = r.Program == "stretcher",
            Matrix = SubstitutionMatrix.Get(r.Matrix),
            MatrixName = r.Matrix,
            GapPenalty = double.Parse(Opt("-gapopen", "10"), CultureInfo.InvariantCulture),
            ExtendPenalty = double.Parse(Opt("-gapextend", "0.5"), CultureInfo.InvariantCulture),
        };
        var a = Headers(r.FastaA)[0];
        var bs = Headers(r.FastaB);
        var writer = new EmbossAlignmentWriter(format, options);
        int k = 0;
        foreach (string aln in r.Alignments.Split('|'))
        {
            var f = aln.Split(':');
            var b = bs[k++];
            var first = new EmbossAlignedSequence(a.Name, f[0], int.Parse(f[2], CultureInfo.InvariantCulture), int.Parse(f[3], CultureInfo.InvariantCulture), Description: a.Desc, Usa: "fasta::a.fa:" + a.Name);
            var second = new EmbossAlignedSequence(b.Name, f[1], int.Parse(f[4], CultureInfo.InvariantCulture), int.Parse(f[5], CultureInfo.InvariantCulture), Description: b.Desc, Usa: "fasta::b.fa:" + b.Name);
            writer.Write(new EmbossPairAlignment(first, second, double.Parse(f[6], CultureInfo.InvariantCulture)));
        }

        return writer.Close();
    }

    private static string FirstDifference(string got, string expected)
    {
        var g = got.Split('\n');
        var e = expected.Split('\n');
        for (int i = 0; i < Math.Max(g.Length, e.Length); i++)
        {
            string a = i < g.Length ? g[i] : "<end>";
            string b = i < e.Length ? e[i] : "<end>";
            if (a != b)
                return $"\n    line {i + 1}: got      '{a}'\n    line {i + 1}: expected '{b}'";
        }

        return "";
    }

    [Test]
    public void EveryFormat_IsByteIdenticalToEmbossOnAllCases()
    {
        var failures = new List<string>();
        foreach (var row in Oracle.Value)
        {
            string got;
            try
            {
                got = Produce(row);
            }
            catch (Exception e) when (e is not AssertionException)
            {
                failures.Add($"{row.Id} {row.Program} {row.Format}: {e.GetType().Name}: {e.Message}");
                continue;
            }

            if (got != row.File)
                failures.Add($"{row.Id} {row.Program} {row.Format} [{string.Join(' ', row.Args)}]{FirstDifference(got, row.File)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {Oracle.Value.Count}: " + string.Join("\n", failures.GroupBy(f => f[(f.IndexOf(' ', f.IndexOf(' ') + 1) + 1)..].Split(' ')[0]).Select(g => g.First()).Take(40)));
    }

    private sealed record SeqOutRow(string Id, string Format, char Type, string Fasta, DateTime Date, string File);

    private static readonly Lazy<List<SeqOutRow>> SeqOutOracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_seqout_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<SeqOutRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new SeqOutRow(f[0], f[1], f[2][0], Unescape(f[3]), DateTime.ParseExact(f[4], "yyyy-MM-dd'T'HH:mm:ss", CultureInfo.InvariantCulture), Unescape(f[5])));
        }

        return list;
    });

    [Test]
    public void SequenceSetFormats_AreByteIdenticalToEmbossSeqretOnSetsOfOneToNineRows()
    {
        var failures = new List<string>();
        foreach (var row in SeqOutOracle.Value)
        {
            var rows = new List<EmbossAlignedSequence>();
            foreach (string block in row.Fasta.Split('>', StringSplitOptions.RemoveEmptyEntries))
            {
                var lines = block.TrimEnd('\n').Split('\n');
                string h = lines[0];
                int sp = h.IndexOf(' ', StringComparison.Ordinal);
                rows.Add(new EmbossAlignedSequence(sp < 0 ? h : h[..sp], string.Concat(lines.Skip(1)), Description: sp < 0 ? "" : h[(sp + 1)..]));
            }

            Assert.That(EmbossAlignmentWriter.TryParseFormat(row.Format, out var format), row.Format);
            string got = EmbossAlignmentWriter.Format(format, rows, new EmbossAlignOptions { RunDate = row.Date, SeqoutType = row.Type });
            if (got != row.File)
                failures.Add($"{row.Id} {row.Format}{FirstDifference(got, row.File)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {SeqOutOracle.Value.Count}: " + string.Join("\n", failures.Take(20)));
        foreach (var fmt in new[] { "fasta", "msf", "clustal", "mega", "meganon", "nexus", "nexusnon", "phylip", "phylipnon", "selex", "treecon" })
            Assert.That(SeqOutOracle.Value.Count(r => r.Format == fmt), Is.GreaterThanOrEqualTo(100), fmt);
    }

    [Test]
    public void NeedleDefaultStyle_HeaderAndBlocks_MatchTheDocumentedLayout()
    {
        // the first example of the oracle written by hand: needle a.fa b.fa -aformat3 pair (EMBOSS 6.6.0)
        var first = new EmbossAlignedSequence("a", "ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT");
        var second = new EmbossAlignedSequence("b", "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT");
        string text = EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, [new EmbossPairAlignment(first, second, 146)], new EmbossAlignOptions
        {
            RunDate = new DateTime(2026, 11, 9, 21, 43, 53, DateTimeKind.Unspecified),
            CommandLine = ["-asequence a.fa", "-bsequence b.fa", "-aformat3 pair", "-auto", "-stdout"],
            Matrix = SubstitutionMatrix.EdnaFull,
            GapPenalty = 10,
            ExtendPenalty = 0.5,
            Global = true,
        });
        Assert.That(text, Does.Contain("# Rundate: Mon  9 Nov 2026 21:43:53\n"));
        Assert.That(text, Does.Contain("# Identity:      32/34 (94.1%)\n"));
        Assert.That(text, Does.Contain("# Score: 146.0\n"));
        Assert.That(text, Does.Contain("a                  1 ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT     33\n                     ||||||||||||||.|||||||| ||||||||||\nb                  1 ACGTACGTAGCTAGATAGCATCGAATCGATCGAT     34\n"));
        Assert.That(text, Does.EndWith("\n#---------------------------------------\n#---------------------------------------\n"));
    }

    [Test]
    public void SequenceStyleFormats_RejectPairOnlyFormats_AndPairFormatsRejectRows()
    {
        var rows = new[] { new EmbossAlignedSequence("a", "AC-T"), new EmbossAlignedSequence("b", "ACGT") };
        Assert.Throws<ArgumentException>(() => EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, rows));
        Assert.Throws<ArgumentException>(() => EmbossAlignmentWriter.Format(EmbossAlignFormat.Sam, rows));
        Assert.Throws<ArgumentException>(() => EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, [new EmbossPairAlignment(rows[0], new EmbossAlignedSequence("b", "ACG"), 1)]));
        Assert.Throws<ArgumentNullException>(() => EmbossAlignmentWriter.Format(EmbossAlignFormat.Pair, (IReadOnlyList<EmbossPairAlignment>)null!));
        var w = new EmbossAlignmentWriter(EmbossAlignFormat.Fasta);
        w.WriteRows(rows);
        w.Close();
        Assert.Throws<InvalidOperationException>(() => w.WriteRows(rows));
    }

    [Test]
    public void Oracle_CoversEveryFormatWithAtLeast100Cases()
    {
        var byFormat = Oracle.Value.GroupBy(r => r.Format).ToDictionary(g => g.Key, g => g.Count());
        foreach (var name in new[] { "unknown", "fasta", "a2m", "msf", "clustal", "mega", "meganon", "nexus", "paup", "nexusnon", "paupnon", "phylip", "phylipnon", "selex", "treecon", "markx0", "markx1", "markx2", "markx3", "markx10", "match", "multiple", "pair", "simple", "sam", "score", "srs", "srspair", "tcoffee" })
            Assert.That(byFormat.GetValueOrDefault(name), Is.GreaterThanOrEqualTo(100), name);
    }

    [Test]
    public void TryParseFormat_ParsesEveryEmbossName_CaseInsensitively()
    {
        foreach (EmbossAlignFormat f in Enum.GetValues<EmbossAlignFormat>())
        {
            Assert.That(EmbossAlignmentWriter.TryParseFormat(EmbossAlignmentWriter.FormatName(f), out var back), Is.True, f.ToString());
            Assert.That(back, Is.EqualTo(f));
            Assert.That(EmbossAlignmentWriter.TryParseFormat(EmbossAlignmentWriter.FormatName(f).ToUpperInvariant(), out back), Is.True);
        }

        Assert.That(EmbossAlignmentWriter.TryParseFormat("bam", out _), Is.False);
        Assert.That(EmbossAlignmentWriter.TryParseFormat(null, out _), Is.False);
    }
}
