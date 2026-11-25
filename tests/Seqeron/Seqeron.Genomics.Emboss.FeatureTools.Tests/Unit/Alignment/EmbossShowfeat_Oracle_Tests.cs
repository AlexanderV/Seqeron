// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S127
using Seqeron.Genomics.Emboss.FeatureTools;
using Seqeron.Genomics.Emboss.SequenceFormats;
using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F260, WP58): <c>showfeat</c> of EMBOSS 6.6.0 (<c>emboss/showfeat.c</c>) through the installed program: <c>showfeat -auto -sequence in.dat [-sformat1 F] -outfile o.sf [-sourcematch ..] [-typematch ..] [-tagmatch ..] [-valuematch ..]
/// [-sort ..] [-joinfeatures] [-html] [-noid] [-nodescription] [-noscale] [-width N] [-collapse] [-noforward] [-noreverse] [-nounknown] [-strand] [-sourceshow] [-position] [-notypeshow] [-tagshow] [-novalueshow] [-stricttags]
/// [-annotation RANGES] [-sbegin1 N] [-send1 N] [-sreverse1]</c>. TestData/Alignment/showfeat_oracle.jsonl.gz, scripts/alignment/generate_showfeat_oracle.py: entries with feature tables in embl, genbank, swiss, nbrf, gff3, refseqp, gff2,
/// EMBL entries with CDS / mRNA / exon features (joins, complements, orders, partial ends), FASTA files without features, damaged files; every option. The output file and the messages of the library are compared byte for byte;
/// the runs the reference cannot run (a fault, a loop, memory errors) are refused by the port.
/// </summary>
[TestFixture]
public class EmbossShowfeat_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 30)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 30)] : "<end>").Replace("\n", "\\n")}')";
    }

    private static EmbossShowfeatRequest Request(string text, List<string> a)
    {
        var q = new List<string>();
        string sm = "*", tm = "*", tgm = "*", vm = "*", sort = "start", ann = "";
        bool jf = false, html = false, idf = true, desc = true, scale = true, coll = false, fw = true, rv = true, unk = true, strand = false, src = false, pos = false, ty = true, tg = false, vl = true, st = false;
        int width = 60;
        for (int i = 0; i < a.Count; i++)
        {
            switch (a[i])
            {
                case "-sourcematch": sm = a[++i]; break;
                case "-typematch": tm = a[++i]; break;
                case "-tagmatch": tgm = a[++i]; break;
                case "-valuematch": vm = a[++i]; break;
                case "-sort": sort = a[++i]; break;
                case "-joinfeatures": jf = true; break;
                case "-html": html = true; break;
                case "-noid": idf = false; break;
                case "-nodescription": desc = false; break;
                case "-noscale": scale = false; break;
                case "-width": width = int.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-collapse": coll = true; break;
                case "-noforward": fw = false; break;
                case "-noreverse": rv = false; break;
                case "-nounknown": unk = false; break;
                case "-strand": strand = true; break;
                case "-sourceshow": src = true; break;
                case "-position": pos = true; break;
                case "-notypeshow": ty = false; break;
                case "-tagshow": tg = true; break;
                case "-novalueshow": vl = false; break;
                case "-stricttags": st = true; break;
                case "-annotation": ann = a[++i]; break;
                default:
                    q.Add(a[i]);
                    if (a[i] is "-sformat1" or "-sbegin1" or "-send1")
                        q.Add(a[++i]);
                    break;
            }
        }

        var quals = EmbossSequenceQualifiers.ParseFile(string.Join(' ', q));
        return new EmbossShowfeatRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, quals)
        {
            SourceMatch = sm, TypeMatch = tm, TagMatch = tgm, ValueMatch = vm, Sort = sort, JoinFeatures = jf, Html = html, Id = idf, Description = desc, Scale = scale, Width = width, Collapse = coll, Forward = fw, Reverse = rv,
            Unknown = unk, Strand = strand, SourceShow = src, Position = pos, TypeShow = ty, TagShow = tg, ValueShow = vl, StrictTags = st, Annotation = ann,
        };
    }

    [Test]
    public void Extractfeat_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.showfeat_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withOutput = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string expected = j.GetProperty("out").GetString()!, err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {spec} {string.Join(' ', args)}]";
            try
            {
                var o = EmbossShowfeatProgram.Showfeat(Request(text, args));
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}");
                else if ((o.Died ? "" : o.Text) != expected)
                    failures.Add($"{id}: text differs at {FirstDifference(o.Text, expected)}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else
                {
                    compared++;
                    if (expected.Length > 0)
                    {
                        withOutput++;
                        formats.Add(fmt);
                    }
                }
            }
            catch (Exception e) when (e is InvalidOperationException or NotSupportedException)
            {
                if (memerr || refusal)
                    refused++;
                else
                    failures.Add($"{id}: refused ({e.Message[..Math.Min(70, e.Message.Length)]}) but the reference ran (exit code {rc})");
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(compared, Is.GreaterThanOrEqualTo(900), $"compared {compared}, refused {refused}");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(880));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(7));
    }
}
