#pragma warning disable S127

using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F263, WP58): <c>twofeat</c> of EMBOSS 6.6.0 (<c>emboss/twofeat.c</c>) through the installed program: <c>twofeat -auto -sequence in.dat -outfile o.rep [-sformat1 F] [-sbegin1 N] [-send1 N] [-asource ..] [-atype ..] [-asense 0|+|-]
/// [-aminscore X] [-amaxscore X] [-atag ..] [-avalue ..] [-bsource ..] ... [-overlap A|O|NO|NW|AW|BW] [-minrange N] [-maxrange N] [-rangetype N|L|R|F] [-sense A|S|O] [-order ..] [-twoout] [-typeout T] [-rformat2 FORMAT] [-raccshow2] ... [-rmaxall2 N] [-rmaxseq2 N]</c>.
/// TestData/Alignment/twofeat_oracle.jsonl.gz, scripts/alignment/generate_twofeat_oracle.py: EMBL entries with CDS / mRNA / exon / gene / misc features (joins, complements, orders, tags, several entries per file), the entries of the reader oracle in 7 formats,
/// damaged files; every option and report format. The report file and the messages (also the fatal error of the opposite sense) are compared byte for byte; the runs the reference cannot run are refused by the port.
/// </summary>
[TestFixture]
public class EmbossTwofeat_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 40)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 40)] : "<end>").Replace("\n", "\\n")}')";
    }

    [Test]
    public void Twofeat_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.twofeat_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withReport = 0, defaultNames = 0, directories = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            var a = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            var cmd = j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string expected = j.GetProperty("out").GetString()!, err = j.GetProperty("err").GetString()!;
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {string.Join(' ', a)}]";
            try
            {
                var q = new List<string>();
                string asource = "*", atype = "*", asense = "0", atag = "*", avalue = "*", bsource = "*", btype = "*", bsense = "0", btag = "*", bvalue = "*", overlap = "A", rangetype = "N", sense = "A", order = "A", typeout = "misc_feature";
                float amin = 0, amax = 0, bmin = 0, bmax = 0;
                int minrange = 0, maxrange = 0;
                bool twoout = false;
                string? outfile = null, rname = null, rext = null, rdir = null;
                var inv = CultureInfo.InvariantCulture;
                for (int i = 0; i < a.Count; i++)
                {
                    switch (a[i])
                    {
                        case "-auto": break;
                        case "-sequence": i++; break;
                        case "-outfile": outfile = a[++i]; break;
                        case "-rname2": rname = a[++i]; break;
                        case "-rextension2": rext = a[++i]; break;
                        case "-rdirectory2": rdir = a[++i]; break;
                        case "-asource": asource = a[++i]; break;
                        case "-atype": atype = a[++i]; break;
                        case "-asense": asense = a[++i]; break;
                        case "-aminscore": amin = float.Parse(a[++i], inv); break;
                        case "-amaxscore": amax = float.Parse(a[++i], inv); break;
                        case "-atag": atag = a[++i]; break;
                        case "-avalue": avalue = a[++i]; break;
                        case "-bsource": bsource = a[++i]; break;
                        case "-btype": btype = a[++i]; break;
                        case "-bsense": bsense = a[++i]; break;
                        case "-bminscore": bmin = float.Parse(a[++i], inv); break;
                        case "-bmaxscore": bmax = float.Parse(a[++i], inv); break;
                        case "-btag": btag = a[++i]; break;
                        case "-bvalue": bvalue = a[++i]; break;
                        case "-overlap": overlap = a[++i]; break;
                        case "-minrange": minrange = int.Parse(a[++i], inv); break;
                        case "-maxrange": maxrange = int.Parse(a[++i], inv); break;
                        case "-rangetype": rangetype = a[++i]; break;
                        case "-sense": sense = a[++i]; break;
                        case "-order": order = a[++i]; break;
                        case "-twoout": twoout = true; break;
                        case "-typeout": typeout = a[++i]; break;
                        default:
                            q.Add(a[i]);
                            if (i + 1 < a.Count && !a[i + 1].StartsWith('-'))
                                q.Add(a[++i]);
                            break;
                    }
                }

                var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport("twofeat", string.Join(' ', q));
                rq = rq with { Name = rname, Extension = rext, Directory = rdir };
                var o = EmbossTwofeatProgram.Twofeat(new EmbossTwofeatRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, inputs[0], rq, cmd, date, outfile ?? "o.rep")
                {
                    DefaultReportFile = outfile is null, ExistingDirectories = ["sub"],
                    ASource = asource, AType = atype, ASense = asense, AMinScore = amin, AMaxScore = amax, ATag = atag, AValue = avalue, BSource = bsource, BType = btype, BSense = bsense, BMinScore = bmin, BMaxScore = bmax, BTag = btag, BValue = bvalue,
                    Overlap = overlap, MinRange = minrange, MaxRange = maxrange, RangeType = rangetype, Sense = sense, Order = order, TwoOut = twoout, TypeOut = typeout,
                });
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                var got = o.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0).ToList();
                var want = err.Split('\n').Where(l => l.Length > 0).ToList();
                string? why = null;
                if ((o.Died ? 1 : 0) != rc)
                    why = $"exit code {(o.Died ? 1 : 0)} expected {rc}";
                else if (!o.Died && o.Report != expected)
                    why = $"report differs at {FirstDifference(o.Report, expected)}";
                else if (!got.SequenceEqual(want))
                    why = $"messages '{string.Join("|", got)}' expected '{string.Join("|", want)}'";
                else if (!o.Died && j.TryGetProperty("outfile", out var of) && of.GetString() is { Length: > 0 } wantFile && o.ReportFile != wantFile)
                    why = $"report file '{o.ReportFile}' expected '{wantFile}'";
                if (why is not null)
                {
                    failures.Add($"{id}: {why}");
                    continue;
                }

                compared++;
                if (!a.Contains("-outfile"))
                    defaultNames++;
                if (a.Contains("-rdirectory2") && expected.Length > 0)
                    directories++;
                if (expected.Length > 0)
                {
                    withReport++;
                    formats.Add(a.Contains("-rformat2") ? a[a.IndexOf("-rformat2") + 1] : "table");
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
        Assert.That(compared, Is.GreaterThanOrEqualTo(955), $"compared {compared}, refused {refused}");
        Assert.That(withReport, Is.GreaterThanOrEqualTo(600));
        Assert.That(defaultNames, Is.GreaterThanOrEqualTo(100), "the runs with the default report file name (-rname2 / -rextension2)");
        Assert.That(directories, Is.GreaterThanOrEqualTo(30), "the runs with -rdirectory2");
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(20));
    }
}
