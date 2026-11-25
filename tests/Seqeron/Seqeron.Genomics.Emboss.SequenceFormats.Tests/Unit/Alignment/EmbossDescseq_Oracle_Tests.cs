// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.SequenceFormats;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3g (B09 F250): <c>descseq</c> of EMBOSS 6.6.0 (<c>emboss/descseq.c</c>: the name and the description of the first sequence replaced or appended to) equals the installed
/// <c>descseq -auto -sequence in.dat [-sformat1 F] [-name N] [-description D] [-append] -outseq o.fa -osformat2 fasta</c> byte for byte on 1 000 runs (TestData/Alignment/descseq_oracle.jsonl.gz,
/// scripts/alignment/generate_descseq_oracle.py): 19 input formats, names and descriptions with blanks, tabs and punctuation, empty values, damaged files.
/// </summary>
[TestFixture]
public class EmbossDescseq_Oracle_Tests
{
    [Test]
    public void Descseq_EqualsTheInstalledProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.descseq_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, masked = 0, refused = 0, died = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string text = j.GetProperty("text").GetString()!;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25;
            string fa = j.GetProperty("fa").GetString()!, err = j.GetProperty("err").GetString()!;
            string? sformat = null;
            string dname = "", ddesc = "";
            bool append = false, lower = false, upper = false;
            for (int i = 0; i < args.Count; i++)
            {
                switch (args[i])
                {
                    case "-sformat1":
                        sformat = args[i + 1];
                        break;
                    case "-name":
                        dname = args[i + 1];
                        break;
                    case "-description":
                        ddesc = args[i + 1];
                        break;
                    case "-append":
                        append = true;
                        break;
                    case "-slower1":
                        lower = true;
                        break;
                    case "-supper1":
                        upper = true;
                        break;
                }
            }

            string id = $"row {j.GetProperty("id")} [{j.GetProperty("fmt").GetString()} {string.Join(' ', args)}]";
            try
            {
                var o = EmbossDescseqPrograms.Descseq(new EmbossDescseqRequest(spec, new Dictionary<string, string> { ["in.dat"] = text }, new EmbossSequenceQualifiers { Format = sformat, Lower = lower, Upper = upper }, dname, ddesc, append));
                if (refusal)
                {
                    failures.Add($"{id}: the reference dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}; messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (o.SequenceText != fa)
                    failures.Add($"{id}: FASTA differs");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (fa.Length > 0 && (dname.Length > 0 || ddesc.Length > 0))
                    masked++;
                if (rc != 0)
                    died++;
            }
            catch (InvalidOperationException) when (refusal)
            {
                refused++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(masked, Is.GreaterThan(700), "runs that changed a name or a description");
        Assert.That(died, Is.GreaterThan(5), "the runs the program refuses");
        _ = refused;
    }
}
