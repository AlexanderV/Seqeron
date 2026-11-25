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
/// B09 WP67 (F309): the USAs <c>db:id</c> of databases that only the data resource catalogue (DRCAT.dat of the EMBOSS data) defines, through seqret and entret of EMBOSS 6.6.0 (<c>queryDbfind</c> -> <c>ajResourceGetDbdata</c>): the database is found
/// in the catalogue when one of its query lines has a format the sequence reader knows (a name or the EDAM number of the format) and the entry is then fetched from the URL of that line. TestData/Alignment/emboss_drcat_usa_oracle.jsonl.gz
/// (scripts/alignment/generate_emboss_drcat_usa_oracle.py): the installed programs on a machine without network access, where a fetch fails ("Unable to read sequence"); with network access some fetches do not return within 60 s
/// (the user's run, handoff/WP62-drcat-lookup/returned), which is not a result the port can have: the port has no network and reports the failed fetch.
/// </summary>
[TestFixture]
public class EmbossDrcat_Oracle_Tests
{
    [Test]
    public void Usas_EqualTheInstalledPrograms_AllRows()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_drcat_usa_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var kinds = new HashSet<string>();
        int rows = 0, found = 0, network = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string spec = args[args.IndexOf("-sequence") + 1];
            string program = j.GetProperty("prog").GetString()!;
            string err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{program} {string.Join(' ', args)}]";
            if (err.StartsWith("Warning: ", StringComparison.Ordinal))
            {
                network++;                                  // a public server (mrs, entrez, das, ensembl) was asked: its answer (database list, HTTP code) is the server's, not the library's
                continue;
            }

            try
            {
                IReadOnlyList<string> messages;
                bool died;
                string text;
                var q = EmbossSequenceQualifiers.ParseFile("");
                if (program == "entret")
                {
                    var e = EmbossEntretProgram.Entret(new EmbossEntretRequest(spec, files, q, false));
                    (messages, died, text) = (e.Messages, e.Died, e.Text);
                }
                else
                {
                    var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, files, q, false, false, "o.fa", "fasta", "o.gff", "gff3", DateTime.UnixEpoch));
                    (messages, died, text) = (o.Messages, o.Died, o.SequenceText);
                }

                string got = messages.Count > 0 ? string.Join("\n", messages) + "\n" : "";
                if ((died ? 1 : 0) != j.GetProperty("rc").GetInt32() || got != err || text != j.GetProperty("out").GetString())
                    failures.Add($"{id}: got '{got.Replace("\n", "|")}' rc {(died ? 1 : 0)} expected '{err.Replace("\n", "|")}' rc {j.GetProperty("rc")}");
                else
                {
                    kinds.Add(err.Split('\n')[0]);
                    if (err.StartsWith("Error: Unable to read sequence", StringComparison.Ordinal))
                        found++;
                }
            }
            catch (Exception e) when (e is InvalidOperationException or NotSupportedException)
            {
                failures.Add($"{id}: refused ({e.Message[..Math.Min(90, e.Message.Length)]}) but the reference ran");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(15)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(800));
        Assert.That(network, Is.LessThan(12), "rows that asked a public server (not compared)");
        Assert.That(found, Is.GreaterThan(20), "rows whose database the catalogue defines (the fetch fails)");
        Assert.That(kinds.Count, Is.GreaterThan(5), "kinds of first messages");
    }

    [Test]
    public void Catalogue_HasTheEntriesOfTheInstalledFile()
    {
        // grep -c '^ID ' DRCAT.dat = 675 entries, 41 IDalt lines; swissprot is an IDalt of UniProtKB_Swiss-Prot, the index keys are lower case
        Assert.Multiple(() =>
        {
            Assert.That(EmbossDrcat.All, Has.Count.EqualTo(675));
            Assert.That(EmbossDrcat.All.Sum(r => r.Idalt.Count), Is.EqualTo(41));
            Assert.That(EmbossDrcat.Lookup("SwissProt")!.Id, Is.EqualTo("UniProtKB_Swiss-Prot"));
            Assert.That(EmbossDrcat.Lookup("embl")!.Query, Has.Count.EqualTo(2));
            Assert.That(EmbossDrcat.Lookup("EMBL")!.Query[1].FormatTerm, Is.EqualTo("1927"));
            Assert.That(EmbossDrcat.Lookup("nosuchdb"), Is.Null);
        });
    }

    /// <summary>The user's run of <c>showdb -only -name</c> (handoff/WP62-drcat-lookup/returned, row 015): "Died: Unknown qualifier -name", status 1. The port refuses a qualifier the library does not know with an <see cref="ArgumentException"/> of the same text.</summary>
    [Test]
    public void Showdb_UnknownQualifier_IsRefusedWithTheMessageOfTheLibrary()
    {
        var e = Assert.Throws<ArgumentException>(() => EmbossProgramRunner.Run("showdb", ["-only", "-name"], new Dictionary<string, string>(), new EmbossRunEnvironment(null, null, null)));
        Assert.That(e!.Message, Does.Contain("Unknown qualifier -name"));
    }

    /// <summary>
    /// The user's run with network access (handoff/WP62-drcat-lookup/returned, rows 001-013: seqret / entret of USAs of the catalogue databases embl, swissprot, uniprot, genbank, refseq, pdb ...): the 6 fetches of seqret that end by themselves fail as they fail on the machine without network
    /// and equal the port; 2 runs of entret used -outseq, which it does not have; the 5 runs that did not end within 60 s (the server accepted the connection and did not answer) are not a result the port can have: it has no network and reports the failed fetch.
    /// </summary>
    [Test]
    public void ReturnedRuns_WithNetwork_EqualThePortOrAreRunsThatNeverEnded()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_drcat_returned_rows.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded data '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int equal = 0, neverEnded = 0, unknownQualifier = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string spec = args[args.IndexOf("-sequence") + 1];
            var q = EmbossSequenceQualifiers.ParseFile("");
            IReadOnlyList<string> messages;
            bool died;
            if (j.GetProperty("err").GetString()!.StartsWith("Died: Unknown qualifier -outseq", StringComparison.Ordinal))
            {
                // entret has no -outseq (a slip of the user's rows.txt): the answer is the command line of the program, which entret (an API, not a program of the ACD runner) does not have
                unknownQualifier++;
                continue;
            }

            if (j.GetProperty("prog").GetString() == "entret")
            {
                var e = EmbossEntretProgram.Entret(new EmbossEntretRequest(spec, new Dictionary<string, string>(), q, false));
                (messages, died) = (e.Messages, e.Died);
            }
            else
            {
                var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, new Dictionary<string, string>(), q, false, false, "o.fa", "fasta", "o.gff", "gff3", DateTime.UnixEpoch));
                (messages, died) = (o.Messages, o.Died);
            }

            string got = messages.Count > 0 ? string.Join("\n", messages) + "\n" : "";
            if (j.GetProperty("timeout").GetBoolean())
            {
                neverEnded++;
                if (!died || !got.Contains("Unable to read sequence", StringComparison.Ordinal))
                    failures.Add($"{spec}: the port should report the failed fetch, got '{got}'");
            }
            else if (got != j.GetProperty("err").GetString() || (died ? 1 : 0) != j.GetProperty("rc").GetInt32())
            {
                failures.Add($"{spec}: got '{got.Replace("\n", "|")}' rc {(died ? 1 : 0)} expected '{j.GetProperty("err").GetString()!.Replace("\n", "|")}' rc {j.GetProperty("rc")}");
            }
            else
            {
                equal++;
            }
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures));
        Assert.That(equal, Is.EqualTo(6), "runs of seqret whose fetch fails");
        Assert.That(unknownQualifier, Is.EqualTo(2), "runs of entret with -outseq, which it does not have");
        Assert.That(neverEnded, Is.EqualTo(5));
    }
}
