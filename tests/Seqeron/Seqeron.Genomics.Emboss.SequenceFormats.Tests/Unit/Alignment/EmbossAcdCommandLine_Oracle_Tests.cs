// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.SequenceFormats;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A6-2 (B09 F177): the resolution of the sequence-qualifier TEXT that follows the sequence arguments of an EMBOSS program (<c>needle a.fa b.fa -sbegin 3 -sre2 ...</c>) against the ACD processing log of the installed EMBOSS 6.6.0
/// (TestData/Alignment/emboss_acd_args_oracle.tsv, scripts/alignment/generate_emboss_acd_args_oracle.py: random texts for the 16 ported programs - numbered / unnumbered names, abbreviations (unique and ambiguous), <c>-</c> / <c>--</c> / <c>/</c> / no dash,
/// any case, <c>=value</c>, <c>-no</c> negations, boolean and integer values the program accepts or refuses). Every row is the program's own statement of what each token set; the port must refuse exactly the texts the program dies on (or that
/// are not sequence qualifiers) and otherwise set the same items to the same final values.
/// </summary>
[TestFixture]
public class EmbossAcdCommandLine_Oracle_Tests
{
    private sealed record Row(int Line, string Program, string Text, string Expected)
    {
        public override string ToString() => $"line {Line} {Program} '{Text}' -> {Expected}";
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_acd_args_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0)
                continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[1], f[2], f[3]));
        }

        return rows;
    });

    /// <summary>What the port makes of the text: <c>ERROR</c> (refused, whatever the reason) or the final values of the sequence qualifiers it set, in the form of the oracle.</summary>
    private static string Resolve(string program, string text, bool checkValues)
    {
        var profile = EmbossAcdProfiles.Get(program);
        var line = new AcdCommandLine(profile);
        try
        {
            line.Run(text.Split(' ', StringSplitOptions.RemoveEmptyEntries));
            if (checkValues)
                _ = EmbossSequenceQualifiers.ParseInputs(program, text);                // values, non-sequence qualifiers, ioffset
        }
        catch (ArgumentException e) when (!e.Message.StartsWith("-ioffset must be 0", StringComparison.Ordinal))        // the program accepts a non-zero offset here and fails reading the file: the port refuses it (a residue string has no file)
        {
            return "ERROR " + e.Message;
        }
        catch (ArgumentException)
        {
            // refused only for -ioffset: the resolution itself is compared
        }

        return string.Join(";", line.Assignments.Select(a => (a.Key.Name, a.Key.PNum, a.Value)).OrderBy(a => a.Name, StringComparer.Ordinal).ThenBy(a => a.PNum)
            .Select(a => $"{a.Name}[{a.PNum}]={a.Value}"));
    }

    private static void Check(string program, int minimum)
    {
        var rows = Oracle.Value.Where(r => r.Program == program).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum), $"oracle rows of {program}");
        var failures = new List<string>();
        foreach (var row in rows)
        {
            bool resolutionOnly = row.Expected.StartsWith('?');                          // the program died later, for a reason that hides the validity of the values
            string expected = row.Expected is "ERROR" or "REFUSED" ? "ERROR" : row.Expected.TrimStart('?');
            string got = Resolve(row.Program, row.Text, !resolutionOnly);
            if (got != expected && !(expected == "ERROR" && got.StartsWith("ERROR ", StringComparison.Ordinal)))
                failures.Add($"{row}\n   port {got}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows.Count(r => r.Expected is not ("ERROR" or "REFUSED") && !r.Expected.StartsWith('?')), Is.GreaterThan(rows.Count / 5), "accepted texts are a real share of the rows");
    }

    [Test] public void Needle_TextResolvesLikeAcd() => Check("needle", 100);
    [Test] public void Water_TextResolvesLikeAcd() => Check("water", 100);
    [Test] public void Stretcher_TextResolvesLikeAcd() => Check("stretcher", 100);
    [Test] public void Matcher_TextResolvesLikeAcd() => Check("matcher", 100);
    [Test] public void Supermatcher_TextResolvesLikeAcd() => Check("supermatcher", 100);
    [Test] public void Needleall_TextResolvesLikeAcd() => Check("needleall", 100);
    [Test] public void Merger_TextResolvesLikeAcd() => Check("merger", 100);
    [Test] public void Megamerger_TextResolvesLikeAcd() => Check("megamerger", 100);
    [Test] public void Seqmatchall_TextResolvesLikeAcd() => Check("seqmatchall", 100);
    [Test] public void Fuzznuc_TextResolvesLikeAcd() => Check("fuzznuc", 100);
    [Test] public void Fuzzpro_TextResolvesLikeAcd() => Check("fuzzpro", 100);
    [Test] public void Fuzztran_TextResolvesLikeAcd() => Check("fuzztran", 100);
    [Test] public void Etandem_TextResolvesLikeAcd() => Check("etandem", 100);
    [Test] public void Equicktandem_TextResolvesLikeAcd() => Check("equicktandem", 100);
    [Test] public void Dotmatcher_TextResolvesLikeAcd() => Check("dotmatcher", 100);
    [Test] public void Distmat_TextResolvesLikeAcd() => Check("distmat", 100);
}
