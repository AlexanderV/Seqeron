using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A6-2 (B09 F177), effect level: the 16 ported EMBOSS programs run on sequence qualifiers written as command-line TEXT after the sequence arguments (<c>needle a.fa b.fa -sbegin 3 -sre2 -nosupper1=...</c>: numbered, unnumbered,
/// abbreviated, <c>=</c> forms, any case; TestData/Alignment/emboss_seqopts_text_oracle.tsv, scripts/alignment/generate_emboss_seqopts_text_oracle.py: 40 rows per program, the program runners of the A5-7 oracle) against what the program
/// REPORTS. The port resolves the text with <see cref="EmbossSequenceQualifiers.ParseInputs"/> and evaluates the row through the same reader / port code as the structured rows (<see cref="EmbossSeqOptsProducer"/>); a text the program refuses must be refused.
/// </summary>
[TestFixture]
public class EmbossSequenceQualifiersText_Oracle_Tests
{
    private sealed record TextRow(int Line, string Program, string Kind, string Params, string[] Seqs, string Text, string Expected);

    private static readonly Lazy<List<TextRow>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_seqopts_text_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<TextRow>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0)
                continue;
            var f = text.Split('\t');
            rows.Add(new TextRow(line, f[1], f[2], f[3], [f[4], f[5], f[6]], f[7], f[8]));
        }

        return rows;
    });

    private static string Format(EmbossSequenceQualifiers q)
        => string.Join('|', new[] { q.Begin, q.End, q.Reverse ? 1 : 0, q.Lower ? 1 : 0, q.Upper ? 1 : 0, q.Circular ? 1 : 0, q.Nucleotide ? 1 : 0, q.Protein ? 1 : 0 });

    private static void Check(string program, int minimum)
    {
        var rows = Oracle.Value.Where(r => r.Program == program).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(minimum), $"oracle rows of {program}");
        var failures = new List<string>();
        int errors = 0, unnumbered = 0;
        foreach (var r in rows)
        {
            string got;
            var shape = new EmbossSeqOptsProducer.Row(r.Line, r.Program, r.Kind, r.Params, r.Seqs, "", "", r.Expected);
            try
            {
                var q = EmbossSequenceQualifiers.ParseInputs(program, r.Text);
                shape = shape with { Quals1 = Format(q[0]), Quals2 = Format(q.Count > 1 ? q[1] : EmbossSequenceQualifiers.None) };
                got = EmbossSeqOptsProducer.Produce(shape);
            }
            catch (ArgumentException)
            {
                got = "ERROR";
            }

            string expected = EmbossSeqOptsProducer.Expected(shape);
            if (got != expected)
                failures.Add($"{program} line {r.Line} '{r.Text}'\n   expected {expected}\n   port     {got}");
            if (expected == "ERROR")
                errors++;
            if (System.Text.RegularExpressions.Regex.IsMatch(r.Text, @"(^|\s)[-/]{0,2}s[a-z]+(\s|$)", System.Text.RegularExpressions.RegexOptions.IgnoreCase))
                unnumbered++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(4)));
        Assert.That(errors, Is.LessThan(rows.Count * 3 / 4), "the refused texts stay a minority");
        Assert.That(unnumbered, Is.GreaterThan(rows.Count / 5), "rows with an unnumbered qualifier");
    }

    [Test] public void Needle_TextLikeEmboss() => Check("needle", 40);
    [Test] public void Water_TextLikeEmboss() => Check("water", 40);
    [Test] public void Stretcher_TextLikeEmboss() => Check("stretcher", 40);
    [Test] public void Matcher_TextLikeEmboss() => Check("matcher", 40);
    [Test] public void Supermatcher_TextLikeEmboss() => Check("supermatcher", 40);
    [Test] public void Needleall_TextLikeEmboss() => Check("needleall", 40);
    [Test] public void Merger_TextLikeEmboss() => Check("merger", 40);
    [Test] public void Megamerger_TextLikeEmboss() => Check("megamerger", 40);
    [Test] public void Seqmatchall_TextLikeEmboss() => Check("seqmatchall", 40);
    [Test] public void Fuzznuc_TextLikeEmboss() => Check("fuzznuc", 40);
    [Test] public void Fuzzpro_TextLikeEmboss() => Check("fuzzpro", 40);
    [Test] public void Fuzztran_TextLikeEmboss() => Check("fuzztran", 40);
    [Test] public void Etandem_TextLikeEmboss() => Check("etandem", 40);
    [Test] public void Equicktandem_TextLikeEmboss() => Check("equicktandem", 40);
    [Test] public void Dotmatcher_TextLikeEmboss() => Check("dotmatcher", 40);
    [Test] public void Distmat_TextLikeEmboss() => Check("distmat", 40);
}
