// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossFeatureTools, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Analysis;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossFeatureTools.Tools;

public partial class EmbossFeatureToolsTools
{
    public static string Qualifier(string name, object value) => "-" + name + " " + string.Format(System.Globalization.CultureInfo.InvariantCulture, "{0}", value);

    // the qualifiers typed after the files (a name with its value is one entry of the report header)
    public static IEnumerable<string> HeaderEntries(string? text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return [];
        var tokens = text.Split([' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
        var entries = new List<string>();
        int i = 0;
        while (i < tokens.Length)
        {
            string t = tokens[i];
            if (t.StartsWith('-') && i + 1 < tokens.Length && (!tokens[i + 1].StartsWith('-') || int.TryParse(tokens[i + 1], out _)))
            {
                entries.Add(t + " " + tokens[i + 1]);
                i += 2;
            }
            else
            {
                entries.Add(t);
                i++;
            }
        }

        return entries;
    }

    /// <summary>The report file text of a program run on one sequence (name + residues become a one-record FASTA file), in <paramref name="reportFormat"/> or with the report qualifiers of the text (<c>-rformat table -raccshow</c>).</summary>
    public static string EmbossReportText(string program, string sequenceName, string sequence, string? sequenceQualifiers, string reportFormat, IEnumerable<string> programEntries, Func<EmbossReportRequest, EmbossProgramReport> run, out string? original, CancellationToken cancellationToken = default)
    {
        string format = (reportFormat ?? "").Trim();
        bool qualifierText = format.StartsWith('-');
        var commandLine = new List<string> { EmbossCommandLineEcho.Parameter("sequence", "sequence.fasta") };
        commandLine.AddRange(programEntries);
        commandLine.AddRange(HeaderEntries(qualifierText ? format : "-rformat " + format));
        commandLine.AddRange(HeaderEntries(sequenceQualifiers));
        var request = new EmbossReportRequest(">" + sequenceName + "\n" + sequence + "\n", sequenceQualifiers, qualifierText ? format : "-rformat " + format, commandLine, DateTime.Now, "sequence.fasta", "stdout", CancellationToken: cancellationToken);
        try
        {
            var result = run(request);
            original = result.OriginalFile;
            return result.Report;
        }
        catch (EmbossReportException e)
        {
            throw new ArgumentException(e.Message, nameof(reportFormat));
        }
    }

    // the USA of the sequence input: the value given, else asis::TEXT, else the file the text is put in
    public static string SequenceSpec(string? sequenceValue, string sequence, string file)
    {
        if (!string.IsNullOrWhiteSpace(sequenceValue))
            return sequenceValue;
        return sequence.TrimStart().StartsWith("asis::", StringComparison.OrdinalIgnoreCase) ? sequence.Trim() : file;
    }

    // the report qualifiers of the text as the program numbers them in the header (-rformat -> -rformat3): the entries of the command line the report header echoes
    public static IEnumerable<string> NumberedReportEntries(string? text, int number)
    {
        foreach (string entry in HeaderEntries(text))
        {
            if (entry.StartsWith("-r", StringComparison.Ordinal) && entry.Length > 2 && char.IsAsciiLetter(entry[2]))
            {
                int blank = entry.IndexOf(' ', StringComparison.Ordinal);
                string name = blank < 0 ? entry : entry[..blank];
                yield return (char.IsAsciiDigit(name[^1]) ? name : name + number.ToString(System.Globalization.CultureInfo.InvariantCulture)) + (blank < 0 ? "" : entry[blank..]);
            }
            else
            {
                yield return entry;
            }
        }
    }

    // The EMBOSS reader's view of a nucleotide input (A4-3): gap characters, digits, blanks and every other non-letter are dropped first (EmbossSequenceReader.Clean); what is left must be IUPAC DNA.
    public static string RequireEmbossIupacDna(string sequence, string paramName) =>
        RequireIupacDna(sequence is null ? sequence! : EmbossSequenceReader.Clean(sequence), paramName);
}
