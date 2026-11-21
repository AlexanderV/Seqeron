using System.Globalization;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Analysis;

/// <summary>One pattern of a pattern list as <c>ajPatternSeqNewList</c> stores it (the report's <c>Pattern_name Mismatch Pattern</c> row).</summary>
/// <param name="Name">The pattern name.</param>
/// <param name="Pattern">The pattern text, as given.</param>
/// <param name="Mismatch">The number of mismatches (<c>-pmismatch</c> or the <c>&lt;mismatch=n&gt;</c> of the pattern file header).</param>
public readonly record struct EmbossProgramPattern(string Name, string Pattern, int Mismatch);

/// <summary>The hits of one pattern of a pattern list, as <see cref="EmbossReports"/> writes them.</summary>
/// <typeparam name="THit">The hit type of the program.</typeparam>
/// <param name="Name">The pattern name.</param>
/// <param name="Pattern">The pattern text.</param>
/// <param name="Mismatches">The mismatches the search allowed.</param>
/// <param name="Hits">The hits of the pattern.</param>
public readonly record struct EmbossPatternHits<THit>(string Name, string Pattern, int Mismatches, IReadOnlyList<THit> Hits);

/// <summary>
/// The pattern list of the ACD type <c>pattern</c> (EMBOSS 6.6.0 <c>ajax/core/ajpat.c</c> <c>ajPatlistSeqRead</c>, read by fuzznuc, fuzzpro and fuzztran): a single pattern, or <c>@file</c> with one pattern per line
/// (<c>-pformat simple</c>) or FASTA style (a <c>&gt;name</c> header, optionally with <c>&lt;mismatch=n&gt;</c>, then the pattern over one or several lines). A5-8c, B09 F192.
/// </summary>
public static class EmbossPatternList
{
    private static readonly Regex MismatchTag = new(@"<mismatch=(\d+)>", RegexOptions.Compiled | RegexOptions.CultureInvariant);

    /// <summary>
    /// Reads a pattern specification. Measured: <c>-pformat fasta</c> is the program's "format 0" like no format (the first line decides: a line starting with <c>&gt;</c> means FASTA style); <c>simple</c> forces one pattern per line.
    /// An empty pattern is dropped, a pattern without a name is called <c>pattern&lt;n&gt;</c> (n = its position in the list).
    /// </summary>
    /// <param name="spec">The pattern, or <c>@file</c>.</param>
    /// <param name="patternName"><c>-pname</c>: the default name (null = none given: "pattern" for a single pattern and for the lines of a simple file, the position-numbered name for unnamed FASTA headers).</param>
    /// <param name="format"><c>-pformat</c>: <c>simple</c>, <c>fasta</c> or null.</param>
    /// <param name="mismatches"><c>-pmismatch</c> (fuzznuc; 0 for the others).</param>
    /// <param name="readFile">Reads the text of a pattern file; required for <c>@file</c>.</param>
    /// <exception cref="ArgumentException">The file cannot be read.</exception>
    public static IReadOnlyList<EmbossProgramPattern> Read(string spec, string? patternName = null, string? format = null, int mismatches = 0, Func<string, string>? readFile = null)
    {
        ArgumentNullException.ThrowIfNull(spec);
        patternName ??= "";
        string namestr = patternName.Length > 0 ? patternName : "pattern";
        if (mismatches < 0)
            mismatches = 0;
        var list = new List<EmbossProgramPattern>();

        if (spec.Length == 0 || spec[0] != '@')
        {
            Add(list, namestr, spec, mismatches);
            return list;
        }

        string path = spec[1..];
        if (readFile is null)
            throw new ArgumentException($"Unable to open pattern file '{path}'", nameof(spec));
        string text = readFile(path);
        var lines = text.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        for (int i = 0; i < lines.Count; i++)
            lines[i] = lines[i].TrimEnd('\r');

        int ifmt = string.Equals(format, "simple", StringComparison.OrdinalIgnoreCase) ? 1 : 0;
        if (ifmt == 0)
            ifmt = lines.Count > 0 && lines[0].StartsWith('>') ? 2 : 1;

        if (ifmt == 1)
        {
            uint npat = 0;
            foreach (var line in lines)
            {
                npat++;
                Add(list, namestr + npat.ToString(CultureInfo.InvariantCulture), line, mismatches);
            }

            return list;
        }

        string name = "";
        var pat = new System.Text.StringBuilder();
        int mismatch = mismatches;
        foreach (var raw in lines)
        {
            string line = raw;
            if (line.Length > 0 && line[0] == '>')
            {
                if (name.Length > 0)
                {
                    Add(list, name, pat.ToString(), mismatch);
                    name = "";
                    pat.Clear();
                    mismatch = mismatches;
                }

                line = line[1..];
                var m = MismatchTag.Match(line);
                if (m.Success)
                {
                    name = m.Groups[1].Value;
                    mismatch = uint.TryParse(name, NumberStyles.None, CultureInfo.InvariantCulture, out uint v) ? (int)v : 0;
                    line = line[..m.Index].TrimEnd(' ', '\t', '\n', '\r', '\f', '\v');
                }

                name = line;
                if (name.Length == 0)
                    name = patternName;
            }
            else
            {
                pat.Append(line);
            }
        }

        if (name.Length == 0)
            name = patternName;
        Add(list, name, pat.ToString(), mismatch);
        return list;
    }

    // ajPatternSeqNewList: an empty pattern is not added; an empty name becomes pattern<n+1>
    private static void Add(List<EmbossProgramPattern> list, string name, string pattern, int mismatch)
    {
        if (pattern.Length == 0)
            return;
        if (name.Length == 0)
            name = "pattern" + (list.Count + 1).ToString(CultureInfo.InvariantCulture);
        list.Add(new EmbossProgramPattern(name, pattern, mismatch));
    }
}
