// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>The program dies on a value of an ACD item (<c>acdSet*</c>: "Died: prog terminated: Bad value for '-name' with -auto defined"), after the messages of the library.</summary>
public sealed class EmbossAcdDieException(IReadOnlyList<string> messages) : Exception("The EMBOSS program died on a command-line value.")
{
    public IReadOnlyList<string> Messages { get; } = messages;
}

/// <summary>
/// The command line of an EMBOSS program resolved by the ACD library (<c>acdArgsParse</c>, <see cref="AcdCommandLine"/>) and read as <c>acdSetInt</c> / <c>acdSetBool</c> / <c>acdSetListSingle</c> read the items:
/// the typed values with the library's messages when a value is invalid (B09 WP61, F280).
/// </summary>
public sealed class EmbossAcdArgs
{
    private readonly AcdProfile _profile;
    private readonly AcdCommandLine _line;
    private readonly string _program;

    private EmbossAcdArgs(string program, AcdProfile profile, AcdCommandLine line)
    {
        _program = program;
        _profile = profile;
        _line = line;
        // ajUtilGetCmdline: the arguments of the run (the oracle runs every program with -auto)
        // acdArgSave: the qualifiers as typed (a value that is no word, or holds one of * ? [ ] { } | ! & ^, in double quotes)
        CommandLine = line.SavedCommandLine;
    }

    /// <summary>The arguments of the run as the library keeps them (<c>ajUtilGetCmdline</c>).</summary>
    public string CommandLine { get; }

    /// <summary>The arguments of the run one entry per argument, as the header of a report prints them.</summary>
    public IReadOnlyList<string> CommandLineEntries => _line.SavedEntries;

    /// <summary>Resolves the tokens of a command line (the qualifiers; the parameters are given as <c>-name value</c>).</summary>
    /// <exception cref="ArgumentException">A token the library would refuse (unknown or ambiguous qualifier, missing value).</exception>
    public static EmbossAcdArgs Parse(string program, IReadOnlyList<string> args)
    {
        var profile = EmbossAcdProfiles.Get(program);
        var line = new AcdCommandLine(profile);
        line.Run(args);
        return new EmbossAcdArgs(program, profile, line);
    }

    public AcdItem? Find(string name) => _profile.List.FirstOrDefault(i => i.Master is null && i.Name == name);

    private AcdItem? FindAssoc(string master, string name, int pnum) => _profile.List.FirstOrDefault(i => i.Master?.Name == master && i.Name == name && (pnum == 0 || i.PNum == pnum));

    /// <summary>The text of an item (a parameter or a qualifier) as the command line gave it; null = not given.</summary>
    public string? Raw(string name) => Find(name) is { } item ? _line.Value(item) : null;

    /// <summary>The text of an associated qualifier (<c>-osformat</c> of the output <c>outseq</c>, number 0: unnumbered); null = not given.</summary>
    public string? RawAssoc(string master, string name, int pnum = 0) => FindAssoc(master, name, pnum) is { } item ? _line.Value(item) : null;

    /// <summary>Dies as the library does when the value of <paramref name="name"/> is bad: the messages, then the death message.</summary>
    public EmbossAcdDieException Die(string name, params string[] messages) => new([.. messages, $"Died: {_program} terminated: Bad value for '-{name}' with -auto defined"]);

    /// <summary><c>acdSetInt</c>: a decimal integer (<c>ajStrToInt</c>).</summary>
    public int Int(string name, int dflt)
    {
        string? raw = Raw(name);
        if (raw is null)
            return dflt;
        if (!AcdCommandLine.TryToInt(raw, out int v))
            throw Die(name, $"Error: Invalid integer value '{raw}'");
        return v;
    }

    /// <summary>
    /// <c>acdSetInt</c> with the attributes <c>minimum</c> / <c>maximum</c> and <c>failrange: N</c> (the usual case of the programs): a value outside the limits is reset to the limit with the warning
    /// "integer value out of range V less than (reset to) L" / "more than (reset to) L"; the default is not checked.
    /// </summary>
    public int IntLimited(string name, int dflt, long min, long max, EmbossFeatureReadContext context)
    {
        string? raw = Raw(name);
        int v = dflt;                                   // acdSetInt checks the default against the limits as well
        if (raw is not null && !AcdCommandLine.TryToInt(raw, out v))
            throw Die(name, $"Error: Invalid integer value '{raw}'");
        if (v < min)
        {
            context.Warn($"integer value out of range {v} less than (reset to) {min}");
            return (int)min;
        }

        if (v > max)
        {
            context.Warn($"integer value out of range {v} more than (reset to) {max}");
            return (int)max;
        }

        return v;
    }

    /// <summary>
    /// <c>acdSetInt</c> with limits that are not failing (<c>failrange: N</c>) and the attributes <c>minimum</c> / <c>maximum</c>, the usual case of the programs: a value below the minimum is reset to it with the warning
    /// "integer value out of range V less than (reset to) L", and the value then checked against the maximum ("... more than (reset to) L"): both can happen for one value (the maximum may be below the minimum: a window longer than the sequence).
    /// B09 WP63, F282.
    /// </summary>
    public int IntClamped(string name, int dflt, long min, long max, EmbossFeatureReadContext context)
    {
        string? raw = Raw(name);
        int v = dflt;
        if (raw is not null && !AcdCommandLine.TryToInt(raw, out v))
            throw Die(name, $"Error: Invalid integer value '{raw}'");
        long value = v;
        if (value < min)
        {
            context.Warn($"integer value out of range {value} less than (reset to) {min}");
            value = min;
        }

        if (value > max)
        {
            context.Warn($"integer value out of range {value} more than (reset to) {max}");
            value = max;
        }

        return (int)value;
    }

    /// <summary>The integer value of an associated qualifier (<c>-rmaxall</c> of the report <c>outfile</c>); null = not given.</summary>
    public int? IntAssoc(string master, string name)
    {
        string? raw = RawAssoc(master, name);
        if (raw is null)
            return null;
        if (!AcdCommandLine.TryToInt(raw, out int v))
            throw Die(name, $"Error: Invalid integer value '{raw}'");
        return v;
    }

    /// <summary>The boolean value of an associated qualifier (<c>-raccshow</c> of the report <c>outfile</c>); null = not given.</summary>
    public bool? BoolAssoc(string master, string name)
    {
        string? raw = RawAssoc(master, name);
        if (raw is null)
            return null;
        if (!AcdCommandLine.TryToBool(raw, out bool v))
            throw Die(name, $"Error: Invalid Y/N value '{raw}'");
        return v;
    }

    /// <summary><c>acdSetFloat</c>.</summary>
    public float Float(string name, float dflt)
    {
        string? raw = Raw(name);
        if (raw is null)
            return dflt;
        if (!float.TryParse(raw, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out float v))
            throw Die(name, $"Error: Invalid floating point value '{raw}'");
        return v;
    }

    /// <summary><c>acdSetBool</c>.</summary>
    public bool Bool(string name, bool dflt)
    {
        string? raw = Raw(name);
        if (raw is null)
            return dflt;
        if (!AcdCommandLine.TryToBool(raw, out bool v))
            throw Die(name, $"Error: Invalid Y/N value '{raw}'");
        return v;
    }

    /// <summary>The string of an item.</summary>
    public string String(string name, string dflt) => Raw(name) ?? dflt;

    /// <summary><c>acdSetListSingle</c> with short codes: the value must be one of <paramref name="codes"/> (case-insensitive).</summary>
    public string Menu(string name, string dflt, params string[] codes)
    {
        string? raw = Raw(name);
        if (raw is null)
            return dflt;
        foreach (string c in codes)
        {
            if (string.Equals(c, raw, StringComparison.OrdinalIgnoreCase))
                return c;
        }

        // a unique abbreviation of a code is accepted
        var starts = codes.Where(c => c.StartsWith(raw, StringComparison.OrdinalIgnoreCase)).ToList();
        if (raw.Length > 0 && starts.Count == 1)
            return starts[0];
        if (raw.Length > 0 && starts.Count > 1)
            throw Die(name, $"Error: '{raw}' is ambiguous ({string.Join(",", starts)})", $"Error: Bad menu option '{raw}'");

        throw Die(name, $"Error: '{raw}' is not a valid menu option", "Accepted short codes are: " + string.Join(",", codes), $"Error: Bad menu option '{raw}'");
    }

    /// <summary>
    /// <c>acdSetList</c> / <c>acdListValue</c>: a list item with several values (<c>-point 2,3</c>): every value of the reply (separated by commas) is matched against the codes and the descriptions of
    /// <paramref name="values"/> (<c>"0:None; 1:Any of the following; 2:Insertions"</c>: delimiter <c>;</c>, code delimiter <c>:</c>): an exact code, an exact description, then a unique prefix; <c>*</c> is every code.
    /// The number of values is limited to <paramref name="min"/> .. <paramref name="max"/>. B09 WP62, F281.
    /// </summary>
    public IReadOnlyList<string> List(string name, string dflt, string values, int min, int max, string delimiter = ";")
    {
        string reply = Raw(name) ?? dflt;
        var lines = new List<(string Code, string Desc)>();
        foreach (string line in values.Split(delimiter))
        {
            int colon = line.IndexOf(':', StringComparison.Ordinal);
            string code = colon < 0 ? line : line[..colon];
            string desc = colon < 0 ? "" : line[(colon + 1)..];
            lines.Add((code.Trim(' ', '\t', '\n', '\r'), desc.Trim(' ', '\t', '\n', '\r')));
        }

        var result = new List<string>();
        var messages = new List<string>();
        bool ok = true;
        foreach (string repstr in AcdTokens(reply, ','))
        {
            if (repstr == "*")
            {
                foreach (var (code, _) in lines)
                    result.Add(code);
                continue;
            }

            int ifound = 0, jfound = 0;
            string hit1 = "", hit2 = "";
            var ambig = new List<string>();
            var valid = new List<string>();
            foreach (var (code, desc) in lines)
            {
                valid.Add(code);
                if (string.Equals(code, repstr, StringComparison.OrdinalIgnoreCase))
                {
                    ifound = 1;
                    hit1 = code;
                    break;
                }

                if (string.Equals(desc, repstr, StringComparison.OrdinalIgnoreCase))
                {
                    jfound = 1;
                    hit2 = code;
                    break;
                }

                if (repstr.Length > 0 && code.StartsWith(repstr, StringComparison.OrdinalIgnoreCase))
                {
                    ifound++;
                    hit1 = code;
                    ambig.Add(code);
                }

                if (repstr.Length > 0 && desc.StartsWith(repstr, StringComparison.OrdinalIgnoreCase))
                {
                    jfound++;
                    hit2 = code;
                    ambig.Add(desc);
                }
            }

            if (ifound == 1)
                result.Add(hit1);
            else if (jfound == 1)
                result.Add(hit2);
            else
            {
                messages.Add(ifound > 0 || jfound > 0 ? $"Error: '{repstr}' is ambiguous ({string.Join(",", ambig)})" : $"Error: '{repstr}' is not a valid menu option");
                if (ifound == 0 && jfound == 0)
                    messages.Add("Accepted short codes are: " + string.Join(",", valid));
                ok = false;
                break;
            }
        }

        if (ok && result.Count < min)
        {
            messages.Add(min <= 1 ? $"Error: Menu needs {min} value" : $"Error: Menu needs {min} values");
            ok = false;
        }

        if (ok && result.Count > max)
        {
            messages.Add(max <= 1 ? $"Error: Menu allows no more than {max} value" : $"Error: Menu allows no more than {max} values");
            ok = false;
        }

        if (!ok)
        {
            messages.Add($"Error: Bad menu option '{reply}'");
            throw Die(name, [.. messages]);
        }

        return result;
    }

    /// <summary><c>ajStrTokenNextParse</c>: the tokens of a text separated by runs of the delimiter (a leading delimiter gives one empty token; a trailing run gives none).</summary>
    internal static IEnumerable<string> AcdTokens(string text, char delimiter)
    {
        int pos = 0;
        while (pos < text.Length)
        {
            int end = text.IndexOf(delimiter, pos);
            if (end < 0)
                end = text.Length;
            yield return text[pos..end];
            pos = end;
            while (pos < text.Length && text[pos] == delimiter)
                pos++;
        }
    }

    /// <summary>The qualifiers of the sequence input <paramref name="master"/> (<c>-sbegin1</c> ...), for a program that reads its sequences from files.</summary>
    public EmbossSequenceQualifiers Sequence(string master)
    {
        var item = Find(master) ?? throw new ArgumentException($"{_program} has no input '{master}'.");
        try
        {
            return EmbossSequenceQualifiers.FromAssociated(item, _line.Value, fileInput: true);
        }
        catch (ArgumentException e)
        {
            throw Die(master, "Error: " + e.Message);
        }
    }

    /// <summary>The output qualifiers of the sequence output <paramref name="master"/> (<c>-osformat -ossingle -osextension -osdirectory -osdbname -osname</c>) and its value.</summary>
    public EmbossOutseqOptions Seqout(string master)
    {
        string? format = RawAssoc(master, "osformat");
        string? single = RawAssoc(master, "ossingle");
        bool isSingle = false;
        if (single is not null && !AcdCommandLine.TryToBool(single, out isSingle))
            throw Die("ossingle", $"Error: Invalid Y/N value '{single}'");

        return new EmbossOutseqOptions(Raw(master) ?? "stdout", format, isSingle, RawAssoc(master, "osdirectory"), RawAssoc(master, "osextension"), RawAssoc(master, "osdbname"), RawAssoc(master, "osname")) { CommandLine = CommandLine, RunDate = EmbossRunEnvironment.Current.RunDate, RefToday = EmbossRunEnvironment.Current.Now };
    }
}
