using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The program dies on a value of an ACD item (<c>acdSet*</c>: "Died: prog terminated: Bad value for '-name' with -auto defined"), after the messages of the library.</summary>
public sealed class EmbossAcdDieException(IReadOnlyList<string> messages) : Exception("The EMBOSS program died on a command-line value.")
{
    public IReadOnlyList<string> Messages { get; } = messages;
}

/// <summary>
/// The command line of an EMBOSS program resolved by the ACD library (<c>acdArgsParse</c>, <see cref="AcdCommandLine"/>) and read as <c>acdSetInt</c> / <c>acdSetBool</c> / <c>acdSetListSingle</c> read the items:
/// the typed values with the library's messages when a value is invalid (B09 WP61, F280).
/// </summary>
internal sealed class EmbossAcdArgs
{
    private readonly AcdProfile _profile;
    private readonly AcdCommandLine _line;
    private readonly string _program;

    private EmbossAcdArgs(string program, AcdProfile profile, AcdCommandLine line, IReadOnlyList<string> args)
    {
        _program = program;
        _profile = profile;
        _line = line;
        // ajUtilGetCmdline: the arguments of the run (the oracle runs every program with -auto)
        // acdArgSave: the qualifiers as typed (a value that is no word, or holds one of * ? [ ] { } | ! & ^, in double quotes)
        CommandLine = "-auto " + string.Join(' ', args.Select(QuoteToken));
    }

    private static string QuoteToken(string t)
    {
        if (t.Length == 0)
            return "\"\"";
        bool word = t[0] == '-' || !(t.Any(char.IsWhiteSpace) || t.IndexOfAny(['*', '?', '[', ']', '{', '}', '|', '!', '&', '^']) >= 0);
        return word ? t : "\"" + t + "\"";
    }

    /// <summary>The arguments of the run as the library keeps them (<c>ajUtilGetCmdline</c>).</summary>
    public string CommandLine { get; }

    /// <summary>Resolves the tokens of a command line (the qualifiers; the parameters are given as <c>-name value</c>).</summary>
    /// <exception cref="ArgumentException">A token the library would refuse (unknown or ambiguous qualifier, missing value).</exception>
    public static EmbossAcdArgs Parse(string program, IReadOnlyList<string> args)
    {
        var profile = EmbossAcdProfiles.Get(program);
        var line = new AcdCommandLine(profile);
        line.Run(args);
        return new EmbossAcdArgs(program, profile, line, args);
    }

    private AcdItem? Find(string name) => _profile.List.FirstOrDefault(i => i.Master is null && i.Name == name);

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
        if (raw is null)
            return dflt;
        if (!AcdCommandLine.TryToInt(raw, out int v))
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
            throw Die(name, $"Error: Invalid boolean value '{raw}'");
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
            throw Die("ossingle", $"Error: Invalid boolean value '{single}'");

        return new EmbossOutseqOptions(Raw(master) ?? "stdout", format, isSingle, RawAssoc(master, "osdirectory"), RawAssoc(master, "osextension"), RawAssoc(master, "osdbname"), RawAssoc(master, "osname")) { CommandLine = CommandLine, RunDate = EmbossRunEnvironment.Current.Now };
    }
}
