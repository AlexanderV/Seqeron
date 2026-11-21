namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// One item of the ACD list of an EMBOSS program as <c>ajacd.c</c> sees it when it resolves a command-line qualifier: a parameter (<c>[-asequence]</c>), a standard / additional / advanced or general qualifier, or an
/// associated qualifier of a master item (<c>-sbegin1</c> belongs to <c>-asequence</c>; <see cref="PNum"/> = the number of the parameter it is tied to, 0 = unnumbered).
/// </summary>
internal sealed class AcdItem(string name, int pnum, bool isParam, string type, AcdItem? master)
{
    public string Name { get; } = name;
    public int PNum { get; } = pnum;
    public bool IsParam { get; } = isParam;
    public string Type { get; } = type;

    /// <summary>The item this one is an associated qualifier of (<c>pa-&gt;Assoc</c>), or null.</summary>
    public AcdItem? Master { get; } = master;

    /// <summary>The associated qualifiers of this item in list order (<c>AssocQuals</c> ... while <c>Assoc</c>).</summary>
    public List<AcdItem> Assoc { get; } = [];

    public bool IsBoolean => Type is "boolean" or "toggle";
}

/// <summary>
/// The qualifier universe of one EMBOSS 6.6.0 program (<see cref="EmbossAcdProfiles"/>: read from <c>&lt;program&gt; -help -verbose</c>): the ACD list the command-line resolver matches a qualifier against.
/// </summary>
internal sealed class AcdProfile
{
    /// <summary>All items in the order of <c>acdList</c>: the general qualifiers, then each parameter / qualifier preceded by its associated qualifiers.</summary>
    public IReadOnlyList<AcdItem> List { get; }

    /// <summary>The number of parameters (<c>acdNParam</c>).</summary>
    public int ParamCount { get; }

    /// <summary>The sequence-input parameters (type sequence / seqall / seqset) in parameter order.</summary>
    public IReadOnlyList<AcdItem> Inputs { get; }

    internal AcdProfile(string program, IReadOnlyList<AcdItem> list)
    {
        Program = program;
        List = list;
        ParamCount = list.Count(i => i.IsParam);
        Inputs = [.. list.Where(i => i.IsParam && i.Type is "sequence" or "seqall" or "seqset" or "seqsetall").OrderBy(i => i.PNum)];
    }

    public string Program { get; }

    /// <summary>Parses the table lines (<c>P name pnum type</c> | <c>Q name type</c> | <c>G name type</c> | <c>A master name pnum type</c>) of <see cref="EmbossAcdProfiles"/>.</summary>
    internal static AcdProfile Parse(string program, string table)
    {
        var general = new List<AcdItem>();
        var main = new List<AcdItem>();
        var byName = new Dictionary<string, AcdItem>(StringComparer.Ordinal);
        var lines = table.Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
        foreach (var line in lines)
        {
            var f = line.Split(' ');
            switch (f[0])
            {
                case "P":
                    var p = new AcdItem(f[1], int.Parse(f[2], System.Globalization.CultureInfo.InvariantCulture), true, f[3], null);
                    main.Add(p);
                    byName[p.Name] = p;
                    break;
                case "Q":
                    var q = new AcdItem(f[1], 0, false, f[2], null);
                    main.Add(q);
                    byName[q.Name] = q;
                    break;
                case "G":
                    general.Add(new AcdItem(f[1], 0, false, f[2], null));
                    break;
            }
        }

        foreach (var line in lines)
        {
            var f = line.Split(' ');
            if (f[0] == "A")
            {
                var master = byName[f[1]];
                var item = new AcdItem(f[2], int.Parse(f[3], System.Globalization.CultureInfo.InvariantCulture), false, f[4], master);
                master.Assoc.Add(item);
            }
        }

        var list = new List<AcdItem>(general);
        foreach (var item in main)
        {
            list.AddRange(item.Assoc);                  // the associated qualifiers precede their master in acdList (measured: fuzznuc "-se" is ambiguous between "send,sequence")
            list.Add(item);
        }

        return new AcdProfile(program, list);
    }
}

internal static partial class EmbossAcdProfiles
{
    private static readonly Dictionary<string, AcdProfile> Cache = new(StringComparer.Ordinal);

    /// <summary>The names of the programs with a profile (the 16 ported programs).</summary>
    internal static IEnumerable<string> Programs => Table.Select(t => t.Program);

    private static readonly string[] SequenceQualifierTable =
    [
        "sbegin integer", "send integer", "sreverse boolean", "sask boolean", "snucleotide boolean", "sprotein boolean", "slower boolean", "supper boolean", "scircular boolean", "squick boolean",
        "sformat string", "iquery string", "ioffset integer", "sdbname string", "sid string", "ufo string", "fformat string", "fopenfile string",
    ];

    /// <summary>A profile with the sequence inputs and the general qualifiers only (no program-specific qualifier): the universe of <see cref="EmbossSequenceQualifiers.Parse(string?, string)"/>, which has no program.</summary>
    internal static AcdProfile SequenceOnly(int inputs)
    {
        var sb = new System.Text.StringBuilder();
        foreach (var g in new[] { "auto", "stdout", "filter", "options", "debug", "verbose", "help", "warning", "error", "fatal", "die", "version" })
            sb.Append("G ").Append(g).Append(" boolean\n");
        var names = inputs == 1 ? new[] { "sequence" } : new[] { "asequence", "bsequence" };
        for (int i = 0; i < inputs; i++)
            sb.Append("P ").Append(names[i]).Append(' ').Append(i + 1).Append(" seqall\n");
        for (int i = 0; i < inputs; i++)
        {
            foreach (var q in SequenceQualifierTable)
            {
                var f = q.Split(' ');
                sb.Append("A ").Append(names[i]).Append(' ').Append(f[0]).Append(' ').Append(i + 1).Append(' ').Append(f[1]).Append('\n');
            }
        }

        return AcdProfile.Parse("sequence", sb.ToString());
    }

    /// <summary>The profile of <paramref name="program"/> (case-sensitive EMBOSS name).</summary>
    /// <exception cref="ArgumentException">No profile for this program.</exception>
    internal static AcdProfile Get(string program)
    {
        lock (Cache)
        {
            if (Cache.TryGetValue(program, out var cached))
                return cached;
            foreach (var (name, items) in Table)
            {
                if (string.Equals(name, program, StringComparison.Ordinal))
                    return Cache[program] = AcdProfile.Parse(name, items);
            }
        }

        throw new ArgumentException($"No ACD profile for the EMBOSS program '{program}'.", nameof(program));
    }
}
