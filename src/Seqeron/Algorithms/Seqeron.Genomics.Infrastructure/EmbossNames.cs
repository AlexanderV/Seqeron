#pragma warning disable S1121, S3358, S3776, S107, S1854, S1066, S2325, S1172

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>A database definition of the EMBOSS configuration (<c>DBNAME name [ attribute: value ... ]</c>): the attributes by name in lower case; the list attributes (<c>field</c>, <c>example</c>, <c>special</c>, <c>edam*</c>, <c>taxon</c>) keep every value.</summary>
internal sealed class EmbossDatabase(string name)
{
    public string Name { get; } = name;

    /// <summary>The single-value attributes (<c>format method type directory filename fields ...</c>).</summary>
    public Dictionary<string, string> Attributes { get; } = [];

    /// <summary>The attributes that can be given more than once.</summary>
    public Dictionary<string, List<string>> Lists { get; } = [];

    public string? Get(string attribute) => Attributes.GetValueOrDefault(attribute);
}

/// <summary>
/// The EMBOSS configuration (<c>ajnam.c</c>: <c>namProcessFile</c> / <c>namListParse</c> of a file such as <c>.embossrc</c> or <c>emboss.default</c>): the words of the file, the variables of <c>SET</c> / <c>ENV</c> (resolved in the values as <c>$name</c>,
/// <c>ajNamResolve</c>), the aliases and the <c>DBNAME</c> definitions (<c>ajNamDatabase</c>, <c>ajNamAliasDatabase</c>, <c>ajNamDbData</c>). A5-8c3h, B09 F256.
/// </summary>
internal sealed class EmbossNames
{
    /// <summary>The attributes of a database definition (<c>namDbAttrs</c>): name and whether it takes a list of values.</summary>
    private static readonly (string Name, bool List)[] DbAttrs =
    [
        ("format", false), ("method", false), ("type", false), ("accession", false), ("app", false), ("appall", false), ("appentry", false), ("appquery", false), ("caseidmatch", false), ("comment", false), ("dbalias", false),
        ("directory", false), ("example", true), ("exclude", false), ("field", true), ("fields", false), ("filters", false), ("filename", false), ("formatall", false), ("formatentry", false), ("formatquery", false),
        ("hasaccession", false), ("httpversion", false), ("identifier", false), ("indexdirectory", false), ("methodall", false), ("methodentry", false), ("methodquery", false), ("namespace", false), ("organisms", false),
        ("proxy", false), ("query", false), ("release", false), ("return", false), ("sequence", false), ("serverversion", false), ("special", true), ("url", false), ("view", false), ("edamdat", true), ("edamfmt", true),
        ("edamid", true), ("edamtpc", true), ("taxon", true),
    ];

    private static readonly Regex VarExp = new(@"^\$([A-Za-z0-9_]+)", RegexOptions.CultureInvariant);

    private readonly Dictionary<string, string> _vars = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, string> _aliases = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, EmbossDatabase> _databases = new(StringComparer.OrdinalIgnoreCase);

    /// <summary>The messages the library printed while reading the definitions (<c>Warning: ...</c>, <c>Error: ...</c>).</summary>
    public List<string> Messages { get; } = [];

    /// <summary>The definitions in the text of the files (in the order the library reads them: later definitions of a name replace earlier ones).</summary>
    public static EmbossNames Parse(params string[] texts)
    {
        var names = new EmbossNames();
        foreach (string text in texts)
            names.ProcessFile(text);
        return names;
    }

    /// <summary><c>ajNamDatabase</c>: the name is a defined database (any case).</summary>
    public bool IsDatabase(string name) => _databases.ContainsKey(name);

    /// <summary><c>ajNamAliasDatabase</c>: the name is a database, or an alias of one (the alias chain is followed).</summary>
    public bool AliasDatabase(ref string name)
    {
        if (_databases.ContainsKey(name))
            return true;
        string tryname = name;
        int guard = 0;
        while (_aliases.TryGetValue(tryname, out string? alias) && guard++ < 100)
        {
            if (_databases.ContainsKey(alias))
            {
                name = alias;
                return true;
            }

            tryname = alias;
        }

        return false;
    }

    /// <summary>The definition of a database (the name as defined, any case); null when there is none.</summary>
    public EmbossDatabase? Database(string name) => _databases.GetValueOrDefault(name);

    /// <summary><c>ajNamGetValueS</c> of a variable (any case).</summary>
    public string? Variable(string name) => _vars.GetValueOrDefault(name);

    /// <summary><c>ajNamResolve</c>: a value that starts with <c>$name</c> gets the value of the variable (an unknown one is dropped with a warning).</summary>
    public string Resolve(string value)
    {
        var m = VarExp.Match(value);
        if (!m.Success)
            return value;
        string rest = value[m.Length..];
        if (_vars.TryGetValue(m.Groups[1].Value, out string? v))
            return v + rest;
        Messages.Add($"Warning: Variable unknown in '{value}'");
        return rest;
    }

    private static int DbAttr(string name)
    {
        string lower = name.ToLowerInvariant();
        int found = 0, j = -1;
        for (int i = 0; i < DbAttrs.Length; i++)
        {
            if (DbAttrs[i].Name == lower)
                return i;
            if (DbAttrs[i].Name.StartsWith(lower, StringComparison.Ordinal))
            {
                found++;
                j = i;
            }
        }

        return found == 1 ? j : -1;
    }

    private static bool IsDbname(string name)
    {
        if (name.Length < 2 || !char.IsAsciiLetter(name[0]))
            return false;
        foreach (char c in name.AsSpan(1))
        {
            if (!char.IsAsciiLetterOrDigit(c) && c != '_')
                return false;
        }

        return true;
    }

    // namProcessFile: the words of the file (blanks split, a quote opens a part that the parse joins again, "]" is a word of its own) and the parse of the words
    private void ProcessFile(string text)
    {
        var words = new List<string>();
        var lines = text.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        var word = new StringBuilder();
        char quote = '\0';
        foreach (string raw in lines)
        {
            string line = raw;
            if (line.EndsWith('\r'))
                line = line[..^1];
            if (line.Length == 0 || line[0] == '#')
                continue;
            foreach (char c in line)
            {
                if (c is ' ' or '\t')
                {
                    if (word.Length > 0)
                    {
                        words.Add(word.ToString());
                        word.Clear();
                    }

                    continue;
                }

                if (c is '\'' or '"')
                {
                    word.Append(c);
                    if (quote != '\0')
                    {
                        if (quote == c)
                            quote = '\0';
                    }
                    else
                    {
                        quote = c;
                    }
                }
                else if (quote == '\0' && word.Length > 0 && c == ']')
                {
                    words.Add(word.ToString());
                    word.Clear();
                    words.Add("]");
                }
                else
                {
                    word.Append(c);
                }
            }

            if (word.Length > 0)
            {
                words.Add(word.ToString());
                word.Clear();
            }
        }

        ListParse(words);
    }

    private enum ParseType { None, Env, Db, Alias }

    // namListParse for the definitions of variables, aliases and databases
    private void ListParse(List<string> words)
    {
        var type = ParseType.None;
        string? name = null, value = null;
        char quoteOpen = '\0', quoteClose = '\0';
        EmbossDatabase? db = null;
        int dbInput = -1;
        bool saveit = false, dbsave = false;
        foreach (string curword0 in words)
        {
            string curword = curword0;
            if (type == ParseType.None)
            {
                string lower = curword.TrimEnd(':').ToLowerInvariant();
                if (lower.Length > 0 && "env".StartsWith(lower, StringComparison.Ordinal) || lower.Length > 0 && "setenv".StartsWith(lower, StringComparison.Ordinal) || lower == "set")
                    type = ParseType.Env;
                else if (lower.Length > 0 && "dbname".StartsWith(lower, StringComparison.Ordinal) || lower == "db")
                    type = ParseType.Db;
                else if (lower.Length > 0 && "alias".StartsWith(lower, StringComparison.Ordinal))
                    type = ParseType.Alias;
                else if (lower.Length > 0 && ("server".StartsWith(lower, StringComparison.Ordinal) || "resource".StartsWith(lower, StringComparison.Ordinal) || "include".StartsWith(lower, StringComparison.Ordinal) || "if".StartsWith(lower, StringComparison.Ordinal) || "ifdef".StartsWith(lower, StringComparison.Ordinal) || "else".StartsWith(lower, StringComparison.Ordinal) || "endif".StartsWith(lower, StringComparison.Ordinal)))
                    throw new NotSupportedException($"The definition type '{curword}' (server, resource, include, if) of the EMBOSS configuration is not read by this port.");
                else
                    Messages.Add($"Error: Invalid definition type '{curword}'");
            }
            else if (quoteOpen != '\0' && curword == "]")
            {
                Messages.Add($"Error: ']' found, unclosed quotes in '{value}'");
                quoteOpen = '\0';
                type = ParseType.None;
            }
            else if (quoteOpen != '\0')
            {
                value = string.Concat(value, " ", curword);
                if (curword.Length > 0 && curword[^1] == quoteClose)
                {
                    value = value[..^1];
                    quoteOpen = quoteClose = '\0';
                    if (type is ParseType.Env or ParseType.Alias)
                        saveit = true;
                    else if (type == ParseType.Db)
                        dbsave = true;
                }
            }
            else if (type is ParseType.Env or ParseType.Alias)
            {
                if (name is not null && value is not null)
                {
                    saveit = true;
                }
                else if (name is not null)
                {
                    if (curword.StartsWith('\''))
                        quoteOpen = quoteClose = '\'';
                    else if (curword.StartsWith('"'))
                        quoteOpen = quoteClose = '"';
                    value = curword;
                    if (quoteOpen != '\0')
                    {
                        value = value[1..];
                        if (value.Length == 0)
                            Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                    }
                    else
                    {
                        saveit = true;
                    }

                    if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                    {
                        quoteOpen = quoteClose = '\0';
                        saveit = true;
                        value = value[..^1];
                    }

                    if (saveit)
                        value = Resolve(value);
                }
                else
                {
                    name = curword;
                }
            }
            else if (type == ParseType.Db)
            {
                if (curword == "[")
                {
                    db = new EmbossDatabase(name ?? "");
                }
                else if (curword == "]")
                {
                    saveit = true;
                }
                else if (name is not null)
                {
                    if (dbInput >= 0)
                    {
                        if (curword.StartsWith('\''))
                            quoteOpen = quoteClose = '\'';
                        else if (curword.StartsWith('"'))
                            quoteOpen = quoteClose = '"';
                        value = curword;
                        if (quoteOpen != '\0')
                        {
                            value = value[1..];
                            if (value.Length == 0)
                                Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                        }
                        else
                        {
                            dbsave = true;
                        }

                        if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                        {
                            quoteOpen = quoteClose = '\0';
                            value = value[..^1];
                            dbsave = true;
                        }

                        if (quoteOpen == '\0')
                            dbsave = true;
                        if (dbsave)
                            value = Resolve(value);
                    }
                    else if (curword.EndsWith(':'))
                    {
                        dbInput = DbAttr(curword.TrimEnd(':'));
                        if (dbInput < 0)
                            Messages.Add($"Warning: .embossrc: bad attribute '{curword.TrimEnd(':').ToLowerInvariant()}' for database '{name}'");
                    }
                    else
                    {
                        Messages.Add($"Warning: .embossrc: unexpected token '{curword}' for database '{name}'");
                    }
                }
                else
                {
                    name = curword;
                    if (!IsDbname(name))
                        Messages.Add($"Error: Invalid database name '{name}'");
                }
            }

            if (dbsave)
            {
                if (dbInput >= 0 && db is not null)
                {
                    var (attrName, isList) = DbAttrs[dbInput];
                    if (isList)
                    {
                        if (!db.Lists.TryGetValue(attrName, out var list))
                            db.Lists[attrName] = list = [];
                        list.Add(value ?? "");
                    }
                    else
                    {
                        if (db.Attributes.TryGetValue(attrName, out string? previous))
                            Messages.Add($"Warning: duplicate db keyword value for '{attrName}' = '{previous}'");
                        db.Attributes[attrName] = value ?? "";
                    }
                }

                value = null;
                dbInput = -1;
                dbsave = false;
            }

            if (saveit)
            {
                if (type == ParseType.Db && name is not null && db is not null)
                    _databases[name] = db;
                else if (type == ParseType.Env && name is not null)
                    _vars[name] = value ?? "";
                else if (type == ParseType.Alias && name is not null)
                    _aliases[name] = value ?? "";
                name = null;
                value = null;
                db = null;
                saveit = false;
                type = ParseType.None;
            }
        }
    }
}
