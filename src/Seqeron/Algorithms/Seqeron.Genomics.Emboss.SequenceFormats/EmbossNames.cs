// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S3776, S107, S1854, S1066, S2325, S1172

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>A database definition of the EMBOSS configuration (<c>DBNAME name [ attribute: value ... ]</c>): the attributes by name in lower case; the list attributes (<c>field</c>, <c>example</c>, <c>special</c>, <c>edam*</c>, <c>taxon</c>) keep every value.</summary>
public sealed class EmbossDatabase(string name)
{
    public string Name { get; } = name;

    /// <summary>The single-value attributes (<c>format method type directory filename fields ...</c>).</summary>
    public Dictionary<string, string> Attributes { get; } = [];

    /// <summary>The attributes that can be given more than once.</summary>
    public Dictionary<string, List<string>> Lists { get; } = [];

    /// <summary>The short name of the definition file (<c>standard</c>, <c>user</c>, or the name of an included file without its directory): <c>NamPEntry.file</c>.</summary>
    public string Defined { get; init; } = "user";

    public string? Get(string attribute) => Attributes.GetValueOrDefault(attribute);
}

/// <summary>
/// The EMBOSS configuration (<c>ajnam.c</c>: <c>namProcessFile</c> / <c>namListParse</c> of a file such as <c>.embossrc</c> or <c>emboss.default</c>): the words of the file, the variables of <c>SET</c> / <c>ENV</c> (resolved in the values as <c>$name</c>,
/// <c>ajNamResolve</c>), the aliases and the <c>DBNAME</c> definitions (<c>ajNamDatabase</c>, <c>ajNamAliasDatabase</c>, <c>ajNamDbData</c>). A5-8c3h, B09 F256.
/// </summary>
public sealed class EmbossNames
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

    // the attributes of a server definition (namSvrAttrs) and of a resource definition (namRsAttrs)
    private static readonly string[] SvrAttrs =
    [
        "method", "type", "cachedirectory", "cachefile", "accession", "app", "appall", "appentry", "appquery", "caseidmatch", "comment", "directory", "field", "filters", "format", "formatall", "formatentry", "formatquery", "hasaccession", "httpversion", "identifier",
        "indexdirectory", "methodall", "methodentry", "methodquery", "proxy", "query", "return", "sequence", "serverversion", "special", "url", "view", "edamdat", "edamfmt", "edamid", "edamtpc", "taxon",
    ];

    private static readonly string[] RsAttrs = ["type", "release", "pagesize", "secpagesize", "cachesize", "seccachesize", "fields", "value"];

    private static readonly Regex VarExp = new(@"^\$([A-Za-z0-9_]+)", RegexOptions.CultureInvariant);

    private readonly Dictionary<string, string> _vars = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, string> _aliases = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, EmbossDatabase> _databases = new(StringComparer.OrdinalIgnoreCase);
    private readonly HashSet<string> _servers = new(StringComparer.OrdinalIgnoreCase);
    private readonly Dictionary<string, Dictionary<string, string>> _resources = new(StringComparer.OrdinalIgnoreCase);
    private Func<string, string?> _files = _ => null;
    private readonly HashSet<string> _included = new(StringComparer.Ordinal);
    private string _fileName = ".embossrc";
    private string _shortName = "user";
    private bool _ifNow = true;
    private readonly Stack<bool> _ifStack = new();

    /// <summary>The messages the library printed while reading the definitions (<c>Warning: ...</c>, <c>Error: ...</c>).</summary>
    public List<string> Messages { get; } = [];

    /// <summary><c>namErrorCount</c>: the errors of the definitions; the program stops ("Died: Error(s) in configuration files") when there are any.</summary>
    public int ErrorCount { get; private set; }

    /// <summary><c>ajNamServer</c>: a defined server.</summary>
    public bool IsServerName(string name) => _servers.Contains(name);

    private int _line = 1;

    // namError: ajErr("File %S line %d: %S")
    private void NamError(string message)
    {
        ErrorCount++;
        Messages.Add($"Error: File {_fileName} line {_line}: {message}");
    }

    /// <summary>The definitions in the text of the files (in the order the library reads them: later definitions of a name replace earlier ones).</summary>
    public static EmbossNames Parse(params string[] texts) => Parse(null, null, texts);

    /// <summary>The definitions in the texts of the files; <paramref name="files"/> gives the files that <c>INCLUDE</c> names (null when there are none), <paramref name="fileName"/> is the name of the files in the messages.</summary>
    public static EmbossNames Parse(Func<string, string?>? files, string? fileName, params string[] texts)
    {
        var names = new EmbossNames { _files = files ?? (_ => null) };
        foreach (string text in texts)
        {
            names._fileName = fileName ?? ".embossrc";
            names.ProcessFile(text);
        }

        return names;
    }

    /// <summary>
    /// The configuration the library always has: the definitions of <c>emboss.standard</c> of the installation (the ontologies, the taxonomy and the data resource catalogue as databases, the public servers; <c>$emboss_standard</c> is the installation
    /// directory) read first, then the texts of the user's files (B09 WP67, F309).
    /// </summary>
    public static EmbossNames ParseWithStandard(Func<string, string?>? files, string? fileName, params string[] texts)
    {
        var names = new EmbossNames { _files = files ?? (_ => null) };
        names._vars["emboss_standard"] = "/usr/share/EMBOSS";
        names._fileName = "emboss.standard";
        names._shortName = "standard";
        names.ProcessFile(EmbossDrcat.ReadResource("emboss.standard.gz"));
        names._shortName = "user";
        foreach (string text in texts)
        {
            names._fileName = fileName ?? ".embossrc";
            names.ProcessFile(text);
        }

        return names;
    }

    /// <summary>The names of the defined databases as defined (<c>ajNamListListDatabases</c>; the caller sorts them).</summary>
    public List<string> DatabaseNames() => [.. _databases.Values.Select(d => d.Name)];

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

    // namProcessFile: the words of the file (blanks split, a quote opens a part that the parse joins again, "]" is a word of its own) and the parse of the words. Returns false when the file had errors.
    private bool ProcessFile(string text)
    {
        int preverrors = ErrorCount;
        var words = new List<string>();
        var counts = new Queue<int>();
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
            counts.Enqueue(words.Count);                                   // the number of words before the line
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

        counts.Enqueue(words.Count);
        ListParse(words, counts);
        return ErrorCount <= preverrors;
    }

    private enum ParseType { None, Env, Svr, Db, Resource, IFile, Alias, If, IfDef, Else, EndIf }

    private ParseType _parseType = ParseType.None;
    private string? _name, _value;

    private static string TypeName(ParseType t) => t switch { ParseType.Env => "SET", ParseType.Svr => "SERVER", ParseType.Db => "DBNAME", ParseType.Resource => "RESOURCE", ParseType.IFile => "INCLUDE", ParseType.Alias => "ALIAS", ParseType.If or ParseType.IfDef => "IF", ParseType.Else => "ELSE", ParseType.EndIf => "ENDIF", _ => "unknown" };

    // ajCharPrefixS(text, word): the text starts with the word
    private static bool PrefixOf(string keyword, string word) => word.Length > 0 && keyword.StartsWith(word, StringComparison.Ordinal);

    private static int AttrIndex(string[] attrs, string name)
    {
        string lower = name.ToLowerInvariant();
        int found = 0, j = -1;
        for (int i = 0; i < attrs.Length; i++)
        {
            if (attrs[i] == lower)
                return i;
            if (attrs[i].StartsWith(lower, StringComparison.Ordinal))
            {
                found++;
                j = i;
            }
        }

        return found == 1 ? j : -1;
    }

    // namRsAttrFieldS: a resource of the type Index has the attributes NAMElen, NAMEpagesize ... for every field it indexes
    private static bool RsAttrField(Dictionary<string, string>? rs, string str)
    {
        if (rs is null || !rs.TryGetValue("type", out string? type) || type != "Index")
            return false;
        string prefix;
        if (str.EndsWith("len", StringComparison.Ordinal))
            prefix = str[..^3];
        else if (str.EndsWith("seccachesize", StringComparison.Ordinal))
            prefix = str[..^12];
        else if (str.EndsWith("cachesize", StringComparison.Ordinal))
            prefix = str[..^9];
        else if (str.EndsWith("secpagesize", StringComparison.Ordinal))
            prefix = str[..^11];
        else if (str.EndsWith("pagesize", StringComparison.Ordinal))
            prefix = str[..^8];
        else
            return false;
        if (prefix.Length == 0)
            return false;
        string fields = rs.TryGetValue("fields", out string? f) ? f : "id,acc,sv,key,des,org,nam,isa,xref";
        return fields.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries).Contains(prefix, StringComparer.Ordinal);
    }

    // ajNamGetValueS: a variable of the definitions (the environment is not part of this port)
    private bool GetValue(string name, out string value)
    {
        value = "";
        string upper = name.StartsWith("EMBOSS_", StringComparison.OrdinalIgnoreCase) ? name.ToUpperInvariant() : "EMBOSS_" + name.ToUpperInvariant();
        if (_vars.TryGetValue(upper, out string? v) || _vars.TryGetValue(name, out v))
        {
            value = v;
            return true;
        }

        return false;
    }

    // namListParse
    private void ListParse(List<string> words, Queue<int> counts)
    {
        var db = (EmbossDatabase?)null;
        Dictionary<string, List<string>>? svLists = null;
        Dictionary<string, string>? svAttrs = null, rsAttrs = null;
        int dbInput = -1, svrInput = -1, rsInput = -1;
        string? rsName = null;
        char quoteOpen = '\0', quoteClose = '\0';
        bool saveit = false, dbsave = false, svrsave = false, rssave = false, saveifblock = false;
        int wordcount = 0, linecount = 0, lineword = 0;
        _name = _value = null;
        _line = 1;
        string savedFile = _fileName;
        foreach (string curword0 in words)
        {
            string curword = curword0;
            while (counts.Count > 0 && lineword < wordcount)
            {
                lineword = counts.Dequeue();
                linecount++;
                _line = linecount - 1;
            }

            wordcount++;
            if (_parseType == ParseType.None)
            {
                string lower = curword.TrimEnd(':').ToLowerInvariant();
                if (PrefixOf("env", lower) || PrefixOf("setenv", lower))
                {
                    _parseType = ParseType.Env;
                }
                else if (PrefixOf("server", lower))
                {
                    _parseType = ParseType.Svr;
                }
                else if (PrefixOf("dbname", lower))
                {
                    _parseType = ParseType.Db;
                }
                else if (PrefixOf("resource", lower))
                {
                    _parseType = ParseType.Resource;
                }
                else if (PrefixOf("alias", lower))
                {
                    _parseType = ParseType.Alias;
                }
                else if (PrefixOf("include", lower))
                {
                    _parseType = ParseType.IFile;
                }
                else if (PrefixOf("if", lower))
                {
                    _parseType = ParseType.If;
                    saveifblock = false;
                }
                else if (PrefixOf("ifdef", lower))
                {
                    _parseType = ParseType.IfDef;
                    saveifblock = false;
                }
                else if (PrefixOf("else", lower))
                {
                    _parseType = ParseType.Else;
                }
                else if (PrefixOf("endif", lower))
                {
                    _parseType = ParseType.EndIf;
                }

                if (_parseType == ParseType.None)
                    NamError($"Invalid definition type '{lower}'");
            }
            else if (quoteOpen != '\0' && curword == "]")
            {
                NamError($"']' found, unclosed quotes in '{_value}'");
                quoteOpen = '\0';
                _parseType = ParseType.None;
            }
            else if (quoteOpen != '\0')
            {
                _value = string.Concat(_value, " ", curword);
                if (curword.Length > 0 && curword[^1] == quoteClose)
                {
                    _value = _value[..^1];
                    quoteOpen = quoteClose = '\0';
                    switch (_parseType)
                    {
                        case ParseType.Env or ParseType.Alias:
                            saveit = true;
                            break;
                        case ParseType.Svr:
                            svrsave = true;
                            break;
                        case ParseType.Db:
                            dbsave = true;
                            break;
                        case ParseType.Resource:
                            rssave = true;
                            break;
                        case ParseType.If or ParseType.IfDef:
                            saveifblock = true;
                            break;
                    }
                }
            }
            else if (_parseType is ParseType.If or ParseType.IfDef)
            {
                if (curword.StartsWith('\''))
                    quoteOpen = quoteClose = '\'';
                else if (curword.StartsWith('"'))
                    quoteOpen = quoteClose = '"';
                _value = curword;
                if (quoteOpen != '\0')
                {
                    _value = _value[1..];
                    if (_value.Length == 0)
                        Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                }
                else
                {
                    saveifblock = true;
                }

                if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                {
                    quoteOpen = quoteClose = '\0';
                    _value = _value[..^1];
                    saveifblock = true;
                }
            }
            else if (_parseType is ParseType.Env or ParseType.Alias)
            {
                if (_name is not null && _value is not null)
                {
                    saveit = true;
                }
                else if (_name is not null)
                {
                    if (curword.StartsWith('\''))
                        quoteOpen = quoteClose = '\'';
                    else if (curword.StartsWith('"'))
                        quoteOpen = quoteClose = '"';
                    _value = curword;
                    if (quoteOpen != '\0')
                    {
                        _value = _value[1..];
                        if (_value.Length == 0)
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
                        _value = _value[..^1];
                    }

                    if (saveit)
                        _value = Resolve(_value);
                }
                else
                {
                    _name = curword;
                }
            }
            else if (_parseType == ParseType.Svr)
            {
                if (curword == "[")
                {
                    svAttrs = [];
                    svLists = [];
                }
                else if (curword == "]")
                {
                    saveit = true;
                }
                else if (_name is not null)
                {
                    if (svrInput >= 0)
                    {
                        if (curword.StartsWith('\''))
                            quoteOpen = quoteClose = '\'';
                        else if (curword.StartsWith('"'))
                            quoteOpen = quoteClose = '"';
                        _value = curword;
                        if (quoteOpen != '\0')
                        {
                            _value = _value[1..];
                            if (_value.Length == 0)
                                Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                        }
                        else
                        {
                            svrsave = true;
                        }

                        if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                        {
                            quoteOpen = quoteClose = '\0';
                            _value = _value[..^1];
                            svrsave = true;
                        }

                        if (quoteOpen == '\0')
                            svrsave = true;
                        if (svrsave)
                            _value = Resolve(_value);
                    }
                    else if (curword.EndsWith(':'))
                    {
                        svrInput = AttrIndex(SvrAttrs, curword.TrimEnd(':'));
                        if (svrInput < 0)
                            Messages.Add($"Warning: {_fileName}: bad attribute '{curword.TrimEnd(':').ToLowerInvariant()}' for server '{_name}'");
                    }
                    else
                    {
                        Messages.Add($"Warning: {_fileName}: unexpected token '{curword}' for server '{_name}'");
                    }
                }
                else
                {
                    _name = curword;
                    if (!IsDbname(_name))
                        Messages.Add($"Error: Invalid server name '{_name}'");                    // ajErr, not namError: the program goes on
                }
            }
            else if (_parseType == ParseType.Db)
            {
                if (curword == "[")
                {
                    db = new EmbossDatabase(_name ?? "") { Defined = _shortName };
                }
                else if (curword == "]")
                {
                    saveit = true;
                }
                else if (_name is not null)
                {
                    if (dbInput >= 0)
                    {
                        if (curword.StartsWith('\''))
                            quoteOpen = quoteClose = '\'';
                        else if (curword.StartsWith('"'))
                            quoteOpen = quoteClose = '"';
                        _value = curword;
                        if (quoteOpen != '\0')
                        {
                            _value = _value[1..];
                            if (_value.Length == 0)
                                Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                        }
                        else
                        {
                            dbsave = true;
                        }

                        if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                        {
                            quoteOpen = quoteClose = '\0';
                            _value = _value[..^1];
                            dbsave = true;
                        }

                        if (quoteOpen == '\0')
                            dbsave = true;
                        if (dbsave)
                            _value = Resolve(_value);
                    }
                    else if (curword.EndsWith(':'))
                    {
                        dbInput = DbAttr(curword.TrimEnd(':'));
                        if (dbInput < 0)
                            Messages.Add($"Warning: {_fileName}: bad attribute '{curword.TrimEnd(':').ToLowerInvariant()}' for database '{_name}'");
                    }
                    else
                    {
                        Messages.Add($"Warning: {_fileName}: unexpected token '{curword}' for database '{_name}'");
                    }
                }
                else
                {
                    _name = curword;
                    if (!IsDbname(_name))
                        Messages.Add($"Error: Invalid database name '{_name}'");
                }
            }
            else if (_parseType == ParseType.Resource)
            {
                if (curword == "[")
                {
                    rsAttrs = [];
                }
                else if (curword == "]")
                {
                    saveit = true;
                }
                else if (_name is not null)
                {
                    if (curword.EndsWith(':'))
                    {
                        string lower = curword.ToLowerInvariant().TrimEnd(':');
                        rsInput = AttrIndex(RsAttrs, lower);
                        rsName = null;
                        if (rsInput < 0)
                        {
                            if (RsAttrField(rsAttrs, lower))
                                rsName = lower;
                            else
                                NamError($"Bad attribute '{lower}' for resource '{_name}'");
                        }
                    }
                    else if (rsName is not null || rsInput >= 0)
                    {
                        if (curword.StartsWith('\''))
                            quoteOpen = quoteClose = '\'';
                        else if (curword.StartsWith('"'))
                            quoteOpen = quoteClose = '"';
                        _value = curword;
                        if (quoteOpen != '\0')
                        {
                            _value = _value[1..];
                            if (_value.Length == 0)
                                Messages.Add($"Error: Bare quote {quoteOpen} found in namListParse");
                        }
                        else
                        {
                            rssave = true;
                        }

                        if (curword.Length > 0 && quoteClose != '\0' && curword[^1] == quoteClose)
                        {
                            quoteOpen = quoteClose = '\0';
                            _value = _value[..^1];
                            rssave = true;
                        }

                        if (quoteOpen == '\0')
                            rssave = true;
                    }

                    if (rssave)
                        _value = Resolve(_value ?? "");
                }
                else
                {
                    _name = curword;
                }
            }
            else if (_parseType == ParseType.IFile)
            {
                _parseType = ParseType.None;
                string saveword = curword;
                string includefn = Resolve(curword);
                if (!_included.Add(includefn))
                {
                    NamError(includefn == saveword ? $"{includefn} already read ... skipping" : $"{saveword} ({includefn}) already read ... skipping");
                }
                else
                {
                    string? incText = _files(includefn);
                    if (incText is null)
                    {
                        NamError(includefn == saveword ? $"Failed to open include file '{includefn}'" : $"Failed to open include file '{saveword}' ({includefn})");
                    }
                    else
                    {
                        string parentFile = _fileName;
                        _fileName = includefn;
                        bool ok = ProcessFile(incText);
                        _fileName = parentFile;
                        _line = linecount - 1;
                        if (!ok)
                            NamError(includefn == saveword ? $"Error(s) found in included file {includefn}" : $"Error(s) found in included file {saveword} ({includefn})");
                    }
                }
            }

            if (svrsave)
            {
                if (svrInput >= 0)
                {
                    string key = SvrAttrs[svrInput];
                    if (key is "field" or "special" or "edamdat" or "edamfmt" or "edamid" or "edamtpc" or "taxon")
                    {
                        if (!svLists!.TryGetValue(key, out var l))
                            svLists[key] = l = [];
                        l.Add(_value ?? "");
                    }
                    else
                    {
                        if (svAttrs!.TryGetValue(key, out string? prev))
                            Messages.Add($"Warning: duplicate server keyword value for '{key}' = '{prev}'");
                        svAttrs[key] = _value ?? "";
                    }
                }

                svrInput = -1;
                svrsave = false;
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
                        list.Add(_value ?? "");
                    }
                    else
                    {
                        if (db.Attributes.TryGetValue(attrName, out string? previous))
                            Messages.Add($"Warning: duplicate db keyword value for '{attrName}' = '{previous}'");
                        db.Attributes[attrName] = _value ?? "";
                    }
                }

                dbInput = -1;
                dbsave = false;
            }

            if (rssave)
            {
                string? key = rsInput >= 0 ? RsAttrs[rsInput] : rsName;
                if (key is not null && rsAttrs is not null)
                    rsAttrs[key] = _value ?? "";
                rsInput = -1;
                rsName = null;
                rssave = false;
            }

            if (_parseType == ParseType.EndIf)
            {
                if (_ifStack.Count == 0)
                {
                    NamError("endif not in if-block");
                }
                else
                {
                    _ifNow = _ifStack.Pop();
                }

                _parseType = ParseType.None;
            }
            else if (_parseType == ParseType.Else)
            {
                _ifNow = !_ifNow;
                _parseType = ParseType.None;
            }

            if (saveifblock)
            {
                if (_parseType is ParseType.If or ParseType.IfDef)
                {
                    string v = Resolve(_value ?? "");
                    _ifStack.Push(_ifNow);
                    if (_parseType == ParseType.If)
                    {
                        string test = GetValue(v, out string tv) ? tv : "";
                        _ifNow = test.Length > 0 && char.ToUpperInvariant(test[0]) is 'Y' or 'T' or '1';
                    }
                    else
                    {
                        _ifNow = GetValue(v, out _);
                    }

                    _parseType = ParseType.None;
                    _value = null;
                }

                saveifblock = false;
            }
            else if (saveit)
            {
                if (_ifNow)
                {
                    switch (_parseType)
                    {
                        case ParseType.Svr when _name is not null:
                            _servers.Add(_name);
                            break;
                        case ParseType.Db when _name is not null && db is not null:
                            _databases[_name] = db;
                            break;
                        case ParseType.Resource when _name is not null:
                            _resources[_name] = rsAttrs ?? [];
                            break;
                        case ParseType.Env when _name is not null:
                            _vars[_name] = _value ?? "";
                            break;
                        case ParseType.Alias when _name is not null:
                            _aliases[_name] = _value ?? "";
                            break;
                    }
                }

                _name = null;
                _value = null;
                db = null;
                svAttrs = null;
                rsAttrs = null;
                saveit = false;
                _parseType = ParseType.None;
            }
        }

        if (_parseType != ParseType.None)
        {
            NamError($"Unexpected end of file in {TypeName(_parseType)} definition");
            _parseType = ParseType.None;
        }

        _included.Clear();
        _fileName = savedFile;
        _ = svLists;
    }
}
