#pragma warning disable S1121, S3358, S3776, S107, S1854, S1066, S2692

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The kind of a query (<c>AJQUERY_*</c> of <c>ajquery.h</c>).</summary>
internal enum EmbossQueryType
{
    Unknown,
    All,
    Query,
    Entry,
}

/// <summary>One field of a query (<c>AjSQueryField</c>): the field name in lower case, the wildcard text and the link operator (<c>AJQLINK_*</c>: init, OR, AND, EOR, NOT, ELSE).</summary>
internal sealed class EmbossQueryField(string field, string wildquery, char link)
{
    public string Field { get; } = field.ToLowerInvariant();

    public string Wildquery { get; } = wildquery;

    /// <summary>I = init, | = OR, &amp; = AND, ^ = EOR, ! = NOT, E = ELSE.</summary>
    public char Link { get; } = link;
}

/// <summary>A query of a text input (<c>AjSQuery</c> of <c>ajquery.c</c>) in the parts a file or a local database query uses.</summary>
internal sealed class EmbossQuery
{
    public List<EmbossQueryField> Fields { get; } = [];

    /// <summary><c>DbFields</c>: the fields of the data type (the sequence fields for a file) or of the database definition.</summary>
    public string DbFields = "";

    public string Qlinks = "", SingleField = "", QryString = "", QryFields = "", Filename = "";

    /// <summary>The database of the query (<c>ajNamDbData</c>): its name as typed, the type (<c>Nucleotide</c> ...), the access method, the data directory, the files (a wildcard) and those excluded, the format of the entries.</summary>
    public string DbName = "", DbType = "", Method = "", Directory = "", Exclude = "", Formatstr = "";
    public bool Wild, CaseId, HasAcc = true;
    public long Fpos;
    public EmbossQueryType QueryType;
    public bool SetQuery;

    /// <summary><c>ajQueryClear</c>.</summary>
    public void Clear()
    {
        Fields.Clear();
        DbFields = Qlinks = SingleField = QryString = QryFields = Filename = "";
        DbName = DbType = Method = Directory = Exclude = Formatstr = "";
        Fpos = 0;
        QueryType = EmbossQueryType.Unknown;
        SetQuery = false;
        Wild = false;
        CaseId = false;
        HasAcc = true;
    }

    /// <summary><c>ajQueryKnownFieldC</c>: the field is one of the fields of the data type or database, or is <c>id</c>, or <c>acc</c> when the database has accession numbers.</summary>
    public bool KnownField(string field)
    {
        foreach (string token in DbFields.Split(['\t', ' ', ',', ';', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries))
        {
            if (string.Equals(token, field, StringComparison.OrdinalIgnoreCase))
                return true;
        }

        if (string.Equals(field, "id", StringComparison.OrdinalIgnoreCase))
            return true;
        return HasAcc && string.Equals(field, "acc", StringComparison.OrdinalIgnoreCase);
    }

    /// <summary><c>ajQueryAddFieldOrC</c> / <c>S</c>: the first field of a query is the initial one.</summary>
    public void AddOr(string field, string wildquery) => Fields.Add(new EmbossQueryField(field, wildquery, Fields.Count == 0 ? 'I' : '|'));

    public void AddAnd(string field, string wildquery) => Fields.Add(new EmbossQueryField(field, wildquery, '&'));

    public void AddNot(string field, string wildquery) => Fields.Add(new EmbossQueryField(field, wildquery, '!'));

    public void AddEor(string field, string wildquery) => Fields.Add(new EmbossQueryField(field, wildquery, '^'));

    public void AddElse(string field, string wildquery) => Fields.Add(new EmbossQueryField(field, wildquery, 'E'));

    /// <summary><c>ajQuerySetWild</c>: a wildcard in the id, or a query with a field that is tested after another.</summary>
    public bool SetWildFlag()
    {
        var wild = new Regex("[*?]", RegexOptions.CultureInvariant);
        foreach (var f in Fields)
        {
            if (f.Field == "id")
            {
                if (wild.IsMatch(f.Wildquery))
                {
                    Wild = true;
                    return true;
                }

                break;
            }
        }

        bool ret = false;
        int ifield = 0;
        foreach (var f in Fields)
        {
            if (ifield++ > 0 && f.Link != 'E')
                return true;
            if (f.Field == "id")
                continue;
            if (wild.IsMatch(f.Wildquery))
            {
                Wild = true;
                return true;
            }

            if (f.Field != "acc")
                ret = true;
        }

        return ret;
    }

    /// <summary><c>ajNamFileQuery</c>: the fields and link operators the sequence reader can test (<c>ajSeqinTypeGetFields</c>, <c>ajSeqinTypeGetQlinks</c>).</summary>
    public void SetFileQuery()
    {
        if (SetQuery)
            return;
        SetQuery = true;
        Qlinks = "|&!^=";
        DbFields = "id acc sv gi des key org";
    }
}

/// <summary>The text input of a program run (<c>AjSTextin</c>): the parts the USA processing of the sequence input reads and writes.</summary>
internal sealed class EmbossTextin
{
    /// <summary><c>Qry</c>: the USA as given.</summary>
    public string Qry = "";

    /// <summary><c>Formatstr</c>: the input format as named by <c>-sformat</c> or the prefix of the USA.</summary>
    public string Formatstr = "";

    /// <summary><c>QryFields</c>: <c>-iquery</c>.</summary>
    public string QryFields = "";

    public string Db = "";

    /// <summary><c>Fpos</c>: <c>-ioffset</c>.</summary>
    public long Fpos;

    public EmbossQuery Query = new();

    /// <summary><c>ajTextinClear</c>: back to the state of a new input except for the format number and the settings of the ACD.</summary>
    public void Clear()
    {
        Db = Qry = Formatstr = QryFields = "";
        Query.Clear();
    }
}

/// <summary>The statics of <c>ajquery.c</c> that keep their value from one query to the next (<c>queryChr</c>, <c>querySvr</c>) and what the parsers need from the environment.</summary>
internal sealed class EmbossQueryEnv
{
    /// <summary>The names of the files (the text of a file by name; null when there is none).</summary>
    public Func<string, string?> Files { get; init; } = _ => null;

    /// <summary><c>ajNamServer</c>: a defined server name.</summary>
    public Func<string, bool> IsServer { get; init; } = _ => false;

    /// <summary><c>queryDbfind</c> for a database name: a defined database (<c>ajNamAliasDatabase</c>).</summary>
    public Func<string, bool> IsDatabase { get; init; } = _ => false;

    /// <summary>The definitions of databases of the EMBOSS configuration (<c>.embossrc</c>); null = none.</summary>
    public EmbossNames? Names { get; init; }

    /// <summary>The names of the files of a directory (the file names without the directory); null when there is no such directory (<c>ajFilebuffNewPathWildExclude</c> lists the directory).</summary>
    public Func<string, IReadOnlyList<string>?> Directory { get; init; } = _ => null;

    internal string QueryChr = "";
    internal string QuerySvr = "";
}

/// <summary>The query-string parsers of <c>ajquery.c</c> (<c>ajQuerystrParseFormat</c>, <c>ajQuerystrParseListfile</c>, <c>ajQuerystrParseRange</c>, <c>ajQuerystrParseRead</c>) for the sequence input. A5-8c3h, B09 F255.</summary>
internal static class EmbossQueryParse
{
    private const RegexOptions Opt = RegexOptions.CultureInvariant;
    private static readonly Regex RegFmt = new("^([A-Za-z0-9-]*)::(.*)$", Opt | RegexOptions.Singleline);
    private static readonly Regex RegSvr = new("^([A-Za-z][A-Za-z0-9_.]+):", Opt | RegexOptions.Singleline);
    private static readonly Regex RegDb = new("^[A-Za-z][A-Za-z0-9_.]*[A-Za-z0-9]$", Opt);
    private static readonly Regex RegQryField = new("^([A-Za-z]+)?(([:]?)(([{]?)([^}]*)(}?)))?$", Opt);
    private static readonly Regex RegQryId = new("^([{]?)([^}]*)(}?)$", Opt);
    private static readonly Regex RegFieldId = new("^(([^{|]+[|])|[^{:]+)[:-]([^{:-]+)([:{]+)([^}]+)(}?)$", Opt);
    private static readonly Regex RegId = new("^(([^|{]+[|])|[^:{%]*[^:{%-])([:%-]?)(({?)([^}]*)(}?))?$", Opt);
    private static readonly Regex RegList = new("^(@|[Ll][Ii][Ss][Tt]:+)(.+)$", Opt);
    private static readonly Regex RegAsis = new("^[Aa][Ss][Ii][Ss]:+(.+)$", Opt);
    private static readonly Regex RegFtp = new("^([Ff][Tt][Pp]:/*.+)$", Opt);
    private static readonly Regex RegHttp = new("^([Hh][Tt][Tt][Pp]:/*.+)$", Opt);
    private static readonly Regex RegRange = new("(.*)[[](-?[0-9]*):(-?[0-9]*)(:([Rr])?)?[]]$", Opt);
    private static readonly Regex RegRefrange = new("(.*):([0-9]*)(:|\\.\\.)([0-9]*)(:([Rr]))?$", Opt);

    private static string Trim(string s) => s.Trim(' ', '\t', '\n');

    // strtol(text, &end, 10) with the check of ajStrToInt / ajStrToLong: the whole text is a number (leading white space and a sign allowed), the value fits a long
    private static bool StrToL(string s, out long value)
    {
        value = 0;
        if (s.Length == 0)
            return false;
        int i = 0;
        while (i < s.Length && s[i] is ' ' or '\t' or '\n' or '\v' or '\f' or '\r')
            i++;
        int start = i;
        if (i < s.Length && s[i] is '+' or '-')
            i++;
        int digits = i;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
            i++;
        if (i == digits || i != s.Length)
            return false;
        return long.TryParse(s.AsSpan(start), System.Globalization.NumberStyles.AllowLeadingSign, System.Globalization.CultureInfo.InvariantCulture, out value);
    }

    /// <summary><c>ajStrToInt</c>: the value of a decimal text, false when it has none (a value out of the range of int is clipped to it, one that does not fit a long gives false and 0).</summary>
    internal static bool ToInt(string s, out int value)
    {
        bool ok = StrToL(s, out long l);
        value = ok ? (int)Math.Clamp(l, int.MinValue, int.MaxValue) : 0;
        return ok;
    }

    /// <summary><c>ajStrToLong</c>.</summary>
    internal static bool ToLong(string s, out long value) => StrToL(s, out value);

    /// <summary>
    /// <c>ajQuerystrParseFormat</c>: strips a <c>format::</c> prefix and sets the format of the input. Returns the prefix (null when there is none, or for an asis / list USA); an unknown format is reported.
    /// </summary>
    public static string? ParseFormat(ref string qry, EmbossTextin textin, Func<string, bool> findFormat, EmbossFeatureReadContext context)
    {
        qry = Trim(qry);
        if (RegList.IsMatch(qry) || RegAsis.IsMatch(qry))
            return null;
        var m = RegFmt.Match(qry);
        if (!m.Success)
        {
            if (textin.Formatstr.Length > 0 && findFormat(textin.Formatstr))
                textin.Query.Formatstr = textin.Formatstr;
            return null;
        }

        string format = m.Groups[1].Value;
        qry = m.Groups[2].Value;
        if (format.Length == 0)
            format = "unknown";
        if (findFormat(format))
        {
            textin.Formatstr = format;
            textin.Query.Formatstr = format;
        }
        else
        {
            context.Error($"Unknown input format '{format}'\n");
        }
        return format;
    }

    /// <summary><c>ajQuerystrParseListfile</c>: strips <c>list::</c> / <c>@</c>; true when it was there.</summary>
    public static bool ParseListfile(ref string qry)
    {
        qry = Trim(qry);
        var m = RegList.Match(qry);
        if (!m.Success)
            return false;
        qry = m.Groups[2].Value;
        return true;
    }

    /// <summary><c>ajQuerystrParseRange</c>: strips <c>[begin:end:r]</c> or <c>:begin:end:r</c> / <c>:begin..end:r</c> and sets the range.</summary>
    public static bool ParseRange(ref string qry, ref int begin, ref int end, ref bool rev)
    {
        var m = RegRange.Match(qry);
        if (m.Success)
        {
            if (m.Groups[2].Length > 0)
                ToInt(m.Groups[2].Value, out begin);
            if (m.Groups[3].Length > 0)
                ToInt(m.Groups[3].Value, out end);
            if (m.Groups[5].Length > 0)
                rev = true;
            qry = m.Groups[1].Value;
            return true;
        }

        m = RegRefrange.Match(qry);
        if (m.Success)
        {
            if (m.Groups[2].Length > 0)
                ToInt(m.Groups[2].Value, out begin);
            if (m.Groups[4].Length > 0)
                ToInt(m.Groups[4].Value, out end);
            if (m.Groups[6].Length > 0)
                rev = true;
            qry = m.Groups[1].Value;
            return true;
        }

        return false;
    }

    /// <summary>The result of <see cref="ParseRead"/>: success, the text of the file the query reads from (null for a database query) and the offset in it.</summary>
    internal sealed record ReadResult(bool Ok, string? Text, string Filename, bool Asis, string? DbName = null)
    {
        /// <summary>The files of a database in reading order (one text for a file).</summary>
        public IReadOnlyList<string>? Texts { get; init; }
    }

    /// <summary>
    /// <c>ajQuerystrParseRead</c> for a file (or asis text) query: finds the file, sets the query fields from <c>file:id</c> / <c>file:field:query</c> / <c>file{...}</c> / <c>-iquery</c> and opens the text.
    /// </summary>
    public static ReadResult ParseRead(ref string pqry, EmbossTextin textin, Func<string, bool> findFormat, EmbossQueryEnv env, EmbossFeatureReadContext context)
    {
        var qry = textin.Query;
        qry.QryString = pqry;
        qry.SingleField = "";
        qry.QryFields = textin.QryFields;
        qry.Fpos = textin.Fpos;
        pqry = Trim(pqry);

        bool isSpecial = false, asis = false;
        var am = RegAsis.Match(pqry);
        if (am.Success)
        {
            asis = true;
            textin.Formatstr = "asis";
            qry.Filename = am.Groups[1].Value;
            isSpecial = true;
        }

        var fm = RegFtp.Match(pqry);
        bool ftp = fm.Success;
        if (ftp)
        {
            qry.Filename = fm.Groups[1].Value;
            isSpecial = true;
        }

        var hm = RegHttp.Match(pqry);
        bool http = hm.Success;
        if (http)
        {
            qry.Filename = hm.Groups[1].Value;
            isSpecial = true;
        }

        bool svrstat = false;
        var sm = RegSvr.Match(pqry);
        if (sm.Success)
        {
            svrstat = true;
            env.QuerySvr = sm.Groups[1].Value;
            if (env.IsServer(env.QuerySvr))
                pqry = pqry[sm.Length..];
            else
            {
                env.QuerySvr = "";
                svrstat = false;
            }
        }

        bool regstat = false, dbstat = false;
        string queryDb = pqry, braceopen = "", braceclose = "";
        if (!isSpecial)
        {
            var m = RegFieldId.Match(pqry);
            regstat = m.Success;
            if (regstat)
            {
                queryDb = m.Groups[1].Value;
                if (DbFind(ref queryDb, env, out dbstat))
                {
                    qry.SingleField = m.Groups[3].Value;
                    braceopen = m.Groups[4].Value;
                    qry.QryFields = m.Groups[5].Value;
                    braceclose = m.Groups[6].Value;
                    braceopen = braceopen.TrimStart(':');
                }
                else
                    regstat = false;
            }

            if (!regstat)
            {
                m = RegId.Match(pqry);
                regstat = m.Success;
                if (regstat)
                {
                    queryDb = m.Groups[1].Value;
                    env.QueryChr = m.Groups[3].Value;
                    braceopen = m.Groups[5].Value;
                    qry.QryFields = m.Groups[6].Value;
                    braceclose = m.Groups[7].Value;
                    regstat = DbFind(ref queryDb, env, out dbstat);
                }
            }
        }

        if (svrstat && !dbstat)
        {
            context.Error($"Server '{env.QuerySvr}:' but no database name in '{pqry}'");
            return new ReadResult(false, null, "", false);
        }

        if (regstat && !dbstat && env.QueryChr == "%")
        {
            ToLong(qry.QryFields, out qry.Fpos);
            regstat = false;
        }

        // no query found, but test -iquery input
        if (!regstat && textin.QryFields.Length > 0)
        {
            var m = RegQryField.Match(textin.QryFields);
            regstat = m.Success;
            if (regstat)
            {
                qry.SingleField = m.Groups[1].Value;
                env.QueryChr = m.Groups[3].Value;
                braceopen = m.Groups[5].Value;
                qry.QryFields = m.Groups[6].Value;
                braceclose = m.Groups[7].Value;
            }
        }

        if (!regstat && textin.QryFields.Length > 0)
        {
            var m = RegQryId.Match(textin.QryFields);
            regstat = m.Success;
            if (regstat)
            {
                braceopen = m.Groups[1].Value;
                qry.QryFields = m.Groups[2].Value;
                braceclose = m.Groups[3].Value;
            }
        }

        if (dbstat)
        {
            // we have a database
            qry.DbName = queryDb;
            if (env.QuerySvr.Length > 0)
                throw new NotSupportedException("A server (server:database) is not read by this port: the data come from a network service.");
            if (!NamDbData(qry, env, context, findFormat))
            {
                qry.QryFields = "";
                return new ReadResult(false, null, "", false);
            }

            if (findFormat(qry.Formatstr))
            {
                textin.Formatstr = qry.Formatstr;
            }
            else
            {
                context.Error($"unknown format '{qry.Formatstr}' for '{pqry}' type '{qry.DbType}'\n");
                return new ReadResult(false, null, "", false);
            }
        }

        qry.QryFields = qry.QryFields.Trim('\t', ' ', '\n', '\r');
        bool inbraces = false;
        if (braceopen.Length > 0)
        {
            if (braceclose.Length > 0)
                inbraces = true;
            else
            {
                context.Error($"Bad query syntax: unclosed braces '{pqry}'");
                qry.QryFields = "";
                return new ReadResult(false, null, "", false);
            }
        }
        else if (braceclose.Length > 0)
        {
            context.Error($"Bad query syntax: unopened braces '{pqry}'");
            qry.QryFields = "";
            return new ReadResult(false, null, "", false);
        }

        if (regstat)
        {
            if (qry.QryFields.Length > 0)
            {
                if (qry.QryFields[0] == '!')
                    qry.QryFields = "*" + qry.QryFields;
                string lastoper = "|";
                int iquery = 0;
                string idstr = "";
                if (qry.SingleField.Length > 0)
                {
                    qry.QueryType = EmbossQueryType.Query;
                    SetQuery(qry, dbstat, env, context);
                    if (!qry.KnownField(qry.SingleField))
                    {
                        context.Error($"Query '{pqry}' query field '{qry.SingleField}' not defined for {Target(qry, dbstat)}");
                        return new ReadResult(false, null, "", false);
                    }

                    string delims = inbraces ? "|&!^ \t\n\r" : "|&!^\t\n\r";
                    var tok = new DelimTokens(qry.QryFields, delims);
                    while (tok.Next(out idstr, out string oper))
                    {
                        if (iquery > 0 && (qry.Qlinks.Length == 0 || qry.Qlinks.IndexOfAny(lastoper.ToCharArray()) < 0))
                        {
                            context.Error($"Query link operator '{lastoper}' not supported by access method '{qry.Method}' for {Target(qry, dbstat)} allowed links are '{qry.Qlinks}'");
                            return new ReadResult(false, null, "", false);
                        }

                        iquery++;
                        AddField(qry, lastoper, qry.SingleField, idstr, context);
                        lastoper = oper.Trim('\t', ' ', '\n', '\r');
                    }
                }
                else if (qry.QryFields.IndexOf(':') > 0)
                {
                    qry.QueryType = EmbossQueryType.Query;
                    SetQuery(qry, dbstat, env, context);
                    var tok = new DelimTokens(qry.QryFields, "|&^!");
                    while (tok.Next(out string qrystr, out string oper))
                    {
                        if (oper.Length > 0 && (qry.Qlinks.Length == 0 || qry.Qlinks.IndexOfAny(oper.ToCharArray()) < 0))
                        {
                            context.Error($"Query link operator '{oper}' not supported by access method '{qry.Method}' for {Target(qry, dbstat)} allowed links are '{qry.Qlinks}'");
                            return new ReadResult(false, null, "", false);
                        }

                        qrystr = qrystr.Trim('\t', ' ', '\n', '\r');
                        int i = qrystr.IndexOf(':');
                        if (i > 0)
                        {
                            qry.SingleField = qrystr[..i];
                            idstr = qrystr[(i + 1)..];
                        }

                        if (!qry.KnownField(qry.SingleField))
                        {
                            context.Error($"Query '{pqry}' query field '{qry.SingleField}' not defined for {Target(qry, dbstat)}");
                            qry.QryFields = "";
                            return new ReadResult(false, null, "", false);
                        }

                        if (iquery > 0 && (qry.Qlinks.Length == 0 || qry.Qlinks.IndexOfAny(lastoper.ToCharArray()) < 0))
                        {
                            context.Error($"Query link operator '{lastoper}' not supported by access method '{qry.Method}' for {Target(qry, dbstat)} allowed links are '{qry.Qlinks}'");
                            return new ReadResult(false, null, "", false);
                        }

                        iquery++;
                        AddField(qry, lastoper, qry.SingleField, idstr, context);
                        lastoper = oper.Trim('\t', ' ', '\n', '\r');
                    }
                }
                else if (qry.QryFields.IndexOfAny([' ', '\t', '\n', '\r', ',', ';', '|']) >= 0)
                {
                    qry.QueryType = EmbossQueryType.Query;
                    SetQuery(qry, dbstat, env, context);
                    foreach (string id in qry.QryFields.Split([' ', '\t', '\n', '\r', ',', ';', '|'], StringSplitOptions.RemoveEmptyEntries))
                        AddIdFields(qry, id);
                }
                else
                {
                    qry.QueryType = EmbossQueryType.Entry;
                    SetQuery(qry, dbstat, env, context);
                    AddIdFields(qry, qry.QryFields);
                }
            }
            else
            {
                qry.QueryType = EmbossQueryType.All;
            }
        }
        else
        {
            qry.QueryType = EmbossQueryType.All;
        }

        if (qry.SetWildFlag() && qry.QueryType == EmbossQueryType.Entry)
            qry.QueryType = EmbossQueryType.Query;

        if (isSpecial)
        {
            if (asis)
                return new ReadResult(true, qry.Filename, "asis", true);
            // ftp and http: the text is read from a server; this port has no network (the files are given by name)
            throw new NotSupportedException("The input is read from a URL (ftp:// or http://): the data are not part of this port.");
        }

        if (dbstat)
        {
            if (qry.Method != "direct")
                throw new NotSupportedException($"The database access method '{qry.Method}' is not implemented by this port (only 'direct', which reads the files of the database as they are).");
            return TextAccessDirect(qry, env, context);
        }

        // a file name that ends with '|' is a command the library runs and reads the output of
        if (queryDb.EndsWith('|'))
            throw new NotSupportedException($"The input name '{queryDb}' ends with '|': the library runs it as a command, which this port does not do.");

        // a file name with a possible query
        string? text = env.Files(queryDb);
        if (text is not null)
        {
            qry.Filename = queryDb;
            if (qry.Fpos != 0)
            {
                // ajTextinAccessOffset: the file is opened and positioned at the byte offset
                qry.QryFields = "";
                return new ReadResult(true, text, queryDb, false);
            }

            return new ReadResult(true, text, queryDb, false);
        }

        context.Error($"Failed to open filename '{queryDb}'");
        return new ReadResult(false, null, "", false);
    }

    // textAccessDirect: the files of the directory that match the file name (a wildcard) and are not excluded, read one after the other (the library lists the directory, so the order is the one of the file system;
    // this port sorts the names)
    private static ReadResult TextAccessDirect(EmbossQuery qry, EmbossQueryEnv env, EmbossFeatureReadContext context)
    {
        if (qry.Filename.Length == 0)
        {
            context.Error("DIRECT access: filename not specified");
            return new ReadResult(false, null, "", false);
        }

        string dir = qry.Directory.Length > 0 ? qry.Directory : ".";
        if (!dir.EndsWith('/'))
            dir += "/";
        var names = env.Directory(dir.TrimEnd('/'));
        if (names is null)
            return new ReadResult(false, null, "", false);
        var accepted = new List<string>();
        foreach (string n in names.Order(StringComparer.Ordinal))
        {
            if (n is "." or "..")
                continue;
            if (FilenameTestExclude(n, qry.Exclude, qry.Filename))
                accepted.Add(n);
        }

        var texts = new List<string>();
        foreach (string n in accepted)
        {
            string? t = env.Files(dir == "./" ? n : dir + n);
            if (t is not null)
                texts.Add(t);
        }

        if (texts.Count == 0)
            return new ReadResult(false, null, "", false);
        return new ReadResult(true, texts[0], qry.Filename, false, qry.DbName) { Texts = texts };
    }

    // ajFilenameTestExclude: the name matches a pattern of the include list and none of the exclude list (patterns separated by blanks, commas, semicolons)
    private static bool FilenameTestExclude(string filename, string exclude, string include)
    {
        bool ret = false;
        char[] sep = [' ', '\t', ',', ';', '\n'];
        foreach (string token in include.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.Match(filename, token))
                ret = true;
        }

        foreach (string token in exclude.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.Match(filename, token))
                ret = false;
        }

        return ret;
    }

    private static readonly string[] DbTypes = ["N", "P", "Nucleotide", "Protein", "Sequence", "Nucfeatures", "Protfeatures", "Features", "Assembly", "Obo", "Refseq", "Resource", "Taxonomy", "Html", "Text", "Url", "Variation", "Xml", "Unknown"];

    // the data type of a database type name (namDbTypes): only the sequence types matter for a sequence input
    private static bool IsSequenceType(string type) => type.ToLowerInvariant() is "n" or "p" or "nucleotide" or "protein" or "sequence";

    // namVarResolve: a value that starts with $name gets the value of the variable (again and again)
    private static string VarResolve(EmbossQueryEnv env, string value)
    {
        var re = new Regex(@"^\$([a-zA-Z0-9_.]+)", RegexOptions.CultureInvariant);
        for (int guard = 0; guard < 100; guard++)
        {
            var m = re.Match(value);
            if (!m.Success)
                break;
            string v = env.Names?.Variable(m.Groups[1].Value) ?? "";
            value = v + value[m.Length..];
        }

        return value;
    }

    // namDbSetAttrStrC: an attribute of the definition, with its default; false when the database has no value and the default is empty
    private static bool DbAttr(EmbossQueryEnv env, EmbossDatabase db, string attribute, ref string target)
    {
        string value;
        if (db.Lists.TryGetValue(attribute, out var list))
        {
            var sb = new StringBuilder();
            char delim = ' ';
            foreach (string item in list)
            {
                foreach (string word in item.Split(' ', StringSplitOptions.RemoveEmptyEntries))
                {
                    if (word == "!")
                        break;
                    if (sb.Length > 0)
                    {
                        sb.Append(delim);
                        delim = ' ';
                    }

                    sb.Append(word);
                }

                delim = ';';
            }

            value = sb.ToString();
        }
        else if (db.Attributes.TryGetValue(attribute, out string? v))
        {
            if (v.Length == 0)
                return false;
            value = v;
        }
        else
        {
            return false;
        }

        target = VarResolve(env, value);
        return true;
    }

    // ajNamDbData for the data type sequence
    private static bool NamDbData(EmbossQuery qry, EmbossQueryEnv env, EmbossFeatureReadContext context, Func<string, bool> findFormat)
    {
        var db = env.Names?.Database(qry.DbName);
        if (db is null)
            return false;
        string liststr = "";
        if (DbAttr(env, db, "type", ref liststr))
        {
            qry.DbType = "";
            foreach (string token0 in liststr.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries))
            {
                string token = token0;
                if (token.Equals("n", StringComparison.OrdinalIgnoreCase))
                    token = "Nucleotide";
                else if (token.Equals("p", StringComparison.OrdinalIgnoreCase))
                    token = "Protein";
                else if (token.Equals("prot", StringComparison.OrdinalIgnoreCase))
                    token = "Protein";
                if (!DbTypes.Contains(token, StringComparer.OrdinalIgnoreCase))
                {
                    context.Error($"Database {qry.DbName} type '{token}' unknown");
                    continue;
                }

                if (IsSequenceType(token))
                {
                    qry.DbType = token;
                    break;
                }
            }

            if (qry.DbType.Length == 0)
            {
                context.Error($"Database {qry.DbName} type(s) '{liststr}' do not match query type 'sequence'");
                return false;
            }
        }

        if (qry.Formatstr.Length > 0)
        {
            if (!FormatTest(qry.Formatstr, qry.DbType))
            {
                context.Error($"Query format '{qry.Formatstr}' does not support datatype '{qry.DbType}'");
                return false;
            }
        }
        else if (DbAttr(env, db, "format", ref liststr))
        {
            bool ok = false;
            foreach (string token in liststr.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries))
            {
                if (FormatTest(token, qry.DbType))
                {
                    ok = true;
                    qry.Formatstr = token;
                }
            }

            if (!ok)
            {
                context.Error($"Database {qry.DbName} format(s) '{liststr}' do not support datatype '{qry.DbType}'");
                return false;
            }
        }

        string tmp = "";
        DbAttr(env, db, "method", ref qry.Method);
        string indexdir = "";
        DbAttr(env, db, "directory", ref indexdir);
        DbAttr(env, db, "indexdirectory", ref qry.Directory);
        DbAttr(env, db, "directory", ref qry.Directory);
        DbAttr(env, db, "exclude", ref qry.Exclude);
        DbAttr(env, db, "filename", ref qry.Filename);
        if (!DbAttr(env, db, "field", ref qry.DbFields))
            DbAttr(env, db, "fields", ref qry.DbFields);
        if (DbAttr(env, db, "caseidmatch", ref tmp))
            qry.CaseId = ToBool(tmp);
        tmp = "";
        if (DbAttr(env, db, "hasaccession", ref tmp))
            qry.HasAcc = ToBool(tmp);
        return true;
    }

    // ajStrToBool: Y, YES, T, TRUE, 1 (any case)
    private static bool ToBool(string s) => s.Length > 0 && (char.ToUpperInvariant(s[0]) is 'Y' or 'T' or '1');

    // namInformatTest: the format is one of the formats of the sequence reader
    private static bool FormatTest(string format, string dbtype) => DbTypes.Contains(dbtype, StringComparer.OrdinalIgnoreCase) && IsSequenceType(dbtype) && EmbossSeqin.FormatNameKnown(format);

    // ajNamDbQuery / ajNamFileQuery for a query that has fields
    private static void SetQuery(EmbossQuery qry, bool dbstat, EmbossQueryEnv env, EmbossFeatureReadContext context)
    {
        if (!dbstat)
        {
            qry.SetFileQuery();
            return;
        }

        if (qry.SetQuery)
            return;
        qry.SetQuery = true;
        var db = env.Names!.Database(qry.DbName);
        if (db is null)
            return;
        string tmp = "";
        if (qry.DbType.Length == 0 && !NamDbData(qry, env, context, _ => true))
            return;
        if (qry.QueryType == EmbossQueryType.All || (qry.QueryType == EmbossQueryType.Unknown && qry.Fields.Count == 0))
        {
            DbAttr(env, db, "methodall", ref qry.Method);
            DbAttr(env, db, "formatall", ref qry.Formatstr);
            qry.QueryType = EmbossQueryType.All;
        }
        else
        {
            DbAttr(env, db, "methodquery", ref qry.Method);
            DbAttr(env, db, "formatquery", ref qry.Formatstr);
            if (qry.QueryType != EmbossQueryType.Query && !qry.SetWildFlag())
            {
                DbAttr(env, db, "methodentry", ref qry.Method);
                DbAttr(env, db, "formatentry", ref qry.Formatstr);
                qry.QueryType = EmbossQueryType.Entry;
            }
            else
            {
                qry.QueryType = EmbossQueryType.Query;
            }
        }

        _ = tmp;
        if (qry.Formatstr.Length == 0)
        {
            context.Error($"No format defined for database '{qry.DbName}'");
            return;
        }

        if (qry.Method.Length == 0)
        {
            context.Error($"No access method for database '{qry.DbName}'");
            return;
        }

        if (qry.Method != "direct")
            throw new NotSupportedException($"The database access method '{qry.Method}' is not implemented by this port (only 'direct', which reads the files of the database as they are).");
        qry.Qlinks = "|&!^=";
    }

    // the end of the messages about a query: "database 'name'" or "datatype 'sequence'"
    private static string Target(EmbossQuery qry, bool dbstat) => dbstat ? $"database '{qry.DbName}'" : "datatype 'sequence'";

    private static void AddIdFields(EmbossQuery qry, string id)
    {
        qry.AddOr("id", id);
        if (qry.HasAcc && qry.KnownField("acc"))
            qry.AddElse("acc", id);
        if (qry.KnownField("sv"))
            qry.AddElse("sv", id);
    }

    private static void AddField(EmbossQuery qry, string oper, string field, string wild, EmbossFeatureReadContext context)
    {
        switch (oper.Length > 0 ? oper[0] : '\0')
        {
            case '|':
                qry.AddOr(field, wild);
                break;
            case '&':
                qry.AddAnd(field, wild);
                break;
            case '!':
                qry.AddNot(field, wild);
                break;
            case '^':
                qry.AddEor(field, wild);
                break;
            default:
                context.Error($"bad query syntax: unknown operator '{oper}'");
                break;
        }
    }

    // queryDbfind for the databases of this port (none are defined): a name that is a file; the library then looks at the data resource catalogue (DRCAT), which this port has not
    private static bool DbFind(ref string db, EmbossQueryEnv env, out bool dbstat)
    {
        dbstat = false;
        bool ret = false;
        if (RegDb.IsMatch(db))
        {
            dbstat = true;
            string alias = db;
            if (env.Names is { } names && names.AliasDatabase(ref alias))
            {
                db = alias;
            }
            else
            {
                if (env.Files(db) is not null)
                {
                    ret = true;
                    dbstat = false;
                }
                else
                {
                    dbstat = false;
                }
            }
        }

        if (dbstat || ret)
            return true;
        if (env.Files(db) is not null)
        {
            ret = true;
            dbstat = false;
        }

        return ret;
    }

    /// <summary><c>ajStrTokenAssignC</c> + <c>ajStrTokenNextParseDelimiters</c>.</summary>
    private sealed class DelimTokens
    {
        private readonly string _s, _delims;
        private int _pos;

        public DelimTokens(string s, string delims)
        {
            _s = s;
            _delims = delims;
            while (_pos < _s.Length && _delims.Contains(_s[_pos], StringComparison.Ordinal))
                _pos++;
        }

        public bool Next(out string token, out string delim)
        {
            token = delim = "";
            if (_pos >= _s.Length)
            {
                _pos = 0;
                while (_pos < _s.Length && _delims.Contains(_s[_pos], StringComparison.Ordinal))
                    _pos++;
                _pos = _s.Length + 1;
                return false;
            }

            int start = _pos;
            while (_pos < _s.Length && !_delims.Contains(_s[_pos], StringComparison.Ordinal))
                _pos++;
            token = _s[start.._pos];
            start = _pos;
            while (_pos < _s.Length && _delims.Contains(_s[_pos], StringComparison.Ordinal))
                _pos++;
            delim = _s[start.._pos];
            return true;
        }
    }
}
