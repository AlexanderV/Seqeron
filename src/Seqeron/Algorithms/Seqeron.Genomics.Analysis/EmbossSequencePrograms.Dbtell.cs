#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The database information programs of EMBOSS 6.6.0 on the standard configuration (<c>emboss.standard</c> of the installation and the data resource catalogue; B09 WP67, F309): <c>dbtell</c> and <c>showdb</c>.</summary>
internal static partial class EmbossSequencePrograms
{
    // the single-value attributes of ajNamDbDetails, in the order of the attribute table (namDbAttrs)
    private static readonly string[] DbDetailAttrs =
    [
        "format", "method", "type", "accession", "app", "appall", "appentry", "appquery", "comment", "dbalias", "directory", "exclude", "fields", "filters", "filename", "formatall", "formatentry", "formatquery", "httpversion",
        "identifier", "indexdirectory", "methodall", "methodentry", "methodquery", "namespace", "organisms", "proxy", "query", "release", "return", "sequence", "serverversion", "url", "view",
    ];

    // every attribute of a definition as dbtell prints it (namDbAttrs order; the list attributes one line per value)
    private static readonly (string Name, bool List)[] DbPrintAttrs =
    [
        ("format", false), ("method", false), ("type", false), ("accession", false), ("app", false), ("appall", false), ("appentry", false), ("appquery", false), ("caseidmatch", false), ("comment", false), ("dbalias", false),
        ("directory", false), ("example", true), ("exclude", false), ("field", true), ("fields", false), ("filters", false), ("filename", false), ("formatall", false), ("formatentry", false), ("formatquery", false),
        ("hasaccession", false), ("httpversion", false), ("identifier", false), ("indexdirectory", false), ("methodall", false), ("methodentry", false), ("methodquery", false), ("namespace", false), ("organisms", false),
        ("proxy", false), ("query", false), ("release", false), ("return", false), ("sequence", false), ("serverversion", false), ("special", true), ("url", false), ("view", false), ("edamdat", true), ("edamfmt", true),
        ("edamid", true), ("edamtpc", true), ("taxon", true),
    ];

    private sealed record DbDetails(string Type, bool Id, bool Qry, bool All, string Comment, string Release, string Methods, string Defined);

    // ajNamDbDetails
    private static DbDetails? DatabaseDetails(EmbossNames names, string name)
    {
        var db = names.Database(name);
        if (db is null)
            return null;
        bool id = false, qry = false, all = false;
        string type = "", comment = "", release = "";
        var methods = new StringBuilder();
        string dbtype = db.Get("type") ?? "";
        foreach (string attr in DbDetailAttrs)
        {
            string? value = db.Get(attr);
            if (string.IsNullOrEmpty(value))
                continue;
            switch (attr)
            {
                case "type":
                    type = value.Equals("N", StringComparison.OrdinalIgnoreCase) ? "Nucleotide" : value.Equals("P", StringComparison.OrdinalIgnoreCase) ? "Protein" : char.ToUpperInvariant(value[0]) + value[1..].ToLowerInvariant();
                    break;
                case "method":
                {
                    int s = EmbossAccessMethods.Scope(value, dbtype);
                    id |= EmbossAccessMethods.HasEntry(s);
                    qry |= EmbossAccessMethods.HasQuery(s);
                    all |= EmbossAccessMethods.HasAll(s);
                    methods.Append(value);
                    break;
                }

                case "methodentry":
                    id |= EmbossAccessMethods.HasEntry(EmbossAccessMethods.Scope(value, dbtype));
                    methods.Append(methods.Length > 0 ? "," : "").Append(value).Append("(id)");
                    break;
                case "methodquery":
                {
                    int s = EmbossAccessMethods.Scope(value, dbtype);
                    id |= EmbossAccessMethods.HasEntry(s);
                    qry |= EmbossAccessMethods.HasQuery(s);
                    methods.Append(methods.Length > 0 ? "," : "").Append(value).Append("(qry)");
                    break;
                }

                case "methodall":
                    all |= EmbossAccessMethods.HasAll(EmbossAccessMethods.Scope(value, dbtype));
                    methods.Append(methods.Length > 0 ? "," : "").Append(value).Append("(all)");
                    break;
                case "comment":
                    comment = value;
                    break;
                case "release":
                    release = value;
                    break;
            }
        }

        if (type.Length == 0)
            throw new NotSupportedException("A database definition without a type (a user definition) is not read by this port.");
        return new DbDetails(type, id, qry, all, comment, release, methods.ToString(), db.Defined);
    }

    private static void RequireStandardConfiguration(IReadOnlyDictionary<string, string> files)
    {
        if (files.ContainsKey(".embossrc") || files.ContainsKey("emboss.default"))
            throw new NotSupportedException("The databases of a user configuration file (.embossrc) are not read by dbtell / showdb of this port: it reports the standard definitions of the installation (the validation and the access method tables of the user-defined types are not ported).");
    }

    // an output file: stdout, or a named file; the directory qualifier is put in front of a name without one
    private static EmbossRunResult WriteOutfile(EmbossAcdArgs a, string name, string text)
    {
        if (name.Length == 0)
            throw a.Die("outfile", "Error: Output file is required");
        string odir = EmbossOutPath.FixDirectory(a.RawAssoc("outfile", "odirectory"));
        if (odir.Length > 0)
            throw new NotSupportedException("The output directory qualifier of dbtell / showdb is not written by this port.");
        if (name == "stdout")
            return new EmbossRunResult(text, "", 0, new Dictionary<string, string>());
        if (name.Contains('/', StringComparison.Ordinal))
            throw new NotSupportedException("An output file in another directory is not written by this port.");
        return new EmbossRunResult("", "", 0, new Dictionary<string, string> { [name] = text });
    }

    /// <summary>
    /// <c>dbtell</c> of EMBOSS 6.6.0 (<c>emboss/dbtell.c</c>): the definition of a database as the configuration holds it (a database defined in <c>emboss.standard</c>: header with the file and the access levels, the attributes; with <c>-full</c> the aliases), or, when
    /// there is no such definition, the line that says a data resource of the catalogue <c>DRCAT</c> has that name. A server (<c>-server</c>) and the databases of a user configuration file are not read.
    /// </summary>
    internal static EmbossRunResult Dbtell(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("dbtell");
        try
        {
            RequireStandardConfiguration(files);
            string dbname = a.String("database", "");
            string svrname = a.String("server", "");
            bool verbose = a.Bool("full", false);
            string outfile = a.Raw("outfile") ?? "stdout";
            if (svrname.Length > 0)
                throw new NotSupportedException("The databases of a server (-server) are not read by this port (the server definitions list their databases from the network).");
            var names = EmbossNames.ParseWithStandard(null, null);
            string truedb = dbname;
            names.AliasDatabase(ref truedb);
            truedb = truedb.ToLowerInvariant();
            bool isalias = dbname != truedb;
            var sb = new StringBuilder();
            string stderr = "";
            var details = DatabaseDetails(names, truedb);
            if (details is not null)
            {
                var db = names.Database(truedb)!;
                sb.Append(isalias ? $"# {dbname} is an alias for {truedb} defined in {details.Defined}\n" : $"# {truedb} is defined in {details.Defined}\n");
                sb.Append($"# access levels id: {(details.Id ? "Yes" : "No")} query: {(details.Qry ? "Yes" : "No")} all: {(details.All ? "Yes" : "No")}\n\n");
                sb.Append($"DBNAME {truedb} [\n");
                foreach (var (name, list) in DbPrintAttrs)
                {
                    IEnumerable<string> values = list ? db.Lists.GetValueOrDefault(name) ?? [] : db.Get(name) is { } v ? [v] : [];
                    foreach (string value in values)
                        sb.Append($"   {name}:{new string(' ', Math.Max(0, 15 - name.Length))}\"{value}\"\n");
                }

                sb.Append("]\n");
            }
            else
            {
                var messages = new List<string>();
                if (EmbossObo.Query($"drcat:{dbname}", "resource", messages).Count > 0)
                    sb.Append($"DBNAME {dbname} [\n   comment: \"defined in DRCAT\"\n]\n");
                stderr = EmbossRunResult.Join(messages);
            }

            _ = verbose;
            return WriteOutfile(a, outfile, sb.ToString()) with { Stderr = stderr };
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }
    }

    // showdbGetFields / showdbGetTaxons / showdbGetAliases / showdbGetExamples: the text of a column and the count
    private static (string Text, int Count) ShowdbList(EmbossDatabase db, string attribute)
    {
        // namDbSetAttrStrC of a list: the words of every value up to a "!" (the first value's words separated by a space, the values by ';')
        var sb = new StringBuilder();
        char delim = ' ';
        foreach (string value in db.Lists.GetValueOrDefault(attribute) ?? [])
        {
            foreach (string word in value.Split(' ', StringSplitOptions.RemoveEmptyEntries))
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

        return (sb.ToString(), 0);
    }

    private static (string Text, int Count) ShowdbFields(EmbossDatabase db)
    {
        string fields = ShowdbList(db, "field").Text;
        if (fields.Length == 0)
            fields = db.Get("fields") ?? "";
        bool hasAcc = db.Get("hasaccession") is not { Length: > 0 } h || h[0] is 'Y' or 'y' or 'T' or 't' or '1';
        var defined = fields.Split(['\t', ' ', ',', ';', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
        int count = 0;
        bool hasid = false, hasacc = false;
        if (!defined.Contains("id", StringComparer.OrdinalIgnoreCase))
        {
            hasid = true;
            count++;
        }

        if (hasAcc && !defined.Contains("acc", StringComparer.OrdinalIgnoreCase))
        {
            hasacc = true;
            count++;
        }

        if (fields.Length == 0)
        {
            if (hasid)
                fields = hasacc ? "id,acc" : "id";
            else
                fields = "-";
        }
        else
        {
            fields = fields.Replace(' ', ',').Replace(';', ',');
            count += 1 + fields.Count(c => c == ',');
        }

        return (fields, count);
    }

    private static (string Text, int Count) ShowdbTaxons(EmbossDatabase db)
    {
        string taxons = ShowdbList(db, "taxon").Text;
        if (taxons.Length == 0)
            return ("-", 0);
        return (taxons.Replace(' ', '_'), 1 + taxons.Count(c => c == ';'));
    }

    private static (string Text, int Count) ShowdbExamples(EmbossDatabase db)
    {
        string examples = ShowdbList(db, "example").Text;
        if (examples.Length == 0)
            return ("-", 0);
        string joined = examples.Replace(' ', ';');
        return (joined, 1 + joined.Count(c => c == ';'));
    }

    /// <summary>
    /// <c>showdb</c> of EMBOSS 6.6.0 (<c>emboss/showdb.c</c>): the table of the configured databases (the standard definitions of the installation: the ontologies, the taxonomy, the data resource catalogue) in groups by type, sorted by definition file and name, with the columns the
    /// qualifiers select (text or <c>-html</c>); one database with <c>-database</c> (a missing one is fatal). The databases of a user configuration file, the aliases they define and the validation of their definitions are not read.
    /// </summary>
    internal static EmbossRunResult Showdb(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("showdb");
        try
        {
            RequireStandardConfiguration(files);
            string onedb = a.String("database", "");
            string outfile = a.Raw("outfile") ?? "stdout";
            bool html = a.Bool("html", false);
            bool showprotein = a.Bool("protein", true), shownucleic = a.Bool("nucleic", true), showobo = a.Bool("obo", true), showfeature = a.Bool("feature", true), showtext = a.Bool("text", true), showxml = a.Bool("xml", true);
            bool showtax = a.Bool("taxonomy", true), showstandard = a.Bool("standard", true), showuser = a.Bool("user", true), showinclude = a.Bool("include", true), showresource = a.Bool("resource", true), showsequence = a.Bool("sequence", true);
            bool full = a.Bool("full", false), only = a.Bool("only", false);
            bool doheader = a.Bool("heading", !only), dotype = a.Bool("type", !only), doid = a.Bool("id", !only), doqry = a.Bool("query", !only), doall = a.Bool("all", !only), docomment = a.Bool("comment", !only);
            bool dofields = a.Bool("fields", full), doaliases = a.Bool("aliases", full), doexamples = a.Bool("examples", full), donumfields = a.Bool("numfields", full), donumaliases = a.Bool("numaliases", full), donumexamples = a.Bool("numexamples", full);
            bool dodefined = a.Bool("defined", full), domethod = a.Bool("methods", full), dotaxons = a.Bool("taxscope", full), doaccess = a.Bool("access", full), dorelease = a.Bool("release", full);
            if (a.Raw("assembly") is { } && !a.Bool("assembly", true))
            {
                // the column of assembly databases is not used by the program (the line is commented out in the source)
            }

            if (!doaccess)
                doid = doqry = doall = false;
            var names = EmbossNames.ParseWithStandard(null, null);
            var sb = new StringBuilder();
            if (html)
                sb.Append("<table border cellpadding=4 bgcolor=\"#FFFFF0\">\n");
            int maxname = 14, maxmethod = 6, maxfield = 6, maxalias = 7, maxexample = 8, maxtax = 6, maxtype = 4, maxdefined = 7, maxrelease = 7;
            void Width(string name, DbDetails d)
            {
                var db = names.Database(name)!;
                maxname = Math.Max(maxname, name.Length);
                maxtype = Math.Max(maxtype, d.Type.Length);
                maxmethod = Math.Max(maxmethod, d.Methods.Length);
                maxdefined = Math.Max(maxdefined, d.Defined.Length);
                maxrelease = Math.Max(maxrelease, d.Release.Length);
                if (dofields)
                    maxfield = Math.Max(maxfield, ShowdbFields(db).Text.Length);
                if (doaliases)
                    maxalias = Math.Max(maxalias, 1);
                if (doexamples)
                    maxexample = Math.Max(maxexample, ShowdbExamples(db).Text.Length);
                if (dotaxons)
                    maxtax = Math.Max(maxtax, ShowdbTaxons(db).Text.Length);
            }

            void Head()
            {
                int imaxfield = maxfield, imaxalias = maxalias, imaxexample = maxexample;
                if (donumfields && dofields)
                    imaxfield += 7;
                if (donumaliases && doaliases)
                    imaxalias += 8;
                if (donumexamples && doexamples)
                    imaxexample += 9;
                int maxlen = new[] { maxname, maxtype, maxmethod, imaxfield, imaxalias, imaxexample, maxtax, maxdefined, maxrelease }.Max();
                string under = new('=', maxlen);
                static string Pad(int n) => n > 0 ? new string(' ', n) : "";
                sb.Append(html ? "<tr><th>Name</th>" : "# Name" + Pad(maxname - 5));
                if (dotype)
                    sb.Append(html ? "<th>Type</th>" : "Type" + Pad(maxtype - 3));
                if (doid)
                    sb.Append(html ? "<th>ID</th>" : "ID  ");
                if (doqry)
                    sb.Append(html ? "<th>Qry</th>" : "Qry ");
                if (doall)
                    sb.Append(html ? "<th>All</th>" : "All ");
                if (domethod)
                    sb.Append(html ? "<th>Method</th>" : "Method" + Pad(maxmethod - 5));
                if (dotaxons)
                    sb.Append(html ? "<th>Taxon</th>" : "Taxon" + Pad(maxtax - 4));
                if (donumfields || dofields)
                    sb.Append(html ? "<th>Fields</th>" : "Fields" + Pad(imaxfield - 5));
                if (donumaliases || doaliases)
                    sb.Append(html ? "<th>Aliases</th>" : "Aliases" + Pad(imaxalias - 6));
                if (donumexamples || doexamples)
                    sb.Append(html ? "<th>Examples</th>" : "Examples" + Pad(imaxexample - 7));
                if (dodefined)
                    sb.Append(html ? "<th>Defined</th>" : "Defined" + Pad(maxdefined - 6));
                if (dorelease)
                    sb.Append(html ? "<th>Release</th>" : "Release" + Pad(maxrelease - 6));
                if (docomment)
                    sb.Append(html ? "<th>Comment</th>" : "Comment");
                if (html)
                {
                    sb.Append("</tr>\n");
                    return;
                }

                sb.Append('\n');
                sb.Append("# ").Append(under.AsSpan(0, Math.Max(0, Math.Min(maxlen, maxname - 2)))).Append(' ');
                if (dotype)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, maxtype))).Append(' ');
                if (doid)
                    sb.Append("==  ");
                if (doqry)
                    sb.Append("=== ");
                if (doall)
                    sb.Append("=== ");
                if (domethod)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, maxmethod))).Append(' ');
                if (dotaxons)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, maxtax))).Append(' ');
                if (donumfields || dofields)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, imaxfield))).Append(' ');
                if (donumaliases || doaliases)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, imaxalias))).Append(' ');
                if (donumexamples || doexamples)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, imaxexample))).Append(' ');
                if (dodefined)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, maxdefined))).Append(' ');
                if (dorelease)
                    sb.Append(under.AsSpan(0, Math.Min(maxlen, maxrelease))).Append(' ');
                if (docomment)
                    sb.Append("=======");
                sb.Append('\n');
            }

            void Row(string name, DbDetails d)
            {
                var db = names.Database(name)!;
                if (html)
                    sb.Append($"<tr><td>{name}</td>");
                else
                    sb.Append(name.PadRight(maxname)).Append(' ');
                if (dotype)
                    sb.Append(html ? $"<td>{d.Type}</td>" : d.Type.Replace(' ', ',').PadRight(maxtype) + " ");
                foreach (var (flag, on) in new[] { (doid, d.Id), (doqry, d.Qry), (doall, d.All) })
                    if (flag)
                        sb.Append(html ? "<td>" : "").Append(on ? "OK  " : "-   ").Append(html ? "</td>" : "");
                if (domethod)
                    sb.Append(html ? $"<td>{d.Methods}</td>" : d.Methods.PadRight(maxmethod) + " ");
                if (dotaxons)
                    sb.Append(html ? "<td>" : "").Append(ShowdbTaxons(db).Text.PadRight(maxtax)).Append(' ').Append(html ? "</td>" : "");
                if (dofields || donumfields)
                {
                    var (text, count) = ShowdbFields(db);
                    sb.Append(html ? "<td>" : "");
                    if (donumfields)
                        sb.Append(count.ToString(System.Globalization.CultureInfo.InvariantCulture).PadLeft(6)).Append(' ');
                    if (dofields)
                        sb.Append(text.PadRight(maxfield)).Append(' ');
                    sb.Append(html ? "</td>" : "");
                }

                if (doaliases || donumaliases)
                {
                    sb.Append(html ? "<td>" : "");
                    if (donumaliases)
                        sb.Append("0".PadLeft(7)).Append(' ');
                    if (doaliases)
                        sb.Append("-".PadRight(maxalias)).Append(' ');
                    sb.Append(html ? "</td>" : "");
                }

                if (doexamples || donumexamples)
                {
                    var (text, count) = ShowdbExamples(db);
                    sb.Append(html ? "<td>" : "");
                    if (donumexamples)
                        sb.Append(count.ToString(System.Globalization.CultureInfo.InvariantCulture).PadLeft(8)).Append(' ');
                    if (doexamples)
                        sb.Append(text.PadRight(maxexample)).Append(' ');
                    sb.Append(html ? "</td>" : "");
                }

                if (dodefined)
                    sb.Append(html ? "<td>" : "").Append(d.Defined.PadRight(maxdefined)).Append(' ').Append(html ? "</td>" : "");
                if (dorelease)
                    sb.Append(html ? "<td>" : "").Append(d.Release.PadRight(maxrelease)).Append(' ').Append(html ? "</td>" : "");
                if (docomment)
                    sb.Append(html ? "<td>" : "").Append(d.Comment.Length > 0 ? d.Comment : "").Append(html ? "</td>" : "");
                sb.Append(html ? "</tr>\n" : "\n");
            }

            if (onedb.Length > 0)
            {
                var d = DatabaseDetails(names, onedb);
                if (d is null)
                {
                    var fatal = EmbossProgramBase.Fatal(ctx, "\n   EMBOSS An error in showdb.c at line 244:\nThe database '" + onedb + "' does not exist");
                    return outfile == "stdout" ? fatal with { Stdout = sb.ToString() } : fatal with { Files = new Dictionary<string, string> { [outfile] = sb.ToString() } };
                }

                Width(onedb, d);
                if (doheader)
                    Head();
                Row(onedb, d);
            }
            else
            {
                var list = names.DatabaseNames();
                list.Sort((x, y) =>
                {
                    int c = string.CompareOrdinal(names.Database(x)!.Defined, names.Database(y)!.Defined);
                    return c != 0 ? c : string.CompareOrdinal(x, y);
                });
                bool Selected(DbDetails d) => !(!showstandard && d.Defined == "standard") && !(!showuser && d.Defined == "user") && !(!showinclude && d.Defined is not ("user" or "standard"));
                foreach (string n in list)
                    if (DatabaseDetails(names, n) is { } d && Selected(d))
                        Width(n, d);
                if (doheader)
                    Head();
                var done = new HashSet<string>(StringComparer.Ordinal);
                void Group(bool show, Func<string, bool> match)
                {
                    foreach (string n in list)
                    {
                        var d = DatabaseDetails(names, n)!;
                        if (!Selected(d) || !show || done.Contains(n) || !match(d.Type))
                            continue;
                        Row(n, d);
                        done.Add(n);
                    }
                }

                static bool Has(string type, string word) => type.Contains(word, StringComparison.Ordinal);
                Group(showprotein, t => Has(t, "Protein") || Has(t, "Protfeatures"));
                Group(shownucleic, t => Has(t, "Nucleotide") || Has(t, "Nucfeatures") || Has(t, "Refseq"));
                Group(showsequence, t => Has(t, "Sequence"));
                Group(showfeature, t => Has(t, "Features") || (shownucleic && Has(t, "Nucfeatures")) || (showprotein && Has(t, "Protfeatures")));
                Group(showtax, t => Has(t, "Taxonomy"));
                Group(showresource, t => Has(t, "Resource"));
                Group(showobo, t => Has(t, "Obo"));
                Group(showxml, t => Has(t, "Xml"));
                Group(showtext, t => Has(t, "Text") || Has(t, "Html"));
            }

            if (html)
                sb.Append("</table>\n");
            return WriteOutfile(a, outfile, sb.ToString());
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }
    }
}
