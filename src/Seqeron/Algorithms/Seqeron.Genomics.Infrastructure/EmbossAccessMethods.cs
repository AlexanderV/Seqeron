namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The access methods of EMBOSS 6.6.0 by data type (<c>ajaxdb/aj*db.c</c>: the lists <c>textAccess</c>, <c>seqAccess</c>, <c>feattabAccess</c>, <c>taxAccess</c>, <c>urlAccess</c>; the lists of the other types are empty) and what each can do: one entry (<c>id</c>), a query (<c>qry</c>), all
/// entries (<c>all</c>): <c>ajTextaccessMethodGetScope</c> and the functions of the types, as <c>namMethod2Scope</c> uses them (B09 WP67, F309: <c>dbtell</c>, <c>showdb</c>).
/// </summary>
internal static class EmbossAccessMethods
{
    private const int Entry = 1, Query = 2, All = 4;

    private static readonly Dictionary<string, int> Text = new(StringComparer.Ordinal)
    {
        ["dbfetch"] = Entry,
        ["wsdbfetch"] = Entry | Query,
        ["ebeye"] = Entry | Query,
        ["emboss"] = Entry | Query | All,
        ["emblcd"] = Entry | Query | All,
        ["entrez"] = Entry | Query | All,
        ["seqhound"] = Entry | Query | All,
        ["mrs"] = Entry | Query,
        ["mrs4"] = Entry | Query,
        ["mrs3"] = Entry | Query,
        ["obda"] = Entry | Query | All,
        ["srs"] = Entry | Query | All,
        ["srsfasta"] = Entry | Query | All,
        ["srswww"] = Entry | Query | All,
        ["url"] = Entry | Query | All,
        ["app"] = Entry | Query | All,
        ["external"] = Entry | Query | All,
        ["direct"] = Entry | Query | All,
    };

    private static readonly Dictionary<string, Dictionary<string, int>> ByType = new(StringComparer.Ordinal)
    {
        ["sequence"] = new(StringComparer.Ordinal) { ["gcg"] = Entry | Query | All, ["embossgcg"] = Entry | Query | All, ["blast"] = Entry | Query | All, ["biomart"] = Entry | Query, ["sql"] = Entry | Query | All, ["ensembl"] = Entry | Query, ["ensemblgenomes"] = Entry | Query, ["das"] = Entry | Query },
        ["features"] = new(StringComparer.Ordinal) { ["das"] = Entry | Query, ["chado"] = Entry, ["gcg"] = Entry | Query | All, ["embossgcg"] = Entry | Query | All },
        ["taxon"] = new(StringComparer.Ordinal) { ["embosstax"] = Entry | Query | All },
        ["url"] = new(StringComparer.Ordinal) { ["urlonly"] = Entry | Query | All },
    };

    /// <summary>The data type of a database type name (<c>namDbTypes</c>, any case); null when unknown.</summary>
    internal static string? DataType(string type) => type.ToLowerInvariant() switch
    {
        "n" or "p" or "nucleotide" or "protein" or "sequence" => "sequence",
        "nucfeatures" or "protfeatures" or "features" => "features",
        "assembly" => "assembly",
        "obo" => "obo",
        "refseq" => "refseq",
        "resource" => "resource",
        "taxonomy" => "taxon",
        "html" or "text" or "unknown" => "text",
        "url" => "url",
        "variation" => "variation",
        "xml" => "xml",
        _ => null,
    };

    /// <summary><c>namMethod2Scope</c> as a number of bits: 1 = entry (<c>id</c>), 2 = query, 4 = all; the text methods first, then those of the data types of <paramref name="dbtype"/> (separated by space, comma or semicolon; an unknown type ends the search).</summary>
    internal static int Scope(string method, string dbtype)
    {
        if (Text.TryGetValue(method, out int scope))
            return scope;
        int result = 0;
        foreach (string token in dbtype.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries))
        {
            string? type = DataType(token);
            if (type is null)
                return result;
            if (ByType.TryGetValue(type, out var table) && table.TryGetValue(method, out int s))
                result |= s;
        }

        return result;
    }

    internal static bool HasEntry(int scope) => (scope & Entry) != 0;

    internal static bool HasQuery(int scope) => (scope & Query) != 0;

    internal static bool HasAll(int scope) => (scope & All) != 0;
}
