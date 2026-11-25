// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S3776, S107, S1854, S1066, S2325, S1172
using System.IO.Compression;
using System.Reflection;
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>A term of a resource (<c>AjPResterm</c>): the identifier and the name of an EDAM or taxonomy term.</summary>
public sealed record EmbossResTerm(string Id, string Name);

/// <summary>A query line of a resource (<c>AjPResquery</c>, the <c>Query</c> line of <c>DRCAT.dat</c>): the datatype, the format, the EDAM term of the format, the identifier terms and the URL with <c>%s</c>.</summary>
public sealed record EmbossResQuery(string Datatype, string Format, string FormatTerm, string Term, string Url);

/// <summary>A cross-reference of a resource (<c>AjPReslink</c>).</summary>
public sealed record EmbossResLink(string Source, string Term);

/// <summary>
/// An entry of the data resource catalogue (<c>AjPResource</c> of <c>ajresource.c</c>, read from <c>DRCAT.dat</c> by <c>resourceinReadDrcat</c> of <c>ajresourceread.c</c>; B09 WP67, F309): the fields the reader keeps and
/// <see cref="Write"/> writes the entry as <c>resourceoutWriteDrcat</c> does. A value <c>None</c> is not stored.
/// </summary>
public sealed class EmbossResource
{
    public string Id = "", Acc = "", Name = "", Desc = "", Url = "", Urllink = "", Urlrest = "", Urlsoap = "";
    public readonly List<string> Idalt = [], Cat = [], Example = [];
    public readonly List<EmbossResTerm> Taxon = [], Edamdat = [], Edamfmt = [], Edamid = [], Edamtpc = [];
    public readonly List<EmbossResLink> Xref = [];
    public readonly List<EmbossResQuery> Query = [];

    /// <summary><c>resourceoutWriteDrcat</c>: every field the reader keeps, a blank line after the entry.</summary>
    public string Write()
    {
        var sb = new StringBuilder();
        sb.Append("ID      ").Append(Id).Append('\n');
        foreach (string a in Idalt)
            sb.Append("IDalt   ").Append(a).Append('\n');
        if (Acc.Length > 0)
            sb.Append("Acc     ").Append(Acc).Append('\n');
        sb.Append("Name    ").Append(Name).Append('\n');
        sb.Append("Desc    ").Append(Desc).Append('\n');
        if (Url.Length > 0)
            sb.Append("URL     ").Append(Url).Append('\n');
        if (Urllink.Length > 0)
            sb.Append("URLlink ").Append(Urllink).Append('\n');
        if (Urlrest.Length > 0)
            sb.Append("URLrest ").Append(Urlrest).Append('\n');
        if (Urlsoap.Length > 0)
            sb.Append("URLsoap ").Append(Urlsoap).Append('\n');
        foreach (string c in Cat)
            sb.Append("Cat     ").Append(c).Append('\n');
        Terms(sb, "Taxon   ", Taxon);
        Terms(sb, "EDAMtpc ", Edamtpc);
        Terms(sb, "EDAMdat ", Edamdat);
        Terms(sb, "EDAMid  ", Edamid);
        Terms(sb, "EDAMfmt ", Edamfmt);
        foreach (var x in Xref)
            sb.Append("Xref    ").Append(x.Source).Append(" | ").Append(x.Term).Append('\n');
        foreach (var q in Query)
            sb.Append("Query    ").Append(q.Datatype).Append(" | ").Append(q.Format).Append(" | ").Append(q.Term).Append(" | ").Append(q.Url).Append('\n');
        foreach (string e in Example)
            sb.Append("Example ").Append(e).Append('\n');
        sb.Append('\n');
        return sb.ToString();

        static void Terms(StringBuilder sb, string tag, List<EmbossResTerm> terms)
        {
            foreach (var t in terms)
                sb.Append(tag).Append(t.Id).Append(" | ").Append(t.Name).Append('\n');
        }
    }
}

/// <summary>The data resource catalogue of the EMBOSS 6.6.0 installation (<c>DRCAT.dat</c> of the data directory), read as <c>ajResourceinNewDrcat</c> reads it through the index of the <c>drcat</c> database (fields <c>id</c> and <c>acc</c>, the identifier and the alternative identifiers).</summary>
public static class EmbossDrcat
{
    private static readonly Lazy<IReadOnlyList<EmbossResource>> LazyAll = new(Load);
    private static readonly Lazy<Dictionary<string, EmbossResource>> LazyIndex = new(BuildIndex);

    /// <summary>Every entry in the order of the file.</summary>
    public static IReadOnlyList<EmbossResource> All => LazyAll.Value;

    /// <summary>The entry whose identifier or alternative identifier is <paramref name="name"/> (the index keys are lower case).</summary>
    public static EmbossResource? Lookup(string name) => LazyIndex.Value.GetValueOrDefault(name.ToLowerInvariant());

    /// <summary>The text of an embedded data file (gzip).</summary>
    internal static string ReadResource(string name) => ResourceCache.GetOrAdd(name, ReadResourceUncached);

    private static readonly System.Collections.Concurrent.ConcurrentDictionary<string, string> ResourceCache = new(StringComparer.Ordinal);

    private static string ReadResourceUncached(string name)
    {
        using var s = typeof(EmbossDrcat).Assembly.GetManifestResourceStream("Seqeron.Genomics.Emboss.SequenceFormats.Data." + name) ?? throw new InvalidOperationException($"Embedded data '{name}' not found.");
        using var gz = new GZipStream(s, CompressionMode.Decompress);
        using var r = new StreamReader(gz, Encoding.Latin1);
        return r.ReadToEnd();
    }

    private static Dictionary<string, EmbossResource> BuildIndex()
    {
        var index = new Dictionary<string, EmbossResource>(StringComparer.Ordinal);
        foreach (var e in All)
        {
            index.TryAdd(e.Id.ToLowerInvariant(), e);
            foreach (string a in e.Idalt)
                index.TryAdd(a.ToLowerInvariant(), e);
        }

        return index;
    }

    private static IReadOnlyList<EmbossResource> Load() => Parse(ReadResource("DRCAT.dat.gz"));

    // ajStrRemoveWhiteExcess: leading and trailing white space removed, runs of white space reduced to one blank
    internal static string WhiteExcess(string s)
    {
        var sb = new StringBuilder();
        bool space = false;
        foreach (char c in s)
        {
            if (c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v')
            {
                space = sb.Length > 0;
                continue;
            }

            if (space)
                sb.Append(' ');
            space = false;
            sb.Append(c);
        }

        return sb.ToString();
    }

    /// <summary><c>resourceinReadDrcat</c> over a whole file: a record starts at a line that begins with <c>ID      </c> and ends at the first line that is empty or at the next <c>ID</c> line.</summary>
    internal static IReadOnlyList<EmbossResource> Parse(string text) => Parse(text, int.MaxValue);

    /// <summary>The first <paramref name="max"/> records of the text from its start (the reader of a query reads one record from the offset of a hit).</summary>
    internal static IReadOnlyList<EmbossResource> Parse(string text, int max)
    {
        var all = new List<EmbossResource>();
        EmbossResource? cur = null;
        int linecnt = 0;
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            string line0 = nl < 0 ? text[pos..] : text[pos..nl];
            pos = nl < 0 ? text.Length : nl + 1;
            string line = line0.TrimEnd('\r');
            if (cur is null)
            {
                if (!line.StartsWith("ID      ", StringComparison.Ordinal))
                    continue;
                if (all.Count >= max)
                    break;
                cur = new EmbossResource { Id = line["ID      ".Length..].Trim() };
                all.Add(cur);
                linecnt = 0;
                continue;
            }

            if (line.Length <= 1)
            {
                cur = null;
                continue;
            }

            linecnt++;
            int sp = line.IndexOfAny([' ', '\t']);
            string name = sp < 0 ? line : line[..sp];
            string rest = sp < 0 ? "" : line[sp..].Trim();
            if (linecnt > 1 && name == "ID")
            {
                cur = new EmbossResource { Id = rest };
                all.Add(cur);
                linecnt = 0;
                continue;
            }

            if (rest == "None")
                continue;
            switch (name)
            {
                case "IDalt": cur.Idalt.Add(rest); break;
                case "Acc": cur.Acc = rest; break;
                case "Name": cur.Name = rest; break;
                case "Desc": cur.Desc = rest; break;
                case "URL": cur.Url = rest; break;
                case "URLlink": cur.Urllink = rest; break;
                case "URLrest": cur.Urlrest = rest; break;
                case "URLsoap": cur.Urlsoap = rest; break;
                case "Cat":
                    foreach (string t in rest.Split(';'))
                        cur.Cat.Add(WhiteExcess(t));
                    break;
                case "EDAMtpc": cur.Edamtpc.Add(Term(rest)); break;
                case "EDAMdat": cur.Edamdat.Add(Term(rest)); break;
                case "EDAMid": cur.Edamid.Add(Term(rest)); break;
                case "EDAMfmt": cur.Edamfmt.Add(Term(rest)); break;
                case "Taxon": cur.Taxon.Add(Term(rest)); break;
                case "Xref":
                {
                    var t = Term(rest);
                    cur.Xref.Add(new EmbossResLink(t.Id, t.Name));
                    break;
                }

                case "Query":
                {
                    string[] parts = rest.Split('|');
                    string datatype = WhiteExcess(parts[0]);
                    string format = parts.Length > 1 ? WhiteExcess(parts[1]) : "";
                    string formatTerm = "";
                    foreach (var f in cur.Edamfmt)
                        if (f.Name == format)
                            formatTerm = f.Id;
                    string term = parts.Length > 2 ? WhiteExcess(parts[2]) : "";
                    string url = parts.Length > 3 ? WhiteExcess(string.Join("|", parts[3..])) : "";
                    cur.Query.Add(new EmbossResQuery(datatype, format, formatTerm, term, url));
                    break;
                }

                case "Example": cur.Example.Add(rest); break;
            }
        }

        return all;
    }

    // the first token up to '|' and the rest, each reduced as ajStrRemoveWhiteExcess does
    private static EmbossResTerm Term(string rest)
    {
        int bar = rest.IndexOf('|');
        return bar < 0 ? new EmbossResTerm(WhiteExcess(rest), "") : new EmbossResTerm(WhiteExcess(rest[..bar]), WhiteExcess(rest[(bar + 1)..]));
    }
}
