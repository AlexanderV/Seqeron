// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The feature vocabularies of EMBOSS 6.6.0 (<c>ajfeat.c</c>: <c>featVocabRead</c>, <c>featVocabReadTags</c>, <c>featVocabReadTypes</c>, <c>featFeatType</c>, <c>featTagName</c>, <c>featTableTypeInternal</c>, <c>ajFeattypeGetExternal</c>,
/// <c>ajFeattagGetNameS</c>, <c>ajFeattagFormat</c>; the <c>featVocabDef</c> table; EMBOSS 6.6.0, kimrutherford/EMBOSS at raw.githubusercontent.com). The tables are built from the embedded data files
/// (<c>EmbossFeatureVocabData.g.cs</c>) by the ported loader, including the aliasing of the <c>name = id</c> lines (the value string of the id is edited IN PLACE and the new name shares it). A5-8c3e, B09 F241.
/// </summary>
public static partial class EmbossFeatureVocab
{
    /// <summary>A mutable string: the <c>AjPStr</c> a table value points to (two keys can share one).</summary>
    public sealed class Box(string v)
    {
        public string V = v;
    }

    /// <summary>An <c>ajTablestrNewCase</c> table: ASCII case-insensitive keys; a <c>put</c> on an existing key replaces the KEY (<c>tableDelStr</c>) and the value; the insertion order is kept for the dumps.</summary>
    public sealed class Table
    {
        private readonly Dictionary<string, (string Key, Box Val)> _map = new(AsciiCaseComparer.Instance);

        public int Count => _map.Count;

        public IEnumerable<(string Key, string Value)> Items => _map.Values.Select(e => (e.Key, e.Val.V));

        /// <summary><c>ajTableFetchS</c>: the value, null when the key is absent.</summary>
        public string? Fetch(string key) => _map.TryGetValue(key, out var e) ? e.Val.V : null;

        public Box? FetchBox(string key) => _map.TryGetValue(key, out var e) ? e.Val : null;

        /// <summary><c>ajTablestrFetchkeyS</c>: the key as it is spelled in the table, null when absent.</summary>
        public string? FetchKey(string key) => _map.TryGetValue(key, out var e) ? e.Key : null;

        /// <summary><c>ajTablePut</c>; true when the key was there (the earlier value is returned by the C function, which the loader only tests for non-null).</summary>
        public bool Put(string key, Box value)
        {
            bool had = _map.ContainsKey(key);
            _map.Remove(key);
            _map[key] = (key, value);
            return had;
        }

        public bool Put(string key, string value) => Put(key, new Box(value));
    }

    private sealed class AsciiCaseComparer : IEqualityComparer<string>
    {
        public static readonly AsciiCaseComparer Instance = new();

        public bool Equals(string? a, string? b)
        {
            if (a is null || b is null)
                return a is null && b is null;
            if (a.Length != b.Length)
                return false;
            for (int i = 0; i < a.Length; i++)
                if (Lower(a[i]) != Lower(b[i]))
                    return false;
            return true;
        }

        public int GetHashCode(string s)
        {
            int h = 17;
            foreach (char c in s)
                h = unchecked(h * 31 + Lower(c));
            return h;
        }

        private static char Lower(char c) => c is >= 'A' and <= 'Z' ? (char)(c + 32) : c;
    }

    /// <summary>One vocabulary: the type table and the tag table (<c>FeatTypeTable*</c> / <c>FeatTagsTable*</c>).</summary>
    public sealed class Vocabulary(Table types, Table tags)
    {
        public readonly Table Types = types;
        public readonly Table Tags = tags;
    }

    private static readonly Dictionary<string, Vocabulary> Loaded = new(StringComparer.Ordinal);
    private static readonly object Gate = new();
    private static Dictionary<string, string>? _files;

    // featVocabDef: the name of a feature format -> the name of the vocabulary of a nucleotide / a protein feature
    private static readonly (string Name, string Truename, bool Nuc, bool Prot)[] VocabDef =
    [
        ("gff3", "gff3", true, false), ("gff3", "gff3protein", false, true), ("gff", "gff3", true, false), ("gff", "gff3protein", false, true), ("gff2", "gff2", true, false), ("gff2", "gff2protein", false, true),
        ("embl", "embl", true, false), ("em", "embl", true, false), ("genbank", "embl", true, false), ("gb", "embl", true, false), ("ddbj", "embl", true, false), ("refseq", "embl", true, false),
        ("refseqp", "refseqp", false, true), ("pir", "pir", false, true), ("nbrf", "pir", false, true), ("sw", "swiss", false, true), ("swissprot", "swiss", false, true),
        ("dasgff", "gff3", true, false), ("dasgff", "gff3protein", false, true), ("debug", "emboss", true, false), ("debug", "protein", false, true),
    ];

    /// <summary>The vocabulary of a true name (<c>embl</c>, <c>emboss</c>, <c>gff2</c>, <c>gff2protein</c>, <c>gff3</c>, <c>gff3protein</c>, <c>pir</c>, <c>protein</c>, <c>refseqp</c>, <c>swiss</c>), read on first use (<c>featVocabInit</c>).</summary>
    internal static Vocabulary Get(string truename)
    {
        lock (Gate)
        {
            if (Loaded.TryGetValue(truename, out var v))
                return v;
            _files ??= DataFiles.ToDictionary(f => f.File, f => f.Text, StringComparer.Ordinal);
            var types = new Table();
            var tags = new Table();
            ReadTags("Etags." + truename, tags, true);
            ReadTypes("Efeatures." + truename, types, tags, true);
            v = new Vocabulary(types, tags);
            Loaded[truename] = v;
            return v;
        }
    }

    /// <summary><c>ajFeatVocabInitNuc</c> / <c>ajFeatVocabInitProt</c>: the vocabulary of a feature FORMAT name for nucleotide or protein features.</summary>
    internal static Vocabulary? ForFormat(string name, bool protein)
    {
        foreach (var d in VocabDef)
        {
            if (protein ? !d.Prot : !d.Nuc)
                continue;
            if (name == d.Truename || name == d.Name)
                return Get(d.Truename);
        }

        return null;
    }

    /// <summary>The internal vocabularies: <c>FeatTypeTableDna</c> (<c>emboss</c>) and <c>FeatTypeTableProtein</c> (<c>protein</c>).</summary>
    public static Vocabulary Dna => Get("emboss");

    public static Vocabulary Protein => Get("protein");

    private static IEnumerable<string> Lines(string name)
    {
        if (!_files!.TryGetValue(name, out var text))
            throw new InvalidOperationException($"Unable to read feature data file '{name}'");
        foreach (string l in text.Split('\n'))
            yield return l.TrimEnd('\r');
    }

    // ajStrRemoveWhiteExcess: leading and trailing white removed, inner runs of white become one blank
    internal static string RemoveWhiteExcess(string s)
    {
        var sb = new StringBuilder();
        bool pendingSpace = false;
        foreach (char c in s)
        {
            if (IsSpace(c))
            {
                pendingSpace = sb.Length > 0;
                continue;
            }

            if (pendingSpace)
                sb.Append(' ');
            pendingSpace = false;
            sb.Append(c);
        }

        return sb.ToString();
    }

    private static bool IsSpace(char c) => c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v';

    private static string FirstWord(string s)
    {
        int i = 0;
        while (i < s.Length && IsSpace(s[i]))
            i++;
        int j = i;
        while (j < s.Length && !IsSpace(s[j]))
            j++;
        return s[i..j];
    }

    private static readonly string[] TagTypes = ["QTEXT", "QWORD", "BTEXT", "TEXT", "VOID", "LIMITED", "QLIMITED", "SPECIAL", "QSPECIAL", "TAGVAL"];

    // featTagName
    private static bool TagName(string line, out string name, out string type, out string rest)
    {
        name = type = rest = "";
        int p = 0;
        while (p < line.Length && IsSpace(line[p]))
            p++;
        int q = p;
        while (p < line.Length && !IsSpace(line[p]))
            p++;
        name = line[q..p];
        while (p < line.Length && IsSpace(line[p]))
            p++;
        if (p >= line.Length)
        {
            name = "";
            return false;
        }

        q = p;
        while (p < line.Length && char.IsAsciiLetter(line[p]))
            p++;
        type = line[q..p];
        if (type != "LIMITED" || type == "QLIMITED")
            return true;
        while (p < line.Length && IsSpace(line[p]))
            p++;
        rest = line[p..];
        return true;
    }

    // featVocabReadTags
    private static void ReadTags(string fname, Table tags, bool recursion)
    {
        int tagscount = tags.Count;
        foreach (string raw in Lines(fname))
        {
            string line = RemoveWhiteExcess(raw);
            if (line.StartsWith('#'))
            {
                if (recursion && line.StartsWith("#include ", StringComparison.OrdinalIgnoreCase))
                {
                    ReadTags(FirstWord(line[9..]), tags, true);
                    tagscount = tags.Count;
                }

                if (recursion && line.StartsWith("#includeonly ", StringComparison.OrdinalIgnoreCase))
                {
                    ReadTags(FirstWord(line[13..]), tags, false);
                    tagscount = tags.Count;
                }
            }
            else if (line.Length > 0 && TagName(line, out string tagname, out string tagtype, out string rest))
            {
                int numtype = Array.IndexOf(TagTypes, tagtype);
                if (numtype < 0)
                    break;
                tagscount++;
                if (tagscount == 1)
                    tags.Put("", tagname);
                var tagstr = new StringBuilder(TagTypes[numtype]).Append(';');
                if (TagTypes[numtype] is "LIMITED" or "QLIMITED")
                {
                    foreach (string token in rest.Split(['"', ',', ' ', '\t'], StringSplitOptions.RemoveEmptyEntries))
                        tagstr.Append(token).Append(';');
                }

                tags.Put(tagname, tagstr.ToString());
            }
        }
    }

    // featFeatType
    private static bool FeatType(string line, out string type, out string intids, out string tag, out string req)
    {
        type = intids = tag = req = "";
        var reqb = new StringBuilder();
        var idsb = new StringBuilder();
        try
        {
            return FeatTypeCore(line, out type, out tag, reqb, idsb);
        }
        finally
        {
            req = reqb.ToString();
            intids = idsb.ToString();
        }
    }

    private static bool FeatTypeCore(string line, out string type, out string tag, StringBuilder req, StringBuilder intids)
    {
        type = tag = "";
        int p = 0;
        while (p < line.Length && IsSpace(line[p]))
            p++;
        if (p >= line.Length)
            return false;
        bool istag = false;
        char C(int i) => i < line.Length ? line[i] : '\0';
        if (C(p) == '/')
        {
            istag = true;
        }
        else if (C(p) is 'm' or 'M' or '-')
        {
            if (C(p + 1) == '/')
            {
                istag = true;
                req.Append(line[p++]);
            }
        }
        else if (C(p) == '+')
        {
            if (C(p + 1) is 'm' or 'M')
            {
                if (C(p + 2) == '/')
                {
                    istag = true;
                    req.Append(line[p++]);
                    req.Append(line[p++]);
                }
            }
            else if (C(p + 1) == '/')
            {
                istag = true;
                req.Append(line[p++]);
            }
        }

        if (istag)
        {
            p++;
            int q = p;
            int n = 0;
            while (p < line.Length && !IsSpace(line[p++]))
                n++;
            tag = line[q..Math.Min(q + n, line.Length)];
            return tag.Length > 0;
        }

        {
            int q = p;
            int n = 0;
            // while(*cp && !isspace(*cp++)): the pointer moves past the first white too
            while (p < line.Length && !IsSpace(line[p++]))
                n++;
            type = line[q..(q + n)];
            if (type.Length == 0)
                return false;
            while (p < line.Length && IsSpace(line[p]))
                p++;
            if (p >= line.Length)
                return true;
            if (line[p] == '=')
                req.Append(line[p++]);
            while (p < line.Length && IsSpace(line[p]))
                p++;
            q = p;
            n = 0;
            while (p < line.Length && !IsSpace(line[p++]))
                n++;
            intids.Append(line, q, n);
            while (p < line.Length)
            {
                while (p < line.Length && IsSpace(line[p]))
                    p++;
                if (p >= line.Length)
                    return true;
                q = p;
                n = 0;
                while (p < line.Length && !IsSpace(line[p++]))
                    n++;
                intids.Append(' ').Append(line, q, n);
            }
        }

        return true;
    }

    // featVocabReadTypes
    private static void ReadTypes(string fname, Table types, Table tags, bool recursion)
    {
        string? savetype = null;
        Box? typtagstr = null;
        bool ismodtype = false, newdefid = false;
        Box? defid = null;
        int filetypecount = 0;
        int typecount = types.Count;
        foreach (string raw in Lines(fname))
        {
            string line = RemoveWhiteExcess(raw);
            if (line.StartsWith('#'))
            {
                if (recursion && line.StartsWith("#include ", StringComparison.OrdinalIgnoreCase))
                {
                    ReadTypes(FirstWord(line[9..]), types, tags, true);
                    typecount = types.Count;
                }

                if (recursion && line.StartsWith("#includeonly ", StringComparison.OrdinalIgnoreCase))
                {
                    ReadTypes(FirstWord(line[13..]), types, tags, false);
                    typecount = types.Count;
                }

                continue;
            }

            if (line.Length == 0 || !FeatType(line, out string type, out string intids, out string tag, out string req))
                continue;
            if (type.Length > 0)
            {
                if (!ismodtype && !string.IsNullOrEmpty(savetype))
                {
                    if (types.Fetch(savetype) is null)
                        types.Put(savetype, typtagstr!);
                    typtagstr = null;
                    savetype = null;
                }

                if (req == "=")
                {
                    ismodtype = true;
                    if (type == ".")
                        type = "";
                    if (types.Fetch(type) is null)
                        types.Put(type, intids);
                    var slot = types.FetchBox(intids) ?? throw new InvalidOperationException($"{fname}: undefined internal ID '{intids}'");
                    int ipos = slot.V.IndexOf(';', StringComparison.Ordinal);
                    if (ipos >= 0)
                        slot.V = slot.V[ipos..];
                    else
                        slot.V = ";" + slot.V;
                    slot.V = type + slot.V;
                    savetype = type;
                    typtagstr = slot;
                }
                else
                {
                    typecount++;
                    filetypecount++;
                    if (typecount == 1)
                    {
                        defid = new Box(type);
                        newdefid = true;
                        if (types.Fetch("") is null)
                            types.Put("", defid);
                    }

                    ismodtype = false;
                    typtagstr = new Box(";");
                    if (intids.Length > 0)
                    {
                        string[] ids = intids.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                        string firstid = ids[0];
                        string storetype = type;
                        bool aliased = false;
                        foreach (string sofaid in ids)
                        {
                            if (newdefid)
                            {
                                defid!.V = sofaid;
                                newdefid = false;
                                defid = null;
                            }

                            // storetype = sofaid is a POINTER to the static token string of ajStrParseWhite: it holds the current token from then on
                            if (aliased)
                                storetype = sofaid;
                            if (!AsciiCaseComparer.Instance.Equals(storetype, firstid) && types.Fetch(storetype) is null)
                            {
                                types.Put(storetype, firstid);
                                aliased = true;
                            }
                        }

                        typtagstr.V = string.Concat(type, typtagstr.V);
                        savetype = firstid;
                    }
                    else
                    {
                        savetype = type;
                    }
                }
            }
            else
            {
                char modchar = req.Length > 0 ? req[0] : '\0';
                char reqchar = req.Length > 0 ? req[^1] : '\0';
                if (modchar == '-')
                {
                    typtagstr!.V = typtagstr.V.Replace(";" + tag + ";", ";", StringComparison.Ordinal);
                }
                else
                {
                    typtagstr!.V = string.Concat(typtagstr.V, reqchar is 'M' or 'm' ? "*" : "", tag, ";");
                }
            }
        }

        if (filetypecount > 0 && savetype is not null && types.Fetch(savetype) is null)
            types.Put(savetype, typtagstr!);
    }

    // ------------------------------------------------------------------------------------------------------------------ lookups

    /// <summary><c>ajFeattagGetNameS</c>: the canonical spelling of a tag, or the default (the first tag read) with <paramref name="known"/> false.</summary>
    internal static string TagNameOf(Table tags, string? tag, out bool known)
    {
        string def = tags.Fetch("") ?? "";
        if (tag is not null)
        {
            string? key = tags.FetchKey(tag);
            if (key is not null)
            {
                known = true;
                return key;
            }
        }

        known = false;
        return def;
    }

    /// <summary><c>ajFeattagFormat</c>: the value type of a canonical tag (the text of its value up to the first semicolon).</summary>
    internal static string TagFormat(Table tags, string name)
    {
        string v = tags.Fetch(name) ?? "";
        int i = v.IndexOf(';', StringComparison.Ordinal);
        return i < 0 ? v : v[..i];
    }

    /// <summary><c>ajFeattagGetLimit</c>: the semicolon-introduced list of valid values (empty when none).</summary>
    internal static string TagLimit(Table tags, string name)
    {
        string v = tags.Fetch(name) ?? "";
        int i = v.IndexOf(';', StringComparison.Ordinal);
        return i < 0 ? "" : v[i..];
    }

    /// <summary><c>featTableTypeInternal</c>: the internal type for a type name of an input format; the default (the key <c>""</c>) when unknown.</summary>
    internal static string? TypeInternal(Table table, string type)
    {
        string? retkey = table.FetchKey(type);
        if (retkey is not null)
        {
            string tmpstr = table.Fetch(retkey) ?? "";
            if (tmpstr.Length == 0 || tmpstr[^1] != ';')
                retkey = table.FetchKey(tmpstr);
            return retkey;
        }

        return table.Fetch("");
    }

    /// <summary><c>ajFeattypeGetExternal</c>: the type name of the vocabulary of an output format for an internal type.</summary>
    public static string TypeExternal(Table table, string type)
    {
        string? retkey = table.FetchKey(type);
        if (retkey is null)
            retkey = table.Fetch("") ?? "";
        string ret = table.Fetch(retkey) ?? "";
        if (ret.Length == 0 || ret[^1] != ';')
        {
            retkey = ret;
            ret = table.Fetch(retkey) ?? "";
        }

        if (ret.Length > 0 && ret[^1] == ';')
        {
            if (ret[0] == ';')
            {
                ret = retkey;
            }
            else
            {
                int i = ret.IndexOf(';', StringComparison.Ordinal);
                string tmp = ret[..i];
                string? k = table.FetchKey(tmp);
                ret = k ?? (table.Fetch(retkey) ?? "");
            }
        }

        return ret;
    }
}
