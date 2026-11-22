#pragma warning disable S107, S3776, S1121, S3358, S1066, S1854, S2325, S3923

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>One file of an EMBL CD-ROM index (<c>TextPCdFile</c>): the bytes, the header (300 bytes: the size of the file, the number of records, the size of a record, ...) and the records after it.</summary>
internal sealed class EmbossCdFile
{
    private readonly byte[] _data;

    public EmbossCdFile(byte[] data)
    {
        _data = data;
        // textCdReadHeader: FileSize (uint), NRecords (uint), RecSize (int16)
        NRecords = _data.Length >= 8 ? BitConverter.ToUInt32(_data, 4) : 0;
        RecSize = _data.Length >= 10 ? BitConverter.ToInt16(_data, 8) : (short)0;
    }

    public uint NRecords { get; }

    public short RecSize { get; }

    // ajReadbinCharTrim: the bytes of the field up to a zero byte, without the blanks at their end
    private string Field(int pos, int size)
    {
        int e = pos;
        int limit = Math.Min(_data.Length, pos + Math.Max(0, size));
        while (e < limit && _data[e] != 0)
            e++;
        string t = pos < _data.Length ? Encoding.Latin1.GetString(_data, pos, Math.Max(0, e - pos)) : "";
        return t.TrimEnd(' ');
    }

    private int Pos(uint ipos) => 300 + (int)(ipos * (uint)RecSize);

    /// <summary><c>textCdIdxName</c>: the name of a record of an index of entries (the record is <c>RecSize - 10</c> characters, then the offset of the entry, the offset of its sequence and the division).</summary>
    public string IdxName(uint ipos) => Field(Pos(ipos), RecSize - 10);

    /// <summary><c>textCdIdxLine</c>.</summary>
    public (string Name, uint AnnOffset, uint SeqOffset, short DivCode) IdxLine(uint ipos)
    {
        int p = Pos(ipos), n = RecSize - 10;
        string name = Field(p, n);
        uint ann = p + n + 4 <= _data.Length ? BitConverter.ToUInt32(_data, p + n) : 0;
        uint sec = p + n + 8 <= _data.Length ? BitConverter.ToUInt32(_data, p + n + 4) : 0;
        short div = p + n + 10 <= _data.Length ? BitConverter.ToInt16(_data, p + n + 8) : (short)0;
        return (name, ann, sec, div);
    }

    /// <summary><c>textCdTrgName</c>: the target of a record of a target file (the record is the number of hits and the first hit, then <c>RecSize - 8</c> characters).</summary>
    public string TrgName(uint ipos) => Field(Pos(ipos) + 8, RecSize - 8);

    /// <summary><c>textCdTrgLine</c>.</summary>
    public (uint NHits, uint FirstHit, string Target) TrgLine(uint ipos)
    {
        int p = Pos(ipos);
        uint nhits = p + 4 <= _data.Length ? BitConverter.ToUInt32(_data, p) : 0;
        uint first = p + 8 <= _data.Length ? BitConverter.ToUInt32(_data, p + 4) : 0;
        return (nhits, first, Field(p + 8, RecSize - 8));
    }

    /// <summary>A record of the file of hits (the number of an entry of the index of entries, 1 based).</summary>
    public int HitAt(uint record)
    {
        int p = Pos(record);
        return p + 4 <= _data.Length ? BitConverter.ToInt32(_data, p) : 0;
    }

    /// <summary>The division file: the name of a division (<c>textCdQryOpen</c> / <c>textCdQryFile</c>): the code (int16) and the name (<c>RecSize - 2</c> characters).</summary>
    public string DivisionName(uint ipos, bool trim)
    {
        int p = Pos(ipos) + 2;
        return trim ? Field(p, RecSize - 2) : RawField(p, RecSize - 2);
    }

    private string RawField(int pos, int size)
    {
        // ajReadbinChar into a C string: the text ends at the first zero byte
        int limit = Math.Min(_data.Length, pos + Math.Max(0, size));
        int e = pos;
        while (e < limit && _data[e] != 0)
            e++;
        return pos < _data.Length ? Encoding.Latin1.GetString(_data, pos, Math.Max(0, e - pos)) : "";
    }
}

/// <summary>
/// The database access method <c>emblcd</c> (<c>textAccessEmblcd</c> of <c>ajtextdb.c</c>): the EMBL CD-ROM index files of the <c>dbi</c> programs (<c>division.lkp</c>, <c>entrynam.idx</c>, and <c>FIELD.trg</c> / <c>FIELD.hit</c> for acnum, seqvn, taxon, keyword, des, gi)
/// are searched for the query and the entries are read from the database files at the offsets of the hits, one entry per hit. A5-8c3h, B09 F276.
/// </summary>
internal static class EmbossCdAccess
{
    private const char Slash = '/';
    private static readonly string[] EmbossFields = ["id", "acc", "sv", "org", "key", "des", "gi"];
    private static readonly string?[] CdFields = [null, "acnum", "seqvn", "taxon", "keyword", "des", "gi"];

    private sealed class Index
    {
        public required EmbossCdFile Division, Entries;
        public required bool[] Skip;
        public required EmbossQuery Query;
        public required EmbossQueryEnv Env;
        public List<EmbossBtHit> Results = [];

        public EmbossCdFile? Open(string name)
        {
            string path = Query.IndexDir.Length == 0 ? name : $"{Query.IndexDir}{Slash}{name}";
            string? text = Env.Files(path);
            return text is null ? null : new EmbossCdFile(Encoding.Latin1.GetBytes(text));
        }

        public bool TrgOpen(string name, out EmbossCdFile? trg, out EmbossCdFile? hit)
        {
            hit = null;
            trg = Open(name + ".trg");
            if (trg is null)
                return false;
            hit = Open(name + ".hit");
            return hit is not null;
        }

        public void AddHit(uint ipos)
        {
            var line = Entries.IdxLine(ipos);
            if (!Skip[line.DivCode - 1])
                Results.Add(new EmbossBtHit((uint)line.DivCode, line.AnnOffset));
        }
    }

    private static int CmpC(string a, string b) => string.CompareOrdinal(a, b);

    // textCdIdxSearch: the position of an entry name (binary search), -1 when there is none
    private static int IdxSearch(Index ix, string entry, out (string Name, uint AnnOffset, uint SeqOffset, short DivCode) line)
    {
        line = default;
        string key = entry.ToUpperInvariant();
        var fil = ix.Entries;
        if (fil.NRecords < 1)
            return -1;
        int ilo = 0, ihi = (int)fil.NRecords - 1, ipos = 0, icmp = 0;
        while (ilo <= ihi)
        {
            ipos = (ilo + ihi) / 2;
            icmp = CmpC(key, fil.IdxName((uint)ipos));
            if (icmp == 0)
                break;
            if (icmp < 0)
                ihi = ipos - 1;
            else
                ilo = ipos + 1;
        }

        if (icmp != 0)
            return -1;
        line = fil.IdxLine((uint)ipos);
        return ipos;
    }

    // textCdTrgSearch: the position of a target (the strange halving of the library), -1 when there is none or it has no hits
    private static int TrgSearch(string entry, EmbossCdFile fp, out (uint NHits, uint FirstHit, string Target) line)
    {
        line = default;
        string key = entry.ToUpperInvariant();
        if (fp.NRecords < 1)
            return -1;
        int ilo = 0, ihi = (int)fp.NRecords;
        int ipos = (ilo + ihi) / 2;
        int icmp = CmpC(key, fp.TrgName((uint)ipos));
        while (icmp != 0)
        {
            if (icmp < 0)
                ihi = ipos;
            else
                ilo = ipos;
            int itry = (ilo + ihi) / 2;
            if (itry == ipos)
                return -1;
            ipos = itry;
            icmp = CmpC(key, fp.TrgName((uint)ipos));
        }

        line = fp.TrgLine((uint)ipos);
        return line.NHits == 0 ? -1 : ipos;
    }

    private static string RemoveWild(string s)
    {
        int i = s.IndexOfAny(['?', '*']);
        return i >= 0 ? s[..i] : s;
    }

    // textCdIdxQuery: the entries whose name matches the wildcard
    private static void IdxQuery(Index ix, string idqry)
    {
        string idstr = idqry.ToUpperInvariant();
        string idpref = RemoveWild(idstr);
        var fil = ix.Entries;
        int jlo, ilo = jlo = 0;
        int jhi, ihi = jhi = (int)fil.NRecords - 1;
        int khi = jhi;
        int ilen = idpref.Length;
        int idlen = fil.RecSize - 10;
        if (idlen < ilen)
            ilen = idlen;
        bool first = true;
        if (ilen > 0)
        {
            while (ilo <= ihi)
            {
                int ipos = (ilo + ihi) / 2;
                string name = Trunc(fil.IdxName((uint)ipos), ilen);
                int icmp = CmpC(idpref, name);
                if (icmp == 0)
                {
                    if (first)
                    {
                        jhi = ihi;
                        first = false;
                        khi = ipos;
                    }

                    jlo = ipos;
                }

                if (icmp > 0)
                    ilo = ipos + 1;
                else
                    ihi = ipos - 1;
            }

            if (first)
                return;
            ilo = jlo + 1;
            ihi = jhi;
            while (ilo <= ihi)
            {
                int ipos = (ilo + ihi) / 2;
                string name = Trunc(fil.IdxName((uint)ipos), ilen);
                int icmp = CmpC(idpref, name);
                if (icmp == 0)
                    khi = ipos;
                if (icmp < 0)
                    ihi = ipos - 1;
                else
                    ilo = ipos + 1;
            }
        }

        for (int i = jlo; i < khi + 1; i++)
        {
            var line = fil.IdxLine((uint)i);
            if (EmbossWild.Match(line.Name, idstr) && !ix.Skip[line.DivCode - 1])
                ix.Results.Add(new EmbossBtHit((uint)line.DivCode, line.AnnOffset));
        }
    }

    private static string Trunc(string s, int n) => s.Length > n ? s[..n] : s;

    // textCdTrgFind: the entries of the targets that match the wildcard
    private static int TrgFind(Index ix, string indexname, string queryName)
    {
        if (!ix.TrgOpen(indexname, out var trgfp, out var hitfp))
            return 0;
        string fdstr = queryName.ToUpperInvariant();
        string fdprefix = RemoveWild(fdstr);
        int b = 0, b2 = 0;
        int t = (int)trgfp!.NRecords - 1;
        int t2 = t, t3 = t;
        int prefixlen = fdprefix.Length;
        int trglen = trgfp.RecSize - 8;
        if (trglen < prefixlen)
            prefixlen = trglen;
        bool first = true;
        if (prefixlen > 0)
        {
            while (b <= t)
            {
                int pos = (t + b) / 2;
                string name = Trunc(trgfp.TrgName((uint)Math.Max(pos, 0)), prefixlen);
                int cmp = CmpC(fdprefix, name);
                if (cmp == 0)
                {
                    if (first)
                    {
                        first = false;
                        t2 = t;
                        t3 = pos;
                    }

                    b2 = pos;
                }

                if (cmp > 0)
                    b = pos + 1;
                else
                    t = pos - 1;
            }

            if (first)
                return 0;
            b = b2 - 1;
            t = t2;
            while (b <= t)
            {
                int pos = (t + b) / 2;
                string name = Trunc(trgfp.TrgName((uint)Math.Max(pos, 0)), prefixlen);
                int cmp = CmpC(fdprefix, name);
                if (cmp == 0)
                    t3 = pos;
                if (cmp < 0)
                    t = pos - 1;
                else
                    b = pos + 1;
            }
        }

        int start = b2, end = t3;
        for (int i = start; i < end + 1; ++i)
        {
            string name = trgfp.TrgName((uint)i);
            if (!EmbossWild.Match(name, fdstr))
                continue;
            var line = trgfp.TrgLine((uint)i);
            for (uint j = 0; j < line.NHits; ++j)
            {
                int k = hitfp!.HitAt(line.FirstHit - 1 + j) - 1;
                ix.AddHit((uint)k);
            }
        }

        return ix.Results.Count;
    }

    // textCdTrgQuery
    private static void TrgQuery(Index ix, string field, string wildqry)
    {
        int ret = 0;
        if (field == "org")
            ret += TrgFind(ix, "taxon", wildqry);
        if (field == "key")
            ret += TrgFind(ix, "keyword", wildqry);
        if (field == "des")
            ret += TrgFind(ix, "des", wildqry);
        if (field == "sv")
            ret += TrgFind(ix, "seqvn", wildqry);
        if (field == "gi")
            ret += TrgFind(ix, "gi", wildqry);
        if (ix.Query.HasAcc && field == "acc")
            ret += TrgFind(ix, "acnum", wildqry);
        _ = ret;
    }

    // textCdQryEntry
    private static bool QryEntry(Index ix)
    {
        int ipos = -1;
        foreach (var fd in ix.Query.Fields)
        {
            if (fd.Link == 'E' && ix.Results.Count > 0)
                continue;
            int i = Array.IndexOf(EmbossFields, fd.Field);
            if (i >= 0)
            {
                if (CdFields[i] is null)
                {
                    ipos = IdxSearch(ix, fd.Wildquery, out var line);
                    if (ipos >= 0 && !ix.Skip[line.DivCode - 1])
                        ix.Results.Add(new EmbossBtHit((uint)line.DivCode, line.AnnOffset));
                }
                else if (ix.TrgOpen(CdFields[i]!, out var trg, out var hit))
                {
                    int trghit = TrgSearch(fd.Wildquery, trg!, out var tl);
                    if (trghit >= 0)
                    {
                        ipos = (int)tl.FirstHit;
                        for (uint ii = 0; ii < tl.NHits; ii++)
                            ix.AddHit((uint)(hit!.HitAt(tl.FirstHit - 1 + ii) - 1));
                    }
                }
            }
            else if (ix.TrgOpen(fd.Field, out var trg, out var hit))
            {
                int trghit = TrgSearch(fd.Wildquery, trg!, out var tl);
                if (trghit >= 0)
                {
                    ipos = (int)tl.FirstHit;
                    for (uint ii = 0; ii < tl.NHits; ii++)
                        ix.AddHit((uint)(hit!.HitAt(tl.FirstHit - 1 + ii) - 1));
                }
            }
        }

        return ipos >= 0 && ix.Results.Count > 0;
    }

    // textCdQryQuery: the hits of every field combined with the link operators (the key of a hit is (file << 4) + offset), then sorted and made unique
    private static bool QryQuery(Index ix)
    {
        var table = new EmbossHitTable(100, true);
        ulong lasthits = 0;
        foreach (var field in ix.Query.Fields)
        {
            if (field.Link == 'E' && lasthits > 0)
                continue;
            ix.Results = [];
            if (field.Field == "id")
                IdxQuery(ix, field.Wildquery);
            else
                TrgQuery(ix, field.Field, field.Wildquery);
            ulong fdhits = (ulong)ix.Results.Count;
            switch (field.Link)
            {
                case 'I':
                case '|':
                    foreach (var h in ix.Results)
                        table.Put(h);
                    break;
                case '&':
                case '^':
                case 'E':
                case '!':
                    {
                        var nt = new EmbossHitTable(fdhits, true);
                        foreach (var h in ix.Results)
                            nt.Put(h);
                        if (field.Link == '&')
                            table.MergeAnd(nt);
                        else if (field.Link == '!')
                            table.MergeNot(nt);
                        else
                            table.MergeEor(nt);
                        break;
                    }
            }

            lasthits = fdhits;
        }

        // ajListSortUnique(textCdEntryCmp): by file then by offset, without duplicates
        var all = table.Values();
        var sorted = all.OrderBy(h => h.Dbno).ThenBy(h => h.Offset).ToList();
        ix.Results = [];
        EmbossBtHit? prev = null;
        foreach (var h in sorted)
        {
            if (prev is not null && prev.Dbno == h.Dbno && prev.Offset == h.Offset)
                continue;
            ix.Results.Add(h);
            prev = h;
        }

        return ix.Results.Count > 0;
    }

    /// <summary><c>textAccessEmblcd</c> for the whole query.</summary>
    public static EmbossQueryParse.ReadResult Access(EmbossQuery qry, EmbossQueryEnv env, EmbossFeatureReadContext context)
    {
        var fail = new EmbossQueryParse.ReadResult(false, null, "", false);
        if (qry.IndexDir.Length == 0)
        {
            context.Error(qry.QueryType == EmbossQueryType.All ? $"no indexdir defined for database {qry.DbName}" : $"no indexdir defined for database '{qry.DbName}'");
            if (qry.QueryType != EmbossQueryType.All)
                context.Warn($"Failed to open index for database '{qry.DbName}'");
            return fail;
        }

        string Path(string name) => qry.IndexDir.Length == 0 ? name : $"{qry.IndexDir}{Slash}{name}";
        string? divText = env.Files(Path("division.lkp"));
        if (divText is null)
        {
            context.Warn($"Cannot open division file in '{qry.IndexDir}' for database '{qry.DbName}'");
            if (qry.QueryType != EmbossQueryType.All)
                context.Warn($"Failed to open index for database '{qry.DbName}'");
            return fail;
        }

        var division = new EmbossCdFile(Encoding.Latin1.GetBytes(divText));
        if (qry.QueryType == EmbossQueryType.All)
        {
            // textCdAll: every file of the database that is included, read one after the other
            var texts = new List<string>();
            for (uint i = 0; i < division.NRecords; i++)
            {
                string name = division.DivisionName(i, true);
                string full = name.Contains('/', StringComparison.Ordinal) || qry.Directory.Length == 0 ? name : qry.Directory.TrimEnd('/') + "/" + name;
                if (!Included(full, qry.Exclude, qry.Filename))
                    continue;
                if (env.Files(full) is { } t)
                    texts.Add(t);
            }

            if (texts.Count == 0)
                return fail;
            qry.QryDone = !qry.CaseId;
            return new EmbossQueryParse.ReadResult(true, texts[0], qry.Filename, false, qry.DbName) { Texts = texts };
        }

        var skip = new bool[Math.Max(1, (int)division.NRecords)];
        for (uint i = 0; i < division.NRecords; i++)
        {
            string name = division.DivisionName(i, true);
            string full = name.Contains('/', StringComparison.Ordinal) || qry.Directory.Length == 0 ? name : qry.Directory.TrimEnd('/') + "/" + name;
            if (!Included(full, qry.Exclude, qry.Filename))
                skip[i] = true;
        }

        string? entText = env.Files(Path("entrynam.idx"));
        if (entText is null)
        {
            context.Error($"Cannot open index file '{Path("entrynam.idx")}' for database '{qry.DbName}'");
            context.Warn($"Failed to open index for database '{qry.DbName}'");
            return fail;
        }

        var ix = new Index { Division = division, Entries = new EmbossCdFile(Encoding.Latin1.GetBytes(entText)), Skip = skip, Query = qry, Env = env };
        bool ok = qry.QueryType == EmbossQueryType.Entry ? QryEntry(ix) : QryQuery(ix);
        _ = ok;
        if (ix.Results.Count == 0)
            return fail;

        // textCdQryNext for every hit: the data file of the division is opened when it changes (a file that is missing ends the read)
        var hitTexts = new List<string>();
        string? current = null;
        uint currentDiv = 0;
        foreach (var h in ix.Results)
        {
            if (h.Dbno != currentDiv)
            {
                currentDiv = h.Dbno;
                string name = division.DivisionName(h.Dbno - 1, false);
                int blank = name.IndexOf(' ', StringComparison.Ordinal);
                string datfile = blank >= 0 ? name[..blank] : name;
                if (datfile.Length == 0)
                {
                    context.Error($"index division file error '{name}'");
                    break;
                }

                string full = datfile.Contains('/', StringComparison.Ordinal) || qry.Directory.Length == 0 ? datfile : qry.Directory.TrimEnd('/') + "/" + datfile;
                current = env.Files(full);
                if (current is null)
                {
                    context.Error($"Cannot open file '{datfile}' for database '{qry.DbName}'");
                    break;
                }
            }

            hitTexts.Add(h.Offset < (ulong)current!.Length ? current[(int)h.Offset..] : "");
        }

        if (hitTexts.Count == 0)
            return fail;
        qry.QryDone = !qry.CaseId;
        return new EmbossQueryParse.ReadResult(true, hitTexts[0], qry.Filename, false, qry.DbName) { HitTexts = hitTexts };
    }

    private static bool Included(string path, string exclude, string include)
    {
        string name = path[(path.LastIndexOf('/') + 1)..];
        bool ret = false;
        char[] sep = [' ', '\t', ',', ';', '\n'];
        foreach (string token in include.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.Match(name, token))
                ret = true;
        }

        foreach (string token in exclude.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.Match(name, token))
                ret = false;
        }

        return ret;
    }
}
