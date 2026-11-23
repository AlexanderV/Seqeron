#pragma warning disable S107, S3776, S1121, S3358, S1066, S1854, S2325

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>An entry of a B+tree index (<c>AjPBtId</c>): the identifier, the number of the database file, how many entries have the identifier (<c>dups</c>: the entries are then in a secondary tree) and the offset in the file.</summary>
internal sealed class EmbossBtId
{
    public string Id = "";
    public uint Dbno, Dups;
    public ulong Offset;
}

/// <summary>A hit of a query (<c>AjPBtHit</c>): the file of the database and the offset of the entry in it.</summary>
internal sealed class EmbossBtHit(uint dbno, ulong offset)
{
    public uint Dbno { get; } = dbno;

    public ulong Offset { get; } = offset;

    /// <summary><c>ajBtreeHitCmp</c>.</summary>
    public static int Compare(EmbossBtHit x, EmbossBtHit y)
    {
        if (x.Dbno != y.Dbno)
            return x.Dbno > y.Dbno ? 1 : -1;
        if (x.Offset != y.Offset)
            return x.Offset > y.Offset ? 1 : -1;
        return 0;
    }

    /// <summary><c>ajBtreeHitHash</c>.</summary>
    public static ulong Hash(EmbossBtHit hit, ulong hashsize) => hashsize == 0 ? 0 : ((hit.Offset + hit.Dbno) >> 2) % hashsize;
}

/// <summary>A keyword of a secondary index (<c>AjPBtPri</c>): the word and the page of the secondary tree that lists its identifiers.</summary>
internal sealed class EmbossBtPri
{
    public string Keyword = "";
    public ulong Treeblock;
}

/// <summary>
/// One B+tree index file of EMBOSS 6.6.0 (<c>ajindex.c</c>: <c>ajBtreeCacheNewReadC</c> with <c>ajBtreeReadParamsC</c>, the page reads and the searches of the read-only index): the file as its bytes (latin-1 text) and its parameter file
/// (<c>NAME.pEXT</c>: pages of <c>Pagesize</c> bytes, <c>Order</c>, <c>Fill</c>, <c>Level</c>, <c>Count</c>, ...). The pages are read as they are in the file (a compressed index keeps its pages tightly packed and its page pointers are the offsets in the
/// compressed file, so the same reads work). A5-8c3h, B09 F275.
/// </summary>
internal sealed class EmbossBtreeCache
{
    private const uint BtLeaf = 4, BtIdBucket = 8, BtPriBucket = 32, BtSecBucket = 64, BtNumBucket = 128, BtSecLeaf = 4096;
    private const int NodeKeyLen = 52, BucketKeyLen = 16;

    private byte[] _data = [];

    // ajFatal: the library stops with the line of its source
    private static EmbossFatalException Fatal(int line, string message) => new($"\n   EMBOSS An error in ajindex.c at line {line}:\n{message}");

    public string FileName = "";
    public bool Secondary, Compressed;
    public uint KwLimit, IdLimit, RefCount, PriPageSize, SecPageSize, Order, Fill, Level, SOrder, SFill;
    public ulong CountUnique, CountAll;
    public ulong SecRootBlock;
    public uint SLevel;

    /// <summary>The extension of the index file of a field (<c>ajBtreeFieldGetExtensionC</c>: <c>btreeFields</c>, a field that is not in the table gets <c>x</c> + its name).</summary>
    public static string FieldExtension(string field) => field switch
    {
        "acc" => "xac", "id" => "xid", "sv" => "xsv", "up" => "xup", "spc" => "xspc", "nam" => "xnm", "key" => "xkw", "des" => "xde", "org" => "xtx", "isa" => "xis", "xref" => "xrf", "ns" => "xns", _ => "x" + field,
    };

    /// <summary><c>ajBtreeFieldGetSecondaryC</c>.</summary>
    public static bool FieldIsSecondary(string extension) => extension switch { "acc" or "id" or "sv" or "up" or "spc" => false, _ => true };

    private static string IndexName(string file, string ext, string idir, string prefix)
    {
        if (idir.Length == 0)
            return $"{file}.{prefix}{ext}";
        return idir.EndsWith('/') ? $"{idir}{file}.{prefix}{ext}" : $"{idir}/{file}.{prefix}{ext}";
    }

    // ajFmtScanS(line, "%*s%u", &value): the value stays as it was when the line has no number
    private static void ScanUint(string line, out uint value)
    {
        value = 0;
        var parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length > 1)
            _ = uint.TryParse(parts[1], out value);
    }

    private static void ScanUlong(string line, out ulong value)
    {
        value = 0;
        var parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length > 1)
            _ = ulong.TryParse(parts[1], out value);
    }

    /// <summary><c>ajBtreeCacheNewReadC</c>: opens the index of a field; null when the parameter file or the index file cannot be read (the library prints "Cannot open param file ...").</summary>
    // the bytes of an index file are made once per text (the data of the installation are the same strings on every call)
    private static readonly System.Runtime.CompilerServices.ConditionalWeakTable<string, byte[]> Bytes = new();

    private static byte[] BytesOf(string text) => Bytes.GetValue(text, t => Encoding.Latin1.GetBytes(t));

    public static EmbossBtreeCache? OpenRead(string file, string ext, string idir, Func<string, string?> files, EmbossFeatureReadContext context)
    {
        string pname = IndexName(file, ext, idir, "p");
        string? ptext = files(pname);
        if (ptext is null)
        {
            context.Error($"Cannot open param file {pname} 2: 'No such file or directory'\n");
            return null;
        }

        var c = new EmbossBtreeCache { FileName = IndexName(file, ext, idir, "") };
        bool secondary = true, setsecond = false;
        ulong pripagecount = 0;
        foreach (string raw in ptext.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (line.StartsWith("Order2 ", StringComparison.Ordinal))
                ScanUint(line, out c.SOrder);
            else if (line.StartsWith("Fill2 ", StringComparison.Ordinal))
                ScanUint(line, out c.SFill);
            if (line.StartsWith("Type ", StringComparison.Ordinal))
            {
                string type = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries).ElementAtOrDefault(1) ?? "";
                if (type == "Secondary")
                    secondary = true;
                else if (type == "Identifier")
                    secondary = false;
                else
                    context.Error($"param file {pname} unknown Type value '{type}'");
                setsecond = true;
            }

            if (line.StartsWith("Compress ", StringComparison.Ordinal))
            {
                string v = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries).ElementAtOrDefault(1) ?? "";
                c.Compressed = v.Length > 0 && char.ToUpperInvariant(v[0]) is 'Y' or 'T' or '1';
            }

            if (line.StartsWith("Pages ", StringComparison.Ordinal))
                ScanUlong(line, out pripagecount);
            if (line.StartsWith("Order ", StringComparison.Ordinal))
                ScanUint(line, out c.Order);
            if (line.StartsWith("Fill ", StringComparison.Ordinal))
                ScanUint(line, out c.Fill);
            if (line.StartsWith("Pagesize ", StringComparison.Ordinal))
                ScanUint(line, out c.PriPageSize);
            if (line.StartsWith("Secpagesize ", StringComparison.Ordinal))
                ScanUint(line, out c.SecPageSize);
            if (line.StartsWith("Level ", StringComparison.Ordinal))
                ScanUint(line, out c.Level);
            if (line.StartsWith("Count ", StringComparison.Ordinal))
                ScanUlong(line, out c.CountUnique);
            if (line.StartsWith("Fullcount ", StringComparison.Ordinal))
                ScanUlong(line, out c.CountAll);
            if (line.StartsWith("Kwlimit ", StringComparison.Ordinal))
                ScanUint(line, out c.KwLimit);
            if (line.StartsWith("Idlimit ", StringComparison.Ordinal))
                ScanUint(line, out c.IdLimit);
            if (line.StartsWith("Reffiles ", StringComparison.Ordinal))
                ScanUint(line, out c.RefCount);
        }

        if (!setsecond)
            secondary = FieldIsSecondary(ext);
        if (c.SecPageSize == 0)
            c.SecPageSize = c.PriPageSize;
        if (secondary && c.IdLimit == 0)
            c.IdLimit = c.KwLimit;
        c.Secondary = secondary;
        _ = pripagecount;
        if (c.PriPageSize == 0 || c.Order < 4 || c.Fill < 4 || c.SOrder < 4 || c.SFill < 4)
        {
            // ajBtreeIdcacheNewC: "cache '%S' pagesize %u order %u too small, increase pagesize"
            if (c.Order < 4)
                context.Error($"cache '{c.FileName}' pagesize {c.PriPageSize} order {c.Order} too small, increase pagesize");
            if (c.Fill < 4)
                context.Error($"cache '{c.FileName}' pagesize {c.PriPageSize} fill {c.Fill} too small, increase pagesize");
            if (c.SOrder < 4)
                context.Error($"cache '{c.FileName}' pagesize {c.PriPageSize} srder {c.SOrder} too small, increase pagesize");
            if (c.SFill < 4)
                context.Error($"cache '{c.FileName}' pagesize {c.SecPageSize} sfill {c.SFill} too small, increase pagesize");
            return null;
        }

        string? text = files(c.FileName);
        if (text is null)
            return null;
        c._data = BytesOf(text);
        return c;
    }

    /// <summary>The page at an offset, as the library reads it (the primary page size; the bytes after the end of the file are zero).</summary>
    private byte[] Page(ulong pos, uint size)
    {
        var buf = new byte[size];
        if (pos < (ulong)_data.Length)
        {
            int n = (int)Math.Min((ulong)size, (ulong)_data.Length - pos);
            Array.Copy(_data, (long)pos, buf, 0, n);
        }

        return buf;
    }

    /// <summary><c>btreePricacheRead</c>.</summary>
    private byte[] PriPage(ulong pos) => Page(pos, PriPageSize);

    /// <summary><c>btreeSeccacheRead</c>.</summary>
    private byte[] SecPage(ulong pos) => Page(pos, SecPageSize);

    private static uint U32(byte[] b, int o) => BitConverter.ToUInt32(b, o);

    private static ulong U64(byte[] b, int o) => BitConverter.ToUInt64(b, o);

    private static uint NodeType(byte[] b) => U32(b, 0);

    private static uint NKeys(byte[] b) => U32(b, 12);

    private static ulong Right(byte[] b) => U64(b, 28);

    /// <summary>The C string at the offset (up to the first zero byte).</summary>
    private static string CStr(byte[] b, int o, int max = -1)
    {
        int e = o;
        int limit = max >= 0 ? Math.Min(b.Length, o + max) : b.Length;
        while (e < limit && b[e] != 0)
            e++;
        return Encoding.Latin1.GetString(b, o, e - o);
    }

    // btreeGetBlockC: the pointer of the first key that is greater than the key (strcmp), else the last pointer
    private ulong GetBlock(byte[] buf, string key)
    {
        uint m = NKeys(buf);
        if (m == 0)
            throw Fatal(6950, $"btreeGetBlockC: No keys in node cache {FileName}");
        int lenptr = NodeKeyLen;
        int keyptr = lenptr + (int)m * 4;
        for (uint i = 0; i < m; ++i)
        {
            uint len = U32(buf, lenptr) + 1;
            string k = CStr(buf, keyptr);
            if (string.CompareOrdinal(k, key) > 0)
                return U64(buf, keyptr + (int)len);
            keyptr += (int)len + 8;
            lenptr += 4;
        }

        return U64(buf, keyptr);
    }

    // btreeGetPointers
    private ulong[] GetPointers(byte[] buf)
    {
        uint nkeys = NKeys(buf);
        if (nkeys == 0)
            throw Fatal(6844, $"btreeGetPointers: No keys in node cache {FileName}");
        var p = new ulong[nkeys + 1];
        int lenptr = NodeKeyLen;
        int keyptr = lenptr + (int)nkeys * 4;
        for (uint i = 0; i < nkeys; ++i)
        {
            uint len = U32(buf, lenptr) + 1;
            keyptr += (int)len;
            p[i] = U64(buf, keyptr);
            keyptr += 8;
            lenptr += 4;
        }

        p[nkeys] = U64(buf, keyptr);
        return p;
    }

    // btreeGetNumPointers
    private ulong[] GetNumPointers(byte[] buf)
    {
        uint nkeys = NKeys(buf);
        if (nkeys == 0)
            throw Fatal(21628, $"GetNumKeys: No keys in node cache {FileName}");
        var p = new ulong[nkeys + 1];
        int pptr = NodeKeyLen + (int)nkeys * 8;
        for (uint i = 0; i <= nkeys; ++i)
        {
            p[i] = U64(buf, pptr);
            pptr += 8;
        }

        return p;
    }

    /// <summary><c>ajStrFmtQuery</c>: lower case; any run of characters that are not letters, digits, <c>*</c>, <c>?</c> or <c>.</c> becomes one underscore.</summary>
    public static string FmtQuery(string s)
    {
        var sb = new StringBuilder();
        int i = 0;
        while (i < s.Length)
        {
            if (!IsKeep(s[i]))
            {
                while (i < s.Length && !IsKeep(s[i]))
                    i++;
                sb.Append('_');
                continue;
            }

            sb.Append(char.ToLowerInvariant(s[i]));
            i++;
        }

        return sb.ToString();
    }

    private static bool IsKeep(char c) => c is (>= 'a' and <= 'z') or (>= 'A' and <= 'Z') or (>= '0' and <= '9') or '*' or '?' or '.';

    // btreeIdentFind: the leaf that holds the key
    private byte[] IdentFind(string key)
    {
        var root = PriPage(0);
        if (Level == 0)
            return root;
        return FindInode(root, key, false);
    }

    // btreePrimaryFindInode / btreePrimaryFindInodeWild (the same descent)
    private byte[] FindInode(byte[] page, string key, bool wild)
    {
        _ = wild;
        var ret = page;
        while (NodeType(ret) != BtLeaf)
            ret = PriPage(GetBlock(ret, key));
        return ret;
    }

    // btreeReadIdbucket (all = false) and btreeIdbucketIdlistAll (all = true): the entries of an id bucket
    private List<EmbossBtId> ReadIdbucket(ulong pagepos, out ulong overflow, bool all = false)
    {
        var list = new List<EmbossBtId>();
        overflow = 0;
        string who = all ? "IdbucketIdlistAll" : "ReadBucket";
        if (pagepos == 0)
            throw all ? Fatal(21449, $"IdbucketIdlistAll: cannot read bucket from root page cache {FileName}") : Fatal(4846, "BucketRead: cannot read bucket from root page");
        var buf = PriPage(pagepos);
        uint nodetype = U32(buf, 0), nentries = U32(buf, 4);
        overflow = U64(buf, 8);
        if (nodetype != BtIdBucket)
            throw all ? Fatal(21462, $"IdbucketIdlistAll: NodeType mismatch. Not bucket ({nodetype}) cache {FileName}") : Fatal(4859, $"ReadBucket: NodeType mismatch. Not bucket ({nodetype})");
        if (nentries > Fill)
            throw Fatal(all ? 21468 : 4864, $"{who}: Bucket too full page: {pagepos} entries: {nentries} max: {Fill} cache {FileName}");
        uint keyskip = 4 + 4 + 8 + RefCount * 8 + 1;
        int kptr = BucketKeyLen;
        int idptr = kptr + (int)nentries * 4;
        for (uint i = 0; i < nentries; ++i)
        {
            uint len = U32(buf, kptr);
            int idlen = (int)(len - keyskip);
            var id = new EmbossBtId { Id = Encoding.Latin1.GetString(buf, idptr, Math.Max(0, idlen)) };
            idptr += idlen + 1;
            id.Dbno = U32(buf, idptr);
            idptr += 4;
            id.Dups = U32(buf, idptr);
            idptr += 4;
            id.Offset = U64(buf, idptr);
            idptr += 8 + (int)RefCount * 8;
            kptr += 4;
            list.Add(id);
        }

        return list;
    }

    // btreeNumbucketBthitlist / btreeNumbucketBtidlist: the entries (offset, reference offsets, file) of a bucket of a secondary tree of numbers
    private void NumBucket(ulong pagepos, string? idname, List<EmbossBtHit>? hits, List<EmbossBtId>? ids)
    {
        bool asHits = hits is not null;
        string who = asHits ? "NumbucketBthitlist" : "NumbucketBtidlist";
        if (pagepos == SecRootBlock)
            throw Fatal(asHits ? 22175 : 22080, $"{who}: cannot read bucket from root page cache {FileName}");
        var buf = SecPage(pagepos);
        uint nodetype = U32(buf, 0), nentries = U32(buf, 4);
        if (nodetype != BtNumBucket && nodetype != BtIdBucket)
            throw Fatal(asHits ? 22187 : 22092, $"{who}: Nodetype mismatch. Not bucket ({nodetype}) cache {FileName}");
        if (nentries > SFill)
            throw Fatal(asHits ? 22193 : 22098, $"{who}: Bucket too full page: {pagepos} entries: {nentries} max: {SFill} cache {FileName}");
        int pptr = BucketKeyLen;
        for (uint i = 0; i < nentries; ++i)
        {
            ulong offset = U64(buf, pptr);
            pptr += 8 + (int)RefCount * 8;
            uint dbno = U32(buf, pptr);
            pptr += 4;
            hits?.Add(new EmbossBtHit(dbno, offset));
            ids?.Add(new EmbossBtId { Id = idname ?? "", Dbno = dbno, Offset = offset });
        }
    }

    // btreeIdentFetchMultiHit / btreeIdentFetchMulti: all the entries of an identifier that occurs several times
    private void FetchMulti(ulong rootblock, string? idname, List<EmbossBtHit>? hits, List<EmbossBtId>? ids)
    {
        var page = SecPage(rootblock);
        SLevel = (uint)Right(page);
        var parray = GetNumPointers(page);
        uint nodetype = NodeType(page);
        while (nodetype != BtSecLeaf && SLevel != 0)
        {
            page = SecPage(parray[0]);
            parray = GetNumPointers(page);
            nodetype = NodeType(page);
        }

        ulong right;
        do
        {
            uint nkeys = NKeys(page);
            for (uint i = 0; i <= nkeys; ++i)
                NumBucket(parray[i], idname, hits, ids);
            right = 0;
            if (SLevel != 0)
            {
                right = Right(page);
                if (right != 0)
                {
                    page = SecPage(right);
                    parray = GetNumPointers(page);
                }
            }
        }
        while (right != 0);
    }

    /// <summary><c>btreeIdentQueryHit</c>: the entry of an identifier (the first of its duplicates).</summary>
    public EmbossBtHit? IdentQueryHit(string key)
    {
        if (CountUnique == 0 && CountAll == 0)
            return null;
        string keystr = FmtQuery(key);
        var page = IdentFind(keystr);
        ulong blockno = GetBlock(page, keystr);
        var ids = ReadIdbucket(blockno, out _);
        foreach (var id in ids)
        {
            if (id.Id == keystr)
                return new EmbossBtHit(id.Dbno, id.Offset);
        }

        return null;
    }

    /// <summary><c>btreeIdentQueryId</c>.</summary>
    private EmbossBtId? IdentQueryId(string key)
    {
        if (CountUnique == 0 && CountAll == 0)
            return null;
        string keystr = FmtQuery(key);
        var page = IdentFind(keystr);
        ulong blockno = GetBlock(page, keystr);
        foreach (var id in ReadIdbucket(blockno, out _))
        {
            if (id.Id == keystr)
                return id;
        }

        return null;
    }

    /// <summary><c>ajBtreeIdentFetchHit</c>: the hits of an identifier (appended to the list).</summary>
    public void IdentFetchHit(string key, List<EmbossBtHit> hitlist, EmbossFeatureReadContext context)
    {
        if (CountUnique == 0 && CountAll == 0)
            return;
        if (RefCount != 0)
            context.Warn($"ajBtreeIdentFetchHit called for cache '{FileName}' with {RefCount} reference files");
        string keystr = FmtQuery(key);
        var page = IdentFind(keystr);
        ulong blockno = GetBlock(page, keystr);
        var ids = ReadIdbucket(blockno, out _);
        EmbossBtId? tid = ids.FirstOrDefault(i => i.Id == keystr);
        if (tid is null)
            return;
        if (tid.Dups == 0)
            hitlist.Add(new EmbossBtHit(tid.Dbno, tid.Offset));
        else
            FetchMulti(tid.Offset, null, hitlist, null);
    }

    // btreeIdbucketIdlistAll: the entries of an id bucket (an identifier with duplicates is replaced by its entries)
    private void IdbucketIdlistAll(ulong pagepos, List<EmbossBtId> idlist)
    {
        foreach (var id in ReadIdbucket(pagepos, out _, true))
        {
            if (id.Dups == 0)
                idlist.Add(id);
            else
                FetchMulti(id.Offset, id.Id, null, idlist);
        }
    }

    // btreeIdleafFetch: the entries of a leaf, sorted by identifier
    private void IdleafFetch(byte[] page, List<EmbossBtId> list)
    {
        uint nkeys = NKeys(page);
        if (nkeys == 0)
            return;
        var parray = GetPointers(page);
        for (uint i = 0; i < nkeys + 1; ++i)
            IdbucketIdlistAll(parray[i], list);
        var sorted = list.OrderBy(x => x.Id, StringComparer.Ordinal).ToList();
        list.Clear();
        list.AddRange(sorted);
    }

    // btreePribucketIdlist
    private void PribucketIdlist(ulong pagepos, List<EmbossBtPri> list)
    {
        if (pagepos == 0)
            throw Fatal(11817, $"PribucketIdlist: cannot read bucket from root page cache {FileName}");
        var buf = PriPage(pagepos);
        uint nodetype = U32(buf, 0), nentries = U32(buf, 4);
        if (nodetype != BtPriBucket)
            throw Fatal(11830, $"PribucketIdlist: NodeType mismatch. Not primary bucket ({nodetype}) cache {FileName}");
        int kptr = BucketKeyLen;
        int code = kptr + (int)nentries * 4;
        for (uint i = 0; i < nentries; ++i)
        {
            uint len = U32(buf, kptr);
            int idlen = (int)len - 8 - 1;
            var pri = new EmbossBtPri { Keyword = Encoding.Latin1.GetString(buf, code, Math.Max(0, idlen)) };
            code += idlen + 1;
            pri.Treeblock = U64(buf, code);
            code += 8;
            kptr += 4;
            list.Add(pri);
        }
    }

    // btreePrileafFetch
    private void PrileafFetch(byte[] page, List<EmbossBtPri> list)
    {
        var parray = GetPointers(page);
        uint nkeys = NKeys(page);
        for (uint i = 0; i < nkeys + 1; ++i)
            PribucketIdlist(parray[i], list);
        var sorted = list.OrderBy(x => x.Keyword, StringComparer.Ordinal).ToList();
        list.Clear();
        list.AddRange(sorted);
    }

    // btreePribucketFindId
    private bool PribucketFindId(ulong pagepos, string id, out ulong treeblock)
    {
        treeblock = 0;
        if (pagepos == 0)
            throw Fatal(13021, $"PribucketFindId: cannot read bucket from root page cache {FileName}");
        var buf = PriPage(pagepos);
        uint nodetype = U32(buf, 0), nentries = U32(buf, 4);
        if (nodetype != BtPriBucket)
            throw Fatal(13034, $"PribucketFindId: NodeType mismatch. Not primary bucket ({nodetype}) cache {FileName}");
        if (nentries > Fill)
            throw Fatal(13040, $"PribucketFindId: Bucket too full  page: {pagepos} entries: {nentries} max: {Fill} cache {FileName}");
        int kptr = BucketKeyLen;
        int code = kptr + (int)nentries * 4;
        for (uint i = 0; i < nentries; ++i)
        {
            uint len = U32(buf, kptr);
            int idlen = (int)len - 8 - 1;
            if (CStr(buf, code) == id)
            {
                treeblock = U64(buf, code + idlen + 1);
                return true;
            }

            code += idlen + 1 + 8;
            kptr += 4;
        }

        return false;
    }

    // btreeSecbucketIdlist: the identifiers of a secondary bucket
    private void SecbucketIdlist(ulong pagepos, List<string> idlist)
    {
        var buf = SecPage(pagepos);
        uint nodetype = U32(buf, 0), nentries = U32(buf, 4);
        if (nodetype != BtSecBucket)
            throw Fatal(13642, $"SecbucketIdlist: NodeType mismatch. Not secondary bucket ({nodetype}) cache {FileName}");
        if (nentries > SFill)
            throw Fatal(13648, $"SecbucketIdlist: Bucket too full page: {pagepos} entries: {nentries} max: {SFill} full cache {FileName}");
        int kptr = BucketKeyLen;
        int code = kptr + (int)nentries * 4;
        for (uint i = 0; i < nentries; ++i)
        {
            uint len = U32(buf, kptr);
            idlist.Add(Encoding.Latin1.GetString(buf, code, (int)len - 1));
            code += (int)len;
            kptr += 4;
        }
    }

    // btreeSecTreeList: the identifiers of the secondary tree at a page
    private List<string> SecTreeList(ulong rootblock)
    {
        var list = new List<string>();
        var page = SecPage(rootblock);
        uint nodetype = NodeType(page);
        if (nodetype == BtSecBucket)
        {
            SecbucketIdlist(rootblock, list);
            return list;
        }

        SLevel = (uint)Right(page);
        var parray = GetPointers(page);
        while (nodetype != BtSecLeaf && SLevel != 0)
        {
            page = SecPage(parray[0]);
            parray = GetPointers(page);
            nodetype = NodeType(page);
        }

        ulong right;
        do
        {
            uint nkeys = NKeys(page);
            for (uint j = 0; j <= nkeys; ++j)
                SecbucketIdlist(parray[j], list);
            right = 0;
            if (SLevel != 0)
            {
                right = Right(page);
                if (right != 0)
                {
                    page = SecPage(right);
                    parray = GetPointers(page);
                }
            }
        }
        while (right != 0);
        return list;
    }

    // btreeReadAllSecLeaves
    private void ReadAllSecLeaves(List<string> list)
    {
        var root = SecPage(SecRootBlock);
        uint nodetype = NodeType(root);
        if (nodetype == BtSecBucket)
        {
            SecbucketIdlist(SecRootBlock, list);
            return;
        }

        ulong right = Right(root);
        SLevel = (uint)right;
        var buf = root;
        ulong secpagepos;
        if (SLevel != 0)
        {
            while (nodetype != BtSecLeaf)
            {
                var parr = GetPointers(buf);
                buf = SecPage(parr[0]);
                nodetype = NodeType(buf);
            }

            secpagepos = GetPointersFirstOfParent(root, buf);
        }
        else
        {
            secpagepos = SecRootBlock;
        }

        var parray = GetPointers(buf);
        uint nkeys = NKeys(buf);
        for (uint i = 0; i < nkeys + 1; ++i)
            SecbucketIdlist(parray[i], list);
        if (SOrder != 0)
        {
            buf = SecPage(secpagepos);
            right = Right(buf);
        }

        while (right != 0 && secpagepos != SecRootBlock)
        {
            secpagepos = right;
            buf = SecPage(secpagepos);
            parray = GetPointers(buf);
            nkeys = NKeys(buf);
            for (uint i = 0; i < nkeys + 1; ++i)
                SecbucketIdlist(parray[i], list);
            right = Right(buf);
        }
    }

    // the pointer parray[0] that btreeReadAllSecLeaves leaves behind after the descent (the page of the first leaf)
    private ulong GetPointersFirstOfParent(byte[] root, byte[] leaf)
    {
        // the descent reads parray[0] of every level; its last value is the page of the leaf
        var buf = root;
        ulong pos = 0;
        while (NodeType(buf) != BtSecLeaf)
        {
            pos = GetPointers(buf)[0];
            buf = SecPage(pos);
        }

        _ = leaf;
        return pos;
    }

    /// <summary><c>ajBtreeKeyFindLen</c>: the page of the secondary tree of a keyword.</summary>
    private bool KeyFindLen(string key, out ulong treeblock)
    {
        treeblock = 0;
        if (CountUnique == 0)
            return false;
        string keystr = FmtQuery(key);
        if (keystr.Length > KwLimit)
            keystr = keystr[..(int)KwLimit];
        var page = IdentFind(keystr);
        ulong blockno = GetBlock(page, keystr);
        return PribucketFindId(blockno, keystr, out treeblock);
    }

    /// <summary><c>ajBtreeKeyFetchHit</c>: the hits of the identifiers of a keyword (<paramref name="idcache"/> is the index of identifiers).</summary>
    public void KeyFetchHit(EmbossBtreeCache idcache, string key, List<EmbossBtHit> hitlist)
    {
        if (!KeyFindLen(key, out ulong treeblock))
            return;
        foreach (string kwid in SecTreeList(treeblock))
        {
            var id = idcache.IdentQueryId(kwid);
            if (id is null)
                continue;
            if (id.Dups == 0)
                hitlist.Add(new EmbossBtHit(id.Dbno, id.Offset));
            else
                idcache.FetchMulti(id.Offset, null, hitlist, null);
        }
    }

    /// <summary><c>btreeKeyFullSearchHit</c>: every entry whose identifier matches the wildcard (the pattern starts with a wildcard).</summary>
    private void KeyFullSearchHit(string key, List<EmbossBtHit> hitlist)
    {
        var page = PriPage(0);
        if (Level != 0)
        {
            while (NodeType(page) != BtLeaf)
                page = PriPage(GetPointers(page)[0]);
        }

        ulong right = 99;
        while (right != 0)
        {
            var list = new List<EmbossBtId>();
            IdleafFetch(page, list);
            foreach (var id in list)
            {
                if (EmbossWild.Match(id.Id, key))
                    hitlist.Add(new EmbossBtHit(id.Dbno, id.Offset));
            }

            right = Right(page);
            if (right != 0)
                page = PriPage(right);
        }

        SortTwoUnique(hitlist);
    }

    /// <summary><c>ajListSortTwoUnique</c> with the file and the offset: sorted, without duplicates.</summary>
    public static void SortTwoUnique(List<EmbossBtHit> list)
    {
        if (list.Count <= 1)
            return;
        var sorted = list.OrderBy(h => h.Dbno).ThenBy(h => h.Offset).ToList();
        list.Clear();
        EmbossBtHit? prev = null;
        foreach (var h in sorted)
        {
            if (prev is not null && EmbossBtHit.Compare(h, prev) == 0)
                continue;
            list.Add(h);
            prev = h;
        }
    }

    /// <summary><c>ajBtreeIdentFetchwildHit</c>: the hits of the identifiers that match the wildcard.</summary>
    public void IdentFetchwildHit(string key, List<EmbossBtHit> hitlist)
    {
        if (CountUnique == 0 && CountAll == 0)
            return;
        string keystr = FmtQuery(key);
        if (keystr.Length > KwLimit)
            keystr = keystr[..(int)KwLimit];
        int p = keystr.IndexOfAny(['*', '?']);
        string prefix;
        if (p >= 0)
        {
            if (p > 0)
            {
                prefix = keystr[..p];
            }
            else
            {
                KeyFullSearchHit(keystr, hitlist);
                return;
            }
        }
        else
        {
            prefix = keystr;
        }

        prefix = FmtQuery(prefix);
        if (prefix.Length > KwLimit)
            prefix = prefix[..(int)KwLimit];
        var page = FindLeafWild(prefix);
        ulong pripagepos = PagePos(page);
        var list = new List<EmbossBtId>();
        IdleafFetch(page, list);
        var queue = new Queue<EmbossBtId>(list);
        EmbossBtId? id = null;
        bool found = false;
        while (queue.Count > 0)
        {
            var cand = queue.Dequeue();
            if (cand.Id.StartsWith(prefix, StringComparison.Ordinal))
            {
                id = cand;
                found = true;
                break;
            }
        }

        if (!found)
        {
            ulong right = Right(page);
            if (right == 0)
                return;
            page = PriPage(right);
            pripagepos = right;
            list = [];
            IdleafFetch(page, list);
            if (list.Count == 0)
                return;
            queue = new Queue<EmbossBtId>(list);
            found = false;
            while (queue.Count > 0)
            {
                var cand = queue.Dequeue();
                if (cand.Id.StartsWith(prefix, StringComparison.Ordinal))
                {
                    id = cand;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
            return;
        var pushed = new List<EmbossBtHit>();                           // ajListPush: every hit goes to the start of the list
        bool finished = false;
        while (!finished)
        {
            if (EmbossWild.Match(id!.Id, keystr))
                pushed.Insert(0, new EmbossBtHit(id.Dbno, id.Offset));
            if (queue.Count == 0)
            {
                var cur = PriPage(pripagepos);
                ulong right = Right(cur);
                if (right == 0)
                {
                    finished = true;
                    continue;
                }

                page = PriPage(right);
                pripagepos = right;
                list = [];
                IdleafFetch(page, list);
                queue = new Queue<EmbossBtId>(list);
                if (queue.Count == 0)
                {
                    finished = true;
                    continue;
                }
            }

            id = queue.Dequeue();
            if (!id.Id.StartsWith(prefix, StringComparison.Ordinal))
                finished = true;
        }

        // the hits are put at the start of the list of the query (ajListPush), then the list is sorted and made unique
        hitlist.InsertRange(0, pushed);
        SortTwoUnique(hitlist);
    }

    private ulong _lastLeafPos;

    private ulong PagePos(byte[] page) => _lastLeafPos;

    // btreePrimaryFetchFindleafWild
    private byte[] FindLeafWild(string key)
    {
        var root = PriPage(0);
        _lastLeafPos = 0;
        if (Level == 0)
            return root;
        var ret = root;
        while (NodeType(ret) != BtLeaf)
        {
            ulong blockno = GetBlock(ret, key);
            ret = PriPage(blockno);
            _lastLeafPos = blockno;
        }

        return ret;
    }

    /// <summary><c>ajBtreeKeyFetchwildHit</c>: the hits of the identifiers of the keywords that match the wildcard.</summary>
    public void KeyFetchwildHit(EmbossBtreeCache idcache, string key, List<EmbossBtHit> hitlist)
    {
        int p = key.IndexOfAny(['*', '?']);
        string prefix;
        if (p >= 0)
        {
            if (p > 0)
            {
                prefix = key[..p];
            }
            else
            {
                KeywordFullSearchHit(key, idcache, hitlist);
                return;
            }
        }
        else
        {
            prefix = key;
        }

        var page = FindLeafWild(prefix);
        ulong pripagepossave = PagePos(page);
        var prilist = new List<EmbossBtPri>();
        PrileafFetch(page, prilist);
        if (prilist.Count == 0)
            return;
        var queue = new Queue<EmbossBtPri>(prilist);
        EmbossBtPri? pri = null;
        bool found = false;
        while (queue.Count > 0)
        {
            var cand = queue.Dequeue();
            if (cand.Keyword.StartsWith(prefix, StringComparison.Ordinal))
            {
                pri = cand;
                found = true;
                break;
            }
        }

        if (!found)
        {
            ulong pagepos = Right(page);
            if (pagepos == 0)
                return;
            page = PriPage(pagepos);
            pripagepossave = pagepos;
            prilist = [];
            PrileafFetch(page, prilist);
            if (prilist.Count == 0)
                return;
            queue = new Queue<EmbossBtPri>(prilist);
            while (queue.Count > 0)
            {
                var cand = queue.Dequeue();
                if (cand.Keyword.StartsWith(prefix, StringComparison.Ordinal))
                {
                    pri = cand;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
            return;
        var strlist = new List<string>();
        bool finished = false;
        while (!finished)
        {
            if (EmbossWild.Match(pri!.Keyword, key))
            {
                SecRootBlock = pri.Treeblock;
                ReadAllSecLeaves(strlist);
            }

            if (queue.Count == 0)
            {
                var cur = PriPage(pripagepossave);
                ulong right = Right(cur);
                if (right == 0)
                {
                    finished = true;
                    continue;
                }

                // the library reads the leaf of the position it left (the page is kept) and then moves to the right neighbour
                var leaf = PriPage(pripagepossave);
                prilist = [];
                PrileafFetch(leaf, prilist);
                pripagepossave = right;
                queue = new Queue<EmbossBtPri>(prilist);
                if (queue.Count == 0)
                {
                    finished = true;
                    continue;
                }
            }

            pri = queue.Dequeue();
            if (!pri.Keyword.StartsWith(prefix, StringComparison.Ordinal))
                finished = true;
        }

        if (strlist.Count > 0)
        {
            // ajListSortUnique(strlist, ajStrVcmp)
            foreach (string id in strlist.Distinct(StringComparer.Ordinal).OrderBy(x => x, StringComparer.Ordinal))
            {
                var btid = idcache.IdentQueryId(id);
                if (btid is not null)
                    hitlist.Add(new EmbossBtHit(btid.Dbno, btid.Offset));
            }
        }
    }

    /// <summary><c>btreeKeywordFullSearchHit</c>: every keyword that matches a wildcard that starts with one.</summary>
    private void KeywordFullSearchHit(string key, EmbossBtreeCache idcache, List<EmbossBtHit> hitlist)
    {
        var strtable = new EmbossStrTable(idcache.CountUnique);
        var page = PriPage(0);
        if (Level != 0)
        {
            while (NodeType(page) != BtLeaf)
                page = PriPage(GetPointers(page)[0]);
        }

        ulong right = 99;
        while (right != 0)
        {
            var list = new List<EmbossBtPri>();
            PrileafFetch(page, list);
            foreach (var pri in list)
            {
                if (!EmbossWild.Match(pri.Keyword, key))
                    continue;
                SecRootBlock = pri.Treeblock;
                var strlist = new List<string>();
                ReadAllSecLeaves(strlist);
                foreach (string id in strlist)
                    strtable.Put(id);
            }

            right = Right(page);
            if (right != 0)
                page = PriPage(right);
        }

        foreach (string id in strtable.Keys())
        {
            var hit = idcache.IdentQueryHit(id);
            if (hit is not null)
                hitlist.Add(hit);
        }
    }
}

/// <summary>The table of strings of <c>ajTablestrNew</c> (a chained hash table with the hash of <c>ajTablestrHash</c>): its keys come back in the order of its buckets.</summary>
internal sealed class EmbossStrTable
{
    private static readonly ulong[] Primes = [3, 7, 13, 31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573, 2097143, 4194301, 8388593, 16777213, 33554393, 67108859, 134217689, 268435399, 536870909, 1073741789, 2147483647, uint.MaxValue];

    private List<string>?[] _buckets;
    private ulong _size;
    private ulong _length;

    public EmbossStrTable(ulong size)
    {
        ulong hint = size / 4;
        int i = 1;
        while (Primes[i] <= hint)
            i++;
        _size = Primes[i - 1];
        _buckets = new List<string>?[_size];
    }

    private static ulong Hash(string key, ulong hashsize)
    {
        ulong hash = 0;
        foreach (char ch in key)
        {
            sbyte c = unchecked((sbyte)ch);
            hash = unchecked(hash * 127 + (ulong)(long)c) % hashsize;
        }

        return hash;
    }

    /// <summary><c>ajTablePut</c> (a key that is there stays where it is).</summary>
    public void Put(string key)
    {
        ulong minsize = _length / 8;
        if (_size < minsize)
            Resize(6 * minsize / 4);
        ulong i = Hash(key, _size);
        var b = _buckets[i] ??= [];
        if (b.Contains(key))
            return;
        b.Insert(0, key);
        _length++;
    }

    private void Resize(ulong hashsize)
    {
        ulong hint = hashsize;
        int i = 1;
        while (Primes[i] <= hint)
            i++;
        ulong iprime = Primes[i - 1];
        if (iprime == _size)
            return;
        var keys = Keys().ToList();
        _size = iprime;
        _buckets = new List<string>?[_size];
        foreach (string k in keys)
        {
            ulong ib = Hash(k, _size);
            (_buckets[ib] ??= []).Insert(0, k);
        }
    }

    /// <summary><c>ajTableToarrayKeys</c>: the keys in the order of the buckets.</summary>
    public IEnumerable<string> Keys()
    {
        foreach (var b in _buckets)
        {
            if (b is null)
                continue;
            foreach (string k in b)
                yield return k;
        }
    }
}
