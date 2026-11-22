#pragma warning disable S107, S3776, S1121, S3358, S1066, S1854, S2325

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The chained hash table of <c>ajtable.c</c> as <c>textEmbossQryQuery</c> uses it to combine the hits of the fields of a query (<c>ajTableNewFunctionLen</c>, <c>ajTablePutClean</c>, <c>ajTableResizeHashsize</c>, <c>ajTableMergeOr</c> / <c>And</c> / <c>Eor</c> / <c>Not</c>,
/// <c>ajTableToarrayValues</c>) with the hash and the comparison of <c>ajBtreeHitHash</c> / <c>ajBtreeHitCmp</c>: the hits come back in the order of the buckets, the newest of a bucket first. A5-8c3h, B09 F275.
/// </summary>
internal sealed class EmbossHitTable
{
    private static readonly ulong[] Primes = [3, 7, 13, 31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573, 2097143, 4194301, 8388593, 16777213, 33554393, 67108859, 134217689, 268435399, 536870909, 1073741789, 2147483647, uint.MaxValue];

    private List<EmbossBtHit>?[] _buckets;
    private readonly bool _ulongKeys;

    public ulong Size { get; private set; }

    public ulong Length { get; private set; }

    /// <summary><c>ajTableNewFunctionLen(size, ...)</c>: the number of buckets is the largest prime that is at most size / 4 (3 at least).</summary>
    /// <param name="size">The number of hits the table is made for.</param>
    /// <param name="ulongKeys">The table of <c>textCdQryQuery</c> (<c>ajTableulongNew</c>): the key of a hit is <c>(file &lt;&lt; 4) + offset</c> (sizeof(ajulong)/2 is the shift of the library), hashed as <c>key &gt;&gt; 2</c>, and two hits are the same when their keys are.</param>
    public EmbossHitTable(ulong size, bool ulongKeys = false)
    {
        _ulongKeys = ulongKeys;
        ulong hint = size / 4;
        int i = 1;
        while (Primes[i] <= hint)
            i++;
        Size = Primes[i - 1];
        _buckets = new List<EmbossBtHit>?[Size];
    }

    /// <summary><c>ajTablePutClean</c>: a hit that is in the table (same file and offset) replaces the one there, in its place.</summary>
    public void Put(EmbossBtHit hit)
    {
        ulong minsize = Length / 8;
        if (Size < minsize)
            ResizeHashsize(6 * minsize / 4);
        ulong i = HashOf(hit, Size);
        var b = _buckets[i] ??= [];
        for (int k = 0; k < b.Count; k++)
        {
            if (Same(hit, b[k]))
            {
                b[k] = hit;
                return;
            }
        }

        b.Insert(0, hit);
        Length++;
    }

    private ulong HashOf(EmbossBtHit h, ulong size) => _ulongKeys ? ((((ulong)h.Dbno << 4) + h.Offset) >> 2) % size : EmbossBtHit.Hash(h, size);

    private bool Same(EmbossBtHit a, EmbossBtHit b) => _ulongKeys ? (((ulong)a.Dbno << 4) + a.Offset) == (((ulong)b.Dbno << 4) + b.Offset) : EmbossBtHit.Compare(a, b) == 0;

    /// <summary><c>ajTableToarrayValues</c>.</summary>
    public List<EmbossBtHit> Values()
    {
        var r = new List<EmbossBtHit>();
        foreach (var b in _buckets)
        {
            if (b is not null)
                r.AddRange(b);
        }

        return r;
    }

    /// <summary><c>ajTableResizeHashsize</c>: the hits are put again, in the order of the old buckets, at the start of their new buckets.</summary>
    public void ResizeHashsize(ulong hashsize)
    {
        int i = 1;
        while (Primes[i] <= hashsize)
            i++;
        ulong iprime = Primes[i - 1];
        if (iprime == Size)
            return;
        var all = Values();
        Size = iprime;
        _buckets = new List<EmbossBtHit>?[Size];
        foreach (var h in all)
        {
            ulong ib = HashOf(h, Size);
            (_buckets[ib] ??= []).Insert(0, h);
        }
    }

    private static void MakeSameSize(EmbossHitTable a, EmbossHitTable b)
    {
        if (b.Size > a.Size)
            a.ResizeHashsize(b.Size);
        else if (a.Size > b.Size)
            b.ResizeHashsize(a.Size);
    }

    private int IndexOf(List<EmbossBtHit>? b, EmbossBtHit h)
    {
        if (b is null)
            return -1;
        for (int k = 0; k < b.Count; k++)
        {
            if (Same(h, b[k]))
                return k;
        }

        return -1;
    }

    // the nodes of table that go to the start of the bucket of this, one after the other (the last one of the bucket ends up first)
    private void Absorb(EmbossHitTable table, ulong i)
    {
        var tb = table._buckets[i];
        if (tb is not null)
        {
            foreach (var q in tb)
            {
                (_buckets[i] ??= []).Insert(0, q);
                table.Length--;
                Length++;
            }
        }

        table._buckets[i] = null;
    }

    /// <summary><c>ajTableMergeOr</c>: this gets the hits of the table that it has not (the table is emptied).</summary>
    public void MergeOr(EmbossHitTable table)
    {
        MakeSameSize(this, table);
        for (ulong i = 0; i < Size; i++)
        {
            var pb = _buckets[i];
            if (pb is not null)
            {
                foreach (var p in pb)
                {
                    int m = IndexOf(table._buckets[i], p);
                    if (m >= 0)
                    {
                        table._buckets[i]!.RemoveAt(m);
                        table.Length--;
                    }
                }
            }

            Absorb(table, i);
        }
    }

    /// <summary><c>ajTableMergeAnd</c>: this keeps the hits that the table has too (the table is emptied).</summary>
    public void MergeAnd(EmbossHitTable table)
    {
        MakeSameSize(this, table);
        for (ulong i = 0; i < Size; i++)
        {
            var pb = _buckets[i];
            if (pb is not null)
            {
                for (int k = 0; k < pb.Count;)
                {
                    if (IndexOf(table._buckets[i], pb[k]) >= 0)
                    {
                        k++;
                    }
                    else
                    {
                        pb.RemoveAt(k);
                        Length--;
                    }
                }
            }

            if (table._buckets[i] is { } tb)
                table.Length -= (ulong)tb.Count;
            table._buckets[i] = null;
        }
    }

    /// <summary><c>ajTableMergeEor</c>: the hits that both have are dropped from both, this gets the rest of the table (the table is emptied).</summary>
    public void MergeEor(EmbossHitTable table)
    {
        MakeSameSize(this, table);
        for (ulong i = 0; i < Size; i++)
        {
            var pb = _buckets[i];
            if (pb is not null)
            {
                for (int k = 0; k < pb.Count;)
                {
                    int m = IndexOf(table._buckets[i], pb[k]);
                    if (m < 0)
                    {
                        k++;
                        continue;
                    }

                    pb.RemoveAt(k);
                    Length--;
                    table._buckets[i]!.RemoveAt(m);
                    table.Length--;
                }
            }

            Absorb(table, i);
        }
    }

    /// <summary><c>ajTableMergeNot</c>: this loses the hits that the table has (the table is emptied).</summary>
    public void MergeNot(EmbossHitTable table)
    {
        MakeSameSize(this, table);
        for (ulong i = 0; i < Size; i++)
        {
            var pb = _buckets[i];
            if (pb is not null)
            {
                for (int k = 0; k < pb.Count;)
                {
                    if (IndexOf(table._buckets[i], pb[k]) >= 0)
                    {
                        pb.RemoveAt(k);
                        Length--;
                    }
                    else
                    {
                        k++;
                    }
                }
            }

            if (table._buckets[i] is { } tb)
                table.Length -= (ulong)tb.Count;
            table._buckets[i] = null;
        }
    }
}

/// <summary>
/// The database access method <c>emboss</c> (<c>textAccessEmboss</c> of <c>ajtextdb.c</c>): the index files of the <c>dbx</c> programs (<c>NAME.ent</c> with the names of the database files, <c>NAME.xid</c> / <c>xac</c> / <c>xsv</c> / <c>xde</c> / <c>xkw</c> / <c>xtx</c> with their
/// <c>.p</c> parameter files) are searched for the query and the entries are read from the database files at the offsets of the hits, one entry per hit, in the order of the hits. A5-8c3h, B09 F275.
/// </summary>
internal static class EmbossIndexAccess
{
    // ajBtreeReadEntriesC: the names of the database files ("directory/name"); the first line that is no comment is skipped by the library (it is the heading of the file)
    private static List<string>? ReadEntries(string file, string indexdir, string directory, Func<string, string?> files)
    {
        uint refcount = 0;
        string fn = indexdir.Length == 0 ? file : indexdir + "/" + file;
        fn += ".ent";
        string? text = files(fn);
        if (text is null)
            throw new EmbossFatalException($"\n   EMBOSS An error in ajindex.c at line 11171:\nCannot open database entries file {fn}");
        var lines = text.Split('\n').ToList();
        if (lines.Count > 0 && lines[^1].Length == 0)
            lines.RemoveAt(lines.Count - 1);
        int idx = 0;
        while (idx < lines.Count)
        {
            string line = lines[idx++].TrimEnd('\r');
            if (line.Length == 0 || line[0] == '#')
                continue;
            if (line.StartsWith("Dual", StringComparison.Ordinal))
                refcount = 1;
            if (line.StartsWith("Reference ", StringComparison.Ordinal))
            {
                var parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length > 1 && int.TryParse(parts[1], out int rc))
                    refcount = (uint)(rc - 1);
            }

            break;
        }

        var result = new List<string>();
        for (; idx < lines.Count; idx++)
        {
            string line = lines[idx].TrimEnd('\r');
            string name = refcount == 0 ? FirstWord(line) : line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries).FirstOrDefault() ?? "";
            result.Add(directory + "/" + name);
        }

        return result;
    }

    // ajFmtScanS(line, "%S", &name): the first word
    private static string FirstWord(string line)
    {
        int a = 0;
        while (a < line.Length && line[a] is ' ' or '\t' or '\n' or '\r' or '\v' or '\f')
            a++;
        int b = a;
        while (b < line.Length && line[b] is not (' ' or '\t' or '\n' or '\r' or '\v' or '\f'))
            b++;
        return line[a..b];
    }

    // ajFilenameTestInclude on the name of a file of the database (without its directory)
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

    /// <summary><c>textAccessEmboss</c> for the whole query: the hits in order, then the text of the file of every hit from its offset.</summary>
    public static EmbossQueryParse.ReadResult Emboss(EmbossQuery qry, EmbossQueryEnv env, EmbossFeatureReadContext context)
    {
        var fail = new EmbossQueryParse.ReadResult(false, null, "", false);
        if (qry.IndexDir.Length == 0)
        {
            context.Error($"no indexdir defined for database '{qry.DbName}'");
            return fail;
        }

        if (qry.Namespace.Length > 0 || qry.Organisms.Length > 0)
            throw new NotSupportedException("The namespace and organisms attributes of a database definition (textEmbossQryNamespace / textEmbossQryOrganisms) are not read by this port.");
        string alias = qry.DbAlias.Length > 0 ? qry.DbAlias : qry.DbName.ToLowerInvariant();
        if (qry.QueryType == EmbossQueryType.All)
        {
            // textEmbossAll: every file of the database that is included, read one after the other
            var texts = new List<string>();
            foreach (string f in ReadEntries(alias, qry.IndexDir, qry.Directory, env.Files)!)
            {
                if (!Included(f, qry.Exclude, qry.Filename))
                    continue;
                if (env.Files(f) is { } t)
                    texts.Add(t);
            }

            if (texts.Count == 0)
                return fail;
            qry.QryDone = !qry.CaseId;
            return new EmbossQueryParse.ReadResult(true, texts[0], qry.Filename, false, qry.DbName) { Texts = texts };
        }

        // textEmbossQryOpen: the id index and the index of every field of the query (the words of the query are lower case)
        var idcache = EmbossBtreeCache.OpenRead(alias, EmbossBtreeCache.FieldExtension("id"), qry.IndexDir, env.Files, context);
        if (idcache is null)
            return fail;
        var dbfiles = ReadEntries(alias, qry.IndexDir, qry.Directory, env.Files)!;
        var caches = new List<EmbossBtreeCache>();
        foreach (var field in qry.Fields)
        {
            var c = EmbossBtreeCache.OpenRead(alias, EmbossBtreeCache.FieldExtension(field.Field), qry.IndexDir, env.Files, context);
            if (c is null)
                return fail;
            caches.Add(c);
        }

        var hits = new List<EmbossBtHit>();
        var wild = new string[qry.Fields.Count];
        for (int i = 0; i < wild.Length; i++)
        {
            qry.Fields[i].Wildquery = qry.Fields[i].Wildquery.ToLowerInvariant();            // textEmbossQryOpen lowers the words of the query in place
            wild[i] = qry.Fields[i].Wildquery;
        }
        if (qry.QueryType == EmbossQueryType.Entry)
        {
            for (int i = 0; i < qry.Fields.Count; i++)
            {
                if (qry.Fields[i].Link == 'E' && hits.Count > 0)
                    continue;
                if (!caches[i].Secondary)
                    caches[i].IdentFetchHit(wild[i], hits, context);
            }
        }
        else
        {
            var table = new EmbossHitTable(100);
            ulong lasthits = 0;
            for (int i = 0; i < qry.Fields.Count; i++)
            {
                var fd = qry.Fields[i];
                if (fd.Link == 'E' && lasthits > 0)
                    continue;
                var list = new List<EmbossBtHit>();
                if (caches[i].Secondary)
                {
                    if (!qry.Wild)
                        caches[i].KeyFetchHit(idcache, wild[i], list);
                    else
                        caches[i].KeyFetchwildHit(idcache, wild[i], list);
                }
                else
                {
                    caches[i].IdentFetchwildHit(wild[i], list);
                }

                ulong fdhits = (ulong)list.Count;
                switch (fd.Link)
                {
                    case 'I':
                        foreach (var h in list)
                            table.Put(h);
                        break;
                    case '|':
                    case '&':
                    case '^':
                    case 'E':
                    case '!':
                        {
                            var nt = new EmbossHitTable(fdhits);
                            foreach (var h in list)
                                nt.Put(h);
                            if (fd.Link == '|')
                                table.MergeOr(nt);
                            else if (fd.Link == '&')
                                table.MergeAnd(nt);
                            else if (fd.Link == '!')
                                table.MergeNot(nt);
                            else
                                table.MergeEor(nt);
                            break;
                        }
                    default:
                        context.Error($"Unexpected query link operator number '{fd.Link}'");
                        break;
                }

                lasthits = fdhits;
            }

            hits.AddRange(table.Values());
        }

        if (hits.Count == 0)
            return fail;

        // textEmbossQryNext for every hit: the files that are excluded are skipped; a file that cannot be opened ends the read
        var skip = new bool[dbfiles.Count];
        if (qry.Exclude.Length > 0)
        {
            for (int i = 0; i < skip.Length; i++)
                skip[i] = !Included(dbfiles[i], qry.Exclude, qry.Filename);
        }

        var hitTexts = new List<string>();
        foreach (var h in hits)
        {
            if (h.Dbno < skip.Length && skip[h.Dbno])
                continue;
            string? t = h.Dbno < dbfiles.Count ? env.Files(dbfiles[(int)h.Dbno]) : null;
            if (t is null)
                break;
            hitTexts.Add(h.Offset < (ulong)t.Length ? t[(int)h.Offset..] : "");
        }

        if (hitTexts.Count == 0)
            return fail;
        qry.QryDone = !qry.CaseId;
        return new EmbossQueryParse.ReadResult(true, hitTexts[0], qry.Filename, false, qry.DbName) { HitTexts = hitTexts };
    }
}
