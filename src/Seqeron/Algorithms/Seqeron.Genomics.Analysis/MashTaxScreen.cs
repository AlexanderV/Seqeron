using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The text of a <c>names.dmp</c> / <c>nodes.dmp</c> pair (the NCBI taxonomy dump) and the optional reference-to-taxonomy mapping file of <c>mash taxscreen -t DIR -m FILE</c>.
/// </summary>
/// <param name="NamesDump">The content of <c>names.dmp</c> (<c>tax_id TAB | TAB name TAB | TAB unique name TAB | TAB name class TAB |</c> per line).</param>
/// <param name="NodesDump">The content of <c>nodes.dmp</c> (<c>tax_id TAB | TAB parent TAB | TAB rank TAB | ...</c> per line).</param>
/// <param name="MappingFile">The content of the mapping file (<c>taxid SEP reference-name</c> per line) or null (<c>-m</c> not given): the taxonomy ID of a reference is that of the file, else the number after the word <c>taxid</c> of its comment.</param>
public sealed record MashTaxonomy(string NamesDump, string NodesDump, string? MappingFile = null);

/// <summary>
/// <c>mash taxscreen</c> (<c>CommandTaxScreen.cpp</c> + <c>taxdb.hpp</c>, marbl/Mash v2.3; B09 A5-10b2, F215): the pass of <c>mash screen</c> over a pool of sequences, the taxonomy of every query reference, the lowest common ancestor of the references that share a hash,
/// and the Kraken-style report of the clade and taxon counts. Reproduced with the program's own containers where they decide the output: the iteration order of its <c>robin_hood::unordered_map</c> (<see cref="RobinHoodMap"/>) fixes the order in which the
/// taxa enter <c>unordered_map&lt;TaxID, TaxCounts&gt;</c>, whose libstdc++ node list (<see cref="LibstdcxxHashMap{T}"/>) is walked while the loop inserts the ancestors (a rehash in the middle changes what the walk visits: reproduced), and the children
/// of a taxon are ordered by <c>std::sort</c> (<c>Seqeron.Genomics.Alignment.LibstdcxxSort</c>, not stable above 16 elements).
/// </summary>
public static class MashTaxScreen
{
    private sealed class Entry(ulong taxId, string rank)
    {
        public readonly ulong TaxId = taxId;
        public readonly string Rank = rank;
        public Entry? Parent;
        public string Name = "";
    }

    private sealed class Counts
    {
        public ulong CladeCount, TaxCount, TaxHashCount, CladeHashCount;
        public List<ulong> Children = [];
    }

    /// <summary>
    /// Runs <c>mash taxscreen [-m FILE] [-t DIR] queries.msh pool...</c>.
    /// </summary>
    /// <param name="queries">The query sketch (what <c>mash sketch</c> wrote; one sketch file, or sequence files sketched with the parameters below).</param>
    /// <param name="pools">The pool files (reads or contigs, nucleotides; translated in 6 frames for an amino-acid query sketch).</param>
    /// <param name="taxonomy">The taxonomy dump and the mapping file.</param>
    /// <param name="k">K-mer size (when the queries are sequences).</param>
    /// <param name="sketchSize">Sketch size (when the queries are sequences).</param>
    /// <param name="seed">Hash seed (when the queries are sequences).</param>
    /// <param name="canonical">False for a <c>-n</c> sketch.</param>
    /// <param name="individual">The sketch was made with <c>-i</c>.</param>
    /// <param name="alphabet">The alphabet of the queries.</param>
    /// <param name="cancellationToken">Polled per chunk and per taxon.</param>
    /// <returns>The report.</returns>
    /// <exception cref="ArgumentException">A dump line cannot be read as the program reads it, a node's parent is missing from the dump (the program leaves the pointer uninitialised: no defined result), or the pool has no record.</exception>
    public static string Run(IReadOnlyList<MashInput> queries, IReadOnlyList<MashInput> pools, MashTaxonomy taxonomy, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, MashAlphabet? alphabet = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentNullException.ThrowIfNull(pools);
        ArgumentNullException.ThrowIfNull(taxonomy);
        ArgumentNullException.ThrowIfNull(taxonomy.NamesDump);
        ArgumentNullException.ThrowIfNull(taxonomy.NodesDump);
        var entries = ReadTaxonomy(taxonomy.NamesDump, taxonomy.NodesDump);
        var scan = MashCommands.ScreenPass(queries, pools, k, sketchSize, seed, canonical, individual, alphabet, cancellationToken);

        // the taxonomy ID of every reference: the mapping file first, then the number after the word "taxid" of its comment
        var referenceTaxIds = new ulong[scan.References.Count];
        if (taxonomy.MappingFile is not null)
        {
            var map = new Dictionary<string, ulong>(StringComparer.Ordinal);
            var s = new CppStream(taxonomy.MappingFile);
            while (s.ReadUInt64(out ulong id))
            {
                s.Ignore(1);
                string name = s.GetLine('\n');
                map.TryAdd(name, id);
            }

            for (int i = 0; i < referenceTaxIds.Length; i++)
            {
                if (map.TryGetValue(scan.References[i].Name, out ulong id))
                    referenceTaxIds[i] = id;
            }
        }

        for (int i = 0; i < referenceTaxIds.Length; i++)
        {
            if (referenceTaxIds[i] != 0)
                continue;
            var words = new CppStream(scan.References[i].Comment);
            ulong taxId = referenceTaxIds[i];
            while (words.ReadWord(out string word))
            {
                if (word == "taxid")
                    words.ReadUInt64(out taxId);
            }

            referenceTaxIds[i] = taxId;
        }

        // for each hash (in the iteration order of the program's robin_hood map) the lowest common ancestor of the taxa of the references that hold it
        var counts = new LibstdcxxHashMap<Counts>(() => new Counts());
        foreach (ulong hash in scan.HashesInProgramOrder)
        {
            cancellationToken.ThrowIfCancellationRequested();
            ulong taxId = 0;
            foreach (int owner in scan.OwnersInProgramOrder(hash))
                taxId = LowestCommonAncestor(entries, referenceTaxIds[owner], taxId);
            var c = counts.GetOrAdd(taxId);
            c.TaxHashCount += 1;
            if (scan.Counts[hash] >= 1)
                c.TaxCount += 1;
        }

        // the clade counts: the loop walks the node list of `counts` and inserts the ancestors on the way
        ulong totalCount = 0;
        for (var it = counts.Begin; it is not null; it = it.Next)
        {
            cancellationToken.ThrowIfCancellationRequested();
            ulong hashCount = it.Value.TaxHashCount;
            ulong count = it.Value.TaxCount;
            totalCount += count;
            entries.TryGetValue(it.Key, out var taxon);
            while (taxon is not null)
            {
                var own = counts.GetOrAdd(taxon.TaxId);
                own.CladeCount += count;
                own.CladeHashCount += hashCount;
                if (taxon.Parent is null)
                    break;
                var children = counts.GetOrAdd(taxon.Parent.TaxId).Children;
                int at = children.BinarySearch(taxon.TaxId);                    // lower_bound over the sorted ids
                if (at < 0)
                    children.Insert(~at, taxon.TaxId);
                taxon = taxon.Parent;
            }
        }

        var sb = new StringBuilder();
        sb.Append("%\thashes\ttaxHashes\thashesDB\ttaxHashesDB\ttaxID\trank\tname\n");
        var zero = counts.Find(0);
        if (zero is not null && (uint)zero.CladeCount > 0)                      // "Should not happen"
            sb.Append(CFormat.Fixed(100 * (uint)zero.CladeCount / (double)totalCount, 4)).Append('\t').Append(CInt(zero.CladeCount)).Append('\t').Append(CInt(zero.TaxCount)).Append("\tno rank\t0\tunclassified\n");
        WriteReport(sb, entries, counts, totalCount, 1, 0, cancellationToken);
        return sb.ToString();
    }

    private static string CInt(ulong v) => ((int)(uint)v).ToString(CultureInfo.InvariantCulture);          // printf("%i") of an unsigned int

    private static void WriteReport(StringBuilder sb, Dictionary<ulong, Entry> entries, LibstdcxxHashMap<Counts> counts, ulong totalCount, ulong taxId, int depth, CancellationToken ct)
    {
        ct.ThrowIfCancellationRequested();
        var it = counts.Find(taxId);
        uint cladeCount = it is null ? 0 : (uint)it.CladeCount;
        if (cladeCount == 0)
            return;
        uint cladeHashCount = (uint)it!.CladeHashCount, taxCount = (uint)it.TaxCount, taxHashCount = (uint)it.TaxHashCount;
        if (!entries.TryGetValue(taxId, out var taxon))
            throw new ArgumentException($"The taxonomy has no entry for the taxon {taxId} that holds hashes: the program dereferences a null pointer here.");
        sb.Append(CFormat.Fixed(100 * cladeCount / (double)totalCount, 4)).Append('\t').Append(CInt(cladeCount)).Append('\t').Append(CInt(taxCount)).Append('\t').Append(CInt(cladeHashCount)).Append('\t').Append(CInt(taxHashCount)).Append('\t')
            .Append(taxon.Rank).Append('\t').Append(taxId.ToString(CultureInfo.InvariantCulture)).Append('\t').Append(' ', 2 * depth).Append(taxon.Name).Append('\n');
        var children = it.Children.ToArray();
        // std::sort(children, [&](int a, int b) { return counts.at(a).cladeCount > counts.at(b).cladeCount; })
        global::Seqeron.Genomics.Alignment.LibstdcxxSort.Sort(children, (a, b) => counts.Find((ulong)(long)(int)a)!.CladeCount > counts.Find((ulong)(long)(int)b)!.CladeCount);
        foreach (ulong child in children)
        {
            if (counts.Find(child) is null)
                break;
            WriteReport(sb, entries, counts, totalCount, child, depth + 1, ct);
        }
    }

    // TaxDB::getLowestCommonAncestor
    private static ulong LowestCommonAncestor(Dictionary<ulong, Entry> entries, ulong a, ulong b)
    {
        if (b == 0)
            return a;
        if (a == 0)
            return b;
        if (!entries.TryGetValue(a, out var ta) || !entries.TryGetValue(b, out var tb))
            return 1;
        var path = new HashSet<Entry>();
        var pta = ta;
        while (pta is not null && pta.TaxId > 1 && pta.Parent is not null)
        {
            if (pta.TaxId == b)
                return b;
            path.Add(pta);
            pta = pta.Parent;
        }

        var ptb = tb;
        while (ptb.TaxId > 0 && ptb.Parent is not null)
        {
            if (path.Contains(ptb))
                return ptb.TaxId;
            ptb = ptb.Parent;
        }

        return 1;
    }

    // TaxDB::TaxDB(namesDump, nodesDump): parseNodesDump, the parent links, parseNamesDump (the stream operations of the C++ code)
    private static Dictionary<ulong, Entry> ReadTaxonomy(string namesDump, string nodesDump)
    {
        var entries = new Dictionary<ulong, Entry>();
        var parents = new List<(Entry Entry, ulong Parent)>();
        var nodes = new CppStream(nodesDump);
        while (nodes.ReadUInt64(out ulong taxId) && nodes.ReadChar(out _) && nodes.ReadUInt64(out ulong parentTaxId) && nodes.ReadChar(out _))
        {
            nodes.Ignore(1);
            string rank = nodes.GetLine('\t');
            if (entries.TryAdd(taxId, new Entry(taxId, rank)))
                parents.Add((entries[taxId], parentTaxId));
            nodes.Ignore(2560, '\n');
        }

        foreach (var (entry, parentId) in parents)
        {
            if (entry.TaxId == parentId)
                continue;
            if (!entries.TryGetValue(parentId, out var parent))
                throw new ArgumentException($"The taxon {entry.TaxId} has the parent {parentId}, which nodes.dmp does not list: the program leaves the parent pointer of the entry uninitialised (undefined behaviour).");
            entry.Parent = parent;
        }

        var names = new CppStream(namesDump);
        while (names.ReadUInt64(out ulong taxId))
        {
            names.Ignore(3);
            string name = names.GetLine('\t');
            names.Ignore(3);
            names.Ignore(256, '|');
            names.Ignore(1);
            string type = names.GetLine('\t');
            if (type == "scientific name" && entries.TryGetValue(taxId, out var entry))
                entry.Name = name;
            names.Ignore(2560, '\n');
        }

        return entries;
    }

    // the std::istream operations the dump readers use (state bits as libstdc++: a failed or exhausted stream extracts nothing more)
    private sealed class CppStream(string text)
    {
        private int _pos;
        private bool _eof, _fail;

        private bool Good => !_eof && !_fail;

        private void SkipSpace()
        {
            while (_pos < text.Length && char.IsWhiteSpace(text[_pos]))
                _pos++;
            if (_pos >= text.Length)
                _eof = true;
        }

        // operator>>(unsigned long long &): optional sign, decimal digits; no digit: failure (the value becomes 0); overflow: the maximum and failure
        public bool ReadUInt64(out ulong value)
        {
            value = 0;
            if (!Good)
            {
                _fail = true;
                return false;
            }

            SkipSpace();
            if (_eof)
            {
                _fail = true;
                return false;
            }

            bool negative = false;
            if (text[_pos] is '+' or '-')
            {
                negative = text[_pos] == '-';
                _pos++;
            }

            int digits = 0;
            bool overflow = false;
            ulong v = 0;
            while (_pos < text.Length && text[_pos] is >= '0' and <= '9')
            {
                ulong d = (ulong)(text[_pos] - '0');
                if (v > (ulong.MaxValue - d) / 10)
                    overflow = true;
                else
                    v = v * 10 + d;
                digits++;
                _pos++;
            }

            if (_pos >= text.Length)
                _eof = true;
            if (digits == 0)
            {
                _fail = true;
                return false;
            }

            if (overflow)
            {
                value = ulong.MaxValue;
                _fail = true;
                return false;
            }

            value = negative ? 0UL - v : v;
            return !_fail;
        }

        // operator>>(std::string &): the next run of non-blank characters
        public bool ReadWord(out string word)
        {
            word = "";
            if (!Good)
            {
                _fail = true;
                return false;
            }

            SkipSpace();
            if (_eof)
            {
                _fail = true;
                return false;
            }

            int start = _pos;
            while (_pos < text.Length && !char.IsWhiteSpace(text[_pos]))
                _pos++;
            if (_pos >= text.Length)
                _eof = true;
            word = text[start.._pos];
            return true;
        }

        // operator>>(char &)
        public bool ReadChar(out char c)
        {
            c = '\0';
            if (!Good)
            {
                _fail = true;
                return false;
            }

            SkipSpace();
            if (_eof)
            {
                _fail = true;
                return false;
            }

            c = text[_pos++];
            if (_pos >= text.Length)
                _eof = true;
            return true;
        }

        public void Ignore(int count) => Ignore(count, null);

        public void Ignore(int count, char? delimiter)
        {
            if (!Good)
            {
                _fail = true;
                return;
            }

            for (int n = 0; n < count; n++)
            {
                if (_pos >= text.Length)
                {
                    _eof = true;
                    return;
                }

                char c = text[_pos++];
                if (delimiter is { } d && c == d)
                    return;
            }
        }

        public string GetLine(char delimiter)
        {
            if (!Good)
            {
                _fail = true;
                return "";
            }

            int start = _pos;
            while (_pos < text.Length && text[_pos] != delimiter)
                _pos++;
            string line = text[start.._pos];
            if (_pos >= text.Length)
            {
                _eof = true;
                if (line.Length == 0)
                    _fail = true;
            }
            else
            {
                _pos++;
            }

            return line;
        }
    }
}

/// <summary>
/// The node list of a libstdc++ <c>std::unordered_map&lt;uint64_t, T&gt;</c> (GCC 13: identity hash, prime bucket counts growing 1, 13, 29, 59, 127, ..., a singly linked list in which a new node goes in front of the first node of its bucket, or at the
/// front of the list when the bucket is empty, and a rehash re-links the nodes in the order of the old list): the iteration order and what an insertion does to a walk that is in progress. B09 F215.
/// </summary>
internal sealed class LibstdcxxHashMap<T>(Func<T> create)
{
    internal sealed class Node(ulong key, T value)
    {
        public readonly ulong Key = key;
        public readonly T Value = value;
        public Node? Next;
    }

    private static long[] BucketCounts => LibstdcxxBuckets.Counts;

    private readonly Node _beforeBegin = new(0, default!);
    private Node?[] _buckets = new Node?[1];
    private int _sizeIndex;
    private long _count;
    private long _nextResize;

    /// <summary>The first node of the list.</summary>
    public Node? Begin => _beforeBegin.Next;

    private static int Bucket(ulong key, int count) => (int)(key % (ulong)count);

    /// <summary>The node of the key or null.</summary>
    public T? Find(ulong key)
    {
        var node = FindNode(key);
        return node is null ? default : node.Value;
    }

    private Node? FindNode(ulong key)
    {
        int bkt = Bucket(key, _buckets.Length);
        var prev = _buckets[bkt];
        if (prev is null)
            return null;
        for (var p = prev.Next; p is not null; p = p.Next)
        {
            if (p.Key == key)
                return p;
            if (p.Next is null || Bucket(p.Next.Key, _buckets.Length) != bkt)
                break;
        }

        return null;
    }

    /// <summary><c>operator[]</c>: the value of the key, inserted (with the default value) when absent.</summary>
    public T GetOrAdd(ulong key)
    {
        var found = FindNode(key);
        if (found is not null)
            return found.Value;
        var node = new Node(key, create());
        // _M_need_rehash(bucket count, element count, 1): the table grows to the next prime of the sequence when the new element exceeds the next resize
        if (_count + 1 > _nextResize)
        {
            _sizeIndex++;
            if (_sizeIndex >= BucketCounts.Length)
                throw new NotSupportedException("The hash table of the taxonomy report grew beyond the sizes this emulation covers.");
            Rehash((int)BucketCounts[_sizeIndex]);
        }

        int bkt = Bucket(key, _buckets.Length);
        if (_buckets[bkt] is { } before)
        {
            node.Next = before.Next;
            before.Next = node;
        }
        else
        {
            node.Next = _beforeBegin.Next;
            _beforeBegin.Next = node;
            if (node.Next is not null)
                _buckets[Bucket(node.Next.Key, _buckets.Length)] = node;
            _buckets[bkt] = _beforeBegin;
        }

        _count++;
        return node.Value;
    }

    private void Rehash(int newCount)
    {
        var buckets = new Node?[newCount];
        var p = _beforeBegin.Next;
        _beforeBegin.Next = null;
        int bbeginBkt = 0;
        while (p is not null)
        {
            var next = p.Next;
            int bkt = Bucket(p.Key, newCount);
            if (buckets[bkt] is null)
            {
                p.Next = _beforeBegin.Next;
                _beforeBegin.Next = p;
                buckets[bkt] = _beforeBegin;
                if (p.Next is not null)
                    buckets[bbeginBkt] = p;
                bbeginBkt = bkt;
            }
            else
            {
                p.Next = buckets[bkt]!.Next;
                buckets[bkt]!.Next = p;
            }

            p = next;
        }

        _buckets = buckets;
        _nextResize = newCount;
    }
}

/// <summary>The bucket counts of a libstdc++ <c>unordered_map</c> as it grows (GCC 13: the first insertion allocates 13, then the next prime of the library's list at least twice the previous count).</summary>
internal static class LibstdcxxBuckets
{
    internal static readonly long[] Counts = [1, 13, 29, 59, 127, 257, 541, 1109, 2357, 5087, 10273, 20753, 42043, 85229, 172933, 351061];
}
