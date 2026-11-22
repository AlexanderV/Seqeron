namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The word table of <c>embword.c</c> (<c>embWordGetTable</c>, <c>embWordPrintTableFI</c>; B09 WP62, F281): an <c>ajTable</c> (<c>ajtable.c</c>) with the hash function <c>wordStrHash</c> (<c>toupper(c) + 31 * hash</c> over the word length) and the comparison
/// <c>ajCharCmpCaseLen</c>: a prime number of buckets from the length of the first sequence, a new key at the head of its bucket, the rehash of <c>ajTableResizeHashsize</c> when the table has more than 8 keys per bucket. The words are printed by the
/// stable merge sort of <c>qsort</c> by decreasing count over the bucket order of the table, so the order of equal counts is the order of the table.
/// </summary>
internal sealed class EmbossWordTable
{
    private static readonly long[] Primes =
    [
        3, 7, 13, 31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573, 2097143, 4194301, 8388593, 16777213, 33554393, 67108859, 134217689, 268435399, 536870909, 1073741789, 2147483647, uint.MaxValue,
    ];

    private sealed class Node(string key, int count)
    {
        public string Key = key;
        public int Count = count;
        public Node? Link;
    }

    private Node?[] _buckets;
    private long _size;
    private long _length;
    private readonly int _wordLength;

    public EmbossWordTable(long sequenceLength, int wordLength)
    {
        _wordLength = wordLength;
        long hint = sequenceLength / 4;
        int i = 1;
        while (Primes[i] <= hint)
            i++;
        _size = Primes[i - 1];
        _buckets = new Node?[_size];
    }

    public long Length => _length;

    private long Hash(string key, long hashsize)
    {
        ulong hashval = 0;
        for (int i = 0; i < _wordLength; i++)
            hashval = char.ToUpperInvariant(key[i]) + 31 * hashval;
        return (long)(hashval % (ulong)hashsize);
    }

    private bool Equal(string a, string b)
    {
        for (int i = 0; i < _wordLength; i++)
        {
            if (char.ToUpperInvariant(a[i]) != char.ToUpperInvariant(b[i]))
                return false;
        }

        return true;
    }

    private Node? Find(string key)
    {
        for (var node = _buckets[Hash(key, _size)]; node is not null; node = node.Link)
        {
            if (Equal(key, node.Key))
                return node;
        }

        return null;
    }

    /// <summary>The count of the word (it was seen); a new word is added with the count 1.</summary>
    public void Add(string text, int start)
    {
        string word = text.Substring(start, _wordLength);
        var rec = Find(word);
        if (rec is not null)
        {
            rec.Count++;
            return;
        }

        // ajTablePut: the table is resized first when it holds more than 8 keys per bucket
        long minsize = _length / 8;
        if (_size < minsize)
            Resize(6 * minsize / 4);
        long i = Hash(word, _size);
        var node = new Node(word, 1) { Link = _buckets[i] };
        _buckets[i] = node;
        _length++;
    }

    // ajTableResizeHashsize
    private void Resize(long hashsize)
    {
        int p = 1;
        while (Primes[p] <= hashsize)
            p++;
        long iprime = Primes[p - 1];
        if (iprime == _size)
            return;
        var nodes = new List<Node>();
        for (long b = 0; b < _size; b++)
        {
            for (var node = _buckets[b]; node is not null; node = node.Link)
                nodes.Add(node);
        }

        _buckets = new Node?[iprime];
        _size = iprime;
        foreach (var node in nodes)
        {
            long i = Hash(node.Key, _size);
            node.Link = _buckets[i];
            _buckets[i] = node;
        }
    }

    /// <summary><c>embWordPrintTableFI</c>: the words with at least <paramref name="mincount"/> occurrences by decreasing count (words of equal count in the order of the table), <c>word TAB count</c>.</summary>
    public string Print(int mincount)
    {
        var values = new List<Node>();
        for (long b = 0; b < _size; b++)
        {
            for (var node = _buckets[b]; node is not null; node = node.Link)
                values.Add(node);
        }

        var sorted = values.OrderByDescending(n => n.Count).ToList();            // stable, as the merge sort of qsort is
        var sb = new System.Text.StringBuilder();
        foreach (var n in sorted)
        {
            if (n.Count < mincount)
                break;
            sb.Append(n.Key).Append('\t').Append(n.Count.ToString(System.Globalization.CultureInfo.InvariantCulture)).Append('\n');
        }

        return sb.ToString();
    }

    /// <summary><c>embWordGetTable</c>: adds the words of the sequence to <paramref name="table"/> (created from the first sequence that is long enough); false when the sequence is too short or has no word without the ambiguity code.</summary>
    public static bool GetTable(ref EmbossWordTable? table, string seq, bool nuc, int wordLength)
    {
        char skip = nuc ? 'N' : 'X';
        int len = seq.Length;
        if (len < wordLength)
            return false;
        table ??= new EmbossWordTable(len, wordLength);
        char At(int p) => p < len ? char.ToUpperInvariant(seq[p]) : '\0';
        int start = 0;
        int i = 0;
        int ilast = len - wordLength;
        int j = 0;
        while (j < wordLength)
        {
            if (At(start + j) == skip)
            {
                while (At(start + j) == skip)
                {
                    i++;
                    start++;
                }

                j = 0;
                if (i > ilast)
                    return false;
            }
            else
            {
                j++;
            }
        }

        j = wordLength - 1;
        while (i <= ilast)
        {
            if (At(start + j) == skip)
            {
                while (At(start + j) == skip)
                {
                    i++;
                    start++;
                }

                i += j;
                start += j;
                continue;
            }

            table.Add(seq, start);
            start++;
            i++;
        }

        return true;
    }
}
