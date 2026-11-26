// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: CyclicSubsetSeed.cc, SubsetSuffixArray*.cc, MultiSequence.cc, lastdb.cc), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Last;

/// <summary>
/// A <c>MultiSequence</c> of LAST: the sequences, each followed by <c>padSize</c> delimiter symbols, in one array of letter codes (the array begins with the pad); the names; the end coordinates of the pads.
/// </summary>
internal sealed class LastMultiSequence
{
    public const int PadSize = 1;

    /// <summary>The code of the delimiter symbol (a blank): the number of letters of the DNA alphabet.</summary>
    public const byte PadCode = 4;

    public List<byte> Seq { get; } = [];

    /// <summary>The coordinate after the pad of every finished sequence (<c>ends</c>; the first entry is the pad at the beginning).</summary>
    public List<int> Ends { get; } = [PadSize];

    public List<string> Names { get; } = [];

    public List<bool> IsReverse { get; } = [];

    public int FinishedSequences => Ends.Count - 1;

    public bool IsFinished => Ends.Count == Names.Count + 1;

    public LastMultiSequence()
    {
        Seq.Add(PadCode);
    }

    public int SeqBeg(int i) => Ends[i];

    public int SeqEnd(int i) => Ends[i + 1] - PadSize;

    public int PadBeg(int i) => Ends[i] - PadSize;

    public int PadEnd(int i) => Ends[i + 1];

    public int SeqLen(int i) => SeqEnd(i) - SeqBeg(i);

    public int PadLen(int i) => PadEnd(i) - PadBeg(i);

    /// <summary><c>whichSequence</c>: <c>upper_bound</c> over the end coordinates.</summary>
    public int WhichSequence(int coordinate)
    {
        int lo = 0, hi = Ends.Count;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (Ends[mid] <= coordinate)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo - 1;
    }

    /// <summary>Starts a new sequence (the name is the first word of the header line; a sequence without a header word has an empty name).</summary>
    public void BeginSequence(string name)
    {
        Names.Add(name);
        IsReverse.Add(false);
    }

    public void FinishSequence()
    {
        for (int i = 0; i < PadSize; i++)
            Seq.Add(PadCode);
        Ends.Add(Seq.Count);
    }

    /// <summary>Removes the finished sequences (<c>reinitForAppending</c>; the port keeps an unfinished sequence outside the buffer).</summary>
    public void ReinitForAppending()
    {
        Seq.Clear();
        Seq.Add(PadCode);
        Names.Clear();
        IsReverse.Clear();
        Ends.Clear();
        Ends.Add(PadSize);
    }

    public string SeqName(int i) => Names[i];

    public char Strand(int i) => IsReverse[i] ? '-' : '+';

    /// <summary><c>reverseComplementOneSequence</c>.</summary>
    public void ReverseComplementOneSequence(int seqNum, byte[] complement)
    {
        int b = SeqBeg(seqNum), e = SeqEnd(seqNum);
        Seq.Reverse(b, e - b);
        for (int i = b; i < e; i++)
            Seq[i] = complement[Seq[i]];
        IsReverse[seqNum] = !IsReverse[seqNum];
    }
}

/// <summary>The seed of <c>lastdb</c> / <c>lastal</c> (<c>CyclicSubsetSeed</c>): per position of the pattern a table letter code -> subset number (255 = not indexed).</summary>
internal sealed class LastSeed
{
    public const int Delimiter = 255;

    private readonly List<byte[]> _maps = [];

    public int Span => _maps.Count;

    public byte[] Map(int depth) => _maps[depth % _maps.Count];

    /// <summary>The default seed of <c>lastdb</c> for DNA (<c>YASS</c>: "1T1001100101", 1 = one subset per letter, T = {A, G} {C, T}, 0 = all letters one subset).</summary>
    public static LastSeed Yass(LastAlphabet alph)
    {
        var seed = new LastSeed();
        foreach (char p in "1T1001100101")
        {
            string[] subsets = p switch { '1' => ["A", "C", "G", "T"], 'T' => ["AG", "CT"], _ => ["ACGT"] };
            var map = new byte[LastAlphabet.Capacity];
            Array.Fill(map, (byte)Delimiter);
            for (int s = 0; s < subsets.Length; s++)
            {
                foreach (char c in subsets[s])
                {
                    map[alph.Encode[c]] = (byte)s;
                    char lower = char.ToLowerInvariant(c);
                    if (lower != c)
                        map[alph.Encode[lower]] = (byte)s;
                }
            }

            seed._maps.Add(map);
        }

        return seed;
    }
}

/// <summary>
/// The index of a database: the start positions of all suffixes whose first letter is in the seed alphabet, sorted by the subset numbers of the letters along the cyclic seed (equal strings in order of position),
/// and the adaptive search of <c>SubsetSuffixArray::match</c>.
/// </summary>
internal sealed class LastSuffixIndex
{
    private readonly byte[] _text;
    private readonly LastSeed _seed;

    public int[] Positions { get; }

    public LastSuffixIndex(byte[] text, int textLength, LastSeed seed)
    {
        _text = text;
        _seed = seed;
        var list = new List<int>();
        byte[] first = seed.Map(0);
        for (int i = 0; i < textLength; i++)
            if (first[text[i]] < LastSeed.Delimiter)
                list.Add(i);
        var keys = list.ToArray();
        // Array.Sort is not stable, so the position breaks ties
        Array.Sort(keys, (a, b) =>
        {
            int c = CompareSuffixes(a, b);
            return c != 0 ? c : a.CompareTo(b);
        });
        Positions = keys;
    }

    private int CompareSuffixes(int a, int b)
    {
        int depth = 0;
        while (true)
        {
            byte[] map = _seed.Map(depth);
            byte x = map[_text[a + depth]];
            byte y = map[_text[b + depth]];
            if (x != y)
                return x.CompareTo(y);
            if (x == LastSeed.Delimiter)
                return 0;
            depth++;
        }
    }

    /// <summary>The subset number of the letter at <paramref name="depth"/> of the suffix at <paramref name="position"/>.</summary>
    private byte SubsetAt(int position, int depth)
    {
        int p = position + depth;
        return p < _text.Length ? _seed.Map(depth)[_text[p]] : (byte)LastSeed.Delimiter;
    }

    /// <summary>
    /// <c>SubsetSuffixArray::match</c>: the smallest depth d >= minDepth for which at most <paramref name="maxHits"/> suffixes match the first d letters of the query (as subsets); the suffixes in index order.
    /// Empty if the query reaches a letter outside the seed alphabet first.
    /// </summary>
    public (int Begin, int End) Match(byte[] query, int queryPos, int maxHits, int minDepth, int maxDepth)
    {
        int beg = 0, end = Positions.Length;
        int depth = 0;
        while (true)
        {
            if (depth >= minDepth && end - beg <= maxHits)
                return (beg, end);
            if (depth >= maxDepth)
                break;
            byte subset = _seed.Map(depth)[query[queryPos + depth]];
            if (subset == LastSeed.Delimiter)
                return (0, 0);
            // narrow [beg, end) to the suffixes whose letter at this depth has this subset
            int lo = beg, hi = end;
            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                if (SubsetAt(Positions[mid], depth) < subset)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            int newBeg = lo;
            hi = end;
            while (lo < hi)
            {
                int mid = lo + (hi - lo) / 2;
                if (SubsetAt(Positions[mid], depth) <= subset)
                    lo = mid + 1;
                else
                    hi = mid;
            }

            beg = newBeg;
            end = lo;
            depth++;
            if (beg == end)
                return (0, 0);
        }

        return (beg, end);
    }
}
