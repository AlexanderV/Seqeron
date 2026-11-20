using Seqeron.Genomics.Infrastructure;
namespace Seqeron.Genomics.Analysis;

/// <summary>
/// One exact word match of <see cref="GenomicAnalyzer.SeqMatchAll"/>: the maximal exact match (at least the word size) of the sequences with indices
/// <see cref="Sequence1"/> &lt; <see cref="Sequence2"/>.
/// </summary>
/// <param name="Sequence1">0-based index of the first sequence of the pair (always the smaller index).</param>
/// <param name="Sequence2">0-based index of the second sequence.</param>
/// <param name="Position1">0-based start in the first sequence (EMBOSS <c>seqmatchall -aformat match</c> shows Position1 + 1).</param>
/// <param name="Position2">0-based start in the second sequence.</param>
/// <param name="Length">Length of the match (the EMBOSS score).</param>
public readonly record struct WordMatch(int Sequence1, int Sequence2, int Position1, int Position2, int Length);

public static partial class GenomicAnalyzer
{
    #region seqmatchall (all-against-all word matches)

    /// <summary>
    /// EMBOSS <c>seqmatchall</c> (6.6.0): the all-against-all comparison of a set of sequences by exact word matches. For every pair i &lt; j with
    /// length(i) &gt; <paramref name="wordSize"/> and length(j) ≥ <paramref name="wordSize"/> every <b>maximal exact match</b> of at least
    /// <paramref name="wordSize"/> symbols is reported (forward strand only; matches may overlap and a position may occur in several matches).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The program builds a table of the words of one sequence and merges runs of consecutive words on a diagonal (<c>embWordBuildMatchTable</c>); the result is exactly the set of
    /// maximal exact matches (left- and right-maximal common substrings with a longer-than-k extension), checked against the binary on 390 random pairs and on the oracle of the
    /// unit tests. It is computed here by the canonical suffix-tree engine of <see cref="FindCommonRegions(string, string, int, CommonRegionMode, StrandSelection, CancellationToken)"/>
    /// (<see cref="CommonRegionMode.MaximalExactMatches"/>): the first sequence is indexed once and every later sequence is a query. Symbols are compared
    /// case-insensitively as literal characters (a protein set works). The program's word table leaves out words around the unknown symbol of the indexed sequence's type (N for a nucleotide sequence, i.e. every character a nucleotide code; X otherwise): its scan only inspects the symbol that enters a window, so after a skipped N a window may still hold an N and match; that scan is reproduced step by step (<c>WordTableMatches</c>), the suffix-tree engine is used when the indexed sequence has no such symbol (R, Y, ... and a protein's N match literally).
    /// Matches are listed per pair in the order of the pairs (0,1), (0,2), …, (1,2), …, then by position in the first sequence and in the second (the program lists them in the
    /// order of its hash table: compare as sets).
    /// </para>
    /// </remarks>
    /// <param name="sequences">The sequences (nucleotide or protein).</param>
    /// <param name="wordSize">Word size (<c>-wordsize</c>, default 4; at least 1).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or an element is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="wordSize"/> &lt; 1.</exception>
    public static IReadOnlyList<WordMatch> SeqMatchAll(IReadOnlyList<string> sequences, int wordSize = 4, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (wordSize < 1)
            throw new ArgumentOutOfRangeException(nameof(wordSize), wordSize, "wordSize must be at least 1 (seqmatchall -wordsize).");
        foreach (var s in sequences)
            ArgumentNullException.ThrowIfNull(s);
        sequences = [.. sequences.Select(EmbossSequenceReader.Clean)];       // the EMBOSS reader drops the gap characters and every other non-letter but * (A4-3)

        var result = new List<WordMatch>();
        for (int i = 0; i < sequences.Count; i++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (sequences[i].Length <= wordSize)
                continue; // seqmatchall: length(seq1) must exceed the word size
            var pair = new PairWordMatcher(sequences[i], wordSize);
            for (int j = i + 1; j < sequences.Count; j++)
            {
                if (sequences[j].Length < wordSize)
                    continue;
                foreach (var m in pair.Matches(sequences[j]).OrderBy(m => m.PositionInText).ThenBy(m => m.PositionInQuery))
                    result.Add(new WordMatch(i, j, m.PositionInText, m.PositionInQuery, m.Length));
            }
        }

        return result;
    }

    /// <summary>
    /// The matches of the words of one sequence (the indexed one, <c>embWordGetTable</c>) against the other sequences (<c>embWordBuildMatchTable</c>): shared by <see cref="SeqMatchAll"/> and
    /// <see cref="MegaMerge"/>.
    /// </summary>
    private sealed class PairWordMatcher
    {
        private readonly string _indexed;
        private readonly string _upper;
        private readonly int _k;
        private readonly char _skip;
        private readonly bool _plain;
        private global::SuffixTree.SuffixTree? _tree;
        private bool[]? _allowed;

        public PairWordMatcher(string indexed, int wordSize)
        {
            _indexed = indexed;
            _k = wordSize;
            // embWordGetTable: the unknown symbol of the first sequence's type (N for a nucleotide sequence, X otherwise) controls which words enter the table
            _skip = EmbossSequenceReader.IsNucleotide(indexed) ? 'N' : 'X';
            _upper = indexed.ToUpperInvariant();
            // The common case (no skip symbol in the indexed sequence, letters only): every word is in the table and the matches are the maximal exact matches of the canonical suffix-tree
            // engine. Otherwise the program's own table scan decides (WordTableMatches).
            _plain = _upper.All(c => c is >= 'A' and <= 'Z' && c != _skip);
        }

        /// <summary>False when the program's word table of the indexed sequence is empty (sequence shorter than the word, or every word starts inside a run of the skip symbol).</summary>
        public bool HasTable => _upper.Length >= _k && (_plain || WordTableStarts(_upper, _k, _skip).Any(a => a));

        public IEnumerable<(int PositionInText, int PositionInQuery, int Length)> Matches(string query)
        {
            if (_plain)
            {
                _tree ??= BuildMaskedTree(_indexed, BuildLetterMask("ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
                return _tree.FindMaximalExactMatches(LetterQuery(query, skip: '\0'), _k);
            }

            _allowed ??= WordTableStarts(_upper, _k, _skip);
            return WordTableMatches(_upper, query.ToUpperInvariant(), _k, _allowed);
        }
    }

    // embWordGetTable's scan: which start positions of the indexed sequence enter the word table. The program only inspects the symbol that enters a window on the right (and the
    // whole first window), and after a skipped symbol restarts the window right behind the run of skipped symbols - so a later window may still contain a skip symbol in its interior.
    private static bool[] WordTableStarts(string seq, int k, char skip)
    {
        int len = seq.Length, ilast = len - k;
        var allowed = new bool[len];
        char At(int p) => p < len ? seq[p] : '\0';
        int i = 0, j = 0;
        while (j < k)
        {
            if (At(i + j) == skip)
            {
                while (At(i + j) == skip)
                    i++;
                j = 0;
                if (i > ilast)
                    return allowed; // no word without the skip symbol: no table
            }
            else
            {
                j++;
            }
        }

        j = k - 1;
        while (i <= ilast)
        {
            if (At(i + j) == skip)
            {
                while (At(i + j) == skip)
                    i++;
                i += j;
                continue;
            }

            allowed[i] = true;
            i++;
        }

        return allowed;
    }

    // embWordBuildMatchTable: every window of the query is looked up in the table of the indexed sequence (literal comparison, no skipping on the query side); consecutive word hits
    // on one diagonal merge into one match of length run + k - 1.
    private static IEnumerable<(int PositionInText, int PositionInQuery, int Length)> WordTableMatches(string a, string b, int k, bool[] allowed)
    {
        var table = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        for (int p = 0; p + k <= a.Length; p++)
        {
            if (!allowed[p])
                continue;
            string word = a.Substring(p, k);
            if (!table.TryGetValue(word, out var list))
                table[word] = list = new List<int>();
            list.Add(p);
        }

        var open = new Dictionary<int, (int StartP, int StartQ, int LastQ)>();
        var done = new List<(int PositionInText, int PositionInQuery, int Length)>();
        for (int q = 0; q + k <= b.Length; q++)
        {
            if (!table.TryGetValue(b.Substring(q, k), out var positions))
                continue;
            foreach (int p in positions)
            {
                int diag = q - p;
                if (open.TryGetValue(diag, out var cur))
                {
                    if (cur.LastQ == q - 1)
                    {
                        open[diag] = cur with { LastQ = q };
                        continue;
                    }

                    done.Add((cur.StartP, cur.StartQ, cur.LastQ - cur.StartQ + k));
                }

                open[diag] = (p, q, q);
            }
        }

        foreach (var cur in open.Values)
            done.Add((cur.StartP, cur.StartQ, cur.LastQ - cur.StartQ + k));
        return done;
    }

    // Upper-case query in which every non-letter and the skipped symbol is the QueryBreak symbol (equal to nothing in a tree).
    private static string LetterQuery(string sequence, char skip)
    {
        var chars = sequence.ToUpperInvariant().ToCharArray();
        for (int i = 0; i < chars.Length; i++)
            if (chars[i] is < 'A' or > 'Z' || chars[i] == skip)
                chars[i] = QueryBreak;
        return new string(chars);
    }


    #endregion
}
