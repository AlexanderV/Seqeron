namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The suffix automaton (DAWG, Blumer et al. 1985; Crochemore 1985; construction as in Crochemore &amp; Rytter 2003 §5.3 and cp-algorithms "Suffix Automaton") of the REVERSE of a text, used for the matching statistics of
/// a query against the text (audit R14-4a, B09 F291): reading the reverse of the query through the automaton of the reverse of the text gives, after the character that is query[i], the length ms(i) of the longest prefix of
/// query[i..] that occurs in the text (Chang &amp; Lawler 1994; Gusfield 1997 §7.8) in O(|text| + |query|) time, together with the automaton state of that prefix. Equal strings are in the same state at the same length, so a
/// (state, length) pair identifies a distinct substring, and the largest end position of a state (propagated up the suffix links) is the lowest start of the string in the text: p = |text| − 1 − maxEnd.
/// Memory: 2 |text| states × (alphabet + 3) ints (80 MB for 2 Mbp; the suffix tree it replaces here took 540 MB).
/// </summary>
/// <remarks>
/// Symbols are byte codes 0 .. alphabet − 1 of the text; a query symbol outside the alphabet (code −1) ends every match. Nothing in this class is specific to nucleotides.
/// </remarks>
internal sealed class MatchingStatisticsAutomaton
{
    private readonly int _alphabet;
    private readonly int _textLength;
    private readonly int[] _next;       // state * alphabet + symbol -> state; 0 (the root, never a target) means "no transition"
    private readonly int[] _link;
    private readonly int[] _len;
    private readonly int[] _maxEnd;     // largest end position of the state in the REVERSED text

    /// <summary>Builds the automaton of the reverse of <paramref name="text"/> (symbol codes 0 .. <paramref name="alphabet"/> − 1).</summary>
    public MatchingStatisticsAutomaton(ReadOnlySpan<byte> text, int alphabet, CancellationToken cancellationToken = default)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(alphabet, 1);
        _alphabet = alphabet;
        _textLength = text.Length;
        int capacity = Math.Max(2, 2 * text.Length);
        _next = new int[checked(capacity * alphabet)];
        _link = new int[capacity];
        _len = new int[capacity];
        _maxEnd = new int[capacity];
        _link[0] = -1;
        _maxEnd[0] = -1;
        int size = 1;
        int last = 0;
        int n = text.Length;
        for (int k = 0; k < n; k++)
        {
            if ((k & 0x3FFF) == 0)
                cancellationToken.ThrowIfCancellationRequested();
            int c = text[n - 1 - k];            // the reversed text
            int cur = size++;
            _len[cur] = _len[last] + 1;
            _maxEnd[cur] = k;
            int p = last;
            while (p != -1 && _next[(p * _alphabet) + c] == 0)
            {
                _next[(p * _alphabet) + c] = cur;
                p = _link[p];
            }

            if (p == -1)
            {
                _link[cur] = 0;
            }
            else
            {
                int q = _next[(p * _alphabet) + c];
                if (_len[p] + 1 == _len[q])
                {
                    _link[cur] = q;
                }
                else
                {
                    int clone = size++;
                    _len[clone] = _len[p] + 1;
                    Array.Copy(_next, q * _alphabet, _next, clone * _alphabet, _alphabet);
                    _link[clone] = _link[q];
                    _maxEnd[clone] = -1;
                    while (p != -1 && _next[(p * _alphabet) + c] == q)
                    {
                        _next[(p * _alphabet) + c] = clone;
                        p = _link[p];
                    }

                    _link[q] = _link[cur] = clone;
                }
            }

            last = cur;
        }

        // endpos(state) = the end positions of its suffix-link children plus its own: the largest end position moves up the links, longest states first (counting sort by len)
        var count = new int[n + 2];
        for (int v = 0; v < size; v++)
            count[_len[v]]++;
        for (int i = 1; i < count.Length; i++)
            count[i] += count[i - 1];
        var order = new int[size];
        for (int v = size - 1; v >= 0; v--)
            order[--count[_len[v]]] = v;
        for (int i = size - 1; i > 0; i--)
        {
            int v = order[i];
            int parent = _link[v];
            if (_maxEnd[v] > _maxEnd[parent])
                _maxEnd[parent] = _maxEnd[v];
        }
    }

    /// <summary>
    /// The matching statistics of <paramref name="query"/>: <c>lengths[i]</c> = the length of the longest prefix of query[i..] that occurs in the text, <c>states[i]</c> = the automaton state of that prefix
    /// (the pair (state, length) identifies the string). A query code of −1 (a symbol that matches nothing) gives 0 there and ends every match across it.
    /// </summary>
    public (int[] Lengths, int[] States) Match(ReadOnlySpan<sbyte> query, CancellationToken cancellationToken = default)
    {
        int m = query.Length;
        var lengths = new int[m];
        var states = new int[m];
        int v = 0, l = 0;
        for (int i = m - 1; i >= 0; i--)
        {
            if ((i & 0x3FFF) == 0)
                cancellationToken.ThrowIfCancellationRequested();
            int c = query[i];
            if (c < 0 || c >= _alphabet)
            {
                v = 0;
                l = 0;
            }
            else
            {
                while (v != 0 && _next[(v * _alphabet) + c] == 0)
                {
                    v = _link[v];
                    l = _len[v];
                }

                int t = _next[(v * _alphabet) + c];
                if (t != 0)
                {
                    v = t;
                    l++;
                }
            }

            lengths[i] = l;
            states[i] = v;
        }

        return (lengths, states);
    }

    /// <summary>The lowest 0-based start in the text of the strings of <paramref name="state"/> (valid for a non-root state).</summary>
    public int FirstOccurrence(int state) => _textLength - 1 - _maxEnd[state];
}
