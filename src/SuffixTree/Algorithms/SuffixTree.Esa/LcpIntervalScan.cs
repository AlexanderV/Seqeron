namespace SuffixTree.Esa;

/// <summary>A forward-only source of the LCP values <c>LCP[1], LCP[2], ..., LCP[n]</c> (<c>LCP[1] = 0</c>; ranks start at 1 because rank 0 is the empty suffix).</summary>
internal interface ILcpStream
{
    /// <summary>The next value; the caller asks exactly <c>n</c> times.</summary>
    int Next();
}

/// <summary>Events of the bottom-up scan of the lcp-interval tree (<see cref="LcpIntervalScan"/>). The sinks of the child-table layouts implement it.</summary>
internal interface IIntervalSink
{
    /// <summary>
    /// <paramref name="k"/> is an l-index of the interval that is open at the scan position; <paramref name="previous"/> is its previous l-index
    /// (the link <c>next[previous] = k</c>), or -1 when <paramref name="k"/> is the first l-index of that interval.
    /// </summary>
    void LIndex(int previous, int k);

    /// <summary>
    /// The internal interval <c>[lb..rb]</c> with value <paramref name="ell"/> and first l-index <paramref name="first"/> is closed (children are
    /// closed before their parent; the root is closed last). <paramref name="lastChild"/> is true when the interval is the last child of its parent
    /// (its right boundary is the parent's) and for the root.
    /// </summary>
    void Close(int ell, int lb, int rb, int first, bool lastChild);
}

/// <summary>An <see cref="ILcpStream"/> over a dense array of <c>n + 1</c> entries (starts at rank 1).</summary>
internal struct ArrayLcpStream(int[] lcp) : ILcpStream
{
    private int _next = 1;

    public int Next() => lcp[_next++];
}

/// <summary>An <see cref="ILcpStream"/> over a <see cref="CompactLcp"/> (starts at rank 1; the escape list is read by a moving cursor).</summary>
internal struct CompactLcpStream : ILcpStream
{
    private CompactLcp.Cursor _cursor;
    private int _next;

    internal CompactLcpStream(CompactLcp lcp)
    {
        _cursor = lcp.CreateCursor();
        _next = 1;
    }

    public int Next() => _cursor.Get(_next++);
}

/// <summary>
/// The sequential scan that every child-table layout is built from (Abouelhoda, Kurtz, Ohlebusch 2004, bottom-up traversal of the lcp-interval
/// tree): it reads the LCP values once, in rank order, keeps a stack of the open intervals (strictly increasing value) and reports the l-indices and
/// the closed intervals to a sink. It needs no suffix array, no text and no random access, only the stream, so a streaming builder (ADR-004, phase
/// PH5) can run it over a sequentially produced LCP file; the in-memory version keeps the stack in RAM, whose height is at most
/// <c>1 + max LCP</c> and at most n (a^n), 16 bytes per entry.
/// </summary>
internal static partial class LcpIntervalScan
{
    private const int InitialCapacity = 64;

    /// <summary>
    /// Scans the <paramref name="n"/> values of <paramref name="lcp"/> and feeds <paramref name="sink"/>. Returns the peak height of the stack (0 for n = 0, which has
    /// no internal interval).
    /// </summary>
    /// <exception cref="ArgumentException">The stream starts with a value other than 0 or contains a negative value (a corrupt LCP array).</exception>
    internal static int Run<TLcp, TSink>(int n, ref TLcp lcp, ref TSink sink)
        where TLcp : struct, ILcpStream
        where TSink : struct, IIntervalSink
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (n == 0) return 0;

        int capacity = InitialCapacity;
        int[] ell = new int[capacity], lb = new int[capacity], first = new int[capacity], last = new int[capacity];
        int sp = 1;                     // the root: value 0, left boundary 0, no l-index yet
        first[0] = -1;
        last[0] = -1;
        int peak = 1;
        for (int i = 1; i <= n + 1; i++)
        {
            int cur = i <= n ? lcp.Next() : -1;
            if (i <= n && (cur < 0 || (i == 1 && cur != 0)))
                throw new ArgumentException($"LCP[{i}] = {cur} is not a valid LCP value (LCP[1] must be 0 and no value may be negative)");

            int left = i - 1;
            while (sp > 0 && ell[sp - 1] > cur)
            {
                sp--;
                int b = lb[sp];
                sink.Close(ell[sp], b, i - 1, first[sp], sp == 0 || cur < ell[sp - 1]);
                left = b;
            }
            if (i > n) break;

            int top = sp - 1;
            if (ell[top] < cur)
            {
                if (sp == capacity)
                {
                    capacity *= 2;
                    Array.Resize(ref ell, capacity);
                    Array.Resize(ref lb, capacity);
                    Array.Resize(ref first, capacity);
                    Array.Resize(ref last, capacity);
                }
                ell[sp] = cur;
                lb[sp] = left;
                first[sp] = i;
                last[sp] = i;
                sp++;
                if (sp > peak) peak = sp;
                sink.LIndex(-1, i);
            }
            else if (first[top] < 0)
            {
                first[top] = i;
                last[top] = i;
                sink.LIndex(-1, i);
            }
            else
            {
                sink.LIndex(last[top], i);
                last[top] = i;
            }
        }
        return peak;
    }
}

/// <summary>A forward-only source of the LCP values <c>LCP[1], ..., LCP[n]</c> at 64 bits (a text above 2^31 symbols has values above <see cref="int.MaxValue"/>).</summary>
internal interface IWideLcpStream
{
    /// <summary>The next value; the caller asks exactly <c>n</c> times.</summary>
    long Next();
}

/// <summary>Events of the bottom-up scan at 64 bits (<see cref="LcpIntervalScan.RunWide{TLcp, TSink}"/>); the meaning of the arguments is that of <see cref="IIntervalSink"/>.</summary>
internal interface IWideIntervalSink
{
    /// <summary><paramref name="k"/> is an l-index of the open interval; <paramref name="previous"/> is its previous l-index, or -1 when <paramref name="k"/> is the first.</summary>
    void LIndex(long previous, long k);

    /// <summary>The internal interval <c>[lb..rb]</c> of value <paramref name="ell"/> and first l-index <paramref name="first"/> is closed (children before their parent, the root last).</summary>
    void Close(long ell, long lb, long rb, long first, bool lastChild);
}

/// <summary>An <see cref="IWideLcpStream"/> over a <see cref="CompactLcp"/> of any kind and size (starts at rank 1).</summary>
internal struct WideCompactLcpStream : IWideLcpStream
{
    private CompactLcp.LongCursor _cursor;
    private long _next;

    internal WideCompactLcpStream(CompactLcp lcp)
    {
        _cursor = lcp.CreateLongCursor();
        _next = 1;
    }

    /// <inheritdoc />
    public long Next() => _cursor.Get(_next++);
}

/// <summary>Counts the intervals and finds the deepest one while the wide scan runs: the numbers of <see cref="StatisticsSink"/> at 64 bits.</summary>
internal struct WideStatisticsSink : IWideIntervalSink
{
    internal long Internal;
    internal long LrsDepth;
    internal long LrsRank;

    /// <inheritdoc />
    public readonly void LIndex(long previous, long k)
    {
    }

    /// <inheritdoc />
    public void Close(long ell, long lb, long rb, long first, bool lastChild)
    {
        if (ell <= 0) return;
        Internal++;
        if (ell > LrsDepth)
        {
            LrsDepth = ell;
            LrsRank = lb;
        }
    }
}

/// <summary>The scan at 64 bits (work package S2.4): <see cref="Run{TLcp, TSink}"/> with every rank, depth and count a <see cref="long"/>, for an LCP array of any size and kind. The 32-bit scan stays for the in-memory builder, whose arrays are limited to 32 bits.</summary>
internal static partial class LcpIntervalScan
{

    /// <summary>
    /// Scans the <paramref name="n"/> values of <paramref name="lcp"/> and feeds <paramref name="sink"/>; returns the peak height of the stack (0 for n = 0). The stack holds one entry per open interval, at most <c>1 + max LCP</c> of them
    /// (<c>a^n</c> has <c>n</c>), 32 bytes each: a text whose stack does not fit an array is refused.
    /// </summary>
    /// <exception cref="ArgumentException">The stream starts with a value other than 0 or contains a negative value (a corrupt LCP array).</exception>
    /// <exception cref="NotSupportedException">The stack of open intervals would hold more entries than an array.</exception>
    internal static long RunWide<TLcp, TSink>(long n, ref TLcp lcp, ref TSink sink)
        where TLcp : struct, IWideLcpStream
        where TSink : struct, IWideIntervalSink
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (n == 0) return 0;

        int capacity = InitialCapacity;
        long[] ell = new long[capacity], lb = new long[capacity], first = new long[capacity], last = new long[capacity];
        int sp = 1;                     // the root: value 0, left boundary 0, no l-index yet
        first[0] = -1;
        last[0] = -1;
        long peak = 1;
        for (long i = 1; i <= n + 1; i++)
        {
            long cur = i <= n ? lcp.Next() : -1;
            if (i <= n && (cur < 0 || (i == 1 && cur != 0)))
                throw new ArgumentException($"LCP[{i}] = {cur} is not a valid LCP value (LCP[1] must be 0 and no value may be negative)");

            long left = i - 1;
            while (sp > 0 && ell[sp - 1] > cur)
            {
                sp--;
                long b = lb[sp];
                sink.Close(ell[sp], b, i - 1, first[sp], sp == 0 || cur < ell[sp - 1]);
                left = b;
            }
            if (i > n) break;

            int top = sp - 1;
            if (ell[top] < cur)
            {
                if (sp == capacity)
                {
                    if (capacity > Array.MaxLength / 2) throw new NotSupportedException($"the stack of the open lcp-intervals holds more than {capacity} entries: a text this deep is not scanned (LcpIntervalScan.RunWide)");
                    capacity *= 2;
                    Array.Resize(ref ell, capacity);
                    Array.Resize(ref lb, capacity);
                    Array.Resize(ref first, capacity);
                    Array.Resize(ref last, capacity);
                }
                ell[sp] = cur;
                lb[sp] = left;
                first[sp] = i;
                last[sp] = i;
                sp++;
                if (sp > peak) peak = sp;
                sink.LIndex(-1, i);
            }
            else if (first[top] < 0)
            {
                first[top] = i;
                last[top] = i;
                sink.LIndex(-1, i);
            }
            else
            {
                sink.LIndex(last[top], i);
                last[top] = i;
            }
        }
        return peak;
    }

    /// <summary>The header statistics of an LCP array of any size and kind (the formulas of <see cref="EsaStatistics"/>), by one wide scan: the numbers that the header of a container records, computed.</summary>
    internal static EsaStatistics StatisticsOf(CompactLcp lcp)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        long n = lcp.Length - 1;
        var stream = new WideCompactLcpStream(lcp);
        var sink = new WideStatisticsSink();
        long peak = RunWide(n, ref stream, ref sink);
        return n == 0
            ? new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0)
            : new EsaStatistics(n, 1 + sink.Internal + (n + 1), sink.Internal, n, n, sink.LrsDepth, sink.LrsRank, peak);
    }
}
