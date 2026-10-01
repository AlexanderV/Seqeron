namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Selects the symmetric-DUST implementation used by
/// <see cref="SequenceComplexity.MaskLowComplexity(string, int, double, char, int, bool, DustEngine)"/> and
/// <see cref="SequenceComplexity.FindLowComplexityIntervals(string, int, double, int, DustEngine)"/>.
/// </summary>
public enum DustEngine
{
    /// <summary>
    /// lh3/sdust <c>sdust_core</c> port (default; the historical behaviour). Every non-ACGT symbol
    /// splits the input into independently scanned ACGT runs. Window ≥ 3, any threshold ≥ 0.
    /// </summary>
    Sdust = 0,

    /// <summary>
    /// NCBI <c>dustmasker</c> parity: a port of <c>CSymDustMasker</c> (algo/dustmask/symdust.cpp) driven by
    /// <c>GetDustMasks_SkipNs</c> (app/dustmasker/dust_mask_app.cpp). IUPAC codes are scanned as bases
    /// (C→1, G→2, T→3, N→pseudo-random 2-bit code from the toolkit's <c>CRandom</c>, every other code→0 = A);
    /// only N runs longer than the window — and leading/trailing N runs of any length — cut the scan, and
    /// those N runs are themselves reported as masked. Restricted to dustmasker's accepted parameter
    /// ranges: window 8–64, level = 10·threshold an integer 2–64, linker 1–32; input must be IUPAC DNA.
    /// </summary>
    Dustmasker = 1,
}

public static partial class SequenceComplexity
{
    // symdust constructor ranges: window_ 8..64, level_ 2..64 (outside them dustmasker silently
    // substitutes the defaults 64 / 20; here they are rejected, as the linker already is — B04 F35).
    private const int DustmaskerMinWindow = 8;
    private const int DustmaskerMaxWindow = 64;
    private const int DustmaskerMinLevel = 2;
    private const int DustmaskerMaxLevel = 64;

    /// <summary>
    /// Parses a DUST engine name, case-insensitively: "sdust" (also null/empty = default) or "dustmasker".
    /// Used by the MCP wrappers' optional <c>engine</c> parameter.
    /// </summary>
    /// <param name="name">Engine name.</param>
    /// <returns>The engine.</returns>
    /// <exception cref="ArgumentException">Thrown for any other name.</exception>
    public static DustEngine ParseDustEngine(string? name) =>
        (name ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "" or "sdust" => DustEngine.Sdust,
            "dustmasker" => DustEngine.Dustmasker,
            _ => throw new ArgumentException($"Unknown DUST engine '{name}' (expected 'sdust' or 'dustmasker').", nameof(name)),
        };

    private static int ValidateDustmaskerParameters(string sequence, int windowSize, double threshold)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, DustmaskerMinWindow);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(windowSize, DustmaskerMaxWindow);
        double level = threshold * 10.0;
        double rounded = Math.Round(level);
        if (Math.Abs(level - rounded) > 1e-9 || rounded < DustmaskerMinLevel || rounded > DustmaskerMaxLevel)
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold,
                "The dustmasker engine takes an integer level = 10·threshold in 2–64 (threshold 0.2–6.4 in steps of 0.1).");
        if (!global::Seqeron.Genomics.Core.SequenceExtensions.IsValidIupacDna(sequence.AsSpan()))
            throw new ArgumentException(
                "The dustmasker engine accepts IUPAC DNA codes only (A C G T R Y S W K M B D H V N, any case).", nameof(sequence));
        return (int)rounded;
    }

    /// <summary>
    /// dustmasker <c>GetDustMasks_SkipNs</c> on an upper-case IUPAC DNA string. Returns 0-based
    /// half-open intervals (dustmasker's closed <c>a - b</c> with End = b + 1), in dustmasker's order.
    /// </summary>
    private static List<(int Start, int End)> FindDustmaskerIntervals(string seq, int window, int level, int linker)
    {
        var closed = SymDustMasker.GetDustMasksSkipNs(seq, (uint)level, window, linker);
        var res = new List<(int Start, int End)>(closed.Count);
        foreach (var (first, second) in closed) res.Add(((int)first, (int)second + 1));
        return res;
    }

    /// <summary>
    /// Line-by-line port of NCBI C++ Toolkit <c>CSymDustMasker</c> (algo/dustmask/symdust.cpp, symdust.hpp)
    /// and of <c>s_FindSegmentWithLongNs</c> / <c>s_InsertMerge</c> / <c>GetDustMasks_SkipNs</c>
    /// (app/dustmasker/dust_mask_app.cpp). US-Government work, public domain. Unsigned C++ arithmetic is kept
    /// where it matters (positions are long, scores uint). Coordinates are closed, as in the toolkit.
    /// </summary>
    private sealed class SymDustMasker
    {
        private const int TripletMask = 0x3F;

        private readonly long _window;
        private readonly long _linker;
        private readonly uint _lowK;
        private readonly uint[] _thresholds;
        private readonly List<Perfect> _p = new();
        private readonly NcbiRandom _random = new(); // convert_t::m_Random (default CRandom = LFG, Reset())

        private SymDustMasker(uint level, long window, long linker)
        {
            _window = window;
            _linker = linker;
            _lowK = level / 5;

            // thresholds_: [0] = 1, [i] = i·level for i = 1 .. window−3 (window − 2 entries).
            _thresholds = new uint[_window - 2];
            _thresholds[0] = 1;
            for (long i = 1; i < _window - 2; ++i) _thresholds[i] = (uint)i * level;
        }

        private struct Perfect
        {
            public long First, Second; // bounds_ (closed, relative to the scan start)
            public uint Score;
            public long Len;

            public Perfect(long first, long second, uint score, long len)
            {
                First = first;
                Second = second;
                Score = score;
                Len = len;
            }
        }

        // CIupac2Ncbi2na_converter: C→1, G→2, T→3, N→m_Random.GetRand() & 3, anything else (A, IUPAC)→0.
        private int Convert(char r) => r switch
        {
            'C' => 1,
            'G' => 2,
            'T' => 3,
            'N' => (int)(_random.GetRand() & 0x3),
            _ => 0,
        };

        /// <summary>dust_mask_app.cpp <c>GetDustMasks_SkipNs</c>.</summary>
        public static List<(long First, long Second)> GetDustMasksSkipNs(string seq, uint level, int window, int linker)
        {
            var duster = new SymDustMasker(level, window, linker);
            var nsRange = FindSegmentWithLongNs(window, seq);

            if (nsRange.Count == 0) return duster.Mask(seq, 0, seq.Length - 1);

            var rv = new List<(long First, long Second)>();
            long seqStart = 0;
            foreach (var itr in nsRange)
            {
                if (itr.First == 0)
                {
                    seqStart = itr.Second + 1;
                    rv.Add(itr);
                    continue;
                }

                var sMask = duster.Mask(seq, seqStart, itr.First - 1);
                if (sMask.Count > 0)
                {
                    InsertMerge(rv, sMask[0], linker);
                    for (int k = 1; k < sMask.Count; k++) rv.Add(sMask[k]);
                    InsertMerge(rv, itr, linker);
                }
                else
                {
                    rv.Add(itr);
                }
                seqStart = itr.Second + 1;
            }

            if (seqStart < seq.Length)
            {
                var sMask = duster.Mask(seq, seqStart, seq.Length - 1);
                if (sMask.Count > 0)
                {
                    InsertMerge(rv, sMask[0], linker);
                    for (int k = 1; k < sMask.Count; k++) rv.Add(sMask[k]);
                }
            }
            return rv;
        }

        // s_FindSegmentWithLongNs: N runs longer than maxNs, plus leading and trailing N runs of any length.
        private static List<(long First, long Second)> FindSegmentWithLongNs(long maxNs, string seq)
        {
            var nsRange = new List<(long First, long Second)>();
            long pos = 0, ns = 0;
            foreach (char c in seq)
            {
                if (c == 'N')
                {
                    ns++;
                }
                else
                {
                    if (ns > 0)
                    {
                        if (ns > maxNs || pos == 0) nsRange.Add((pos, pos + ns - 1));
                        pos += ns;
                        ns = 0;
                    }
                    pos++;
                }
            }
            if (ns > 0) nsRange.Add((pos, pos + ns - 1));
            return nsRange;
        }

        // s_InsertMerge: joins only when the previous end + linker equals the new start exactly.
        private static void InsertMerge(List<(long First, long Second)> list, (long First, long Second) newMask, long linker)
        {
            if (list.Count > 0 && list[^1].Second + linker == newMask.First)
            {
                list[^1] = (list[^1].First, newMask.Second);
                return;
            }
            list.Add(newMask);
        }

        /// <summary>symdust.cpp <c>operator()(seq, start, stop)</c>.</summary>
        private List<(long First, long Second)> Mask(string seq, long start, long stop)
        {
            var res = new List<(long First, long Second)>();
            if (seq.Length == 0) return res;
            if (stop >= seq.Length) stop = seq.Length - 1;
            if (start > stop) start = stop;

            while (stop > 2 + start) // there must be at least one triplet
            {
                _p.Clear();
                var w = new Triplets(_window, _lowK, _p, _thresholds);

                char c1 = seq[(int)start], c2 = seq[(int)start + 1];
                int hi = Convert(c1);
                int t = (hi << 2) + Convert(c2);

                long it = start + w.Stop + 2;
                bool done = false;
                while (!done && it <= stop)
                {
                    SaveMaskedRegions(res, w.Start, start);

                    // shift the window
                    t = ((t << 2) & TripletMask) + (Convert(seq[(int)it]) & 0x3);
                    ++it;

                    if (w.ShiftWindow(t))
                    {
                        if (w.NeedsProcessing()) w.FindPerfect();
                    }
                    else
                    {
                        while (it <= stop)
                        {
                            SaveMaskedRegions(res, w.Start, start);
                            t = ((t << 2) & TripletMask) + (Convert(seq[(int)it]) & 0x3);

                            if (w.ShiftWindow(t))
                            {
                                done = true;
                                break;
                            }
                            ++it;
                        }
                    }
                }

                // append the rest of the perfect intervals to the result
                long wstart = w.Start;
                while (_p.Count > 0)
                {
                    SaveMaskedRegions(res, wstart, start);
                    ++wstart;
                }

                if (w.Start > 0) start += w.Start;
                else break;
            }

            return res;
        }

        private void SaveMaskedRegions(List<(long First, long Second)> res, long wstart, long start)
        {
            if (_p.Count == 0) return;

            var b = _p[^1];
            if (b.First < wstart)
            {
                var b1 = (First: b.First + start, Second: b.Second + start);
                if (res.Count > 0)
                {
                    long s = res[^1].Second;
                    if (s + _linker >= b1.First) res[^1] = (res[^1].First, Math.Max(s, b1.Second));
                    else res.Add(b1);
                }
                else
                {
                    res.Add(b1);
                }

                while (_p.Count > 0 && _p[^1].First < wstart) _p.RemoveAt(_p.Count - 1);
            }
        }

        /// <summary>symdust <c>triplets</c>: the window's triplet deque (index 0 = newest) and its counts.</summary>
        private sealed class Triplets
        {
            private readonly int[] _buf;    // circular deque storage
            private int _head, _count;      // _head = index of the front (newest) element
            private readonly long _maxSize;
            private readonly uint _lowK;
            private readonly List<Perfect> _p;
            private readonly uint[] _thresholds;
            private readonly int[] _cw = new int[64];
            private readonly int[] _cv = new int[64];
            private readonly int[] _counts = new int[64];
            private uint _rw, _rv, _numDiff;
            private long _l;

            public Triplets(long window, uint lowK, List<Perfect> perfectList, uint[] thresholds)
            {
                _maxSize = window - 2;
                _buf = new int[_maxSize + 1];
                _lowK = lowK;
                _p = perfectList;
                _thresholds = thresholds;
            }

            public long Start { get; private set; }
            public long Stop { get; private set; }

            private int this[long index] => _buf[(int)((_head + index) % _buf.Length)];

            private void PushFront(int t)
            {
                _head = (_head - 1 + _buf.Length) % _buf.Length;
                _buf[_head] = t;
                _count++;
            }

            private int PopBack()
            {
                int v = _buf[(_head + _count - 1) % _buf.Length];
                _count--;
                return v;
            }

            private static void AddTripletInfo(ref uint r, int[] c, int t) { r += (uint)c[t]; ++c[t]; }

            private static void RemTripletInfo(ref uint r, int[] c, int t) { --c[t]; r -= (uint)c[t]; }

            public bool NeedsProcessing()
            {
                long count = Stop - _l;
                return count < _count && 10 * _rw > _thresholds[count];
            }

            private bool ShiftHigh(int t)
            {
                int s = PopBack();
                RemTripletInfo(ref _rw, _cw, s);
                if (_cw[s] == 0) --_numDiff;
                ++Start;

                PushFront(t);
                if (_cw[t] == 0) ++_numDiff;
                AddTripletInfo(ref _rw, _cw, t);
                ++Stop;

                if (_numDiff <= 1)
                {
                    _p.Insert(0, new Perfect(Start, Stop + 1, 0, 0));
                    return false;
                }
                return true;
            }

            public bool ShiftWindow(int t)
            {
                // shift the left end of the window, if necessary
                if (_count >= _maxSize)
                {
                    if (_numDiff <= 1) return ShiftHigh(t);

                    int s = PopBack();
                    RemTripletInfo(ref _rw, _cw, s);
                    if (_cw[s] == 0) --_numDiff;

                    if (_l == Start)
                    {
                        ++_l;
                        RemTripletInfo(ref _rv, _cv, s);
                    }

                    ++Start;
                }

                PushFront(t);
                if (_cw[t] == 0) ++_numDiff;
                AddTripletInfo(ref _rw, _cw, t);
                AddTripletInfo(ref _rv, _cv, t);

                if (_cv[t] > _lowK)
                {
                    long off = _count - (_l - Start) - 1;
                    int u;
                    do
                    {
                        u = this[off];
                        RemTripletInfo(ref _rv, _cv, u);
                        ++_l;
                        off--;
                    } while (u != t);
                }

                ++Stop;

                if (_count >= _maxSize && _numDiff <= 1)
                {
                    _p.Clear();
                    _p.Insert(0, new Perfect(Start, Stop + 1, 0, 0));
                    return false;
                }
                return true;
            }

            public void FindPerfect()
            {
                long count = Stop - _l; // the suffix length
                Array.Copy(_cv, _counts, 64);

                uint score = _rv;
                int perfectIter = 0;
                uint maxPerfectScore = 0;
                long maxLen = 0;
                long pos = _l - 1; // skipping the suffix

                for (long idx = count; idx < _count; ++idx, ++count, --pos)
                {
                    int tt = this[idx];
                    int cnt = _counts[tt];
                    AddTripletInfo(ref score, _counts, tt);

                    if (cnt > 0 && score * 10 > _thresholds[count])
                    {
                        // found the candidate for the perfect interval: get the max score for the
                        // existing perfect intervals within the current suffix
                        while (perfectIter < _p.Count && pos <= _p[perfectIter].First)
                        {
                            var pi = _p[perfectIter];
                            if (maxPerfectScore == 0 || maxLen * pi.Score > maxPerfectScore * pi.Len)
                            {
                                maxPerfectScore = pi.Score;
                                maxLen = pi.Len;
                            }
                            ++perfectIter;
                        }

                        // check if the current suffix score is at least as large
                        if (maxPerfectScore == 0 || score * maxLen >= maxPerfectScore * count)
                        {
                            maxPerfectScore = score;
                            maxLen = count;
                            _p.Insert(perfectIter, new Perfect(pos, Stop + 1, maxPerfectScore, count));
                        }
                    }
                }
            }
        }

        /// <summary>
        /// NCBI C++ Toolkit <c>CRandom</c> (util/random_gen.cpp/.hpp), default LFG method: lagged Fibonacci,
        /// lags 33/13, modulus 2^32 on the state, <c>GetRand() = x_GetRand32Bits() &gt;&gt; 1</c>, seeded by
        /// <c>Reset()</c> with the fixed <c>sm_State</c> table (so dustmasker's N substitution is deterministic).
        /// </summary>
        private sealed class NcbiRandom
        {
            private const int StateSize = 33;
            private const int StateOffset = 12;

            private static readonly uint[] SmState =
            {
                0xd53f1852, 0xdfc78b83, 0x4f256096, 0x0e643df7,
                0x82c359bf, 0xc7794dfa, 0xd5e9ffaa, 0x2c8cb64a,
                0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589,
                0xa5853646, 0xe71576e2, 0x0dae30df, 0xb09ce711,
                0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8,
                0xd788d765, 0x67498087, 0x9d7aba26, 0x261351d4,
                0x411ee7ea, 0x0393a263, 0x2c5a5835, 0xc115fcd8,
                0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98,
                0x6e37bd55,
            };

            private readonly uint[] _state = (uint[])SmState.Clone();
            private int _rj = StateOffset;
            private int _rk = StateSize - 1;

            public uint GetRand()
            {
                uint r = unchecked(_state[_rk] + _state[_rj--]);
                _state[_rk--] = r;
                if (_rk < 0) _rk = StateSize - 1;
                else if (_rj < 0) _rj = StateSize - 1;
                return r >> 1; // discard the least-random bit
            }
        }
    }
}
