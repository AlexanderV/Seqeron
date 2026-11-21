namespace Seqeron.Genomics.Alignment;

/// <summary>The alignment mode of SPOA (<c>spoa::AlignmentType</c>): local (Smith-Waterman), global (Needleman-Wunsch) or semi-global (overlap).</summary>
public enum SpoaAlignmentType
{
    /// <summary>Smith-Waterman.</summary>
    Sw = 0,

    /// <summary>Needleman-Wunsch.</summary>
    Nw = 1,

    /// <summary>Overlap (semi-global, <c>-l 2</c>).</summary>
    Ov = 2,
}

/// <summary>The gap model of SPOA (<c>spoa::AlignmentSubtype</c>).</summary>
public enum SpoaAlignmentSubtype
{
    /// <summary>g * i.</summary>
    Linear,

    /// <summary>g + (i - 1) * e.</summary>
    Affine,

    /// <summary>The better of two affine functions.</summary>
    Convex,
}

/// <summary>
/// The partial order alignment engine of SPOA 4.0.8 as the installed program runs it: the SIMD engine (<c>simd_alignment_engine_implementation.hpp</c>) emulated lane by lane - the int16 or int32 variant
/// chosen by the worst-case score, vectors of <c>registerBits / (16 or 32)</c> lanes whose padding lanes carry <c>padding_penalty</c> and take part in the row maxima, wrapping lane arithmetic, the shifted
/// prefix maxima of <c>_mmxxx_prefix_max</c> and the traceback of the program (including <c>_mmxxx_max_value</c>, which starts at 0, and <c>_mmxxx_index_of</c>). The sequential <c>SisdAlignmentEngine</c> differs from it
/// on semi-global alignments and on rows whose maximum is in a padding lane, so it is not used.
/// </summary>
public sealed class SpoaEngine
{
    private readonly int _m, _n, _g, _e, _q, _c;
    private readonly int _registerBits;

    /// <summary>The alignment mode.</summary>
    public SpoaAlignmentType Type { get; }

    /// <summary>The gap model derived from the penalties.</summary>
    public SpoaAlignmentSubtype Subtype { get; }

    private SpoaEngine(SpoaAlignmentType type, SpoaAlignmentSubtype subtype, int m, int n, int g, int e, int q, int c, int registerBits)
    {
        Type = type;
        Subtype = subtype;
        _m = m;
        _n = n;
        _g = g;
        _e = e;
        _q = q;
        _c = c;
        _registerBits = registerBits;
    }

    /// <summary><c>AlignmentEngine::Create</c>: the parameters are <c>std::int8_t</c> values; a linear model sets e = g, an affine model sets q = g and c = e.</summary>
    /// <param name="registerBits">256 = AVX2 (default; what the installed program runs on a CPU with AVX2), 128 = SSE4.1 / SSE2.</param>
    /// <exception cref="ArgumentException"><c>std::invalid_argument</c> of SPOA.</exception>
    public static SpoaEngine Create(int type, int m, int n, int g, int e, int q, int c, int registerBits = 256)
    {
        if (type is not (0 or 1 or 2))
            throw new ArgumentException("[spoa::AlignmentEngine::Create] error: invalid alignment type!");
        if (g > 0 || q > 0)
            throw new ArgumentException("[spoa::AlignmentEngine::Create] error: gap opening penalty must be non-positive!");
        if (e > 0 || c > 0)
            throw new ArgumentException("[spoa::AlignmentEngine::Create] error: gap extension penalty must be non-positive!");
        SpoaAlignmentSubtype subtype;
        if (g >= e)
            subtype = SpoaAlignmentSubtype.Linear;
        else if (g <= q || e >= c)
            subtype = SpoaAlignmentSubtype.Affine;
        else
            subtype = SpoaAlignmentSubtype.Convex;
        if (subtype == SpoaAlignmentSubtype.Linear)
        {
            e = g;
        }
        else if (subtype == SpoaAlignmentSubtype.Affine)
        {
            q = g;
            c = e;
        }

        return new SpoaEngine((SpoaAlignmentType)type, subtype, m, n, g, e, q, c, registerBits);
    }

    /// <summary><c>WorstCaseAlignmentScore</c>.</summary>
    private long WorstCaseAlignmentScore(long i, long j)
    {
        long GapScore(long len) => len == 0 ? 0 : Math.Min(_g + ((len - 1) * _e), _q + ((len - 1) * _c));
        return Math.Min(-1 * ((_m * Math.Min(i, j)) + GapScore(Math.Abs(i - j))), GapScore(i) + GapScore(j));
    }

    /// <summary><c>Prealloc</c>: for an input without any residue the program's allocation size underflows (std::length_error).</summary>
    public static bool PreallocFails(uint maxSequenceLen) => maxSequenceLen == 0;

    /// <summary>Aligns a sequence to the graph (<c>AlignmentEngine::Align</c>); returns the (node id, sequence position) pairs and the score.</summary>
    /// <exception cref="ArgumentException">The score may overflow (<c>std::invalid_argument</c>).</exception>
    public List<(int Node, int Position)> Align(string sequence, SpoaGraph graph, out int score)
    {
        score = 0;
        uint sequenceLen = (uint)sequence.Length;
        if (graph.Nodes.Count == 0 || sequenceLen == 0)
            return [];
        long worst = WorstCaseAlignmentScore(sequenceLen + 8L, graph.Nodes.Count);
        if (worst < int.MinValue + 1024L)
            throw new ArgumentException("[spoa::SimdAlignmentEngine::Align] error: possible overflow!");
        bool wide = worst < short.MinValue + 1024L;
        var run = new Run(this, sequence, graph, wide);
        return run.Execute(out score);
    }

    private sealed class Run
    {
        private readonly SpoaEngine _s;
        private readonly string _sequence;
        private readonly SpoaGraph _graph;
        private readonly bool _wide;
        private readonly int _lanes;
        private readonly int _logLanes;
        private readonly int _negInf;
        private readonly long _normalWidth;
        private readonly int _width;
        private readonly int _height;
        private readonly int[] _rank;
        private readonly int[] _firstColumn;
        private int[] _h = [], _f = [], _e = [], _o = [], _q = [];
        private int[] _profile = [];
        private int[][] _masks = [];

        public Run(SpoaEngine engine, string sequence, SpoaGraph graph, bool wide)
        {
            _s = engine;
            _sequence = sequence;
            _graph = graph;
            _wide = wide;
            _lanes = engine._registerBits / (wide ? 32 : 16);
            _logLanes = (int)Math.Log2(_lanes);
            _negInf = wide ? int.MinValue + 1024 : short.MinValue + 1024;
            _normalWidth = sequence.Length;
            _width = (int)Math.Ceiling((double)sequence.Length / _lanes);
            _height = graph.Nodes.Count + 1;
            _rank = new int[graph.Nodes.Count];
            _firstColumn = new int[_height * 3];
        }

        // ------------------------------------------------------------------ lane arithmetic

        private int Wrap(int v) => _wide ? v : (short)v;

        private int[] Add(int[] a, int[] b)
        {
            var r = new int[_lanes];
            for (int i = 0; i < _lanes; i++)
                r[i] = Wrap(a[i] + b[i]);
            return r;
        }

        private int[] Sub(int[] a, int[] b)
        {
            var r = new int[_lanes];
            for (int i = 0; i < _lanes; i++)
                r[i] = Wrap(a[i] - b[i]);
            return r;
        }

        private int[] Max(int[] a, int[] b)
        {
            var r = new int[_lanes];
            for (int i = 0; i < _lanes; i++)
                r[i] = Math.Max(a[i], b[i]);
            return r;
        }

        private int[] Or(int[] a, int[] b)
        {
            var r = new int[_lanes];
            for (int i = 0; i < _lanes; i++)
                r[i] = Wrap(a[i] | b[i]);
            return r;
        }

        private int[] Set1(int v)
        {
            var r = new int[_lanes];
            Array.Fill(r, Wrap(v));
            return r;
        }

        /// <summary>Shift towards higher lanes (<c>_mmxxx_slli_si</c> by <paramref name="count"/> lanes), zero fill.</summary>
        private int[] Slli(int[] a, int count = 1)
        {
            var r = new int[_lanes];
            for (int i = count; i < _lanes; i++)
                r[i] = a[i - count];
            return r;
        }

        /// <summary>The last lane moved to lane 0 (<c>_mmxxx_srli_si(a, kRSS)</c>), zero elsewhere.</summary>
        private int[] SrliLast(int[] a)
        {
            var r = new int[_lanes];
            r[0] = a[_lanes - 1];
            return r;
        }

        private void PrefixMax(int[] a, int[][] penalties, int offset)
        {
            for (int k = 0; k < _logLanes; k++)
            {
                var shifted = Or(_masks[k], Slli(Add(a, penalties[offset + k]), 1 << k));
                var m = Max(a, shifted);
                Array.Copy(m, a, _lanes);
            }
        }

        private int[] Cell(int[] matrix, int row, int j) => matrix.AsSpan(((row * _width) + j) * _lanes, _lanes).ToArray();

        private void Store(int[] matrix, int row, int j, int[] v) => Array.Copy(v, 0, matrix, ((row * _width) + j) * _lanes, _lanes);

        private int At(int[] matrix, int row, long col) => matrix[(row * _width * _lanes) + (int)col];

        // ------------------------------------------------------------------ initialization

        private void Initialize()
        {
            var s = _s;
            int paddingPenalty = -1 * Math.Max(Math.Max(Math.Abs(s._m), Math.Abs(s._n)), Math.Max(Math.Abs(s._g), Math.Abs(s._q)));
            int codes = (int)_graph.NumCodes;
            _profile = new int[codes * _width * _lanes];
            for (int i = 0; i < codes; i++)
            {
                char c = _graph.Decoder((uint)i);
                for (int j = 0; j < _width; j++)
                {
                    for (int k = 0; k < _lanes; k++)
                    {
                        long col = ((long)j * _lanes) + k;
                        int value = paddingPenalty;
                        if (col < _normalWidth)
                            value = c == _sequence[(int)col] ? s._m : s._n;
                        _profile[(((i * _width) + j) * _lanes) + k] = value;
                    }
                }
            }

            var rankToNode = _graph.RankToNode;
            for (int i = 0; i < rankToNode.Count; i++)
                _rank[rankToNode[i].Id] = i;
            int size = _height * _width * _lanes;
            _h = new int[size];
            if (s.Subtype != SpoaAlignmentSubtype.Linear)
            {
                _f = new int[size];
                _e = new int[size];
            }

            if (s.Subtype == SpoaAlignmentSubtype.Convex)
            {
                _o = new int[size];
                _q = new int[size];
            }

            int h = _height;
            switch (s.Subtype)
            {
                case SpoaAlignmentSubtype.Convex:
                    for (int j = 0; j < _width; j++)
                    {
                        Store(_o, 0, j, Set1(_negInf));
                        var q = new int[_lanes];
                        for (int k = 0; k < _lanes; k++)
                            q[k] = Wrap(s._q + (((j * _lanes) + k) * s._c));
                        Store(_q, 0, j, q);
                    }

                    _firstColumn[2 * h] = 0;
                    for (int i = 1; i < h; i++)
                    {
                        var edges = rankToNode[i - 1].InEdges;
                        int penalty = edges.Count == 0 ? s._q - s._c : _negInf;
                        foreach (var it in edges)
                            penalty = Math.Max(penalty, _firstColumn[(2 * h) + _rank[it.Tail.Id] + 1]);
                        _firstColumn[(2 * h) + i] = penalty + s._c;
                    }

                    break;
            }

            if (s.Subtype != SpoaAlignmentSubtype.Linear)
            {
                {
                    for (int j = 0; j < _width; j++)
                    {
                        Store(_f, 0, j, Set1(_negInf));
                        var ev = new int[_lanes];
                        for (int k = 0; k < _lanes; k++)
                            ev[k] = Wrap(s._g + (((j * _lanes) + k) * s._e));
                        Store(_e, 0, j, ev);
                    }

                    _firstColumn[h] = 0;
                    for (int i = 1; i < h; i++)
                    {
                        var edges = rankToNode[i - 1].InEdges;
                        int penalty = edges.Count == 0 ? s._g - s._e : _negInf;
                        foreach (var it in edges)
                            penalty = Math.Max(penalty, _firstColumn[h + _rank[it.Tail.Id] + 1]);
                        _firstColumn[h + i] = penalty + s._e;
                    }
                }
            }

            switch (s.Type)
            {
                case SpoaAlignmentType.Sw:
                    for (int j = 0; j < _width; j++)
                        Store(_h, 0, j, new int[_lanes]);
                    for (int i = 0; i < h; i++)
                        _firstColumn[i] = 0;
                    break;
                case SpoaAlignmentType.Nw:
                    switch (s.Subtype)
                    {
                        case SpoaAlignmentSubtype.Convex:
                            for (int i = 0; i < h; i++)
                                _firstColumn[i] = Math.Max(_firstColumn[h + i], _firstColumn[(2 * h) + i]);
                            for (int j = 0; j < _width; j++)
                                Store(_h, 0, j, Max(Cell(_e, 0, j), Cell(_q, 0, j)));
                            break;
                        case SpoaAlignmentSubtype.Affine:
                            for (int i = 0; i < h; i++)
                                _firstColumn[i] = _firstColumn[h + i];
                            for (int j = 0; j < _width; j++)
                                Store(_h, 0, j, Cell(_e, 0, j));
                            break;
                        default:
                            _firstColumn[0] = 0;
                            for (int i = 1; i < h; i++)
                            {
                                var edges = rankToNode[i - 1].InEdges;
                                int penalty = edges.Count == 0 ? 0 : _negInf;
                                foreach (var it in edges)
                                    penalty = Math.Max(penalty, _firstColumn[_rank[it.Tail.Id] + 1]);
                                _firstColumn[i] = penalty + s._g;
                            }

                            LinearBorder();
                            break;
                    }

                    break;
                default:
                    switch (s.Subtype)
                    {
                        case SpoaAlignmentSubtype.Convex:
                            for (int j = 0; j < _width; j++)
                                Store(_h, 0, j, Max(Cell(_e, 0, j), Cell(_q, 0, j)));
                            break;
                        case SpoaAlignmentSubtype.Affine:
                            for (int j = 0; j < _width; j++)
                                Store(_h, 0, j, Cell(_e, 0, j));
                            break;
                        default:
                            LinearBorder();
                            break;
                    }

                    for (int i = 0; i < h; i++)
                        _firstColumn[i] = 0;
                    break;
            }
        }

        private void LinearBorder()
        {
            for (int j = 0; j < _width; j++)
            {
                var v = new int[_lanes];
                for (int k = 0; k < _lanes; k++)
                    v[k] = Wrap(_s._g + (((j * _lanes) + k) * _s._g));
                Store(_h, 0, j, v);
            }
        }

        private void InitMasks(bool convex)
        {
            int count = convex ? 2 * _logLanes : _logLanes;
            _masks = new int[_logLanes + 1][];
            var unpacked = new int[_lanes];
            for (int i = 0, j = 0; i < _lanes && j < _logLanes; i++)
            {
                unpacked[i] = _negInf;
                if ((i & (i + 1)) == 0)
                    _masks[j++] = (int[])unpacked.Clone();
            }

            var inf = Set1(_negInf);
            _masks[_logLanes] = Slli(inf);
            _ = count;
        }

        // ------------------------------------------------------------------ the three gap models

        public List<(int Node, int Position)> Execute(out int score)
        {
            score = 0;
            Initialize();
            return _s.Subtype switch
            {
                SpoaAlignmentSubtype.Linear => Linear(out score),
                SpoaAlignmentSubtype.Affine => Affine(out score),
                _ => Convex(out score),
            };
        }

        private (int MaxI, int MaxJ, int MaxScore) Finish(int maxScore, int maxI)
        {
            var node = maxI >= 1 ? _graph.RankToNode[maxI - 1] : null;
            int maxJ = -1;
            if (_s.Type == SpoaAlignmentType.Sw)
            {
                maxJ = IndexOf(maxI, maxScore);
            }
            else if (_s.Type == SpoaAlignmentType.Ov)
            {
                maxJ = node!.OutEdges.Count == 0 ? IndexOf(maxI, maxScore) : (int)_normalWidth - 1;
            }
            else
            {
                maxJ = (int)_normalWidth - 1;
            }

            return (maxI, maxJ, maxScore);
        }

        private int IndexOf(int row, int value)
        {
            for (int i = 0; i < _width; i++)
            {
                for (int j = 0; j < _lanes; j++)
                {
                    if (_h[(((row * _width) + i) * _lanes) + j] == value)
                        return (i * _lanes) + j;
                }
            }

            return -1;
        }

        /// <summary><c>_mmxxx_max_value</c>: the maximum of the lanes, starting at 0.</summary>
        private static int MaxValue(int[] v)
        {
            int max = 0;
            foreach (int x in v)
                max = Math.Max(max, x);
            return max;
        }

        private void RowMaximum(SpoaGraph.Node it, int i, int[] score, int[] hLast, ref int maxScore, ref int maxI)
        {
            if (_s.Type == SpoaAlignmentType.Sw)
            {
                int row = MaxValue(score);
                if (maxScore < row)
                {
                    maxScore = row;
                    maxI = i;
                }
            }
            else if (_s.Type == SpoaAlignmentType.Ov)
            {
                if (it.OutEdges.Count == 0)
                {
                    int row = MaxValue(score);
                    if (maxScore < row)
                    {
                        maxScore = row;
                        maxI = i;
                    }
                }
            }
            else if (it.OutEdges.Count == 0)
            {
                int row = hLast[(int)((_normalWidth - 1) % _lanes)];
                if (maxScore < row)
                {
                    maxScore = row;
                    maxI = i;
                }
            }
        }

        private List<(int Node, int Position)> Linear(out int score)
        {
            score = 0;
            var s = _s;
            InitMasks(false);
            var penalties = new int[_logLanes][];
            penalties[0] = Set1(s._g);
            for (int i = 1; i < _logLanes; i++)
                penalties[i] = Add(penalties[i - 1], penalties[i - 1]);
            int maxScore = s.Type == SpoaAlignmentType.Sw ? 0 : _negInf;
            int maxI = -1;
            var zeroes = new int[_lanes];
            var g = Set1(s._g);
            foreach (var it in _graph.RankToNode)
            {
                int i = _rank[it.Id] + 1;
                int predI = it.InEdges.Count == 0 ? 0 : _rank[it.InEdges[0].Tail.Id] + 1;
                var x = SrliLast(Set1(_firstColumn[predI]));
                for (int j = 0; j < _width; j++)
                {
                    var hp = Cell(_h, predI, j);
                    var t1 = SrliLast(hp);
                    var hr = Or(Slli(hp), x);
                    x = t1;
                    Store(_h, i, j, Max(Add(hr, Profile(it, j)), Add(hp, g)));
                }

                for (int p = 1; p < it.InEdges.Count; p++)
                {
                    predI = _rank[it.InEdges[p].Tail.Id] + 1;
                    x = SrliLast(Set1(_firstColumn[predI]));
                    for (int j = 0; j < _width; j++)
                    {
                        var hp = Cell(_h, predI, j);
                        var t1 = SrliLast(hp);
                        var m = Or(Slli(hp), x);
                        x = t1;
                        Store(_h, i, j, Max(Cell(_h, i, j), Max(Add(m, Profile(it, j)), Add(hp, g))));
                    }
                }

                var scoreVec = Set1(_negInf);
                x = SrliLast(Add(Set1(_firstColumn[i]), g));
                int[] last = zeroes;
                for (int j = 0; j < _width; j++)
                {
                    var hr = Max(Cell(_h, i, j), Or(x, _masks[_logLanes]));
                    PrefixMax(hr, penalties, 0);
                    x = SrliLast(Add(hr, g));
                    if (s.Type == SpoaAlignmentType.Sw)
                        hr = Max(hr, zeroes);
                    scoreVec = Max(scoreVec, hr);
                    Store(_h, i, j, hr);
                    last = hr;
                }

                RowMaximum(it, i, scoreVec, last, ref maxScore, ref maxI);
            }

            if (maxI == -1)
                return [];
            score = maxScore;
            var (mi, mj, _) = Finish(maxScore, maxI);
            return Backtrack(mi, mj);
        }

        private int[] Profile(SpoaGraph.Node it, int j) => _profile.AsSpan(((((int)it.Code * _width) + j) * _lanes), _lanes).ToArray();

        private List<(int Node, int Position)> Affine(out int score)
        {
            score = 0;
            var s = _s;
            InitMasks(false);
            var penalties = new int[_logLanes][];
            penalties[0] = Set1(s._e);
            for (int i = 1; i < _logLanes; i++)
                penalties[i] = Add(penalties[i - 1], penalties[i - 1]);
            int maxScore = s.Type == SpoaAlignmentType.Sw ? 0 : _negInf;
            int maxI = -1;
            var zeroes = new int[_lanes];
            var g = Set1(s._g - s._e);
            var e = Set1(s._e);
            foreach (var it in _graph.RankToNode)
            {
                int i = _rank[it.Id] + 1;
                int predI = it.InEdges.Count == 0 ? 0 : _rank[it.InEdges[0].Tail.Id] + 1;
                var x = SrliLast(Set1(_firstColumn[predI]));
                for (int j = 0; j < _width; j++)
                {
                    var hp = Cell(_h, predI, j);
                    Store(_f, i, j, Add(Max(Add(hp, g), Cell(_f, predI, j)), e));
                    Store(_h, i, j, Add(Or(Slli(hp), x), Profile(it, j)));
                    x = SrliLast(hp);
                }

                for (int p = 1; p < it.InEdges.Count; p++)
                {
                    predI = _rank[it.InEdges[p].Tail.Id] + 1;
                    x = SrliLast(Set1(_firstColumn[predI]));
                    for (int j = 0; j < _width; j++)
                    {
                        var hp = Cell(_h, predI, j);
                        Store(_f, i, j, Max(Cell(_f, i, j), Add(Max(Add(hp, g), Cell(_f, predI, j)), e)));
                        Store(_h, i, j, Max(Cell(_h, i, j), Add(Or(Slli(hp), x), Profile(it, j))));
                        x = SrliLast(hp);
                    }
                }

                var scoreVec = zeroes;
                x = Set1(_firstColumn[i]);
                int[] last = zeroes;
                for (int j = 0; j < _width; j++)
                {
                    var hr = Max(Cell(_h, i, j), Cell(_f, i, j));
                    var ev = Add(Add(Or(Slli(hr), SrliLast(x)), g), e);
                    PrefixMax(ev, penalties, 0);
                    Store(_e, i, j, ev);
                    hr = Max(hr, ev);
                    x = Max(hr, Sub(ev, g));
                    if (s.Type == SpoaAlignmentType.Sw)
                        hr = Max(hr, zeroes);
                    scoreVec = Max(scoreVec, hr);
                    Store(_h, i, j, hr);
                    last = hr;
                }

                RowMaximum(it, i, scoreVec, last, ref maxScore, ref maxI);
            }

            if (maxI == -1)
                return [];
            score = maxScore;
            var (mi, mj, _) = Finish(maxScore, maxI);
            return Backtrack(mi, mj);
        }

        private List<(int Node, int Position)> Convex(out int score)
        {
            score = 0;
            var s = _s;
            InitMasks(true);
            var penalties = new int[2 * _logLanes][];
            penalties[0] = Set1(s._e);
            for (int i = 1; i < _logLanes; i++)
                penalties[i] = Add(penalties[i - 1], penalties[i - 1]);
            penalties[_logLanes] = Set1(s._c);
            for (int i = _logLanes + 1; i < 2 * _logLanes; i++)
                penalties[i] = Add(penalties[i - 1], penalties[i - 1]);
            int maxScore = s.Type == SpoaAlignmentType.Sw ? 0 : _negInf;
            int maxI = -1;
            var zeroes = new int[_lanes];
            var g = Set1(s._g - s._e);
            var e = Set1(s._e);
            var q = Set1(s._q - s._c);
            var c = Set1(s._c);
            foreach (var it in _graph.RankToNode)
            {
                int i = _rank[it.Id] + 1;
                int predI = it.InEdges.Count == 0 ? 0 : _rank[it.InEdges[0].Tail.Id] + 1;
                var x = SrliLast(Set1(_firstColumn[predI]));
                for (int j = 0; j < _width; j++)
                {
                    var hp = Cell(_h, predI, j);
                    Store(_f, i, j, Add(Max(Add(hp, g), Cell(_f, predI, j)), e));
                    Store(_o, i, j, Add(Max(Add(hp, q), Cell(_o, predI, j)), c));
                    Store(_h, i, j, Add(Or(Slli(hp), x), Profile(it, j)));
                    x = SrliLast(hp);
                }

                for (int p = 1; p < it.InEdges.Count; p++)
                {
                    predI = _rank[it.InEdges[p].Tail.Id] + 1;
                    x = SrliLast(Set1(_firstColumn[predI]));
                    for (int j = 0; j < _width; j++)
                    {
                        var hp = Cell(_h, predI, j);
                        Store(_f, i, j, Max(Cell(_f, i, j), Add(Max(Add(hp, g), Cell(_f, predI, j)), e)));
                        Store(_o, i, j, Max(Cell(_o, i, j), Add(Max(Add(hp, q), Cell(_o, predI, j)), c)));
                        Store(_h, i, j, Max(Cell(_h, i, j), Add(Or(Slli(hp), x), Profile(it, j))));
                        x = SrliLast(hp);
                    }
                }

                x = Set1(_firstColumn[i]);
                var y = Set1(_firstColumn[i]);
                var scoreVec = zeroes;
                int[] last = zeroes;
                for (int j = 0; j < _width; j++)
                {
                    var hr = Max(Cell(_h, i, j), Max(Cell(_f, i, j), Cell(_o, i, j)));
                    var ev = Add(Add(Or(Slli(hr), SrliLast(x)), g), e);
                    PrefixMax(ev, penalties, 0);
                    Store(_e, i, j, ev);
                    var qv = Add(Add(Or(Slli(hr), SrliLast(y)), q), c);
                    PrefixMax(qv, penalties, _logLanes);
                    Store(_q, i, j, qv);
                    hr = Max(hr, Max(ev, qv));
                    x = Max(hr, Sub(ev, g));
                    y = Max(hr, Sub(qv, q));
                    if (s.Type == SpoaAlignmentType.Sw)
                        hr = Max(hr, zeroes);
                    scoreVec = Max(scoreVec, hr);
                    Store(_h, i, j, hr);
                    last = hr;
                }

                RowMaximum(it, i, scoreVec, last, ref maxScore, ref maxI);
            }

            if (maxI == -1)
                return [];
            score = maxScore;
            var (mi, mj, _) = Finish(maxScore, maxI);
            return Backtrack(mi, mj);
        }

        // ------------------------------------------------------------------ traceback

        /// <summary>The value of a matrix at (row, column) for the traceback; column -1 is the first column of the matrix.</summary>
        private int V(int[] matrix, int row, long col) => col < 0 ? _firstColumn[row] : At(matrix, row, col);

        private List<(int Node, int Position)> Backtrack(int maxI, int maxJ)
        {
            var s = _s;
            var rankToNode = _graph.RankToNode;
            var alignment = new List<(int Node, int Position)>();
            int i = maxI, j = maxJ;
            int prevI = 0, prevJ = 0;
            int h = _height;
            while (true)
            {
                if (j == -1 || i == 0)
                    break;
                var it = rankToNode[i - 1];
                var predecessors = it.InEdges.Count == 0 ? [0] : it.InEdges.Select(x => _rank[x.Tail.Id] + 1).ToList();
                int hij = At(_h, i, j);
                if (s.Type == SpoaAlignmentType.Sw && hij == 0)
                    break;
                int prof = _profile[((int)it.Code * _width * _lanes) + j];
                bool found = false, extendLeft = false, extendUp = false;
                foreach (int p in predecessors)
                {
                    if (hij == V(_h, p, j - 1) + prof)
                    {
                        prevI = p;
                        prevJ = j - 1;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    foreach (int p in predecessors)
                    {
                        if (UpMatches(hij, p, j, out extendUp))
                        {
                            prevI = p;
                            prevJ = j;
                            found = true;
                            break;
                        }
                    }
                }

                if (!found && LeftMatches(hij, i, j, out extendLeft))
                {
                    prevI = i;
                    prevJ = j - 1;
                }

                alignment.Add((i == prevI ? -1 : (int)rankToNode[i - 1].Id, j == prevJ ? -1 : j));
                i = prevI;
                j = prevJ;
                if (extendLeft)
                {
                    while (true)
                    {
                        alignment.Add((-1, j));
                        j--;
                        if (j == -1 || V(_e, i, j) + s._e != V(_e, i, (long)j + 1) || (s.Subtype == SpoaAlignmentSubtype.Convex && V(_q, i, j) + s._c != V(_q, i, (long)j + 1)))
                            break;
                    }
                }
                else if (extendUp)
                {
                    ExtendUp(ref i, j, ref prevI, alignment);
                }
            }

            if (s.Type == SpoaAlignmentType.Nw)
            {
                while (i == 0 && j != -1)
                {
                    alignment.Add((-1, j));
                    j--;
                }

                while (i != 0 && j == -1)
                {
                    alignment.Add(((int)rankToNode[i - 1].Id, -1));
                    var it = rankToNode[i - 1];
                    if (it.InEdges.Count == 0)
                    {
                        i = 0;
                    }
                    else
                    {
                        foreach (var jt in it.InEdges)
                        {
                            int predI = _rank[jt.Tail.Id] + 1;
                            if (FirstColumnStep(i, predI, h))
                            {
                                i = predI;
                                break;
                            }
                        }
                    }
                }
            }

            alignment.Reverse();
            return alignment;
        }

        private bool FirstColumnStep(int i, int predI, int h)
        {
            var s = _s;
            if (s.Subtype == SpoaAlignmentSubtype.Linear)
                return _firstColumn[i] == _firstColumn[predI] + s._g;
            if (s.Subtype == SpoaAlignmentSubtype.Affine)
                return _firstColumn[i] == _firstColumn[predI] + s._e;
            return _firstColumn[h + i] == _firstColumn[h + predI] + s._e || _firstColumn[(2 * h) + i] == _firstColumn[(2 * h) + predI] + s._c;
        }

        private bool UpMatches(int hij, int p, int j, out bool extendUp)
        {
            var s = _s;
            extendUp = false;
            int hp = At(_h, p, j);
            if (s.Subtype == SpoaAlignmentSubtype.Linear)
                return hij == hp + s._g;
            if (s.Subtype == SpoaAlignmentSubtype.Affine)
            {
                if (hij == At(_f, p, j) + s._e)
                {
                    extendUp = true;
                    return true;
                }

                return hij == hp + s._g;
            }

            if (hij == At(_f, p, j) + s._e)
            {
                extendUp = true;
                return true;
            }

            if (hij == hp + s._g)
                return true;
            if (hij == At(_o, p, j) + s._c)
            {
                extendUp = true;
                return true;
            }

            return hij == hp + s._q;
        }

        private bool LeftMatches(int hij, int i, int j, out bool extendLeft)
        {
            var s = _s;
            extendLeft = false;
            int hl = V(_h, i, j - 1);
            if (s.Subtype == SpoaAlignmentSubtype.Linear)
                return hij == hl + s._g;
            if (hij == EAt(i, j - 1) + s._e)
            {
                extendLeft = true;
                return true;
            }

            if (hij == hl + s._g)
                return true;
            if (s.Subtype == SpoaAlignmentSubtype.Convex)
            {
                if (hij == QAt(i, j - 1) + s._c)
                {
                    extendLeft = true;
                    return true;
                }

                return hij == hl + s._q;
            }

            return false;
        }

        private int EAt(int row, long col) => V(_e, row, col);

        private int QAt(int row, long col) => V(_q, row, col);

        /// <summary>The <c>extend_up</c> loop: moves up through the F / O gaps, emitting (node, -1).</summary>
        private void ExtendUp(ref int i, int j, ref int outerPrevI, List<(int Node, int Position)> alignment)
        {
            var s = _s;
            var rankToNode = _graph.RankToNode;
            while (true)
            {
                var preds = rankToNode[i - 1].InEdges.Select(x => _rank[x.Tail.Id] + 1).ToList();
                int prevI = 0;
                bool stop;
                if (s.Subtype == SpoaAlignmentSubtype.Affine)
                {
                    stop = false;
                    foreach (int p in preds)
                    {
                        stop = At(_f, i, j) == At(_h, p, j) + s._g;
                        if (stop || At(_f, i, j) == At(_f, p, j) + s._e)
                        {
                            prevI = p;
                            break;
                        }
                    }
                }
                else
                {
                    stop = true;
                    foreach (int p in preds)
                    {
                        if (At(_f, i, j) == At(_f, p, j) + s._e || At(_o, i, j) == At(_o, p, j) + s._c)
                        {
                            prevI = p;
                            stop = false;
                            break;
                        }
                    }

                    if (stop)
                    {
                        foreach (int p in preds)
                        {
                            if (At(_f, i, j) == At(_h, p, j) + s._g || At(_o, i, j) == At(_h, p, j) + s._q)
                            {
                                prevI = p;
                                break;
                            }
                        }
                    }
                }

                alignment.Add(((int)rankToNode[i - 1].Id, -1));
                i = prevI;
                outerPrevI = prevI;
                if (stop || i == 0)
                    break;
            }
        }
    }
}
