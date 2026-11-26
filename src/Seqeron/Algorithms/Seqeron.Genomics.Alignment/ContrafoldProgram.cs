#pragma warning disable S3776, S107, S138, S1541, S1854, S3358, CA1814, S2368, S125, S1450, S2234
using System.Globalization;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// CONTRAfold 2.02 (Do, Woods, Batzoglou 2006, Bioinformatics 22:e90-e98) as <c>contrafold predict SEQUENCE.fa --posteriors CUTOFF OUT</c> runs it with the default parameters of the program (<c>contrafold.params.complementary</c>,
/// <c>RunPredictionMode&lt;float&gt;</c>: the whole inference is carried out in single precision with the polynomial approximations <c>Fast_Exp</c> and <c>Fast_LogExpPlusOne</c> of <c>LogSpace.hpp</c>):
/// the inside and outside algorithms of the conditional log-linear model (<c>ComputeInside</c>, <c>ComputeOutside</c>) and the posterior base-pairing probabilities (<c>ComputePosterior</c>) in the sparse BPSEQ text (<c>PrintSparseBPSEQ</c>).
/// This is the program that <c>mafft --contrafold</c> (<c>contrafoldwrap</c>) runs for every sequence.
/// <para>Licence: CONTRAfold is distributed under the BSD 3-clause licence (Copyright (c) 2006, Chuong Do; manual.tex "License (BSD)"): the notice is kept in <c>THIRD-PARTY-NOTICES.md</c>; this port is part of the MIT assembly.
/// Source read: contrafold_v2_02.tar.gz (<c>InferenceEngine.ipp</c>, <c>InferenceEngine.hpp</c>, <c>LogSpace.hpp</c>, <c>Config.hpp</c>, <c>ComputationEngine.ipp</c>, <c>SparseMatrix.ipp</c>, <c>SStruct.cpp</c>, <c>ParameterManager.ipp</c>, <c>MakeDefaults.pl</c>).</para>
/// </summary>
public static partial class ContrafoldProgram
{
    private const int Alphabet = 4; // ACGU; every other letter is the "unknown" index 4
    private const int MaxHairpin = 30, MaxBulge = 30, MaxInternal = 30, MaxInternalSymmetric = 15, MaxInternalAsymmetry = 28, MaxInternalExplicit = 4, MaxSingle = 30;
    private const float NegInf = (float)-2e20;
    private const float NegInfHalf = (float)(-2e20 / 2);
    private const float LogAddLimit = (float)11.8624794162;

    private static readonly Lazy<Dictionary<string, float>> Defaults = new(() =>
    {
        var d = new Dictionary<string, float>(StringComparer.Ordinal);
        foreach (var line in DefaultComplementaryValues.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            int sp = line.IndexOf(' ');
            d[line[..sp]] = (float)double.Parse(line[(sp + 1)..], CultureInfo.InvariantCulture); // the generated table holds double literals that the program converts to float
        }
        return d;
    });

    /// <summary>
    /// <c>contrafold predict file --posteriors cutoff out</c> for the single sequence of the file: the text of the output file, one line per nucleotide, <c>"i letter j:p j:p ..."</c> (the pairs (i, j), j &gt; i, whose posterior probability is at least
    /// <paramref name="cutoff"/>; the probabilities in the default format of an output stream: six significant digits).
    /// </summary>
    /// <param name="sequence">The sequence (letters; A, C, G, U / T as in the program: other letters, including T, are the unknown symbol that pairs with nothing).</param>
    /// <param name="cutoff">The posterior cutoff (<c>--posteriors 0.01</c> in mafft).</param>
    /// <exception cref="ArgumentException">The sequence is empty ("Zero-length sequence read." of the program).</exception>
    public static string PredictPosteriors(string sequence, double cutoff = 0.01)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0)
            throw new ArgumentException("Zero-length sequence read.", nameof(sequence));
        var engine = new Engine(sequence);
        engine.ComputeInside();
        engine.ComputeOutside();
        engine.ComputePosterior();
        return engine.Sparse((float)cutoff, sequence);
    }

    private static float FastExp(float x)
    {
        if (x < (float)-2.4915033807)
        {
            if (x < (float)-5.8622823336)
            {
                if (x < (float)-9.91152)
                    return 0f;
                return (((float)0.0000803850 * x + (float)0.0021627428) * x + (float)0.0194708555) * x + (float)0.0588080014;
            }
            if (x < (float)-3.8396630909)
                return (((float)0.0013889414 * x + (float)0.0244676474) * x + (float)0.1471290604) * x + (float)0.3042757740;
            return (((float)0.0072335607 * x + (float)0.0906002677) * x + (float)0.3983111356) * x + (float)0.6245959221;
        }
        if (x < (float)-0.6725053211)
        {
            if (x < (float)-1.4805375919)
                return (((float)0.0232410351 * x + (float)0.2085645908) * x + (float)0.6906367911) * x + (float)0.8682322329;
            return (((float)0.0573782771 * x + (float)0.3580258429) * x + (float)0.9121133217) * x + (float)0.9793091728;
        }
        if (x < 0f)
            return (((float)0.1199175927 * x + (float)0.4815668234) * x + (float)0.9975991939) * x + (float)0.9999505077;
        return x > (float)46.052 ? (float)1e20 : MathF.Exp(x);
    }

    private static float FastLogExpPlusOne(float x)
    {
        if (x < (float)3.3792499610)
        {
            if (x < (float)1.6320158198)
            {
                if (x < (float)0.6615367791)
                    return (((float)-0.0065591595 * x + (float)0.1276442762) * x + (float)0.4996554598) * x + (float)0.6931542306;
                return (((float)-0.0155157557 * x + (float)0.1446775699) * x + (float)0.4882939746) * x + (float)0.6958092989;
            }
            if (x < (float)2.4912588184)
                return (((float)-0.0128909247 * x + (float)0.1301028251) * x + (float)0.5150398748) * x + (float)0.6795585882;
            return (((float)-0.0072142647 * x + (float)0.0877540853) * x + (float)0.6208708362) * x + (float)0.5909675829;
        }
        if (x < (float)5.7890710412)
        {
            if (x < (float)4.4261691294)
                return (((float)-0.0031455354 * x + (float)0.0467229449) * x + (float)0.7592532310) * x + (float)0.4348794399;
            return (((float)-0.0010110698 * x + (float)0.0185943421) * x + (float)0.8831730747) * x + (float)0.2523695427;
        }
        if (x < (float)7.8162726752)
            return (((float)-0.0001962780 * x + (float)0.0046084408) * x + (float)0.9634431978) * x + (float)0.0983148903;
        return (((float)-0.0000113994 * x + (float)0.0003734731) * x + (float)0.9959107193) * x + (float)0.0149855051;
    }

    /// <summary><c>Fast_LogPlusEquals(float &amp;x, float y)</c>.</summary>
    private static void LogPlusEquals(ref float x, float y)
    {
        if (x < y)
            (x, y) = (y, x);
        if (y > NegInfHalf && x - y < LogAddLimit)
            x = FastLogExpPlusOne(x - y) + y;
    }

    internal sealed class Engine
    {
        private readonly int _l, _size;
        private readonly int[] _s, _offset;
        private readonly bool[] _allowPaired;
        // parameters (float, [..] indices as in the program; the index 4 = unknown symbol = 0)
        private readonly float[,] _basePair = new float[5, 5];
        private readonly float[,,,] _terminalMismatch = new float[5, 5, 5, 5];
        private readonly float[] _hairpinCache = new float[MaxHairpin + 1];
        private readonly float[,] _singleCache = new float[MaxSingle + 1, MaxSingle + 1];
        private readonly float[] _bulge0x1 = new float[5], _bulge1x0 = new float[5];
        private readonly float[,] _internal1x1 = new float[5, 5];
        private readonly float[,,,] _helixStacking = new float[5, 5, 5, 5];
        private readonly float[,] _helixClosing = new float[5, 5];
        private readonly float[,,] _dangleLeft = new float[5, 5, 5], _dangleRight = new float[5, 5, 5];
        private readonly float _multiBase, _multiUnpaired, _multiPaired, _externalUnpaired, _externalPaired;
        private float[] _f5i = null!, _fci = null!, _fmi = null!, _fm1i = null!, _f5o = null!, _fco = null!, _fmo = null!, _fm1o = null!, _posterior = null!;

        /// <summary>The log partition coefficient from the inside pass (<c>ComputeLogPartitionCoefficient</c>: <c>F5i[L]</c>).</summary>
        internal float InsideLogPartition => _f5i[_l];

        /// <summary>The same quantity from the outside pass (<c>F5o[0]</c>): the two agree when the recursions are consistent.</summary>
        internal float OutsideLogPartition => _f5o[0];

        /// <summary>The posterior base-pairing probability of the positions i &lt; j (1-based), after <see cref="ComputePosterior"/>.</summary>
        internal float Posterior(int i, int j) => _posterior[_offset[i] + j];

        internal Engine(string sequence)
        {
            _l = sequence.Length;
            _size = (_l + 1) * (_l + 2) / 2;
            _s = new int[_l + 1];
            _s[0] = Alphabet;
            for (int i = 1; i <= _l; i++)
                _s[i] = SymbolOf(sequence[i - 1]);
            _offset = new int[_l + 1];
            for (int i = 0; i <= _l; i++)
                _offset[i] = i * (_l + 1 + _l + 1 - i - 1) / 2;
            _allowPaired = new bool[_size];
            for (int i = 1; i <= _l; i++)
                for (int j = i + 1; j <= _l; j++)
                    _allowPaired[_offset[i] + j] = IsComplementary(_s[i], _s[j]);

            var v = Defaults.Value;
            const string A = "ACGU";
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    _basePair[i, j] = v[Symmetric("base_pair_", A[i].ToString() + A[j], A[j].ToString() + A[i])];
            for (int i1 = 0; i1 < 4; i1++)
                for (int j1 = 0; j1 < 4; j1++)
                    for (int i2 = 0; i2 < 4; i2++)
                        for (int j2 = 0; j2 < 4; j2++)
                            _terminalMismatch[i1, j1, i2, j2] = v["terminal_mismatch_" + A[i1] + A[j1] + A[i2] + A[j2]];
            var hairpin = new float[MaxHairpin + 1];
            for (int i = 0; i <= MaxHairpin; i++)
                hairpin[i] = v["hairpin_length_at_least_" + i.ToString(CultureInfo.InvariantCulture)];
            var explicitInternal = new float[MaxInternalExplicit + 1, MaxInternalExplicit + 1];
            for (int i = 1; i <= MaxInternalExplicit; i++)
                for (int j = 1; j <= MaxInternalExplicit; j++)
                    explicitInternal[i, j] = v["internal_explicit_" + Math.Min(i, j).ToString(CultureInfo.InvariantCulture) + "_" + Math.Max(i, j).ToString(CultureInfo.InvariantCulture)];
            var bulge = new float[MaxBulge + 1];
            for (int i = 1; i <= MaxBulge; i++)
                bulge[i] = v["bulge_length_at_least_" + i.ToString(CultureInfo.InvariantCulture)];
            var internalLength = new float[MaxInternal + 1];
            for (int i = 2; i <= MaxInternal; i++)
                internalLength[i] = v["internal_length_at_least_" + i.ToString(CultureInfo.InvariantCulture)];
            var symmetric = new float[MaxInternalSymmetric + 1];
            for (int i = 1; i <= MaxInternalSymmetric; i++)
                symmetric[i] = v["internal_symmetric_length_at_least_" + i.ToString(CultureInfo.InvariantCulture)];
            var asymmetry = new float[MaxInternalAsymmetry + 1];
            for (int i = 1; i <= MaxInternalAsymmetry; i++)
                asymmetry[i] = v["internal_asymmetry_at_least_" + i.ToString(CultureInfo.InvariantCulture)];
            for (int i = 0; i < 4; i++)
            {
                _bulge0x1[i] = v["bulge_0x1_nucleotides_" + A[i]];
                _bulge1x0[i] = _bulge0x1[i];
            }
            for (int i1 = 0; i1 < 4; i1++)
                for (int i2 = 0; i2 < 4; i2++)
                    _internal1x1[i1, i2] = v[Symmetric("internal_1x1_nucleotides_", A[i1].ToString() + A[i2], A[i2].ToString() + A[i1])];
            for (int i1 = 0; i1 < 4; i1++)
                for (int j1 = 0; j1 < 4; j1++)
                    for (int i2 = 0; i2 < 4; i2++)
                        for (int j2 = 0; j2 < 4; j2++)
                            _helixStacking[i1, j1, i2, j2] = v[Symmetric("helix_stacking_", A[i1].ToString() + A[j1] + A[i2] + A[j2], A[j2].ToString() + A[i2] + A[j1] + A[i1])];
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    _helixClosing[i, j] = v["helix_closing_" + A[i] + A[j]];
            _multiBase = v["multi_base"];
            _multiUnpaired = v["multi_unpaired"];
            _multiPaired = v["multi_paired"];
            for (int i1 = 0; i1 < 4; i1++)
                for (int j1 = 0; j1 < 4; j1++)
                    for (int k = 0; k < 4; k++)
                    {
                        _dangleLeft[i1, j1, k] = v["dangle_left_" + A[i1] + A[j1] + A[k]];
                        _dangleRight[i1, j1, k] = v["dangle_right_" + A[i1] + A[j1] + A[k]];
                    }
            _externalUnpaired = v["external_unpaired"];
            _externalPaired = v["external_paired"];

            // InitializeCache: the cumulative length scores
            _hairpinCache[0] = hairpin[0];
            for (int i = 1; i <= MaxHairpin; i++)
                _hairpinCache[i] = _hairpinCache[i - 1] + hairpin[i];
            var tb = Cumulative(bulge);
            var ti = Cumulative(internalLength);
            var ts = Cumulative(symmetric);
            var ta = Cumulative(asymmetry);
            for (int l1 = 0; l1 <= MaxSingle; l1++)
            {
                for (int l2 = 0; l1 + l2 <= MaxSingle; l2++)
                {
                    float c = 0f;
                    if (l1 != 0 || l2 != 0)
                    {
                        if (l1 == 0 || l2 == 0)
                        {
                            c += tb[Math.Min(MaxBulge, l1 + l2)];
                        }
                        else
                        {
                            if (l1 <= MaxInternalExplicit && l2 <= MaxInternalExplicit)
                                c += explicitInternal[l1, l2];
                            c += ti[Math.Min(MaxInternal, l1 + l2)];
                            if (l1 == l2)
                                c += ts[Math.Min(MaxInternalSymmetric, l1)];
                            c += ta[Math.Min(MaxInternalAsymmetry, Math.Abs(l1 - l2))];
                        }
                    }
                    _singleCache[l1, l2] = c;
                }
            }
        }

        private static float[] Cumulative(float[] at)
        {
            var r = new float[at.Length];
            r[0] = at[0];
            for (int i = 1; i < at.Length; i++)
                r[i] = r[i - 1] + at[i];
            return r;
        }

        private static string Symmetric(string prefix, string a, string b) => string.CompareOrdinal(prefix + a, prefix + b) < 0 ? prefix + a : prefix + b;

        private static int SymbolOf(char c) => char.ToUpperInvariant(c) switch { 'A' => 0, 'C' => 1, 'G' => 2, 'U' => 3, _ => Alphabet };

        private static bool IsComplementary(int a, int b) => (a, b) is (0, 3) or (3, 0) or (2, 3) or (3, 2) or (1, 2) or (2, 1);

        // ------------------------------------------------------------------ the scores of the inference engine (ScoreJunctionA ... ScoreSingleNucleotides)

        private float ScoreJunctionA(int i, int j) =>
            0f
            + _helixClosing[_s[i], _s[j + 1]]
            + (i < _l ? _dangleLeft[_s[i], _s[j + 1], _s[i + 1]] : 0f)
            + (j > 0 ? _dangleRight[_s[i], _s[j + 1], _s[j]] : 0f);

        private float ScoreJunctionB(int i, int j) => 0f + _helixClosing[_s[i], _s[j + 1]] + _terminalMismatch[_s[i], _s[j + 1], _s[i + 1], _s[j]];

        private float ScoreBasePair(int i, int j) => 0f + _basePair[_s[i], _s[j]];

        private float ScoreHairpin(int i, int j) => 0f + ScoreJunctionB(i, j) + _hairpinCache[Math.Min(j - i, MaxHairpin)];

        private float ScoreHelixStacking(int i, int j) => _helixStacking[_s[i], _s[j], _s[i + 1], _s[j - 1]];

        private float ScoreSingleNucleotides(int i, int j, int p, int q)
        {
            int l1 = p - i, l2 = j - q;
            return 0f + 0f + 0f
                + (l1 == 0 && l2 == 1 ? _bulge0x1[_s[j]] : 0f)
                + (l1 == 1 && l2 == 0 ? _bulge1x0[_s[i + 1]] : 0f)
                + (l1 == 1 && l2 == 1 ? _internal1x1[_s[i + 1], _s[j]] : 0f);
        }

        private float MultiUnpaired => _multiUnpaired + 0f;

        private float ExternalUnpaired => _externalUnpaired + 0f;

        // ------------------------------------------------------------------ ComputeInside

        internal void ComputeInside()
        {
            int L = _l;
            var off = _offset;
            _f5i = Fill(L + 1);
            _fci = Fill(_size);
            _fmi = Fill(_size);
            _fm1i = Fill(_size);
            for (int i = L; i >= 0; i--)
            {
                for (int j = i; j <= L; j++)
                {
                    float fm2i = NegInf;
                    if (i + 2 <= j)
                    {
                        for (int k = i + 1; k < j; k++)
                            LogPlusEquals(ref fm2i, _fm1i[off[i] + k] + _fmi[off[k] + j]);
                    }
                    if (0 < i && j < L && _allowPaired[off[i] + j + 1])
                    {
                        float sumI = NegInf;
                        LogPlusEquals(ref sumI, ScoreHairpin(i, j));
                        {
                            float scoreHelix = i + 2 <= j ? ScoreBasePair(i + 1, j) + ScoreHelixStacking(i, j + 1) : 0f;
                            float scoreOther = ScoreJunctionB(i, j);
                            for (int p = i; p <= Math.Min(i + MaxSingle, j); p++)
                            {
                                int qMin = Math.Max(p + 2, p - i + j - MaxSingle);
                                int baseIdx = off[p + 1] - 1;
                                for (int q = j; q >= qMin; q--)
                                {
                                    if (!_allowPaired[off[p + 1] + q])
                                        continue;
                                    float score = p == i && q == j
                                        ? scoreHelix + _fci[baseIdx + q]
                                        : scoreOther + _singleCache[p - i, j - q] + _fci[baseIdx + q] + ScoreBasePair(p + 1, q) + ScoreJunctionB(q, p) + ScoreSingleNucleotides(i, j, p, q);
                                    LogPlusEquals(ref sumI, score);
                                }
                            }
                        }
                        LogPlusEquals(ref sumI, fm2i + ScoreJunctionA(i, j) + _multiPaired + _multiBase);
                        _fci[off[i] + j] = sumI;
                    }
                    if (0 < i && i + 2 <= j && j < L)
                    {
                        float sumI = NegInf;
                        if (_allowPaired[off[i + 1] + j])
                            LogPlusEquals(ref sumI, _fci[off[i + 1] + j - 1] + ScoreJunctionA(j, i) + _multiPaired + ScoreBasePair(i + 1, j));
                        LogPlusEquals(ref sumI, _fm1i[off[i + 1] + j] + MultiUnpaired);
                        _fm1i[off[i] + j] = sumI;
                    }
                    if (0 < i && i + 2 <= j && j < L)
                    {
                        float sumI = NegInf;
                        LogPlusEquals(ref sumI, fm2i);
                        LogPlusEquals(ref sumI, _fmi[off[i] + j - 1] + MultiUnpaired);
                        LogPlusEquals(ref sumI, _fm1i[off[i] + j]);
                        _fmi[off[i] + j] = sumI;
                    }
                }
            }
            _f5i[0] = 0f;
            for (int j = 1; j <= L; j++)
            {
                float sumI = NegInf;
                LogPlusEquals(ref sumI, _f5i[j - 1] + ExternalUnpaired);
                for (int k = 0; k < j; k++)
                    if (_allowPaired[off[k + 1] + j])
                        LogPlusEquals(ref sumI, _f5i[k] + _fci[off[k + 1] + j - 1] + _externalPaired + ScoreBasePair(k + 1, j) + ScoreJunctionA(j, k));
                _f5i[j] = sumI;
            }
        }

        private static float[] Fill(int n)
        {
            var a = new float[n];
            Array.Fill(a, NegInf);
            return a;
        }

        // ------------------------------------------------------------------ ComputeOutside

        internal void ComputeOutside()
        {
            int L = _l;
            var off = _offset;
            _f5o = Fill(L + 1);
            _fco = Fill(_size);
            _fmo = Fill(_size);
            _fm1o = Fill(_size);
            _f5o[L] = 0f;
            for (int j = L; j >= 1; j--)
            {
                LogPlusEquals(ref _f5o[j - 1], _f5o[j] + ExternalUnpaired);
                for (int k = 0; k < j; k++)
                {
                    if (_allowPaired[off[k + 1] + j])
                    {
                        float temp = _f5o[j] + _externalPaired + ScoreBasePair(k + 1, j) + ScoreJunctionA(j, k);
                        LogPlusEquals(ref _f5o[k], temp + _fci[off[k + 1] + j - 1]);
                        LogPlusEquals(ref _fco[off[k + 1] + j - 1], temp + _f5i[k]);
                    }
                }
            }
            for (int i = 0; i <= L; i++)
            {
                for (int j = L; j >= i; j--)
                {
                    float fm2o = NegInf;
                    if (0 < i && i + 2 <= j && j < L)
                    {
                        LogPlusEquals(ref fm2o, _fmo[off[i] + j]);
                        LogPlusEquals(ref _fmo[off[i] + j - 1], _fmo[off[i] + j] + MultiUnpaired);
                        LogPlusEquals(ref _fm1o[off[i] + j], _fmo[off[i] + j]);
                    }
                    if (0 < i && i + 2 <= j && j < L)
                    {
                        if (_allowPaired[off[i + 1] + j])
                            LogPlusEquals(ref _fco[off[i + 1] + j - 1], _fm1o[off[i] + j] + ScoreJunctionA(j, i) + _multiPaired + ScoreBasePair(i + 1, j));
                        LogPlusEquals(ref _fm1o[off[i + 1] + j], _fm1o[off[i] + j] + MultiUnpaired);
                    }
                    if (0 < i && j < L && _allowPaired[off[i] + j + 1])
                    {
                        {
                            float scoreHelix = i + 2 <= j ? _fco[off[i] + j] + ScoreBasePair(i + 1, j) + ScoreHelixStacking(i, j + 1) : 0f;
                            float scoreOther = _fco[off[i] + j] + ScoreJunctionB(i, j);
                            for (int p = i; p <= Math.Min(i + MaxSingle, j); p++)
                            {
                                int qMin = Math.Max(p + 2, p - i + j - MaxSingle);
                                int baseIdx = off[p + 1] - 1;
                                for (int q = j; q >= qMin; q--)
                                {
                                    if (!_allowPaired[off[p + 1] + q])
                                        continue;
                                    float value = p == i && q == j
                                        ? scoreHelix
                                        : scoreOther + _singleCache[p - i, j - q] + ScoreBasePair(p + 1, q) + ScoreJunctionB(q, p) + ScoreSingleNucleotides(i, j, p, q);
                                    LogPlusEquals(ref _fco[baseIdx + q], value);
                                }
                            }
                        }
                        LogPlusEquals(ref fm2o, _fco[off[i] + j] + ScoreJunctionA(i, j) + _multiPaired + _multiBase);
                    }
                    if (i + 2 <= j)
                    {
                        for (int k = i + 1; k < j; k++)
                        {
                            LogPlusEquals(ref _fm1o[off[i] + k], fm2o + _fmi[off[k] + j]);
                            LogPlusEquals(ref _fmo[off[k] + j], fm2o + _fm1i[off[i] + k]);
                        }
                    }
                }
            }
        }

        // ------------------------------------------------------------------ ComputePosterior

        internal void ComputePosterior()
        {
            int L = _l;
            var off = _offset;
            _posterior = new float[_size];
            float Z = _f5i[L];
            for (int i = L; i >= 0; i--)
            {
                for (int j = i; j <= L; j++)
                {
                    if (0 < i && j < L && _allowPaired[off[i] + j + 1])
                    {
                        float outside = _fco[off[i] + j] - Z;
                        float scoreHelix = i + 2 <= j ? outside + ScoreBasePair(i + 1, j) + ScoreHelixStacking(i, j + 1) : 0f;
                        float scoreOther = outside + ScoreJunctionB(i, j);
                        for (int p = i; p <= Math.Min(i + MaxSingle, j); p++)
                        {
                            int qMin = Math.Max(p + 2, p - i + j - MaxSingle);
                            int baseIdx = off[p + 1] - 1;
                            for (int q = j; q >= qMin; q--)
                            {
                                if (!_allowPaired[off[p + 1] + q])
                                    continue;
                                _posterior[off[p + 1] + q] += FastExp(p == i && q == j
                                    ? scoreHelix + _fci[baseIdx + q]
                                    : scoreOther + _singleCache[p - i, j - q] + _fci[baseIdx + q] + ScoreBasePair(p + 1, q) + ScoreJunctionB(q, p) + ScoreSingleNucleotides(i, j, p, q));
                            }
                        }
                    }
                    if (0 < i && i + 2 <= j && j < L && _allowPaired[off[i + 1] + j])
                        _posterior[off[i + 1] + j] += FastExp(_fm1o[off[i] + j] + _fci[off[i + 1] + j - 1] + ScoreJunctionA(j, i) + _multiPaired + ScoreBasePair(i + 1, j) - Z);
                }
            }
            for (int j = 1; j <= L; j++)
            {
                float outside = _f5o[j] - Z;
                for (int k = 0; k < j; k++)
                    if (_allowPaired[off[k + 1] + j])
                        _posterior[off[k + 1] + j] += FastExp(outside + _f5i[k] + _fci[off[k + 1] + j - 1] + _externalPaired + ScoreBasePair(k + 1, j) + ScoreJunctionA(j, k));
            }
            for (int i = 1; i <= L; i++)
                for (int j = i + 1; j <= L; j++)
                    _posterior[off[i] + j] = Math.Clamp(_posterior[off[i] + j], 0f, 1f);
        }

        /// <summary><c>GetPosterior(cutoff)</c> and <c>SparseMatrix::PrintSparseBPSEQ</c>: the pairs at least as probable as the cutoff, the values in the default format of an output stream (six significant digits).</summary>
        internal string Sparse(float cutoff, string sequence)
        {
            var sb = new System.Text.StringBuilder();
            for (int i = 1; i <= _l; i++)
            {
                sb.Append(i.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(sequence[i - 1]);
                for (int j = i; j <= _l; j++)
                {
                    float v = _posterior[_offset[i] + j] >= cutoff ? _posterior[_offset[i] + j] : 0f;
                    if (v != 0f)
                        sb.Append(' ').Append(j.ToString(CultureInfo.InvariantCulture)).Append(':').Append(((double)v).ToString("G6", CultureInfo.InvariantCulture));
                }
                sb.Append('\n');
            }
            return sb.ToString();
        }
    }
}
