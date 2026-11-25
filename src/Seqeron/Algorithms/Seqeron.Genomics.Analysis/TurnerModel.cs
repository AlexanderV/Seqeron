using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The Turner 2004 nearest-neighbour energy model of <see cref="RnaSecondaryStructure"/> as one immutable,
/// array-backed parameter set at a given temperature — the single place where every loop term is evaluated
/// (MFE, traceback, structure evaluation, partition function and the public single-term functions all call it).
/// </summary>
/// <remarks>
/// <para>
/// <b>Temperature.</b> The free-energy tables are the NNDB 37 °C values (<see cref="Turner2004Parameters"/>,
/// <see cref="Turner2004InteriorLoopTables"/>). At another temperature every tabulated entry is rescaled with
/// its Turner 2004 enthalpy exactly as ViennaRNA does (<c>params.c</c> <c>RESCALE_dG</c>):
/// <c>ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15 K</c>, the ΔH being ViennaRNA 2.7.2's
/// (<see cref="Turner2004Enthalpies"/>); the logarithmic loop-size extrapolations scale their
/// 1.75·R·T / 1.08·T/T₃₇ prefactor (<c>lxc = lxc37·T/310.15</c>); the bulge degeneracy term is −RT·ln(states);
/// NNDB-only terms without a published enthalpy (special G·U closure −2.2, all-C hairpins, special-C bulge −0.9,
/// 5'GGUC/3'CUGG three-stack −4.12) are held at their 37 °C value. At 37 °C the model is bit-for-bit the NNDB
/// table model (per-loop results rounded to 0.01 kcal/mol as tabulated); at other temperatures the rescaled
/// values are used unrounded (ViennaRNA's partition-function convention, <c>pf_smooth</c>).
/// </para>
/// <para>
/// Index conventions are ViennaRNA's: pair type 1..7 = CG GC GU UG AU UA non-standard
/// (<see cref="Turner2004InteriorLoopTables.PairType"/>), base code 0..4 = N A C G U. A nucleotide that is not
/// A/C/G/U/T is the unknown base N whose mismatch / dangle entries are the maximum over the four bases.
/// </para>
/// </remarks>
internal sealed class TurnerModel
{
    /// <summary>"Impossible" energy (a pair that cannot be formed).</summary>
    public const double Inf = 1e18;

    /// <summary>ViennaRNA <c>GASCONST</c> in cal/(mol·K) (<c>params/constants.h</c>).</summary>
    public const double GasConstantCalPerMolK = 1.98717;

    /// <summary>Temperature of the NNDB / ViennaRNA parameter measurement, 37 °C = 310.15 K (<c>Tmeasure</c>).</summary>
    public const double ReferenceTemperatureKelvin = 310.15;

    // NNDB hairpin first-mismatch bonuses (hairpin-mismatch-parameters), additive on the terminal mismatch.
    private const double UuGaMismatchBonus = -0.9;
    private const double GgMismatchBonus = -0.8;
    // NNDB-only constants (no ViennaRNA counterpart, no published enthalpy → held at 37 °C).
    private const double SpecialGuClosureBonus = -2.2;
    private const double AllCPenaltyPerNt = 0.3;
    private const double AllCPenaltyConstantGt3 = 1.6;
    private const double AllCPenalty3Nt = 1.5;
    private const double BulgeSpecialCBonus = -0.9;

    /// <summary>Multibranch closing term a (ViennaRNA <c>ML_closing37</c> 930 dcal), kcal/mol at 37 °C.</summary>
    public const double MlClosing37 = 9.3;

    /// <summary>Multibranch per-branch term c (<c>ML_intern37</c> −90 dcal), kcal/mol at 37 °C.</summary>
    public const double MlIntern37 = -0.9;

    /// <summary>Multibranch per-unpaired term b (<c>ML_BASE37</c> 0), kcal/mol at 37 °C.</summary>
    public const double MlBase37 = 0.0;

    /// <summary>Special 5'GGUC/3'CUGG three-stack total (NNDB GU parameters, note b).</summary>
    public const double SpecialGgucCugg3Stack = -4.12;

    // NNDB internal-loop log extrapolation 1.08·ln(n/6) (n > 30) and hairpin / bulge 1.75·R·T·ln(n/9 | n/6),
    // both at 37 °C as published (R = 1.987), scaled by T/310.15 like ViennaRNA's lxc.
    private const double InternalLoopLxc37 = 1.08;
    private const double HairpinBulgeLxc37 = 1.75 * 1.987 * ReferenceTemperatureKelvin / 1000.0;
    private const double MaxNinio = 3.0;

    private const int MaxTabulatedLoop = 30;

    private static readonly char[] BaseLetters = { 'N', 'A', 'C', 'G', 'U' };
    private static readonly byte[] PairTypeTable = BuildPairTypeTable();
    private static readonly (char Five, char Three)[] TypeLetters =
        { (' ', ' '), ('C', 'G'), ('G', 'C'), ('G', 'U'), ('U', 'G'), ('A', 'U'), ('U', 'A'), (' ', ' ') };

    private static readonly object CacheLock = new();
    private static readonly Dictionary<double, TurnerModel> Cache = new();

    /// <summary>The 37 °C NNDB model.</summary>
    public static TurnerModel At37 { get; } = new(ReferenceTemperatureKelvin);

    /// <summary>Temperature in Kelvin.</summary>
    public double Temperature { get; }

    /// <summary>RT in kcal/mol (ViennaRNA <c>kT</c>).</summary>
    public double RT { get; }

    /// <summary>T / 310.15 K (ViennaRNA <c>tempf</c>).</summary>
    public double TempFactor { get; }

    /// <summary>True at exactly 37 °C: loop terms are rounded to 0.01 kcal/mol like the NNDB tables.</summary>
    public bool IsReference { get; }

    /// <summary>Terminal AU/GU penalty (per helix end closed by A·U or G·U).</summary>
    public double TerminalAu { get; }

    /// <summary>Multibranch loop closing term a (ML_closing).</summary>
    public double MlClosing { get; }

    /// <summary>Multibranch per-branch term c (ML_intern).</summary>
    public double MlIntern { get; }

    /// <summary>Multibranch per-unpaired term b (ML_BASE), 0 in the Turner 2004 set.</summary>
    public double MlBase { get; }

    private readonly double ninio;
    private readonly double lxcHairpinBulge;
    private readonly double lxcInternal;
    private readonly double[] stack = new double[64];          // [type*8 + type2], both read 5'→3'
    private readonly double[] hairpinInit = new double[MaxTabulatedLoop + 1];
    private readonly double[] bulgeInit = new double[MaxTabulatedLoop + 1];
    private readonly double[] internalInit = new double[MaxTabulatedLoop + 1];
    private readonly double[] hairpinMismatch = new double[8 * 25];   // [type][x=s[i+1]][y=s[j-1]], terminal AU/GU not included
    private readonly double[] hairpinMismatchAu = new double[8 * 25];  // same, terminal AU/GU folded in (ViennaRNA mismatchH)
    private readonly double[] extMismatch = new double[8 * 25];       // [type(i,j)][s[i-1]][s[j+1]]
    private readonly double[] multiMismatch = new double[8 * 25];
    private readonly double[] dangle5 = new double[8 * 5];            // [type(i,j)][s[i-1]]
    private readonly double[] dangle3 = new double[8 * 5];            // [type(i,j)][s[j+1]]
    private readonly double[] int11 = new double[8 * 8 * 25];
    private readonly double[] int21 = new double[8 * 8 * 125];
    private readonly double[] int22 = new double[8 * 8 * 625];
    private readonly double[] mismatchInterior = new double[8 * 25];
    private readonly double[] mismatch1n = new double[8 * 25];
    private readonly double[] mismatch23 = new double[8 * 25];
    private readonly Dictionary<int, double> specialHairpin = new();

    private TurnerModel(double kelvin)
    {
        Temperature = kelvin;
        RT = GasConstantCalPerMolK * kelvin / 1000.0;
        TempFactor = kelvin / ReferenceTemperatureKelvin;
        IsReference = kelvin == ReferenceTemperatureKelvin;
        double f = TempFactor;

        // ΔG(T) = ΔH − (ΔH − ΔG37)·f, applied only to defined enthalpies (ViennaRNA RESCALE_dG).
        double Rs(double dG, int dHdcal)
            => IsReference || dHdcal >= Turner2004Enthalpies.Inf ? dG : dHdcal / 100.0 - (dHdcal / 100.0 - dG) * f;

        TerminalAu = Rs(Turner2004Parameters.TerminalAuGuPenalty, Turner2004Enthalpies.TerminalAudH);
        MlClosing = Rs(MlClosing37, Turner2004Enthalpies.MlClosingdH);
        MlIntern = Rs(MlIntern37, Turner2004Enthalpies.MlInterndH);
        MlBase = Rs(MlBase37, Turner2004Enthalpies.MlBasedH);
        ninio = Rs(0.6, Turner2004Enthalpies.NiniodH);
        lxcHairpinBulge = HairpinBulgeLxc37 * f;
        lxcInternal = InternalLoopLxc37 * f;

        // Stack: ours keyed "xy/zw" = 5'-x y-3' / 3'-z w-5' (outer pair x·z, inner pair y·w); ViennaRNA stack[type][rtype(type2)].
        var stackdH = Turner2004Enthalpies.StackdH;
        for (int t = 1; t <= 6; t++)
            for (int t2 = 1; t2 <= 6; t2++)
            {
                var (a, b) = TypeLetters[t];
                var (c, d) = TypeLetters[t2];
                double dg = Turner2004Parameters.StackingEnergies[$"{a}{c}/{b}{d}"];
                int rt2 = PairType(CodeOf(d), CodeOf(c));
                stack[t * 8 + t2] = Rs(dg, stackdH[t * 8 + rt2]);
            }

        for (int n = 3; n <= MaxTabulatedLoop; n++)
            hairpinInit[n] = Rs(Turner2004Parameters.HairpinLoopInitiation[n], Turner2004Enthalpies.HairpindH[n]);
        for (int n = 1; n <= MaxTabulatedLoop; n++)
            bulgeInit[n] = Rs(Turner2004Parameters.BulgeLoopInitiation[n], Turner2004Enthalpies.BulgedH[n]);
        for (int n = 4; n <= MaxTabulatedLoop; n++)
            internalInit[n] = Rs(Turner2004Parameters.InternalLoopInitiation[n], Turner2004Enthalpies.InternalLoopdH[n]);

        BuildMismatchTables();
        BuildDangleTables();
        BuildInteriorTables();
        BuildSpecialHairpins();
    }

    /// <summary>The model at <paramref name="kelvin"/> (cached; exactly the 37 °C model at 310.15 K).</summary>
    public static TurnerModel ForTemperature(double kelvin)
    {
        if (kelvin == ReferenceTemperatureKelvin) return At37;
        lock (CacheLock)
        {
            if (!Cache.TryGetValue(kelvin, out var m))
            {
                if (Cache.Count > 16) Cache.Clear();
                Cache[kelvin] = m = new TurnerModel(kelvin);
            }
            return m;
        }
    }

    // ---------------------------------------------------------------- encoding / pairing

    private static byte[] BuildPairTypeTable()
    {
        var t = new byte[25];
        for (int a = 0; a < 5; a++)
            for (int b = 0; b < 5; b++)
                t[a * 5 + b] = (byte)Turner2004InteriorLoopTables.PairType(BaseLetters[a], BaseLetters[b]);
        return t;
    }

    /// <summary>ViennaRNA base code of a nucleotide (A 1, C 2, G 3, U/T 4, anything else 0), case-insensitive.</summary>
    public static byte CodeOf(char c) => (byte)Turner2004InteriorLoopTables.BaseCode(c);

    /// <summary>Encodes a sequence (T read as U) into base codes.</summary>
    public static byte[] Encode(string seq)
    {
        var codes = new byte[seq.Length];
        for (int i = 0; i < codes.Length; i++) codes[i] = CodeOf(seq[i]);
        return codes;
    }

    /// <summary>ViennaRNA pair type of the ordered pair of base codes (1..6 canonical, 7 non-standard).</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    public static int PairType(int five, int three) => PairTypeTable[five * 5 + three];

    /// <summary>True for a canonical pair (Watson–Crick or G·U).</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    public static bool CanPair(int five, int three) => PairTypeTable[five * 5 + three] != 7;

    /// <summary>True when the pair type is A·U, U·A, G·U or U·G (terminal AU/GU penalty).</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    public static bool IsAuGu(int type) => type is >= 3 and <= 6;

    private double Round(double x) => IsReference ? Math.Round(x, 2) : x;

    // ---------------------------------------------------------------- table construction

    private void BuildMismatchTables()
    {
        var hdH = Turner2004Enthalpies.MismatchHdH;
        var edH = Turner2004Enthalpies.MismatchExtdH;
        var mdH = Turner2004Enthalpies.MismatchMdH;
        double f = TempFactor;
        double Rs(double dG, int dH) => IsReference || dH >= Turner2004Enthalpies.Inf ? dG : dH / 100.0 - (dH / 100.0 - dG) * f;

        for (int t = 1; t <= 6; t++)
        {
            var (c5, c3) = TypeLetters[t];
            // dG37 over A/C/G/U first, then the N (index 0) rows/columns = maximum over the substitutions.
            var hp = new double[25];
            var tm = new double[25];
            for (int x = 1; x <= 4; x++)
                for (int y = 1; y <= 4; y++)
                {
                    char a = BaseLetters[x], b = BaseLetters[y];
                    double e = Turner2004Parameters.TerminalMismatchEnergies[$"{c5}{a}{b}{c3}"];
                    tm[x * 5 + y] = e;
                    if ((a == 'U' && b == 'U') || (a == 'G' && b == 'A')) e += UuGaMismatchBonus;
                    if (a == 'G' && b == 'G') e += GgMismatchBonus;
                    hp[x * 5 + y] = e;
                }
            FillUnknownMax(hp);
            // Ext/ML mismatch stacks seen from the loop: closing pair read (s[j], s[i]) = (c5', c3') of the key, so the
            // table of pair type (i,j) = (b_i, b_j) holds TM["{b_j}{s[j+1]}{s[i-1]}{b_i}"] at [s[i-1]][s[j+1]].
            var (bi, bj) = TypeLetters[t];
            var ext = new double[25];
            for (int n5 = 1; n5 <= 4; n5++)
                for (int n3 = 1; n3 <= 4; n3++)
                    ext[n5 * 5 + n3] = Turner2004Parameters.TerminalMismatchEnergies[$"{bj}{BaseLetters[n3]}{BaseLetters[n5]}{bi}"];
            FillUnknownMax(ext);
            for (int k = 0; k < 25; k++)
            {
                // ViennaRNA folds the terminal AU/GU penalty of a hairpin (size ≥ 4) into mismatchH, with its own enthalpy.
                double hpAu = hp[k] + (IsAuGu(t) ? Turner2004Parameters.TerminalAuGuPenalty : 0.0);
                hairpinMismatchAu[t * 25 + k] = Rs(hpAu, hdH[t * 25 + k]);
                double auPart = IsAuGu(t) ? TerminalAu : 0.0;
                hairpinMismatch[t * 25 + k] = IsReference ? hp[k] : hairpinMismatchAu[t * 25 + k] - auPart;
                // Exterior / multibranch terminal mismatches never destabilise (ViennaRNA params.c: (mm > 0) ? 0 : mm).
                extMismatch[t * 25 + k] = Math.Min(0.0, Rs(ext[k], edH[t * 25 + k]));
                multiMismatch[t * 25 + k] = Math.Min(0.0, Rs(ext[k], mdH[t * 25 + k]));
            }
        }
    }

    // Replaces row 0 / column 0 (unknown base N) with the maximum over the A/C/G/U substitutions.
    private static void FillUnknownMax(double[] m)
    {
        for (int x = 1; x <= 4; x++)
        {
            double mx = double.NegativeInfinity;
            for (int y = 1; y <= 4; y++) mx = Math.Max(mx, m[x * 5 + y]);
            m[x * 5] = mx;
        }
        for (int y = 0; y <= 4; y++)
        {
            double mx = double.NegativeInfinity;
            for (int x = 1; x <= 4; x++) mx = Math.Max(mx, m[x * 5 + y]);
            m[y] = mx;
        }
    }

    private void BuildDangleTables()
    {
        var d5dH = Turner2004Enthalpies.Dangle5dH;
        var d3dH = Turner2004Enthalpies.Dangle3dH;
        double f = TempFactor;
        // Dangling ends never destabilise (ViennaRNA params.c: "must be <= 0"): a rescaled positive value is set to 0.
        double Rs(double dG, int dH) => IsReference ? dG : Math.Min(0.0, dH / 100.0 - (dH / 100.0 - dG) * f);
        for (int t = 1; t <= 6; t++)
        {
            var (bi, bj) = TypeLetters[t];
            double mx5 = double.NegativeInfinity, mx3 = double.NegativeInfinity;
            for (int x = 1; x <= 4; x++)
            {
                double e5 = Turner2004Parameters.DanglingEnd5Energies[$"{bj}{BaseLetters[x]}{bi}"];
                double e3 = Turner2004Parameters.DanglingEnd3Energies[$"{bj}{BaseLetters[x]}{bi}"];
                dangle5[t * 5 + x] = Rs(e5, d5dH[t * 5 + x]);
                dangle3[t * 5 + x] = Rs(e3, d3dH[t * 5 + x]);
                mx5 = Math.Max(mx5, e5);
                mx3 = Math.Max(mx3, e3);
            }
            dangle5[t * 5] = Rs(mx5, d5dH[t * 5]);
            dangle3[t * 5] = Rs(mx3, d3dH[t * 5]);
        }
    }

    private void BuildInteriorTables()
    {
        double f = TempFactor;
        double R(int dG37, int dH) => IsReference ? dG37 / 100.0 : (dH / 100.0 - (dH / 100.0 - dG37 / 100.0) * f);
        var i11 = Turner2004Enthalpies.Int11dH;
        var i21 = Turner2004Enthalpies.Int21dH;
        var i22 = Turner2004Enthalpies.Int22dH;
        for (int t = 1; t <= 7; t++)
            for (int t2 = 1; t2 <= 7; t2++)
            {
                for (int a = 0; a < 5; a++)
                    for (int b = 0; b < 5; b++)
                    {
                        int ix = ((t * 8 + t2) * 5 + a) * 5 + b;
                        int dH = i11[ix];
                        int dG = Turner2004InteriorLoopTables.Int11(t, t2, a, b);
                        int11[ix] = dH >= Turner2004Enthalpies.Inf ? dG / 100.0 : R(dG, dH);
                        for (int c = 0; c < 5; c++)
                        {
                            int ix3 = (((t * 8 + t2) * 5 + a) * 5 + b) * 5 + c;
                            int dH3 = i21[ix3];
                            int dG3 = Turner2004InteriorLoopTables.Int21(t, t2, a, b, c);
                            int21[ix3] = dH3 >= Turner2004Enthalpies.Inf ? dG3 / 100.0 : R(dG3, dH3);
                            for (int d = 0; d < 5; d++)
                            {
                                int ix4 = ((((t * 8 + t2) * 5 + a) * 5 + b) * 5 + c) * 5 + d;
                                int dH4 = i22[ix4];
                                int dG4 = Turner2004InteriorLoopTables.Int22(t, t2, a, b, c, d);
                                int22[ix4] = dH4 >= Turner2004Enthalpies.Inf ? dG4 / 100.0 : R(dG4, dH4);
                            }
                        }
                    }
            }
        var mI = Turner2004Enthalpies.MismatchIdH;
        var m1 = Turner2004Enthalpies.Mismatch1nIdH;
        var m23 = Turner2004Enthalpies.Mismatch23IdH;
        for (int t = 1; t <= 7; t++)
            for (int x = 0; x < 5; x++)
                for (int y = 0; y < 5; y++)
                {
                    int ix = (t * 5 + x) * 5 + y;
                    int o = (t * 5 + x) * 5 + y;
                    mismatchInterior[o] = IsReference ? Turner2004InteriorLoopTables.MismatchInterior(t, x, y) / 100.0
                        : R(Turner2004InteriorLoopTables.MismatchInterior(t, x, y), mI[ix]);
                    mismatch1n[o] = IsReference ? Turner2004InteriorLoopTables.Mismatch1nInterior(t, x, y) / 100.0
                        : R(Turner2004InteriorLoopTables.Mismatch1nInterior(t, x, y), m1[ix]);
                    mismatch23[o] = IsReference ? Turner2004InteriorLoopTables.Mismatch23Interior(t, x, y) / 100.0
                        : R(Turner2004InteriorLoopTables.Mismatch23Interior(t, x, y), m23[ix]);
                }
    }

    private static int SpecialKey(ReadOnlySpan<byte> codes)
    {
        int k = 0;
        foreach (byte c in codes) k = k * 5 + c;
        return codes.Length * 1_000_000 + k;
    }

    private void BuildSpecialHairpins()
    {
        var dHByLoop = new Dictionary<string, int>();
        for (int k = 0; k < Turner2004Enthalpies.Triloops.Length; k++) dHByLoop[Turner2004Enthalpies.Triloops[k]] = Turner2004Enthalpies.TriloopsdH[k];
        for (int k = 0; k < Turner2004Enthalpies.Tetraloops.Length; k++) dHByLoop[Turner2004Enthalpies.Tetraloops[k]] = Turner2004Enthalpies.TetraloopsdH[k];
        for (int k = 0; k < Turner2004Enthalpies.Hexaloops.Length; k++) dHByLoop[Turner2004Enthalpies.Hexaloops[k]] = Turner2004Enthalpies.HexaloopsdH[k];
        foreach (var (loop, dG) in Turner2004Parameters.SpecialHairpinLoops)
        {
            int dH = dHByLoop[loop]; // every NNDB special loop is a ViennaRNA special loop
            var codes = new byte[loop.Length];
            for (int i = 0; i < codes.Length; i++) codes[i] = CodeOf(loop[i]);
            specialHairpin[SpecialKey(codes)] = IsReference ? dG : dH / 100.0 - (dH / 100.0 - dG) * TempFactor;
        }
    }

    // ---------------------------------------------------------------- loop terms

    /// <summary>
    /// Stacking energy of the pair of types (<paramref name="type"/> outer, <paramref name="type2"/> inner, both read
    /// 5'→3'); 0 when either is not a canonical pair.
    /// </summary>
    public double Stack(int type, int type2) => stack[type * 8 + type2];

    /// <summary>
    /// Hairpin loop free energy of <c>closing5' loop closing3'</c> (<paramref name="hairpin"/> = base codes of the
    /// closing pair and the loop in between). By default the closing pair's terminal AU/GU penalty is NOT included
    /// (NNDB counts it in the helix); with <paramref name="includeTerminalAu"/> it is, the way ViennaRNA folds it into
    /// <c>mismatchH</c> (size ≥ 4, with that table's enthalpy) or adds it (size 3). A special tri/tetra/hexaloop value is
    /// the measured loop total and contains the closure already (<paramref name="isSpecial"/>). Loops shorter than 3 nt
    /// return the prohibitive 100.
    /// </summary>
    public double HairpinLoop(ReadOnlySpan<byte> hairpin, bool specialGuClosure, out bool isSpecial, bool includeTerminalAu = false)
    {
        isSpecial = false;
        int size = hairpin.Length - 2;
        if (size >= 3 && (size == 3 || size == 4 || size == 6) && specialHairpin.Count > 0)
        {
            bool known = true;
            foreach (byte c in hairpin) if (c == 0) { known = false; break; }
            if (known && specialHairpin.TryGetValue(SpecialKey(hairpin), out double sp))
            {
                isSpecial = true;
                return Round(sp);
            }
        }

        double e;
        if (size >= 3 && size <= MaxTabulatedLoop) e = hairpinInit[size];
        else if (size > MaxTabulatedLoop) e = hairpinInit[9] + lxcHairpinBulge * Math.Log(size / 9.0); // NNDB: ΔG(9) + 1.75RT·ln(n/9)
        else return 100.0; // loops < 3 nt are sterically impossible

        byte c5 = hairpin[0], c3 = hairpin[^1];
        if (size >= 4)
        {
            int type = PairType(c5, c3);
            e += (includeTerminalAu ? hairpinMismatchAu : hairpinMismatch)[type * 25 + hairpin[1] * 5 + hairpin[^2]];
            if (specialGuClosure && c5 == 3 && c3 == 4) e += SpecialGuClosureBonus;
        }

        bool allC = true;
        for (int k = 1; k <= size; k++) if (hairpin[k] != 2) { allC = false; break; }
        if (allC) e += size == 3 ? AllCPenalty3Nt : AllCPenaltyPerNt * size + AllCPenaltyConstantGt3;

        if (includeTerminalAu && size == 3 && IsAuGu(PairType(c5, c3))) e += TerminalAu;
        return Round(e);
    }

    /// <summary>Hairpin loop closed by (i,j) inside the sequence <paramref name="s"/>, incl. the closing pair's terminal AU/GU penalty.</summary>
    public double HairpinTerm(byte[] s, int i, int j)
    {
        bool specialGu = s[i] == 3 && s[j] == 4 && i >= 2 && s[i - 1] == 3 && s[i - 2] == 3;
        return HairpinLoop(s.AsSpan(i, j - i + 1), specialGu, out _, includeTerminalAu: true);
    }

    /// <summary>Bulge loop of <paramref name="size"/> nucleotides (NNDB bulge rules; special C, degeneracy −RT·ln(states)).</summary>
    public double Bulge(int size, int bulged, int outer5, int outer3, int inner5, int inner3,
        int numStates = 1, bool? bulgeOn5PrimeStrand = null)
    {
        double e = size <= MaxTabulatedLoop ? bulgeInit[size] : bulgeInit[6] + lxcHairpinBulge * Math.Log(size / 6.0);
        int t1 = PairType(outer5, outer3), t2 = PairType(inner5, inner3);
        if (size == 1)
        {
            e += Stack(t1, t2); // stacking as if the bulge were absent
            bool on5 = outer5 == 2 || inner5 == 2;
            bool on3 = outer3 == 2 || inner3 == 2;
            bool adjacentC = bulgeOn5PrimeStrand switch { true => on5, false => on3, null => on5 && on3 };
            if (bulged == 2 && adjacentC) e += BulgeSpecialCBonus;
            if (numStates > 1) e -= RT * Math.Log(numStates);
        }
        else
        {
            if (IsAuGu(t1)) e += TerminalAu;
            if (IsAuGu(t2)) e += TerminalAu;
        }
        return Round(e);
    }

    private double InternalInit(int n)
        => n <= MaxTabulatedLoop ? internalInit[n] : internalInit[6] + lxcInternal * Math.Log(n / 6.0);

    private double Asymmetry(int d) => Math.Min(MaxNinio, ninio * d);

    /// <summary>
    /// Degree-2 loop with n1 / n2 unpaired nucleotides (ViennaRNA <c>vrna_E_internal</c>): stack (0,0), bulge (one side
    /// empty) or interior loop. <paramref name="outer5"/>/<paramref name="outer3"/> = bases of the outer pair (i,j),
    /// <paramref name="inner5"/>/<paramref name="inner3"/> = bases of the inner pair (p,q); the four mismatch bases are
    /// s[i+1], s[j−1], s[p−1], s[q+1].
    /// </summary>
    public double InternalLoop(int n1, int n2, int outer5, int outer3, int inner5, int inner3,
        int mm5Outer, int mm3Outer, int mm5Inner, int mm3Inner, int bulgeStates = 1)
    {
        int nl = Math.Max(n1, n2), ns = Math.Min(n1, n2);
        if (nl == 0) return Stack(PairType(outer5, outer3), PairType(inner5, inner3));
        if (ns == 0)
            return Bulge(nl, n1 > 0 ? mm5Outer : mm3Outer, outer5, outer3, inner5, inner3, bulgeStates, n1 > 0);

        int type = PairType(outer5, outer3);
        int type2 = PairType(inner3, inner5);
        int si1 = mm5Outer, sj1 = mm3Outer, sp1 = mm5Inner, sq1 = mm3Inner;

        if (ns == 1)
        {
            if (nl == 1) return int11[((type * 8 + type2) * 5 + si1) * 5 + sj1];
            if (nl == 2)
                return n1 == 1
                    ? int21[(((type * 8 + type2) * 5 + si1) * 5 + sq1) * 5 + sj1]
                    : int21[(((type2 * 8 + type) * 5 + sq1) * 5 + si1) * 5 + sp1];
            return Round(InternalInit(nl + 1) + Asymmetry(nl - ns)
                + mismatch1n[(type * 5 + si1) * 5 + sj1] + mismatch1n[(type2 * 5 + sq1) * 5 + sp1]);
        }
        if (ns == 2 && nl == 2)
            return int22[((((type * 8 + type2) * 5 + si1) * 5 + sp1) * 5 + sq1) * 5 + sj1];
        if (ns == 2 && nl == 3)
            return Round(InternalInit(5) + ninio
                + mismatch23[(type * 5 + si1) * 5 + sj1] + mismatch23[(type2 * 5 + sq1) * 5 + sp1]);
        return Round(InternalInit(nl + ns) + Asymmetry(nl - ns)
            + mismatchInterior[(type * 5 + si1) * 5 + sj1] + mismatchInterior[(type2 * 5 + sq1) * 5 + sp1]);
    }

    /// <summary>Loop closed by (i,j) with inner pair (ip,jp) inside <paramref name="s"/>, including the n = 1 bulge degeneracy.</summary>
    public double InteriorTerm(byte[] s, int i, int j, int ip, int jp)
    {
        int n1 = ip - i - 1, n2 = j - jp - 1;
        int states = 1;
        if (n1 + n2 == 1)
        {
            // NNDB bulge loops, Example 1: identical neighbouring bases give equivalent bulge positions.
            int bulged = n1 == 1 ? s[i + 1] : s[j - 1];
            if (n1 == 1)
                states += (s[i] == bulged ? 1 : 0) + (s[ip] == bulged ? 1 : 0);
            else
                states += (s[j] == bulged ? 1 : 0) + (s[jp] == bulged ? 1 : 0);
        }
        return InternalLoop(n1, n2, s[i], s[j], s[ip], s[jp],
            i + 1 < s.Length ? s[i + 1] : 0, j > 0 ? s[j - 1] : 0, ip > 0 ? s[ip - 1] : 0, jp + 1 < s.Length ? s[jp + 1] : 0, states);
    }

    /// <summary>True when (i,j) closes the special 5'GGUC/3'CUGG three-stack.</summary>
    public static bool IsSpecialGgucCuggStack(byte[] s, int i, int j)
        => j - i >= 7 && s[i] == 3 && s[j] == 2 && s[i + 1] == 3 && s[j - 1] == 4
           && s[i + 2] == 4 && s[j - 2] == 3 && s[i + 3] == 2 && s[j - 3] == 3;

    /// <summary>Helix end (i,j) facing the exterior loop, d2 dangles (ViennaRNA <c>vrna_E_exterior_stem</c>).</summary>
    public double ExteriorStem(byte[] s, int i, int j)
    {
        int type = PairType(s[i], s[j]);
        double e = IsAuGu(type) ? TerminalAu : 0.0;
        bool has5 = i > 0, has3 = j < s.Length - 1;
        if (has5 && has3) e += extMismatch[type * 25 + s[i - 1] * 5 + s[j + 1]];
        else if (has5) e += dangle5[type * 5 + s[i - 1]];
        else if (has3) e += dangle3[type * 5 + s[j + 1]];
        return e;
    }

    /// <summary>Branch helix (i,j) of a multibranch loop, d2 (ViennaRNA <c>E_MLstem(type, s[i−1], s[j+1])</c>).</summary>
    public double MultiloopStem(byte[] s, int i, int j)
    {
        int type = PairType(s[i], s[j]);
        return MlIntern + (IsAuGu(type) ? TerminalAu : 0.0) + multiMismatch[type * 25 + s[i - 1] * 5 + s[j + 1]];
    }

    /// <summary>Closing pair (i,j) of a multibranch loop, d2 (ViennaRNA <c>E_MLstem(rtype, s[j−1], s[i+1]) + ML_closing</c>).</summary>
    public double MultiloopClosing(byte[] s, int i, int j)
    {
        int type = PairType(s[j], s[i]);
        return MlClosing + MlIntern + (IsAuGu(type) ? TerminalAu : 0.0) + multiMismatch[type * 25 + s[j - 1] * 5 + s[i + 1]];
    }

    // ---------------------------------------------------------------- public-API helpers (case/T-normalised chars)

    /// <summary>Terminal mismatch (exterior / multibranch) of the helix end <c>closing5 mismatch5 / closing3 mismatch3</c>.</summary>
    public double TerminalMismatch(char closing5, char closing3, char mismatch5, char mismatch3)
    {
        int type = PairType(CodeOf(closing3), CodeOf(closing5));
        if (type == 7) return 0.0;
        // table of pair (i,j) = (closing3, closing5) at [s[i-1] = mismatch3][s[j+1] = mismatch5]
        return multiMismatch[type * 25 + CodeOf(mismatch3) * 5 + CodeOf(mismatch5)];
    }

    /// <summary>Single dangling base on the helix end <c>closing5 / closing3</c>.</summary>
    public double Dangle(char closing5, char closing3, char dangling, bool is3Prime)
    {
        int type = PairType(CodeOf(closing3), CodeOf(closing5));
        if (type == 7) return 0.0;
        return (is3Prime ? dangle3 : dangle5)[type * 5 + CodeOf(dangling)];
    }
}
