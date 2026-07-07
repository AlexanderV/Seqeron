// PROTMOTIF-HMM-001 (08_DIFFERENTIAL, strategy BRUTE) — differential test of the Plan7 local-mode
// Forward/Backward DP against an INDEPENDENT explicit-path-enumeration oracle.
//
// WHY THIS EXISTS. The mutation survivors in Plan7ProfileHmm.RunForwardBackward live in the
// Backward / posterior-decoding half of the DP (bM/bI/bD/bN/bB/bE/bC/bJ). Those cells are NOT
// observable through the public scores or domain envelopes, which RE-CONVERGE under index/threshold
// mutation — so exact-output tests cannot distinguish a corrupted Backward index (they are
// "equivalent mutants" only with respect to the public API). The engine therefore exposes the raw
// Forward/Backward matrices to this test (DebugLocalForwardBackward, internal via IVT).
//
// THE ORACLE. We recompute every Forward and Backward cell from FIRST PRINCIPLES by explicitly
// enumerating Plan7 state paths — NOT by re-running the DP recurrence. Forward F[i][s] is the
// log-sum over all path PREFIXES that end in state s having emitted i residues; Backward B[i][s] is
// the log-sum over all path SUFFIXES from s (i emitted) that emit residues i+1..L and reach T. This
// is a genuinely independent computation (path enumeration vs. dynamic programming), so a mutated
// index in the recurrence yields a cell that diverges from the enumerated truth -> killed.
//
// SELF-VALIDATION. The oracle's own total log-sum is asserted equal to the engine's ForwardNats AND
// to the public LocalForwardScore, so an oracle bug cannot silently masquerade as a passing test.
//
// Topology / parameters are read verbatim from the engine snapshot (LocalEntryLn occupancy weights,
// _transitionLn, _matchEmissionLn, the static HMMER amino background, and the N/J/C length config),
// because those are the MODEL (ground truth), not the algorithm under test. Sources: Eddy (2011)
// PLoS Comput Biol 7:e1002195; HMMER generic_fwdback.c (p7_GForward/p7_GBackward, local esc=0);
// Durbin et al. (1998) §3.2 (Forward/Backward).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
[Category("PROTMOTIF-HMM-001")]
public class Plan7ProfileHmm_ForwardBackwardDifferential_Tests
{
    private const double NegInf = double.NegativeInfinity;
    private const double Tol = 1e-11;

    // ---- amino alphabet & small HMM builders (parameters fully under test control) ----

    private const string Alpha = "ACDEFGHIKLMNPQRSTVWY";

    private static double[] EmissionLn(double pA, double pC)
    {
        var row = new double[20];
        for (int i = 0; i < 20; i++)
        {
            double p = Alpha[i] == 'A' ? pA : Alpha[i] == 'C' ? pC : 1e-9;
            row[i] = Math.Log(p);
        }
        return row;
    }

    private static double[] TransLn(double mm, double mi, double md, double im, double ii, double dm, double dd)
    {
        double L(double p) => double.IsNaN(p) ? NegInf : Math.Log(p);
        return new[] { L(mm), L(mi), L(md), L(im), L(ii), L(dm), L(dd) };
    }

    // A toy HMM bundled with the parameter arrays it was built from, so the glocal oracle can read
    // them (the glocal Forward/Viterbi uses the INSTANCE background passed here, not the static one).
    private sealed record ToyHmm(Plan7ProfileHmm Hmm, int M, double[][] Match, double[][] Insert, double[][] Trans, double[] Bg);

    // 2-node toy HMM over {A,C} (identical to the canonical hand HMM, so ViterbiScore("AC") == 0.5187…).
    private static ToyHmm BuildHmm2()
    {
        var bg = EmissionLn(0.6, 0.4);
        var match = new double[3][];
        var insert = new double[3][];
        var trans = new double[3][];
        insert[0] = EmissionLn(0.6, 0.4);
        trans[0] = TransLn(0.9, 0.05, 0.05, 0.5, 0.5, 1.0, double.NaN);
        match[1] = EmissionLn(0.7, 0.3); insert[1] = EmissionLn(0.6, 0.4);
        trans[1] = TransLn(0.8, 0.1, 0.1, 0.5, 0.5, 1.0, double.NaN);
        match[2] = EmissionLn(0.2, 0.8); insert[2] = EmissionLn(0.6, 0.4);
        trans[2] = TransLn(1.0, double.NaN, double.NaN, 1.0, 0.5, 1.0, double.NaN);
        return new ToyHmm(Plan7ProfileHmm.FromLogParameters("toy2", 2, match, insert, trans, bg), 2, match, insert, trans, bg);
    }

    // 1-node toy HMM — the m==1 edge case (the `m >= 1` delete-init guard is equivalent for m>=2).
    private static ToyHmm BuildHmm1()
    {
        var bg = EmissionLn(0.6, 0.4);
        var match = new double[2][];
        var insert = new double[2][];
        var trans = new double[2][];
        insert[0] = EmissionLn(0.6, 0.4);
        trans[0] = TransLn(0.9, 0.05, 0.05, 0.5, 0.5, 1.0, double.NaN);
        match[1] = EmissionLn(0.7, 0.3); insert[1] = EmissionLn(0.6, 0.4);
        trans[1] = TransLn(1.0, double.NaN, double.NaN, 1.0, 0.5, 1.0, double.NaN);
        return new ToyHmm(Plan7ProfileHmm.FromLogParameters("toy1", 1, match, insert, trans, bg), 1, match, insert, trans, bg);
    }

    // 3-node toy HMM — exercises D2/D3, M3, I2 and longer delete chains.
    private static ToyHmm BuildHmm3()
    {
        var bg = EmissionLn(0.6, 0.4);
        var match = new double[4][];
        var insert = new double[4][];
        var trans = new double[4][];
        insert[0] = EmissionLn(0.6, 0.4);
        trans[0] = TransLn(0.85, 0.1, 0.05, 0.5, 0.5, 1.0, double.NaN);
        match[1] = EmissionLn(0.7, 0.3); insert[1] = EmissionLn(0.6, 0.4);
        trans[1] = TransLn(0.75, 0.15, 0.1, 0.5, 0.5, 0.6, 0.4);
        match[2] = EmissionLn(0.3, 0.7); insert[2] = EmissionLn(0.6, 0.4);
        trans[2] = TransLn(0.7, 0.2, 0.1, 0.5, 0.5, 0.6, 0.4);
        match[3] = EmissionLn(0.5, 0.5); insert[3] = EmissionLn(0.6, 0.4);
        trans[3] = TransLn(1.0, double.NaN, double.NaN, 1.0, 0.5, 1.0, double.NaN);
        return new ToyHmm(Plan7ProfileHmm.FromLogParameters("toy3", 3, match, insert, trans, bg), 3, match, insert, trans, bg);
    }

    // ---- independent enumeration oracle ----

    private sealed class Oracle
    {
        private readonly int _m;
        private readonly int _len;
        private readonly int[] _dsq;            // 1..len
        private readonly double[] _bm;          // local entry, index 1..m
        private readonly double[][] _trans;     // [0..m][7]
        private readonly double[][] _match;     // [0..m][20]
        private readonly double[] _bg;          // [20]
        private readonly double _xLoop, _xMove, _eMove, _eLoop;
        private readonly Dictionary<(char, int, int), double> _fwd = new();

        public Oracle(Plan7ProfileHmm.DebugForwardBackwardSnapshot s, int m)
        {
            _m = m;
            _len = s.Dsq.Length - 1;
            _dsq = s.Dsq;
            _bm = s.LocalEntryLn;
            _trans = s.TransitionLn;
            _match = s.MatchEmissionLn;
            _bg = s.BackgroundLn;
            _xLoop = s.XLoop; _xMove = s.XMove; _eMove = s.EMove; _eLoop = s.ELoop;
            ForwardDfs('N', 0, 0, 0.0);
        }

        private static double LogSumExp(double a, double b)
        {
            if (double.IsNegativeInfinity(a)) return b;
            if (double.IsNegativeInfinity(b)) return a;
            double hi = Math.Max(a, b), lo = Math.Min(a, b);
            return hi + Math.Log(1.0 + Math.Exp(lo - hi));
        }

        // Emission log-odds added when ENTERING an emitting state v of node k at residue position i.
        // Only match states emit a non-zero log-odds (insert/N/J/C emit background -> 0).
        private double Emit(char v, int k, int i)
        {
            if (v != 'M') return 0.0;
            int x = _dsq[i];
            if (x < 0) return 0.0;
            double me = _match[k][x];
            return double.IsNegativeInfinity(me) ? NegInf : me - _bg[x];
        }

        // Successors of (kind,k): (vKind, vK, transitionLogOdds, emits).
        private IEnumerable<(char V, int Vk, double L, bool Emits)> Succ(char kind, int k)
        {
            switch (kind)
            {
                case 'N':
                    yield return ('N', 0, _xLoop, true);
                    yield return ('B', 0, _xMove, false);
                    break;
                case 'B':
                    for (int j = 1; j <= _m; j++) yield return ('M', j, _bm[j], true);
                    break;
                case 'M':
                    if (k < _m)
                    {
                        yield return ('M', k + 1, _trans[k][0], true);  // TMM
                        yield return ('I', k, _trans[k][1], true);      // TMI
                        yield return ('D', k + 1, _trans[k][2], false); // TMD
                    }
                    yield return ('E', 0, 0.0, false);                  // local exit esc=0
                    break;
                case 'I':
                    yield return ('M', k + 1, _trans[k][3], true);      // TIM
                    yield return ('I', k, _trans[k][4], true);          // TII
                    break;
                case 'D':
                    if (k < _m)
                    {
                        yield return ('M', k + 1, _trans[k][5], true);  // TDM
                        yield return ('D', k + 1, _trans[k][6], false); // TDD
                    }
                    yield return ('E', 0, 0.0, false);
                    break;
                case 'E':
                    yield return ('C', 0, _eMove, false);
                    yield return ('J', 0, _eLoop, false);
                    break;
                case 'J':
                    yield return ('J', 0, _xLoop, true);
                    yield return ('B', 0, _xMove, false);
                    break;
                case 'C':
                    yield return ('C', 0, _xLoop, true);
                    yield return ('T', 0, _xMove, false);
                    break;
                case 'T':
                    break;
            }
        }

        private void ForwardDfs(char kind, int k, int i, double pp)
        {
            var key = (kind, k, i);
            _fwd[key] = _fwd.TryGetValue(key, out var cur) ? LogSumExp(cur, pp) : pp;
            if (kind == 'T') return;
            foreach (var (v, vk, l, emits) in Succ(kind, k))
            {
                if (double.IsNegativeInfinity(l)) continue;
                if (emits)
                {
                    if (i + 1 > _len) continue;
                    double e = Emit(v, vk, i + 1);
                    if (double.IsNegativeInfinity(e)) continue;
                    ForwardDfs(v, vk, i + 1, pp + l + e);
                }
                else
                {
                    ForwardDfs(v, vk, i, pp + l);
                }
            }
        }

        public double Forward(char kind, int k, int i)
            => _fwd.TryGetValue((kind, k, i), out var v) ? v : NegInf;

        public double Total => Forward('T', 0, _len);

        // Suffix log-sum from (kind,k) having emitted i residues, emitting i+1..len and reaching T.
        public double Backward(char kind, int k, int i)
        {
            if (kind == 'T') return i == _len ? 0.0 : NegInf;
            double acc = NegInf;
            foreach (var (v, vk, l, emits) in Succ(kind, k))
            {
                if (double.IsNegativeInfinity(l)) continue;
                if (emits)
                {
                    if (i + 1 > _len) continue;
                    double e = Emit(v, vk, i + 1);
                    if (double.IsNegativeInfinity(e)) continue;
                    acc = LogSumExp(acc, l + e + Backward(v, vk, i + 1));
                }
                else
                {
                    acc = LogSumExp(acc, l + Backward(v, vk, i));
                }
            }
            return acc;
        }
    }

    private static void AssertCell(double expected, double actual, string what)
    {
        if (double.IsNegativeInfinity(expected) && double.IsNegativeInfinity(actual)) return;
        Assert.That(actual, Is.EqualTo(expected).Within(Tol),
            $"{what}: expected {expected} (enumeration), got {actual} (engine)");
    }

    public static IEnumerable<TestCaseData> Cases()
    {
        yield return new TestCaseData(2, "AC", false).SetName("Hmm2_AC_unihit");
        yield return new TestCaseData(2, "AC", true).SetName("Hmm2_AC_multihit");
        yield return new TestCaseData(2, "ACA", true).SetName("Hmm2_ACA_multihit");
        yield return new TestCaseData(2, "ACAC", true).SetName("Hmm2_ACAC_multihit");
        yield return new TestCaseData(3, "AC", true).SetName("Hmm3_AC_multihit");
        yield return new TestCaseData(3, "ACA", false).SetName("Hmm3_ACA_unihit");
        yield return new TestCaseData(3, "ACAC", true).SetName("Hmm3_ACAC_multihit");
    }

    [TestCaseSource(nameof(Cases))]
    public void ForwardBackwardMatrices_MatchIndependentPathEnumeration(int nodes, string seq, bool multihit)
    {
        var hmm = (nodes == 2 ? BuildHmm2() : BuildHmm3()).Hmm;
        var snap = hmm.DebugLocalForwardBackward(seq, multihit);
        var fb = snap.Fb;
        int m = nodes, len = seq.Length;
        var oracle = new Oracle(snap, m);

        // (0) self-validation: the oracle's own total must equal the engine's Forward total, and (in
        // multihit mode, which is the LocalForwardScore configuration) the public bit-path score too.
        AssertCell(oracle.Total, fb.ForwardNats, "oracle.Total vs engine.ForwardNats");
        if (multihit)
            AssertCell(fb.ForwardNats, hmm.LocalForwardScore(seq), "engine.ForwardNats vs LocalForwardScore");

        // (1) Forward cells, rows 1..len.
        for (int i = 1; i <= len; i++)
        {
            AssertCell(oracle.Forward('N', 0, i), fb.FN[i], $"FN[{i}]");
            AssertCell(oracle.Forward('B', 0, i), fb.FB[i], $"FB[{i}]");
            AssertCell(oracle.Forward('E', 0, i), fb.FE[i], $"FE[{i}]");
            AssertCell(oracle.Forward('C', 0, i), fb.FC[i], $"FC[{i}]");
            AssertCell(oracle.Forward('J', 0, i), fb.FJ[i], $"FJ[{i}]");
            for (int k = 1; k <= m; k++)
            {
                AssertCell(oracle.Forward('M', k, i), fb.FM[i][k], $"FM[{i}][{k}]");
                AssertCell(oracle.Forward('D', k, i), fb.FD[i][k], $"FD[{i}][{k}]");
                if (k < m) AssertCell(oracle.Forward('I', k, i), fb.FI[i][k], $"FI[{i}][{k}]");
            }
        }

        // (2) Backward cells, rows 1..len — the survivor-bearing half.
        for (int i = 1; i <= len; i++)
        {
            AssertCell(oracle.Backward('N', 0, i), fb.BN[i], $"BN[{i}]");
            AssertCell(oracle.Backward('B', 0, i), fb.BB[i], $"BB[{i}]");
            AssertCell(oracle.Backward('E', 0, i), fb.BE[i], $"BE[{i}]");
            AssertCell(oracle.Backward('C', 0, i), fb.BC[i], $"BC[{i}]");
            AssertCell(oracle.Backward('J', 0, i), fb.BJ[i], $"BJ[{i}]");
            for (int k = 1; k <= m; k++)
            {
                AssertCell(oracle.Backward('M', k, i), fb.BM[i][k], $"BM[{i}][{k}]");
                AssertCell(oracle.Backward('D', k, i), fb.BD[i][k], $"BD[{i}][{k}]");
                if (k < m) AssertCell(oracle.Backward('I', k, i), fb.BI[i][k], $"BI[{i}][{k}]");
            }
        }

        // (3) Forward-Backward consistency: at every row i, the posterior over the emitting states that
        // produced residue i sums to 1 — i.e. log-sum of F'[i][s]+B[i][s] over those states equals the
        // total (a recurrence-free law; cf. generic_decoding.c p7_GDomainDecoding). For the N/J/C
        // special states the forward partner of "emitted residue i" is F*[i-1]+XLoop (the self-loop that
        // emitted residue i), NOT F*[i] — FC[i]/FJ[i] also carry the silent E->C / E->J entry mass.
        for (int i = 1; i <= len; i++)
        {
            double acc = NegInf;
            void Add(double f, double b) { if (!double.IsNegativeInfinity(f) && !double.IsNegativeInfinity(b)) acc = LseTest(acc, f + b); }
            Add(fb.FN[i - 1] + fb.XLoop, fb.BN[i]);
            Add(fb.FC[i - 1] + fb.XLoop, fb.BC[i]);
            Add(fb.FJ[i - 1] + fb.XLoop, fb.BJ[i]);
            for (int k = 1; k <= m; k++)
            {
                Add(fb.FM[i][k], fb.BM[i][k]);
                if (k < m) Add(fb.FI[i][k], fb.BI[i][k]);
            }
            AssertCell(fb.ForwardNats, acc, $"FB-consistency row {i}");
        }
    }

    // ---- glocal Viterbi/Forward differential (a SEPARATE DP from the local Forward/Backward) ----
    //
    // ViterbiScore/ForwardScore are HMMER glocal full-length scores: the path must enter at the model
    // BEGIN (B->M1 or B->D1, NO local entry at an arbitrary M_k), traverse every node, and exit only
    // from M_m or D_m (no N/C/J flanking, no multihit). We enumerate every such path and compare the
    // log-sum (Forward) and max (Viterbi) to the public scores. Glocal emission log-odds use the
    // INSTANCE background passed to the model (not the static HMMER background), so the oracle reads it
    // from the ToyHmm. Independent of the glocal DP recurrence -> kills its index/guard mutants.

    private static int Res(char c)
    {
        int i = Alpha.IndexOf(char.ToUpperInvariant(c));
        return i; // -1 if outside the 20-letter amino alphabet
    }

    private static double GlocalEmit(ToyHmm t, char v, int k, int a)
    {
        if (a < 0) return 0.0;
        double e = v == 'M' ? t.Match[k][a] : v == 'I' ? t.Insert[k][a] : double.NaN;
        if (double.IsNaN(e)) return 0.0; // non-emitting destination
        return double.IsNegativeInfinity(e) ? NegInf : e - t.Bg[a];
    }

    private static IEnumerable<(char V, int Vk, double L, bool Emits)> GlocalSucc(ToyHmm t, char kind, int k)
    {
        int m = t.M;
        switch (kind)
        {
            case 'S':
                // Absolute start (i=0). The all-delete entry B->D1 is available ONLY here — RunGlocal
                // sets prevD[1] at i=0 and never re-derives D1 inside the loop (line ~353: fromM=-inf
                // for node 1). The match entry, by contrast, floats (handled via B below).
                yield return ('B', 0, 0.0, false);            // enter the (floating) begin
                yield return ('D', 1, t.Trans[0][2], false);  // all-delete start (no free prefix)
                break;
            case 'B':
                // The BEGIN is score 0 at EVERY column (RunGlocal line ~335: fromM = 0 for node 1),
                // so M1 may be entered at any residue — the prefix floats free. Zero-cost B->B skip.
                yield return ('B', 0, 0.0, true);             // skip a leading residue (free)
                yield return ('M', 1, t.Trans[0][0], true);   // B->M1 (TMM slot of BEGIN node)
                break;
            case 'M':
                if (k < m)
                {
                    yield return ('M', k + 1, t.Trans[k][0], true);
                    yield return ('D', k + 1, t.Trans[k][2], false);
                }
                yield return ('I', k, t.Trans[k][1], true);
                break;
            case 'I':
                if (k < m) yield return ('M', k + 1, t.Trans[k][3], true);
                yield return ('I', k, t.Trans[k][4], true);
                break;
            case 'D':
                if (k < m)
                {
                    yield return ('M', k + 1, t.Trans[k][5], true);
                    yield return ('D', k + 1, t.Trans[k][6], false);
                }
                break;
        }
    }

    private static (double Fwd, double Vit) GlocalEnumerate(ToyHmm t, string seq)
    {
        int m = t.M, n = seq.Length;
        var dsq = new int[n + 1];
        for (int i = 1; i <= n; i++) dsq[i] = Res(seq[i - 1]);
        var paths = new List<double>();
        void Walk(char kind, int k, int i, double acc)
        {
            if ((kind == 'M' || kind == 'D') && k == m && i == n) paths.Add(acc); // M_m/D_m -> E (log 0)
            foreach (var (v, vk, l, emits) in GlocalSucc(t, kind, k))
            {
                if (double.IsNegativeInfinity(l)) continue;
                if (emits)
                {
                    if (i + 1 > n) continue;
                    double e = GlocalEmit(t, v, vk, dsq[i + 1]);
                    if (double.IsNegativeInfinity(e)) continue;
                    Walk(v, vk, i + 1, acc + l + e);
                }
                else Walk(v, vk, i, acc + l);
            }
        }
        Walk('S', 0, 0, 0.0);
        double fwd = NegInf, vit = NegInf;
        foreach (var p in paths) { fwd = LseTest(fwd, p); vit = Math.Max(vit, p); }
        return (fwd, vit);
    }

    public static IEnumerable<TestCaseData> GlocalCases()
    {
        // Empty sequence exercises ONLY the i=0 init: the all-delete path and NewRow boundary. This is
        // where the `m >= 1` delete guard (equivalent for m>=2), the `j <= m` D-chain bound, and the
        // NewRow `j <= m` fill (a missed index m leaks 0.0 instead of -inf) become observable.
        yield return new TestCaseData(1, "").SetName("Glocal_Hmm1_empty");
        yield return new TestCaseData(2, "").SetName("Glocal_Hmm2_empty");
        yield return new TestCaseData(3, "").SetName("Glocal_Hmm3_empty");
        yield return new TestCaseData(1, "A").SetName("Glocal_Hmm1_A");
        yield return new TestCaseData(1, "C").SetName("Glocal_Hmm1_C");
        yield return new TestCaseData(1, "AC").SetName("Glocal_Hmm1_AC");
        yield return new TestCaseData(2, "AC").SetName("Glocal_Hmm2_AC");
        yield return new TestCaseData(2, "C").SetName("Glocal_Hmm2_C");        // forces B->D1 entry
        yield return new TestCaseData(2, "AAC").SetName("Glocal_Hmm2_AAC");    // forces an insert
        yield return new TestCaseData(2, "ACAC").SetName("Glocal_Hmm2_ACAC");  // two inserts
        yield return new TestCaseData(3, "AC").SetName("Glocal_Hmm3_AC");      // a delete in the middle
        yield return new TestCaseData(3, "C").SetName("Glocal_Hmm3_C");        // two deletes
        yield return new TestCaseData(3, "ACA").SetName("Glocal_Hmm3_ACA");
        yield return new TestCaseData(3, "ACAC").SetName("Glocal_Hmm3_ACAC");
    }

    [TestCaseSource(nameof(GlocalCases))]
    public void GlocalForwardViterbi_MatchIndependentEnumeration(int nodes, string seq)
    {
        var t = nodes == 1 ? BuildHmm1() : nodes == 2 ? BuildHmm2() : BuildHmm3();
        var (fwd, vit) = GlocalEnumerate(t, seq);

        // Independent anchor: the 2-node "AC" optimum is the exact Durbin derivation 0.5187937934151676.
        if (nodes == 2 && seq == "AC")
            AssertCell(0.5187937934151676, vit, "Durbin anchor Viterbi(AC)");

        AssertCell(vit, t.Hmm.ViterbiScore(seq), $"ViterbiScore({seq})");
        AssertCell(fwd, t.Hmm.ForwardScore(seq), $"ForwardScore({seq})");
        Assert.That(fwd, Is.GreaterThanOrEqualTo(vit - 1e-9), "Forward >= Viterbi");
    }

    private static double LseTest(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double hi = Math.Max(a, b), lo = Math.Min(a, b);
        return hi + Math.Log(1.0 + Math.Exp(lo - hi));
    }
}
