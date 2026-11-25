namespace Seqeron.Genomics.Analysis;

public static partial class RnaSecondaryStructure
{
    #region Pseudoknot Prediction (pknotsRG canonical simple recursive pseudoknots, exact O(n^4) DP)

    // The pknotsRG grammar (Reeder & Giegerich 2004, BMC Bioinformatics 5:104, PMC514697) over the shared Turner 2004
    // model. Pseudoknot-specific parameters are the published pknotsRG ones (Energy.lhs, github.com/jensreeder/pknotsRG):
    //   pkinit = 900 (new pseudoknot: 9.0 kcal/mol), npp = 30 (every unpaired nucleotide of a pseudoknot loop: 0.3),
    //   pbp = 0 (base pair inside the knot: 0.0), pkml = +600 (a knot used as a branch of a multiloop / of another
    //   knot's loop: 6.0, Algebras.lhs energyalgext). Loop, stack and junction energies are NOT re-tabulated: every
    //   nested term is the shared Turner 2004 (ViennaRNA d2) TurnerModel term, see KnotLayer.
    private const double Pseudoknot_InitiationPenalty = 9.0;       // pkinit
    private const double Pseudoknot_UnpairedLoopPenalty = 0.3;     // npp, per unpaired nucleotide of a knot loop
    private const double Pseudoknot_ComponentPenalty = 6.0;        // pkml, knot as a multiloop / knot-loop component
    private const int Pseudoknot_MinHelixPairs = 2;                // canonical helices need >= 2 stacked pairs
    private const int Pseudoknot_MinSpan = 11;                     // 2*(2+2) + 1 (explicit) + 0 (middle) + 2 (explicit)

    /// <summary>
    /// A pseudoknot prediction: the folded sequence, the dot-bracket annotation (nested pairs <c>()</c>, crossing
    /// helices on the further bracket families of <see cref="ToDotBracket"/>), the full base-pair set and the
    /// total free energy.
    /// </summary>
    /// <param name="Sequence">The folded sequence (upper-cased; T read as U).</param>
    /// <param name="DotBracket">
    /// Dot-bracket of the pair set as rendered by <see cref="ToDotBracket"/> (ViennaRNA <c>vrna_db_from_ptable</c>
    /// families). Without a pseudoknot this is the plain MFE dot-bracket (single bracket family).
    /// </param>
    /// <param name="BasePairs">All base pairs as 0-based (5' &lt; 3') tuples, sorted by 5' position.</param>
    /// <param name="FreeEnergy">ΔG° (kcal/mol) of the returned structure under the pknotsRG-over-Turner-2004 model.</param>
    /// <param name="HasPseudoknot">True iff the returned structure contains a crossing (pseudoknotted) helix.</param>
    public readonly record struct PseudoknotStructure(
        string Sequence,
        string DotBracket,
        IReadOnlyList<(int Position1, int Position2)> BasePairs,
        double FreeEnergy,
        bool HasPseudoknot);

    /// <summary>
    /// Minimum-free-energy structure that may contain AT MOST ONE canonical simple recursive H-type pseudoknot
    /// (Reeder &amp; Giegerich 2004, class csr-PK), the knot standing in the exterior loop and its three loops folding
    /// pseudoknot-free. This is the exact optimum of that search space (the pknotsRG O(n⁴) dynamic programme with the
    /// knot count limited to one), not a heuristic: the rest of the sequence (before, after and inside the knot loops)
    /// folds by the shared Turner 2004 engine. The result is never worse than <see cref="CalculateMfeStructure"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Knot = a·u·b·v·a'·w·b' (5'→3'), helix a pairs a', b pairs b'. Canonization (the paper's rules): helices have
    /// equal strand lengths and no bulges (rule 1), a·a' is the maximal helix starting at its outer pair (rule 2),
    /// b·b' likewise but shortened where the maximal helices would compete for the bases of v (rule 3, boundary
    /// taken from b as in pknotsRG <c>getbeta</c>); each helix has ≥ 2 pairs; one nucleotide after a and two before b'
    /// stay unpaired (pknotsRG <c>front</c>/<c>back</c>), v may be empty.
    /// </para>
    /// <para>
    /// Energy: 9.0 (initiation) + Turner 2004 stacks of both helices (incl. the GGUC/CUGG three-stack) + 0.3 per
    /// unpaired nucleotide of u, v, w (the three bridging nucleotides included) + d2 mismatch/terminal-AU of the four
    /// helix ends + each stem component of a loop scored like an exterior-loop stem. See
    /// docs/algorithms/RnaStructure/Pseudoknot_Prediction.md for the full formula and the deviations from pknotsRG
    /// (Turner 1999 parameters, relaxed 1-nt-bulge helices, coaxial stacking).
    /// </para>
    /// <para>Time O(n⁴) (four free boundaries i, l, k, j), space O(n²).</para>
    /// </remarks>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <param name="temperature">Folding temperature in kelvin for the Turner terms (default 310.15 = 37 °C).</param>
    /// <returns>
    /// The best structure; for null/empty/too-short input, an empty pseudoknot-free structure
    /// (no pairs, all dots, ΔG = 0).
    /// </returns>
    public static PseudoknotStructure PredictStructurePseudoknot(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
        => FoldWithKnots(rnaSequence, minLoopSize, temperature, nestedKnots: false);

    /// <summary>
    /// Full pknotsRG grammar: MFE structure over ALL canonical simple recursive pseudoknots, including several knots
    /// side by side, knots as branches of multiloops and knots inside the loops u, v, w of other knots
    /// (Reeder &amp; Giegerich 2004: "the unpaired strands u, v, w … fold internally in an arbitrary way, including simple
    /// recursive pseudoknots"). Exact optimum of that space; never worse than <see cref="PredictStructurePseudoknot"/>
    /// or <see cref="CalculateMfeStructure"/>. A knot used as a multiloop branch or as a component of a knot loop pays
    /// the pknotsRG <c>pkml</c> penalty 6.0 (in place of the multiloop branch term).
    /// </summary>
    /// <remarks>Same energy model, canonization and complexity as <see cref="PredictStructurePseudoknot"/>.</remarks>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <param name="temperature">Folding temperature in kelvin for the Turner terms (default 310.15 = 37 °C).</param>
    public static PseudoknotStructure PredictStructurePseudoknotRecursive(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
        => FoldWithKnots(rnaSequence, minLoopSize, temperature, nestedKnots: true);

    private static PseudoknotStructure FoldWithKnots(string rnaSequence, int minLoopSize, double temperature, bool nestedKnots)
    {
        if (minLoopSize < MinHairpinLoop) minLoopSize = MinHairpinLoop;
        var model = ModelFor(temperature);

        string seq = rnaSequence is null ? "" : rnaSequence.ToUpperInvariant().Replace('T', 'U');
        int n = seq.Length;
        if (n < Pseudoknot_MinSpan)
        {
            var plain = CalculateMfeStructure(seq, minLoopSize, temperature);
            return new PseudoknotStructure(plain.Sequence, plain.DotBracket, plain.BasePairs, plain.FreeEnergy, HasPseudoknot: false);
        }

        byte[] s = TurnerModel.Encode(seq);
        double[] vBuf = new double[n * n];
        double[] wmBuf = new double[n * n];
        double[] w = new double[n];
        var knots = new KnotLayer(model, s, n, minLoopSize, nestedKnots);

        FillDp(model, s, n, minLoopSize, vBuf, wmBuf, w, knots);
        var pairs = TracebackMfe(model, s, n, minLoopSize, vBuf, wmBuf, w, knots);
        pairs.Sort((a, b) => a.Position1.CompareTo(b.Position1));

        double energy = RoundEnergy(model, nestedKnots ? w[n - 1] : knots.AtMostOneKnot[n - 1]);
        return new PseudoknotStructure(seq, ToDotBracket(n, pairs), pairs, energy, HasPseudoknot: knots.KnotsTraced > 0);
    }

    /// <summary>
    /// The pseudoknot layer of the Zuker–Stiegler fill (<see cref="FillDp"/>): per interval [i,j] the best canonical
    /// knot <c>K</c> spanning exactly [i,j] and the best fold <c>C</c> of a knot LOOP segment [i,j] (a sequence of stem
    /// components and unpaired nucleotides). Every nested energy is the shared <see cref="TurnerModel"/> term.
    /// </summary>
    private sealed class KnotLayer
    {
        private const double INF = TurnerModel.Inf;
        private readonly TurnerModel _m;
        private readonly byte[] _s;
        private readonly int _n;
        private readonly int[] _maxPairs;       // consecutive pairs (p,q),(p+1,q-1),… (rule 2: maximal helix)
        private readonly double[] _helixFull;   // Turner energy of that maximal helix

        /// <summary>True: knots may sit in multiloops and in knot loops (full grammar); false: one knot, exterior loop.</summary>
        internal bool NestedKnots { get; }
        /// <summary>K(i,j): best knot spanning exactly [i,j] (INF if none), without the multiloop/loop component penalty.</summary>
        internal double[] K { get; }
        /// <summary>C(i,j): best knot-loop segment fold of [i,j]; C of an empty interval is 0.</summary>
        internal double[] C { get; }
        /// <summary>Single-knot mode: best prefix [0..j] fold with at most one knot.</summary>
        internal double[] AtMostOneKnot { get; }
        internal int KnotsTraced { get; private set; }

        internal KnotLayer(TurnerModel m, byte[] s, int n, int minLoop, bool nestedKnots)
        {
            _m = m; _s = s; _n = n; NestedKnots = nestedKnots;
            K = new double[n * n];
            C = new double[n * n];
            AtMostOneKnot = new double[n];
            Array.Fill(K, INF);

            _maxPairs = new int[n * n];
            for (int span = 2; span <= n; span++)
                for (int p = 0; p + span - 1 < n; p++)
                {
                    int q = p + span - 1;
                    if (!TurnerModel.CanPair(s[p], s[q])) continue;
                    _maxPairs[p * n + q] = 1 + (p + 1 < q - 1 ? _maxPairs[(p + 1) * n + (q - 1)] : 0);
                }
            _helixFull = new double[n * n];
            for (int p = 0; p < n; p++)
                for (int q = p + 1; q < n; q++)
                    if (_maxPairs[p * n + q] >= Pseudoknot_MinHelixPairs)
                        _helixFull[p * n + q] = HelixEnergy(p, q, _maxPairs[p * n + q]);

            // Segments too short for any pair are all-unpaired (FillDp starts at span minLoop+2).
            for (int span = 1; span < minLoop + 2 && span <= n; span++)
                for (int i = 0; i + span - 1 < n; i++)
                    C[i * n + i + span - 1] = Pseudoknot_UnpairedLoopPenalty * span;
        }

        private double Cget(int i, int j) => j < i ? 0.0 : C[i * _n + j];

        // Turner 2004 energy of the bulge-free helix of `len` pairs starting at (p,q): consecutive stacks, with the
        // special GGUC/CUGG three-stack as the alternative FillDp offers for V(i,j).
        private double HelixEnergy(int p, int q, int len)
        {
            Span<double> e = len <= 64 ? stackalloc double[len + 1] : new double[len + 1];
            e[len] = 0; e[len - 1] = 0;
            for (int t = len - 2; t >= 0; t--)
            {
                double best = _m.InteriorTerm(_s, p + t, q - t, p + t + 1, q - t - 1) + e[t + 1];
                if (t + 3 <= len - 1 && TurnerModel.IsSpecialGgucCuggStack(_s, p + t, q - t))
                    best = Math.Min(best, SpecialGGUC_CUGG_3Stack + e[t + 3]);
                e[t] = best;
            }
            return e[0];
        }

        /// <summary>Fills K(i,j) and C(i,j) once V(i,j) is known (spans are visited in increasing order).</summary>
        internal void FillCell(int i, int j, double[] vBuf)
        {
            int n = _n, ij = i * n + j;
            if (j - i + 1 >= Pseudoknot_MinSpan) K[ij] = ComputeKnot(i, j, out _, out _, out _);

            double best = Pseudoknot_UnpairedLoopPenalty + Cget(i + 1, j);
            for (int q = i; q <= j; q++)
            {
                double v = vBuf[i * n + q];
                if (v < INF)
                    best = Math.Min(best, v + _m.ExteriorStem(_s, i, q) + Cget(q + 1, j));
                if (NestedKnots && K[i * n + q] < INF)
                    best = Math.Min(best, K[i * n + q] + Pseudoknot_ComponentPenalty + Cget(q + 1, j));
            }
            C[ij] = best;
        }

        /// <summary>min over p of W(p−1) + K(p,j): a knot closing the exterior prefix [0..j].</summary>
        internal double BestExteriorKnot(int j, double[] w)
        {
            double best = INF;
            for (int p = 0; p + Pseudoknot_MinSpan - 1 <= j; p++)
            {
                double k = K[p * _n + j];
                if (k < INF) best = Math.Min(best, (p > 0 ? w[p - 1] : 0.0) + k);
            }
            return best;
        }

        /// <summary>Single-knot mode: the prefix fold with at most one knot, AtMostOneKnot[j].</summary>
        internal void FillExterior(int j, double[] vBuf, double[] w)
        {
            if (NestedKnots) return;
            int n = _n;
            double best = w[j];
            if (j > 0) best = Math.Min(best, AtMostOneKnot[j - 1]);
            for (int p = 0; p <= j; p++)
            {
                double v = vBuf[p * n + j];
                if (v < INF)
                    best = Math.Min(best, (p > 0 ? AtMostOneKnot[p - 1] : 0.0) + v + _m.ExteriorStem(_s, p, j));
            }
            best = Math.Min(best, BestExteriorKnot(j, w));
            AtMostOneKnot[j] = best;
        }

        // The knot spanning [i,j]: a = [i..i+α−1] with a' = [k−α+1..k] (maximal helix of (i,k)), b = [l..l+β−1] with
        // b' = [j−β+1..j] (maximal helix of (l,j), shortened so that b stays left of a' and w keeps ≥ 2 nt).
        // Loops: u' = [i+α+1..l−1] after the explicit unpaired nt i+α, v = [l+β..k−α], w' = [k+1..j−β−2] before the two
        // explicit unpaired nts j−β−1, j−β.
        private double ComputeKnot(int i, int j, out int bestK, out int bestL, out int bestBeta)
        {
            int n = _n;
            double best = INF;
            bestK = bestL = bestBeta = -1;
            double constant = Pseudoknot_InitiationPenalty + 3 * Pseudoknot_UnpairedLoopPenalty;

            for (int k = i + 5; k <= j - 4; k++)
            {
                int alpha = _maxPairs[i * n + k];
                if (alpha < Pseudoknot_MinHelixPairs) continue;
                int a5Last = i + alpha - 1, a3First = k - alpha + 1;
                double eA = constant + _helixFull[i * n + k] + _m.ExteriorStem(_s, i, k)
                            + _m.InnerStemEnd(_s, a5Last, a3First);

                for (int l = i + alpha + 1; l <= k - alpha - 1; l++)
                {
                    int beta0 = _maxPairs[l * n + j];
                    if (beta0 < Pseudoknot_MinHelixPairs) continue;
                    int beta = Math.Min(beta0, Math.Min(k - alpha - l + 1, j - k - 2));
                    if (beta < Pseudoknot_MinHelixPairs) continue;

                    double e = eA + (beta == beta0 ? _helixFull[l * n + j] : HelixEnergy(l, j, beta))
                               + _m.ExteriorStem(_s, l, j) + _m.InnerStemEnd(_s, l + beta - 1, j - beta + 1)
                               + Cget(i + alpha + 1, l - 1) + Cget(l + beta, k - alpha) + Cget(k + 1, j - beta - 2);
                    if (e < best) { best = e; bestK = k; bestL = l; bestBeta = beta; }
                }
            }
            return best;
        }

        /// <summary>Traceback of the work items 'X' (prefix with ≤ 1 knot), 'K' (knot) and 'C' (knot-loop segment).</summary>
        internal void Trace(char kind, int i, int j, Stack<(char Kind, int I, int J)> stack,
            List<(int Position1, int Position2)> pairs, double[] vBuf, double[] w)
        {
            const double EPS = Dp_TraceEps;
            int n = _n;
            if (kind == 'X')
            {
                if (j < 0) return;
                double target = AtMostOneKnot[j];
                if (Math.Abs(target - w[j]) <= EPS) { stack.Push(('W', 0, j)); return; }
                if (j > 0 && Math.Abs(target - AtMostOneKnot[j - 1]) <= EPS) { stack.Push(('X', 0, j - 1)); return; }
                for (int p = 0; p <= j; p++)
                {
                    double v = vBuf[p * n + j];
                    if (v < INF && Math.Abs(target - ((p > 0 ? AtMostOneKnot[p - 1] : 0.0) + v + _m.ExteriorStem(_s, p, j))) <= EPS)
                    {
                        if (p > 0) stack.Push(('X', 0, p - 1));
                        stack.Push(('V', p, j));
                        return;
                    }
                }
                for (int p = 0; p + Pseudoknot_MinSpan - 1 <= j; p++)
                {
                    double k = K[p * n + j];
                    if (k < INF && Math.Abs(target - ((p > 0 ? w[p - 1] : 0.0) + k)) <= EPS)
                    {
                        if (p > 0) stack.Push(('W', 0, p - 1));
                        stack.Push(('K', p, j));
                        return;
                    }
                }
                return;
            }

            if (kind == 'C')
            {
                if (j < i) return;
                double target = C[i * n + j];
                if (Math.Abs(target - (Pseudoknot_UnpairedLoopPenalty + Cget(i + 1, j))) <= EPS)
                {
                    PushLoop(stack, i + 1, j);
                    return;
                }
                for (int q = i; q <= j; q++)
                {
                    double v = vBuf[i * n + q];
                    if (v < INF && Math.Abs(target - (v + _m.ExteriorStem(_s, i, q) + Cget(q + 1, j))) <= EPS)
                    {
                        stack.Push(('V', i, q));
                        PushLoop(stack, q + 1, j);
                        return;
                    }
                    if (NestedKnots && K[i * n + q] < INF
                        && Math.Abs(target - (K[i * n + q] + Pseudoknot_ComponentPenalty + Cget(q + 1, j))) <= EPS)
                    {
                        stack.Push(('K', i, q));
                        PushLoop(stack, q + 1, j);
                        return;
                    }
                }
                return;
            }

            // kind == 'K'
            {
                ComputeKnot(i, j, out int k, out int l, out int beta);
                int alpha = _maxPairs[i * n + k];
                for (int t = 0; t < alpha; t++) pairs.Add((i + t, k - t));
                for (int t = 0; t < beta; t++) pairs.Add((l + t, j - t));
                PushLoop(stack, i + alpha + 1, l - 1);
                PushLoop(stack, l + beta, k - alpha);
                PushLoop(stack, k + 1, j - beta - 2);
                KnotsTraced++;
            }
        }

        private static void PushLoop(Stack<(char Kind, int I, int J)> stack, int i, int j)
        {
            if (j >= i) stack.Push(('C', i, j));
        }
    }

    #endregion
}
