namespace Seqeron.Genomics.Analysis;

public static partial class RnaSecondaryStructure
{
    #region Partition function (McCaskill 1990 on the Turner 2004 model)

    /// <summary>
    /// Result of the McCaskill (1990) equilibrium partition-function computation on the Turner 2004
    /// nearest-neighbour model — the counterpart of ViennaRNA 2.7.2 <c>fc.pf()</c> + <c>fc.bpp()</c> +
    /// <c>fc.mfe()</c> + <c>fc.pr_structure()</c> + <c>fc.mean_bp_distance()</c> + <c>fc.centroid()</c>.
    /// </summary>
    /// <param name="PartitionFunction">
    /// Z = Σ_S exp(−E(S)/RT) over every pseudoknot-free structure S of the sequence (the open chain, E = 0, is one term,
    /// so Z ≥ 1). It is <c>+∞</c> when it exceeds the range of a double (≈ e^709; very long or very cold ensembles) —
    /// use <see cref="LogPartitionFunction"/> / <see cref="EnsembleFreeEnergy"/>, which are always finite.
    /// </param>
    /// <param name="BasePairProbabilities">
    /// P[(i,j)] = (Σ_{S ∋ (i,j)} exp(−E(S)/RT)) / Z, the equilibrium probability that 0-based positions i &lt; j pair
    /// (McCaskill's outside recursion; ViennaRNA <c>bpp()</c>). Only pairs with P &gt; 0 are listed.
    /// </param>
    /// <param name="EnsembleFreeEnergy">
    /// G = −RT·ln Z (kcal/mol), ViennaRNA's <c>pf()</c> energy; G ≤ MFE.
    /// </param>
    /// <param name="LogPartitionFunction">ln Z (finite).</param>
    /// <param name="UnpairedProbabilities">
    /// p_unpaired[i] = 1 − Σ_j P[(i,j)] over both i&lt;j and j&lt;i (one entry per position), so
    /// p_unpaired[i] + Σ_j P(i,j) = 1.
    /// </param>
    /// <param name="MfeStructure">Dot-bracket of the minimum-free-energy structure (<see cref="CalculateMfeStructure"/>).</param>
    /// <param name="MfeFreeEnergy">Its free energy (kcal/mol, ViennaRNA <c>mfe()</c>).</param>
    /// <param name="MfeProbability">
    /// Boltzmann probability of the MFE structure within the ensemble, exp((G − E_mfe)/RT) (ViennaRNA <c>pr_structure(mfe_struct)</c>).
    /// </param>
    /// <param name="EnsembleDiversity">
    /// Mean base-pair distance of the ensemble ⟨d⟩ = 2·Σ_{i&lt;j} P(i,j)(1 − P(i,j)) (ViennaRNA <c>mean_bp_distance()</c>).
    /// </param>
    /// <param name="CentroidStructure">
    /// Centroid structure: every pair with P &gt; 0.5 (minimal mean base-pair distance to the ensemble; ViennaRNA <c>centroid()</c>).
    /// </param>
    /// <param name="Temperature">Absolute temperature (K) of the computation.</param>
    public readonly record struct PartitionFunctionResult(
        double PartitionFunction,
        IReadOnlyDictionary<(int I, int J), double> BasePairProbabilities,
        double EnsembleFreeEnergy,
        double LogPartitionFunction,
        IReadOnlyList<double> UnpairedProbabilities,
        string MfeStructure,
        double MfeFreeEnergy,
        double MfeProbability,
        double EnsembleDiversity,
        string CentroidStructure,
        double Temperature);

    /// <summary>
    /// Result of the McCaskill (1990) partition function evaluated on the SAME Turner-2004
    /// nearest-neighbour energy model used by <see cref="CalculateMinimumFreeEnergy"/>:
    /// the Boltzmann-weighted sum over the full ensemble of pseudoknot-free secondary
    /// structures, the equilibrium base-pair probabilities, the per-base unpaired
    /// probabilities, and the ensemble free energy.
    /// </summary>
    /// <param name="PartitionFunction">
    /// Z = Σ_S exp(−E(S)/RT) over every admissible structure S, E(S) scored with the Turner-2004
    /// model (McCaskill JS 1990, Biopolymers 29:1105; Lorenz et al. 2011, ViennaRNA). Z ≥ 1
    /// (the open chain, E = 0, is always one term).
    /// </param>
    /// <param name="BasePairProbabilities">
    /// P[(i,j)] = (Σ_{S ∋ (i,j)} exp(−E(S)/RT)) / Z, the equilibrium probability that 0-based
    /// positions i &lt; j pair, via the inside/outside recursion. Only pairs with P &gt; 0 are listed.
    /// </param>
    /// <param name="UnpairedProbabilities">
    /// p_unpaired[i] = 1 − Σ_j P[(i,j)] (summed over both i&lt;j and j&lt;i), the equilibrium
    /// probability that 0-based position i is unpaired (Lorenz et al. 2011). One entry per position.
    /// </param>
    /// <param name="EnsembleFreeEnergy">
    /// ΔG_ensemble = −RT·ln Z (kcal/mol). By the laws of statistical mechanics this is a lower
    /// bound on the MFE returned by <see cref="CalculateMinimumFreeEnergy"/> for the same input.
    /// </param>
    public readonly record struct UnpairedProbabilityResult(
        double PartitionFunction,
        IReadOnlyDictionary<(int I, int J), double> BasePairProbabilities,
        IReadOnlyList<double> UnpairedProbabilities,
        double EnsembleFreeEnergy);

    // ViennaRNA model default sfact (model.c): the partition-function scale is chosen so that the largest Boltzmann
    // factor, exp(−MFE/RT), is divided out by sfact·MFE spread evenly over the nucleotides.
    private const double PfScaleFactor = 1.07;

    // Scaled inside matrices (ViennaRNA pf_scale): every quantity for a stretch of k nucleotides is stored multiplied
    // by us^k (us = 1/pf_scale ≤ 1), which keeps all intermediate values of long sequences inside the double range.
    private sealed class PfInside
    {
        public required int N;
        public required double Us;
        public required double[] Usp;   // Usp[k] = Us^k
        public required double[] V;     // V̂(i,j)  pair (i,j) closes the structure, i*n+j
        public required double[] Wm;    // ŴM(i,j) multiloop segment with ≥ 1 branch
        public required double[] Wm1;   // ŴM1(i,j) multiloop segment with exactly 1 branch starting at i
        public required double[] Wl;    // Wl[j] = Ŵ of the prefix of length j (positions 0..j−1), Wl[0] = 1
        public double ZScaled => Wl[N];
    }

    /// <summary>
    /// Computes the McCaskill (1990) equilibrium partition function Z of the full ensemble of pseudoknot-free
    /// secondary structures on the Turner 2004 nearest-neighbour model — the SAME energy model as
    /// <see cref="CalculateMinimumFreeEnergy"/> (ViennaRNA 2.7.2 <c>RNAfold -p</c>, d2 dangles) — together with the
    /// equilibrium base-pair probabilities, unpaired probabilities, ensemble free energy, MFE structure and its
    /// probability, ensemble diversity and centroid structure.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Every place the MFE recursion takes a <c>min</c> and ADDS a loop energy ΔG, the partition function takes a
    /// <c>sum</c> and MULTIPLIES by exp(−ΔG/RT) of the same <see cref="TurnerModel"/> loop energy. The decomposition is
    /// unambiguous (every structure has exactly one derivation; multiloop branches go through a single-branch matrix as
    /// in ViennaRNA's <c>qm1</c>), so Z = Σ_S exp(−E(S)/RT) over exactly the structures the MFE folder explores. Base-pair
    /// probabilities use McCaskill's inside/outside recursion on the same grammar (O(n³) time, O(n²) memory).
    /// </para>
    /// <para>
    /// Numerical stability: all matrices are scaled per nucleotide (ViennaRNA <c>pf_scale</c> from the MFE,
    /// <c>sfact</c> = 1.07), so n = 2000 and beyond neither overflows nor underflows; Z itself is reported as a double
    /// (+∞ beyond ≈ e^709) next to its logarithm. The partition function agrees with ViennaRNA <c>fc.pf()</c> up to the
    /// documented parameter rounding of ViennaRNA's default file (NNDB values are used here).
    /// </para>
    /// <para>
    /// Temperature: <paramref name="temperature"/> rescales every Turner 2004 parameter with its enthalpy exactly like
    /// ViennaRNA (<c>md.temperature</c>, see <see cref="TurnerModel"/>), and RT = 1.98717e-3·T kcal/mol — not just RT.
    /// </para>
    /// </remarks>
    /// <param name="sequence">RNA (or DNA; T read as U) sequence; case-insensitive. Characters other than A/C/G/U/T (N, IUPAC, gaps) never pair.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not a positive finite number.</exception>
    public static PartitionFunctionResult CalculatePartitionFunction(
        string sequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var model = ModelFor(temperature);
        if (minLoopSize < 3) minLoopSize = 3;

        string seq = sequence.ToUpperInvariant().Replace('T', 'U');
        int n = seq.Length;
        var mfe = CalculateMfeStructure(seq, minLoopSize, temperature);
        var noPairs = new Dictionary<(int, int), double>();

        if (n < minLoopSize + 2)
        {
            // No pair can form: only the open chain (Z = 1), every base unpaired.
            var allUnpaired = new double[n];
            Array.Fill(allUnpaired, 1.0);
            return new PartitionFunctionResult(1.0, noPairs, 0.0, 0.0, allUnpaired,
                mfe.DotBracket, mfe.FreeEnergy, 1.0, 0.0, new string('.', n), temperature);
        }

        byte[] s = TurnerModel.Encode(seq);
        double rt = model.RT;

        // Scale from the MFE; a pathological ensemble that still leaves the double range is retried with
        // the unscaled-energy estimate before giving up (never silently NaN).
        PfInside? inside = null;
        double logUs = 0;
        foreach (double sfact in new[] { PfScaleFactor, 1.0, 0.9, 1.2 })
        {
            logUs = sfact * mfe.FreeEnergy / (rt * n);
            inside = FillPartitionInside(model, s, n, minLoopSize, Math.Exp(logUs));
            if (double.IsFinite(inside.ZScaled) && inside.ZScaled > 0) break;
            inside = null;
        }
        if (inside is null)
            return GroundStateEnsemble(mfe, n, rt, temperature);

        double logZ = Math.Log(inside.ZScaled) - n * logUs;
        double z = inside.Usp[n] > 1e-290 && double.IsFinite(inside.ZScaled / inside.Usp[n]) ? inside.ZScaled / inside.Usp[n] : Math.Exp(logZ);
        double g = -rt * logZ + 0.0; // + 0.0 turns −0 into +0 for Z = 1

        double[] prob = ComputePairProbabilities(model, s, n, minLoopSize, inside);

        var bpp = new Dictionary<(int I, int J), double>();
        var paired = new double[n];
        double diversity = 0.0;
        var centroid = new char[n];
        Array.Fill(centroid, '.');
        for (int i = 0; i < n; i++)
            for (int j = i + minLoopSize + 1; j < n; j++)
            {
                double p = prob[i * n + j];
                if (!(p > 0)) continue;
                if (p > 1) p = 1; // round-off guard
                bpp[(i, j)] = p;
                paired[i] += p;
                paired[j] += p;
                diversity += p * (1 - p);
                if (p > 0.5) { centroid[i] = '('; centroid[j] = ')'; }
            }

        var unpaired = new double[n];
        for (int i = 0; i < n; i++) unpaired[i] = Math.Clamp(1.0 - paired[i], 0.0, 1.0);

        double mfeProbability = Math.Min(1.0, Math.Exp((g - mfe.FreeEnergy) / rt));

        return new PartitionFunctionResult(z, bpp, g, logZ, unpaired, mfe.DotBracket, mfe.FreeEnergy,
            mfeProbability, 2.0 * diversity, new string(centroid), temperature);
    }

    // T → 0 limit: when a single loop's Boltzmann weight already exceeds the double range (RT below ≈ 0.01 kcal/mol,
    // i.e. a few Kelvin) the ensemble is the ground state — G = MFE, ln Z = −MFE/RT, the MFE structure has probability 1
    // (ties between degenerate ground states are not resolved).
    private static PartitionFunctionResult GroundStateEnsemble(MfeStructure mfe, int n, double rt, double temperature)
    {
        double logZ = -mfe.FreeEnergy / rt;
        double z = logZ < 709.0 ? Math.Exp(logZ) : double.PositiveInfinity;
        var bpp = new Dictionary<(int I, int J), double>();
        var unpaired = new double[n];
        Array.Fill(unpaired, 1.0);
        foreach (var (a, b) in mfe.BasePairs)
        {
            bpp[(a, b)] = 1.0;
            unpaired[a] = 0.0;
            unpaired[b] = 0.0;
        }
        return new PartitionFunctionResult(z, bpp, mfe.FreeEnergy, logZ, unpaired, mfe.DotBracket,
            mfe.FreeEnergy, 1.0, 0.0, mfe.DotBracket, temperature);
    }

    /// <summary>
    /// Calculates the Boltzmann probability of a single structure within the ensemble:
    /// p(S) = exp(−E_S/RT) / exp(−E_ensemble/RT) = exp((E_ensemble − E_S)/RT), where E_ensemble = −RT·ln Z is the
    /// ensemble free energy (ViennaRNA <c>pr_structure</c> / <c>pr_energy</c>: <c>exp((dG − e)/kT)</c>), with
    /// RT = 1.98717e-3·T kcal/mol (ViennaRNA GASCONST). Evaluated in the exponent so that large |E| cannot overflow.
    /// </summary>
    /// <remarks>
    /// Source: McCaskill JS (1990) Biopolymers 29:1105-1119, p(P|S) = Z⁻¹·exp(−βE(P)); ViennaRNA 2.7.2
    /// <c>probabilities/probs_structures.c</c>.
    /// </remarks>
    /// <param name="structureEnergy">Free energy of the structure (kcal/mol).</param>
    /// <param name="ensembleEnergy">Ensemble free energy −RT·ln Z (kcal/mol), e.g. <see cref="PartitionFunctionResult.EnsembleFreeEnergy"/>.</param>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not a positive finite number.</exception>
    public static double CalculateStructureProbability(
        double structureEnergy, double ensembleEnergy, double temperature = DefaultTemperatureKelvin)
    {
        if (!(temperature > 0) || double.IsInfinity(temperature))
            throw new ArgumentOutOfRangeException(nameof(temperature), "Temperature must be a positive finite number (Kelvin).");

        double rt = TurnerModel.GasConstantCalPerMolK * temperature / 1000.0; // kcal/mol
        return Math.Exp((ensembleEnergy - structureEnergy) / rt);
    }

    /// <summary>
    /// The McCaskill (1990) partition function on the Turner-2004 model reduced to the quantities of
    /// <see cref="UnpairedProbabilityResult"/>; identical to the corresponding fields of
    /// <see cref="CalculatePartitionFunction"/> (same engine).
    /// </summary>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not a positive finite number.</exception>
    public static UnpairedProbabilityResult CalculateUnpairedProbabilities(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        var r = CalculatePartitionFunction(rnaSequence ?? string.Empty, minLoopSize, temperature);
        return new UnpairedProbabilityResult(r.PartitionFunction, r.BasePairProbabilities, r.UnpairedProbabilities, r.EnsembleFreeEnergy);
    }

    /// <summary>
    /// TargetScan/Agarwal-2015 structural accessibility (SA): the equilibrium probability that a
    /// contiguous window of <paramref name="windowLength"/> nucleotides ending at
    /// <paramref name="windowEnd"/> (0-based, inclusive) is entirely unpaired, computed EXACTLY
    /// from the Turner-2004 McCaskill partition function as Z_open(window) / Z, where
    /// Z_open(window) is the partition function over structures in which no base pair is incident
    /// to any position of the window (every base in the window unpaired).
    /// </summary>
    /// <remarks>
    /// RNAplfold (<c>RNAplfold -L 40 -W 80 -u 20</c> in <c>targetscan_70_context_scores.pl</c>,
    /// <c>runRNAplfold_all_UTRs</c>) reports, in row <c>i</c> column <c>L</c> of the <c>_lunp</c>
    /// file, the probability that the length-<c>L</c> stretch ENDING at position <c>i</c> is
    /// unpaired (RNAplfold man page; Bernhart et al. 2006). TargetScan's <c>getSA_contribution</c>
    /// reads column index 13 of that row (the 14th unpaired-probability column → <c>L = 14</c>)
    /// at the row 7 nt downstream of the seed-match start, i.e. the unpaired probability of the
    /// 14-nt region centred on the match to miRNA nucleotides 7–8 (Agarwal et al. 2015, eLife
    /// 4:e05005, Fig 4A: "the log10 value of the unpaired probability for a 14-nt region centered
    /// on the match to miRNA nucleotides 7 and 8"). The joint region-unpaired probability is the
    /// exact McCaskill ensemble quantity Z_open/Z.
    /// </remarks>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="windowEnd">0-based inclusive index of the last nucleotide of the unpaired window.</param>
    /// <param name="windowLength">Length of the unpaired window (TargetScan uses 14).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3).</param>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <returns>P(window entirely unpaired) ∈ [0,1]; 1.0 when no pair can form at all.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> ≤ 0, or the window does not fit in the sequence.</exception>
    public static double CalculateRegionUnpairedProbability(
        string rnaSequence, int windowEnd, int windowLength,
        int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        var model = ModelFor(temperature);
        if (windowLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(windowLength), "Window length must be positive.");
        if (minLoopSize < 3) minLoopSize = 3;

        string seq = string.IsNullOrEmpty(rnaSequence)
            ? string.Empty
            : rnaSequence.ToUpperInvariant().Replace('T', 'U');
        int n = seq.Length;

        int windowStart = windowEnd - windowLength + 1;
        if (windowStart < 0 || windowEnd >= n)
            throw new ArgumentOutOfRangeException(nameof(windowEnd),
                "The unpaired window must lie entirely within the sequence.");

        if (n < minLoopSize + 2)
            return 1.0; // no pair can form anywhere → the region is unpaired with probability 1.

        byte[] s = TurnerModel.Encode(seq);

        // Both partition functions use the same scale, so it cancels in the ratio. ViennaRNA's default scale for an
        // unknown MFE is −185 + (T−37)·7.27 cal/(mol·nt) (params.c, vrna_exp_params_rescale), which suffices here.
        double ePerNt = (-185.0 + (temperature - TurnerModel.ReferenceTemperatureKelvin) * 7.27) / 1000.0;
        double us = Math.Min(1.0, Math.Exp(PfScaleFactor * ePerNt / model.RT));

        var all = FillPartitionInside(model, s, n, minLoopSize, us);
        double z = all.ZScaled;
        if (!double.IsFinite(z) || z <= 0)
            return 0.0;

        // Z_open: the same partition function with EVERY base in the window forbidden from pairing.
        var open = FillPartitionInside(model, s, n, minLoopSize, us, windowStart, windowEnd);
        double pu = open.ZScaled / z;
        if (pu < 0) pu = 0;
        if (pu > 1) pu = 1; // round-off guard; Z_open ≤ Z by construction.
        return pu;
    }

    /// <summary>
    /// Fills the scaled inside matrices V̂/ŴM/ŴM1/Ŵ — the Boltzmann-domain mirror of <see cref="FillDp"/>. With
    /// <paramref name="forbiddenStart"/> ≤ <paramref name="forbiddenEnd"/> every pair incident to that window is excluded.
    /// </summary>
    /// <remarks>
    /// Unambiguous decomposition of the SAME energy model as FillDp (d2 dangles, linear multiloop, ML_base = 0 ⇒ an
    /// unpaired multiloop nucleotide has weight 1), each term carrying the scale us per nucleotide it spans:
    /// <code>
    ///   ŴM1(i,j) = ŴM1(i,j−1)·us + V̂(i,j)·e^{−Mstem(i,j)/RT}                (one branch starting at i)
    ///   ŴM(i,j)  = Σ_k [ us^(k−i) + ŴM(i,k−1) ]·ŴM1(k,j)                      (≥ 1 branch)
    ///   V̂(i,j)  = e^{−H}·us^len + Σ e^{−I}·V̂(p,q)·us^(unpaired+2) + GGUC corr. + e^{−Mclose}·us²·Σ_k ŴM(i+1,k−1)·ŴM1(k,j−1)
    ///   Ŵ(j)     = Ŵ(j−1)·us + Σ_i Ŵ(i−1)·V̂(i,j)·e^{−Ext(i,j)/RT}
    /// </code>
    /// </remarks>
    private static PfInside FillPartitionInside(
        TurnerModel m, byte[] s, int n, int minLoopSize, double us,
        int forbiddenStart = -1, int forbiddenEnd = -1)
    {
        const int MAXLOOP = Dp_MAXLOOP;
        double rt = m.RT;
        bool Forbidden(int p) => forbiddenStart <= forbiddenEnd && p >= forbiddenStart && p <= forbiddenEnd;

        var usp = new double[n + 3];
        usp[0] = 1.0;
        for (int k = 1; k < usp.Length; k++) usp[k] = usp[k - 1] * us;

        var v = new double[n * n];
        var wm = new double[n * n];
        var wm1 = new double[n * n];
        double ggucWeight = Math.Exp(-SpecialGGUC_CUGG_3Stack / rt);

        for (int span = minLoopSize + 2; span <= n; span++)
        {
            for (int i = 0; i <= n - span; i++)
            {
                int j = i + span - 1;
                int ij = i * n + j;

                if (TurnerModel.CanPair(s[i], s[j]) && !Forbidden(i) && !Forbidden(j))
                {
                    // Hairpin loop: all loop bases [i+1..j−1] unpaired.
                    double acc = Math.Exp(-m.HairpinTerm(s, i, j) / rt) * usp[span];

                    // Stack / bulge / interior loop with inner pair (ip,jp).
                    for (int ip = i + 1; ip < j - 1 && ip - i - 1 <= MAXLOOP; ip++)
                    {
                        int jpMin = Math.Max(ip + 1, j - (MAXLOOP - (ip - i - 1)) - 1);
                        for (int jp = j - 1; jp >= jpMin; jp--)
                        {
                            double inner = v[ip * n + jp];
                            if (inner <= 0) continue;
                            double e = m.InteriorTerm(s, i, j, ip, jp);
                            if (e >= Dp_INF) continue;
                            acc += Math.Exp(-e / rt) * inner * usp[(ip - i) + (j - jp)];
                        }
                    }

                    // Special GGUC/CUGG three-stack: replaces the three individual stacks (as FillDp does).
                    if (TurnerModel.IsSpecialGgucCuggStack(s, i, j) && v[(i + 1) * n + j - 1] > 0 && v[(i + 2) * n + j - 2] > 0)
                    {
                        double inner3 = v[(i + 3) * n + j - 3];
                        if (inner3 > 0)
                            acc += (ggucWeight - Math.Exp(-GgucSplitEnergy(m, s, i, j) / rt)) * inner3 * usp[6];
                    }

                    // Multiloop closed by (i,j): ≥ 2 branches inside.
                    double ml = 0.0;
                    for (int k = i + 2; k < j - 1; k++)
                    {
                        double right = wm1[k * n + j - 1];
                        if (right <= 0) continue;
                        double left = wm[(i + 1) * n + k - 1];
                        if (left > 0) ml += left * right;
                    }
                    if (ml > 0)
                        acc += Math.Exp(-m.MultiloopClosing(s, i, j) / rt) * ml * usp[2];

                    v[ij] = acc;
                }

                // ŴM1 / ŴM: multiloop parts, only meaningful for 0 < i, j < n−1 (neighbours exist).
                if (i > 0 && j < n - 1)
                {
                    double w1 = j - 1 >= i + minLoopSize + 1 ? wm1[ij - 1] * us : 0.0;
                    if (v[ij] > 0) w1 += v[ij] * Math.Exp(-m.MultiloopStem(s, i, j) / rt);
                    wm1[ij] = w1;

                    double w = 0.0;
                    for (int k = i; k <= j - minLoopSize - 1; k++)
                    {
                        double branch = wm1[k * n + j];
                        if (branch <= 0) continue;
                        double before = usp[k - i]; // Dp_ML_base = 0 ⇒ the unpaired prefix has weight 1 (× scale)
                        if (k > i) before += wm[i * n + k - 1];
                        w += before * branch;
                    }
                    wm[ij] = w;
                }
            }
        }

        // Ŵ(j): exterior-loop partition function of the prefix, Wl[j] for positions 0..j−1.
        var wl = new double[n + 1];
        wl[0] = 1.0;
        for (int j = 0; j < n; j++)
        {
            double wj = wl[j] * us; // j unpaired externally
            for (int i = 0; i + minLoopSize + 1 <= j; i++)
            {
                double vij = v[i * n + j];
                if (vij <= 0) continue;
                wj += wl[i] * vij * Math.Exp(-m.ExteriorStem(s, i, j) / rt);
            }
            wl[j + 1] = wj;
        }

        return new PfInside { N = n, Us = us, Usp = usp, V = v, Wm = wm, Wm1 = wm1, Wl = wl };
    }

    // Energy of the three individual stacks (i,j)→(i+1,j−1)→(i+2,j−2)→(i+3,j−3) that the GGUC/CUGG special replaces.
    private static double GgucSplitEnergy(TurnerModel m, byte[] s, int i, int j)
        => m.InteriorTerm(s, i, j, i + 1, j - 1)
         + m.InteriorTerm(s, i + 1, j - 1, i + 2, j - 2)
         + m.InteriorTerm(s, i + 2, j - 2, i + 3, j - 3);

    /// <summary>
    /// McCaskill's outside recursion on the grammar of <see cref="FillPartitionInside"/>: returns P(i,j) at index i*n+j
    /// for every admissible pair (Σ of the Boltzmann weights of all structures containing the pair, divided by Z).
    /// </summary>
    /// <remarks>
    /// Ōv(i,j) is the weight of everything outside a pair (i,j) given that it is formed: the exterior loop
    /// (Ŵ(i−1)·e^{−Ext}·Ŵr(j+1)), every enclosing interior loop (k,l), the GGUC correction, and (i,j) as a multiloop branch via
    /// ŴM1/ŴM (outside values Ōwm1, Ōwm; Ōg(i,j) = Σ_{jj≥j} Ōwm1(i,jj)·us^(jj−j)). Cells are processed by decreasing span, so
    /// every outside value is known before it is used; P(i,j) = V̂(i,j)·Ōv(i,j)/Ẑ.
    /// </remarks>
    private static double[] ComputePairProbabilities(TurnerModel m, byte[] s, int n, int minLoopSize, PfInside f)
    {
        const int MAXLOOP = Dp_MAXLOOP;
        double rt = m.RT;
        double us = f.Us;
        double[] usp = f.Usp, v = f.V, wm = f.Wm, wm1 = f.Wm1, wl = f.Wl;
        double zs = f.ZScaled;
        double ggucWeight = Math.Exp(-SpecialGGUC_CUGG_3Stack / rt);

        // Ŵr(i): scaled partition function of the suffix starting at i (positions i..n−1).
        var wr = new double[n + 2];
        wr[n] = 1.0;
        for (int i = n - 1; i >= 0; i--)
        {
            double acc = wr[i + 1] * us;
            for (int j = i + minLoopSize + 1; j < n; j++)
            {
                double vij = v[i * n + j];
                if (vij <= 0) continue;
                acc += vij * Math.Exp(-m.ExteriorStem(s, i, j) / rt) * wr[j + 1];
            }
            wr[i] = acc;
        }

        var vo = new double[n * n];     // Ōv
        var wmo = new double[n * n];    // Ōwm
        var wm1o = new double[n * n];   // Ōwm1
        var go = new double[n * n];     // Ōg
        var mlOut = new double[n * n];  // Ōv·e^{−Mclose}·us²: the pair (i,j) as the closing pair of a multiloop
        var prob = new double[n * n];

        for (int span = n; span >= minLoopSize + 2; span--)
        {
            for (int i = 0; i <= n - span; i++)
            {
                int j = i + span - 1;
                int ij = i * n + j;
                bool inMl = i > 0 && j < n - 1;

                if (inMl)
                {
                    // Ōwm(i,j): ŴM(i,j) as the left part of a longer ŴM(i,j2), or of a multiloop closed by (i−1,j0).
                    double t = 0.0;
                    for (int j2 = j + minLoopSize + 2; j2 <= n - 2; j2++)
                        t += wmo[i * n + j2] * wm1[(j + 1) * n + j2];
                    for (int j0 = j + minLoopSize + 3; j0 <= n - 1; j0++)
                        t += mlOut[(i - 1) * n + j0] * wm1[(j + 1) * n + j0 - 1];
                    wmo[ij] = t;

                    // Ōwm1(i,j): ŴM1(i,j) as the branch of ŴM(i2,j) (i2 ≤ i), or of a multiloop closed by (i0,j+1).
                    t = 0.0;
                    for (int i2 = 1; i2 <= i; i2++)
                    {
                        double wo = wmo[i2 * n + j];
                        if (wo == 0.0) continue;
                        t += wo * ((i > i2 ? wm[i2 * n + i - 1] : 0.0) + usp[i - i2]);
                    }
                    for (int i0 = 0; i0 <= i - minLoopSize - 3; i0++)
                        t += mlOut[i0 * n + j + 1] * wm[(i0 + 1) * n + i - 1];
                    wm1o[ij] = t;
                    go[ij] = t + (j + 1 <= n - 2 ? us * go[ij + 1] : 0.0);
                }

                if (v[ij] > 0)
                {
                    // (a) (i,j) faces the exterior loop.
                    double o = wl[i] * Math.Exp(-m.ExteriorStem(s, i, j) / rt) * wr[j + 1];

                    // (b) (i,j) is the inner pair of a stack / bulge / interior loop closed by (k,l).
                    for (int k = i - 1; k >= 0 && i - k - 1 <= MAXLOOP; k--)
                    {
                        int lMax = Math.Min(n - 1, j + 1 + (MAXLOOP - (i - k - 1)));
                        for (int l = j + 1; l <= lMax; l++)
                        {
                            double outer = vo[k * n + l];
                            if (outer <= 0) continue;
                            double e = m.InteriorTerm(s, k, l, i, j);
                            if (e >= Dp_INF) continue;
                            o += outer * Math.Exp(-e / rt) * usp[(i - k) + (l - j)];
                        }
                    }

                    // (c) fourth pair of a special GGUC/CUGG three-stack closed by (i−3,j+3).
                    if (i >= 3 && j + 3 <= n - 1 && vo[(i - 3) * n + j + 3] > 0
                        && TurnerModel.IsSpecialGgucCuggStack(s, i - 3, j + 3)
                        && v[(i - 2) * n + j + 2] > 0 && v[(i - 1) * n + j + 1] > 0)
                        o += vo[(i - 3) * n + j + 3] * (ggucWeight - Math.Exp(-GgucSplitEnergy(m, s, i - 3, j + 3) / rt)) * usp[6];

                    // (d) (i,j) is a branch of a multiloop.
                    if (inMl)
                        o += go[ij] * Math.Exp(-m.MultiloopStem(s, i, j) / rt);

                    vo[ij] = o;
                    mlOut[ij] = o * Math.Exp(-m.MultiloopClosing(s, i, j) / rt) * usp[2];
                    prob[ij] += v[ij] * o / zs;

                    // GGUC correction also belongs to the two middle pairs of the three-stack closed by (i,j).
                    if (TurnerModel.IsSpecialGgucCuggStack(s, i, j) && v[(i + 1) * n + j - 1] > 0 && v[(i + 2) * n + j - 2] > 0)
                    {
                        double inner3 = v[(i + 3) * n + j - 3];
                        if (inner3 > 0)
                        {
                            double corr = o * (ggucWeight - Math.Exp(-GgucSplitEnergy(m, s, i, j) / rt)) * inner3 * usp[6] / zs;
                            prob[(i + 1) * n + j - 1] += corr;
                            prob[(i + 2) * n + j - 2] += corr;
                        }
                    }
                }
            }
        }

        return prob;
    }

    #endregion
}
