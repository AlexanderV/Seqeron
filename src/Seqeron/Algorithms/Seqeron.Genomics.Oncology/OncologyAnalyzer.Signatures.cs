namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region SBS-96 Trinucleotide Context Catalog

    /// <summary>
    /// Number of single-base-substitution channels in the SBS-96 classification: six pyrimidine substitution
    /// subtypes × four 5' bases × four 3' bases. Source: Alexandrov et al. (2013), Nature 500:415-421 — "96
    /// possible mutation types (6 types of substitution * 4 types of 5' base * 4 types of 3' base)"; COSMIC
    /// SBS96 (https://cancer.sanger.ac.uk/signatures/sbs/sbs96/). Value = 96.
    /// </summary>
    public const int Sbs96ChannelCount = 96;

    /// <summary>
    /// The four 5'/3' flanking bases of an SBS-96 trinucleotide context, in canonical alphabetical order
    /// (A, C, G, T). Source: COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019, BMC Genomics
    /// 20:685) — each substitution has "sixteen possible trinucleotide (4 types of 5' base * 4 types of 3' base)".
    /// </summary>
    private static readonly char[] ContextBases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// The six single-base substitutions of the SBS-6 / SBS-96 classification, each referred to by the
    /// pyrimidine (C or T) of the mutated Watson-Crick base pair, in canonical order. Source: COSMIC SBS96 —
    /// "C>A, C>G, C>T, T>A, T>C, and T>G"; Alexandrov et al. (2013); Bergstrom et al. (2019).
    /// </summary>
    private static readonly (char Ref, char Alt)[] PyrimidineSubstitutions =
    {
        ('C', 'A'), ('C', 'G'), ('C', 'T'), ('T', 'A'), ('T', 'C'), ('T', 'G')
    };

    /// <summary>
    /// Classifies a single-base substitution into its SBS-96 trinucleotide-context channel, folded onto the
    /// pyrimidine strand. The mutated base sits in the centre of the trinucleotide; the channel is rendered as
    /// <c>5'[REF&gt;ALT]3'</c> (e.g. <c>A[C&gt;A]A</c>). Each substitution is referred to by the pyrimidine of
    /// the mutated Watson-Crick base pair: when the reference (mutated) base is a purine (A or G), the
    /// trinucleotide context and the substitution are reverse-complemented onto the pyrimidine strand before
    /// counting (e.g. a G&gt;T at 5'-T G A-3' folds to <c>T[C&gt;A]A</c>). Source: Alexandrov et al. (2013),
    /// Nature 500:415-421; COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019) — "using the purine
    /// base of the Watson-Crick base-pair for classifying mutation types will require taking the reverse
    /// complement sequence". Complement map A↔T, C↔G (Watson-Crick base pairing).
    /// </summary>
    /// <param name="fivePrime">Reference base immediately 5' of the mutated base (A/C/G/T, case-insensitive).</param>
    /// <param name="referenceBase">The mutated (reference) base (A/C/G/T, case-insensitive).</param>
    /// <param name="alternateBase">The substituted (alternate) base (A/C/G/T, case-insensitive, ≠ reference).</param>
    /// <param name="threePrime">Reference base immediately 3' of the mutated base (A/C/G/T, case-insensitive).</param>
    /// <returns>The SBS-96 channel label in the form <c>5'[REF&gt;ALT]3'</c> with a pyrimidine reference base.</returns>
    /// <exception cref="ArgumentException">A base is not A/C/G/T, or the reference and alternate bases are equal.</exception>
    public static string ClassifySbsContext(char fivePrime, char referenceBase, char alternateBase, char threePrime)
    {
        char five = NormalizeBase(fivePrime, nameof(fivePrime));
        char reference = NormalizeBase(referenceBase, nameof(referenceBase));
        char alternate = NormalizeBase(alternateBase, nameof(alternateBase));
        char three = NormalizeBase(threePrime, nameof(threePrime));

        if (reference == alternate)
        {
            throw new ArgumentException(
                $"A substitution requires reference ≠ alternate (got '{reference}' = '{alternate}').",
                nameof(alternateBase));
        }

        // Fold purine-reference substitutions onto the pyrimidine strand by reverse-complementing the
        // trinucleotide context AND the substitution (SigProfiler / COSMIC). For a pyrimidine reference
        // (C or T) the mutation is already on the pyrimidine strand and is kept as-is.
        if (reference is 'A' or 'G')
        {
            char foldedFive = Complement(three);   // 3' neighbour becomes the 5' neighbour after reversal
            char foldedThree = Complement(five);   // 5' neighbour becomes the 3' neighbour after reversal
            reference = Complement(reference);
            alternate = Complement(alternate);
            five = foldedFive;
            three = foldedThree;
        }

        return $"{five}[{reference}>{alternate}]{three}";
    }

    /// <summary>
    /// Enumerates all 96 canonical SBS-96 channel labels (<c>5'[REF&gt;ALT]3'</c>), in deterministic
    /// substitution-major order: the six pyrimidine substitutions (C&gt;A, C&gt;G, C&gt;T, T&gt;A, T&gt;C,
    /// T&gt;G), then 5' base (A,C,G,T), then 3' base (A,C,G,T). Source: COSMIC SBS96; Alexandrov et al. (2013)
    /// — 6 × 4 × 4 = 96. The ordering is a presentation convention and does not affect per-variant classification.
    /// </summary>
    /// <returns>The 96 distinct channel labels.</returns>
    public static IReadOnlyList<string> EnumerateSbs96Channels()
    {
        var channels = new List<string>(Sbs96ChannelCount);
        foreach (var (reference, alternate) in PyrimidineSubstitutions)
        {
            foreach (char five in ContextBases)
            {
                foreach (char three in ContextBases)
                {
                    channels.Add($"{five}[{reference}>{alternate}]{three}");
                }
            }
        }

        return channels;
    }

    /// <summary>
    /// Builds the SBS-96 mutational catalog (the 96-channel spectrum) from a collection of single-base
    /// substitutions, each given as its 5' base, reference (mutated) base, alternate base, and 3' base.
    /// Every variant is classified via <see cref="ClassifySbsContext"/> (with pyrimidine-strand folding) and
    /// tallied into its channel. All 96 channels are present in the result, including those with a zero count,
    /// so the spectrum has a fixed shape. The sum of the counts equals the number of input variants (each
    /// classifiable variant contributes exactly one count — the catalog is a partition). Source: Alexandrov
    /// et al. (2013); COSMIC SBS96; SigProfilerMatrixGenerator (Bergstrom et al. 2019).
    /// </summary>
    /// <param name="variants">
    /// SBS variants as (FivePrime, Reference, Alternate, ThreePrime) tuples; each must be a valid single-base
    /// substitution (A/C/G/T bases, reference ≠ alternate).
    /// </param>
    /// <returns>
    /// A dictionary keyed by all 96 channel labels mapping to the count of variants in each channel.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException">Any variant has an invalid base or reference == alternate.</exception>
    public static IReadOnlyDictionary<string, int> Build96ContextCatalog(
        IEnumerable<(char FivePrime, char Reference, char Alternate, char ThreePrime)> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        // Initialise all 96 channels to zero so the spectrum always has the full, fixed shape.
        var catalog = new Dictionary<string, int>(Sbs96ChannelCount, StringComparer.Ordinal);
        foreach (string channel in EnumerateSbs96Channels())
        {
            catalog[channel] = 0;
        }

        foreach (var (fivePrime, reference, alternate, threePrime) in variants)
        {
            string channel = ClassifySbsContext(fivePrime, reference, alternate, threePrime);
            catalog[channel]++;
        }

        return catalog;
    }

    /// <summary>
    /// Returns the Watson-Crick complement of a DNA base (A↔T, C↔G). Source: complementary base pairing,
    /// adenine pairs with thymine and cytosine pairs with guanine.
    /// </summary>
    private static char Complement(char baseChar) => baseChar switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"'{baseChar}' is not a DNA base (A/C/G/T).", nameof(baseChar))
    };

    /// <summary>
    /// Validates and upper-cases a single DNA base, rejecting anything that is not A/C/G/T.
    /// </summary>
    private static char NormalizeBase(char baseChar, string paramName)
    {
        char upper = char.ToUpperInvariant(baseChar);
        if (upper is not ('A' or 'C' or 'G' or 'T'))
        {
            throw new ArgumentException($"'{baseChar}' is not a DNA base (A/C/G/T).", paramName);
        }

        return upper;
    }

    #endregion


    #region Signature Fitting / Refitting (ONCO-SIG-002)

    /// <summary>
    /// Convergence tolerance ε for the Lawson-Hanson active-set NNLS main loop: the iteration stops when the
    /// largest gradient component over the inactive set R is ≤ ε. Source: Lawson C.L. &amp; Hanson R.J. (1974),
    /// <i>Solving Least Squares Problems</i>, Ch. 23 — the active-set algorithm terminates when
    /// max(w_R) ≤ ε (https://en.wikipedia.org/wiki/Non-negative_least_squares). A small positive value
    /// (1e-12) makes the stop effectively exact for the well-conditioned signature matrices used here.
    /// </summary>
    private const double NnlsTolerance = 1e-12;

    /// <summary>
    /// Maximum number of outer (active-set growth) iterations for the NNLS solver. The Lawson-Hanson method
    /// adds at most one index per outer iteration and is guaranteed to terminate in a finite number of steps;
    /// this cap (a small multiple of the number of signatures) is a safety bound against floating-point
    /// non-termination. Source: Lawson &amp; Hanson (1974), Ch. 23 (finite termination of the active-set method).
    /// </summary>
    private const int NnlsMaxOuterIterationsPerSignature = 30;

    /// <summary>
    /// The result of fitting (refitting) an observed mutational catalog to a set of reference signatures.
    /// </summary>
    /// <param name="Exposures">
    /// The fitted non-negative contribution (weight) of each reference signature, in signature order — the
    /// NNLS solution x of min‖S·x − d‖² with x ≥ 0 (Blokzijl et al. 2018; Lawson &amp; Hanson 1974).
    /// </param>
    /// <param name="NormalizedExposures">
    /// <paramref name="Exposures"/> divided by their sum, so they sum to 1 when the total is positive (all
    /// zero otherwise) — the proportion of mutations attributed to each signature (Rosenthal et al. 2016).
    /// </param>
    /// <param name="Reconstruction">The reconstructed catalog S·x (Rosenthal et al. 2016, R = S·W).</param>
    /// <param name="ReconstructionCosineSimilarity">
    /// Cosine similarity between the observed catalog d and the reconstruction S·x — the reconstruction
    /// quality measure (Blokzijl et al. 2018); 1.0 means the catalog is exactly representable.
    /// </param>
    public readonly record struct SignatureFitResult(
        IReadOnlyList<double> Exposures,
        IReadOnlyList<double> NormalizedExposures,
        IReadOnlyList<double> Reconstruction,
        double ReconstructionCosineSimilarity);

    /// <summary>
    /// Computes the cosine similarity between two equal-length non-negative vectors:
    /// <code>sim(A,B) = Σᵢ AᵢBᵢ / ( √(Σᵢ Aᵢ²) · √(Σᵢ Bᵢ²) )</code>
    /// i.e. the dot product divided by the product of the Euclidean norms — the cosine of the angle between
    /// the two vectors. The value lies in [0, 1] for non-negative inputs (0 = independent / orthogonal,
    /// 1 = identical direction) and is invariant to positive scaling of either vector. Source: Blokzijl et al.
    /// (2018), <i>Genome Medicine</i> 10:33 (§ "Mutational profile similarity"); Pan &amp; Wang (2020), iMutSig.
    /// When either vector has zero Euclidean norm the cosine is undefined (division by zero); this method
    /// returns 0.0 for that degenerate case (no shared direction).
    /// </summary>
    /// <param name="a">First vector (e.g. a mutational profile / 96-channel catalog).</param>
    /// <param name="b">Second vector of the same length as <paramref name="a"/>.</param>
    /// <returns>Cosine similarity in [0, 1]; 0.0 when either vector is all-zero.</returns>
    /// <exception cref="ArgumentNullException">Either argument is null.</exception>
    /// <exception cref="ArgumentException">The two vectors have different lengths, or are empty.</exception>
    public static double CosineSimilarity(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        if (a.Count == 0 || b.Count == 0)
        {
            throw new ArgumentException("Cosine similarity is undefined for empty vectors.", nameof(a));
        }

        if (a.Count != b.Count)
        {
            throw new ArgumentException(
                $"Vectors must have the same length (got {a.Count} and {b.Count}).", nameof(b));
        }

        double dot = 0.0;
        double normASquared = 0.0;
        double normBSquared = 0.0;
        for (int i = 0; i < a.Count; i++)
        {
            dot += a[i] * b[i];
            normASquared += a[i] * a[i];
            normBSquared += b[i] * b[i];
        }

        if (normASquared == 0.0 || normBSquared == 0.0)
        {
            // Undefined (zero-norm vector has no direction); treated as no shared direction.
            return 0.0;
        }

        return dot / (Math.Sqrt(normASquared) * Math.Sqrt(normBSquared));
    }

    /// <summary>
    /// Reconstructs a mutational catalog from reference signatures and their exposures: the matrix-vector
    /// product S·x, where each entry is Σⱼ signatures[j][channel] · exposures[j]. Source: Rosenthal et al.
    /// (2016), <i>Genome Biology</i> 17:31 — the reconstructed profile is S·W.
    /// </summary>
    /// <param name="signatures">
    /// Reference signatures as a list of equal-length channel vectors (one vector per signature; element
    /// <c>signatures[j][k]</c> is the weight of signature j in channel k).
    /// </param>
    /// <param name="exposures">The per-signature exposure (weight) of each signature, same count as signatures.</param>
    /// <returns>The reconstructed catalog vector of length equal to the signature channel count.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">No signatures, ragged signatures, or count mismatch with exposures.</exception>
    public static IReadOnlyList<double> ReconstructCatalog(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> exposures)
    {
        int channelCount = ValidateSignatures(signatures);
        ArgumentNullException.ThrowIfNull(exposures);

        if (exposures.Count != signatures.Count)
        {
            throw new ArgumentException(
                $"Exposure count ({exposures.Count}) must equal signature count ({signatures.Count}).",
                nameof(exposures));
        }

        var reconstruction = new double[channelCount];
        for (int j = 0; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j];
            double weight = exposures[j];
            for (int k = 0; k < channelCount; k++)
            {
                reconstruction[k] += signature[k] * weight;
            }
        }

        return reconstruction;
    }

    /// <summary>
    /// Fits (refits) an observed mutational catalog to a set of caller-supplied reference signatures by solving
    /// the non-negative least-squares problem
    /// <code>minₓ ‖ S·x − d ‖₂²,  subject to x ≥ 0</code>
    /// where the columns of S are the reference signatures, d is the observed catalog, and x is the fitted
    /// per-signature exposure (contribution) vector. The decomposition is the standard signature-refitting
    /// model: the catalog is projected onto the non-negative cone spanned by the reference signatures
    /// (Blokzijl et al. 2018). Reference signature profiles are <b>not</b> hardcoded — they are supplied by the
    /// caller (e.g. COSMIC SBS profiles); this method only performs the fit. The NNLS problem is solved with
    /// the Lawson-Hanson active-set algorithm (Lawson &amp; Hanson 1974). The result also exposes the
    /// proportion-normalised exposures (Rosenthal et al. 2016), the reconstruction S·x, and the cosine
    /// similarity between d and S·x as a reconstruction-quality measure (Blokzijl et al. 2018).
    /// </summary>
    /// <param name="catalog">The observed mutational catalog d (e.g. a 96-channel SBS spectrum), non-negative.</param>
    /// <param name="signatures">Reference signatures as equal-length channel vectors (one per signature).</param>
    /// <returns>The fit result: exposures, normalised exposures, reconstruction, and reconstruction cosine.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">
    /// No signatures, ragged signatures, or the catalog length differs from the signature channel count.
    /// </exception>
    public static SignatureFitResult FitSignatures(
        IReadOnlyList<double> catalog,
        IReadOnlyList<IReadOnlyList<double>> signatures)
    {
        ArgumentNullException.ThrowIfNull(catalog);
        int channelCount = ValidateSignatures(signatures);

        if (catalog.Count != channelCount)
        {
            throw new ArgumentException(
                $"Catalog length ({catalog.Count}) must equal the signature channel count ({channelCount}).",
                nameof(catalog));
        }

        double[] exposures = SolveNonNegativeLeastSquares(signatures, catalog, channelCount);

        IReadOnlyList<double> reconstruction = ReconstructCatalog(signatures, exposures);
        double[] normalized = NormalizeExposures(exposures);
        double reconstructionCosine = CosineSimilarity(catalog, reconstruction);

        return new SignatureFitResult(exposures, normalized, reconstruction, reconstructionCosine);
    }

    /// <summary>
    /// Normalises exposures into proportions that sum to 1 (each divided by the total), or all zeros when the
    /// total is zero. Source: Rosenthal et al. (2016) — "the weights W are normalized between 0 and 1".
    /// </summary>
    private static double[] NormalizeExposures(double[] exposures)
    {
        double sum = 0.0;
        foreach (double e in exposures)
        {
            sum += e;
        }

        var normalized = new double[exposures.Length];
        if (sum > 0.0)
        {
            for (int i = 0; i < exposures.Length; i++)
            {
                normalized[i] = exposures[i] / sum;
            }
        }

        return normalized;
    }

    /// <summary>
    /// Solves minₓ ‖ S·x − d ‖₂² subject to x ≥ 0 with the Lawson-Hanson active-set algorithm.
    /// Source: Lawson C.L. &amp; Hanson R.J. (1974), <i>Solving Least Squares Problems</i>, Ch. 23
    /// (https://en.wikipedia.org/wiki/Non-negative_least_squares). Index set P holds the passive (free,
    /// possibly non-zero) variables; R holds the active (clamped-to-zero) variables; the gradient
    /// w = Sᵀ(d − S·x) selects the next variable to free.
    /// </summary>
    /// <param name="signatures">Signature matrix S (column j = signatures[j]).</param>
    /// <param name="catalog">Observed vector d.</param>
    /// <param name="channelCount">Number of channels (rows of S, length of d).</param>
    /// <returns>The NNLS solution x (length = number of signatures).</returns>
    private static double[] SolveNonNegativeLeastSquares(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> catalog,
        int channelCount)
    {
        int n = signatures.Count;
        var x = new double[n];
        bool[] passive = new bool[n]; // true => index in P, false => in R

        int maxOuter = n * NnlsMaxOuterIterationsPerSignature;
        int outer = 0;

        while (outer++ < maxOuter)
        {
            // w = Sᵀ(d − S·x); only inactive (R) components matter for selection.
            double[] residual = ComputeResidual(signatures, x, catalog, channelCount);
            double[] gradient = ComputeGradient(signatures, residual);

            int j = -1;
            double maxGradient = NnlsTolerance;
            for (int i = 0; i < n; i++)
            {
                if (!passive[i] && gradient[i] > maxGradient)
                {
                    maxGradient = gradient[i];
                    j = i;
                }
            }

            if (j < 0)
            {
                // R empty or max(w_R) ≤ ε — KKT conditions satisfied.
                break;
            }

            passive[j] = true;

            // Inner loop: solve the unconstrained LS on P; if any becomes ≤ 0, take the bounded step.
            int innerGuard = 0;
            while (innerGuard++ <= n)
            {
                double[] s = SolveLeastSquaresOnPassiveSet(signatures, catalog, passive, channelCount, n);

                double minPassive = double.PositiveInfinity;
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && s[i] < minPassive)
                    {
                        minPassive = s[i];
                    }
                }

                if (minPassive > 0.0)
                {
                    x = s;
                    break;
                }

                // α = min over i in P with s_i ≤ 0 of x_i / (x_i − s_i).
                double alpha = double.PositiveInfinity;
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && s[i] <= 0.0)
                    {
                        double denom = x[i] - s[i];
                        if (denom != 0.0)
                        {
                            double candidate = x[i] / denom;
                            if (candidate < alpha)
                            {
                                alpha = candidate;
                            }
                        }
                    }
                }

                if (double.IsPositiveInfinity(alpha))
                {
                    // Numerical safeguard: no feasible step (should not occur for valid inputs).
                    x = s;
                    break;
                }

                for (int i = 0; i < n; i++)
                {
                    x[i] += alpha * (s[i] - x[i]);
                }

                // Move indices with x_i ≤ 0 from P back to R.
                for (int i = 0; i < n; i++)
                {
                    if (passive[i] && x[i] <= 0.0)
                    {
                        x[i] = 0.0;
                        passive[i] = false;
                    }
                }
            }
        }

        return x;
    }

    /// <summary>Computes the residual d − S·x.</summary>
    private static double[] ComputeResidual(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        double[] x,
        IReadOnlyList<double> catalog,
        int channelCount)
    {
        var residual = new double[channelCount];
        for (int k = 0; k < channelCount; k++)
        {
            residual[k] = catalog[k];
        }

        for (int j = 0; j < signatures.Count; j++)
        {
            double weight = x[j];
            if (weight == 0.0)
            {
                continue;
            }

            IReadOnlyList<double> signature = signatures[j];
            for (int k = 0; k < channelCount; k++)
            {
                residual[k] -= signature[k] * weight;
            }
        }

        return residual;
    }

    /// <summary>Computes the gradient w = Sᵀ·residual.</summary>
    private static double[] ComputeGradient(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        double[] residual)
    {
        var gradient = new double[signatures.Count];
        for (int j = 0; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j];
            double sum = 0.0;
            for (int k = 0; k < residual.Length; k++)
            {
                sum += signature[k] * residual[k];
            }

            gradient[j] = sum;
        }

        return gradient;
    }

    /// <summary>
    /// Solves the unconstrained least-squares problem restricted to the passive set P:
    /// s_P = ((S_P)ᵀ S_P)⁻¹ (S_P)ᵀ d, with s_R = 0, via the normal equations solved by Gaussian elimination.
    /// Source: Lawson &amp; Hanson (1974), Ch. 23.
    /// </summary>
    private static double[] SolveLeastSquaresOnPassiveSet(
        IReadOnlyList<IReadOnlyList<double>> signatures,
        IReadOnlyList<double> catalog,
        bool[] passive,
        int channelCount,
        int n)
    {
        var passiveIndices = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (passive[i])
            {
                passiveIndices.Add(i);
            }
        }

        int p = passiveIndices.Count;
        var s = new double[n];
        if (p == 0)
        {
            return s;
        }

        // Normal equations: (S_Pᵀ S_P) z = S_Pᵀ d.
        var ata = new double[p, p];
        var atb = new double[p];
        for (int a = 0; a < p; a++)
        {
            IReadOnlyList<double> sigA = signatures[passiveIndices[a]];
            double rhs = 0.0;
            for (int k = 0; k < channelCount; k++)
            {
                rhs += sigA[k] * catalog[k];
            }

            atb[a] = rhs;

            for (int b = 0; b < p; b++)
            {
                IReadOnlyList<double> sigB = signatures[passiveIndices[b]];
                double dot = 0.0;
                for (int k = 0; k < channelCount; k++)
                {
                    dot += sigA[k] * sigB[k];
                }

                ata[a, b] = dot;
            }
        }

        double[] z = SolveLinearSystem(ata, atb, p);
        for (int a = 0; a < p; a++)
        {
            s[passiveIndices[a]] = z[a];
        }

        return s;
    }

    /// <summary>
    /// Solves the dense linear system M·z = rhs (M is p×p, symmetric positive semi-definite here) by Gaussian
    /// elimination with partial pivoting. Standard direct method (CLRS, §28; Numerical Recipes §2.1).
    /// </summary>
    private static double[] SolveLinearSystem(double[,] matrix, double[] rhs, int p)
    {
        // Work on copies so the inputs are not mutated.
        var m = new double[p, p];
        var b = new double[p];
        for (int i = 0; i < p; i++)
        {
            b[i] = rhs[i];
            for (int k = 0; k < p; k++)
            {
                m[i, k] = matrix[i, k];
            }
        }

        for (int col = 0; col < p; col++)
        {
            // Partial pivot: largest magnitude in this column at or below the diagonal.
            int pivot = col;
            double best = Math.Abs(m[col, col]);
            for (int row = col + 1; row < p; row++)
            {
                double magnitude = Math.Abs(m[row, col]);
                if (magnitude > best)
                {
                    best = magnitude;
                    pivot = row;
                }
            }

            if (pivot != col)
            {
                for (int k = 0; k < p; k++)
                {
                    (m[col, k], m[pivot, k]) = (m[pivot, k], m[col, k]);
                }

                (b[col], b[pivot]) = (b[pivot], b[col]);
            }

            double diagonal = m[col, col];
            if (diagonal == 0.0)
            {
                // Singular column (collinear signatures); leave this component at 0.
                continue;
            }

            for (int row = col + 1; row < p; row++)
            {
                double factor = m[row, col] / diagonal;
                if (factor == 0.0)
                {
                    continue;
                }

                for (int k = col; k < p; k++)
                {
                    m[row, k] -= factor * m[col, k];
                }

                b[row] -= factor * b[col];
            }
        }

        // Back-substitution.
        var z = new double[p];
        for (int row = p - 1; row >= 0; row--)
        {
            double sum = b[row];
            for (int k = row + 1; k < p; k++)
            {
                sum -= m[row, k] * z[k];
            }

            double diagonal = m[row, row];
            z[row] = diagonal == 0.0 ? 0.0 : sum / diagonal;
        }

        return z;
    }

    /// <summary>
    /// Validates a signature matrix (non-null, at least one signature, each signature non-null and equal-length,
    /// non-empty) and returns the common channel count.
    /// </summary>
    private static int ValidateSignatures(IReadOnlyList<IReadOnlyList<double>> signatures)
    {
        ArgumentNullException.ThrowIfNull(signatures);

        if (signatures.Count == 0)
        {
            throw new ArgumentException("At least one reference signature is required.", nameof(signatures));
        }

        IReadOnlyList<double> first = signatures[0]
            ?? throw new ArgumentException("Signature vectors cannot be null.", nameof(signatures));
        int channelCount = first.Count;
        if (channelCount == 0)
        {
            throw new ArgumentException("Signature vectors cannot be empty.", nameof(signatures));
        }

        for (int j = 1; j < signatures.Count; j++)
        {
            IReadOnlyList<double> signature = signatures[j]
                ?? throw new ArgumentException("Signature vectors cannot be null.", nameof(signatures));
            if (signature.Count != channelCount)
            {
                throw new ArgumentException(
                    $"All signatures must have the same length (signature 0 has {channelCount}, " +
                    $"signature {j} has {signature.Count}).",
                    nameof(signatures));
            }
        }

        return channelCount;
    }

    #endregion


    #region De-novo Signature Extraction via NMF (ONCO-SIG-002)

    /// <summary>
    /// Maximum number of multiplicative-update iterations for the de-novo NMF signature extraction. Source:
    /// Lee &amp; Seung (2001), <i>Algorithms for Non-negative Matrix Factorization</i>, NIPS 13 — the
    /// multiplicative updates are "applied iteratively until W and H converge"
    /// (https://en.wikipedia.org/wiki/Non-negative_matrix_factorization). NMF is non-convex, so a finite cap is
    /// a safety bound; this default mirrors the iteration budgets used by SigProfiler-style extractors.
    /// </summary>
    public const int DefaultNmfMaxIterations = 10_000;

    /// <summary>
    /// Default relative-improvement convergence tolerance for the NMF objective: iteration stops when the
    /// per-iteration decrease of the Frobenius residual ‖V − WH‖²_F, relative to the previous value, drops below
    /// this threshold. Source: Lee &amp; Seung (2001), Theorem 1 — the objective is monotonically non-increasing,
    /// so a small relative-change stop is a valid convergence test. Value = 1e-10.
    /// </summary>
    public const double DefaultNmfTolerance = 1e-10;

    /// <summary>
    /// Default RNG seed for the non-negative random initialisation of the NMF factors. Fixed so that, for a
    /// given count matrix V, rank k, iteration budget and tolerance, the extracted signatures and exposures are
    /// reproducible (NMF is non-convex and initialisation-dependent). Mirrors the fixed-seed convention used by
    /// <see cref="DefaultBootstrapSeed"/>. Value = 42.
    /// </summary>
    public const int DefaultNmfSeed = 42;

    /// <summary>
    /// A small additive floor on the multiplicative-update denominators (and on the random initialisation) to
    /// avoid 0/0 when a row or column of a factor collapses to zero. Source: standard regularisation of the
    /// Lee &amp; Seung multiplicative updates whose denominators (WᵀWH, WHHᵀ) can vanish
    /// (https://en.wikipedia.org/wiki/Non-negative_matrix_factorization). Value = 1e-12, far below any
    /// meaningful mutation count.
    /// </summary>
    private const double NmfEpsilon = 1e-12;

    /// <summary>
    /// The result of de-novo NMF signature extraction from a mutation-count matrix V ≈ W·H at a caller-specified
    /// rank k (Lee &amp; Seung 2001; Alexandrov et al. 2013).
    /// </summary>
    /// <param name="Signatures">
    /// The extracted signatures W as a list of k channel vectors (one vector of length <c>ChannelCount</c> per
    /// signature). Each signature is L1-normalised so its channel weights sum to 1 — a probability distribution
    /// across the mutation channels, per the COSMIC / SigProfiler convention (Alexandrov et al. 2020).
    /// </param>
    /// <param name="Exposures">
    /// The exposures H as a k × samples matrix (<c>Exposures[j][s]</c> = activity of signature j in sample s).
    /// Non-negative; absorbs the per-signature scale removed by L1-normalising <paramref name="Signatures"/>.
    /// </param>
    /// <param name="FinalResidual">
    /// The final Frobenius reconstruction residual ‖V − W·H‖²_F (squared) after the last iteration.
    /// </param>
    /// <param name="Iterations">The number of multiplicative-update iterations actually performed.</param>
    /// <param name="ObjectiveHistory">
    /// The Frobenius residual ‖V − W·H‖²_F recorded after each iteration (length = <see cref="Iterations"/>),
    /// which is monotonically non-increasing (Lee &amp; Seung 2001, Theorem 1).
    /// </param>
    public readonly record struct SignatureExtractionResult(
        IReadOnlyList<IReadOnlyList<double>> Signatures,
        IReadOnlyList<IReadOnlyList<double>> Exposures,
        double FinalResidual,
        int Iterations,
        IReadOnlyList<double> ObjectiveHistory);

    /// <summary>
    /// De-novo mutational-signature extraction by Non-negative Matrix Factorization (NMF). Given a non-negative
    /// mutation-count matrix V (channels × samples) and a caller-specified rank k, factorises
    /// <code>V ≈ W·H,  W ≥ 0  (channels × k),  H ≥ 0  (k × samples)</code>
    /// where the columns of W are the de-novo signatures and H holds their per-sample exposures. Unlike
    /// <see cref="FitSignatures"/> (which refits exposures against caller-supplied reference signatures), the
    /// signatures here are <b>discovered from the data</b> — no reference profiles are required. The factors are
    /// found with the Lee &amp; Seung (2001) multiplicative update rules for the squared Euclidean (Frobenius)
    /// objective ‖V − W·H‖²_F:
    /// <code>
    /// H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)
    /// W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)
    /// </code>
    /// iterated until the relative decrease of the objective falls below <paramref name="tolerance"/> or
    /// <paramref name="maxIterations"/> is reached. Each extracted signature column of W is then L1-normalised to
    /// sum to 1 (a probability distribution over the channels), with the removed scale absorbed into H, per the
    /// COSMIC / SigProfiler convention (Alexandrov et al. 2013, 2020). NMF is non-convex, so the factorisation is
    /// a local optimum dependent on the (seeded, deterministic) non-negative random initialisation.
    /// </summary>
    /// <param name="countMatrix">
    /// The mutation-count matrix V as a list of rows, one per channel (e.g. 96 SBS channels); each row is a
    /// vector over the samples (<c>countMatrix[channel][sample]</c>). All entries must be finite and ≥ 0. Every
    /// row must have the same length (the sample count), and there must be at least one sample.
    /// </param>
    /// <param name="rank">The number of signatures k to extract; 1 ≤ k ≤ channel count.</param>
    /// <param name="maxIterations">Maximum multiplicative-update iterations (&gt; 0).</param>
    /// <param name="tolerance">Relative-improvement convergence tolerance (≥ 0).</param>
    /// <param name="seed">RNG seed for the non-negative random initialisation (reproducibility).</param>
    /// <returns>The extracted signatures W (L1-normalised), exposures H, residual, iteration count and history.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="countMatrix"/> or any row is null.</exception>
    /// <exception cref="ArgumentException">
    /// Empty matrix, zero samples, ragged rows, a negative or non-finite entry, rank &lt; 1, rank &gt; channel
    /// count, maxIterations ≤ 0, or tolerance &lt; 0.
    /// </exception>
    public static SignatureExtractionResult ExtractSignatures(
        IReadOnlyList<IReadOnlyList<double>> countMatrix,
        int rank,
        int maxIterations = DefaultNmfMaxIterations,
        double tolerance = DefaultNmfTolerance,
        int seed = DefaultNmfSeed)
        => ExtractSignatures(countMatrix, rank, NmfObjective.Frobenius, maxIterations, tolerance, seed);

    /// <summary>
    /// Selects the NMF objective (and hence the multiplicative-update variant) used by
    /// <see cref="ExtractSignatures(IReadOnlyList{IReadOnlyList{double}}, int, NmfObjective, int, double, int)"/>.
    /// Both variants are from Lee &amp; Seung (2001) and both are monotonically non-increasing in their respective
    /// objective (Theorems 1 and 2).
    /// </summary>
    public enum NmfObjective
    {
        /// <summary>
        /// Squared-Euclidean (Frobenius) objective ‖V − WH‖²_F with the Theorem-1 multiplicative updates
        /// <c>H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)</c>, <c>W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)</c>. Lee &amp; Seung (2001), Theorem 1.
        /// </summary>
        Frobenius = 0,

        /// <summary>
        /// Generalized Kullback-Leibler (Poisson) divergence
        /// <c>D(V‖WH) = Σ ( V·log(V/WH) − V + WH )</c> with the Theorem-2 multiplicative updates
        /// <c>H_aμ ← H_aμ · (Σ_i W_ia V_iμ/(WH)_iμ) / (Σ_i W_ia)</c> and
        /// <c>W_ia ← W_ia · (Σ_μ H_aμ V_iμ/(WH)_iμ) / (Σ_μ H_aμ)</c>. This is the objective the SigProfiler
        /// mutational-signature extractor actually optimises (Alexandrov et al. 2013; Islam et al. 2022).
        /// Lee &amp; Seung (2001), Theorem 2.
        /// </summary>
        KullbackLeibler = 1,
    }

    /// <summary>
    /// De-novo mutational-signature extraction by NMF (V ≈ W·H) at a caller-specified rank k, allowing the caller
    /// to pick the objective: the squared-Euclidean (Frobenius) variant (Lee &amp; Seung 2001, Theorem 1) or the
    /// generalized Kullback-Leibler / Poisson variant (Theorem 2) — the latter being the objective SigProfiler
    /// uses for mutational signatures (Alexandrov et al. 2013; Islam et al. 2022). The chosen objective is
    /// monotonically non-increasing across iterations; iteration stops on relative-improvement
    /// <paramref name="tolerance"/> or <paramref name="maxIterations"/>. Signatures (columns of W) are then
    /// L1-normalised to sum to 1 with the scale absorbed into H (COSMIC convention).
    /// <see cref="SignatureExtractionResult.FinalResidual"/> and
    /// <see cref="SignatureExtractionResult.ObjectiveHistory"/> hold the value of the <em>selected</em>
    /// objective (Frobenius residual ‖V − WH‖²_F, or KL divergence D(V‖WH)).
    /// </summary>
    /// <param name="countMatrix">Mutation-count matrix V (channels × samples); rows are channels.</param>
    /// <param name="rank">Number of signatures k to extract; 1 ≤ k ≤ channel count.</param>
    /// <param name="objective">Which Lee &amp; Seung objective / update variant to optimise.</param>
    /// <param name="maxIterations">Maximum multiplicative-update iterations (&gt; 0).</param>
    /// <param name="tolerance">Relative-improvement convergence tolerance (≥ 0).</param>
    /// <param name="seed">RNG seed for the non-negative random initialisation (reproducibility).</param>
    /// <returns>The extracted signatures W (L1-normalised), exposures H, objective value, iteration count, history.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="countMatrix"/> or any row is null.</exception>
    /// <exception cref="ArgumentException">Invalid matrix, rank, maxIterations or tolerance (see other overload).</exception>
    public static SignatureExtractionResult ExtractSignatures(
        IReadOnlyList<IReadOnlyList<double>> countMatrix,
        int rank,
        NmfObjective objective,
        int maxIterations = DefaultNmfMaxIterations,
        double tolerance = DefaultNmfTolerance,
        int seed = DefaultNmfSeed)
    {
        double[][] v = ValidateCountMatrix(countMatrix, out int channelCount, out int sampleCount);
        ValidateRankAndStop(rank, channelCount, maxIterations, tolerance);

        (double[][] w, double[][] h, double finalObjective, int iterations, List<double> history) =
            RunNmf(v, rank, objective, maxIterations, tolerance, seed, channelCount, sampleCount);

        IReadOnlyList<IReadOnlyList<double>> signatures = TransposeColumnsToSignatures(w, channelCount, rank);
        IReadOnlyList<IReadOnlyList<double>> exposures = RowsToReadOnly(h);

        return new SignatureExtractionResult(signatures, exposures, finalObjective, iterations, history);
    }

    /// <summary>
    /// Validates rank and stop-criterion parameters shared by the extraction entry points.
    /// </summary>
    private static void ValidateRankAndStop(int rank, int channelCount, int maxIterations, double tolerance)
    {
        if (rank < 1)
        {
            throw new ArgumentException($"Rank k must be ≥ 1 (got {rank}).", nameof(rank));
        }

        if (rank > channelCount)
        {
            throw new ArgumentException(
                $"Rank k ({rank}) cannot exceed the channel count ({channelCount}).", nameof(rank));
        }

        if (maxIterations <= 0)
        {
            throw new ArgumentException($"maxIterations must be > 0 (got {maxIterations}).", nameof(maxIterations));
        }

        if (tolerance < 0)
        {
            throw new ArgumentException($"tolerance must be ≥ 0 (got {tolerance}).", nameof(tolerance));
        }
    }

    /// <summary>
    /// Core multiplicative-update NMF loop shared by both objectives. Returns the (L1-column-normalised) factors,
    /// the final value of the selected objective, the iteration count, and the per-iteration objective history.
    /// The objective is monotonically non-increasing (Lee &amp; Seung 2001, Theorems 1 &amp; 2).
    /// </summary>
    private static (double[][] W, double[][] H, double FinalObjective, int Iterations, List<double> History) RunNmf(
        double[][] v, int rank, NmfObjective objective, int maxIterations, double tolerance, int seed,
        int channelCount, int sampleCount)
    {
        // Non-negative random initialisation (Lee & Seung do not prescribe one; uniform (0,1] is standard).
        var rng = new Random(seed);
        double[][] w = InitializeNonNegativeFactor(rng, channelCount, rank);   // channels × k
        double[][] h = InitializeNonNegativeFactor(rng, rank, sampleCount);    // k × samples

        var history = new List<double>(maxIterations);
        double previousObjective = double.PositiveInfinity;
        int iterations = 0;

        for (int iter = 0; iter < maxIterations; iter++)
        {
            if (objective == NmfObjective.KullbackLeibler)
            {
                // H_aμ ← H_aμ · (Σ_i W_ia V_iμ/(WH)_iμ) / (Σ_i W_ia)   — Lee & Seung (2001), Theorem 2.
                UpdateHKl(w, h, v, channelCount, rank, sampleCount);
                // W_ia ← W_ia · (Σ_μ H_aμ V_iμ/(WH)_iμ) / (Σ_μ H_aμ)   — Lee & Seung (2001), Theorem 2.
                UpdateWKl(w, h, v, channelCount, rank, sampleCount);
            }
            else
            {
                // H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)   — Lee & Seung (2001), Theorem 1.
                UpdateH(w, h, v, channelCount, rank, sampleCount);
                // W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)   — Lee & Seung (2001), Theorem 1.
                UpdateW(w, h, v, channelCount, rank, sampleCount);
            }

            iterations = iter + 1;
            double currentObjective = objective == NmfObjective.KullbackLeibler
                ? KullbackLeiblerDivergence(w, h, v, channelCount, rank, sampleCount)
                : FrobeniusResidualSquared(w, h, v, channelCount, rank, sampleCount);
            history.Add(currentObjective);

            // Relative-improvement stop. The objective is monotonically non-increasing (Theorems 1 & 2), so
            // previousObjective ≥ currentObjective; the decrease is non-negative.
            double decrease = previousObjective - currentObjective;
            double denominator = previousObjective > NmfEpsilon ? previousObjective : 1.0;
            if (!double.IsInfinity(previousObjective) && decrease / denominator < tolerance)
                break;

            previousObjective = currentObjective;
        }

        // L1-normalise each signature column of W so its channel weights sum to 1, absorbing the scale into the
        // corresponding row of H (Alexandrov et al. 2013/2020; COSMIC SBS — signatures are probability
        // distributions over the channels). This fixes the NMF scaling ambiguity without changing W·H.
        NormalizeSignatureColumns(w, h, channelCount, rank, sampleCount);

        double finalObjective = objective == NmfObjective.KullbackLeibler
            ? KullbackLeiblerDivergence(w, h, v, channelCount, rank, sampleCount)
            : FrobeniusResidualSquared(w, h, v, channelCount, rank, sampleCount);

        return (w, h, finalObjective, iterations, history);
    }

    /// <summary>
    /// Validates the count matrix V (non-null, non-empty, rectangular, finite, non-negative) and returns it as a
    /// dense jagged array of rows, with the channel and sample counts.
    /// </summary>
    private static double[][] ValidateCountMatrix(
        IReadOnlyList<IReadOnlyList<double>> countMatrix, out int channelCount, out int sampleCount)
    {
        ArgumentNullException.ThrowIfNull(countMatrix);

        channelCount = countMatrix.Count;
        if (channelCount == 0)
        {
            throw new ArgumentException("The count matrix must have at least one channel (row).", nameof(countMatrix));
        }

        IReadOnlyList<double> firstRow = countMatrix[0]
            ?? throw new ArgumentException("Count-matrix rows cannot be null.", nameof(countMatrix));
        sampleCount = firstRow.Count;
        if (sampleCount == 0)
        {
            throw new ArgumentException("The count matrix must have at least one sample (column).", nameof(countMatrix));
        }

        var v = new double[channelCount][];
        for (int c = 0; c < channelCount; c++)
        {
            IReadOnlyList<double> row = countMatrix[c]
                ?? throw new ArgumentException("Count-matrix rows cannot be null.", nameof(countMatrix));
            if (row.Count != sampleCount)
            {
                throw new ArgumentException(
                    $"All channel rows must have the same sample count (row 0 has {sampleCount}, " +
                    $"row {c} has {row.Count}).",
                    nameof(countMatrix));
            }

            var dense = new double[sampleCount];
            for (int s = 0; s < sampleCount; s++)
            {
                double value = row[s];
                if (double.IsNaN(value) || double.IsInfinity(value) || value < 0)
                {
                    throw new ArgumentException(
                        $"Count-matrix entries must be finite and ≥ 0 (row {c}, sample {s} = {value}).",
                        nameof(countMatrix));
                }

                dense[s] = value;
            }

            v[c] = dense;
        }

        return v;
    }

    /// <summary>
    /// Builds a non-negative factor of the given shape with entries drawn uniformly from (0, 1], floored by
    /// <see cref="NmfEpsilon"/> so no entry is exactly zero (a zero row/column cannot recover under the
    /// multiplicative updates). Source: standard non-negative random initialisation for Lee &amp; Seung NMF.
    /// </summary>
    private static double[][] InitializeNonNegativeFactor(Random rng, int rows, int cols)
    {
        var factor = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            var rowValues = new double[cols];
            for (int j = 0; j < cols; j++)
            {
                rowValues[j] = rng.NextDouble() + NmfEpsilon;
            }

            factor[i] = rowValues;
        }

        return factor;
    }

    /// <summary>
    /// Multiplicative update of H: <c>H ← H ⊙ (Wᵀ V) ⊘ (Wᵀ W H)</c>. Source: Lee &amp; Seung (2001), Theorem 1
    /// (Euclidean objective). The denominator is floored by <see cref="NmfEpsilon"/> to avoid 0/0.
    /// </summary>
    private static void UpdateH(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        // numerator = Wᵀ V  (k × samples); denominator = Wᵀ (W H)  (k × samples).
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        for (int a = 0; a < k; a++)
        {
            for (int s = 0; s < samples; s++)
            {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int c = 0; c < channels; c++)
                {
                    numerator += w[c][a] * v[c][s];
                    denominator += w[c][a] * wh[c][s];
                }

                if (denominator < NmfEpsilon)
                {
                    denominator = NmfEpsilon;
                }

                h[a][s] *= numerator / denominator;
            }
        }
    }

    /// <summary>
    /// Multiplicative update of W: <c>W ← W ⊙ (V Hᵀ) ⊘ (W H Hᵀ)</c>. Source: Lee &amp; Seung (2001), Theorem 1
    /// (Euclidean objective). The denominator is floored by <see cref="NmfEpsilon"/> to avoid 0/0.
    /// </summary>
    private static void UpdateW(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        // numerator = V Hᵀ  (channels × k); denominator = (W H) Hᵀ  (channels × k).
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        for (int c = 0; c < channels; c++)
        {
            for (int a = 0; a < k; a++)
            {
                double numerator = 0.0;
                double denominator = 0.0;
                for (int s = 0; s < samples; s++)
                {
                    numerator += v[c][s] * h[a][s];
                    denominator += wh[c][s] * h[a][s];
                }

                if (denominator < NmfEpsilon)
                {
                    denominator = NmfEpsilon;
                }

                w[c][a] *= numerator / denominator;
            }
        }
    }

    /// <summary>
    /// Multiplicative KL/Poisson update of H:
    /// <c>H_aμ ← H_aμ · (Σ_i W_ia V_iμ/(WH)_iμ) / (Σ_i W_ia)</c>. Source: Lee &amp; Seung (2001), Theorem 2
    /// (divergence objective). The reconstruction (WH)_iμ is floored by <see cref="NmfEpsilon"/> to avoid V/0,
    /// and the column-sum denominator Σ_i W_ia is floored to avoid 0/0.
    /// </summary>
    private static void UpdateHKl(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        for (int a = 0; a < k; a++)
        {
            // Denominator Σ_i W_ia is independent of the sample μ.
            double columnSum = 0.0;
            for (int c = 0; c < channels; c++)
            {
                columnSum += w[c][a];
            }

            if (columnSum < NmfEpsilon)
            {
                columnSum = NmfEpsilon;
            }

            for (int s = 0; s < samples; s++)
            {
                double numerator = 0.0;
                for (int c = 0; c < channels; c++)
                {
                    double reconstructed = wh[c][s] < NmfEpsilon ? NmfEpsilon : wh[c][s];
                    numerator += w[c][a] * v[c][s] / reconstructed;
                }

                h[a][s] *= numerator / columnSum;
            }
        }
    }

    /// <summary>
    /// Multiplicative KL/Poisson update of W:
    /// <c>W_ia ← W_ia · (Σ_μ H_aμ V_iμ/(WH)_iμ) / (Σ_μ H_aμ)</c>. Source: Lee &amp; Seung (2001), Theorem 2
    /// (divergence objective). The reconstruction (WH)_iμ is floored by <see cref="NmfEpsilon"/> to avoid V/0,
    /// and the row-sum denominator Σ_μ H_aμ is floored to avoid 0/0.
    /// </summary>
    private static void UpdateWKl(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        double[][] wh = MultiplyWh(w, h, channels, k, samples);

        // Denominator Σ_μ H_aμ is independent of the channel i.
        var rowSum = new double[k];
        for (int a = 0; a < k; a++)
        {
            double sum = 0.0;
            for (int s = 0; s < samples; s++)
            {
                sum += h[a][s];
            }

            rowSum[a] = sum < NmfEpsilon ? NmfEpsilon : sum;
        }

        for (int c = 0; c < channels; c++)
        {
            for (int a = 0; a < k; a++)
            {
                double numerator = 0.0;
                for (int s = 0; s < samples; s++)
                {
                    double reconstructed = wh[c][s] < NmfEpsilon ? NmfEpsilon : wh[c][s];
                    numerator += h[a][s] * v[c][s] / reconstructed;
                }

                w[c][a] *= numerator / rowSum[a];
            }
        }
    }

    /// <summary>
    /// The generalized Kullback-Leibler divergence <c>D(V‖WH) = Σ_iμ ( V_iμ log(V_iμ/(WH)_iμ) − V_iμ + (WH)_iμ )</c>.
    /// Source: Lee &amp; Seung (2001), the divergence objective. The term V·log(V/WH) is taken as 0 when V = 0
    /// (limit x log x → 0); (WH)_iμ is floored by <see cref="NmfEpsilon"/> to avoid log(·/0).
    /// </summary>
    private static double KullbackLeiblerDivergence(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        double sum = 0.0;
        for (int c = 0; c < channels; c++)
        {
            for (int s = 0; s < samples; s++)
            {
                double reconstructed = 0.0;
                for (int a = 0; a < k; a++)
                {
                    reconstructed += w[c][a] * h[a][s];
                }

                if (reconstructed < NmfEpsilon)
                {
                    reconstructed = NmfEpsilon;
                }

                double observed = v[c][s];
                if (observed > 0.0)
                {
                    sum += observed * Math.Log(observed / reconstructed) - observed + reconstructed;
                }
                else
                {
                    // V = 0: the log term vanishes (limit x log x → 0); only +(WH) remains.
                    sum += reconstructed;
                }
            }
        }

        return sum;
    }

    /// <summary>
    /// Computes the dense product W·H (channels × samples) of the current factors.
    /// </summary>
    private static double[][] MultiplyWh(double[][] w, double[][] h, int channels, int k, int samples)
    {
        var wh = new double[channels][];
        for (int c = 0; c < channels; c++)
        {
            var row = new double[samples];
            for (int s = 0; s < samples; s++)
            {
                double sum = 0.0;
                for (int a = 0; a < k; a++)
                {
                    sum += w[c][a] * h[a][s];
                }

                row[s] = sum;
            }

            wh[c] = row;
        }

        return wh;
    }

    /// <summary>
    /// The squared Frobenius reconstruction residual ‖V − W·H‖²_F = Σ (V − W·H)². Source: Lee &amp; Seung (2001),
    /// the Euclidean objective F(W,H) = ‖V − WH‖²_F.
    /// </summary>
    private static double FrobeniusResidualSquared(double[][] w, double[][] h, double[][] v, int channels, int k, int samples)
    {
        double sum = 0.0;
        for (int c = 0; c < channels; c++)
        {
            for (int s = 0; s < samples; s++)
            {
                double reconstructed = 0.0;
                for (int a = 0; a < k; a++)
                {
                    reconstructed += w[c][a] * h[a][s];
                }

                double diff = v[c][s] - reconstructed;
                sum += diff * diff;
            }
        }

        return sum;
    }

    /// <summary>
    /// L1-normalises each signature column of W so its channel weights sum to 1, absorbing the removed scale into
    /// the matching row of H so that W·H is unchanged. Source: Alexandrov et al. (2013/2020); COSMIC SBS — each
    /// signature is a probability distribution over the channels. A column that sums to zero is left as zeros.
    /// </summary>
    private static void NormalizeSignatureColumns(double[][] w, double[][] h, int channels, int k, int samples)
    {
        for (int a = 0; a < k; a++)
        {
            double columnSum = 0.0;
            for (int c = 0; c < channels; c++)
            {
                columnSum += w[c][a];
            }

            if (columnSum <= NmfEpsilon)
            {
                continue;
            }

            for (int c = 0; c < channels; c++)
            {
                w[c][a] /= columnSum;
            }

            // Absorb the scale into H so W·H is invariant: H[a][·] *= columnSum.
            for (int s = 0; s < samples; s++)
            {
                h[a][s] *= columnSum;
            }
        }
    }

    /// <summary>
    /// Converts the channels × k factor W into k signature vectors (column a of W → signature a of length
    /// <paramref name="channels"/>), the per-signature channel-vector layout used elsewhere in this class.
    /// </summary>
    private static IReadOnlyList<IReadOnlyList<double>> TransposeColumnsToSignatures(double[][] w, int channels, int k)
    {
        var signatures = new List<IReadOnlyList<double>>(k);
        for (int a = 0; a < k; a++)
        {
            var signature = new double[channels];
            for (int c = 0; c < channels; c++)
            {
                signature[c] = w[c][a];
            }

            signatures.Add(signature);
        }

        return signatures;
    }

    /// <summary>
    /// Wraps the rows of H as a read-only list of read-only rows.
    /// </summary>
    private static IReadOnlyList<IReadOnlyList<double>> RowsToReadOnly(double[][] h)
    {
        var rows = new List<IReadOnlyList<double>>(h.Length);
        foreach (double[] row in h)
        {
            rows.Add(row);
        }

        return rows;
    }

    // ---- Automatic rank / model-stability selection (Brunet 2004; Alexandrov 2013 / SigProfiler) ----

    /// <summary>
    /// Default number of independent NMF runs (random restarts) per candidate rank used by
    /// <see cref="SelectRank"/>. Source: Brunet et al. (2004) build the consensus matrix as the average of the
    /// connectivity matrices over "multiple runs"; SigProfilerExtractor uses partition clustering of replicate
    /// factorizations (Islam et al. 2022). A modest default keeps deterministic tests fast; callers raise it for
    /// production. Value = 20.
    /// </summary>
    public const int DefaultRankSelectionRuns = 20;

    /// <summary>
    /// Default average-stability acceptance threshold for <see cref="SelectRank"/>. Source: SigProfilerExtractor
    /// considers a solution stable when "signatures have an average stability above 0.80" (Islam et al. 2022;
    /// default <c>stability=0.8</c>). Value = 0.80.
    /// </summary>
    public const double DefaultRankStabilityThreshold = 0.80;

    /// <summary>
    /// Per-candidate-rank diagnostics produced by <see cref="SelectRank"/>, making the rank choice auditable
    /// (Brunet et al. 2004; Alexandrov et al. 2013).
    /// </summary>
    /// <param name="Rank">The candidate rank k.</param>
    /// <param name="CopheneticCorrelation">
    /// The cophenetic correlation coefficient of the consensus matrix over the runs — the Pearson correlation
    /// between the consensus-induced sample distances and their cophenetic distances from average-linkage
    /// hierarchical clustering (Brunet et al. 2004). 1.0 = perfectly stable clustering; the first rank where it
    /// "begins to fall" is the Brunet rank estimate.
    /// </param>
    /// <param name="AverageStability">
    /// The mean over signatures of the per-signature average silhouette width across the runs (cosine distance
    /// in signature space). Source: Alexandrov et al. (2013) / SigProfilerExtractor "silhouette value of the
    /// cluster corresponding to that signature" (Islam et al. 2022); silhouette per Rousseeuw (1987).
    /// </param>
    /// <param name="MinimumStability">The minimum per-signature silhouette stability over the k clusters.</param>
    /// <param name="MeanReconstructionError">
    /// The mean over runs of the final reconstruction error (selected objective; Frobenius residual ‖V−WH‖²_F or
    /// KL divergence) — lower is a better data fit.
    /// </param>
    public readonly record struct RankStability(
        int Rank,
        double CopheneticCorrelation,
        double AverageStability,
        double MinimumStability,
        double MeanReconstructionError);

    /// <summary>
    /// The result of automatic NMF rank selection over a candidate range (Brunet et al. 2004; Alexandrov et al.
    /// 2013 / SigProfilerExtractor, Islam et al. 2022).
    /// </summary>
    /// <param name="SelectedRank">
    /// The chosen number of signatures: the largest candidate rank whose average signature stability is at least
    /// <c>stabilityThreshold</c> and whose minimum per-signature stability is at least <c>minStability</c>
    /// (SigProfiler: "average stability above 0.80 with no individual signature having stability below 0.20").
    /// If no rank meets the threshold, the rank with the highest average stability is returned.
    /// </param>
    /// <param name="PerRank">Per-candidate-rank diagnostics (cophenetic, stability, error), in ascending rank order.</param>
    public readonly record struct RankSelectionResult(
        int SelectedRank,
        IReadOnlyList<RankStability> PerRank);

    /// <summary>
    /// Automatically selects the number of de-novo signatures (NMF rank) for a mutation-count matrix by the
    /// SigProfiler/Brunet model-stability approach. For each candidate rank k in <paramref name="minRank"/>..
    /// <paramref name="maxRank"/>, NMF is run <paramref name="runs"/> times from a <b>fixed, deterministic seed
    /// sequence</b>; the runs are summarised by (a) the cophenetic correlation coefficient of the consensus
    /// matrix (Brunet et al. 2004), (b) the average signature <em>stability</em> measured as the per-signature
    /// average silhouette width across runs (Alexandrov et al. 2013; Islam et al. 2022), and (c) the mean
    /// reconstruction error. The selected rank is the largest k whose average stability is ≥
    /// <paramref name="stabilityThreshold"/> and whose minimum per-signature stability is ≥
    /// <paramref name="minStability"/> (SigProfiler default 0.80 / 0.20); if none qualifies, the highest-average-
    /// stability rank is returned. All per-rank diagnostics are returned so the choice is auditable.
    /// </summary>
    /// <param name="countMatrix">Mutation-count matrix V (channels × samples); rows are channels.</param>
    /// <param name="minRank">Smallest candidate rank (≥ 1).</param>
    /// <param name="maxRank">Largest candidate rank (≥ <paramref name="minRank"/>, ≤ channel count).</param>
    /// <param name="objective">NMF objective for the per-run extractions.</param>
    /// <param name="runs">Independent NMF runs per candidate rank (≥ 1).</param>
    /// <param name="stabilityThreshold">Minimum acceptable average signature stability (default 0.80).</param>
    /// <param name="minStability">Minimum acceptable per-signature stability (default 0.20).</param>
    /// <param name="maxIterations">Maximum multiplicative-update iterations per run (&gt; 0).</param>
    /// <param name="tolerance">Relative-improvement convergence tolerance (≥ 0).</param>
    /// <param name="seed">Base RNG seed; run r at rank k uses a deterministic derived seed for reproducibility.</param>
    /// <returns>The selected rank plus per-rank stability/error diagnostics.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="countMatrix"/> or a row is null.</exception>
    /// <exception cref="ArgumentException">
    /// Invalid matrix; minRank &lt; 1; maxRank &lt; minRank; maxRank &gt; channel count; runs &lt; 1; bad iteration
    /// / tolerance; or a stability threshold outside [0, 1].
    /// </exception>
    public static RankSelectionResult SelectRank(
        IReadOnlyList<IReadOnlyList<double>> countMatrix,
        int minRank,
        int maxRank,
        NmfObjective objective = NmfObjective.KullbackLeibler,
        int runs = DefaultRankSelectionRuns,
        double stabilityThreshold = DefaultRankStabilityThreshold,
        double minStability = DefaultMinSignatureStability,
        int maxIterations = DefaultNmfMaxIterations,
        double tolerance = DefaultNmfTolerance,
        int seed = DefaultNmfSeed)
    {
        double[][] v = ValidateCountMatrix(countMatrix, out int channelCount, out int sampleCount);

        if (minRank < 1)
        {
            throw new ArgumentException($"minRank must be ≥ 1 (got {minRank}).", nameof(minRank));
        }

        if (maxRank < minRank)
        {
            throw new ArgumentException(
                $"maxRank ({maxRank}) must be ≥ minRank ({minRank}).", nameof(maxRank));
        }

        if (maxRank > channelCount)
        {
            throw new ArgumentException(
                $"maxRank ({maxRank}) cannot exceed the channel count ({channelCount}).", nameof(maxRank));
        }

        if (runs < 1)
        {
            throw new ArgumentException($"runs must be ≥ 1 (got {runs}).", nameof(runs));
        }

        if (stabilityThreshold < 0.0 || stabilityThreshold > 1.0)
        {
            throw new ArgumentException(
                $"stabilityThreshold must be in [0, 1] (got {stabilityThreshold}).", nameof(stabilityThreshold));
        }

        if (minStability < 0.0 || minStability > 1.0)
        {
            throw new ArgumentException(
                $"minStability must be in [0, 1] (got {minStability}).", nameof(minStability));
        }

        ValidateRankAndStop(minRank, channelCount, maxIterations, tolerance);

        var perRank = new List<RankStability>(maxRank - minRank + 1);
        for (int k = minRank; k <= maxRank; k++)
        {
            // Collect the signature sets and per-sample cluster assignments of all runs at this rank.
            var runSignatures = new List<double[][]>(runs); // each: k signatures × channels
            var connectivity = new double[sampleCount][];   // accumulated consensus (sum of connectivity)
            for (int s = 0; s < sampleCount; s++)
            {
                connectivity[s] = new double[sampleCount];
            }

            double errorSum = 0.0;
            for (int r = 0; r < runs; r++)
            {
                // Deterministic per-(rank, run) seed so the whole procedure is reproducible.
                int runSeed = DerivedSeed(seed, k, r);
                (double[][] w, double[][] h, double finalObjective, _, _) =
                    RunNmf(v, k, objective, maxIterations, tolerance, runSeed, channelCount, sampleCount);

                errorSum += finalObjective;

                // Signatures of this run as k channel-vectors (column a of W).
                var sigs = new double[k][];
                for (int a = 0; a < k; a++)
                {
                    sigs[a] = new double[channelCount];
                    for (int c = 0; c < channelCount; c++)
                    {
                        sigs[a][c] = w[c][a];
                    }
                }

                runSignatures.Add(sigs);

                // Connectivity: samples assigned to the same metagene (argmax over H column) are connected.
                int[] assignment = AssignSamplesToClusters(h, k, sampleCount);
                for (int i = 0; i < sampleCount; i++)
                {
                    for (int j = 0; j < sampleCount; j++)
                    {
                        if (assignment[i] == assignment[j])
                        {
                            connectivity[i][j] += 1.0;
                        }
                    }
                }
            }

            // Consensus matrix = average connectivity over runs (Brunet 2004).
            for (int i = 0; i < sampleCount; i++)
            {
                for (int j = 0; j < sampleCount; j++)
                {
                    connectivity[i][j] /= runs;
                }
            }

            double cophenetic = CopheneticCorrelation(connectivity, sampleCount);
            (double avgStability, double minStab) = SignatureStability(runSignatures, k, channelCount);
            double meanError = errorSum / runs;

            perRank.Add(new RankStability(k, cophenetic, avgStability, minStab, meanError));
        }

        int selected = ChooseRank(perRank, stabilityThreshold, minStability);
        return new RankSelectionResult(selected, perRank);
    }

    /// <summary>
    /// Default minimum acceptable per-signature stability for <see cref="SelectRank"/>. Source:
    /// SigProfilerExtractor requires "no individual signature having stability below 0.20" (Islam et al. 2022;
    /// default <c>min_stability=0.2</c>). Value = 0.20.
    /// </summary>
    public const double DefaultMinSignatureStability = 0.20;

    /// <summary>
    /// Deterministic per-(rank, run) seed derivation so the whole rank-selection procedure is reproducible while
    /// each run still gets a distinct initialisation. Combines the base seed, rank and run index.
    /// </summary>
    private static int DerivedSeed(int baseSeed, int rank, int run)
    {
        unchecked
        {
            int hash = baseSeed;
            hash = (hash * 397) ^ rank;
            hash = (hash * 397) ^ run;
            return hash;
        }
    }

    /// <summary>
    /// Assigns each sample (column of H) to the cluster of its largest exposure (metagene), per Brunet et al.
    /// (2004): "Sample assignment is determined by its largest metagene expression value."
    /// </summary>
    private static int[] AssignSamplesToClusters(double[][] h, int k, int samples)
    {
        var assignment = new int[samples];
        for (int s = 0; s < samples; s++)
        {
            int best = 0;
            double bestValue = h[0][s];
            for (int a = 1; a < k; a++)
            {
                if (h[a][s] > bestValue)
                {
                    bestValue = h[a][s];
                    best = a;
                }
            }

            assignment[s] = best;
        }

        return assignment;
    }

    /// <summary>
    /// The cophenetic correlation coefficient of the consensus matrix (Brunet et al. 2004): the Pearson
    /// correlation between the sample distances induced by the consensus matrix (distance = 1 − consensus) and
    /// the cophenetic distances from average-linkage hierarchical clustering on those distances. Returns 1.0 when
    /// every consensus entry equals 1 (perfectly stable clustering — e.g. rank 1), the value defined by Brunet.
    /// </summary>
    private static double CopheneticCorrelation(double[][] consensus, int samples)
    {
        if (samples < 2)
        {
            return 1.0;
        }

        // Distance induced by the consensus matrix: d = 1 − consensus (consensus is a similarity in [0, 1]).
        var distance = new double[samples][];
        for (int i = 0; i < samples; i++)
        {
            distance[i] = new double[samples];
            for (int j = 0; j < samples; j++)
            {
                distance[i][j] = 1.0 - consensus[i][j];
            }
        }

        double[][] cophenetic = AverageLinkageCopheneticDistances(distance, samples);

        // Pearson correlation over the strictly-upper-triangular pairs (i < j).
        double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0, sumYY = 0.0;
        int n = 0;
        for (int i = 0; i < samples; i++)
        {
            for (int j = i + 1; j < samples; j++)
            {
                double x = distance[i][j];
                double y = cophenetic[i][j];
                sumX += x;
                sumY += y;
                sumXY += x * y;
                sumXX += x * x;
                sumYY += y * y;
                n++;
            }
        }

        double meanX = sumX / n;
        double meanY = sumY / n;
        double covXY = sumXY / n - meanX * meanY;
        double varX = sumXX / n - meanX * meanX;
        double varY = sumYY / n - meanY * meanY;

        // Degenerate: a constant distance vector (e.g. all-ones consensus ⇒ all-zero distance) has zero
        // variance; Brunet's "perfect consensus" case maps to cophenetic = 1.
        if (varX <= NmfEpsilon || varY <= NmfEpsilon)
        {
            return 1.0;
        }

        return covXY / Math.Sqrt(varX * varY);
    }

    /// <summary>
    /// Computes the cophenetic distance matrix from average-linkage (UPGMA) agglomerative clustering of the given
    /// distance matrix. The cophenetic distance between two samples is the linkage distance at which their
    /// clusters first merge (Brunet et al. 2004 use average linkage by default).
    /// </summary>
    private static double[][] AverageLinkageCopheneticDistances(double[][] distance, int samples)
    {
        // Active clusters: each starts as a singleton {i}. clusterMembers holds the sample indices per cluster.
        var members = new List<List<int>>(samples);
        var active = new List<int>(samples);
        for (int i = 0; i < samples; i++)
        {
            members.Add(new List<int> { i });
            active.Add(i);
        }

        // Pairwise average-linkage distances between active clusters, indexed by cluster id.
        var clusterDistance = new Dictionary<(int, int), double>();
        for (int a = 0; a < active.Count; a++)
        {
            for (int b = a + 1; b < active.Count; b++)
            {
                clusterDistance[(active[a], active[b])] = distance[active[a]][active[b]];
            }
        }

        var cophenetic = new double[samples][];
        for (int i = 0; i < samples; i++)
        {
            cophenetic[i] = new double[samples];
        }

        int nextId = samples;
        while (active.Count > 1)
        {
            // Find the closest pair of active clusters.
            double best = double.PositiveInfinity;
            int bi = 0, bj = 1;
            for (int a = 0; a < active.Count; a++)
            {
                for (int b = a + 1; b < active.Count; b++)
                {
                    double d = clusterDistance[Key(active[a], active[b])];
                    if (d < best)
                    {
                        best = d;
                        bi = a;
                        bj = b;
                    }
                }
            }

            int ci = active[bi];
            int cj = active[bj];

            // The merge height 'best' is the cophenetic distance between every member of ci and every member of cj.
            foreach (int p in members[ci])
            {
                foreach (int q in members[cj])
                {
                    cophenetic[p][q] = best;
                    cophenetic[q][p] = best;
                }
            }

            // Merge cj into a new cluster id; average-linkage distance to every other active cluster.
            var merged = new List<int>(members[ci]);
            merged.AddRange(members[cj]);
            members.Add(merged);
            int newId = nextId++;

            foreach (int other in active)
            {
                if (other == ci || other == cj)
                {
                    continue;
                }

                // Average linkage: mean over all member-pairs (size-weighted average of the two sub-distances).
                double dCi = clusterDistance[Key(ci, other)];
                double dCj = clusterDistance[Key(cj, other)];
                double weighted = (members[ci].Count * dCi + members[cj].Count * dCj)
                                  / (members[ci].Count + members[cj].Count);
                clusterDistance[Key(newId, other)] = weighted;
            }

            active.RemoveAt(bj);
            active.RemoveAt(bi);
            active.Add(newId);
        }

        return cophenetic;

        static (int, int) Key(int a, int b) => a < b ? (a, b) : (b, a);
    }

    /// <summary>
    /// Computes signature stability across the NMF runs as the per-signature average silhouette width
    /// (Alexandrov et al. 2013; Islam et al. 2022 "silhouette value of the cluster corresponding to that
    /// signature"; silhouette per Rousseeuw 1987), clustering the <c>runs × k</c> extracted signatures into k
    /// clusters by greedy cosine matching against the first run's signatures as the reference partition. Distance
    /// is cosine distance (1 − cosine similarity). Returns (mean over signatures of average silhouette,
    /// minimum over signatures). Each signature is L1-normalised, so cosine distance is well defined.
    /// </summary>
    private static (double Average, double Minimum) SignatureStability(
        IReadOnlyList<double[][]> runSignatures, int k, int channels)
    {
        int runs = runSignatures.Count;
        if (runs < 2 || k < 1)
        {
            // With a single run there is no cross-run dispersion to measure; treat as perfectly stable.
            return (1.0, 1.0);
        }

        // Flatten all signatures and assign each to the cluster (reference signature of run 0) it is closest to
        // by cosine, requiring a one-to-one match within each run (Hungarian reduced to greedy per run).
        int total = runs * k;
        var points = new double[total][];
        var cluster = new int[total];
        var reference = runSignatures[0];

        int idx = 0;
        for (int r = 0; r < runs; r++)
        {
            double[][] sigs = runSignatures[r];
            bool[] taken = new bool[k];
            // Greedily assign each signature of this run to its best unused reference cluster.
            for (int a = 0; a < k; a++)
            {
                points[idx] = sigs[a];
                int bestCluster = -1;
                double bestSim = double.NegativeInfinity;
                for (int c = 0; c < k; c++)
                {
                    if (taken[c])
                    {
                        continue;
                    }

                    double sim = CosineSimilarity(sigs[a], reference[c]);
                    if (sim > bestSim)
                    {
                        bestSim = sim;
                        bestCluster = c;
                    }
                }

                taken[bestCluster] = true;
                cluster[idx] = bestCluster;
                idx++;
            }
        }

        // Per-point silhouette with cosine distance.
        var clusterSilhouetteSum = new double[k];
        var clusterSize = new int[k];
        for (int i = 0; i < total; i++)
        {
            clusterSize[cluster[i]]++;
        }

        for (int i = 0; i < total; i++)
        {
            int ci = cluster[i];

            // a(i): mean intra-cluster distance.
            double aSum = 0.0;
            int aCount = 0;
            // b(i): min over other clusters of mean distance to that cluster.
            var otherSum = new double[k];
            var otherCount = new int[k];
            for (int j = 0; j < total; j++)
            {
                if (j == i)
                {
                    continue;
                }

                double d = 1.0 - CosineSimilarity(points[i], points[j]);
                if (cluster[j] == ci)
                {
                    aSum += d;
                    aCount++;
                }
                else
                {
                    otherSum[cluster[j]] += d;
                    otherCount[cluster[j]]++;
                }
            }

            double a = aCount > 0 ? aSum / aCount : 0.0;
            double b = double.PositiveInfinity;
            for (int c = 0; c < k; c++)
            {
                if (c == ci || otherCount[c] == 0)
                {
                    continue;
                }

                double mean = otherSum[c] / otherCount[c];
                if (mean < b)
                {
                    b = mean;
                }
            }

            double silhouette;
            if (double.IsPositiveInfinity(b))
            {
                // Only one cluster present: silhouette is conventionally 0 (Rousseeuw 1987).
                silhouette = 0.0;
            }
            else
            {
                double denom = Math.Max(a, b);
                silhouette = denom <= NmfEpsilon ? 0.0 : (b - a) / denom;
            }

            clusterSilhouetteSum[ci] += silhouette;
        }

        double avg = 0.0;
        double min = double.PositiveInfinity;
        int nonEmpty = 0;
        for (int c = 0; c < k; c++)
        {
            if (clusterSize[c] == 0)
            {
                continue;
            }

            double clusterStability = clusterSilhouetteSum[c] / clusterSize[c];
            avg += clusterStability;
            min = Math.Min(min, clusterStability);
            nonEmpty++;
        }

        avg = nonEmpty > 0 ? avg / nonEmpty : 0.0;
        if (double.IsPositiveInfinity(min))
        {
            min = 0.0;
        }

        return (avg, min);
    }

    /// <summary>
    /// Chooses the rank: the largest candidate whose average stability ≥ <paramref name="stabilityThreshold"/>
    /// and minimum per-signature stability ≥ <paramref name="minStability"/> (SigProfiler 0.80 / 0.20); if none
    /// qualifies, the rank with the highest average stability (Alexandrov 2013 stability/error trade-off).
    /// </summary>
    private static int ChooseRank(
        IReadOnlyList<RankStability> perRank, double stabilityThreshold, double minStability)
    {
        int selected = -1;
        for (int i = 0; i < perRank.Count; i++)
        {
            RankStability rs = perRank[i];
            if (rs.AverageStability >= stabilityThreshold && rs.MinimumStability >= minStability)
            {
                selected = rs.Rank; // keep the largest qualifying rank
            }
        }

        if (selected >= 0)
        {
            return selected;
        }

        // Fallback: highest average stability (ties → smallest rank, the most parsimonious model).
        RankStability best = perRank[0];
        for (int i = 1; i < perRank.Count; i++)
        {
            if (perRank[i].AverageStability > best.AverageStability)
            {
                best = perRank[i];
            }
        }

        return best.Rank;
    }

    // ---- COSMIC / reference matching by cosine similarity (Alexandrov 2013/2020; Islam et al. 2022) ----

    /// <summary>
    /// The best-matching reference signature for one extracted signature, by cosine similarity
    /// (Alexandrov et al. 2013/2020; Islam et al. 2022 — de-novo signatures are matched to references by
    /// maximising cosine similarity).
    /// </summary>
    /// <param name="ExtractedIndex">Index of the extracted signature this match describes.</param>
    /// <param name="ReferenceIndex">Index of the closest reference signature.</param>
    /// <param name="CosineSimilarity">
    /// The cosine similarity in [0, 1] between the extracted signature and the closest reference (1 = identical
    /// direction). Cosine is scale-invariant, so a positively-scaled copy of a reference matches it with 1.0.
    /// </param>
    public readonly record struct SignatureMatch(
        int ExtractedIndex,
        int ReferenceIndex,
        double CosineSimilarity);

    /// <summary>
    /// Matches each extracted (de-novo) signature to its closest reference signature by cosine similarity,
    /// labelling each extracted signature with the index of the best-matching reference and the cosine value
    /// (Alexandrov et al. 2013/2020; Islam et al. 2022). The reference set is <b>caller-supplied</b> (e.g. COSMIC
    /// SBS profiles) — no profiles are hardcoded. This is the per-signature reduction of SigProfiler's Hungarian
    /// "maximise total cosine" pairing: each extracted signature is labelled with its single nearest reference.
    /// </summary>
    /// <param name="extractedSignatures">
    /// The de-novo signatures to label (e.g. <see cref="SignatureExtractionResult.Signatures"/>); one channel
    /// vector per signature.
    /// </param>
    /// <param name="referenceSignatures">
    /// The reference signatures to match against (e.g. COSMIC SBS); one channel vector per reference, each of the
    /// same channel count as the extracted signatures.
    /// </param>
    /// <returns>One <see cref="SignatureMatch"/> per extracted signature, in extracted order.</returns>
    /// <exception cref="ArgumentNullException">Any argument or signature vector is null.</exception>
    /// <exception cref="ArgumentException">
    /// Either set is empty, signatures are ragged, or the extracted and reference channel counts differ.
    /// </exception>
    public static IReadOnlyList<SignatureMatch> MatchToReferenceSignatures(
        IReadOnlyList<IReadOnlyList<double>> extractedSignatures,
        IReadOnlyList<IReadOnlyList<double>> referenceSignatures)
    {
        int extractedChannels = ValidateSignatures(extractedSignatures);
        int referenceChannels = ValidateSignatures(referenceSignatures);

        if (extractedChannels != referenceChannels)
        {
            throw new ArgumentException(
                $"Extracted signatures have {extractedChannels} channels but reference signatures have "
                + $"{referenceChannels}; they must match.", nameof(referenceSignatures));
        }

        var matches = new List<SignatureMatch>(extractedSignatures.Count);
        for (int e = 0; e < extractedSignatures.Count; e++)
        {
            IReadOnlyList<double> extracted = extractedSignatures[e];
            int bestRef = 0;
            double bestCosine = CosineSimilarity(extracted, referenceSignatures[0]);
            for (int r = 1; r < referenceSignatures.Count; r++)
            {
                double cosine = CosineSimilarity(extracted, referenceSignatures[r]);
                if (cosine > bestCosine)
                {
                    bestCosine = cosine;
                    bestRef = r;
                }
            }

            matches.Add(new SignatureMatch(e, bestRef, bestCosine));
        }

        return matches;
    }

    #endregion


    #region Signature Exposure Bootstrap Confidence Intervals (ONCO-SIG-003)

    /// <summary>
    /// Default number of bootstrap replicates for exposure confidence-interval estimation. Source:
    /// Senkin S. (2021), MSA — confidence intervals are derived from "1000 bootstrap variations"
    /// (https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/); Wang S. et al., sigminer
    /// <c>sig_fit_bootstrap</c> recommends "Bootstrap replicates &gt;= 100"
    /// (https://shixiangwang.github.io/sigminer-doc/sigfit.html). Default = 1000.
    /// </summary>
    public const int DefaultBootstrapReplicates = 1000;

    /// <summary>
    /// Default RNG seed for the multinomial resampling step. Fixed so that, for a given catalog,
    /// signatures, replicate count and confidence level, the returned intervals are reproducible
    /// (deterministic test/clinical re-runs). Mirrors the fixed-seed convention used by
    /// <see cref="Seqeron.Genomics.Phylogenetics"/> bootstrap. Value = 42.
    /// </summary>
    public const int DefaultBootstrapSeed = 42;

    /// <summary>
    /// Default two-sided confidence level for the percentile bootstrap interval. Source: Efron B. (1979),
    /// percentile method — a 95% interval is the [2.5%, 97.5%] percentiles of the bootstrap distribution
    /// (Senkin 2021, MSA: "95% confidence intervals … taking [2.5%, 97.5%] percentiles of the resulting
    /// bootstrap activities", https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/). Default = 0.95.
    /// </summary>
    public const double DefaultBootstrapConfidence = 0.95;

    /// <summary>
    /// Selects how each bootstrap replicate of the mutational catalog is resampled in
    /// <see cref="BootstrapExposures"/>. Both schemes are described by Senkin (2021), MSA
    /// (<i>BMC Bioinformatics</i> 22:540, https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/), which notes that
    /// "the conditional distribution of a vector of independent Poisson variables is equivalent to multinomial
    /// distribution" — the two differ only in whether the total mutational burden N is held fixed.
    /// </summary>
    public enum BootstrapResampling
    {
        /// <summary>
        /// Fixed-N multinomial resample: the observed total N = Σ catalog mutations is redistributed across
        /// channels with probabilities pₖ = catalogₖ / N, so every replicate has exactly N mutations.
        /// This is the sigminer <c>sig_fit_bootstrap</c> scheme,
        /// <c>sample(K, total, replace=TRUE, prob=catalog/sum(catalog))</c>
        /// (https://github.com/ShixiangWang/sigminer/blob/master/R/sig_fit_bootstrap.R). Default — preserves the
        /// historical behaviour of <see cref="BootstrapExposures"/> byte-for-byte.
        /// </summary>
        Multinomial = 0,

        /// <summary>
        /// Poisson resample (Senkin 2021, MSA): each channel count is drawn independently as
        /// Poisson(observedₖ), so "for any given mutation category … the distribution of bootstrapped mutation
        /// counts follows a Poisson distribution" and "the total mutational burden is no longer fixed". This is
        /// the Poisson noise variant of the MSA parametric bootstrap, where the variance of each channel equals
        /// its mean (the observed count).
        /// </summary>
        Poisson = 1,
    }

    /// <summary>
    /// A bootstrap confidence interval for one signature's exposure (activity), produced by
    /// <see cref="BootstrapExposures"/>.
    /// </summary>
    /// <param name="PointEstimate">
    /// The exposure for this signature from the NNLS fit of the <i>observed</i> (un-resampled) catalog —
    /// the point estimate the interval is centred on (Senkin 2021; Huang et al. 2018).
    /// </param>
    /// <param name="Mean">The mean of this signature's exposure across the bootstrap replicates.</param>
    /// <param name="Lower">
    /// Lower bound of the percentile interval — the (½(1−c))·100-th percentile of the replicate exposures
    /// (2.5th percentile for c = 0.95). Source: Efron (1979) percentile method; Senkin (2021).
    /// </param>
    /// <param name="Upper">
    /// Upper bound — the (1−½(1−c))·100-th percentile of the replicate exposures (97.5th for c = 0.95).
    /// </param>
    /// <param name="Confidence">The two-sided confidence level c the interval was computed at.</param>
    public readonly record struct ExposureConfidenceInterval(
        double PointEstimate,
        double Mean,
        double Lower,
        double Upper,
        double Confidence);

    /// <summary>
    /// Estimates per-signature exposure confidence intervals by the parametric bootstrap: the observed integer
    /// mutational catalog is repeatedly resampled, each resampled catalog is refit to the reference signatures
    /// by NNLS (<see cref="FitSignatures(IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>), and a
    /// two-sided percentile confidence interval is taken per signature from the resulting bootstrap exposure
    /// distribution. The point estimate is the NNLS exposure of the un-resampled observed catalog. The
    /// <paramref name="resampling"/> parameter selects the resampling scheme:
    /// <list type="bullet">
    /// <item><see cref="BootstrapResampling.Multinomial"/> (default) — each replicate is a draw of
    /// N = Σ catalog mutations from the multinomial distribution with per-channel probabilities
    /// pₖ = catalogₖ / N (fixed total N).</item>
    /// <item><see cref="BootstrapResampling.Poisson"/> — each channel count is drawn independently as
    /// Poisson(observedₖ); the total burden is no longer fixed (Senkin 2021, MSA Poisson variant).</item>
    /// </list>
    /// <para>
    /// Sources: Huang X., Wojtowicz D., Przytycka T.M. (2018), <i>Bioinformatics</i> 34(2):330–337 — bootstrap
    /// resampling of the mutation-count vector to assess decomposition confidence; Senkin S. (2021), MSA,
    /// <i>BMC Bioinformatics</i> 22:540 — "mutations are accumulated following Poisson distributions for each
    /// mutation class", "drawing counts from independent binomial distributions, so that the total mutational
    /// burden is no longer fixed … for any given mutation category … the distribution of bootstrapped mutation
    /// counts follows a Poisson distribution", "the conditional distribution of a vector of independent Poisson
    /// variables is equivalent to multinomial distribution", "95% confidence intervals … taking [2.5%, 97.5%]
    /// percentiles"; Wang S. et al., sigminer <c>sig_fit_bootstrap</c> — resample via
    /// <c>sample(K, total, replace=TRUE, prob=catalog/sum(catalog))</c> (a multinomial draw) then refit;
    /// Efron B. (1979) percentile interval. Reference signature profiles are caller-supplied (not fabricated).
    /// </para>
    /// </summary>
    /// <param name="catalog">
    /// The observed mutational catalog as non-negative integer per-channel mutation counts (e.g. a 96-channel
    /// SBS spectrum). Counts (not proportions) are required because each resample needs the per-channel totals.
    /// </param>
    /// <param name="signatures">Reference signatures as equal-length channel vectors (one per signature).</param>
    /// <param name="replicates">Number of bootstrap replicates (≥ 1; default <see cref="DefaultBootstrapReplicates"/>).</param>
    /// <param name="confidence">Two-sided confidence level in (0, 1); default <see cref="DefaultBootstrapConfidence"/>.</param>
    /// <param name="seed">RNG seed for the resampling; fixed value makes results reproducible.</param>
    /// <param name="resampling">
    /// The bootstrap resampling scheme; <see cref="BootstrapResampling.Multinomial"/> (the historical default,
    /// fixed N) or <see cref="BootstrapResampling.Poisson"/> (Senkin 2021 Poisson variant).
    /// </param>
    /// <returns>One <see cref="ExposureConfidenceInterval"/> per signature, in signature order.</returns>
    /// <exception cref="ArgumentNullException">Any argument (or a signature vector) is null.</exception>
    /// <exception cref="ArgumentException">
    /// No signatures, ragged signatures, catalog length ≠ channel count, or a negative catalog count.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// replicates &lt; 1, confidence outside (0, 1), or an unrecognised <paramref name="resampling"/> value.
    /// </exception>
    public static IReadOnlyList<ExposureConfidenceInterval> BootstrapExposures(
        IReadOnlyList<int> catalog,
        IReadOnlyList<IReadOnlyList<double>> signatures,
        int replicates = DefaultBootstrapReplicates,
        double confidence = DefaultBootstrapConfidence,
        int seed = DefaultBootstrapSeed,
        BootstrapResampling resampling = BootstrapResampling.Multinomial)
    {
        ArgumentNullException.ThrowIfNull(catalog);
        int channelCount = ValidateSignatures(signatures);

        if (catalog.Count != channelCount)
        {
            throw new ArgumentException(
                $"Catalog length ({catalog.Count}) must equal the signature channel count ({channelCount}).",
                nameof(catalog));
        }

        if (replicates < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(replicates), replicates, "At least one bootstrap replicate is required.");
        }

        if (double.IsNaN(confidence) || confidence <= 0.0 || confidence >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(confidence), confidence, "Confidence must be in the open interval (0, 1).");
        }

        if (resampling != BootstrapResampling.Multinomial && resampling != BootstrapResampling.Poisson)
        {
            throw new ArgumentOutOfRangeException(
                nameof(resampling), resampling, "Unrecognised bootstrap resampling scheme.");
        }

        // Observed counts and total N = Σ catalog (the multinomial sample size).
        long total = 0;
        var observed = new double[channelCount];
        for (int k = 0; k < channelCount; k++)
        {
            int count = catalog[k];
            if (count < 0)
            {
                throw new ArgumentException("Catalog counts cannot be negative.", nameof(catalog));
            }

            observed[k] = count;
            total += count;
        }

        int signatureCount = signatures.Count;

        // Point estimate: NNLS exposures of the observed (un-resampled) catalog (Senkin 2021; Huang 2018).
        IReadOnlyList<double> pointEstimate = FitSignatures(observed, signatures).Exposures;

        // Per-signature bootstrap exposure distributions.
        var replicateExposures = new double[signatureCount][];
        for (int j = 0; j < signatureCount; j++)
        {
            replicateExposures[j] = new double[replicates];
        }

        var random = new Random(seed);
        var resampled = new double[channelCount];
        for (int rep = 0; rep < replicates; rep++)
        {
            if (resampling == BootstrapResampling.Poisson)
            {
                PoissonResample(observed, random, resampled);
            }
            else
            {
                MultinomialResample(observed, total, random, resampled);
            }

            IReadOnlyList<double> exposures = FitSignatures(resampled, signatures).Exposures;
            for (int j = 0; j < signatureCount; j++)
            {
                replicateExposures[j][rep] = exposures[j];
            }
        }

        double lowerProbability = (1.0 - confidence) / 2.0;
        double upperProbability = 1.0 - lowerProbability;

        var intervals = new ExposureConfidenceInterval[signatureCount];
        for (int j = 0; j < signatureCount; j++)
        {
            double[] distribution = replicateExposures[j];
            double mean = Mean(distribution);
            double lower = Percentile(distribution, lowerProbability);
            double upper = Percentile(distribution, upperProbability);
            intervals[j] = new ExposureConfidenceInterval(
                pointEstimate[j], mean, lower, upper, confidence);
        }

        return intervals;
    }

    /// <summary>
    /// Draws one multinomial resample of the catalog: N = <paramref name="total"/> mutations are assigned to
    /// channels with probabilities pₖ = observedₖ / N, written into <paramref name="destination"/>. Implements
    /// the standard sequential conditional-binomial construction of a multinomial draw: channel k receives
    /// Binomial(remaining, pₖ / Σ_{i≥k} pᵢ) of the mutations not yet assigned. Source: Senkin (2021), MSA
    /// (multinomial/Poisson resampling); sigminer <c>sig_fit_bootstrap</c>
    /// (<c>sample(..., replace=TRUE, prob=catalog/sum(catalog))</c>). When N = 0 the resample is all zeros.
    /// </summary>
    private static void MultinomialResample(double[] observed, long total, Random random, double[] destination)
    {
        int channelCount = observed.Length;
        Array.Clear(destination, 0, channelCount);

        if (total == 0)
        {
            return;
        }

        long remaining = total;
        double remainingProbabilityMass = total; // Σ observed = N (un-normalized probability mass).

        for (int k = 0; k < channelCount && remaining > 0; k++)
        {
            double weight = observed[k];
            if (weight > 0.0)
            {
                // Conditional probability of falling in channel k among the not-yet-assigned channels.
                double p = remainingProbabilityMass > 0.0 ? weight / remainingProbabilityMass : 0.0;
                if (p >= 1.0)
                {
                    // Last channel with any mass: it takes all remaining draws.
                    destination[k] = remaining;
                    remaining = 0;
                }
                else
                {
                    long drawn = SampleBinomial(remaining, p, random);
                    destination[k] = drawn;
                    remaining -= drawn;
                }
            }

            remainingProbabilityMass -= weight;
        }
    }

    /// <summary>
    /// Draws one Poisson resample of the catalog: each channel k is independently assigned
    /// Poisson(observedₖ) mutations, written into <paramref name="destination"/>. Implements the Senkin (2021)
    /// MSA Poisson-noise variant — "for any given mutation category … the distribution of bootstrapped mutation
    /// counts follows a Poisson distribution" with mean equal to the observed count (variance = mean) — so the
    /// total mutational burden is not fixed across replicates. Source: Senkin S. (2021), MSA,
    /// <i>BMC Bioinformatics</i> 22:540 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8567580/). A channel whose
    /// observed count is 0 has mean 0 and therefore always resamples to 0.
    /// </summary>
    private static void PoissonResample(double[] observed, Random random, double[] destination)
    {
        int channelCount = observed.Length;
        for (int k = 0; k < channelCount; k++)
        {
            destination[k] = SamplePoisson(observed[k], random);
        }
    }

    /// <summary>
    /// Samples a Poisson(lambda) variate by Knuth's multiplication-of-uniforms algorithm: draw uniforms until
    /// their running product falls below e^(−lambda); the number of draws minus one is the variate. Source:
    /// Knuth D.E., <i>The Art of Computer Programming</i>, Vol. 2 (Seminumerical Algorithms), §3.4.1 — the
    /// standard exact algorithm for generating Poisson deviates. lambda ≤ 0 returns 0 (degenerate mean-0 case).
    /// </summary>
    private static long SamplePoisson(double lambda, Random random)
    {
        if (lambda <= 0.0)
        {
            return 0;
        }

        double limit = Math.Exp(-lambda);
        long count = 0;
        double product = random.NextDouble();
        while (product > limit)
        {
            count++;
            product *= random.NextDouble();
        }

        return count;
    }

    /// <summary>
    /// Samples a Binomial(n, p) variate by summing n independent Bernoulli(p) trials. Exact for the
    /// per-channel sample sizes encountered in mutational catalogs. Source: standard definition of the
    /// Binomial distribution as the number of successes in n independent Bernoulli(p) trials.
    /// </summary>
    private static long SampleBinomial(long n, double p, Random random)
    {
        if (p <= 0.0)
        {
            return 0;
        }

        if (p >= 1.0)
        {
            return n;
        }

        long successes = 0;
        for (long i = 0; i < n; i++)
        {
            if (random.NextDouble() < p)
            {
                successes++;
            }
        }

        return successes;
    }

    /// <summary>Arithmetic mean of a non-empty array.</summary>
    private static double Mean(double[] values)
    {
        double sum = 0.0;
        foreach (double value in values)
        {
            sum += value;
        }

        return sum / values.Length;
    }

    /// <summary>
    /// Computes the empirical percentile (quantile) of <paramref name="values"/> at probability
    /// <paramref name="probability"/> ∈ [0, 1] using linear interpolation between order statistics on the
    /// 0-based rank h = probability·(n − 1) (the "linear interpolation of the modes for the order statistics"
    /// /R type-7 / NumPy default convention): result = x₍⌊h⌋₎ + (h − ⌊h⌋)·(x₍⌊h⌋₊₁₎ − x₍⌊h⌋₎) over the sorted
    /// values. Source: Hyndman R.J. &amp; Fan Y. (1996), <i>The American Statistician</i> 50(4):361–365
    /// (sample-quantile type 7); used to realise the Efron (1979) percentile bootstrap interval.
    /// </summary>
    private static double Percentile(double[] values, double probability)
    {
        int n = values.Length;
        if (n == 1)
        {
            return values[0];
        }

        var sorted = (double[])values.Clone();
        Array.Sort(sorted);

        double rank = probability * (n - 1);
        int lowerIndex = (int)Math.Floor(rank);
        int upperIndex = (int)Math.Ceiling(rank);
        double fraction = rank - lowerIndex;

        if (lowerIndex == upperIndex)
        {
            return sorted[lowerIndex];
        }

        return sorted[lowerIndex] + fraction * (sorted[upperIndex] - sorted[lowerIndex]);
    }

    #endregion


    #region Mutational Process Classification (ONCO-SIG-004)

    /// <summary>
    /// Minimum normalized relative contribution for a single mutational signature to be reported as
    /// present/active. A signature whose contribution falls below this fraction is excluded (set to zero).
    /// Source: Rosenthal R. et al. (2016), deconstructSigs, <i>Genome Biology</i> 17:31 — "the weights W
    /// are normalized between 0 and 1 and any signature with Wᵢ &lt; 6% is excluded"
    /// (https://doi.org/10.1186/s13059-016-0893-4); reference implementation <c>whichSignatures.R</c>
    /// declares <c>signature.cutoff = 0.06</c> and applies <c>weights[weights &lt; signature.cutoff] &lt;- 0</c>.
    /// The comparison is strict less-than, so a contribution of exactly 0.06 is retained. Value = 0.06.
    /// </summary>
    public const double DefaultSignatureContributionCutoff = 0.06;

    /// <summary>
    /// A recognized mutational process (mutagenic aetiology) inferred from active COSMIC SBS signatures.
    /// Aetiology assignments are from the COSMIC SBS catalogue (https://cancer.sanger.ac.uk/signatures/sbs/;
    /// Alexandrov et al. 2020, <i>Nature</i> 578:94–101).
    /// </summary>
    public enum MutationalProcess
    {
        /// <summary>Signature label not mapped to any recognized COSMIC aetiology.</summary>
        Unknown = 0,

        /// <summary>
        /// Aging / clock-like mutagenesis. COSMIC SBS1 ("Spontaneous deamination of 5-methylcytosine
        /// (clock-like signature)") and SBS5 ("Unknown (clock-like signature)").
        /// </summary>
        Aging,

        /// <summary>
        /// APOBEC cytidine-deaminase activity. COSMIC SBS2 and SBS13 ("Activity of APOBEC family of
        /// cytidine deaminases").
        /// </summary>
        Apobec,

        /// <summary>Tobacco smoking. COSMIC SBS4 ("Tobacco smoking").</summary>
        TobaccoSmoking,

        /// <summary>Ultraviolet light exposure. COSMIC SBS7a/7b/7c/7d ("Ultraviolet light exposure").</summary>
        UltravioletLight,

        /// <summary>
        /// Defective DNA mismatch repair. COSMIC SBS6, SBS15, SBS26 ("Defective DNA mismatch repair") and
        /// SBS20 ("Concurrent POLD1 mutations and defective DNA mismatch repair").
        /// </summary>
        MismatchRepairDeficiency,
    }

    /// <summary>
    /// COSMIC SBS signature label → mutational process map, taken verbatim from the COSMIC SBS catalogue
    /// proposed-aetiology strings (https://cancer.sanger.ac.uk/signatures/sbs/; Alexandrov et al. 2020).
    /// Only the five canonical processes named for ONCO-SIG-004 are mapped: Aging (SBS1, SBS5),
    /// APOBEC (SBS2, SBS13), Tobacco smoking (SBS4), UV (SBS7a–d), MMR deficiency (SBS6, SBS15, SBS20, SBS26).
    /// Labels are matched case-insensitively. Signatures outside this map resolve to
    /// <see cref="MutationalProcess.Unknown"/>.
    /// </summary>
    private static readonly IReadOnlyDictionary<string, MutationalProcess> SbsToProcess =
        new Dictionary<string, MutationalProcess>(StringComparer.OrdinalIgnoreCase)
        {
            ["SBS1"] = MutationalProcess.Aging,
            ["SBS5"] = MutationalProcess.Aging,
            ["SBS2"] = MutationalProcess.Apobec,
            ["SBS13"] = MutationalProcess.Apobec,
            ["SBS4"] = MutationalProcess.TobaccoSmoking,
            ["SBS7a"] = MutationalProcess.UltravioletLight,
            ["SBS7b"] = MutationalProcess.UltravioletLight,
            ["SBS7c"] = MutationalProcess.UltravioletLight,
            ["SBS7d"] = MutationalProcess.UltravioletLight,
            ["SBS6"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS15"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS20"] = MutationalProcess.MismatchRepairDeficiency,
            ["SBS26"] = MutationalProcess.MismatchRepairDeficiency,
        };

    /// <summary>
    /// Resolves a COSMIC SBS signature label (e.g. <c>"SBS2"</c>, <c>"SBS7a"</c>) to its mutational process
    /// per the COSMIC proposed-aetiology assignments. Matching is case-insensitive. Unmapped or
    /// unknown-aetiology labels return <see cref="MutationalProcess.Unknown"/>. Source: COSMIC SBS catalogue
    /// (https://cancer.sanger.ac.uk/signatures/sbs/; Alexandrov et al. 2020, <i>Nature</i> 578:94–101).
    /// </summary>
    /// <param name="signatureLabel">A COSMIC SBS signature label.</param>
    /// <returns>The mapped <see cref="MutationalProcess"/>, or <see cref="MutationalProcess.Unknown"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="signatureLabel"/> is null.</exception>
    public static MutationalProcess GetMutationalProcess(string signatureLabel)
    {
        ArgumentNullException.ThrowIfNull(signatureLabel);
        return SbsToProcess.TryGetValue(signatureLabel, out MutationalProcess process)
            ? process
            : MutationalProcess.Unknown;
    }

    /// <summary>
    /// The aggregated activity of one mutational process within a tumour sample, produced by
    /// <see cref="ClassifyMutationalProcess"/>.
    /// </summary>
    /// <param name="Process">The mutational process (aetiology).</param>
    /// <param name="Contribution">
    /// The summed normalized relative contribution of all of this process's signatures that survived the
    /// per-signature cutoff (a fraction in [0, 1]). Source: additive deconstructSigs weights (Rosenthal 2016).
    /// </param>
    public readonly record struct ProcessActivity(MutationalProcess Process, double Contribution);

    /// <summary>
    /// The result of classifying which mutational processes are active in a tumour sample.
    /// </summary>
    /// <param name="ActiveProcesses">
    /// The active processes (those with at least one surviving signature), in descending contribution order
    /// then by process for deterministic ordering. Empty when no signature survives the cutoff.
    /// </param>
    /// <param name="DominantProcess">
    /// The active process with the largest aggregated contribution, or <see cref="MutationalProcess.Unknown"/>
    /// when no process is active.
    /// </param>
    public readonly record struct MutationalProcessClassification(
        IReadOnlyList<ProcessActivity> ActiveProcesses,
        MutationalProcess DominantProcess);

    /// <summary>
    /// Classifies the mutational processes active in a tumour from per-signature exposures (activities).
    /// <para>
    /// The exposures are first converted to normalized relative contributions (each divided by their sum,
    /// so they sum to 1) per deconstructSigs ("the weights W are normalized between 0 and 1"). A signature
    /// is then declared present/active only when its normalized contribution is at least
    /// <paramref name="contributionCutoff"/> (default 6%); contributions below the cutoff are dropped, so
    /// the surviving contributions can sum to less than 1 (Rosenthal et al. 2016, deconstructSigs:
    /// "any signature with Wᵢ &lt; 6% is excluded"; <c>weights[weights &lt; signature.cutoff] &lt;- 0</c>).
    /// Each surviving signature is mapped to its mutational process via the COSMIC SBS aetiology map
    /// (<see cref="GetMutationalProcess"/>), and the per-signature contributions are summed per process.
    /// The dominant process is the active process with the largest aggregated contribution.
    /// </para>
    /// <para>
    /// Reference signature <i>profiles</i> are caller-supplied elsewhere (e.g. COSMIC SBS matrices used by
    /// <see cref="FitSignatures(IReadOnlyList{double}, IReadOnlyList{IReadOnlyList{double}})"/>); this method
    /// consumes only the resulting exposures and their COSMIC labels.
    /// </para>
    /// Sources: Rosenthal R. et al. (2016), <i>Genome Biology</i> 17:31
    /// (https://doi.org/10.1186/s13059-016-0893-4); COSMIC SBS catalogue
    /// (https://cancer.sanger.ac.uk/signatures/sbs/); Alexandrov L.B. et al. (2020), <i>Nature</i> 578:94–101.
    /// </summary>
    /// <param name="exposures">
    /// Per-signature exposures as (COSMIC SBS label, non-negative activity) pairs. Activities are arbitrary
    /// non-negative magnitudes (mutation counts or proportions); they are normalized internally.
    /// </param>
    /// <param name="contributionCutoff">
    /// Minimum normalized relative contribution (strict lower bound is a contribution that is &lt; this value)
    /// for a signature to be active; default <see cref="DefaultSignatureContributionCutoff"/> (0.06).
    /// </param>
    /// <returns>The active processes (descending contribution) and the dominant process.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="exposures"/> (or a label) is null.</exception>
    /// <exception cref="ArgumentException">An exposure is negative.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="contributionCutoff"/> is outside [0, 1).</exception>
    public static MutationalProcessClassification ClassifyMutationalProcess(
        IReadOnlyList<(string Signature, double Exposure)> exposures,
        double contributionCutoff = DefaultSignatureContributionCutoff)
    {
        ArgumentNullException.ThrowIfNull(exposures);

        if (double.IsNaN(contributionCutoff) || contributionCutoff < 0.0 || contributionCutoff >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(contributionCutoff), contributionCutoff,
                "The contribution cutoff must be in the half-open interval [0, 1).");
        }

        double total = 0.0;
        for (int i = 0; i < exposures.Count; i++)
        {
            (string signature, double exposure) = exposures[i];
            ArgumentNullException.ThrowIfNull(signature);
            if (double.IsNaN(exposure) || exposure < 0.0)
            {
                throw new ArgumentException("Exposures must be non-negative.", nameof(exposures));
            }

            total += exposure;
        }

        var empty = new MutationalProcessClassification(
            Array.Empty<ProcessActivity>(), MutationalProcess.Unknown);

        // Σ exposure = 0 ⇒ normalization undefined ⇒ no active processes (INV-5).
        if (total <= 0.0)
        {
            return empty;
        }

        // Per-signature normalized contribution, cutoff, then aggregate surviving contributions by process.
        var byProcess = new Dictionary<MutationalProcess, double>();
        foreach ((string signature, double exposure) in exposures)
        {
            double contribution = exposure / total;
            if (contribution < contributionCutoff)
            {
                continue; // deconstructSigs: weights[weights < signature.cutoff] <- 0
            }

            MutationalProcess process = GetMutationalProcess(signature);
            if (process == MutationalProcess.Unknown)
            {
                continue; // unmapped-aetiology signatures contribute to no recognized process (COSMIC)
            }

            byProcess.TryGetValue(process, out double accumulated);
            byProcess[process] = accumulated + contribution;
        }

        if (byProcess.Count == 0)
        {
            return empty;
        }

        // Deterministic order: descending contribution, then by process enum for ties.
        var active = byProcess
            .Select(kv => new ProcessActivity(kv.Key, kv.Value))
            .OrderByDescending(p => p.Contribution)
            .ThenBy(p => p.Process)
            .ToArray();

        return new MutationalProcessClassification(active, active[0].Process);
    }

    #endregion


    #region Tumor Gene Expression Outlier / Signature Score (ONCO-EXPR-001)

    /// <summary>
    /// Direction of an expression outlier relative to the reference cohort.
    /// </summary>
    public enum ExpressionDirection
    {
        /// <summary>z &gt; +threshold — expression is elevated versus the reference cohort (overexpressed).</summary>
        Over,

        /// <summary>z &lt; −threshold — expression is reduced versus the reference cohort (underexpressed).</summary>
        Under,
    }

    /// <summary>
    /// A single gene whose sample expression is an outlier relative to its reference cohort.
    /// </summary>
    /// <param name="Gene">Gene identifier.</param>
    /// <param name="ZScore">The expression z-score z = (value − μ)/σ of the gene in this sample.</param>
    /// <param name="Direction">Whether the gene is over- or under-expressed (sign of the z-score).</param>
    public readonly record struct ExpressionOutlier(string Gene, double ZScore, ExpressionDirection Direction);

    /// <summary>
    /// Computes the expression z-score of a single sample value relative to a reference cohort:
    /// z = (value − μ) / σ, where μ is the arithmetic mean and σ is the <b>sample</b> standard deviation
    /// (divisor n − 1) of the cohort.
    /// </summary>
    /// <remarks>
    /// Source: cBioPortal mRNA z-score normalization specification
    /// (https://docs.cbioportal.org/z-score-normalization-script/) — z = "(r - mu)/sigma where r is the raw
    /// expression value, and mu and sigma are the mean and standard deviation". The reference implementation
    /// <c>NormalizeExpressionLevels.java</c> (cbioportal-core) computes σ with divisor (n − 1)
    /// (<c>std = std/(double)(v.length-1)</c>), i.e. the sample standard deviation, and aborts when σ = 0.
    /// </remarks>
    /// <param name="value">The sample expression value (raw <c>r</c>) on a normalization scale.</param>
    /// <param name="referenceCohort">Reference expression values for this gene; at least two values required.</param>
    /// <returns>The z-score of <paramref name="value"/> relative to the cohort.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="referenceCohort"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// The cohort has fewer than two values (sample standard deviation is undefined), or the cohort has a
    /// standard deviation of 0 (no defined z-score; mirrors the reference implementation's fatal error).
    /// </exception>
    public static double CalculateExpressionZScore(double value, IReadOnlyList<double> referenceCohort)
    {
        ArgumentNullException.ThrowIfNull(referenceCohort);

        int n = referenceCohort.Count;

        // Sample standard deviation (divisor n − 1) is undefined for fewer than two observations.
        if (n < 2)
        {
            throw new ArgumentException(
                "Reference cohort must contain at least two values to compute a sample standard deviation.",
                nameof(referenceCohort));
        }

        // Zero-spread (all values identical) is the documented σ = 0 error case
        // (NormalizeExpressionLevels.java fatal error). Detect it robustly by the
        // cohort's RANGE rather than from the computed SD: with a constant cohort,
        // accumulating the mean incurs a ≤1-ULP rounding error, which leaves a tiny
        // non-zero Σ(rᵢ − μ)² and an SD on the order of 1e-17 — the `sd == 0` test
        // alone would then be bypassed and emit a spurious ~1e20 z-score. The
        // max == min check is exact for a no-spread cohort and immune to that drift.
        double first = referenceCohort[0];
        bool anyDifferent = false;
        for (int i = 1; i < n; i++)
        {
            if (referenceCohort[i] != first)
            {
                anyDifferent = true;
                break;
            }
        }

        if (!anyDifferent)
        {
            throw new ArgumentException(
                "Cannot compute a z-score relative to a reference cohort with a standard deviation of 0.",
                nameof(referenceCohort));
        }

        double mean = 0.0;
        for (int i = 0; i < n; i++)
        {
            mean += referenceCohort[i];
        }

        mean /= n;

        double sumSquaredDeviations = 0.0;
        for (int i = 0; i < n; i++)
        {
            double d = referenceCohort[i] - mean;
            sumSquaredDeviations += d * d;
        }

        // Sample standard deviation: divisor (n − 1) per NormalizeExpressionLevels.java std().
        double sd = Math.Sqrt(sumSquaredDeviations / (n - 1));

        if (sd == 0.0)
        {
            throw new ArgumentException(
                "Cannot compute a z-score relative to a reference cohort with a standard deviation of 0.",
                nameof(referenceCohort));
        }

        return (value - mean) / sd;
    }

    /// <summary>
    /// Identifies genes whose sample expression is an outlier relative to per-gene reference cohorts, using
    /// the z-score rule z &gt; +threshold (overexpressed) or z &lt; −threshold (underexpressed).
    /// </summary>
    /// <remarks>
    /// Source: cBioPortal FAQ (https://docs.cbioportal.org/user-guide/faq/) — "samples with expression
    /// z-scores &gt;2 or &lt;-2 in any queried genes are considered altered." The default threshold is 2.0
    /// (<see cref="DefaultExpressionOutlierThreshold"/>); the comparison is strict, so |z| = threshold is not
    /// reported as an outlier.
    /// </remarks>
    /// <param name="sampleExpression">The sample's expression value per gene.</param>
    /// <param name="referenceCohorts">Per-gene reference cohort values; must contain every gene in the sample.</param>
    /// <param name="threshold">Absolute z-score threshold (must be positive). Default 2.0.</param>
    /// <returns>The outlier genes, in the iteration order of <paramref name="sampleExpression"/>.</returns>
    /// <exception cref="ArgumentNullException">Either dictionary is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> is not positive.</exception>
    /// <exception cref="ArgumentException">A sampled gene has no reference cohort, or a cohort is degenerate.</exception>
    public static IReadOnlyList<ExpressionOutlier> IdentifyOutlierGenes(
        IReadOnlyDictionary<string, double> sampleExpression,
        IReadOnlyDictionary<string, IReadOnlyList<double>> referenceCohorts,
        double threshold = DefaultExpressionOutlierThreshold)
    {
        ArgumentNullException.ThrowIfNull(sampleExpression);
        ArgumentNullException.ThrowIfNull(referenceCohorts);

        if (threshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(threshold), threshold, "Outlier threshold must be positive.");
        }

        var outliers = new List<ExpressionOutlier>();

        foreach (KeyValuePair<string, double> entry in sampleExpression)
        {
            if (!referenceCohorts.TryGetValue(entry.Key, out IReadOnlyList<double>? cohort))
            {
                throw new ArgumentException(
                    $"No reference cohort supplied for gene '{entry.Key}'.", nameof(referenceCohorts));
            }

            double z = CalculateExpressionZScore(entry.Value, cohort);

            // Strict thresholds: z > +t ⇒ over, z < −t ⇒ under; |z| = t is not an outlier (cBioPortal FAQ).
            if (z > threshold)
            {
                outliers.Add(new ExpressionOutlier(entry.Key, z, ExpressionDirection.Over));
            }
            else if (z < -threshold)
            {
                outliers.Add(new ExpressionOutlier(entry.Key, z, ExpressionDirection.Under));
            }
        }

        return outliers;
    }

    /// <summary>
    /// Computes the combined z-score (gene-signature / pathway activity score) over a set of member-gene
    /// z-scores: a = (Σᵢ zᵢ) / √k, where k is the number of member genes.
    /// </summary>
    /// <remarks>
    /// Source: Lee E. et al. (2008), "Inferring Pathway Activity toward Precise Disease Classification",
    /// PLoS Comput Biol 4(11):e1000217 (https://doi.org/10.1371/journal.pcbi.1000217) — the per-gene
    /// z-scores of a gene set are "averaged into a combined z-score … the square root of the number of member
    /// genes is used in the denominator to stabilize the variance of the mean." Corroborated by the GSVA
    /// "combined z-score" method. Member z-scores are caller-supplied (e.g. from
    /// <see cref="CalculateExpressionZScore"/>); the signature gene set is caller-defined.
    /// </remarks>
    /// <param name="memberZScores">The per-gene z-scores of the signature's member genes (k ≥ 1).</param>
    /// <returns>The combined z-score activity a = (Σ z) / √k.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="memberZScores"/> is null.</exception>
    /// <exception cref="ArgumentException">The signature is empty (k = 0; the score is undefined).</exception>
    public static double CalculateSignatureScore(IReadOnlyList<double> memberZScores)
    {
        ArgumentNullException.ThrowIfNull(memberZScores);

        int k = memberZScores.Count;
        if (k == 0)
        {
            throw new ArgumentException(
                "Signature must contain at least one member gene z-score.", nameof(memberZScores));
        }

        double sum = 0.0;
        for (int i = 0; i < k; i++)
        {
            sum += memberZScores[i];
        }

        return sum / Math.Sqrt(k);
    }

    #endregion
}
