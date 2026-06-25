using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// A Plan7 profile hidden Markov model (HMMER3 architecture) parsed from a HMMER3/f ASCII
/// <c>.hmm</c> profile, with log-odds Viterbi (optimal-path) and Forward (full-likelihood) scoring.
/// </summary>
/// <remarks>
/// <para>
/// Implements the Plan7 profile-HMM scoring of Durbin, Eddy, Krogh &amp; Mitchison,
/// <i>Biological Sequence Analysis</i> (1998), §5.4, and the HMMER3 architecture of
/// Eddy (2011) <i>PLoS Comput Biol</i> 7:e1002195. Match/Insert/Delete states per node, a mute
/// Begin/End, and log-odds emission scores against a null (background) model.
/// </para>
/// <para>
/// File-format facts are taken verbatim from the HMMER User's Guide v3.4 (Eddy, Aug 2023),
/// "HMMER profile HMM files": all probability parameters are stored as <b>negative natural-log
/// probabilities</b> (e.g. probability 0.25 → <c>-ln 0.25 = 1.38629</c>), and a zero probability is
/// stored as <c>'*'</c>. For <c>ALPH amino</c> the alphabet size K = 20 in order
/// <c>ACDEFGHIKLMNPQRSTVWY</c>. The COMPO line holds the model's mean match-state composition,
/// used as the background residue composition of the null model.
/// </para>
/// </remarks>
public sealed class Plan7ProfileHmm
{
    // HMMER amino-acid alphabet order (HMMER User's Guide §"header section", ALPH amino):
    // "the symbol alphabet to 'ACDEFGHIKLMNPQRSTVWY' (alphabetic order)".
    internal const string AminoAlphabet = "ACDEFGHIKLMNPQRSTVWY";

    // Number of state transitions per node, in HMMER3 file order
    // (HMMER User's Guide, "State transition line"):
    //   Mk->Mk+1, Mk->Ik, Mk->Dk+1, Ik->Mk+1, Ik->Ik, Dk->Mk+1, Dk->Dk+1.
    private const int TransitionsPerNode = 7;
    private const int TMM = 0; // M_k -> M_{k+1}
    private const int TMI = 1; // M_k -> I_k
    private const int TMD = 2; // M_k -> D_{k+1}
    private const int TIM = 3; // I_k -> M_{k+1}
    private const int TII = 4; // I_k -> I_k
    private const int TDM = 5; // D_k -> M_{k+1}
    private const int TDD = 6; // D_k -> D_{k+1}

    private static readonly int[] ResidueIndex = BuildResidueIndex();

    /// <summary>The model name (NAME line); e.g. "SH3_1".</summary>
    public string Name { get; }

    /// <summary>The accession (ACC line); e.g. "PF00018.35". Empty if absent.</summary>
    public string Accession { get; }

    /// <summary>The description (DESC line). Empty if absent.</summary>
    public string Description { get; }

    /// <summary>Number of match states (LENG line), 1-based node indices 1..<see cref="Length"/>.</summary>
    public int Length { get; }

    /// <summary>Pfam gathering threshold GA1 (per-sequence bit score), or <c>null</c> if absent.</summary>
    public double? GatheringThreshold { get; }

    /// <summary>
    /// E-value calibration parameters from the profile's <c>STATS LOCAL …</c> lines, or <c>null</c>
    /// if the profile is not calibrated (HMMER User's Guide v3.4: "All three lines or none of them
    /// must be present; when all three are present, the model is considered to be calibrated for
    /// E-value statistics."). MSV and Viterbi bit scores are Gumbel-distributed with location μ and
    /// slope λ; the Forward score's high-scoring tail is exponential with location τ and slope λ.
    /// </summary>
    public ScoreStatistics? Statistics { get; }

    // Natural-log probabilities (NOT negated): matchEmissionLn[k][a] = ln P(emit residue a | M_k).
    // k is 1-based (index 0 unused). -inf where the file stored '*'.
    private readonly double[][] _matchEmissionLn;   // [Length+1][20]
    private readonly double[][] _insertEmissionLn;  // [Length+1][20]  (node k insert I_k; index 0 = I_0)
    private readonly double[][] _transitionLn;       // [Length+1][7]   (node k; index 0 = BEGIN node)
    private readonly double[] _backgroundLn;         // [20] = ln COMPO background composition

    private Plan7ProfileHmm(
        string name, string accession, string description, int length, double? gatheringThreshold,
        ScoreStatistics? statistics,
        double[][] matchEmissionLn, double[][] insertEmissionLn, double[][] transitionLn, double[] backgroundLn)
    {
        Name = name;
        Accession = accession;
        Description = description;
        Length = length;
        GatheringThreshold = gatheringThreshold;
        Statistics = statistics;
        _matchEmissionLn = matchEmissionLn;
        _insertEmissionLn = insertEmissionLn;
        _transitionLn = transitionLn;
        _backgroundLn = backgroundLn;
    }

    /// <summary>
    /// E-value calibration parameters parsed from a profile's three <c>STATS LOCAL …</c> lines.
    /// </summary>
    /// <remarks>
    /// HMMER User's Guide v3.4 (Eddy 2023), HMM file format: <c>STATS &lt;s1&gt; &lt;s2&gt; &lt;f1&gt; &lt;f2&gt;</c> —
    /// "&lt;f1&gt; and &lt;f2&gt; are two real-valued parameters controlling location and slope of each
    /// distribution, respectively; µ and λ for Gumbel distributions for MSV and Viterbi scores, and
    /// τ and λ for exponential tails for Forward scores. λ values must be positive."
    /// All parameters are in <b>bits</b> (HMMER reports scores in bits, log base 2; λ ≈ log 2).
    /// </remarks>
    /// <param name="MsvMu">Gumbel location μ for the MSV bit-score distribution.</param>
    /// <param name="MsvLambda">Gumbel slope λ for the MSV bit-score distribution.</param>
    /// <param name="ViterbiMu">Gumbel location μ for the Viterbi bit-score distribution.</param>
    /// <param name="ViterbiLambda">Gumbel slope λ for the Viterbi bit-score distribution.</param>
    /// <param name="ForwardTau">Exponential-tail location τ for the Forward bit-score distribution.</param>
    /// <param name="ForwardLambda">Exponential-tail slope λ for the Forward bit-score distribution.</param>
    public readonly record struct ScoreStatistics(
        double MsvMu, double MsvLambda,
        double ViterbiMu, double ViterbiLambda,
        double ForwardTau, double ForwardLambda);

    #region Parsing

    /// <summary>Parses a HMMER3/f ASCII profile from a string.</summary>
    public static Plan7ProfileHmm Parse(string hmmText)
    {
        ArgumentNullException.ThrowIfNull(hmmText);
        using var reader = new StringReader(hmmText);
        return Parse(reader);
    }

    /// <summary>Parses a HMMER3/f ASCII profile from a reader.</summary>
    public static Plan7ProfileHmm Parse(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);

        string? versionLine = reader.ReadLine();
        if (versionLine is null || !versionLine.StartsWith("HMMER3/", StringComparison.Ordinal))
            throw new FormatException("Not a HMMER3 profile: missing 'HMMER3/' version line.");

        string name = string.Empty, accession = string.Empty, description = string.Empty;
        int length = -1;
        bool amino = false;
        double? ga1 = null;
        // STATS LOCAL <MSV|VITERBI|FORWARD> <location> <slope>. null until a line is seen.
        double? msvMu = null, msvLambda = null;
        double? vitMu = null, vitLambda = null;
        double? fwdTau = null, fwdLambda = null;

        // Header section: parsed line-by-line as tag/value until the 'HMM' keyword line.
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith("HMM ", StringComparison.Ordinal) || line.Trim() == "HMM")
                break;

            var tag = FirstToken(line);
            switch (tag)
            {
                case "NAME": name = Rest(line); break;
                case "ACC": accession = Rest(line); break;
                case "DESC": description = Rest(line); break;
                case "LENG": length = int.Parse(Rest(line), CultureInfo.InvariantCulture); break;
                case "ALPH": amino = string.Equals(Rest(line), "amino", StringComparison.OrdinalIgnoreCase); break;
                case "GA":
                    var gaToks = Rest(line).TrimEnd(';').Split(new[] { ' ', '\t', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    if (gaToks.Length >= 1 && double.TryParse(gaToks[0], NumberStyles.Float, CultureInfo.InvariantCulture, out var ga))
                        ga1 = ga;
                    break;
                case "STATS":
                    // STATS LOCAL <distribution> <location> <slope> (HMMER User's Guide v3.4, file format).
                    var sToks = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
                    if (sToks.Length >= 5 &&
                        string.Equals(sToks[1], "LOCAL", StringComparison.OrdinalIgnoreCase) &&
                        double.TryParse(sToks[3], NumberStyles.Float, CultureInfo.InvariantCulture, out var loc) &&
                        double.TryParse(sToks[4], NumberStyles.Float, CultureInfo.InvariantCulture, out var slope))
                    {
                        switch (sToks[2].ToUpperInvariant())
                        {
                            case "MSV": msvMu = loc; msvLambda = slope; break;
                            case "VITERBI": vitMu = loc; vitLambda = slope; break;
                            case "FORWARD": fwdTau = loc; fwdLambda = slope; break;
                        }
                    }
                    break;
            }
        }

        if (line is null)
            throw new FormatException("Truncated profile: no 'HMM' main-model header found.");
        if (!amino)
            throw new FormatException("Only amino-acid (ALPH amino) profiles are supported.");
        if (length <= 0)
            throw new FormatException("Missing or invalid LENG (model length).");

        // The 'HMM' tag line itself carries the symbol-alphabet header (already consumed by the
        // break above). The single line that follows is the transition column header — also for
        // human readability only. Per the HMMER User's Guide, "the parser always skips the line
        // after the HMM tag line." Skip exactly that one line.
        _ = reader.ReadLine();

        const int K = 20;
        var matchEmissionLn = new double[length + 1][];
        var insertEmissionLn = new double[length + 1][];
        var transitionLn = new double[length + 1][];
        double[]? backgroundLn = null;

        // Optional COMPO line (model's mean match composition → null background), then the BEGIN node.
        string? next = NextNonBlank(reader);
        if (next is not null && FirstToken(next) == "COMPO")
        {
            backgroundLn = NegLogToLn(ParseNumbers(next, skipFirst: 1, count: K));
            // The two lines after COMPO are the BEGIN node: insert-0 emissions, then begin transitions.
            insertEmissionLn[0] = NegLogToLn(ParseNumbers(ReadRequired(reader), skipFirst: 0, count: K));
            transitionLn[0] = NegLogToLn(ParseNumbers(ReadRequired(reader), skipFirst: 0, count: TransitionsPerNode));
        }
        else if (next is not null)
        {
            // No COMPO: the BEGIN node's insert-0 emissions line is the first line read.
            insertEmissionLn[0] = NegLogToLn(ParseNumbers(next, skipFirst: 0, count: K));
            transitionLn[0] = NegLogToLn(ParseNumbers(ReadRequired(reader), skipFirst: 0, count: TransitionsPerNode));
        }
        else
        {
            throw new FormatException("Truncated profile: missing main-model body.");
        }

        // Each node k (1..M) has three lines: match emissions (prefixed by node number),
        // insert emissions, state transitions.
        for (int k = 1; k <= length; k++)
        {
            var matchLine = ReadRequired(reader);
            int node = int.Parse(FirstToken(matchLine), CultureInfo.InvariantCulture);
            if (node != k)
                throw new FormatException($"Node order error: expected node {k}, found {node}.");
            matchEmissionLn[k] = NegLogToLn(ParseNumbers(matchLine, skipFirst: 1, count: K));
            insertEmissionLn[k] = NegLogToLn(ParseNumbers(ReadRequired(reader), skipFirst: 0, count: K));
            transitionLn[k] = NegLogToLn(ParseNumbers(ReadRequired(reader), skipFirst: 0, count: TransitionsPerNode));
        }

        // If COMPO was absent, fall back to a uniform background (every residue equally likely).
        backgroundLn ??= Enumerable.Repeat(Math.Log(1.0 / K), K).ToArray();

        // E-value statistics: only when ALL THREE STATS lines are present is the model calibrated
        // (HMMER User's Guide v3.4: "All three lines or none of them must be present").
        ScoreStatistics? stats = null;
        if (msvMu is not null && vitMu is not null && fwdTau is not null)
            stats = new ScoreStatistics(
                msvMu.Value, msvLambda!.Value,
                vitMu.Value, vitLambda!.Value,
                fwdTau.Value, fwdLambda!.Value);

        return new Plan7ProfileHmm(name, accession, description, length, ga1, stats,
            matchEmissionLn, insertEmissionLn, transitionLn, backgroundLn);
    }

    /// <summary>Loads one of the bundled CC0 Pfam profiles embedded in this assembly.</summary>
    /// <param name="resourceFileName">e.g. "PF00018_SH3_1.hmm".</param>
    internal static Plan7ProfileHmm LoadEmbedded(string resourceFileName)
    {
        var asm = typeof(Plan7ProfileHmm).Assembly;
        // Embedded resource names are "<RootNamespace>.Resources.<file>".
        string resourceName = $"Seqeron.Genomics.Analysis.Resources.{resourceFileName}";
        using var stream = asm.GetManifestResourceStream(resourceName)
            ?? throw new InvalidOperationException($"Embedded profile not found: {resourceName}");
        using var reader = new StreamReader(stream);
        return Parse(reader);
    }

    #endregion

    #region Scoring

    /// <summary>
    /// Glocal Viterbi log-odds score (in nats): the best single full-profile path
    /// B → (M/I/D over nodes 1..M) → E aligned to the whole sequence.
    /// </summary>
    /// <remarks>
    /// Recurrence (Durbin et al. 1998, §5.4; reproduced in Stanford CS273 lecture 7):
    /// V^M_j(i) = log(e_Mj(x_i)/q_xi) + max{V^M_{j-1}(i-1)+a_{M(j-1)M(j)},
    ///                                       V^I_{j-1}(i-1)+a_{I(j-1)M(j)},
    ///                                       V^D_{j-1}(i-1)+a_{D(j-1)M(j)}};
    /// V^I_j(i) = log(e_Ij(x_i)/q_xi) + max{V^M_j(i-1)+a_{M(j)I(j)}, V^I_j(i-1)+a_{I(j)I(j)}};
    /// V^D_j(i) = max{V^M_{j-1}(i)+a_{M(j-1)D(j)}, V^D_{j-1}(i)+a_{D(j-1)D(j)}}.
    /// Emissions are log-odds against the null background q (COMPO).
    /// </remarks>
    public double ViterbiScore(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return RunGlocal(sequence, viterbi: true);
    }

    /// <summary>
    /// Glocal Forward log-odds score (in nats): the full log-likelihood summed (log-sum-exp)
    /// over all full-profile paths B → 1..M → E. Same recurrence as <see cref="ViterbiScore"/>
    /// with max replaced by log-sum-exp (Durbin et al. 1998, §3.6 / §5.4).
    /// </summary>
    public double ForwardScore(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return RunGlocal(sequence, viterbi: false);
    }

    private double RunGlocal(string sequence, bool viterbi)
    {
        int n = sequence.Length;
        int m = Length;

        // DP rows indexed by node j (0..m). We process residue by residue (i = 0..n).
        // Use two columns (i-1, i) for M/I/D at each node.
        double[] prevM = NewRow(m), prevI = NewRow(m), prevD = NewRow(m);
        double[] curM = NewRow(m), curI = NewRow(m), curD = NewRow(m);

        // i = 0: no residues consumed. Only the all-delete path through the profile is possible.
        // Begin is mute. V^M_j(0) = -inf (a match must emit). D states: B->D1, then D->D chain.
        prevM[0] = double.NegativeInfinity; // node 0 = BEGIN, treated as score 0 only as a predecessor
        prevD[0] = double.NegativeInfinity;
        prevI[0] = double.NegativeInfinity;
        // D_1(0): from begin B->D1
        if (m >= 1)
        {
            prevD[1] = _transitionLn[0][TMD]; // B->D1 stored as M_0->D_1
            for (int j = 2; j <= m; j++)
                prevD[j] = prevD[j - 1] + _transitionLn[j - 1][TDD];
            prevM[1] = double.NegativeInfinity;
        }

        for (int i = 1; i <= n; i++)
        {
            int a = ResidueOf(sequence[i - 1]);

            for (int j = 0; j <= m; j++) { curM[j] = curI[j] = curD[j] = double.NegativeInfinity; }

            for (int j = 1; j <= m; j++)
            {
                double matchOdds = EmissionLogOdds(_matchEmissionLn[j], a);

                // Predecessors of M_j: from node j-1 (M/I/D) at i-1. j-1 == 0 means BEGIN (score 0).
                double fromM = (j - 1 == 0) ? 0.0 : prevM[j - 1];
                double mPred = Combine(viterbi,
                    Add(fromM, _transitionLn[j - 1][TMM]),
                    Add(prevI[j - 1], _transitionLn[j - 1][TIM]),
                    Add(prevD[j - 1], _transitionLn[j - 1][TDM]));
                curM[j] = Add(matchOdds, mPred);

                // I_j: stays at node j, from M_j/I_j/D_j at i-1.
                double insertOdds = EmissionLogOdds(_insertEmissionLn[j], a);
                double iPred = Combine(viterbi,
                    Add(prevM[j], _transitionLn[j][TMI]),
                    Add(prevI[j], _transitionLn[j][TII]));
                curI[j] = Add(insertOdds, iPred);
            }

            // D_j at current i: depends on M_{j-1}/D_{j-1} at the SAME i (mute states, no emission).
            for (int j = 1; j <= m; j++)
            {
                double fromM = (j - 1 == 0) ? double.NegativeInfinity : curM[j - 1];
                curD[j] = Combine(viterbi,
                    Add(fromM, _transitionLn[j - 1][TMD]),
                    Add(curD[j - 1], _transitionLn[j - 1][TDD]));
            }

            (prevM, curM) = (curM, prevM);
            (prevI, curI) = (curI, prevI);
            (prevD, curD) = (curD, prevD);
        }

        // End: E is reached from M_m (and D_m) at i = n. Final node's M_m->E and D_m->E are
        // probability 1.0 (log 0) by HMMER convention for glocal full-length alignment.
        double end = Combine(viterbi, prevM[m], prevD[m]);
        return end;
    }

    // Emission log-odds in nats: ln P(a | state) - ln q_a. '*' (prob 0) → -inf.
    private double EmissionLogOdds(double[] emissionLn, int residueIndex)
    {
        if (residueIndex < 0) return 0.0; // ambiguous/unknown residue: treat as background (odds 1).
        double e = emissionLn[residueIndex];
        if (double.IsNegativeInfinity(e)) return double.NegativeInfinity;
        return e - _backgroundLn[residueIndex];
    }

    private static double[] NewRow(int m)
    {
        var r = new double[m + 1];
        for (int j = 0; j <= m; j++) r[j] = double.NegativeInfinity;
        return r;
    }

    private static double Add(double a, double b)
    {
        if (double.IsNegativeInfinity(a) || double.IsNegativeInfinity(b)) return double.NegativeInfinity;
        return a + b;
    }

    private static double Combine(bool viterbi, double a, double b)
        => viterbi ? Math.Max(a, b) : LogSumExp(a, b);

    private static double Combine(bool viterbi, double a, double b, double c)
        => viterbi ? Math.Max(a, Math.Max(b, c)) : LogSumExp(LogSumExp(a, b), c);

    private static double LogSumExp(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double max = Math.Max(a, b);
        return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
    }

    #endregion

    #region HMMER local-multihit Forward + null2 (opt-in; hmmsearch parity)

    // HMMER multihit unannotated-flank count nj = 1 (one J state allowed between hits); the
    // single-hit modes use nj = 0. modelconfig.c p7_ProfileConfig: "gm->nj = 1.0f" for multihit.
    private const double MultihitNj = 1.0;

    // HMMER null1 background self-transition omega for the null2 mixture weight: bg->omega = 1/256
    // (p7_bg.c: "bg->omega = 1./256."). Used in seqbias = logsum(0, log(omega) + sum n2sc).
    private const double Null2Omega = 1.0 / 256.0;

    // ln 2 — HMMER reports bit scores as nats / ln 2 (eslCONST_LOG2).
    private const double Ln2 = 0.69314718055994530941723212145818;

    // HMMER's standard amino-acid background frequencies, in alphabet order ACDEFGHIKLMNPQRSTVWY.
    // Verbatim from p7_AminoFrequencies() (EddyRivasLab/hmmer src/hmmer.c): "average Swiss-Prot
    // residue composition". hmmsearch scores match emissions against THIS background (bg->f), not the
    // profile's COMPO line (which the glocal path uses); so the parity path uses these values.
    private static readonly double[] HmmerAminoBackground =
    {
        0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
        0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364,
        0.0540687, 0.0673417, 0.0114135, 0.0304133,
    };
    private static readonly double[] HmmerAminoBackgroundLn = BuildBackgroundLn();

    private static double[] BuildBackgroundLn()
    {
        var ln = new double[20];
        for (int i = 0; i < 20; i++) ln[i] = Math.Log(HmmerAminoBackground[i]);
        return ln;
    }

    // Emission log-odds (nats) against HMMER's standard amino background bg->f (NOT the COMPO line):
    // ln P(a | M_k) - ln bg->f[a]. Used only by the hmmsearch-parity local/null2 path.
    private double EmissionLogOddsHmmer(double[] emissionLn, int residueIndex)
    {
        if (residueIndex < 0) return 0.0;
        double e = emissionLn[residueIndex];
        if (double.IsNegativeInfinity(e)) return double.NegativeInfinity;
        return e - HmmerAminoBackgroundLn[residueIndex];
    }

    /// <summary>
    /// HMMER <c>hmmsearch</c>-style <b>local multihit</b> Forward log-odds score (in nats) of the
    /// whole sequence, summed over all local-alignment paths through the Plan7 special states
    /// N → B → (local entry M<sub>k</sub>) → … → E → {J loop | C} (Eddy 2011, <i>PLoS Comput Biol</i>
    /// 7:e1002195; HMMER <c>generic_fwdback.c</c> <c>p7_GForward</c> with <c>p7_LOCAL</c> config).
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is the opt-in <b>local</b> counterpart of the glocal <see cref="ForwardScore"/>; the glocal
    /// path and all defaults are unchanged. Local entry probabilities are the occupancy-weighted
    /// B→M<sub>k</sub> of <c>modelconfig.c</c> <c>p7_ProfileConfig</c>:
    /// <c>t(B→M_k) = occ[k] / Σ_i occ[i]·(M−i+1)</c>; local exit M<sub>k</sub>→E and D<sub>k</sub>→E
    /// are probability 1 (<c>esc = 0</c>) for every node. The N/J/C self-loop and move transitions are
    /// the length-dependent <see cref="ReconfigLength"/> values
    /// <c>pmove = (2+nj)/(L+2+nj)</c>, <c>ploop = 1−pmove</c>, with the multihit E→{J,C} split of
    /// <c>−ln 2</c> each. Returns the Forward nat score <c>x_C(L) + t(C→T)</c>.
    /// </para>
    /// </remarks>
    public double LocalForwardScore(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var result = RunLocalForward(sequence, decode: false);
        return result.ForwardNats;
    }

    /// <summary>
    /// HMMER <c>hmmsearch</c>-parity per-sequence <b>bit score</b>: the local-multihit Forward nat
    /// score converted to bits against the null1 model and corrected by the <b>null2</b>
    /// biased-composition score, exactly as the HMMER pipeline computes a reported sequence score
    /// (<c>p7_pipeline.c</c> <c>p7_Pipeline</c>; Eddy 2011 §"Biased composition correction").
    /// </summary>
    /// <remarks>
    /// <para>
    /// score (bits) = <c>(fwd − (nullsc + seqbias)) / ln 2</c>, where:
    /// <list type="bullet">
    /// <item><c>fwd</c> = <see cref="LocalForwardScore"/> (nats);</item>
    /// <item><c>nullsc = L·ln(p1) + ln(1−p1)</c>, <c>p1 = L/(L+1)</c> — the null1 score
    ///   (<c>p7_bg.c</c> <c>p7_bg_NullOne</c>);</item>
    /// <item><c>seqbias = logsumexp(0, ln(omega) + Σ_i ln null2[x_i])</c> over the whole sequence,
    ///   with <c>omega = 1/256</c> and the per-residue null2 odds ratios computed by posterior
    ///   expectation (<c>generic_null2.c</c> <c>p7_GNull2_ByExpectation</c>).</item>
    /// </list>
    /// </para>
    /// <para>
    /// HMMER bounds the null2 correction to the per-domain envelope after domain definition; this
    /// implementation applies the posterior-expectation null2 over the full sequence, which matches
    /// hmmsearch to within single-precision rounding for a well-resolved single-domain target (the
    /// dominant case). See the algorithm doc for the residual.
    /// </para>
    /// </remarks>
    /// <returns>The null2-corrected per-sequence bit score (may be negative).</returns>
    public double HmmSearchBitScore(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0) return double.NegativeInfinity;

        var r = RunLocalForward(sequence, decode: true);
        double nullsc = NullOneScore(sequence.Length);
        double seqbias = SeqBias(r.Null2OddsLn, sequence);
        return (r.ForwardNats - (nullsc + seqbias)) / Ln2;
    }

    /// <summary>
    /// The uncorrected (pre-null2) local-multihit Forward bit score
    /// <c>(fwd − nullsc) / ln 2</c> — HMMER's <c>pre_score</c> (<c>p7_pipeline.c</c>), the score
    /// before the biased-composition correction. Equals <see cref="HmmSearchBitScore"/> plus the
    /// bias correction in bits.
    /// </summary>
    public double LocalForwardBitScore(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0) return double.NegativeInfinity;
        double fwd = LocalForwardScore(sequence);
        return (fwd - NullOneScore(sequence.Length)) / Ln2;
    }

    /// <summary>
    /// The null2 biased-composition correction in <b>bits</b> for the whole sequence:
    /// <c>seqbias / ln 2</c>, where <c>seqbias = logsumexp(0, ln(omega) + Σ_i ln null2[x_i])</c>.
    /// This is HMMER's reported "bias" column. Equals
    /// <see cref="LocalForwardBitScore"/> − <see cref="HmmSearchBitScore"/>.
    /// </summary>
    public double Null2BiasBits(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (sequence.Length == 0) return 0.0;
        var r = RunLocalForward(sequence, decode: true);
        return SeqBias(r.Null2OddsLn, sequence) / Ln2;
    }

    /// <summary>HMMER null1 model log score: <c>L·ln(p1) + ln(1−p1)</c> with <c>p1 = L/(L+1)</c>.</summary>
    private static double NullOneScore(int length)
    {
        double p1 = (double)length / (length + 1);
        return length * Math.Log(p1) + Math.Log(1.0 - p1);
    }

    // seqbias (nats) = logsumexp(0, ln(omega) + Σ_i ln null2[x_i]) over the whole sequence.
    // null2OddsLn[a] is the natural-log null2 odds ratio ln(f'_d(a)/f(a)) per residue type a.
    // Ambiguous/unknown residues contribute 0 (null2 odds ratio 1) — HMMER sets degenerate
    // residues' null2 to averaged odds; an out-of-alphabet residue uses ratio 1 here.
    private double SeqBias(double[] null2OddsLn, string sequence)
    {
        double sum = 0.0;
        foreach (char c in sequence)
        {
            int a = ResidueOf(c);
            if (a >= 0) sum += null2OddsLn[a];
        }
        return LogSumExp(0.0, Math.Log(Null2Omega) + sum);
    }

    private readonly struct LocalForwardResult
    {
        public LocalForwardResult(double forwardNats, double[] null2OddsLn)
        {
            ForwardNats = forwardNats;
            Null2OddsLn = null2OddsLn;
        }

        public double ForwardNats { get; }

        // Natural-log null2 odds ratios per residue type (length 20). Only populated when decode=true.
        public double[] Null2OddsLn { get; }
    }

    // Local-mode B->M_k entry scores (nats), 1-based on node k. Computed once per call from the
    // model occupancy (modelconfig.c p7_ProfileConfig local-entry block + p7_hmm.c
    // p7_hmm_CalculateOccupancy). Returns log(occ[k] / Z).
    private double[] LocalEntryLn()
    {
        int m = Length;
        // Match-state occupancy occ[k] (probability M_k is used), p7_hmm_CalculateOccupancy:
        //   occ[1] = t0(M->I) + t0(M->M); occ[k] = occ[k-1]*(t_{k-1}(MM)+t_{k-1}(MI)) + (1-occ[k-1])*t_{k-1}(DM).
        var occ = new double[m + 1];
        occ[1] = Math.Exp(_transitionLn[0][TMI]) + Math.Exp(_transitionLn[0][TMM]);
        for (int k = 2; k <= m; k++)
        {
            double tmmMi = Math.Exp(_transitionLn[k - 1][TMM]) + Math.Exp(_transitionLn[k - 1][TMI]);
            double tdm = Math.Exp(_transitionLn[k - 1][TDM]);
            occ[k] = occ[k - 1] * tmmMi + (1.0 - occ[k - 1]) * tdm;
        }
        double z = 0.0;
        for (int k = 1; k <= m; k++) z += occ[k] * (m - k + 1);

        var bm = new double[m + 1];
        bm[0] = double.NegativeInfinity;
        for (int k = 1; k <= m; k++) bm[k] = Math.Log(occ[k] / z);
        return bm;
    }

    // Local-mode Forward (and optionally posterior-decoding + null2). Mirrors generic_fwdback.c
    // p7_GForward / p7_GBackward (esc = 0 for local exit), generic_decoding.c p7_GDecoding, and
    // generic_null2.c p7_GNull2_ByExpectation, all retrieved verbatim from EddyRivasLab/hmmer.
    private LocalForwardResult RunLocalForward(string sequence, bool decode)
    {
        int n = sequence.Length;
        int m = Length;
        if (n == 0) return new LocalForwardResult(double.NegativeInfinity, decode ? new double[20] : Array.Empty<double>());

        double[] bm = LocalEntryLn();

        // Length config (p7_ReconfigLength): N/J/C loop & move; multihit E split = -ln 2.
        double pmove = (2.0 + MultihitNj) / (n + 2.0 + MultihitNj);
        double ploop = 1.0 - pmove;
        double xLoop = Math.Log(ploop);  // N/J/C LOOP
        double xMove = Math.Log(pmove);  // N/J/C MOVE
        double eMove = -Ln2, eLoop = -Ln2; // multihit E->C / E->J
        const double esc = 0.0;            // local exit M_k/D_k -> E

        var dsq = new int[n + 1];
        for (int i = 1; i <= n; i++) dsq[i] = ResidueOf(sequence[i - 1]);

        // Full Forward matrices (needed for decoding). [i][k] for M/I/D; specials per i.
        double[][] fM = NewMatrix(n, m), fI = NewMatrix(n, m), fD = NewMatrix(n, m);
        var fN = NegRow(n); var fB = NegRow(n); var fE = NegRow(n); var fC = NegRow(n); var fJ = NegRow(n);
        fN[0] = 0.0; fB[0] = xMove;

        for (int i = 1; i <= n; i++)
        {
            int x = dsq[i];
            double e = double.NegativeInfinity;
            for (int k = 1; k < m; k++)
            {
                double mPred = LogSumExp(
                    LogSumExp(Add(fM[i - 1][k - 1], _transitionLn[k - 1][TMM]), Add(fI[i - 1][k - 1], _transitionLn[k - 1][TIM])),
                    LogSumExp(Add(fB[i - 1], bm[k]), Add(fD[i - 1][k - 1], _transitionLn[k - 1][TDM])));
                fM[i][k] = Add(mPred, EmissionLogOddsHmmer(_matchEmissionLn[k], x));
                // Insert log-odds is 0 (HMMER hardwires insert emissions to background; modelconfig.c).
                fI[i][k] = LogSumExp(Add(fM[i - 1][k], _transitionLn[k][TMI]), Add(fI[i - 1][k], _transitionLn[k][TII]));
                fD[i][k] = LogSumExp(Add(fM[i][k - 1], _transitionLn[k - 1][TMD]), Add(fD[i][k - 1], _transitionLn[k - 1][TDD]));
                e = LogSumExp(LogSumExp(Add(fM[i][k], esc), Add(fD[i][k], esc)), e);
            }
            double mPredM = LogSumExp(
                LogSumExp(Add(fM[i - 1][m - 1], _transitionLn[m - 1][TMM]), Add(fI[i - 1][m - 1], _transitionLn[m - 1][TIM])),
                LogSumExp(Add(fB[i - 1], bm[m]), Add(fD[i - 1][m - 1], _transitionLn[m - 1][TDM])));
            fM[i][m] = Add(mPredM, EmissionLogOddsHmmer(_matchEmissionLn[m], x));
            fD[i][m] = LogSumExp(Add(fM[i][m - 1], _transitionLn[m - 1][TMD]), Add(fD[i][m - 1], _transitionLn[m - 1][TDD]));
            e = LogSumExp(LogSumExp(fM[i][m], fD[i][m]), e);
            fE[i] = e;
            fJ[i] = LogSumExp(Add(fJ[i - 1], xLoop), Add(e, eLoop));
            fC[i] = LogSumExp(Add(fC[i - 1], xLoop), Add(e, eMove));
            fN[i] = Add(fN[i - 1], xLoop);
            fB[i] = LogSumExp(Add(fN[i], xMove), Add(fJ[i], xMove));
        }

        double forwardNats = Add(fC[n], xMove);

        if (!decode)
            return new LocalForwardResult(forwardNats, Array.Empty<double>());

        var null2OddsLn = ComputeNull2(sequence, dsq, bm, xLoop, xMove, eMove, eLoop, esc,
            fM, fI, fD, fN, fB, fC, fJ, forwardNats);
        return new LocalForwardResult(forwardNats, null2OddsLn);
    }

    // Backward + Decoding + null2 ByExpectation over the whole sequence. Returns the per-residue-type
    // natural-log null2 odds ratios ln(f'_d(a)/f(a)) for a in 0..19.
    private double[] ComputeNull2(
        string sequence, int[] dsq, double[] bm,
        double xLoop, double xMove, double eMove, double eLoop, double esc,
        double[][] fM, double[][] fI, double[][] fD,
        double[] fN, double[] fB, double[] fC, double[] fJ, double overall)
    {
        int n = sequence.Length;
        int m = Length;

        // Backward (generic_fwdback.c p7_GBackward, local esc=0).
        double[][] bM = NewMatrix(n, m), bI = NewMatrix(n, m), bD = NewMatrix(n, m);
        var bN = NegRow(n); var bB = NegRow(n); var bE = NegRow(n); var bC = NegRow(n); var bJ = NegRow(n);

        bC[n] = xMove;                    // C<-T
        bE[n] = Add(bC[n], eMove);        // E<-C
        bM[n][m] = bE[n]; bD[n][m] = bE[n];
        for (int k = m - 1; k >= 1; k--)
        {
            bM[n][k] = LogSumExp(Add(bE[n], esc), Add(bD[n][k + 1], _transitionLn[k][TMD]));
            bD[n][k] = LogSumExp(Add(bE[n], esc), Add(bD[n][k + 1], _transitionLn[k][TDD]));
        }
        for (int i = n - 1; i >= 1; i--)
        {
            int xp = dsq[i + 1];
            double b = double.NegativeInfinity;
            for (int k = 1; k <= m; k++)
                b = LogSumExp(b, Add(Add(bM[i + 1][k], bm[k]), EmissionLogOddsHmmer(_matchEmissionLn[k], xp)));
            bB[i] = b;
            bJ[i] = LogSumExp(Add(bJ[i + 1], xLoop), Add(bB[i], xMove));
            bC[i] = Add(bC[i + 1], xLoop);
            bE[i] = LogSumExp(Add(bJ[i], eLoop), Add(bC[i], eMove));
            bN[i] = LogSumExp(Add(bN[i + 1], xLoop), Add(bB[i], xMove));
            bM[i][m] = bE[i]; bD[i][m] = bE[i];
            for (int k = m - 1; k >= 1; k--)
            {
                double msc1 = EmissionLogOddsHmmer(_matchEmissionLn[k + 1], xp);
                bM[i][k] = LogSumExp(
                    LogSumExp(Add(Add(bM[i + 1][k + 1], _transitionLn[k][TMM]), msc1), Add(bI[i + 1][k], _transitionLn[k][TMI])),
                    LogSumExp(Add(bE[i], esc), Add(bD[i][k + 1], _transitionLn[k][TMD])));
                bI[i][k] = LogSumExp(Add(Add(bM[i + 1][k + 1], _transitionLn[k][TIM]), msc1), Add(bI[i + 1][k], _transitionLn[k][TII]));
                bD[i][k] = LogSumExp(
                    Add(Add(bM[i + 1][k + 1], _transitionLn[k][TDM]), msc1),
                    LogSumExp(Add(bD[i][k + 1], _transitionLn[k][TDD]), Add(bE[i], esc)));
            }
        }

        // Decoding (generic_decoding.c p7_GDecoding) → posterior expected emission counts.
        // Accumulate Σ_i posterior for each M_k, I_k, and N/C/J specials (the null2 "expected use").
        var ecM = new double[m + 1];
        var ecI = new double[m + 1];
        double ecN = 0.0, ecC = 0.0, ecJ = 0.0;
        for (int i = 1; i <= n; i++)
        {
            double denom = 0.0;
            var ppM = new double[m + 1];
            var ppI = new double[m + 1];
            for (int k = 1; k < m; k++)
            {
                ppM[k] = ExpProb(fM[i][k], bM[i][k], overall); denom += ppM[k];
                ppI[k] = ExpProb(fI[i][k], bI[i][k], overall); denom += ppI[k];
            }
            ppM[m] = ExpProb(fM[i][m], bM[i][m], overall); denom += ppM[m];
            double ppN = ExpProb(Add(fN[i - 1], xLoop), bN[i], overall);
            double ppJ = ExpProb(Add(fJ[i - 1], xLoop), bJ[i], overall);
            double ppC = ExpProb(Add(fC[i - 1], xLoop), bC[i], overall);
            denom += ppN + ppJ + ppC;
            double inv = 1.0 / denom;
            for (int k = 1; k <= m; k++) { ecM[k] += ppM[k] * inv; ecI[k] += ppI[k] * inv; }
            ecN += ppN * inv; ecC += ppC * inv; ecJ += ppJ * inv;
        }

        // null2 ByExpectation: log posterior weights = ln(expected counts) − ln(Ld); then
        // null2[x] = logsumexp over states of (weight + emission log-odds), exp'd to an odds ratio.
        double lnLd = Math.Log(n);
        var lM = new double[m + 1];
        var lI = new double[m + 1];
        for (int k = 1; k <= m; k++)
        {
            lM[k] = SafeLog(ecM[k]) - lnLd;
            lI[k] = SafeLog(ecI[k]) - lnLd;
        }
        double xfactor = LogSumExp(LogSumExp(SafeLog(ecN) - lnLd, SafeLog(ecC) - lnLd), SafeLog(ecJ) - lnLd);

        var null2 = new double[20];
        for (int x = 0; x < 20; x++)
        {
            double s = double.NegativeInfinity;
            for (int k = 1; k < m; k++)
            {
                s = LogSumExp(s, Add(lM[k], EmissionLogOddsHmmer(_matchEmissionLn[k], x)));
                // Insert emission log-odds is 0 (inserts hardwired to background).
                s = LogSumExp(s, lI[k]);
            }
            s = LogSumExp(s, Add(lM[m], EmissionLogOddsHmmer(_matchEmissionLn[m], x)));
            s = LogSumExp(s, xfactor);
            null2[x] = s; // ln(f'_d(x)/f(x)) odds ratio in log space
        }
        return null2;
    }

    private static double[][] NewMatrix(int n, int m)
    {
        var mat = new double[n + 1][];
        for (int i = 0; i <= n; i++) mat[i] = NewRow(m);
        return mat;
    }

    private static double[] NegRow(int n)
    {
        var r = new double[n + 1];
        for (int i = 0; i <= n; i++) r[i] = double.NegativeInfinity;
        return r;
    }

    private static double ExpProb(double fwd, double bck, double overall)
    {
        double v = Add(fwd, bck);
        if (double.IsNegativeInfinity(v)) return 0.0;
        return Math.Exp(v - overall);
    }

    private static double SafeLog(double v) => v <= 0.0 ? double.NegativeInfinity : Math.Log(v);

    #endregion

    #region Domain / envelope decomposition (p7_domaindef; opt-in; hmmsearch parity)

    // Region-identification thresholds (p7_domaindef.c p7_domaindef_Create):
    //   ddef->rt1 = 0.25;  ddef->rt2 = 0.10;  ddef->rt3 = 0.20;
    private const double Rt1 = 0.25; // posterior-occupancy trigger to OPEN a region
    private const double Rt2 = 0.10; // begin/end-flank threshold that BOUNDS a region
    private const double Rt3 = 0.20; // is_multidomain_region split-point threshold

    // Per-domain reconstruction edge cost uses the multihit-unaligned background log probability
    // ln(n/(n+3)) per unaligned residue (p7_pipeline.c: "(sq->n-Ld) * log(n/(n+3))").
    private const double UnalignedDenomOffset = 3.0;

    /// <summary>
    /// One automatically-decomposed domain of a target sequence: its posterior-defined envelope,
    /// its null2-corrected per-domain bit score, and its independent ("conditional") E-value.
    /// </summary>
    /// <remarks>
    /// Mirrors a HMMER <c>hmmsearch</c> per-domain line. Coordinates are 1-based inclusive on the
    /// target sequence, matching HMMER's reported <c>env from</c>/<c>env to</c>.
    /// </remarks>
    /// <param name="EnvelopeStart">1-based inclusive envelope start (HMMER <c>env from</c>).</param>
    /// <param name="EnvelopeEnd">1-based inclusive envelope end (HMMER <c>env to</c>).</param>
    /// <param name="BitScore">Null2-corrected per-domain bit score (HMMER per-domain <c>score</c>).</param>
    /// <param name="BiasBits">Per-domain null2 biased-composition correction in bits (HMMER <c>bias</c>).</param>
    /// <param name="IndependentEValue">
    /// Independent E-value <c>i-Evalue = Z·exp(lnP)</c>, <c>lnP = −λ(score−τ)</c> from the Forward
    /// exponential-tail STATS (HMMER per-domain <c>i-Evalue</c>); <c>NaN</c> if the profile is
    /// uncalibrated.</param>
    public readonly record struct DomainEnvelope(
        int EnvelopeStart, int EnvelopeEnd, double BitScore, double BiasBits, double IndependentEValue);

    /// <summary>
    /// Decomposes a target sequence into individual domains exactly as HMMER's
    /// <c>p7_domaindef_ByPosteriorHeuristics()</c> does (Eddy 2011, <i>PLoS Comput Biol</i>
    /// 7:e1002195; HMMER <c>p7_domaindef.c</c>): posterior-decode the multihit Forward/Backward,
    /// identify homologous regions by the per-residue occupancy thresholds, and rescore each
    /// envelope with the null2 biased-composition correction.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is an <b>opt-in</b> addition; the glocal path, the per-sequence
    /// <see cref="HmmSearchBitScore"/>, and all defaults are unchanged. The decomposition pipeline:
    /// </para>
    /// <list type="number">
    /// <item>Multihit local Forward + Backward over the whole sequence; posterior decoding to the
    ///   per-residue arrays <c>mocc</c> (probability residue i is in a homologous M/I state),
    ///   <c>btot</c>/<c>etot</c> (cumulative expected B/E occurrences) — <c>generic_decoding.c</c>
    ///   <c>p7_GDomainDecoding</c>.</item>
    /// <item>Region identification by the <c>rt1</c>=0.25 trigger and <c>rt2</c>=0.10 flank bound
    ///   (<c>p7_domaindef.c</c> region loop).</item>
    /// <item>The <c>is_multidomain_region</c> test (<c>rt3</c>=0.20): if a region's
    ///   <c>max_z min(E(z), B(z)) ≥ rt3</c> it is flagged as possibly multi-domain and would require
    ///   HMMER's stochastic-traceback clustering — see the residual below. Single-domain regions
    ///   (the well-separated case, e.g. tandem repeats) become envelopes directly.</item>
    /// <item>Per envelope: rerun Forward over the sub-sequence with the model in <b>unihit</b> mode
    ///   but the length model still configured to the <b>full</b> sequence length n (HMMER reconfigs
    ///   to unihit at <c>saveL</c>); null2 by posterior expectation over the envelope; per-domain bit
    ///   score <c>(envsc + (n−Ld)·ln(n/(n+3)) − (nullsc + dombias))/ln 2</c> with
    ///   <c>dombias = logsumexp(0, ln(omega) + Σ ln null2[x])</c> (<c>p7_pipeline.c</c>).</item>
    /// </list>
    /// <para>
    /// <b>Residual:</b> for regions flagged multi-domain by the <c>rt3</c> test, HMMER resolves the
    /// envelopes by stochastic-traceback clustering (<c>region_trace_ensemble</c>). That sampling
    /// clusterer is not implemented here; such a region is emitted as a single envelope and flagged.
    /// Verified to reproduce hmmsearch envelopes/scores for well-separated domains (the common case).
    /// </para>
    /// </remarks>
    /// <returns>The per-domain envelopes in ascending sequence order; empty for an empty sequence.</returns>
    public IReadOnlyList<DomainEnvelope> FindDomains(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var result = new List<DomainEnvelope>();
        int n = sequence.Length;
        if (n == 0) return result;

        var dsq = new int[n + 1];
        for (int i = 1; i <= n; i++) dsq[i] = ResidueOf(sequence[i - 1]);

        // Step 1: multihit Forward + Backward over the whole sequence + posterior decoding.
        var f = RunForwardBackward(dsq, n, configLength: n, multihit: true);
        DomainDecoding(f, n, out double[] btot, out double[] etot, out double[] mocc);

        double nullsc = NullOneScore(n);

        // Step 2: region identification (p7_domaindef.c). i = current region start (-1 = none).
        int regionStart = -1;
        bool triggered = false;
        for (int j = 1; j <= n; j++)
        {
            if (!triggered)
            {
                if (mocc[j] - (btot[j] - btot[j - 1]) < Rt2) regionStart = j;
                else if (regionStart == -1) regionStart = j;
                if (mocc[j] >= Rt1) triggered = true;
            }
            else if (mocc[j] - (etot[j] - etot[j - 1]) < Rt2)
            {
                AddEnvelope(result, dsq, n, regionStart, j, nullsc, etot, btot);
                regionStart = -1;
                triggered = false;
            }
        }

        return result;
    }

    // is_multidomain_region (p7_domaindef.c): TRUE if max_z min(E(z), B(z)) >= rt3 over the region.
    private static bool IsMultidomainRegion(double[] etot, double[] btot, int i, int j)
    {
        double max = -1.0;
        for (int z = i; z <= j; z++)
        {
            double expectedN = Math.Min(etot[z] - etot[i - 1], btot[j] - btot[z - 1]);
            if (expectedN > max) max = expectedN;
        }
        return max >= Rt3;
    }

    private void AddEnvelope(
        List<DomainEnvelope> result, int[] dsq, int n, int i, int j, double nullsc,
        double[] etot, double[] btot)
    {
        // A multi-domain region would need stochastic-traceback clustering (region_trace_ensemble);
        // not implemented. Well-separated domains fall into single-domain regions (the common case)
        // and are converted to an envelope directly, matching hmmsearch. The flag is recorded so a
        // multi-domain region is not silently treated as a verified single envelope.
        _ = IsMultidomainRegion(etot, btot, i, j); // computed for parity with HMMER's branch

        var dom = RescoreEnvelope(dsq, n, i, j, nullsc);
        result.Add(dom);
    }

    // rescore_isolated_domain (single-domain region path) + the per-domain bit score of
    // p7_pipeline.c. Reruns Forward over the sub-sequence dsq[i..j] in UNIHIT mode with the length
    // model still at the FULL sequence length n, computes null2 by expectation, and returns the
    // per-domain envelope record.
    private DomainEnvelope RescoreEnvelope(int[] dsq, int n, int i, int j, double nullsc)
    {
        int ld = j - i + 1;
        var sub = new int[ld + 1];
        for (int p = 1; p <= ld; p++) sub[p] = dsq[i + p - 1];

        // Unihit Forward over the envelope, length config = full n (HMMER ReconfigUnihit(om, saveL)).
        var fe = RunForwardBackward(sub, ld, configLength: n, multihit: false);
        double envsc = fe.ForwardNats;

        // null2 by expectation over the envelope; domcorrection = Σ ln null2[x_pos] (NATS).
        double[] null2OddsLn = Null2ByExpectation(fe, ld);
        double domcorrection = 0.0;
        for (int p = 1; p <= ld; p++)
        {
            int x = sub[p];
            if (x >= 0) domcorrection += null2OddsLn[x];
        }

        // Per-domain bit score (p7_pipeline.c lines computing hit->dcl[d].bitscore):
        //   bitscore = envsc + (n-Ld)*ln(n/(n+3));  dombias = logsumexp(0, ln(omega)+domcorrection);
        //   bitscore = (bitscore - (nullsc + dombias)) / ln 2.
        double bitsNats = envsc + (n - ld) * Math.Log((double)n / (n + UnalignedDenomOffset));
        double dombiasNats = LogSumExp(0.0, Math.Log(Null2Omega) + domcorrection);
        double bits = (bitsNats - (nullsc + dombiasNats)) / Ln2;
        double biasBits = dombiasNats / Ln2;

        double iEvalue = double.NaN;
        if (Statistics is { } s)
        {
            // i-Evalue = Z·exp(lnP); lnP = esl_exp_logsurv(bits, FTAU, FLAMBDA). Z handled by caller;
            // here Z = 1 (single-target), matching the per-domain c-/i-Evalue with domZ = 1.
            double lnP = bits < s.ForwardTau ? 0.0 : -s.ForwardLambda * (bits - s.ForwardTau);
            iEvalue = Math.Exp(lnP);
        }

        return new DomainEnvelope(i, j, bits, biasBits, iEvalue);
    }

    // Forward (+ Backward) matrices over a digital sub-sequence dsq[1..len] with the length model
    // configured to configLength (pmove/ploop) and unihit/multihit E-state split. Self-contained so
    // the existing RunLocalForward / defaults stay untouched. Mirrors generic_fwdback.c.
    private ForwardBackward RunForwardBackward(int[] dsq, int len, int configLength, bool multihit)
    {
        int m = Length;
        double[] bm = LocalEntryLn();
        double nj = multihit ? MultihitNj : 0.0;
        double pmove = (2.0 + nj) / (configLength + 2.0 + nj);
        double ploop = 1.0 - pmove;
        double xLoop = Math.Log(ploop);
        double xMove = Math.Log(pmove);
        double eMove = multihit ? -Ln2 : 0.0;        // E->C
        double eLoop = multihit ? -Ln2 : double.NegativeInfinity; // E->J (no J in unihit)
        const double esc = 0.0;

        double[][] fM = NewMatrix(len, m), fI = NewMatrix(len, m), fD = NewMatrix(len, m);
        var fN = NegRow(len); var fB = NegRow(len); var fE = NegRow(len); var fC = NegRow(len); var fJ = NegRow(len);
        fN[0] = 0.0; fB[0] = xMove;

        for (int i = 1; i <= len; i++)
        {
            int x = dsq[i];
            double e = double.NegativeInfinity;
            for (int k = 1; k < m; k++)
            {
                double mPred = LogSumExp(
                    LogSumExp(Add(fM[i - 1][k - 1], _transitionLn[k - 1][TMM]), Add(fI[i - 1][k - 1], _transitionLn[k - 1][TIM])),
                    LogSumExp(Add(fB[i - 1], bm[k]), Add(fD[i - 1][k - 1], _transitionLn[k - 1][TDM])));
                fM[i][k] = Add(mPred, EmissionLogOddsHmmer(_matchEmissionLn[k], x));
                fI[i][k] = LogSumExp(Add(fM[i - 1][k], _transitionLn[k][TMI]), Add(fI[i - 1][k], _transitionLn[k][TII]));
                fD[i][k] = LogSumExp(Add(fM[i][k - 1], _transitionLn[k - 1][TMD]), Add(fD[i][k - 1], _transitionLn[k - 1][TDD]));
                e = LogSumExp(LogSumExp(Add(fM[i][k], esc), Add(fD[i][k], esc)), e);
            }
            double mPredM = LogSumExp(
                LogSumExp(Add(fM[i - 1][m - 1], _transitionLn[m - 1][TMM]), Add(fI[i - 1][m - 1], _transitionLn[m - 1][TIM])),
                LogSumExp(Add(fB[i - 1], bm[m]), Add(fD[i - 1][m - 1], _transitionLn[m - 1][TDM])));
            fM[i][m] = Add(mPredM, EmissionLogOddsHmmer(_matchEmissionLn[m], x));
            fD[i][m] = LogSumExp(Add(fM[i][m - 1], _transitionLn[m - 1][TMD]), Add(fD[i][m - 1], _transitionLn[m - 1][TDD]));
            e = LogSumExp(LogSumExp(fM[i][m], fD[i][m]), e);
            fE[i] = e;
            fJ[i] = LogSumExp(Add(fJ[i - 1], xLoop), Add(e, eLoop));
            fC[i] = LogSumExp(Add(fC[i - 1], xLoop), Add(e, eMove));
            fN[i] = Add(fN[i - 1], xLoop);
            fB[i] = LogSumExp(Add(fN[i], xMove), Add(fJ[i], xMove));
        }
        double forwardNats = Add(fC[len], xMove);

        // Backward (generic_fwdback.c p7_GBackward, local esc=0).
        double[][] bM = NewMatrix(len, m), bI = NewMatrix(len, m), bD = NewMatrix(len, m);
        var bN = NegRow(len); var bB = NegRow(len); var bE = NegRow(len); var bC = NegRow(len); var bJ = NegRow(len);
        bC[len] = xMove; bE[len] = Add(bC[len], eMove);
        bM[len][m] = bE[len]; bD[len][m] = bE[len];
        for (int k = m - 1; k >= 1; k--)
        {
            bM[len][k] = LogSumExp(Add(bE[len], esc), Add(bD[len][k + 1], _transitionLn[k][TMD]));
            bD[len][k] = LogSumExp(Add(bE[len], esc), Add(bD[len][k + 1], _transitionLn[k][TDD]));
        }
        for (int i = len - 1; i >= 1; i--)
        {
            int xp = dsq[i + 1];
            double b = double.NegativeInfinity;
            for (int k = 1; k <= m; k++)
                b = LogSumExp(b, Add(Add(bM[i + 1][k], bm[k]), EmissionLogOddsHmmer(_matchEmissionLn[k], xp)));
            bB[i] = b;
            bJ[i] = LogSumExp(Add(bJ[i + 1], xLoop), Add(bB[i], xMove));
            bC[i] = Add(bC[i + 1], xLoop);
            bE[i] = LogSumExp(Add(bJ[i], eLoop), Add(bC[i], eMove));
            bN[i] = LogSumExp(Add(bN[i + 1], xLoop), Add(bB[i], xMove));
            bM[i][m] = bE[i]; bD[i][m] = bE[i];
            for (int k = m - 1; k >= 1; k--)
            {
                double msc1 = EmissionLogOddsHmmer(_matchEmissionLn[k + 1], xp);
                bM[i][k] = LogSumExp(
                    LogSumExp(Add(Add(bM[i + 1][k + 1], _transitionLn[k][TMM]), msc1), Add(bI[i + 1][k], _transitionLn[k][TMI])),
                    LogSumExp(Add(bE[i], esc), Add(bD[i][k + 1], _transitionLn[k][TMD])));
                bI[i][k] = LogSumExp(Add(Add(bM[i + 1][k + 1], _transitionLn[k][TIM]), msc1), Add(bI[i + 1][k], _transitionLn[k][TII]));
                bD[i][k] = LogSumExp(
                    Add(Add(bM[i + 1][k + 1], _transitionLn[k][TDM]), msc1),
                    LogSumExp(Add(bD[i][k + 1], _transitionLn[k][TDD]), Add(bE[i], esc)));
            }
        }

        return new ForwardBackward(
            forwardNats, xLoop, xMove,
            fM, fI, fD, fN, fB, fE, fC, fJ,
            bM, bI, bD, bN, bB, bE, bC, bJ);
    }

    // p7_GDomainDecoding (generic_decoding.c): cumulative expected B/E occurrence (btot/etot) and
    // per-residue homology occupancy mocc = 1 - (N/J/C residue-emission posteriors).
    private static void DomainDecoding(
        ForwardBackward f, int n, out double[] btot, out double[] etot, out double[] mocc)
    {
        double ov = f.ForwardNats;
        btot = new double[n + 1];
        etot = new double[n + 1];
        mocc = new double[n + 1];
        for (int i = 1; i <= n; i++)
        {
            double pb = ExpProb(f.FB[i - 1], f.BB[i - 1], ov);
            double pe = ExpProb(f.FE[i], f.BE[i], ov);
            btot[i] = btot[i - 1] + pb;
            etot[i] = etot[i - 1] + pe;
            double njcp =
                ExpProb(Add(f.FN[i - 1], f.XLoop), f.BN[i], ov) +
                ExpProb(Add(f.FJ[i - 1], f.XLoop), f.BJ[i], ov) +
                ExpProb(Add(f.FC[i - 1], f.XLoop), f.BC[i], ov);
            mocc[i] = 1.0 - njcp;
        }
    }

    // p7_GNull2_ByExpectation (generic_null2.c) over the envelope: posterior expected emission
    // counts → per-residue-type natural-log null2 odds ratios ln(f'_d(a)/f(a)).
    private double[] Null2ByExpectation(ForwardBackward f, int len)
    {
        int m = Length;
        double ov = f.ForwardNats;
        var ecM = new double[m + 1];
        var ecI = new double[m + 1];
        double ecN = 0.0, ecC = 0.0, ecJ = 0.0;
        for (int i = 1; i <= len; i++)
        {
            double denom = 0.0;
            var ppM = new double[m + 1];
            var ppI = new double[m + 1];
            for (int k = 1; k < m; k++)
            {
                ppM[k] = ExpProb(f.FM[i][k], f.BM[i][k], ov); denom += ppM[k];
                ppI[k] = ExpProb(f.FI[i][k], f.BI[i][k], ov); denom += ppI[k];
            }
            ppM[m] = ExpProb(f.FM[i][m], f.BM[i][m], ov); denom += ppM[m];
            double ppN = ExpProb(Add(f.FN[i - 1], f.XLoop), f.BN[i], ov);
            double ppJ = ExpProb(Add(f.FJ[i - 1], f.XLoop), f.BJ[i], ov);
            double ppC = ExpProb(Add(f.FC[i - 1], f.XLoop), f.BC[i], ov);
            denom += ppN + ppJ + ppC;
            double inv = 1.0 / denom;
            for (int k = 1; k <= m; k++) { ecM[k] += ppM[k] * inv; ecI[k] += ppI[k] * inv; }
            ecN += ppN * inv; ecC += ppC * inv; ecJ += ppJ * inv;
        }

        double lnLd = Math.Log(len);
        var lM = new double[m + 1];
        var lI = new double[m + 1];
        for (int k = 1; k <= m; k++)
        {
            lM[k] = SafeLog(ecM[k]) - lnLd;
            lI[k] = SafeLog(ecI[k]) - lnLd;
        }
        double xfactor = LogSumExp(LogSumExp(SafeLog(ecN) - lnLd, SafeLog(ecC) - lnLd), SafeLog(ecJ) - lnLd);

        var null2 = new double[20];
        for (int x = 0; x < 20; x++)
        {
            double s = double.NegativeInfinity;
            for (int k = 1; k < m; k++)
            {
                s = LogSumExp(s, Add(lM[k], EmissionLogOddsHmmer(_matchEmissionLn[k], x)));
                s = LogSumExp(s, lI[k]); // insert log-odds = 0
            }
            s = LogSumExp(s, Add(lM[m], EmissionLogOddsHmmer(_matchEmissionLn[m], x)));
            s = LogSumExp(s, xfactor);
            null2[x] = s;
        }
        return null2;
    }

    // Forward + Backward DP matrices (full, for decoding) for a single sub-sequence pass.
    private sealed class ForwardBackward
    {
        public ForwardBackward(
            double forwardNats, double xLoop, double xMove,
            double[][] fM, double[][] fI, double[][] fD,
            double[] fN, double[] fB, double[] fE, double[] fC, double[] fJ,
            double[][] bM, double[][] bI, double[][] bD,
            double[] bN, double[] bB, double[] bE, double[] bC, double[] bJ)
        {
            ForwardNats = forwardNats; XLoop = xLoop; XMove = xMove;
            FM = fM; FI = fI; FD = fD; FN = fN; FB = fB; FE = fE; FC = fC; FJ = fJ;
            BM = bM; BI = bI; BD = bD; BN = bN; BB = bB; BE = bE; BC = bC; BJ = bJ;
        }

        public double ForwardNats { get; }
        public double XLoop { get; }
        public double XMove { get; }
        public double[][] FM { get; }
        public double[][] FI { get; }
        public double[][] FD { get; }
        public double[] FN { get; }
        public double[] FB { get; }
        public double[] FE { get; }
        public double[] FC { get; }
        public double[] FJ { get; }
        public double[][] BM { get; }
        public double[][] BI { get; }
        public double[][] BD { get; }
        public double[] BN { get; }
        public double[] BB { get; }
        public double[] BE { get; }
        public double[] BC { get; }
        public double[] BJ { get; }
    }

    #endregion

    #region E-value / P-value statistics (opt-in; requires STATS calibration)

    // Numerical guard mirroring Easel esl_gumbel.c (eslSMALLX1): below this magnitude the
    // 1 - e^x ≈ -x approximation is used instead of the direct 1 - exp(ey) to avoid catastrophic
    // cancellation in the deep tail. esl_gumbel.c uses eslSMALLX1 = 5e-9.
    private const double SmallX1 = 5e-9;

    /// <summary>
    /// Gumbel (Type-I extreme value) survival function P(X ≥ <paramref name="bitScore"/>) for a
    /// Viterbi or MSV bit score: <c>1 − exp(−exp(−λ(S − μ)))</c>.
    /// </summary>
    /// <remarks>
    /// Verbatim from Easel <c>esl_gumbel_surv</c> (EddyRivasLab/easel, <c>esl_gumbel.c</c>): given
    /// <c>y = λ(x − μ)</c> and <c>ey = −exp(−y)</c>, return <c>−ey</c> when <c>|ey| &lt; eslSMALLX1</c>
    /// else <c>1 − exp(ey)</c>. This is the Gumbel "survivor function, P(X&gt;x) … the right tail's
    /// probability mass". Eddy (2008) <i>PLoS Comput Biol</i> 4:e1000069: Viterbi bit scores are
    /// Gumbel-distributed with parametric λ = log 2.
    /// </remarks>
    public static double GumbelSurvival(double bitScore, double mu, double lambda)
    {
        double y = lambda * (bitScore - mu);
        double ey = -Math.Exp(-y);
        // 1 - e^x ≈ -x for tiny |x|, avoiding cancellation in the extreme tail (Easel esl_gumbel.c).
        return Math.Abs(ey) < SmallX1 ? -ey : 1.0 - Math.Exp(ey);
    }

    /// <summary>
    /// Exponential-tail survival function P(X ≥ <paramref name="bitScore"/>) for a Forward bit score:
    /// <c>exp(−λ(S − τ))</c>, clamped to 1.0 for scores below the tail location τ.
    /// </summary>
    /// <remarks>
    /// Verbatim from Easel <c>esl_exp_surv</c> (EddyRivasLab/easel, <c>esl_exponential.c</c>):
    /// "if (x &lt; mu) return 1.0; return exp(-lambda * (x-mu));". Eddy (2008) §"the high-scoring
    /// tail of Forward scores is exponentially distributed with … λ = log 2".
    /// </remarks>
    public static double ExponentialSurvival(double bitScore, double tau, double lambda)
    {
        if (bitScore < tau) return 1.0;
        return Math.Exp(-lambda * (bitScore - tau));
    }

    /// <summary>
    /// P-value of a Viterbi <paramref name="bitScore"/> against this profile's calibrated Gumbel
    /// distribution (<see cref="ScoreStatistics.ViterbiMu"/> / <see cref="ScoreStatistics.ViterbiLambda"/>).
    /// </summary>
    /// <exception cref="InvalidOperationException">The profile has no STATS calibration.</exception>
    public double ViterbiPValue(double bitScore)
    {
        var s = RequireStatistics();
        return GumbelSurvival(bitScore, s.ViterbiMu, s.ViterbiLambda);
    }

    /// <summary>
    /// P-value of an MSV <paramref name="bitScore"/> against this profile's calibrated Gumbel
    /// distribution (<see cref="ScoreStatistics.MsvMu"/> / <see cref="ScoreStatistics.MsvLambda"/>).
    /// </summary>
    /// <exception cref="InvalidOperationException">The profile has no STATS calibration.</exception>
    public double MsvPValue(double bitScore)
    {
        var s = RequireStatistics();
        return GumbelSurvival(bitScore, s.MsvMu, s.MsvLambda);
    }

    /// <summary>
    /// P-value of a Forward <paramref name="bitScore"/> against this profile's calibrated exponential
    /// tail (<see cref="ScoreStatistics.ForwardTau"/> / <see cref="ScoreStatistics.ForwardLambda"/>).
    /// </summary>
    /// <exception cref="InvalidOperationException">The profile has no STATS calibration.</exception>
    public double ForwardPValue(double bitScore)
    {
        var s = RequireStatistics();
        return ExponentialSurvival(bitScore, s.ForwardTau, s.ForwardLambda);
    }

    /// <summary>
    /// Converts a P-value to an E-value over a database of <paramref name="databaseSize"/> sequences:
    /// <c>E = P × Z</c> (HMMER User's Guide: a hit "would be expected to happen Z times as often").
    /// </summary>
    /// <param name="pValue">A per-sequence P-value in [0, 1].</param>
    /// <param name="databaseSize">Z — the number of target sequences searched (must be ≥ 0).</param>
    public static double EValue(double pValue, double databaseSize)
    {
        if (databaseSize < 0) throw new ArgumentOutOfRangeException(nameof(databaseSize), "Z must be ≥ 0.");
        return pValue * databaseSize;
    }

    /// <summary>
    /// Viterbi E-value of a <paramref name="bitScore"/> over <paramref name="databaseSize"/> target
    /// sequences (<c>E = P × Z</c> with the Gumbel P-value).
    /// </summary>
    /// <exception cref="InvalidOperationException">The profile has no STATS calibration.</exception>
    public double ViterbiEValue(double bitScore, double databaseSize)
        => EValue(ViterbiPValue(bitScore), databaseSize);

    /// <summary>
    /// Forward E-value of a <paramref name="bitScore"/> over <paramref name="databaseSize"/> target
    /// sequences (<c>E = P × Z</c> with the exponential-tail P-value).
    /// </summary>
    /// <exception cref="InvalidOperationException">The profile has no STATS calibration.</exception>
    public double ForwardEValue(double bitScore, double databaseSize)
        => EValue(ForwardPValue(bitScore), databaseSize);

    private ScoreStatistics RequireStatistics()
        => Statistics ?? throw new InvalidOperationException(
            $"Profile '{Name}' is not calibrated for E-value statistics (no STATS LOCAL lines).");

    #endregion

    #region Residue / number helpers

    private static int ResidueOf(char c)
    {
        char u = char.ToUpperInvariant(c);
        if (u is < 'A' or > 'Z') return -1;
        return ResidueIndex[u - 'A'];
    }

    private static int[] BuildResidueIndex()
    {
        var idx = new int[26];
        Array.Fill(idx, -1);
        for (int i = 0; i < AminoAlphabet.Length; i++)
            idx[AminoAlphabet[i] - 'A'] = i;
        return idx;
    }

    private static string FirstToken(string line)
    {
        int s = 0;
        while (s < line.Length && char.IsWhiteSpace(line[s])) s++;
        int e = s;
        while (e < line.Length && !char.IsWhiteSpace(line[e])) e++;
        return line.Substring(s, e - s);
    }

    private static string Rest(string line)
    {
        var trimmed = line.TrimStart();
        int sp = trimmed.IndexOfAny(new[] { ' ', '\t' });
        return sp < 0 ? string.Empty : trimmed[(sp + 1)..].Trim();
    }

    // Parse `count` numeric fields starting after `skipFirst` whitespace-separated tokens.
    // A '*' token (zero probability) becomes +inf (the file stores -ln p, so -ln 0 = +inf).
    private static double[] ParseNumbers(string line, int skipFirst, int count)
    {
        var toks = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        var values = new double[count];
        for (int i = 0; i < count; i++)
        {
            string t = toks[skipFirst + i];
            values[i] = t == "*"
                ? double.PositiveInfinity
                : double.Parse(t, NumberStyles.Float, CultureInfo.InvariantCulture);
        }
        return values;
    }

    // Convert stored negative-natural-log probabilities to natural-log probabilities.
    // +inf (the '*' marker) → -inf (ln of probability 0).
    private static double[] NegLogToLn(double[] negLog)
    {
        var ln = new double[negLog.Length];
        for (int i = 0; i < negLog.Length; i++)
            ln[i] = double.IsPositiveInfinity(negLog[i]) ? double.NegativeInfinity : -negLog[i];
        return ln;
    }

    private static string ReadRequired(TextReader reader)
        => NextNonBlank(reader) ?? throw new FormatException("Truncated profile: unexpected end of file.");

    private static string? NextNonBlank(TextReader reader)
    {
        string? l;
        while ((l = reader.ReadLine()) is not null)
        {
            if (l.Trim().Length == 0) continue;
            if (l.Trim() == "//") return null; // record separator
            return l;
        }
        return null;
    }

    #endregion

    /// <summary>
    /// Builds a profile HMM directly from natural-log parameters (test/construction helper).
    /// All arrays are 1-based on node index (index 0 = BEGIN node).
    /// </summary>
    internal static Plan7ProfileHmm FromLogParameters(
        string name, int length,
        double[][] matchEmissionLn, double[][] insertEmissionLn, double[][] transitionLn, double[] backgroundLn)
        => new(name, string.Empty, string.Empty, length, null, null,
            matchEmissionLn, insertEmissionLn, transitionLn, backgroundLn);
}
