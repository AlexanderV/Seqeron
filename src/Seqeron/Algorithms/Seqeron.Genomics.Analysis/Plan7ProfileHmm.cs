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

    // Natural-log probabilities (NOT negated): matchEmissionLn[k][a] = ln P(emit residue a | M_k).
    // k is 1-based (index 0 unused). -inf where the file stored '*'.
    private readonly double[][] _matchEmissionLn;   // [Length+1][20]
    private readonly double[][] _insertEmissionLn;  // [Length+1][20]  (node k insert I_k; index 0 = I_0)
    private readonly double[][] _transitionLn;       // [Length+1][7]   (node k; index 0 = BEGIN node)
    private readonly double[] _backgroundLn;         // [20] = ln COMPO background composition

    private Plan7ProfileHmm(
        string name, string accession, string description, int length, double? gatheringThreshold,
        double[][] matchEmissionLn, double[][] insertEmissionLn, double[][] transitionLn, double[] backgroundLn)
    {
        Name = name;
        Accession = accession;
        Description = description;
        Length = length;
        GatheringThreshold = gatheringThreshold;
        _matchEmissionLn = matchEmissionLn;
        _insertEmissionLn = insertEmissionLn;
        _transitionLn = transitionLn;
        _backgroundLn = backgroundLn;
    }

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

        return new Plan7ProfileHmm(name, accession, description, length, ga1,
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
        => new(name, string.Empty, string.Empty, length, null,
            matchEmissionLn, insertEmissionLn, transitionLn, backgroundLn);
}
