using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides algorithms for protein motif finding, domain prediction, and pattern matching.
/// </summary>
public static class ProteinMotifFinder
{
    #region Records and Types

    /// <summary>
    /// Represents a found protein motif.
    /// </summary>
    public readonly record struct MotifMatch(
        int Start,
        int End,
        string Sequence,
        string MotifName,
        string Pattern,
        double Score,
        double EValue);

    /// <summary>
    /// Represents a protein domain.
    /// </summary>
    public readonly record struct ProteinDomain(
        string Name,
        string Accession,
        int Start,
        int End,
        double Score,
        string Description);

    /// <summary>
    /// A profile-HMM domain hit carrying the Viterbi bit <see cref="Score"/> and the HMMER
    /// E-value derived from the profile's <c>STATS LOCAL VITERBI</c> Gumbel calibration.
    /// </summary>
    /// <param name="Name">Domain family name (e.g. "SH3").</param>
    /// <param name="Accession">Pfam accession (e.g. "PF00018").</param>
    /// <param name="Start">0-based start of the scored span (sequence start for the glocal score).</param>
    /// <param name="End">0-based inclusive end of the scored span.</param>
    /// <param name="Score">Viterbi log-odds score in bits.</param>
    /// <param name="EValue">Viterbi E-value E = P·Z (P from the profile's Gumbel STATS, Z = database size).</param>
    /// <param name="Description">Domain description.</param>
    public readonly record struct ProteinDomainHit(
        string Name,
        string Accession,
        int Start,
        int End,
        double Score,
        double EValue,
        string Description);

    /// <summary>
    /// Represents a PROSITE-style pattern.
    /// </summary>
    public readonly record struct PrositePattern(
        string Accession,
        string Name,
        string Pattern,
        string RegexPattern,
        string Description);

    /// <summary>
    /// Signal-peptide cleavage-site prediction result (von Heijne 1986 weight-matrix method).
    /// </summary>
    /// <param name="CleavagePosition">
    /// 1-based position of the first residue of the mature protein; cleavage occurs between
    /// residue <c>CleavagePosition − 1</c> (position −1) and <c>CleavagePosition</c> (position +1).
    /// </param>
    /// <param name="Score">
    /// von Heijne weight-matrix score at the best site: sum of log-odds weights over positions
    /// −13..+2 (natural log of observed/expected residue frequencies).
    /// </param>
    /// <param name="SignalSequence">Predicted signal peptide (residues 1..<c>CleavagePosition − 1</c>).</param>
    /// <param name="WindowSequence">The 15-residue scoring window (positions −13..+2) at the best site.</param>
    /// <param name="IsLikelySignalPeptide"><c>true</c> when <see cref="Score"/> ≥ the von Heijne acceptance threshold (3.5).</param>
    public readonly record struct SignalPeptide(
        int CleavagePosition,
        double Score,
        string SignalSequence,
        string WindowSequence,
        bool IsLikelySignalPeptide);

    #endregion

    #region Common Motif Patterns (PROSITE-style)

    /// <summary>
    /// Common protein motifs in PROSITE format.
    /// PROSITE patterns verified against https://prosite.expasy.org/ (Hulo et al. 2007).
    /// Non-PROSITE patterns verified against published literature — see inline references.
    /// </summary>
    public static readonly IReadOnlyDictionary<string, PrositePattern> CommonMotifs = new Dictionary<string, PrositePattern>
    {
        // PROSITE PS00001: N-glycosylation site — https://prosite.expasy.org/PS00001
        ["PS00001"] = new("PS00001", "ASN_GLYCOSYLATION", "N-{P}-[ST]-{P}", @"N[^P][ST][^P]", "N-glycosylation site"),

        // PROSITE PS00005: Protein kinase C phosphorylation site — https://prosite.expasy.org/PS00005
        ["PS00005"] = new("PS00005", "PKC_PHOSPHO_SITE", "[ST]-x-[RK]", @"[ST].[RK]", "Protein kinase C phosphorylation site"),

        // PROSITE PS00006: Casein kinase II phosphorylation site — https://prosite.expasy.org/PS00006
        ["PS00006"] = new("PS00006", "CK2_PHOSPHO_SITE", "[ST]-x(2)-[DE]", @"[ST].{2}[DE]", "Casein kinase II phosphorylation site"),

        // PROSITE PS00004: cAMP/cGMP-dependent phosphorylation site — https://prosite.expasy.org/PS00004
        ["PS00004"] = new("PS00004", "CAMP_PHOSPHO_SITE", "[RK](2)-x-[ST]", @"[RK]{2}.[ST]", "cAMP-dependent phosphorylation site"),

        // PROSITE PS00007: Tyrosine kinase phosphorylation site 1 — https://prosite.expasy.org/PS00007
        ["PS00007"] = new("PS00007", "TYR_PHOSPHO_SITE_1", "[RK]-x(2)-[DE]-x(3)-Y", @"[RK].{2}[DE].{3}Y", "Tyrosine kinase phosphorylation site 1"),

        // PROSITE PS00008: N-myristoylation site — https://prosite.expasy.org/PS00008
        ["PS00008"] = new("PS00008", "MYRISTYL", "G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}", @"G[^EDRKHPFYW].{2}[STAGCN][^P]", "N-myristoylation site"),

        // PROSITE PS00009: Amidation site — https://prosite.expasy.org/PS00009
        ["PS00009"] = new("PS00009", "AMIDATION", "x-G-[RK]-[RK]", @".G[RK][RK]", "Amidation site"),

        // PROSITE PS00016: Cell attachment sequence (RGD) — https://prosite.expasy.org/PS00016
        ["PS00016"] = new("PS00016", "RGD", "R-G-D", @"RGD", "Cell attachment sequence"),

        // PROSITE PS00017: ATP/GTP-binding site motif A (P-loop) — https://prosite.expasy.org/PS00017
        ["PS00017"] = new("PS00017", "ATP_GTP_A", "[AG]-x(4)-G-K-[ST]", @"[AG].{4}GK[ST]", "ATP/GTP-binding site motif A (P-loop)"),

        // PROSITE PS00018: EF-hand calcium-binding domain — https://prosite.expasy.org/PS00018
        ["PS00018"] = new("PS00018", "EF_HAND_1", "D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]",
            @"D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]", "EF-hand calcium-binding domain"),

        // PROSITE PS00028: Zinc finger C2H2 type — https://prosite.expasy.org/PS00028
        ["PS00028"] = new("PS00028", "ZINC_FINGER_C2H2_1", "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H",
            @"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H", "Zinc finger C2H2 type"),

        // PROSITE PS00029: Leucine zipper pattern — https://prosite.expasy.org/PS00029
        ["PS00029"] = new("PS00029", "LEUCINE_ZIPPER", "L-x(6)-L-x(6)-L-x(6)-L", @"L.{6}L.{6}L.{6}L", "Leucine zipper pattern"),

        // NLS — Monopartite nuclear localization signal
        // Consensus: K-K/R-X-K/R (Chelsky et al. 1989; Dingwall & Laskey 1991)
        // Position 1 is K only per Chelsky consensus; confirmed by Wikipedia "Nuclear localization sequence"
        ["NLS1"] = new("NLS1", "NLS_MONOPARTITE", "K-[KR]-x-[KR]", @"K[KR].[KR]", "Monopartite nuclear localization signal"),

        // NES — Nuclear export signal (CRM1-dependent)
        // General NES consensus: Φ1-x(2,3)-Φ2-x(2,3)-Φ3-x-Φ4 where Φ = L,I,V,F,M
        // la Cour et al. (2004) Protein Eng Des Sel 17:527-536
        // All Φ positions allow full hydrophobic set per paper definition
        ["NES1"] = new("NES1", "NES", "[LIVFM]-x(2,3)-[LIVFM]-x(2,3)-[LIVFM]-x-[LIVFM]", @"[LIVFM].{2,3}[LIVFM].{2,3}[LIVFM].[LIVFM]", "Nuclear export signal"),

        // SIM — SUMO-interacting motif (type b)
        // Consensus: Ψ-x-Ψ-Ψ where Ψ = V/I/L (Hecker et al. 2006, JBC 281:16117-27)
        ["SIM1"] = new("SIM1", "SUMO_INTERACTION", "[VIL]-x-[VIL]-[VIL]", @"[VIL].[VIL][VIL]", "SUMO interaction motif"),

        // WW domain binding motif — PY motif
        // Consensus: PPxY (Chen & Sudol 1995, PNAS 92:7819-23)
        ["WW1"] = new("WW1", "WW_BINDING", "P-P-x-Y", @"PP.Y", "WW domain binding motif (PY motif)"),

        // SH3 domain binding motif class I
        // Consensus: +xxPxxP where + = R/K (Mayer 2001, J Cell Sci 114:1253-63)
        ["SH3_1"] = new("SH3_1", "SH3_BINDING_1", "[RK]-x(2)-P-x(2)-P", @"[RK].{2}P.{2}P", "SH3 domain binding motif class I"),
    };

    #endregion

    #region Motif Finding

    /// <summary>
    /// Finds all occurrences of common motifs in a protein sequence.
    /// </summary>
    public static IEnumerable<MotifMatch> FindCommonMotifs(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        foreach (var motif in CommonMotifs.Values)
        {
            foreach (var match in FindMotifByPattern(upper, motif.RegexPattern, motif.Name, motif.Accession))
            {
                yield return match;
            }
        }
    }

    /// <summary>
    /// Upper bound on time spent matching a single user-supplied regex pattern. A pathological
    /// (catastrophic-backtracking) pattern that exceeds this budget is treated as a non-matching
    /// failure (no matches) rather than being allowed to hang the scan.
    /// </summary>
    private static readonly TimeSpan RegexMatchTimeout = TimeSpan.FromSeconds(2);

    /// <summary>
    /// Finds all occurrences of a specific pattern in a protein sequence.
    /// Uses lookahead-based matching to discover overlapping occurrences,
    /// consistent with PROSITE ScanProsite behavior (De Castro et al. 2006).
    /// </summary>
    public static IEnumerable<MotifMatch> FindMotifByPattern(
        string proteinSequence,
        string regexPattern,
        string motifName = "Custom",
        string patternId = "")
    {
        if (string.IsNullOrEmpty(proteinSequence) || string.IsNullOrEmpty(regexPattern))
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        Regex regex;
        try
        {
            // Lookahead wrapper enables overlapping match discovery.
            // A bounded match timeout protects against catastrophic-backtracking patterns
            // (e.g. "(A+)+B" over a long homopolymer): a valid-but-pathological regex compiles
            // fine but could otherwise blow up exponentially at match time and hang the scan.
            // RegexMatchTimeoutException is swallowed below — consistent with the documented
            // "invalid/pathological pattern → no matches, never a hang" contract.
            regex = new Regex("(?=(" + regexPattern + "))", RegexOptions.IgnoreCase, RegexMatchTimeout);
        }
        catch
        {
            yield break;
        }

        // Materialize the matches under a guard so a backtracking timeout (thrown during the
        // lazy regex walk) is swallowed rather than leaking out of this enumerator.
        List<Match> matches;
        try
        {
            matches = regex.Matches(upper).Cast<Match>().ToList();
        }
        catch (RegexMatchTimeoutException)
        {
            yield break;
        }

        foreach (Match match in matches)
        {
            var captured = match.Groups[1];

            // A motif occurrence must span at least one residue. A zero-width capture
            // (e.g. a pattern like ".{0}", "$", "(?=...)", or "A?" matching the empty
            // string) would otherwise yield a degenerate MotifMatch with an out-of-bounds
            // Start (== sequence length) and End < Start at the end-of-string position,
            // violating the 0 ≤ Start ≤ End ≤ n−1 coordinate invariant. Skip such captures.
            if (captured.Length == 0)
                continue;

            double score = CalculateMotifScore(captured.Value, regexPattern);

            yield return new MotifMatch(
                Start: match.Index,
                End: match.Index + captured.Length - 1,
                Sequence: captured.Value,
                MotifName: motifName,
                Pattern: patternId,
                Score: score,
                EValue: CalculateEValue(proteinSequence.Length, captured.Length, score));
        }
    }

    /// <summary>
    /// Finds motif using PROSITE pattern syntax.
    /// </summary>
    public static IEnumerable<MotifMatch> FindMotifByProsite(
        string proteinSequence,
        string prositePattern,
        string motifName = "Custom")
    {
        string regexPattern = ConvertPrositeToRegex(prositePattern);
        return FindMotifByPattern(proteinSequence, regexPattern, motifName, prositePattern);
    }

    /// <summary>
    /// Converts PROSITE pattern to regex.
    /// </summary>
    public static string ConvertPrositeToRegex(string prositePattern)
    {
        if (string.IsNullOrEmpty(prositePattern))
            return "";

        var sb = new StringBuilder();
        int i = 0;

        while (i < prositePattern.Length)
        {
            char c = prositePattern[i];

            if (c == '-')
            {
                // Separator, skip
                i++;
            }
            else if (c == 'x')
            {
                // Any amino acid
                if (i + 1 < prositePattern.Length && prositePattern[i + 1] == '(')
                {
                    // x(n) or x(n,m)
                    int end = prositePattern.IndexOf(')', i);
                    if (end > i)
                    {
                        string range = prositePattern.Substring(i + 2, end - i - 2);
                        sb.Append(".{" + range + "}");
                        i = end + 1;
                    }
                    else
                    {
                        sb.Append('.');
                        i++;
                    }
                }
                else
                {
                    sb.Append('.');
                    i++;
                }
            }
            else if (c == '[')
            {
                // Character class - may contain '>' (C-terminus) or '<' (N-terminus)
                // Per PROSITE User Manual §IV.E: "In some rare cases (e.g. PS00267
                // or PS00539), '>' can also occur inside square brackets for the
                // C-terminal element. 'F-[GSTV]-P-R-L-[G>]' means that either
                // 'F-[GSTV]-P-R-L-G' or 'F-[GSTV]-P-R-L>' are considered."
                int end = prositePattern.IndexOf(']', i);
                if (end > i)
                {
                    string content = prositePattern.Substring(i + 1, end - i - 1);
                    bool hasCterm = content.Contains('>');
                    bool hasNterm = content.Contains('<');

                    if (hasCterm || hasNterm)
                    {
                        string letters = content
                            .Replace(">", "")
                            .Replace("<", "");

                        if (hasCterm && letters.Length > 0)
                        {
                            sb.Append("(?:");
                            sb.Append(letters.Length == 1 ? letters : "[" + letters + "]");
                            sb.Append("|$)");
                        }
                        else if (hasNterm && letters.Length > 0)
                        {
                            sb.Append("(?:^|");
                            sb.Append(letters.Length == 1 ? letters : "[" + letters + "]");
                            sb.Append(')');
                        }
                        else if (hasCterm)
                        {
                            sb.Append('$');
                        }
                        else
                        {
                            sb.Append('^');
                        }
                    }
                    else
                    {
                        sb.Append(prositePattern.AsSpan(i, end - i + 1));
                    }
                    i = end + 1;
                }
                else
                {
                    sb.Append(c);
                    i++;
                }
            }
            else if (c == '{')
            {
                // Exclusion class - convert to [^...]
                int end = prositePattern.IndexOf('}', i);
                if (end > i)
                {
                    string excluded = prositePattern.Substring(i + 1, end - i - 1);
                    sb.Append("[^" + excluded + "]");
                    i = end + 1;
                }
                else
                {
                    sb.Append(c);
                    i++;
                }
            }
            else if (c == '<')
            {
                // N-terminus
                sb.Append('^');
                i++;
            }
            else if (c == '>')
            {
                // C-terminus
                sb.Append('$');
                i++;
            }
            else if (c == '(')
            {
                // Repetition count after amino acid
                int end = prositePattern.IndexOf(')', i);
                if (end > i)
                {
                    string range = prositePattern.Substring(i + 1, end - i - 1);
                    sb.Append("{" + range + "}");
                    i = end + 1;
                }
                else
                {
                    i++;
                }
            }
            else if (char.IsLetter(c))
            {
                // Single amino acid
                sb.Append(char.ToUpperInvariant(c));
                i++;
            }
            else if (c == '.')
            {
                // PROSITE User Manual: "A period ends the pattern"
                break;
            }
            else
            {
                // Unsupported construct. The supported PROSITE PA-line grammar is handled by
                // the branches above ('-', 'x', '[...]', '{...}', '<', '>', '(...)', letters,
                // '.'). Any other character — notably the extended ScanProsite *query*
                // metacharacter '*' (Kleene star, e.g. '<{C}*>'), or stray '?'/'+' — would
                // otherwise fall through and be silently dropped, mis-parsing the pattern.
                // Mirror the "reject, don't silently drop" policy used by the Newick parser
                // and throw rather than producing a deceptively-valid regex.
                throw new FormatException(
                    $"Unsupported PROSITE construct: '{c}' at position {i} in pattern " +
                    $"\"{prositePattern}\". The PROSITE→regex converter supports only the " +
                    "standard PA-line grammar (residue letters, x, [...], {...}, <, >, " +
                    "(n)/(m,n) repetition, '-' separators and the '.' terminator); the " +
                    "extended ScanProsite query metacharacter '*' (Kleene star) is not supported.");
            }
        }

        return sb.ToString();
    }

    #endregion

    #region Signal Peptide Prediction

    // --- von Heijne (1986) weight-matrix parameters --------------------------------------------
    // Source: von Heijne G. "A new method for predicting signal sequence cleavage sites."
    //         Nucleic Acids Res. 14:4683-4690 (1986). https://doi.org/10.1093/nar/14.11.4683
    // Reference implementation: EMBOSS 6.6.0 `sigcleave` (emboss/sigcleave.c) with data files
    //         data/Esig.euk (161 eukaryotic signal peptides) and data/Esig.pro (36 prokaryotic).
    //         https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html

    /// <summary>
    /// Number of residue columns in the weight matrix: positions −13..+2 inclusive (von Heijne 1986).
    /// </summary>
    private const int WeightMatrixColumns = 15;

    /// <summary>
    /// Window offset of the first scoring column relative to the candidate +1 residue.
    /// Column 0 = position −13; EMBOSS <c>sigcleave</c> uses <c>pval = −13</c>.
    /// </summary>
    private const int FirstColumnOffset = -13;

    /// <summary>
    /// Matrix column index of position −3 (the conserved small-residue site of the (−3,−1) rule).
    /// EMBOSS treats a zero count here as a strong penalty (<see cref="ZeroCountPenaltyPseudocount"/>).
    /// </summary>
    private const int MinusThreeColumn = 10;

    /// <summary>Matrix column index of position −1 (the residue immediately N-terminal to the cleavage site).</summary>
    private const int MinusOneColumn = 12;

    /// <summary>
    /// Pseudocount substituted for a zero count at columns −3 / −1, producing a large negative
    /// log-odds weight (EMBOSS <c>sigcleave.c</c>: <c>mat[i][j] = 1.0e-10</c> for j==10 || j==12).
    /// </summary>
    private const double ZeroCountPenaltyPseudocount = 1.0e-10;

    /// <summary>
    /// Pseudocount substituted for a zero count at all other columns
    /// (EMBOSS <c>sigcleave.c</c>: <c>mat[i][j] = 1.0</c>).
    /// </summary>
    private const double ZeroCountDefaultPseudocount = 1.0;

    /// <summary>
    /// Default minimum acceptance weight. EMBOSS recommends <c>-minweight ≥ 3.5</c>: at this level the
    /// method identifies ≈95% of signal peptides and rejects ≈95% of non-signal peptides
    /// (von Heijne 1986; EMBOSS <c>sigcleave</c> documentation).
    /// </summary>
    public const double DefaultMinWeight = 3.5;

    /// <summary>
    /// Eukaryotic amino-acid count matrix (161 aligned signal peptides), rows in alphabetical
    /// one-letter order A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y; columns −13..+2 then the
    /// per-residue background "Expect" count. Verbatim from EMBOSS 6.6.0 data/Esig.euk
    /// (von Heijne 1986).
    /// </summary>
    private static readonly int[][] EukaryoticCounts =
    {
        //         -13 -12 -11 -10  -9  -8  -7  -6  -5  -4  -3  -2  -1  +1  +2   (Expect is held separately)
        /* A */ new[] { 16, 13, 14, 15, 20, 18, 18, 17, 25, 15, 47,  6, 80, 18,  6 },
        /* C */ new[] {  3,  6,  9,  7,  9, 14,  6,  8,  5,  6, 19,  3,  9,  8,  3 },
        /* D */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  5,  3,  0,  5,  0, 10, 11 },
        /* E */ new[] {  0,  0,  0,  1,  0,  0,  0,  0,  3,  7,  0,  7,  0, 13, 14 },
        /* F */ new[] { 13,  9, 11, 11,  6,  7, 18, 13,  4,  5,  0, 13,  0,  6,  4 },
        /* G */ new[] {  4,  4,  3,  6,  3, 13,  3,  2, 19, 34,  5,  7, 39, 10,  7 },
        /* H */ new[] {  0,  0,  0,  0,  0,  1,  1,  0,  5,  0,  0,  6,  0,  4,  2 },
        /* I */ new[] { 15, 15,  8,  6, 11,  5,  4,  8,  5,  1, 10,  5,  0,  8,  7 },
        /* K */ new[] {  0,  0,  0,  1,  0,  0,  1,  0,  0,  4,  0,  2,  0, 11,  9 },
        /* L */ new[] { 71, 68, 72, 79, 78, 45, 64, 49, 10, 23,  8, 20,  1,  8,  4 },
        /* M */ new[] {  0,  3,  7,  4,  1,  6,  2,  2,  0,  0,  0,  1,  0,  1,  2 },
        /* N */ new[] {  0,  1,  0,  1,  1,  0,  0,  0,  3,  3,  0, 10,  0,  4,  7 },
        /* P */ new[] {  2,  0,  2,  0,  0,  4,  1,  8, 20, 14,  0,  1,  3,  0, 22 },
        /* Q */ new[] {  0,  0,  0,  1,  0,  6,  1,  0, 10,  8,  0, 18,  3, 19, 10 },
        /* R */ new[] {  2,  0,  0,  0,  0,  1,  0,  0,  7,  4,  0, 15,  0, 12,  9 },
        /* S */ new[] {  9,  3,  8,  6, 13, 10, 15, 16, 26, 11, 23, 17, 20, 15, 10 },
        /* T */ new[] {  2, 10,  5,  4,  5, 13,  7,  7, 12,  6, 17,  8,  6,  3, 10 },
        /* V */ new[] { 20, 25, 15, 18, 13, 15, 11, 27,  0, 12, 32,  3,  0,  8, 17 },
        /* W */ new[] {  4,  3,  3,  1,  1,  2,  6,  3,  1,  3,  0,  9,  0,  2,  0 },
        /* Y */ new[] {  0,  1,  4,  0,  0,  1,  3,  1,  1,  2,  0,  5,  0,  1,  7 },
    };

    /// <summary>Eukaryotic background "Expect" counts (one per residue), from EMBOSS data/Esig.euk.</summary>
    private static readonly double[] EukaryoticExpect =
    {
        /* A */ 14.5, /* C */ 4.5, /* D */ 8.9, /* E */ 10.0, /* F */ 5.6,
        /* G */ 12.1, /* H */ 3.4, /* I */ 7.4, /* K */ 11.3, /* L */ 12.1,
        /* M */ 2.7,  /* N */ 7.1, /* P */ 7.4, /* Q */ 6.3,  /* R */ 7.6,
        /* S */ 11.4, /* T */ 9.7, /* V */ 11.1,/* W */ 1.8,  /* Y */ 5.6,
    };

    /// <summary>
    /// Prokaryotic amino-acid count matrix (36 aligned signal peptides); same row/column layout
    /// as <see cref="EukaryoticCounts"/>. Verbatim from EMBOSS 6.6.0 data/Esig.pro (von Heijne 1986).
    /// </summary>
    private static readonly int[][] ProkaryoticCounts =
    {
        /* A */ new[] { 10,  8,  8,  9,  6,  7,  5,  6,  7,  7, 24,  2, 31, 18,  4 },
        /* C */ new[] {  1,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0 },
        /* D */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  8 },
        /* E */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  4,  8 },
        /* F */ new[] {  2,  4,  3,  4,  1,  1,  8,  0,  4,  1,  0,  7,  0,  1,  0 },
        /* G */ new[] {  4,  2,  2,  2,  3,  5,  2,  4,  2,  2,  0,  2,  2,  1,  0 },
        /* H */ new[] {  0,  0,  1,  0,  0,  0,  0,  1,  1,  0,  0,  7,  0,  1,  0 },
        /* I */ new[] {  3,  1,  5,  1,  5,  0,  1,  3,  0,  0,  0,  0,  0,  0,  2 },
        /* K */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  3,  0 },
        /* L */ new[] {  8, 11,  9,  8,  9, 13,  1,  0,  2,  2,  1,  2,  0,  0,  1 },
        /* M */ new[] {  0,  2,  1,  1,  3,  2,  3,  0,  1,  2,  0,  4,  0,  0,  1 },
        /* N */ new[] {  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  3,  0,  1,  4 },
        /* P */ new[] {  0,  1,  1,  1,  1,  1,  2,  3,  5,  2,  0,  0,  0,  0,  5 },
        /* Q */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  3,  0,  0,  1 },
        /* R */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
        /* S */ new[] {  1,  0,  1,  4,  4,  1,  5, 15,  5,  8,  5,  2,  2,  0,  0 },
        /* T */ new[] {  2,  0,  4,  2,  2,  2,  2,  2,  5,  1,  3,  0,  1,  1,  2 },
        /* V */ new[] {  5,  7,  1,  3,  1,  4,  7,  0,  0,  4,  3,  0,  0,  2,  0 },
        /* W */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
        /* Y */ new[] {  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0 },
    };

    /// <summary>Prokaryotic background "Expect" counts, from EMBOSS data/Esig.pro.</summary>
    private static readonly double[] ProkaryoticExpect =
    {
        /* A */ 3.2, /* C */ 1.0, /* D */ 2.0, /* E */ 2.2, /* F */ 1.3,
        /* G */ 2.7, /* H */ 0.8, /* I */ 1.7, /* K */ 2.5, /* L */ 2.7,
        /* M */ 0.6, /* N */ 1.6, /* P */ 1.7, /* Q */ 1.4, /* R */ 1.7,
        /* S */ 2.6, /* T */ 2.2, /* V */ 2.5, /* W */ 0.4, /* Y */ 1.3,
    };

    /// <summary>Residue rows of the weight matrices, in matrix order (index = matrix row).</summary>
    private const string MatrixResidues = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// Log-odds weight matrix [residueRow][column] = ln(count/expect), built once from
    /// <see cref="EukaryoticCounts"/> with EMBOSS pseudocount handling.
    /// </summary>
    private static readonly double[][] EukaryoticWeights = BuildWeightMatrix(EukaryoticCounts, EukaryoticExpect);

    /// <summary>Log-odds weight matrix built from <see cref="ProkaryoticCounts"/>.</summary>
    private static readonly double[][] ProkaryoticWeights = BuildWeightMatrix(ProkaryoticCounts, ProkaryoticExpect);

    /// <summary>
    /// Predicts the signal-peptide cleavage site using von Heijne's (1986) position-specific
    /// weight-matrix method, scoring identically to the EMBOSS <c>sigcleave</c> reference
    /// implementation: at each candidate site the score is the sum of log-odds weights
    /// (ln of observed/expected residue frequency) over positions −13..+2, and the highest-scoring
    /// site is returned.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence (one-letter codes; case-insensitive).</param>
    /// <param name="prokaryote">
    /// When <c>true</c>, scores against the 36-sequence prokaryotic matrix; otherwise the default
    /// 161-sequence eukaryotic matrix. Source: EMBOSS <c>sigcleave -prokaryote</c>.
    /// </param>
    /// <param name="minWeight">
    /// Acceptance threshold for <see cref="SignalPeptide.IsLikelySignalPeptide"/>; defaults to
    /// <see cref="DefaultMinWeight"/> (3.5).
    /// </param>
    /// <returns>
    /// The best-scoring cleavage site, or <c>null</c> when the sequence is null/empty or too short
    /// to contain a full scoring window (fewer than <see cref="WeightMatrixColumns"/> residues).
    /// </returns>
    public static SignalPeptide? PredictSignalPeptide(
        string proteinSequence,
        bool prokaryote = false,
        double minWeight = DefaultMinWeight)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < WeightMatrixColumns)
            return null;

        string upper = proteinSequence.ToUpperInvariant();
        double[][] weights = prokaryote ? ProkaryoticWeights : EukaryoticWeights;

        double bestScore = double.NegativeInfinity;
        int bestSite = -1;

        // Candidate site i (0-based) is the +1 residue (first residue of the mature protein).
        // The scoring window spans matrix columns 0..14 → sequence positions i+FirstColumnOffset..i+1.
        // EMBOSS sigcleave.c: weight = Σ matrix[aa(j)][ic], j = i−13..i+1, ic = 0..14.
        for (int i = 0; i < upper.Length; i++)
        {
            double weight = 0.0;
            for (int col = 0; col < WeightMatrixColumns; col++)
            {
                int j = i + FirstColumnOffset + col;
                if (j < 0 || j >= upper.Length)
                    continue;

                int row = MatrixResidues.IndexOf(upper[j]);
                if (row < 0)
                    continue; // non-standard residue (X, B, Z, *) contributes nothing, as in EMBOSS

                weight += weights[row][col];
            }

            if (weight > bestScore)
            {
                bestScore = weight;
                bestSite = i;
            }
        }

        if (bestSite < 0)
            return null;

        // Cleavage is between residue (bestSite − 1) and bestSite; mature protein starts at bestSite.
        int cleavagePosition = bestSite + 1; // 1-based mature start
        int windowStart = Math.Max(0, bestSite + FirstColumnOffset);
        int windowEnd = Math.Min(upper.Length - 1, bestSite + 1); // inclusive (position +2)

        return new SignalPeptide(
            CleavagePosition: cleavagePosition,
            Score: bestScore,
            SignalSequence: upper.Substring(0, bestSite),
            WindowSequence: upper.Substring(windowStart, windowEnd - windowStart + 1),
            IsLikelySignalPeptide: bestScore >= minWeight);
    }

    /// <summary>
    /// Builds the log-odds weight matrix from raw counts using the EMBOSS <c>sigcleave</c> transform:
    /// a zero count becomes <see cref="ZeroCountPenaltyPseudocount"/> at columns −3/−1 (a strong
    /// penalty) and <see cref="ZeroCountDefaultPseudocount"/> elsewhere, then each cell is
    /// <c>ln(count / expect)</c>. Source: EMBOSS 6.6.0 emboss/sigcleave.c <c>sigcleave_readSig</c>.
    /// </summary>
    private static double[][] BuildWeightMatrix(int[][] counts, double[] expect)
    {
        var weights = new double[counts.Length][];
        for (int row = 0; row < counts.Length; row++)
        {
            weights[row] = new double[WeightMatrixColumns];
            for (int col = 0; col < WeightMatrixColumns; col++)
            {
                double value = counts[row][col];
                if (value == 0.0)
                {
                    value = (col == MinusThreeColumn || col == MinusOneColumn)
                        ? ZeroCountPenaltyPseudocount
                        : ZeroCountDefaultPseudocount;
                }

                weights[row][col] = Math.Log(value / expect[row]);
            }
        }

        return weights;
    }

    #endregion

    #region Transmembrane Prediction

    // --- Kyte & Doolittle (1982) hydropathy parameters -----------------------------------------
    // Source: Kyte J, Doolittle RF. "A simple method for displaying the hydropathic character of a
    //         protein." J Mol Biol. 157(1):105-132 (1982). https://doi.org/10.1016/0022-2836(82)90515-0
    // TM-detection parameters (window 19, threshold 1.6) verified against the Davidson College
    //         Kyte-Doolittle background page (citing Kyte & Doolittle 1982):
    //         https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
    // Windowing (arithmetic mean of the window's per-residue values) matches Biopython 1.x
    //         Bio.SeqUtils.ProtParam.protein_scale with edge weight 1.0:
    //         https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py

    /// <summary>
    /// Default sliding-window width for transmembrane-helix detection. Kyte &amp; Doolittle (1982)
    /// found a 19-residue window optimal for identifying membrane-spanning segments.
    /// </summary>
    private const int DefaultTransmembraneWindow = 19;

    /// <summary>
    /// Default hydropathy threshold: a window whose mean Kyte-Doolittle score exceeds 1.6 is
    /// considered part of a transmembrane segment (Kyte &amp; Doolittle 1982; Davidson background page).
    /// </summary>
    private const double DefaultTransmembraneThreshold = 1.6;

    /// <summary>
    /// Minimum span (residues) a candidate region must cover to be reported as a transmembrane helix.
    /// A single α-helix needs ≈18–21 residues to cross the ≈30 Å lipid bilayer, matching the
    /// 19-residue scanning window (each residue rises ≈1.5 Å along the helix axis). Set to the
    /// window width so that any region containing at least one above-threshold window qualifies.
    /// </summary>
    private const int MinTransmembraneHelixLength = DefaultTransmembraneWindow;

    /// <summary>
    /// Predicts transmembrane α-helices with the Kyte &amp; Doolittle (1982) hydropathy method:
    /// the mean hydropathy is computed over a sliding window of <paramref name="windowSize"/>
    /// residues, and each maximal run of windows whose mean exceeds <paramref name="threshold"/>
    /// is reported as a candidate transmembrane segment.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence (one-letter codes; case-insensitive).</param>
    /// <param name="windowSize">Sliding-window width; defaults to 19 (Kyte &amp; Doolittle 1982).</param>
    /// <param name="threshold">Mean-hydropathy cutoff; defaults to 1.6 (Kyte &amp; Doolittle 1982).</param>
    /// <returns>
    /// One tuple per predicted segment: 0-based inclusive <c>Start</c>/<c>End</c> residue indices and
    /// the peak (maximum) window mean within the segment. Empty when the sequence is null/empty or
    /// shorter than <paramref name="windowSize"/>.
    /// </returns>
    public static IEnumerable<(int Start, int End, double Score)> PredictTransmembraneHelices(
        string proteinSequence,
        int windowSize = DefaultTransmembraneWindow,
        double threshold = DefaultTransmembraneThreshold)
    {
        if (string.IsNullOrEmpty(proteinSequence) || windowSize <= 0 || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        // Calculate hydropathy profile
        var hydropathy = CalculateHydropathyProfile(upper, windowSize);

        // Find regions above threshold
        int? regionStart = null;
        double maxScore = 0;

        for (int i = 0; i < hydropathy.Count; i++)
        {
            if (hydropathy[i] >= threshold)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    maxScore = hydropathy[i];
                }
                else
                {
                    maxScore = Math.Max(maxScore, hydropathy[i]);
                }
            }
            else if (regionStart != null)
            {
                // End of region: the last above-threshold window starts at (i - 1) and therefore
                // covers residues up to (i - 1) + windowSize - 1. The reported End is that last
                // covered residue (0-based inclusive), clamped to the last residue index.
                int start = regionStart.Value;
                int lastCoveredResidue = i - 1 + windowSize - 1;

                if (lastCoveredResidue - start + 1 >= MinTransmembraneHelixLength)
                {
                    yield return (start, Math.Min(lastCoveredResidue, upper.Length - 1), maxScore);
                }

                regionStart = null;
                maxScore = 0;
            }
        }

        // Handle last region
        if (regionStart != null)
        {
            int start = regionStart.Value;
            int lastCoveredResidue = hydropathy.Count - 1 + windowSize - 1;

            if (lastCoveredResidue - start + 1 >= MinTransmembraneHelixLength)
            {
                yield return (start, Math.Min(lastCoveredResidue, upper.Length - 1), maxScore);
            }
        }
    }

    /// <summary>
    /// Kyte-Doolittle hydropathy scale (Kyte &amp; Doolittle 1982, J Mol Biol 157:105-132).
    /// Values verified against QIAGEN CLC Genomics Workbench "Hydrophobicity scales" reference.
    /// </summary>
    private static readonly Dictionary<char, double> HydropathyScale = new()
    {
        {'A', 1.8}, {'R', -4.5}, {'N', -3.5}, {'D', -3.5}, {'C', 2.5},
        {'Q', -3.5}, {'E', -3.5}, {'G', -0.4}, {'H', -3.2}, {'I', 4.5},
        {'L', 3.8}, {'K', -3.9}, {'M', 1.9}, {'F', 2.8}, {'P', -1.6},
        {'S', -0.8}, {'T', -0.7}, {'W', -0.9}, {'Y', -1.3}, {'V', 4.2}
    };

    /// <summary>
    /// Computes the Kyte-Doolittle sliding-window hydropathy profile as the arithmetic mean of the
    /// per-residue scores in each window (Biopython <c>protein_scale</c> with edge weight 1.0).
    /// Non-standard residues (X, B, Z, *) carry no scale value and are excluded from the mean.
    /// </summary>
    private static List<double> CalculateHydropathyProfile(string sequence, int windowSize)
    {
        var profile = new List<double>();

        for (int i = 0; i <= sequence.Length - windowSize; i++)
        {
            double sum = 0;
            int count = 0;

            for (int j = i; j < i + windowSize; j++)
            {
                if (HydropathyScale.TryGetValue(sequence[j], out double value))
                {
                    sum += value;
                    count++;
                }
            }

            profile.Add(count > 0 ? sum / count : 0);
        }

        return profile;
    }

    #endregion

    #region Disorder Prediction

    /// <summary>
    /// Predicts intrinsically disordered regions.
    /// </summary>
    public static IEnumerable<(int Start, int End, double Score)> PredictDisorderedRegions(
        string proteinSequence,
        int windowSize = 21,
        double threshold = 0.5)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        // Disorder-promoting amino acids
        var disorderPropensity = new Dictionary<char, double>
        {
            {'A', 0.06}, {'R', 0.18}, {'N', 0.23}, {'D', 0.19}, {'C', -0.02},
            {'Q', 0.23}, {'E', 0.24}, {'G', 0.17}, {'H', 0.10}, {'I', -0.49},
            {'L', -0.34}, {'K', 0.26}, {'M', -0.23}, {'F', -0.42}, {'P', 0.41},
            {'S', 0.14}, {'T', 0.04}, {'W', -0.49}, {'Y', -0.31}, {'V', -0.39}
        };

        var profile = new List<double>();

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double sum = 0;
            int count = 0;

            for (int j = i; j < i + windowSize; j++)
            {
                if (disorderPropensity.TryGetValue(upper[j], out double value))
                {
                    sum += value;
                    count++;
                }
            }

            double avgPropensity = count > 0 ? sum / count : 0;
            // Normalize to 0-1 range
            double normalized = (avgPropensity + 0.5) / 1.0;
            profile.Add(Math.Max(0, Math.Min(1, normalized)));
        }

        // Find regions above threshold
        int? regionStart = null;
        double maxScore = 0;

        for (int i = 0; i < profile.Count; i++)
        {
            if (profile[i] >= threshold)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    maxScore = profile[i];
                }
                else
                {
                    maxScore = Math.Max(maxScore, profile[i]);
                }
            }
            else if (regionStart != null)
            {
                int start = regionStart.Value;
                int end = i - 1 + windowSize / 2;

                if (end - start >= 10)
                {
                    yield return (start, end, maxScore);
                }

                regionStart = null;
                maxScore = 0;
            }
        }

        if (regionStart != null)
        {
            int start = regionStart.Value;
            int end = profile.Count - 1 + windowSize / 2;

            if (end - start >= 10)
            {
                yield return (start, Math.Min(end, upper.Length - 1), maxScore);
            }
        }
    }

    #endregion

    #region Coiled-Coil Prediction

    // --- Coiled-coil heptad-repeat constants (all source-traceable) ---

    // Heptad repeat is denoted (abcdefg)n; the hydrophobic core positions are a and d.
    // Lupas A, Van Dyke M, Stock J (1991) Science 252:1162-1164; Mason JM, Arndt KM (2004)
    // ChemBioChem 5(2):170-176.
    private const int HeptadLength = 7;          // a,b,c,d,e,f,g
    private const int HeptadPositionA = 0;        // core position a
    private const int HeptadPositionD = 3;        // core position d

    // Hydrophobic-core residues at positions a and d: "isoleucine, leucine, or valine"
    // (Wikipedia "Coiled coil", citing Mason & Arndt 2004).
    private static readonly char[] CoiledCoilCoreResidues = { 'I', 'L', 'V' };

    // Default scoring window = 28 residues (4 heptads) and 7 candidate registers per Lupas (1991):
    // "a window can be assigned seven different heptad repeat frames ... 28 positions in a gliding window".
    private const int DefaultCoiledCoilWindow = 28;

    // Default threshold = 0.5: a window is coiled-coil if MORE THAN HALF (>= 50%) of its a/d positions
    // are occupied by hydrophobic-core residues ("predominantly hydrophobic at a and d", Mason & Arndt 2004).
    private const double DefaultCoiledCoilThreshold = 0.5;

    // Minimum reported region length = 21 residues (3 heptads): naturally occurring coiled coils are
    // built from multiple heptads, (abcdefg)1-(abcdefg)2-(abcdefg)3... (Mason & Arndt 2004).
    private const int MinCoiledCoilRegion = 3 * HeptadLength;

    /// <summary>
    /// Predicts coiled-coil regions by heptad-repeat analysis. For every residue window the score is the
    /// fraction of its heptad a/d core positions occupied by a hydrophobic core residue (I, L or V),
    /// maximised over the seven possible heptad registers (frames); contiguous windows scoring at or above
    /// <paramref name="threshold"/> form a region. This is the fully-specified heptad a/d occupancy model
    /// (Lupas 1991; Mason &amp; Arndt 2004); it deliberately does NOT use the COILS position-specific
    /// scoring matrix, whose weights were not available from authoritative sources.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence (case-insensitive).</param>
    /// <param name="windowSize">Sliding-window length in residues (default 28 = 4 heptads).</param>
    /// <param name="threshold">Minimum a/d hydrophobic-core occupancy fraction in [0,1] (default 0.5).</param>
    /// <returns>Non-overlapping regions (Start, End inclusive 0-based; Score = peak occupancy fraction).</returns>
    public static IEnumerable<(int Start, int End, double Score)> PredictCoiledCoils(
        string proteinSequence,
        int windowSize = DefaultCoiledCoilWindow,
        double threshold = DefaultCoiledCoilThreshold)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();
        var coreResidues = new HashSet<char>(CoiledCoilCoreResidues);

        int windowCount = upper.Length - windowSize + 1;
        var profile = new double[windowCount];
        for (int i = 0; i < windowCount; i++)
        {
            profile[i] = BestHeptadOccupancy(upper, i, windowSize, coreResidues);
        }

        // Group contiguous above-threshold windows into regions; a window at profile index i covers
        // residues [i, i + windowSize - 1]. Region score = peak window occupancy in the run.
        int? regionStart = null;
        double peak = 0;

        for (int i = 0; i < profile.Length; i++)
        {
            if (profile[i] >= threshold)
            {
                if (regionStart is null)
                {
                    regionStart = i;
                    peak = profile[i];
                }
                else
                {
                    peak = Math.Max(peak, profile[i]);
                }
            }
            else if (regionStart is not null)
            {
                var region = BuildRegion(regionStart.Value, i - 1, windowSize, peak);
                if (region.HasValue)
                    yield return region.Value;
                regionStart = null;
                peak = 0;
            }
        }

        if (regionStart is not null)
        {
            var region = BuildRegion(regionStart.Value, profile.Length - 1, windowSize, peak);
            if (region.HasValue)
                yield return region.Value;
        }
    }

    /// <summary>
    /// Returns the maximum, over the seven heptad registers, of the fraction of a/d core positions in the
    /// window starting at <paramref name="windowStart"/> that are occupied by a hydrophobic core residue.
    /// </summary>
    private static double BestHeptadOccupancy(
        string upper, int windowStart, int windowSize, HashSet<char> coreResidues)
    {
        double best = 0;
        for (int register = 0; register < HeptadLength; register++)
        {
            int coreCount = 0;
            int hydrophobicCount = 0;
            for (int j = 0; j < windowSize; j++)
            {
                int index = windowStart + j;
                int heptadPosition = Mod(index - register, HeptadLength);
                if (heptadPosition != HeptadPositionA && heptadPosition != HeptadPositionD)
                    continue;

                coreCount++;
                if (coreResidues.Contains(upper[index]))
                    hydrophobicCount++;
            }

            if (coreCount == 0)
                continue;

            double occupancy = (double)hydrophobicCount / coreCount;
            if (occupancy > best)
                best = occupancy;
        }

        return best;
    }

    /// <summary>
    /// Maps a run of above-threshold window indices [<paramref name="firstWindow"/>,
    /// <paramref name="lastWindow"/>] to a residue region, keeping it only if it spans at least
    /// <see cref="MinCoiledCoilRegion"/> residues.
    /// </summary>
    private static (int Start, int End, double Score)? BuildRegion(
        int firstWindow, int lastWindow, int windowSize, double peak)
    {
        int start = firstWindow;
        int end = lastWindow + windowSize - 1;
        if (end - start + 1 < MinCoiledCoilRegion)
            return null;
        return (start, end, peak);
    }

    /// <summary>Non-negative modulo (C# % can be negative for negative operands).</summary>
    private static int Mod(int value, int modulus)
    {
        int result = value % modulus;
        return result < 0 ? result + modulus : result;
    }

    #endregion

    #region Low Complexity Regions

    /// <summary>
    /// Default SEG sliding-window length W. Source: NCBI SEG man page ("Trigger window length
    /// [Default 12]") and NCBI blast_seg.c (<c>kSegWindow = 12</c>).
    /// </summary>
    private const int DefaultSegWindow = 12;

    /// <summary>
    /// Default SEG trigger complexity K1 (locut), in bits/residue. A window with complexity ≤ K1
    /// triggers a low-complexity segment. Source: NCBI SEG man page ("Trigger complexity [Default 2.2]")
    /// and blast_seg.c (<c>kSegLocut = 2.2</c>).
    /// </summary>
    private const double DefaultSegTriggerComplexity = 2.2;

    /// <summary>
    /// Default SEG extension complexity K2 (hicut), in bits/residue. A triggered segment is extended
    /// over neighbouring windows whose complexity is ≤ K2. Source: NCBI SEG man page
    /// ("Extension complexity [Default 2.5]") and blast_seg.c (<c>kSegHicut = 2.5</c>).
    /// </summary>
    private const double DefaultSegExtensionComplexity = 2.5;

    /// <summary>
    /// Finds low-complexity regions in a protein sequence using the SEG algorithm of
    /// Wootton &amp; Federhen (1993). Local complexity is the Shannon entropy of the residue
    /// composition in a sliding window, measured in bits per residue:
    /// K = −Σ pᵢ·log₂(pᵢ), with maximum log₂(20) ≈ 4.322 for the 20 amino acids.
    /// A window with complexity ≤ <paramref name="triggerComplexity"/> (K1) triggers a
    /// low-complexity segment, which is then extended over adjacent windows whose complexity
    /// is ≤ <paramref name="extensionComplexity"/> (K2).
    /// References: Wootton &amp; Federhen (1993) Comput. Chem. 17:149–163, eq. (3);
    /// NCBI SEG (ncbi-seg) man page; NCBI blast_seg.c (s_Entropy).
    /// </summary>
    /// <param name="proteinSequence">Protein sequence (case-insensitive).</param>
    /// <param name="windowSize">Sliding-window length W (default 12).</param>
    /// <param name="triggerComplexity">Trigger complexity K1 in bits/residue (default 2.2).</param>
    /// <param name="extensionComplexity">Extension complexity K2 in bits/residue (default 2.5).</param>
    /// <returns>
    /// 0-based, inclusive low-complexity regions as <c>(Start, End, Complexity)</c>, where
    /// <c>Complexity</c> is the minimum window complexity (bits/residue) observed inside the region.
    /// Empty if the sequence is null/empty or shorter than the window.
    /// </returns>
    public static IEnumerable<(int Start, int End, double Complexity)> FindLowComplexityRegions(
        string proteinSequence,
        int windowSize = DefaultSegWindow,
        double triggerComplexity = DefaultSegTriggerComplexity,
        double extensionComplexity = DefaultSegExtensionComplexity)
    {
        if (windowSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(windowSize), "Window size must be positive.");

        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();
        int windowCount = upper.Length - windowSize + 1;

        // Per-window SEG complexity (bits/residue) for every window start position.
        var windowComplexity = new double[windowCount];
        for (int i = 0; i < windowCount; i++)
            windowComplexity[i] = CalculateSegComplexity(upper.AsSpan(i, windowSize));

        // SEG pass 1: find trigger windows (complexity ≤ K1); pass 2: extend over adjacent
        // windows with complexity ≤ K2. Windows form a contiguous run; the residue span of the
        // run is [firstStart, lastStart + W - 1] (0-based inclusive).
        int i2 = 0;
        while (i2 < windowCount)
        {
            if (windowComplexity[i2] > extensionComplexity)
            {
                i2++;
                continue;
            }

            // Start of a maximal run of windows with complexity ≤ K2.
            int runStart = i2;
            double minComplexity = windowComplexity[i2];
            bool triggered = windowComplexity[i2] <= triggerComplexity;
            int runEnd = i2;
            while (runEnd + 1 < windowCount && windowComplexity[runEnd + 1] <= extensionComplexity)
            {
                runEnd++;
                if (windowComplexity[runEnd] < minComplexity)
                    minComplexity = windowComplexity[runEnd];
                if (windowComplexity[runEnd] <= triggerComplexity)
                    triggered = true;
            }

            // Only emit the run if at least one window actually triggered (complexity ≤ K1).
            if (triggered)
                yield return (runStart, runEnd + windowSize - 1, minComplexity);

            i2 = runEnd + 1;
        }
    }

    /// <summary>
    /// Computes the SEG local complexity of a residue window as Shannon entropy in bits/residue:
    /// K = −Σ pᵢ·log₂(pᵢ), where pᵢ is the fraction of the window occupied by residue type i.
    /// A homopolymer window yields 0; a window of W distinct residues yields log₂(W).
    /// Source: Wootton &amp; Federhen (1993) eq. (3); NCBI blast_seg.c s_Entropy; SeqComplex `ce`.
    /// </summary>
    private static double CalculateSegComplexity(ReadOnlySpan<char> window)
    {
        Span<int> counts = stackalloc int[char.MaxValue + 1];
        foreach (char c in window)
            counts[c]++;

        double length = window.Length;
        double entropy = 0.0;
        foreach (char c in window)
        {
            // Visit each distinct residue once: tally on the first occurrence only.
            if (counts[c] == 0)
                continue;
            double p = counts[c] / length;
            entropy -= p * Math.Log2(p);
            counts[c] = 0; // mark consumed so duplicates are not double-counted
        }

        return entropy;
    }

    #endregion

    #region Scoring Functions

    /// <summary>
    /// Number of standard amino acids used as the alphabet size for scoring.
    /// </summary>
    private const int AminoAcidAlphabetSize = 20;

    /// <summary>
    /// Calculates motif score as total information content (IC) in bits.
    /// IC per position = log2(alphabet_size / allowed_count).
    /// Reference: Schneider TD, Stephens RM (1990). "Sequence logos: a new way to
    /// display consensus sequences." Nucleic Acids Res 18:6097-6100.
    /// </summary>
    private static double CalculateMotifScore(string matchSequence, string pattern)
    {
        var allowedCounts = ParseRegexAllowedCounts(pattern);

        double totalIC = 0;
        foreach (int allowed in allowedCounts)
        {
            if (allowed > 0 && allowed < AminoAcidAlphabetSize)
                totalIC += Math.Log2((double)AminoAcidAlphabetSize / allowed);
        }

        return totalIC;
    }

    /// <summary>
    /// Calculates E-value: expected number of random matches in a sequence.
    /// E = (N - L + 1) × P_random, where P_random = 2^(-IC_total) = ∏(k_i / 20).
    /// This is a direct combinatorial probability under uniform amino acid distribution;
    /// IC definition per Schneider & Stephens (1990) Nucleic Acids Res 18:6097-6100.
    /// </summary>
    private static double CalculateEValue(int sequenceLength, int motifLength, double score)
    {
        int searchSpace = Math.Max(1, sequenceLength - motifLength + 1);
        double probability = Math.Pow(2.0, -score);
        return searchSpace * probability;
    }

    /// <summary>
    /// Parses a regex pattern to determine the number of allowed amino acids per position.
    /// Used by CalculateMotifScore for information-content-based scoring.
    /// Handles: single letters (1), character classes [ABC] (N), negated classes [^ABC] (20-N),
    /// any position '.' (20), and quantifiers {n} / {n,m} (minimum repeat count).
    /// </summary>
    private static List<int> ParseRegexAllowedCounts(string regexPattern)
    {
        var counts = new List<int>();
        int i = 0;

        while (i < regexPattern.Length)
        {
            int allowed;

            if (regexPattern[i] == '[')
            {
                i++; // skip [
                bool negated = i < regexPattern.Length && regexPattern[i] == '^';
                if (negated) i++; // skip ^

                int letterCount = 0;
                while (i < regexPattern.Length && regexPattern[i] != ']')
                {
                    letterCount++;
                    i++;
                }
                if (i < regexPattern.Length) i++; // skip ]

                allowed = negated
                    ? AminoAcidAlphabetSize - letterCount
                    : letterCount;
                allowed = Math.Clamp(allowed, 1, AminoAcidAlphabetSize);
            }
            else if (regexPattern[i] == '.')
            {
                i++;
                allowed = AminoAcidAlphabetSize;
            }
            else if (char.IsLetter(regexPattern[i]))
            {
                i++;
                allowed = 1;
            }
            else
            {
                // Skip non-pattern characters (e.g., parentheses from lookahead wrappers)
                i++;
                continue;
            }

            // Check for quantifier: {n} or {n,m}
            int repeat = 1;
            if (i < regexPattern.Length && regexPattern[i] == '{')
            {
                i++; // skip {
                var num = new StringBuilder();
                while (i < regexPattern.Length && regexPattern[i] != ',' && regexPattern[i] != '}')
                {
                    num.Append(regexPattern[i]);
                    i++;
                }
                if (int.TryParse(num.ToString(), out int n))
                    repeat = n;

                // Skip to closing }
                while (i < regexPattern.Length && regexPattern[i] != '}')
                    i++;
                if (i < regexPattern.Length) i++; // skip }
            }

            for (int j = 0; j < repeat; j++)
                counts.Add(allowed);
        }

        return counts;
    }

    #endregion

    #region Domain Finding

    // --- Exact PROSITE PATTERN signatures (deterministic, regex-like, citable) ---
    //
    // The following three signatures are reproduced VERBATIM from official PROSITE PATTERN
    // entries and translated to .NET regex using the PROSITE pattern syntax rules
    // (PROSITE/ScanProsite user manual: 'x' → any residue, [..] → allowed set,
    // {..} → excluded set, '-' separators dropped, x(n) → repetition):
    // https://prosite.expasy.org/scanprosite/scanprosite_doc.html
    //
    // These are deterministic patterns, NOT trained HMM profiles, so they are reproduced exactly.

    /// <summary>
    /// Zinc finger C2H2 PROSITE PATTERN PS00028 (Pfam PF00096): C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.
    /// Source: PROSITE PS00028 (https://prosite.expasy.org/PS00028); Krishna SS et al. (2003) NAR 31:532–550.
    /// </summary>
    private const string ZincFingerC2H2Pattern = @"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";

    /// <summary>
    /// WD-repeats signature PROSITE PATTERN PS00678 (WD_REPEATS_1; Pfam PF00400), verbatim:
    /// [LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN].
    /// Translated: x → '.', x(2) → '.{2}', {P} → '[^P]', {DP} → '[^DP]', '-' separators dropped.
    /// Source: PROSITE PS00678 (https://prosite.expasy.org/PS00678); Neer EJ et al. (1994) Nature 371:297–300.
    /// </summary>
    private const string Wd40RepeatPattern =
        @"[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]";

    /// <summary>
    /// Walker A / P-loop (ATP/GTP-binding) PROSITE PATTERN PS00017 (Pfam PF00069 kinase ATP-binding site):
    /// [AG]-x(4)-G-K-[ST]. Source: PROSITE PS00017 (https://prosite.expasy.org/PS00017);
    /// Walker JE et al. (1982) EMBO J 1:945–951.
    /// </summary>
    private const string WalkerAPattern = @"[AG].{4}GK[ST]";

    /// <summary>
    /// Finds common protein domains using deterministic PROSITE PATTERN signatures.
    /// </summary>
    /// <remarks>
    /// Only domains with an EXACT PROSITE pattern are detected here: zinc finger C2H2 (PS00028),
    /// WD-repeats (PS00678), and the Walker A / P-loop ATP-binding site (PS00017). Domains whose
    /// only PROSITE signature is a weight-matrix PROFILE — SH3 (PS50002) and PDZ (PS50106) — have
    /// no deterministic pattern and are intentionally NOT detected here: a profile/HMM is a trained
    /// model that cannot be reproduced as an exact regex without fabricating a signature.
    /// </remarks>
    public static IEnumerable<ProteinDomain> FindDomains(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            yield break;

        // Zinc finger C2H2 — exact PROSITE pattern PS00028.
        var zincFingers = FindMotifByPattern(proteinSequence,
            ZincFingerC2H2Pattern, "Zinc Finger C2H2", "PF00096");
        foreach (var zf in zincFingers)
        {
            yield return new ProteinDomain("Zinc Finger C2H2", "PF00096",
                zf.Start, zf.End, zf.Score, "Zinc finger, C2H2 type");
        }

        // WD40 repeats — exact PROSITE pattern PS00678 (WD_REPEATS_1).
        var wd40 = FindMotifByPattern(proteinSequence,
            Wd40RepeatPattern, "WD40 Repeat", "PF00400");
        foreach (var wd in wd40)
        {
            yield return new ProteinDomain("WD40 Repeat", "PF00400",
                wd.Start, wd.End, wd.Score, "WD40/YVTN repeat-like-containing domain");
        }

        // Kinase ATP-binding / P-loop — exact PROSITE pattern PS00017 (Walker A).
        var kinase = FindMotifByPattern(proteinSequence,
            WalkerAPattern, "Protein Kinase", "PF00069");
        foreach (var k in kinase)
        {
            yield return new ProteinDomain("Protein Kinase ATP-binding", "PF00069",
                k.Start, k.End, k.Score, "Protein kinase domain, ATP-binding site");
        }
    }

    #endregion

    #region Profile-HMM Domain Detection (Plan7; opt-in)

    // Bundled CC0 Pfam profile HMMs (see Resources/README.md). These domains have NO deterministic
    // PROSITE pattern — they are trained profile HMMs — so they are detected with the Plan7 engine,
    // not the regex FindDomains path. Provenance + CC0 licence documented in the Evidence artifact.
    private static readonly (string Resource, string Name, string Accession, string Description)[] BundledProfiles =
    {
        ("PF00018_SH3_1.hmm", "SH3", "PF00018", "SH3 domain"),
        ("PF00595_PDZ.hmm",   "PDZ", "PF00595", "PDZ domain"),
        ("PF00400_WD40.hmm",  "WD40", "PF00400", "WD domain, G-beta repeat"),
    };

    private static readonly Lazy<Plan7ProfileHmm[]> LazyBundledHmms = new(() =>
        BundledProfiles.Select(p => Plan7ProfileHmm.LoadEmbedded(p.Resource)).ToArray());

    /// <summary>
    /// Detects SH3, PDZ and WD40 domains by scoring the protein against bundled Pfam profile HMMs
    /// (PF00018, PF00595, PF00400) with the Plan7 Viterbi log-odds algorithm.
    /// </summary>
    /// <remarks>
    /// This is an <b>opt-in</b> path, independent of <see cref="FindDomains"/> (which uses exact
    /// PROSITE patterns and is unchanged). A domain is reported when its Viterbi log-odds score in
    /// bits is at least <paramref name="minBitScore"/>. The reported <c>Score</c> is the Viterbi
    /// bit score; <c>Start</c>/<c>End</c> span the whole input (the glocal full-profile score is a
    /// whole-sequence quantity, not a sub-alignment envelope — see the algorithm doc residual).
    /// </remarks>
    /// <param name="proteinSequence">Amino-acid sequence.</param>
    /// <param name="minBitScore">Minimum Viterbi bit score to report a domain. Default 10 bits.</param>
    public static IEnumerable<ProteinDomain> FindDomainsByHmm(string proteinSequence, double minBitScore = DefaultHmmMinBitScore)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            yield break;

        var hmms = LazyBundledHmms.Value;
        for (int p = 0; p < BundledProfiles.Length; p++)
        {
            double bits = ScoreBits(hmms[p], proteinSequence);
            if (bits >= minBitScore)
            {
                var meta = BundledProfiles[p];
                yield return new ProteinDomain(meta.Name, meta.Accession,
                    0, proteinSequence.Length - 1, bits, meta.Description);
            }
        }
    }

    /// <summary>
    /// Scores a protein against a single bundled Pfam profile HMM, returning the Plan7 Viterbi
    /// log-odds score in <b>bits</b> (natural-log score / ln 2). Opt-in; PROSITE path unchanged.
    /// </summary>
    /// <param name="accession">Pfam accession: "PF00018" (SH3), "PF00595" (PDZ) or "PF00400" (WD40).</param>
    public static double ScoreDomainHmm(string proteinSequence, string accession)
    {
        ArgumentNullException.ThrowIfNull(proteinSequence);
        ArgumentNullException.ThrowIfNull(accession);
        int idx = Array.FindIndex(BundledProfiles, p => p.Accession == accession);
        if (idx < 0)
            throw new ArgumentException($"Unknown bundled Pfam profile: '{accession}'.", nameof(accession));
        return ScoreBits(LazyBundledHmms.Value[idx], proteinSequence);
    }

    // ln 2, for converting natural-log (nat) log-odds scores to bits.
    private const double LogOddsBitsPerNat = 0.69314718055994530941723212145818; // ln 2
    // Default reporting threshold (bits) for FindDomainsByHmm.
    private const double DefaultHmmMinBitScore = 10.0;
    // Default database size Z for E-value reporting (single-target search).
    private const double DefaultDatabaseSize = 1.0;

    private static double ScoreBits(Plan7ProfileHmm hmm, string sequence)
    {
        double nats = hmm.ViterbiScore(sequence);
        if (double.IsNegativeInfinity(nats)) return double.NegativeInfinity;
        return nats / LogOddsBitsPerNat;
    }

    /// <summary>
    /// Detects SH3/PDZ/WD40 domains and reports each hit's Viterbi bit score together with its
    /// HMMER E-value, derived from the bundled profile's <c>STATS LOCAL VITERBI</c> Gumbel
    /// calibration (E = P·Z, P = <c>1 − exp(−exp(−λ(S − μ)))</c>).
    /// </summary>
    /// <remarks>
    /// Opt-in; the exact-PROSITE <see cref="FindDomains"/> path and the bit-score-only
    /// <see cref="FindDomainsByHmm(string, double)"/> overload are unchanged. The reported E-value is
    /// computed from the profile's stored STATS parameters applied to the engine's Viterbi bit score.
    /// Exact <c>hmmsearch</c>-reported-E-value parity additionally requires HMMER's full pipeline
    /// (MSV/bias prefilters and the null2 biased-composition correction); see the algorithm doc.
    /// </remarks>
    /// <param name="proteinSequence">Amino-acid sequence.</param>
    /// <param name="databaseSize">Z — number of target sequences searched. Default 1.</param>
    /// <param name="minBitScore">Minimum Viterbi bit score to report a domain. Default 10 bits.</param>
    public static IEnumerable<ProteinDomainHit> FindDomainHitsByHmm(
        string proteinSequence, double databaseSize = DefaultDatabaseSize, double minBitScore = DefaultHmmMinBitScore)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            yield break;

        var hmms = LazyBundledHmms.Value;
        for (int p = 0; p < BundledProfiles.Length; p++)
        {
            double bits = ScoreBits(hmms[p], proteinSequence);
            if (bits >= minBitScore)
            {
                var meta = BundledProfiles[p];
                double eValue = hmms[p].ViterbiEValue(bits, databaseSize);
                yield return new ProteinDomainHit(meta.Name, meta.Accession,
                    0, proteinSequence.Length - 1, bits, eValue, meta.Description);
            }
        }
    }

    /// <summary>
    /// Scores a protein against a single bundled Pfam profile HMM and returns its Viterbi bit score
    /// and HMMER E-value (E = P·Z) from the profile's Gumbel <c>STATS LOCAL VITERBI</c> calibration.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence.</param>
    /// <param name="accession">Pfam accession: "PF00018" (SH3), "PF00595" (PDZ) or "PF00400" (WD40).</param>
    /// <param name="databaseSize">Z — number of target sequences searched. Default 1.</param>
    public static (double BitScore, double EValue) ScoreDomainHmmEValue(
        string proteinSequence, string accession, double databaseSize = DefaultDatabaseSize)
    {
        ArgumentNullException.ThrowIfNull(proteinSequence);
        ArgumentNullException.ThrowIfNull(accession);
        int idx = Array.FindIndex(BundledProfiles, p => p.Accession == accession);
        if (idx < 0)
            throw new ArgumentException($"Unknown bundled Pfam profile: '{accession}'.", nameof(accession));
        var hmm = LazyBundledHmms.Value[idx];
        double bits = ScoreBits(hmm, proteinSequence);
        double eValue = hmm.ViterbiEValue(bits, databaseSize);
        return (bits, eValue);
    }

    /// <summary>
    /// One automatically-decomposed Pfam domain hit of a multi-domain protein: the family, its
    /// posterior-defined envelope, and its null2-corrected per-domain bit score / independent E-value.
    /// </summary>
    /// <param name="Name">Domain family name (e.g. "WD40").</param>
    /// <param name="Accession">Pfam accession (e.g. "PF00400").</param>
    /// <param name="EnvelopeStart">1-based inclusive envelope start (HMMER <c>env from</c>).</param>
    /// <param name="EnvelopeEnd">1-based inclusive envelope end (HMMER <c>env to</c>).</param>
    /// <param name="BitScore">Null2-corrected per-domain bit score (HMMER per-domain <c>score</c>).</param>
    /// <param name="BiasBits">Per-domain null2 biased-composition correction in bits (HMMER <c>bias</c>).</param>
    /// <param name="IndependentEValue">Independent E-value <c>i-Evalue = Z·exp(lnP)</c> (HMMER <c>i-Evalue</c>).</param>
    /// <param name="Description">Domain description.</param>
    public readonly record struct DomainEnvelopeHit(
        string Name,
        string Accession,
        int EnvelopeStart,
        int EnvelopeEnd,
        double BitScore,
        double BiasBits,
        double IndependentEValue,
        string Description);

    /// <summary>
    /// Decomposes a protein into individual SH3/PDZ/WD40 domains by HMMER's automatic
    /// <c>hmmsearch</c> domain/envelope decomposition (<c>p7_domaindef</c>): posterior region
    /// identification + per-envelope null2-corrected scoring. A multi-domain target (e.g. a multi-WD40
    /// β-propeller) is automatically split into one hit per domain, each with its own envelope and score.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is an <b>opt-in</b> path, independent of <see cref="FindDomains"/>,
    /// <see cref="FindDomainsByHmm(string, double)"/> and <see cref="FindDomainHitsByHmm"/>, all of
    /// which are unchanged. Each bundled profile is scored against the whole sequence; every envelope
    /// whose per-domain bit score is at least <paramref name="minBitScore"/> is reported as a separate
    /// <see cref="DomainEnvelopeHit"/>. See <see cref="Plan7ProfileHmm.FindDomains"/> for the algorithm
    /// and the stochastic-clustering residual for closely-overlapping domains.
    /// </para>
    /// </remarks>
    /// <param name="proteinSequence">Amino-acid sequence.</param>
    /// <param name="minBitScore">Minimum per-domain bit score to report an envelope. Default 10 bits.</param>
    /// <returns>The per-domain envelope hits across all bundled families, in ascending start order.</returns>
    public static IReadOnlyList<DomainEnvelopeHit> FindDomainEnvelopes(
        string proteinSequence, double minBitScore = DefaultHmmMinBitScore)
    {
        var hits = new List<DomainEnvelopeHit>();
        if (string.IsNullOrEmpty(proteinSequence))
            return hits;

        var hmms = LazyBundledHmms.Value;
        for (int p = 0; p < BundledProfiles.Length; p++)
        {
            var meta = BundledProfiles[p];
            foreach (var env in hmms[p].FindDomains(proteinSequence))
            {
                if (env.BitScore < minBitScore) continue;
                hits.Add(new DomainEnvelopeHit(meta.Name, meta.Accession,
                    env.EnvelopeStart, env.EnvelopeEnd, env.BitScore, env.BiasBits,
                    env.IndependentEValue, meta.Description));
            }
        }

        hits.Sort((a, b) => a.EnvelopeStart.CompareTo(b.EnvelopeStart));
        return hits;
    }

    /// <summary>
    /// Decomposes a protein into individual domains against one named bundled Pfam profile HMM
    /// (HMMER <c>hmmsearch</c> domain/envelope decomposition; <see cref="FindDomainEnvelopes"/>).
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence.</param>
    /// <param name="accession">Pfam accession: "PF00018" (SH3), "PF00595" (PDZ) or "PF00400" (WD40).</param>
    /// <returns>The per-domain envelopes for that family, in ascending sequence order.</returns>
    public static IReadOnlyList<Plan7ProfileHmm.DomainEnvelope> FindDomainEnvelopes(
        string proteinSequence, string accession)
    {
        ArgumentNullException.ThrowIfNull(proteinSequence);
        ArgumentNullException.ThrowIfNull(accession);
        int idx = Array.FindIndex(BundledProfiles, p => p.Accession == accession);
        if (idx < 0)
            throw new ArgumentException($"Unknown bundled Pfam profile: '{accession}'.", nameof(accession));
        return LazyBundledHmms.Value[idx].FindDomains(proteinSequence);
    }

    #endregion
}
