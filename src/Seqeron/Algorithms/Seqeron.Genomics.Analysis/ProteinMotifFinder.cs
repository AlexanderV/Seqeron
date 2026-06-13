using System;
using System.Collections.Generic;
using System.Linq;
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
            // Lookahead wrapper enables overlapping match discovery
            regex = new Regex("(?=(" + regexPattern + "))", RegexOptions.IgnoreCase);
        }
        catch
        {
            yield break;
        }

        var matches = regex.Matches(upper);
        foreach (Match match in matches)
        {
            var captured = match.Groups[1];
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
                        sb.Append(prositePattern.Substring(i, end - i + 1));
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

    /// <summary>
    /// Predicts transmembrane helices using hydropathy analysis.
    /// </summary>
    public static IEnumerable<(int Start, int End, double Score)> PredictTransmembraneHelices(
        string proteinSequence,
        int windowSize = 19,
        double threshold = 1.6)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
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
                // End of region
                int start = regionStart.Value;
                int end = i - 1 + windowSize;

                if (end - start >= 15) // Minimum TM helix length
                {
                    yield return (start, Math.Min(end, upper.Length - 1), maxScore);
                }

                regionStart = null;
                maxScore = 0;
            }
        }

        // Handle last region
        if (regionStart != null)
        {
            int start = regionStart.Value;
            int end = hydropathy.Count - 1 + windowSize;

            if (end - start >= 15)
            {
                yield return (start, Math.Min(end, upper.Length - 1), maxScore);
            }
        }
    }

    /// <summary>
    /// Kyte-Doolittle hydropathy scale.
    /// </summary>
    private static readonly Dictionary<char, double> HydropathyScale = new()
    {
        {'A', 1.8}, {'R', -4.5}, {'N', -3.5}, {'D', -3.5}, {'C', 2.5},
        {'Q', -3.5}, {'E', -3.5}, {'G', -0.4}, {'H', -3.2}, {'I', 4.5},
        {'L', 3.8}, {'K', -3.9}, {'M', 1.9}, {'F', 2.8}, {'P', -1.6},
        {'S', -0.8}, {'T', -0.7}, {'W', -0.9}, {'Y', -1.3}, {'V', 4.2}
    };

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

    /// <summary>
    /// Predicts coiled-coil regions using heptad repeat analysis.
    /// </summary>
    public static IEnumerable<(int Start, int End, double Score)> PredictCoiledCoils(
        string proteinSequence,
        int windowSize = 28,
        double threshold = 0.5)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        // Coiled-coil scoring based on heptad positions (a-g)
        // Positions a and d favor hydrophobic residues
        var positionWeights = new Dictionary<int, Dictionary<char, double>>
        {
            // Position a (hydrophobic)
            { 0, new Dictionary<char, double> { {'L', 0.9}, {'I', 0.8}, {'V', 0.7}, {'M', 0.6}, {'A', 0.4} } },
            // Position d (hydrophobic)
            { 3, new Dictionary<char, double> { {'L', 0.9}, {'I', 0.8}, {'V', 0.7}, {'M', 0.6}, {'A', 0.4} } },
            // Position e (charged, negative)
            { 4, new Dictionary<char, double> { {'E', 0.8}, {'D', 0.6}, {'Q', 0.4}, {'N', 0.4} } },
            // Position g (charged, positive)
            { 6, new Dictionary<char, double> { {'K', 0.8}, {'R', 0.7}, {'E', 0.5}, {'Q', 0.4} } }
        };

        var profile = new List<double>();

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double score = 0;
            int count = 0;

            for (int j = 0; j < windowSize; j++)
            {
                int pos = j % 7;
                char aa = upper[i + j];

                if (positionWeights.TryGetValue(pos, out var weights))
                {
                    if (weights.TryGetValue(aa, out double weight))
                    {
                        score += weight;
                    }
                }
                count++;
            }

            profile.Add(score / (windowSize / 7 * 4));
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
                int end = i - 1 + windowSize;

                if (end - start >= 21) // Minimum 3 heptads
                {
                    yield return (start, Math.Min(end, upper.Length - 1), maxScore);
                }

                regionStart = null;
                maxScore = 0;
            }
        }

        if (regionStart != null)
        {
            int start = regionStart.Value;
            int end = profile.Count - 1 + windowSize;

            if (end - start >= 21)
            {
                yield return (start, Math.Min(end, upper.Length - 1), maxScore);
            }
        }
    }

    #endregion

    #region Low Complexity Regions

    /// <summary>
    /// Finds low complexity regions (compositionally biased).
    /// </summary>
    public static IEnumerable<(int Start, int End, char DominantAa, double Frequency)> FindLowComplexityRegions(
        string proteinSequence,
        int windowSize = 12,
        double threshold = 0.4)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < windowSize)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        int? regionStart = null;
        char dominantAa = ' ';
        double maxFreq = 0;

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            string window = upper.Substring(i, windowSize);
            var composition = window.GroupBy(c => c)
                .OrderByDescending(g => g.Count())
                .First();

            double freq = (double)composition.Count() / windowSize;

            if (freq >= threshold)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    dominantAa = composition.Key;
                    maxFreq = freq;
                }
                else
                {
                    if (freq > maxFreq)
                    {
                        dominantAa = composition.Key;
                        maxFreq = freq;
                    }
                }
            }
            else if (regionStart != null)
            {
                yield return (regionStart.Value, i - 1 + windowSize, dominantAa, maxFreq);
                regionStart = null;
                maxFreq = 0;
            }
        }

        if (regionStart != null)
        {
            yield return (regionStart.Value, upper.Length - 1, dominantAa, maxFreq);
        }
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
                string num = "";
                while (i < regexPattern.Length && regexPattern[i] != ',' && regexPattern[i] != '}')
                {
                    num += regexPattern[i];
                    i++;
                }
                if (int.TryParse(num, out int n))
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

    /// <summary>
    /// Finds common protein domains using signature patterns.
    /// </summary>
    public static IEnumerable<ProteinDomain> FindDomains(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            yield break;

        // Check for zinc finger domains
        var zincFingers = FindMotifByPattern(proteinSequence,
            @"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H", "Zinc Finger C2H2", "PF00096");
        foreach (var zf in zincFingers)
        {
            yield return new ProteinDomain("Zinc Finger C2H2", "PF00096",
                zf.Start, zf.End, zf.Score, "Zinc finger, C2H2 type");
        }

        // Check for WD40 repeats
        var wd40 = FindMotifByPattern(proteinSequence,
            @"[LIVMFYWC].{5,12}[WF]D", "WD40 Repeat", "PF00400");
        foreach (var wd in wd40)
        {
            yield return new ProteinDomain("WD40 Repeat", "PF00400",
                wd.Start, wd.End, wd.Score, "WD40/YVTN repeat-like-containing domain");
        }

        // Check for SH3 domain signature
        var sh3 = FindMotifByPattern(proteinSequence,
            @"[LIVMF].{2}[GA]W[FYW].{5,8}[LIVMF]", "SH3", "PF00018");
        foreach (var s in sh3)
        {
            yield return new ProteinDomain("SH3", "PF00018",
                s.Start, s.End, s.Score, "SH3 domain");
        }

        // Check for PDZ domain
        var pdz = FindMotifByPattern(proteinSequence,
            @"[LIVMF][ST][LIVMF].{2}G[LIVMF].{3,4}[LIVMF].{2}[DEN]", "PDZ", "PF00595");
        foreach (var p in pdz)
        {
            yield return new ProteinDomain("PDZ", "PF00595",
                p.Start, p.End, p.Score, "PDZ domain");
        }

        // Check for kinase domain ATP-binding
        var kinase = FindMotifByPattern(proteinSequence,
            @"[AG].{4}GK[ST]", "Protein Kinase", "PF00069");
        foreach (var k in kinase)
        {
            yield return new ProteinDomain("Protein Kinase ATP-binding", "PF00069",
                k.Start, k.End, k.Score, "Protein kinase domain, ATP-binding site");
        }
    }

    #endregion
}
