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
    /// Signal peptide prediction result.
    /// </summary>
    public readonly record struct SignalPeptide(
        int CleavagePosition,
        string NRegion,
        string HRegion,
        string CRegion,
        double Score,
        double Probability);

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

        // PROSITE PS00007: Tyrosine kinase phosphorylation site — https://prosite.expasy.org/PS00007
        ["PS00007"] = new("PS00007", "TYR_PHOSPHO_SITE", "[RK]-x(2)-[DE]-x(3)-Y", @"[RK].{2}[DE].{3}Y", "Tyrosine kinase phosphorylation site"),

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
        ["NLS1"] = new("NLS1", "NLS_MONOPARTITE", "[KR]-[KR]-x-[KR]", @"[KR][KR].[KR]", "Monopartite nuclear localization signal"),

        // NES — Nuclear export signal (CRM1-dependent)
        // Simplified from la Cour et al. (2004) Nucleic Acids Res 32:W142-5
        // General NES consensus: Φ1-x(2,3)-Φ2-x(2,3)-Φ3-x-Φ4 where Φ = hydrophobic
        // L-enriched form used here reduces false positives
        ["NES1"] = new("NES1", "NES", "L-x(2,3)-[LIVFM]-x(2,3)-L-x-[LI]", @"L.{2,3}[LIVFM].{2,3}L.[LI]", "Nuclear export signal"),

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
                // Character class
                int end = prositePattern.IndexOf(']', i);
                if (end > i)
                {
                    sb.Append(prositePattern.Substring(i, end - i + 1));
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
            else
            {
                i++;
            }
        }

        return sb.ToString();
    }

    #endregion

    #region Signal Peptide Prediction

    /// <summary>
    /// Predicts signal peptide cleavage site.
    /// </summary>
    public static SignalPeptide? PredictSignalPeptide(string proteinSequence, int maxLength = 70)
    {
        if (string.IsNullOrEmpty(proteinSequence) || proteinSequence.Length < 15)
            return null;

        string upper = proteinSequence.ToUpperInvariant();
        int searchLength = Math.Min(maxLength, upper.Length);

        // Signal peptide has three regions:
        // N-region (1-5 aa, positively charged)
        // H-region (7-15 aa, hydrophobic)
        // C-region (3-7 aa, polar, ends with small amino acids A, G, S)

        double bestScore = 0;
        int bestCleavage = -1;
        string bestN = "", bestH = "", bestC = "";

        // Scan for cleavage sites
        for (int cleavage = 15; cleavage <= Math.Min(35, searchLength); cleavage++)
        {
            // Check -3, -1 rule (small amino acids at positions -3 and -1)
            char m3 = upper[cleavage - 3];
            char m1 = upper[cleavage - 1];

            if (!IsSmallAminoAcid(m3) || !IsSmallAminoAcid(m1))
                continue;

            // Score the regions
            string nRegion = upper.Substring(0, Math.Min(5, cleavage / 3));
            string hRegion = upper.Substring(nRegion.Length, cleavage - nRegion.Length - 5);
            string cRegion = upper.Substring(cleavage - 5, 5);

            if (hRegion.Length < 5)
                continue;

            double nScore = ScoreNRegion(nRegion);
            double hScore = ScoreHydrophobicRegion(hRegion);
            double cScore = ScoreCRegion(cRegion);

            double totalScore = nScore * 0.2 + hScore * 0.5 + cScore * 0.3;

            if (totalScore > bestScore)
            {
                bestScore = totalScore;
                bestCleavage = cleavage;
                bestN = nRegion;
                bestH = hRegion;
                bestC = cRegion;
            }
        }

        if (bestCleavage < 0 || bestScore < 0.4)
            return null;

        return new SignalPeptide(
            CleavagePosition: bestCleavage,
            NRegion: bestN,
            HRegion: bestH,
            CRegion: bestC,
            Score: bestScore,
            Probability: Math.Min(1.0, bestScore * 1.2));
    }

    private static bool IsSmallAminoAcid(char aa)
    {
        return "AGSTC".Contains(aa);
    }

    private static double ScoreNRegion(string region)
    {
        // N-region should have positive charges (K, R)
        int positiveCharges = region.Count(c => c == 'K' || c == 'R');
        return Math.Min(1.0, positiveCharges * 0.4);
    }

    private static double ScoreHydrophobicRegion(string region)
    {
        // H-region should be hydrophobic
        int hydrophobic = region.Count(c => "AILMFVW".Contains(c));
        double ratio = (double)hydrophobic / region.Length;
        return ratio > 0.5 ? 1.0 : ratio * 2;
    }

    private static double ScoreCRegion(string region)
    {
        // C-region: polar, small residues at -3, -1
        int smallPolar = region.Count(c => "AGSTN".Contains(c));
        return Math.Min(1.0, smallPolar * 0.25);
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

        // Minimum non-zero score for fully unconstrained patterns
        return totalIC > 0 ? totalIC : 0.01;
    }

    /// <summary>
    /// Calculates E-value: expected number of random matches in a sequence.
    /// E = (N - L + 1) × P_random, where P_random = 2^(-IC_total).
    /// Reference: Altschul SF et al. (1990). "Basic Local Alignment Search Tool."
    /// J Mol Biol 215:403-410.
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
