using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for RNA splice site prediction and intron/exon analysis.
/// </summary>
public static class SpliceSitePredictor
{
    #region Records and Types

    /// <summary>
    /// Represents a predicted splice site.
    /// </summary>
    public readonly record struct SpliceSite(
        int Position,
        SpliceSiteType Type,
        string Motif,
        double Score,
        double Confidence);

    /// <summary>
    /// Result of explicit branch-point detection for a 3' acceptor site,
    /// per the human <c>yUnAy</c> consensus (Gao et al. 2008).
    /// </summary>
    /// <param name="Found">True when a candidate branch point scoring at or above the
    /// requested minimum was located in the upstream search window.</param>
    /// <param name="BranchPointPosition">0-based index of the branch-point adenosine
    /// (motif position 0) in the input sequence, or -1 when none was found.</param>
    /// <param name="DistanceFromAg">Distance in nucleotides from the branch-point
    /// adenosine to the AG (number of bases between the branch point and the splice
    /// site); 0 when none was found.</param>
    /// <param name="Motif">The 5-nt <c>yUnAy</c> motif (positions -3..+1) at the
    /// detected branch point, or the empty string when none was found.</param>
    /// <param name="Score">Conservation-weighted match score in [0, 1] of the best
    /// candidate (1.0 = perfect yUnAy consensus); 0 when none was found.</param>
    /// <param name="PolypyrimidineTractFraction">Pyrimidine (C/U) fraction of the
    /// tract between the branch point and the AG (the 4-24 nt downstream window).</param>
    public readonly record struct BranchPointResult(
        bool Found,
        int BranchPointPosition,
        int DistanceFromAg,
        string Motif,
        double Score,
        double PolypyrimidineTractFraction);

    /// <summary>
    /// Types of splice sites.
    /// </summary>
    public enum SpliceSiteType
    {
        Donor,      // 5' splice site (GT/GU)
        Acceptor,   // 3' splice site (AG)
        Branch,     // Branch point (A)
        U12Donor,   // Minor spliceosome donor (AT/AU)
        U12Acceptor // Minor spliceosome acceptor (AC)
    }

    /// <summary>
    /// Represents a predicted intron.
    /// </summary>
    public readonly record struct Intron(
        int Start,
        int End,
        int Length,
        SpliceSite DonorSite,
        SpliceSite AcceptorSite,
        SpliceSite? BranchPoint,
        string Sequence,
        IntronType Type,
        double Score);

    /// <summary>
    /// Types of introns.
    /// </summary>
    public enum IntronType
    {
        U2,     // Major spliceosome (GT-AG)
        U12,    // Minor spliceosome (AT-AC)
        GcAg,   // GC-AG variant
        Unknown
    }

    /// <summary>
    /// Represents an exon.
    /// </summary>
    public readonly record struct Exon(
        int Start,
        int End,
        int Length,
        ExonType Type,
        int? Phase,
        string Sequence);

    /// <summary>
    /// Types of exons.
    /// </summary>
    public enum ExonType
    {
        Initial,    // First exon (contains start codon)
        Internal,   // Internal exon
        Terminal,   // Last exon (contains stop codon)
        Single      // Single exon gene
    }

    /// <summary>
    /// Complete gene structure prediction.
    /// </summary>
    public readonly record struct GeneStructure(
        IReadOnlyList<Exon> Exons,
        IReadOnlyList<Intron> Introns,
        string SplicedSequence,
        double OverallScore);

    #endregion

    #region Position Weight Matrices

    // Donor site IUPAC consensus: MAG|GURAGU (M=A/C, R=A/G)
    // Positions -3 to +5 relative to exon-intron boundary (9 positions total).
    // Binary weights: 1.0 = matches consensus nucleotide(s), 0.0 = no match.
    // Sources: Shapiro & Senapathy (1987) NAR 15:7155; Mount (1982) NAR 10:459;
    //          Burge et al. (1999) The RNA World, CSHL Press.
    private static readonly Dictionary<int, Dictionary<char, double>> DonorPwm = new()
    {
        { -3, new Dictionary<char, double> { {'A', 1.0}, {'C', 1.0}, {'G', 0.0}, {'U', 0.0} } }, // M (A or C)
        { -2, new Dictionary<char, double> { {'A', 1.0}, {'C', 0.0}, {'G', 0.0}, {'U', 0.0} } }, // A
        { -1, new Dictionary<char, double> { {'A', 0.0}, {'C', 0.0}, {'G', 1.0}, {'U', 0.0} } }, // G
        { 0,  new Dictionary<char, double> { {'A', 0.0}, {'C', 0.0}, {'G', 1.0}, {'U', 0.0} } }, // G (invariant)
        { 1,  new Dictionary<char, double> { {'A', 0.0}, {'C', 0.0}, {'G', 0.0}, {'U', 1.0} } }, // U (invariant)
        { 2,  new Dictionary<char, double> { {'A', 1.0}, {'C', 0.0}, {'G', 1.0}, {'U', 0.0} } }, // R (A or G)
        { 3,  new Dictionary<char, double> { {'A', 1.0}, {'C', 0.0}, {'G', 0.0}, {'U', 0.0} } }, // A
        { 4,  new Dictionary<char, double> { {'A', 0.0}, {'C', 0.0}, {'G', 1.0}, {'U', 0.0} } }, // G
        { 5,  new Dictionary<char, double> { {'A', 0.0}, {'C', 0.0}, {'G', 0.0}, {'U', 1.0} } }  // U
    };

    // Consensus acceptor site: (Y)nNCAG|G
    // Position -15 to +1 relative to splice site
    private static readonly Dictionary<int, Dictionary<char, double>> AcceptorPwm = new()
    {
        { -15, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.30}, {'G', 0.10}, {'U', 0.50} } },
        { -10, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.30}, {'G', 0.10}, {'U', 0.50} } },
        { -5,  new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { -4,  new Dictionary<char, double> { {'A', 0.05}, {'C', 0.40}, {'G', 0.05}, {'U', 0.50} } },
        { -3,  new Dictionary<char, double> { {'A', 0.05}, {'C', 0.70}, {'G', 0.05}, {'U', 0.20} } },
        { -2,  new Dictionary<char, double> { {'A', 1.00}, {'C', 0.00}, {'G', 0.00}, {'U', 0.00} } }, // A
        { -1,  new Dictionary<char, double> { {'A', 0.00}, {'C', 0.00}, {'G', 1.00}, {'U', 0.00} } }, // G
        { 0,   new Dictionary<char, double> { {'A', 0.20}, {'C', 0.15}, {'G', 0.50}, {'U', 0.15} } }
    };

    // --- Human branch-point consensus (Gao et al. 2008) ---
    // Consensus BPS = yUnAy at motif positions -3..+1, with the branch-point
    // adenosine at position 0. Per-position conservation (lariat RT-PCR, n=181
    // confirmed branch sites): pos -3 pyrimidine C+U = 47.0%+32.0% = 79.0%;
    // pos -2 U = 74.6%; pos 0 A = 92.3%; pos +1 pyrimidine C+U = 33.1%+42.0% = 75.1%;
    // pos -1 = n (any, uninformative). Source: Gao K, Masuda A, Matsuura T, Ohno K
    // (2008) "Human branch point consensus sequence is yUnAy", Nucleic Acids Res
    // 36(7):2257-2267, Table 1 / Figure 2 (DOI 10.1093/nar/gkn073).
    private const double BpPos3PyrimidineConservation = 0.790;  // pos -3 y (C+U)
    private const double BpPos2UracilConservation = 0.746;      // pos -2 U
    private const double BpBranchAdenosineConservation = 0.923; // pos  0 A (branch point)
    private const double BpPos1PyrimidineConservation = 0.751;  // pos +1 y (C+U)

    // Branch-point search window upstream of the 3' acceptor AG. Gao et al. (2008)
    // report 83% of branch sites at positions -34..-21 relative to the 3' end of the
    // intron (median -26, mean -27.7 +/- 7.6); Mercer et al. (2015, Genome Res
    // 25:290-303) report most branch sites 19-35 nt from the 3'ss; the spliceosome
    // literature gives an outer envelope of ~18-50 nt. The default window uses the
    // conservative 18-40 nt envelope that brackets the Gao -34..-21 core.
    private const int BranchPointMinDistanceFromAg = 18; // nt upstream of AG (inclusive)
    private const int BranchPointMaxDistanceFromAg = 40; // nt upstream of AG (inclusive)

    // Polypyrimidine tract between the branch point and the AG spans 4-24 nt
    // downstream of the branch point (Gao et al. 2008).
    private const int BranchPointPptMinSpan = 4;  // nt downstream of branch point
    private const int BranchPointPptMaxSpan = 24; // nt downstream of branch point

    // Branch point consensus: YNYURAC (Y=C/U, R=A/G, N=any)
    // Typically 18-40 nt upstream of 3' splice site
    private static readonly Dictionary<int, Dictionary<char, double>> BranchPointPwm = new()
    {
        { 0, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { 1, new Dictionary<char, double> { {'A', 0.25}, {'C', 0.25}, {'G', 0.25}, {'U', 0.25} } },
        { 2, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { 3, new Dictionary<char, double> { {'A', 0.05}, {'C', 0.05}, {'G', 0.05}, {'U', 0.85} } },
        { 4, new Dictionary<char, double> { {'A', 0.60}, {'C', 0.05}, {'G', 0.30}, {'U', 0.05} } },
        { 5, new Dictionary<char, double> { {'A', 0.95}, {'C', 0.02}, {'G', 0.02}, {'U', 0.01} } }, // Branch A
        { 6, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.60}, {'G', 0.10}, {'U', 0.20} } }
    };

    #endregion

    #region Splice Site Prediction

    /// <summary>
    /// Finds all potential donor (5') splice sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindDonorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 6)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 0; i <= upper.Length - 6; i++)
        {
            // Check for canonical GT/GU
            if (upper[i] == 'G' && upper[i + 1] == 'U')
            {
                double score = ScoreDonorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // Non-canonical GC
            else if (includeNonCanonical && upper[i] == 'G' && upper[i + 1] == 'C')
            {
                double score = ScoreDonorSite(upper, i); // GC naturally scores lower (position +1 mismatches U consensus)
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // U12 spliceosome AT/AU
            else if (includeNonCanonical && upper[i] == 'A' && (upper[i + 1] == 'U' || upper[i + 1] == 'T'))
            {
                double score = ScoreU12DonorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.U12Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Finds all potential acceptor (3') splice sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindAcceptorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 20)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 15; i < upper.Length - 1; i++)
        {
            // Check for canonical AG
            if (upper[i] == 'A' && upper[i + 1] == 'G')
            {
                double score = ScoreAcceptorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i + 1, // Position of G in AG (last intron nucleotide)
                        Type: SpliceSiteType.Acceptor,
                        Motif: GetMotifContext(upper, i, 15, 2),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // U12 spliceosome AC
            else if (includeNonCanonical && upper[i] == 'A' && upper[i + 1] == 'C')
            {
                double score = ScoreU12AcceptorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i + 1,
                        Type: SpliceSiteType.U12Acceptor,
                        Motif: GetMotifContext(upper, i, 15, 2),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Detects the explicit branch point upstream of a 3' acceptor AG using the human
    /// <c>yUnAy</c> branch-point consensus (Gao et al. 2008). This is an opt-in,
    /// additive analysis; the default <see cref="FindAcceptorSites"/> scoring
    /// (PWM consensus + polypyrimidine-tract fraction) is unchanged.
    /// </summary>
    /// <remarks>
    /// The branch-point adenosine sits at motif position 0 of the consensus
    /// <c>yUnAy</c> (positions -3..+1). It is searched in the
    /// <see cref="BranchPointMinDistanceFromAg"/>..<see cref="BranchPointMaxDistanceFromAg"/>
    /// (default 18-40) nt window upstream of the AG, bracketing the -34..-21 core
    /// reported by Gao et al. (2008). Each candidate is scored by a conservation-weighted
    /// fraction over the four informative positions (-3 pyrimidine, -2 U, 0 A, +1
    /// pyrimidine); position -1 is uninformative ('n'). The polypyrimidine-tract fraction
    /// between the branch point and the AG is reported alongside the best candidate.
    /// </remarks>
    /// <param name="sequence">The pre-mRNA / intron sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="acceptorAgPosition">0-based index of the terminal <c>G</c> of the 3'
    /// acceptor AG (i.e. the <see cref="SpliceSite.Position"/> reported by
    /// <see cref="FindAcceptorSites"/>).</param>
    /// <param name="minScore">Minimum conservation-weighted score in [0, 1] for a
    /// candidate to count as found. Defaults to 0.5.</param>
    /// <returns>The best branch-point candidate in the upstream window, or a
    /// <see cref="BranchPointResult"/> with <c>Found = false</c> when none qualifies.</returns>
    public static BranchPointResult FindAcceptorBranchPoint(
        string sequence,
        int acceptorAgPosition,
        double minScore = 0.5)
    {
        var notFound = new BranchPointResult(false, -1, 0, string.Empty, 0, 0);

        if (string.IsNullOrEmpty(sequence) || acceptorAgPosition <= 0 ||
            acceptorAgPosition >= sequence.Length)
            return notFound;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // The reported acceptor Position is the terminal G of the AG; the splice
        // junction (3' end of the intron) lies just after it. Distance is measured
        // from the branch-point adenosine to that 3' end.
        int agEnd = acceptorAgPosition; // index of the G of AG (last intronic nt)

        BranchPointResult best = notFound;

        // Scan the upstream window. distance = (agEnd - branchA), so branchA = agEnd - distance.
        for (int distance = BranchPointMinDistanceFromAg;
             distance <= BranchPointMaxDistanceFromAg;
             distance++)
        {
            int branchA = agEnd - distance;
            // Need motif positions -3..+1 → indices [branchA-3, branchA+1] in bounds.
            if (branchA - 3 < 0 || branchA + 1 >= upper.Length)
                continue;

            double score = ScoreBranchPointConsensus(upper, branchA);
            if (score < minScore || score <= best.Score)
                continue;

            string motif = upper.Substring(branchA - 3, 5); // yUnAy (positions -3..+1)
            double pptFraction = ComputeBranchPointPptFraction(upper, branchA, agEnd);

            best = new BranchPointResult(
                Found: true,
                BranchPointPosition: branchA,
                DistanceFromAg: distance,
                Motif: motif,
                Score: score,
                PolypyrimidineTractFraction: pptFraction);
        }

        return best;
    }

    /// <summary>
    /// Finds branch point sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindBranchPoints(
        string sequence,
        int searchStart = 0,
        int searchEnd = -1,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 7)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        int end = searchEnd < 0 ? upper.Length - 7 : Math.Min(searchEnd, upper.Length - 7);

        for (int i = Math.Max(0, searchStart); i <= end; i++)
        {
            double score = ScoreBranchPoint(upper, i);
            if (score >= minScore)
            {
                yield return new SpliceSite(
                    Position: i + 5, // Position of branch A
                    Type: SpliceSiteType.Branch,
                    Motif: upper.Substring(i, 7),
                    Score: score,
                    Confidence: CalculateConfidence(score, 0.5, 1.0));
            }
        }
    }

    #endregion

    #region Scoring Functions

    private static double ScoreDonorSite(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        foreach (var (offset, weights) in DonorPwm)
        {
            int pos = position + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += weight;
                    count++;
                }
            }
        }

        return count > 0 ? score / count : 0;
    }

    private static double ScoreAcceptorSite(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        // Score polypyrimidine tract
        int pptScore = 0;
        for (int i = position - 15; i < position - 3; i++)
        {
            if (i >= 0 && i < sequence.Length)
            {
                if (sequence[i] == 'C' || sequence[i] == 'U')
                    pptScore++;
            }
        }

        score += pptScore / 12.0 * 2; // PPT contribution

        // Score AG and context — PWM offsets are relative to splice site (position + 2)
        foreach (var (offset, weights) in AcceptorPwm)
        {
            int pos = position + 2 + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += Math.Log2(weight / 0.25 + 0.01);
                    count++;
                }
            }
        }

        double normalized = (score / (count + 1) + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    /// <summary>
    /// Scores a candidate branch point against the human <c>yUnAy</c> consensus
    /// (Gao et al. 2008) as a conservation-weighted match fraction in [0, 1].
    /// <paramref name="branchA"/> is the 0-based index of the candidate branch-point
    /// adenosine (motif position 0). Position -1 ('n') is uninformative and ignored.
    /// A perfect canonical branch point (y at -3, U at -2, A at 0, y at +1) scores 1.0.
    /// </summary>
    private static double ScoreBranchPointConsensus(string sequence, int branchA)
    {
        // Caller guarantees [branchA-3, branchA+1] are in bounds.
        char pos3 = sequence[branchA - 3]; // y (C/U)
        char pos2 = sequence[branchA - 2]; // U
        char pos0 = sequence[branchA];     // A (branch point)
        char pos1 = sequence[branchA + 1]; // y (C/U)

        double matched = 0;
        if (pos3 == 'C' || pos3 == 'U') matched += BpPos3PyrimidineConservation;
        if (pos2 == 'U') matched += BpPos2UracilConservation;
        if (pos0 == 'A') matched += BpBranchAdenosineConservation;
        if (pos1 == 'C' || pos1 == 'U') matched += BpPos1PyrimidineConservation;

        const double maxScore = BpPos3PyrimidineConservation + BpPos2UracilConservation +
                                BpBranchAdenosineConservation + BpPos1PyrimidineConservation;
        return matched / maxScore;
    }

    /// <summary>
    /// Pyrimidine (C/U) fraction of the polypyrimidine tract between a branch point and
    /// the AG. The tract spans 4-24 nt downstream of the branch point (Gao et al. 2008);
    /// this measures the bases strictly between the branch-point adenosine and the AG,
    /// clamped to that 4-24 nt span.
    /// </summary>
    private static double ComputeBranchPointPptFraction(string sequence, int branchA, int agEnd)
    {
        // The AG occupies [agEnd-1, agEnd]; the tract lies between the branch point
        // and the A of the AG. Start one base after the branch point.
        int tractStart = branchA + BranchPointPptMinSpan;
        int tractEnd = Math.Min(branchA + BranchPointPptMaxSpan, agEnd - 2); // up to base before A of AG
        if (tractEnd < tractStart)
        {
            tractStart = branchA + 1;
            tractEnd = agEnd - 2;
        }

        int count = 0;
        int total = 0;
        for (int i = tractStart; i <= tractEnd; i++)
        {
            if (i < 0 || i >= sequence.Length)
                continue;
            total++;
            if (sequence[i] == 'C' || sequence[i] == 'U')
                count++;
        }

        return total > 0 ? (double)count / total : 0;
    }

    private static double ScoreBranchPoint(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        foreach (var (offset, weights) in BranchPointPwm)
        {
            int pos = position + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += Math.Log2(weight / 0.25 + 0.01);
                    count++;
                }
            }
        }

        double normalized = (score / count + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    private static double ScoreU12DonorSite(string sequence, int position)
    {
        // Simplified U12 scoring - /ATATCC/ consensus
        string motif = position + 6 <= sequence.Length
            ? sequence.Substring(position, 6)
            : "";

        int matches = 0;
        string consensus = "AUAUCC";
        for (int i = 0; i < Math.Min(motif.Length, consensus.Length); i++)
        {
            if (motif[i] == consensus[i])
                matches++;
        }

        return matches / 6.0;
    }

    private static double ScoreU12AcceptorSite(string sequence, int position)
    {
        // U12-type 3' splice site consensus: YCCAC (Hall & Padgett 1994, Jackson 1991)
        // position = index of A in the AC dinucleotide
        // Expected upstream: sequence[pos-3]=Y, sequence[pos-2]=C, sequence[pos-1]=C
        if (position < 3 || position + 1 >= sequence.Length)
            return 0;

        if (sequence[position] != 'A' || sequence[position + 1] != 'C')
            return 0;

        double score = 0;

        // YCCAC consensus: check C at position-1 (second C of YCC)
        if (sequence[position - 1] == 'C')
            score += 1.0;

        // YCCAC consensus: check C at position-2 (first C of YCC)
        if (position - 2 >= 0 && sequence[position - 2] == 'C')
            score += 1.0;

        // YCCAC consensus: check Y (pyrimidine) at position-3
        if (position - 3 >= 0 && (sequence[position - 3] == 'C' || sequence[position - 3] == 'U'))
            score += 0.5;

        // Polypyrimidine tract upstream
        int pptCount = 0;
        int pptTotal = 0;
        for (int i = position - 15; i < position - 3; i++)
        {
            if (i >= 0 && i < sequence.Length)
            {
                pptTotal++;
                if (sequence[i] == 'C' || sequence[i] == 'U')
                    pptCount++;
            }
        }

        double pptFraction = pptTotal > 0 ? (double)pptCount / pptTotal : 0;
        score += pptFraction;

        // Max possible: 1.0 (C at -1) + 1.0 (C at -2) + 0.5 (Y at -3) + 1.0 (PPT) = 3.5
        double normalized = score / 3.5;
        return Math.Max(0, Math.Min(1, normalized));
    }

    #endregion

    #region Intron Prediction

    /// <summary>
    /// Predicts introns by pairing donor and acceptor sites.
    /// </summary>
    public static IEnumerable<Intron> PredictIntrons(
        string sequence,
        int minIntronLength = 60,
        int maxIntronLength = 100000,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        var donors = FindDonorSites(upper, minScore * 0.8, true).ToList();
        var acceptors = FindAcceptorSites(upper, minScore * 0.8, true).ToList();

        foreach (var donor in donors)
        {
            foreach (var acceptor in acceptors)
            {
                int intronLength = acceptor.Position - donor.Position + 1;
                if (intronLength < minIntronLength || intronLength > maxIntronLength)
                    continue;

                // Determine intron type
                var type = DetermineIntronType(donor, acceptor);

                // Find branch point
                int searchStart = Math.Max(0, acceptor.Position - 50);
                int searchEnd = acceptor.Position - 18;
                var branchPoints = FindBranchPoints(upper, searchStart, searchEnd, 0.4).ToList();
                var bestBranch = branchPoints.OrderByDescending(b => b.Score).FirstOrDefault();

                double combinedScore = bestBranch.Score > 0
                    ? (donor.Score + acceptor.Score + bestBranch.Score) / 3
                    : (donor.Score + acceptor.Score) / 2;

                if (combinedScore >= minScore)
                {
                    yield return new Intron(
                        Start: donor.Position,
                        End: acceptor.Position,
                        Length: intronLength,
                        DonorSite: donor,
                        AcceptorSite: acceptor,
                        BranchPoint: bestBranch.Score > 0 ? bestBranch : null,
                        Sequence: upper.Substring(donor.Position, intronLength),
                        Type: type,
                        Score: combinedScore);
                }
            }
        }
    }

    private static IntronType DetermineIntronType(SpliceSite donor, SpliceSite acceptor)
    {
        if (donor.Type == SpliceSiteType.U12Donor || acceptor.Type == SpliceSiteType.U12Acceptor)
            return IntronType.U12;

        string donorMotif = donor.Motif.ToUpperInvariant();
        // Motif = GetMotifContext(seq, pos, 3, 6): donor dinucleotide at offset min(3, pos)
        int offset = Math.Min(3, donor.Position);
        if (donorMotif.Length >= offset + 2)
        {
            string dinuc = donorMotif.Substring(offset, 2);
            if (dinuc == "GC")
                return IntronType.GcAg;
            if (dinuc == "GU" || dinuc == "GT")
                return IntronType.U2;
        }

        return IntronType.Unknown;
    }

    #endregion

    #region Gene Structure Prediction

    /// <summary>
    /// Predicts exon/intron structure of a gene.
    /// </summary>
    public static GeneStructure PredictGeneStructure(
        string sequence,
        int minExonLength = 30,
        int minIntronLength = 60,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new GeneStructure(
                new List<Exon>(),
                new List<Intron>(),
                "",
                0);
        }

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Get all introns
        var introns = PredictIntrons(upper, minIntronLength, 100000, minScore)
            .OrderByDescending(i => i.Score)
            .ToList();

        // Select non-overlapping introns greedily
        var selectedIntrons = SelectNonOverlappingIntrons(introns);

        // Derive exons from intron positions
        var exons = DeriveExons(upper, selectedIntrons, minExonLength);

        // Generate spliced sequence from the SAME exon set that DeriveExons reports,
        // so that splicedSequence == concat(reportedExons) by construction and the
        // coverage invariants stay internally consistent (INV-3 / INV-4). Sub-threshold
        // inter-intron regions filtered by minExonLength are excluded from BOTH the exon
        // list AND the spliced product, preserving the documented filtering semantics.
        string splicedSequence = GenerateSplicedSequence(exons);

        double overallScore = selectedIntrons.Count > 0
            ? selectedIntrons.Average(i => i.Score)
            : 0;

        return new GeneStructure(
            Exons: exons,
            Introns: selectedIntrons,
            SplicedSequence: splicedSequence,
            OverallScore: overallScore);
    }

    private static List<Intron> SelectNonOverlappingIntrons(List<Intron> introns)
    {
        var selected = new List<Intron>();
        var used = new HashSet<int>();

        foreach (var intron in introns)
        {
            bool overlaps = false;
            for (int pos = intron.Start; pos <= intron.End; pos++)
            {
                if (used.Contains(pos))
                {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps)
            {
                selected.Add(intron);
                for (int pos = intron.Start; pos <= intron.End; pos++)
                    used.Add(pos);
            }
        }

        return selected.OrderBy(i => i.Start).ToList();
    }

    private static List<Exon> DeriveExons(string sequence, List<Intron> introns, int minExonLength)
    {
        var exons = new List<Exon>();

        if (introns.Count == 0)
        {
            // Single exon gene
            exons.Add(new Exon(
                Start: 0,
                End: sequence.Length - 1,
                Length: sequence.Length,
                Type: ExonType.Single,
                Phase: 0,
                Sequence: sequence));
            return exons;
        }

        int currentPos = 0;

        for (int i = 0; i < introns.Count; i++)
        {
            var intron = introns[i];
            int exonEnd = intron.Start - 1;

            if (exonEnd - currentPos + 1 >= minExonLength)
            {
                var exonType = i == 0 ? ExonType.Initial : ExonType.Internal;
                exons.Add(new Exon(
                    Start: currentPos,
                    End: exonEnd,
                    Length: exonEnd - currentPos + 1,
                    Type: exonType,
                    Phase: CalculatePhase(exons),
                    Sequence: sequence.Substring(currentPos, exonEnd - currentPos + 1)));
            }

            currentPos = intron.End + 1;
        }

        // Terminal exon
        if (sequence.Length - currentPos >= minExonLength)
        {
            exons.Add(new Exon(
                Start: currentPos,
                End: sequence.Length - 1,
                Length: sequence.Length - currentPos,
                Type: ExonType.Terminal,
                Phase: CalculatePhase(exons),
                Sequence: sequence.Substring(currentPos)));
        }

        return exons;
    }

    private static int CalculatePhase(List<Exon> previousExons)
    {
        int totalLength = previousExons.Sum(e => e.Length);
        return totalLength % 3;
    }

    /// <summary>
    /// Builds the spliced (mature mRNA) sequence by concatenating the exon
    /// sequences exactly as reported by <see cref="DeriveExons"/>. Using the same
    /// exon set guarantees splicedSequence == concat(reportedExons) (INV-4) and
    /// keeps coverage internally consistent with the reported exons (INV-3),
    /// including sub-threshold inter-intron regions that DeriveExons drops.
    /// </summary>
    private static string GenerateSplicedSequence(List<Exon> exons)
    {
        if (exons.Count == 0)
            return string.Empty;

        var sb = new StringBuilder();
        foreach (var exon in exons.OrderBy(e => e.Start))
            sb.Append(exon.Sequence);

        return sb.ToString();
    }

    #endregion

    #region Alternative Splicing

    /// <summary>
    /// Detects potential alternative splicing patterns.
    /// </summary>
    public static IEnumerable<(string Type, int Position, string Description)> DetectAlternativeSplicing(
        string sequence,
        double minScore = 0.4)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        var donors = FindDonorSites(sequence, minScore).ToList();
        var acceptors = FindAcceptorSites(sequence, minScore).ToList();

        // Exon skipping: multiple introns with shared boundaries
        for (int i = 0; i < donors.Count; i++)
        {
            int validAcceptors = acceptors.Count(a => a.Position > donors[i].Position + 60);
            if (validAcceptors > 1)
            {
                yield return ("ExonSkipping", donors[i].Position,
                    $"Potential exon skipping from donor at {donors[i].Position}");
            }
        }

        // Alternative 5' splice sites (multiple donors before same acceptor)
        var donorGroups = donors.GroupBy(d => d.Position / 50).Where(g => g.Count() > 1);
        foreach (var group in donorGroups)
        {
            yield return ("Alt5SS", group.First().Position,
                $"Alternative 5' splice sites at positions {string.Join(", ", group.Select(d => d.Position))}");
        }

        // Alternative 3' splice sites (multiple acceptors after same donor)
        var acceptorGroups = acceptors.GroupBy(a => a.Position / 50).Where(g => g.Count() > 1);
        foreach (var group in acceptorGroups)
        {
            yield return ("Alt3SS", group.First().Position,
                $"Alternative 3' splice sites at positions {string.Join(", ", group.Select(a => a.Position))}");
        }
    }

    /// <summary>
    /// Finds retained intron candidates (introns that might be retained in some transcripts).
    /// </summary>
    public static IEnumerable<Intron> FindRetainedIntronCandidates(
        string sequence,
        double minScore = 0.5)
    {
        var introns = PredictIntrons(sequence, 60, 500, minScore).ToList();

        // Short introns with moderate scores are candidates for retention
        foreach (var intron in introns.Where(i => i.Length < 500 && i.Score < 0.8))
        {
            yield return intron;
        }
    }

    #endregion

    #region Utility Methods

    private static string GetMotifContext(string sequence, int position, int upstream, int downstream)
    {
        int start = Math.Max(0, position - upstream);
        int end = Math.Min(sequence.Length, position + downstream);
        return sequence.Substring(start, end - start);
    }

    private static double CalculateConfidence(double score, double minExpected, double maxExpected)
    {
        return Math.Max(0, Math.Min(1, (score - minExpected) / (maxExpected - minExpected)));
    }

    /// <summary>
    /// Calculates MaxEntScan-like score for splice sites.
    /// </summary>
    public static double CalculateMaxEntScore(string motif, SpliceSiteType type)
    {
        if (string.IsNullOrEmpty(motif))
            return 0;

        string upper = motif.ToUpperInvariant().Replace('T', 'U');

        if (type == SpliceSiteType.Donor)
        {
            // Score 9-mer around GT
            double score = 0;
            for (int i = 0; i < upper.Length && i < 9; i++)
            {
                int offset = i - 3;
                if (DonorPwm.TryGetValue(offset, out var weights))
                {
                    if (weights.TryGetValue(upper[i], out double w))
                        score += Math.Log2(w + 0.01);
                }
            }
            return score;
        }
        else if (type == SpliceSiteType.Acceptor)
        {
            double score = 0;
            for (int i = 0; i < upper.Length; i++)
            {
                int offset = i - 15;
                if (AcceptorPwm.TryGetValue(offset, out var weights))
                {
                    if (weights.TryGetValue(upper[i], out double w))
                        score += Math.Log2(w + 0.01);
                }
            }
            return score;
        }

        return 0;
    }

    /// <summary>
    /// Checks if a sequence position is within a coding region (simple heuristic).
    /// </summary>
    public static bool IsWithinCodingRegion(string sequence, int position, int frame = 0)
    {
        if (position < 0 || position >= sequence.Length)
            return false;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Look for start codon upstream
        int searchStart = Math.Max(0, position - 300);
        for (int i = searchStart; i <= position - 3; i++)
        {
            if (upper.Substring(i, 3) == "AUG")
            {
                // Check if position is in frame
                return (position - i) % 3 == frame;
            }
        }

        return false;
    }

    #endregion
}
