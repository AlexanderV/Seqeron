using System.Globalization;
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
        U12Donor,   // Minor spliceosome donor (/RTATCCTY: AT-AC or GT-AG U12 introns)
        U12Acceptor // Minor spliceosome acceptor (AC or AG terminus; U12 branch site UUCCUUAAC)
    }

    /// <summary>
    /// Scoring model for <see cref="FindDonorSites(string, double, bool, DonorScoringModel, ShapiroSenapathyCategory)"/>.
    /// </summary>
    public enum DonorScoringModel
    {
        /// <summary>Shapiro &amp; Senapathy (1987) consensus-value score on positions -2..+6,
        /// reported in [0, 1] (the published 0-100 score divided by 100). Default.</summary>
        ShapiroSenapathy,

        /// <summary>Yeo &amp; Burge (2004) MaxEntScan <c>score5ss</c> in bits on the 9-nt window
        /// -3..+6 (<see cref="ScoreDonorMaxEnt"/>); <c>minScore</c> is then a threshold in bits.</summary>
        MaxEntScan
    }

    /// <summary>
    /// Scoring model for <see cref="FindAcceptorSites(string, double, bool, AcceptorScoringModel, ShapiroSenapathyCategory)"/>.
    /// </summary>
    public enum AcceptorScoringModel
    {
        /// <summary>Shapiro &amp; Senapathy (1987) 3' splice-site score on positions -14..+1
        /// (pyrimidine region best 8 of 10 + CAGG), reported in [0, 1] (published 0-100 / 100). Default.</summary>
        ShapiroSenapathy,

        /// <summary>Yeo &amp; Burge (2004) MaxEntScan <c>score3ss</c> in bits on the 23-nt window
        /// -20..+3 (<see cref="ScoreAcceptorMaxEnt"/>); <c>minScore</c> is then a threshold in bits.</summary>
        MaxEntScan
    }

    /// <summary>
    /// Branch-point consensus model for <see cref="FindBranchPoints(string, int, int, double, BranchPointModel)"/>.
    /// </summary>
    public enum BranchPointModel
    {
        /// <summary>Human <c>yUnAy</c> (Gao et al. 2008), conservation-weighted. Default.</summary>
        HumanYUnAy,

        /// <summary><i>S. cerevisiae</i> <c>UACUAAC</c> (Langford &amp; Gallwitz 1983), match fraction.</summary>
        YeastUACUAAC
    }

    /// <summary>
    /// GENBANK category subtable of Shapiro &amp; Senapathy (1987) Table 1 used for the
    /// donor (5' splice site) and acceptor (3' splice site) percentages. The paper scores a sequence with the subtable of
    /// its organism; <see cref="Primate"/> is the default (human sequences).
    /// </summary>
    public enum ShapiroSenapathyCategory
    {
        /// <summary>Primate (542 5' ss, 542 3' ss).</summary>
        Primate,
        /// <summary>All seven GENBANK categories pooled (1446 5' ss, 1432 3' ss).</summary>
        Total,
        /// <summary>Rodent (393 5' ss, 378 3' ss).</summary>
        Rodent,
        /// <summary>Other mammals (65 5' ss, 67 3' ss).</summary>
        OtherMammal,
        /// <summary>Other vertebrates (209 5' ss, 202 3' ss).</summary>
        OtherVertebrate,
        /// <summary>Invertebrates (80 5' ss, 86 3' ss).</summary>
        Invertebrate,
        /// <summary>Plants (111 5' ss, 110 3' ss).</summary>
        Plant,
        /// <summary>Viruses (46 5' ss, 47 3' ss).</summary>
        Viral
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

    // --- Shapiro & Senapathy (1987) 5' splice-site (donor) score ---
    // Source: Shapiro MB, Senapathy P (1987) "RNA splice junctions of different classes of
    // eukaryotes: sequence statistics and functional implications in gene expression."
    // Nucleic Acids Res 15(17):7155-7174 (DOI 10.1093/nar/15.17.7155).
    // Table 1 (p. 7159): percentage of each nucleotide at 5' ss positions -3..+6 (exon -3..-1,
    // intron +1..+6; +1/+2 = the GT), one subtable per GENBANK category. Rows below are
    // -3,-2,-1,+1,+2,+3,+4,+5,+6; columns A, C, G, T(=U). Transcribed from the scanned paper.
    // Methods (p. 7157-7158): every subsequence of length 8 (positions -2..+6) is scored as
    //   score = 100 (t - mint) / (maxt - mint)
    // where t is the sum of the percentages of the observed nucleotides and mint / maxt are
    // the sums of the lowest / highest percentages of the eight positions (worked example,
    // primate AAGTGAGT: t = 509, mint = 47, maxt = 595, score 84.3). Position -3 is tabulated
    // but not part of the scored subsequence. This library reports score / 100 in [0, 1].
    private static readonly IReadOnlyDictionary<ShapiroSenapathyCategory, int[,]> ShapiroSenapathyDonorTables =
        new Dictionary<ShapiroSenapathyCategory, int[,]>
        {
            [ShapiroSenapathyCategory.Total] = new int[,]
            {
                { 36, 36, 18, 11 }, { 56, 13, 14, 17 }, { 11, 3, 78, 8 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 60, 3, 35, 3 }, { 69, 9, 12, 11 }, { 7, 6, 82, 5 }, { 14, 14, 19, 53 },
            },
            [ShapiroSenapathyCategory.Primate] = new int[,]
            {
                { 32, 37, 19, 12 }, { 58, 13, 15, 15 }, { 10, 4, 78, 8 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 57, 2, 39, 2 }, { 71, 8, 12, 9 }, { 5, 6, 84, 5 }, { 16, 15, 22, 47 },
            },
            [ShapiroSenapathyCategory.Rodent] = new int[,]
            {
                { 38, 35, 16, 11 }, { 56, 14, 14, 16 }, { 11, 3, 79, 7 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 55, 4, 37, 4 }, { 73, 7, 11, 8 }, { 7, 6, 82, 5 }, { 14, 14, 19, 54 },
            },
            [ShapiroSenapathyCategory.OtherMammal] = new int[,]
            {
                { 38, 43, 11, 8 }, { 52, 8, 22, 18 }, { 8, 0, 86, 6 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 42, 0, 52, 6 }, { 65, 5, 22, 9 }, { 6, 5, 89, 0 }, { 13, 19, 13, 56 },
            },
            [ShapiroSenapathyCategory.OtherVertebrate] = new int[,]
            {
                { 35, 36, 19, 11 }, { 57, 12, 15, 15 }, { 11, 3, 77, 9 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 68, 1, 29, 2 }, { 71, 8, 11, 11 }, { 10, 6, 79, 5 }, { 11, 14, 26, 50 },
            },
            [ShapiroSenapathyCategory.Invertebrate] = new int[,]
            {
                { 38, 28, 25, 10 }, { 63, 13, 10, 15 }, { 8, 5, 80, 8 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 75, 0, 21, 4 }, { 76, 11, 7, 6 }, { 7, 4, 82, 7 }, { 17, 13, 4, 66 },
            },
            [ShapiroSenapathyCategory.Plant] = new int[,]
            {
                { 41, 34, 17, 8 }, { 47, 13, 10, 31 }, { 9, 4, 77, 10 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 64, 7, 19, 10 }, { 35, 19, 10, 36 }, { 10, 6, 72, 12 }, { 10, 17, 8, 64 },
            },
            [ShapiroSenapathyCategory.Viral] = new int[,]
            {
                { 35, 28, 22, 15 }, { 54, 26, 9, 11 }, { 35, 7, 54, 4 }, { 0, 0, 100, 0 }, { 0, 0, 0, 100 },
                { 83, 0, 17, 0 }, { 70, 11, 15, 4 }, { 15, 2, 74, 9 }, { 11, 9, 13, 67 },
            },
        };

    // Table rows used by the S&S 5' score: positions -2..+6 (rows 1..8; row 0 = -3 unused).
    private const int ShapiroSenapathyFirstScoredRow = 1;
    private const int ShapiroSenapathyRowCount = 9;
    // Offset of the scored subsequence start (-2) relative to the G of GT (+1).
    private const int ShapiroSenapathyUpstream = 2;
    // Offset of the last scored nucleotide (+6) relative to the G of GT (+1).
    private const int ShapiroSenapathyDownstream = 5;

    // U12-type 5' splice-site consensus /RTATCCTY at intron positions +1..+8 (R = A/G, Y = C/T):
    // Sharp PA, Burge CB (1997) Cell 91:875-879; Burge CB, Tuschl T, Sharp PA (1999) The RNA
    // World 2nd ed., CSHL Press; Hall SL, Padgett RA (1994) J Mol Biol 239:357 (ATATCCTT, AT-AC);
    // AT-AC and GT-AG termini both occur in U12-type introns (Dietrich, Incorvaia & Padgett 1997
    // Mol Cell 1:151). Stored in RNA alphabet; each entry lists the nucleotides that match.
    private static readonly string[] U12DonorConsensus = { "AG", "U", "A", "U", "C", "C", "U", "CU" };

    // --- Shapiro & Senapathy (1987) 3' splice-site (acceptor) score ---
    // Source: Shapiro & Senapathy (1987) NAR 15:7155, Table 1 (p. 7159): percentage of each
    // nucleotide at 3' ss positions -14..+1 (intron -14..-1, -2/-1 = the AG; +1 = first exon
    // nucleotide), one subtable per GENBANK category. Rows below are -14..-1,+1; columns
    // A, C, G, T(=U). Transcribed from the scanned paper.
    // Methods (p. 7158): left of the NCAGG consensus the 10 pyrimidine locations (-14..-5) are
    // scored by the best 8 of 10 percentages; the CAGG part (-3,-2,-1,+1) by its 4 percentages;
    //   score = 100 ((t1 - l1)/(h1 - l1) + (t2 - l2)/(h2 - l2)) / 2
    // with t1 = sum of the best 8 of the 10 observed percentages at -14..-5, l1 (h1) = sum of the
    // lowest (highest) 8 of the 10 per-position minima (maxima), t2 / l2 / h2 = observed / minimum
    // / maximum percentages summed over -3,-2,-1,+1; position -4 (N) is not scored. "Any sequence
    // containing an AG in those 10 locations [is] eliminated from consideration." Reproduces the
    // published 3' ss scores of Table 5 (40/41 at 1 dp). This library reports score / 100.
    private static readonly IReadOnlyDictionary<ShapiroSenapathyCategory, int[,]> ShapiroSenapathyAcceptorTables =
        new Dictionary<ShapiroSenapathyCategory, int[,]>
        {
            [ShapiroSenapathyCategory.Total] = new int[,]
            {
                { 11, 30, 13, 46 }, { 12, 32, 12, 44 }, { 9, 27, 11, 53 }, { 8, 29, 11, 51 }, { 11, 30, 10, 49 },
                { 11, 32, 11, 46 }, { 10, 38, 9, 44 }, { 10, 38, 8, 44 }, { 8, 38, 7, 47 }, { 8, 35, 7, 49 },
                { 25, 27, 24, 24 }, { 4, 72, 1, 23 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 26, 14, 50, 9 },
            },
            [ShapiroSenapathyCategory.Primate] = new int[,]
            {
                { 9, 31, 15, 45 }, { 9, 33, 13, 45 }, { 7, 31, 11, 51 }, { 7, 35, 7, 51 }, { 10, 35, 7, 47 },
                { 10, 35, 11, 44 }, { 7, 43, 7, 42 }, { 9, 41, 8, 42 }, { 6, 39, 6, 48 }, { 6, 40, 8, 46 },
                { 23, 29, 23, 24 }, { 3, 74, 1, 22 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 28, 13, 49, 10 },
            },
            [ShapiroSenapathyCategory.Rodent] = new int[,]
            {
                { 9, 35, 10, 46 }, { 11, 36, 11, 42 }, { 5, 29, 10, 56 }, { 6, 28, 15, 51 }, { 7, 32, 10, 51 },
                { 9, 32, 11, 48 }, { 7, 40, 6, 47 }, { 6, 41, 6, 47 }, { 9, 43, 7, 42 }, { 9, 39, 6, 46 },
                { 28, 30, 16, 26 }, { 4, 73, 0, 23 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 28, 14, 49, 9 },
            },
            [ShapiroSenapathyCategory.OtherMammal] = new int[,]
            {
                { 8, 29, 3, 60 }, { 5, 40, 14, 41 }, { 5, 29, 11, 56 }, { 0, 27, 10, 63 }, { 10, 43, 6, 40 },
                { 9, 18, 9, 64 }, { 3, 45, 12, 40 }, { 1, 51, 4, 43 }, { 10, 51, 4, 34 }, { 7, 52, 7, 33 },
                { 33, 25, 19, 22 }, { 0, 79, 0, 21 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 13, 30, 51, 6 },
            },
            [ShapiroSenapathyCategory.OtherVertebrate] = new int[,]
            {
                { 9, 35, 13, 43 }, { 10, 33, 10, 47 }, { 8, 24, 10, 59 }, { 8, 29, 9, 53 }, { 9, 28, 10, 53 },
                { 5, 35, 13, 47 }, { 8, 39, 10, 43 }, { 10, 41, 7, 43 }, { 5, 39, 6, 50 }, { 7, 27, 8, 59 },
                { 21, 27, 36, 17 }, { 5, 74, 0, 22 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 21, 17, 53, 8 },
            },
            [ShapiroSenapathyCategory.Invertebrate] = new int[,]
            {
                { 25, 14, 19, 42 }, { 25, 21, 12, 42 }, { 24, 13, 13, 49 }, { 19, 21, 19, 41 }, { 16, 18, 11, 55 },
                { 26, 41, 7, 26 }, { 24, 17, 12, 47 }, { 22, 33, 3, 42 }, { 12, 25, 0, 63 }, { 8, 20, 3, 70 },
                { 16, 16, 26, 42 }, { 1, 76, 0, 22 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 38, 9, 41, 12 },
            },
            [ShapiroSenapathyCategory.Plant] = new int[,]
            {
                { 21, 17, 12, 50 }, { 21, 15, 15, 49 }, { 21, 22, 12, 45 }, { 18, 19, 13, 50 }, { 21, 15, 25, 40 },
                { 20, 21, 15, 44 }, { 27, 15, 15, 44 }, { 17, 11, 25, 47 }, { 14, 21, 18, 47 }, { 14, 19, 11, 56 },
                { 27, 14, 35, 24 }, { 4, 62, 0, 35 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 21, 14, 57, 8 },
            },
            [ShapiroSenapathyCategory.Viral] = new int[,]
            {
                { 2, 26, 30, 41 }, { 17, 30, 15, 38 }, { 21, 28, 15, 36 }, { 15, 15, 15, 55 }, { 21, 17, 13, 49 },
                { 19, 17, 2, 62 }, { 15, 26, 11, 49 }, { 17, 26, 9, 49 }, { 6, 32, 6, 55 }, { 11, 34, 4, 51 },
                { 38, 26, 26, 11 }, { 11, 60, 6, 23 }, { 100, 0, 0, 0 }, { 0, 0, 100, 0 }, { 15, 4, 66, 15 },
            },
        };

    // S&S 3' ss layout: 15 rows -14..-1,+1; rows 0..9 = pyrimidine region -14..-5 (best 8 of 10),
    // row 10 = -4 (N, unscored), rows 11..14 = C,A,G (-3..-1) and +1.
    private const int ShapiroSenapathyAcceptorPyrimidineRows = 10;
    private const int ShapiroSenapathyAcceptorPyrimidineBest = 8;
    private const int ShapiroSenapathyAcceptorCaggFirstRow = 11;
    private const int ShapiroSenapathyAcceptorRowCount = 15;
    // Offset of the A of AG (-2) from the first scored nucleotide (-14).
    private const int ShapiroSenapathyAcceptorUpstream = 12;

    // --- U12-type 3' splice site: branch-site recognition ---
    // U12-type introns lack a polypyrimidine tract; their 3' ss is recognised through the highly
    // conserved branch site, which lies close to the 3' end (Turunen et al. 2013 WIREs RNA 4:61,
    // optimum 11-13 nt; Dietrich et al. 2001 RNA 7:1845). Consensus UUCCUUAAC with the branch
    // adenosine at the second A (Brock, Dietrich & Padgett 2008 RNA 14:2430: NOL1 intron F, fully
    // consensus branch site -UUCCUUAAC-, branch A90). Branch-site-to-acceptor distance window
    // 8-21 nt: Levine & Durbin 2001 NAR 29:4006 (human genome U12 scan). AT-AC (AC terminus) and
    // GT-AG (AG terminus) U12 introns (Dietrich, Incorvaia & Padgett 1997 Mol Cell 1:151).
    // Stored in RNA alphabet; branch A at index 7.
    private const string U12BranchSiteConsensus = "UUCCUUAAC";
    private const int U12BranchSiteBranchIndex = 7;
    private const int U12BranchToAcceptorMinDistance = 8;
    private const int U12BranchToAcceptorMaxDistance = 21;

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

    // Yeast (S. cerevisiae) branch-point sequence UACUAAC, strictly conserved, branch A = the
    // last A (index 5) (Langford & Gallwitz 1983 Cell 33:519; Gao et al. 2008 Introduction).
    private const string YeastBranchPointConsensus = "UACUAAC";
    private const int YeastBranchPointBranchIndex = 5;

    // Human yUnAy motif (positions -3..+1): length and index of the branch A.
    private const int HumanBranchPointMotifLength = 5;
    private const int HumanBranchPointBranchIndex = 3;

    #endregion

    #region Splice Site Prediction

    /// <summary>
    /// Finds all potential donor (5') splice sites in a sequence, scored with the
    /// Shapiro &amp; Senapathy (1987) primate 5' splice-site table.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="FindDonorSites(string, double, bool, DonorScoringModel, ShapiroSenapathyCategory)"/>
    /// with <see cref="DonorScoringModel.ShapiroSenapathy"/> and <see cref="ShapiroSenapathyCategory.Primate"/>.
    /// </remarks>
    public static IEnumerable<SpliceSite> FindDonorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false) =>
        FindDonorSites(sequence, minScore, includeNonCanonical,
            DonorScoringModel.ShapiroSenapathy, ShapiroSenapathyCategory.Primate);

    /// <summary>
    /// Finds all potential donor (5') splice sites in a sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The input is upper-cased and read as RNA (T = U). A candidate is a <c>GU</c> (and, when
    /// <paramref name="includeNonCanonical"/> is true, a <c>GC</c>) dinucleotide; its
    /// <see cref="SpliceSite.Position"/> is the 0-based index of the <c>G</c>, i.e. the first
    /// intron nucleotide (+1); the exon-intron boundary lies immediately before it.
    /// <see cref="SpliceSite.Motif"/> is the -3..+6 context (3 exon + 6 intron nt, truncated at
    /// the sequence ends).
    /// </para>
    /// <para>
    /// <b><see cref="DonorScoringModel.ShapiroSenapathy"/> (default):</b> the Shapiro &amp;
    /// Senapathy (1987) score of the 8-nt subsequence -2..+6,
    /// <c>(t - mint) / (maxt - mint)</c> with <c>t</c> the sum of the Table 1 percentages of the
    /// observed nucleotides (subtable <paramref name="category"/>), i.e. the published 0-100
    /// score / 100 (e.g. primate <c>AAGTGAGT</c> → 0.843). A candidate whose -2..+6 window runs
    /// off the sequence or contains a non-A/C/G/T/U character is not scored (no site).
    /// <see cref="SpliceSite.Confidence"/> = clamp((score - 0.5) / 0.5, 0, 1).
    /// </para>
    /// <para>
    /// <b><see cref="DonorScoringModel.MaxEntScan"/>:</b> the Yeo &amp; Burge (2004) MaxEntScan
    /// <c>score5ss</c> (<see cref="ScoreDonorMaxEnt"/>) of the 9-nt window -3..+6, in bits;
    /// <paramref name="minScore"/> is a threshold in bits. A window that runs off the sequence
    /// or contains a non-A/C/G/T/U character is skipped. <see cref="SpliceSite.Confidence"/> =
    /// 2^s / (1 + 2^s), the posterior probability of a true site at prior odds 1 (score5 is the
    /// log2 likelihood ratio site/background).
    /// </para>
    /// <para>
    /// <b>U12-type donors</b> (<paramref name="includeNonCanonical"/>, Shapiro-Senapathy model
    /// only): every <c>AU</c> or <c>GU</c> start is also scored against the U12 5' splice-site
    /// consensus <c>/RUAUCCUY</c> (+1..+8; Sharp &amp; Burge 1997, Burge et al. 1999) as the
    /// fraction of matching positions, and reported as <see cref="SpliceSiteType.U12Donor"/>
    /// when the score reaches <paramref name="minScore"/>. A <c>GU</c> site can therefore be
    /// reported both as a U2 <see cref="SpliceSiteType.Donor"/> and as a U12 candidate.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA or RNA sequence (case-insensitive); null/empty yields no sites.</param>
    /// <param name="minScore">Minimum score to report a site (model units).</param>
    /// <param name="includeNonCanonical">Also report GC donors and U12-type donors.</param>
    /// <param name="model">Scoring model for U2-type (GU/GC) donors.</param>
    /// <param name="category">Shapiro &amp; Senapathy (1987) Table 1 subtable (ignored by MaxEntScan).</param>
    public static IEnumerable<SpliceSite> FindDonorSites(
        string sequence,
        double minScore,
        bool includeNonCanonical,
        DonorScoringModel model,
        ShapiroSenapathyCategory category = ShapiroSenapathyCategory.Primate)
    {
        if (!Enum.IsDefined(model))
            throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown donor scoring model.");
        if (!ShapiroSenapathyDonorTables.ContainsKey(category))
            throw new ArgumentOutOfRangeException(nameof(category), category, "Unknown Shapiro & Senapathy category.");

        return FindDonorSitesIterator(sequence, minScore, includeNonCanonical, model, category);
    }

    private static IEnumerable<SpliceSite> FindDonorSitesIterator(
        string sequence,
        double minScore,
        bool includeNonCanonical,
        DonorScoringModel model,
        ShapiroSenapathyCategory category)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 0; i + 1 < upper.Length; i++)
        {
            char first = upper[i];
            char second = upper[i + 1];

            // U2-type candidate: canonical GU, or (opt-in) non-canonical GC.
            if (first == 'G' && (second == 'U' || (includeNonCanonical && second == 'C')))
            {
                double? score = model == DonorScoringModel.MaxEntScan
                    ? ScoreDonorSiteMaxEnt(upper, i)
                    : ScoreDonorSiteShapiroSenapathy(upper, i, category);

                if (score is double s && s >= minScore)
                {
                    double confidence = model == DonorScoringModel.MaxEntScan
                        ? 1.0 / (1.0 + Math.Pow(2.0, -s))
                        : CalculateConfidence(s, 0.5, 1.0);

                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: s,
                        Confidence: confidence);
                }
            }

            // U12-type candidate: /RU (AU or GU) start, scored against /RUAUCCUY.
            if (includeNonCanonical && model == DonorScoringModel.ShapiroSenapathy &&
                (first == 'A' || first == 'G') && second == 'U')
            {
                double? score = ScoreU12DonorSite(upper, i);
                if (score is double s && s >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.U12Donor,
                        Motif: GetMotifContext(upper, i, 3, 8),
                        Score: s,
                        Confidence: CalculateConfidence(s, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Finds all potential acceptor (3') splice sites in a sequence, scored with the
    /// Shapiro &amp; Senapathy (1987) primate 3' splice-site table.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="FindAcceptorSites(string, double, bool, AcceptorScoringModel, ShapiroSenapathyCategory)"/>
    /// with <see cref="AcceptorScoringModel.ShapiroSenapathy"/> and <see cref="ShapiroSenapathyCategory.Primate"/>.
    /// </remarks>
    public static IEnumerable<SpliceSite> FindAcceptorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false) =>
        FindAcceptorSites(sequence, minScore, includeNonCanonical,
            AcceptorScoringModel.ShapiroSenapathy, ShapiroSenapathyCategory.Primate);

    /// <summary>
    /// Finds all potential acceptor (3') splice sites in a sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The input is upper-cased and read as RNA (T = U). A candidate is an <c>AG</c> dinucleotide
    /// (intron positions -2/-1); its <see cref="SpliceSite.Position"/> is the 0-based index of the
    /// <c>G</c>, i.e. the last intron nucleotide; the intron-exon boundary lies immediately after it.
    /// <see cref="SpliceSite.Motif"/> is the -20..+3 context (20 intron + 3 exon nt, the MaxEntScan
    /// score3 window, truncated at the sequence ends).
    /// </para>
    /// <para>
    /// <b><see cref="AcceptorScoringModel.ShapiroSenapathy"/> (default):</b> the Shapiro &amp;
    /// Senapathy (1987) 3' ss score of the 15-nt subsequence -14..+1 (Table 1 subtable
    /// <paramref name="category"/>): <c>((t1 - l1)/(h1 - l1) + (t2 - l2)/(h2 - l2)) / 2</c>, where
    /// the pyrimidine region -14..-5 is scored by the best 8 of its 10 percentages and the
    /// <c>CAGG</c> part by the percentages at -3, -2, -1, +1 (position -4 unscored); i.e. the
    /// published 0-100 score / 100 (e.g. primate <c>TTATTCTACTCCAG|G</c> → 0.958). A candidate with
    /// an <c>AG</c> inside -14..-5 is eliminated (as in the paper), and one whose -14..+1 window runs
    /// off the sequence or contains a non-A/C/G/T/U character is not scored (no site).
    /// <see cref="SpliceSite.Confidence"/> = clamp((score - 0.5) / 0.5, 0, 1).
    /// </para>
    /// <para>
    /// <b><see cref="AcceptorScoringModel.MaxEntScan"/>:</b> the Yeo &amp; Burge (2004) MaxEntScan
    /// <c>score3ss</c> (<see cref="ScoreAcceptorMaxEnt"/>) of the 23-nt window -20..+3, in bits;
    /// <paramref name="minScore"/> is a threshold in bits. A window that runs off the sequence or
    /// contains a non-A/C/G/T/U character is skipped. <see cref="SpliceSite.Confidence"/> =
    /// 2^s / (1 + 2^s).
    /// </para>
    /// <para>
    /// <b>U12-type acceptors</b> (<paramref name="includeNonCanonical"/>, Shapiro-Senapathy model
    /// only): U12-type introns have no polypyrimidine tract and are recognised through their branch
    /// site. Every <c>AC</c> (AT-AC) or <c>AG</c> (GT-AG U12) terminus is scored as the best
    /// fraction of matching positions to the U12 branch-site consensus <c>UUCCUUAAC</c> (branch A
    /// = second A; Brock et al. 2008) over all placements whose branch A lies 8-21 nt upstream of
    /// the last intron nucleotide (Levine &amp; Durbin 2001), and reported as
    /// <see cref="SpliceSiteType.U12Acceptor"/> when it reaches <paramref name="minScore"/>.
    /// Placements running off the sequence or containing a non-A/C/G/T/U character are skipped.
    /// An <c>AG</c> can therefore be reported both as a U2 <see cref="SpliceSiteType.Acceptor"/> and
    /// as a U12 candidate.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA or RNA sequence (case-insensitive); null/empty yields no sites.</param>
    /// <param name="minScore">Minimum score to report a site (model units).</param>
    /// <param name="includeNonCanonical">Also report U12-type (AC or AG) acceptors.</param>
    /// <param name="model">Scoring model for U2-type (AG) acceptors.</param>
    /// <param name="category">Shapiro &amp; Senapathy (1987) Table 1 subtable (ignored by MaxEntScan).</param>
    public static IEnumerable<SpliceSite> FindAcceptorSites(
        string sequence,
        double minScore,
        bool includeNonCanonical,
        AcceptorScoringModel model,
        ShapiroSenapathyCategory category = ShapiroSenapathyCategory.Primate)
    {
        if (!Enum.IsDefined(model))
            throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown acceptor scoring model.");
        if (!ShapiroSenapathyAcceptorTables.ContainsKey(category))
            throw new ArgumentOutOfRangeException(nameof(category), category, "Unknown Shapiro & Senapathy category.");

        return FindAcceptorSitesIterator(sequence, minScore, includeNonCanonical, model, category);
    }

    private static IEnumerable<SpliceSite> FindAcceptorSitesIterator(
        string sequence,
        double minScore,
        bool includeNonCanonical,
        AcceptorScoringModel model,
        ShapiroSenapathyCategory category)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 0; i + 1 < upper.Length; i++)
        {
            if (upper[i] != 'A')
                continue;
            char second = upper[i + 1];

            // U2-type candidate: AG (A at -2, G at -1).
            if (second == 'G')
            {
                double? score = model == AcceptorScoringModel.MaxEntScan
                    ? ScoreAcceptorSiteMaxEnt(upper, i)
                    : ScoreAcceptorSiteShapiroSenapathy(upper, i, category);

                if (score is double s && s >= minScore)
                {
                    double confidence = model == AcceptorScoringModel.MaxEntScan
                        ? 1.0 / (1.0 + Math.Pow(2.0, -s))
                        : CalculateConfidence(s, 0.5, 1.0);

                    yield return new SpliceSite(
                        Position: i + 1, // G of AG (last intron nucleotide)
                        Type: SpliceSiteType.Acceptor,
                        Motif: GetMotifContext(upper, i, 18, 5),
                        Score: s,
                        Confidence: confidence);
                }
            }

            // U12-type candidate: AC (AT-AC) or AG (GT-AG) terminus, scored by its branch site.
            if (includeNonCanonical && model == AcceptorScoringModel.ShapiroSenapathy &&
                (second == 'C' || second == 'G'))
            {
                double? score = ScoreU12AcceptorSite(upper, i + 1);
                if (score is double s && s >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i + 1,
                        Type: SpliceSiteType.U12Acceptor,
                        Motif: GetMotifContext(upper, i, 18, 5),
                        Score: s,
                        Confidence: CalculateConfidence(s, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Detects the explicit branch point upstream of a 3' acceptor AG using the human
    /// <c>yUnAy</c> branch-point consensus (Gao et al. 2008). This is an opt-in,
    /// additive analysis; the <see cref="FindAcceptorSites(string, double, bool)"/> score
    /// (Shapiro &amp; Senapathy 1987) has no branch-point term.
    /// </summary>
    /// <remarks>
    /// The branch-point adenosine sits at motif position 0 of the consensus
    /// <c>yUnAy</c> (positions -3..+1). It is searched in the
    /// <see cref="BranchPointMinDistanceFromAg"/>..<see cref="BranchPointMaxDistanceFromAg"/>
    /// (default 18-40) nt window upstream of the AG, bracketing the -34..-21 core
    /// reported by Gao et al. (2008). Each candidate is scored by a conservation-weighted
    /// fraction over the four informative positions (-3 pyrimidine, -2 U, 0 A, +1
    /// pyrimidine); position -1 is uninformative ('n'). A placement whose 5-nt motif contains a
    /// non-A/C/G/T/U character is not scored. On ties the candidate nearest the AG wins (Gao et al.
    /// 2008: downstream branch sites are used more often). The pyrimidine fraction of the
    /// polypyrimidine tract (branch point +4..+24, Gao et al. 2008, ending before the AG) is
    /// reported alongside the best candidate.
    /// </remarks>
    /// <param name="sequence">The pre-mRNA / intron sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="acceptorAgPosition">0-based index of the terminal <c>G</c> of the 3'
    /// acceptor AG (i.e. the <see cref="SpliceSite.Position"/> reported by
    /// <see cref="FindAcceptorSites(string, double, bool)"/>).</param>
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
            if (branchA - 3 < 0 || branchA + 1 >= upper.Length ||
                !IsAcgu(upper, branchA - 3, HumanBranchPointMotifLength))
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
    /// Finds branch-point candidates in a sequence with the human <c>yUnAy</c> consensus
    /// (Gao et al. 2008).
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="FindBranchPoints(string, int, int, double, BranchPointModel)"/> with
    /// <see cref="BranchPointModel.HumanYUnAy"/>.
    /// </remarks>
    public static IEnumerable<SpliceSite> FindBranchPoints(
        string sequence,
        int searchStart = 0,
        int searchEnd = -1,
        double minScore = 0.5) =>
        FindBranchPoints(sequence, searchStart, searchEnd, minScore, BranchPointModel.HumanYUnAy);

    /// <summary>
    /// Finds branch-point candidates in a sequence.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Every motif placement whose first nucleotide index lies in
    /// [<paramref name="searchStart"/>, <paramref name="searchEnd"/>] (inclusive; -1 = to the end)
    /// and which fits inside the sequence is scored; placements containing a non-A/C/G/T/U
    /// character are skipped. <see cref="SpliceSite.Position"/> is the 0-based index of the branch
    /// adenosine position, <see cref="SpliceSite.Motif"/> the scored motif (RNA alphabet).
    /// </para>
    /// <para>
    /// <b><see cref="BranchPointModel.HumanYUnAy"/> (default):</b> 5-nt motif <c>yUnAy</c>
    /// (-3..+1, branch A at 0), scored as the conservation-weighted match fraction of Gao et al.
    /// (2008) Table 4 (y@-3 0.790, U@-2 0.746, A@0 0.923, y@+1 0.751; -1 uninformative); 1.0 for a
    /// perfect <c>yUnAy</c>.
    /// </para>
    /// <para>
    /// <b><see cref="BranchPointModel.YeastUACUAAC"/>:</b> 7-nt motif, fraction of positions matching
    /// the strictly conserved <i>S. cerevisiae</i> branch-point sequence <c>UACUAAC</c> (branch A =
    /// sixth nucleotide; Langford &amp; Gallwitz 1983).
    /// </para>
    /// <see cref="SpliceSite.Confidence"/> = clamp((score - 0.5) / 0.5, 0, 1).
    /// </remarks>
    public static IEnumerable<SpliceSite> FindBranchPoints(
        string sequence,
        int searchStart,
        int searchEnd,
        double minScore,
        BranchPointModel model)
    {
        if (!Enum.IsDefined(model))
            throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown branch-point model.");

        return FindBranchPointsIterator(sequence, searchStart, searchEnd, minScore, model);
    }

    private static IEnumerable<SpliceSite> FindBranchPointsIterator(
        string sequence,
        int searchStart,
        int searchEnd,
        double minScore,
        BranchPointModel model)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        bool yeast = model == BranchPointModel.YeastUACUAAC;
        int motifLength = yeast ? YeastBranchPointConsensus.Length : HumanBranchPointMotifLength;
        int branchIndex = yeast ? YeastBranchPointBranchIndex : HumanBranchPointBranchIndex;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        int lastStart = upper.Length - motifLength;
        int end = searchEnd < 0 ? lastStart : Math.Min(searchEnd, lastStart);

        for (int i = Math.Max(0, searchStart); i <= end; i++)
        {
            if (!IsAcgu(upper, i, motifLength))
                continue;

            double score = yeast
                ? ConsensusMatchFraction(upper, i, YeastBranchPointConsensus)
                : ScoreBranchPointConsensus(upper, i + branchIndex);

            if (score >= minScore)
            {
                yield return new SpliceSite(
                    Position: i + branchIndex,
                    Type: SpliceSiteType.Branch,
                    Motif: upper.Substring(i, motifLength),
                    Score: score,
                    Confidence: CalculateConfidence(score, 0.5, 1.0));
            }
        }
    }

    #endregion

    #region Scoring Functions

    /// <summary>
    /// Shapiro &amp; Senapathy (1987) 5' splice-site score in [0, 1] of the 8-nt subsequence
    /// -2..+6 around the donor whose +1 nucleotide is at <paramref name="position"/>;
    /// null when the subsequence runs off the sequence or contains a non-A/C/G/U character.
    /// <paramref name="sequence"/> must be upper-case RNA.
    /// </summary>
    private static double? ScoreDonorSiteShapiroSenapathy(
        string sequence, int position, ShapiroSenapathyCategory category)
    {
        int start = position - ShapiroSenapathyUpstream;
        int end = position + ShapiroSenapathyDownstream;
        if (start < 0 || end >= sequence.Length)
            return null;

        int[,] table = ShapiroSenapathyDonorTables[category];
        int t = 0, minT = 0, maxT = 0;
        for (int row = ShapiroSenapathyFirstScoredRow; row < ShapiroSenapathyRowCount; row++)
        {
            int column = NucleotideColumn(sequence[start + row - ShapiroSenapathyFirstScoredRow]);
            if (column < 0)
                return null;

            t += table[row, column];
            int lo = int.MaxValue, hi = int.MinValue;
            for (int c = 0; c < 4; c++)
            {
                lo = Math.Min(lo, table[row, c]);
                hi = Math.Max(hi, table[row, c]);
            }
            minT += lo;
            maxT += hi;
        }

        return (double)(t - minT) / (maxT - minT);
    }

    /// <summary>
    /// MaxEntScan score5 (bits) of the 9-nt window -3..+6 around the donor whose +1 nucleotide
    /// is at <paramref name="position"/>; null when the window runs off the sequence or
    /// contains a non-A/C/G/U character. <paramref name="sequence"/> must be upper-case RNA.
    /// </summary>
    private static double? ScoreDonorSiteMaxEnt(string sequence, int position)
    {
        int start = position - MaxEntDonorGtStart;
        if (start < 0 || start + MaxEntDonorWindowLength > sequence.Length)
            return null;

        for (int k = start; k < start + MaxEntDonorWindowLength; k++)
        {
            if (NucleotideColumn(sequence[k]) < 0)
                return null;
        }

        return ScoreDonorMaxEnt(sequence.Substring(start, MaxEntDonorWindowLength));
    }

    /// <summary>Column of an upper-case RNA/DNA nucleotide in an A,C,G,T(U) table; -1 otherwise.</summary>
    private static int NucleotideColumn(char c) => c switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'U' or 'T' => 3,
        _ => -1
    };

    /// <summary>
    /// Shapiro &amp; Senapathy (1987) 3' splice-site score in [0, 1] of the 15-nt subsequence
    /// -14..+1 around the acceptor whose A of AG (-2) is at <paramref name="agIndex"/>; null when the
    /// subsequence runs off the sequence, contains a non-A/C/G/U character, or has an AG inside the
    /// pyrimidine region -14..-5 (eliminated in the paper). <paramref name="sequence"/> must be
    /// upper-case RNA.
    /// </summary>
    private static double? ScoreAcceptorSiteShapiroSenapathy(
        string sequence, int agIndex, ShapiroSenapathyCategory category)
    {
        int start = agIndex - ShapiroSenapathyAcceptorUpstream;
        if (start < 0 || start + ShapiroSenapathyAcceptorRowCount > sequence.Length ||
            !IsAcgu(sequence, start, ShapiroSenapathyAcceptorRowCount))
            return null;

        for (int k = start; k < start + ShapiroSenapathyAcceptorPyrimidineRows - 1; k++)
        {
            if (sequence[k] == 'A' && sequence[k + 1] == 'G')
                return null;
        }

        int[,] table = ShapiroSenapathyAcceptorTables[category];

        // Pyrimidine region: best 8 of the 10 observed percentages; bounds from the 8 lowest
        // per-position minima / 8 highest per-position maxima.
        var observed = new int[ShapiroSenapathyAcceptorPyrimidineRows];
        var minima = new int[ShapiroSenapathyAcceptorPyrimidineRows];
        var maxima = new int[ShapiroSenapathyAcceptorPyrimidineRows];
        for (int row = 0; row < ShapiroSenapathyAcceptorPyrimidineRows; row++)
        {
            observed[row] = table[row, NucleotideColumn(sequence[start + row])];
            (minima[row], maxima[row]) = RowRange(table, row);
        }

        int t1 = observed.OrderByDescending(v => v).Take(ShapiroSenapathyAcceptorPyrimidineBest).Sum();
        int l1 = minima.OrderBy(v => v).Take(ShapiroSenapathyAcceptorPyrimidineBest).Sum();
        int h1 = maxima.OrderByDescending(v => v).Take(ShapiroSenapathyAcceptorPyrimidineBest).Sum();

        // CAGG part: -3, -2, -1, +1.
        int t2 = 0, l2 = 0, h2 = 0;
        for (int row = ShapiroSenapathyAcceptorCaggFirstRow; row < ShapiroSenapathyAcceptorRowCount; row++)
        {
            t2 += table[row, NucleotideColumn(sequence[start + row])];
            var (lo, hi) = RowRange(table, row);
            l2 += lo;
            h2 += hi;
        }

        return ((double)(t1 - l1) / (h1 - l1) + (double)(t2 - l2) / (h2 - l2)) / 2;
    }

    /// <summary>Lowest and highest percentage of a Shapiro &amp; Senapathy table row.</summary>
    private static (int Min, int Max) RowRange(int[,] table, int row)
    {
        int lo = int.MaxValue, hi = int.MinValue;
        for (int c = 0; c < 4; c++)
        {
            lo = Math.Min(lo, table[row, c]);
            hi = Math.Max(hi, table[row, c]);
        }
        return (lo, hi);
    }

    /// <summary>
    /// MaxEntScan score3 (bits) of the 23-nt window -20..+3 around the acceptor whose A of AG is at
    /// <paramref name="agIndex"/>; null when the window runs off the sequence or contains a
    /// non-A/C/G/U character. <paramref name="sequence"/> must be upper-case RNA.
    /// </summary>
    private static double? ScoreAcceptorSiteMaxEnt(string sequence, int agIndex)
    {
        int start = agIndex - MaxEntAcceptorAgStart;
        if (start < 0 || start + MaxEntAcceptorWindowLength > sequence.Length ||
            !IsAcgu(sequence, start, MaxEntAcceptorWindowLength))
            return null;

        return ScoreAcceptorMaxEnt(sequence.Substring(start, MaxEntAcceptorWindowLength));
    }

    /// <summary>True when <paramref name="length"/> characters from <paramref name="start"/> are
    /// all A/C/G/T/U (upper case).</summary>
    private static bool IsAcgu(string sequence, int start, int length)
    {
        for (int k = start; k < start + length; k++)
        {
            if (NucleotideColumn(sequence[k]) < 0)
                return false;
        }
        return true;
    }

    /// <summary>Fraction of positions of <paramref name="consensus"/> matched at
    /// <paramref name="start"/> (caller guarantees bounds).</summary>
    private static double ConsensusMatchFraction(string sequence, int start, string consensus)
    {
        int matches = 0;
        for (int k = 0; k < consensus.Length; k++)
        {
            if (sequence[start + k] == consensus[k])
                matches++;
        }
        return (double)matches / consensus.Length;
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

    /// <summary>
    /// Fraction of the 8 intron positions +1..+8 matching the U12-type 5' splice-site
    /// consensus /RUAUCCUY (Sharp &amp; Burge 1997); null when the 8-nt window runs off the
    /// sequence. <paramref name="sequence"/> must be upper-case RNA.
    /// </summary>
    private static double? ScoreU12DonorSite(string sequence, int position)
    {
        if (position < 0 || position + U12DonorConsensus.Length > sequence.Length)
            return null;

        int matches = 0;
        for (int k = 0; k < U12DonorConsensus.Length; k++)
        {
            if (U12DonorConsensus[k].Contains(sequence[position + k]))
                matches++;
        }

        return (double)matches / U12DonorConsensus.Length;
    }

    /// <summary>
    /// U12-type acceptor score: best fraction of positions matching the U12 branch-site consensus
    /// <c>UUCCUUAAC</c> over the placements whose branch A lies
    /// <see cref="U12BranchToAcceptorMinDistance"/>..<see cref="U12BranchToAcceptorMaxDistance"/> nt
    /// upstream of <paramref name="lastIntronIndex"/> (the terminal C/G of the intron); placements
    /// off the sequence or with a non-A/C/G/U character are skipped; null when none is scorable.
    /// </summary>
    private static double? ScoreU12AcceptorSite(string sequence, int lastIntronIndex)
    {
        double? best = null;
        for (int d = U12BranchToAcceptorMinDistance; d <= U12BranchToAcceptorMaxDistance; d++)
        {
            int start = lastIntronIndex - d - U12BranchSiteBranchIndex;
            if (start < 0 || start + U12BranchSiteConsensus.Length > sequence.Length ||
                !IsAcgu(sequence, start, U12BranchSiteConsensus.Length))
                continue;

            double score = ConsensusMatchFraction(sequence, start, U12BranchSiteConsensus);
            if (best is null || score > best)
                best = score;
        }
        return best;
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
    /// Calculates a splice-site score for a motif of the given type.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Donor:</b> the Yeo &amp; Burge (2004) MaxEntScan <c>score5ss</c> maximum-entropy score in bits,
    /// computed by <see cref="ScoreDonorMaxEnt"/> (the motif must be the 9-nt window: 3 exonic + 6
    /// intronic nt, GT at 1-based positions 4-5; e.g. <c>CAGGTAAGT</c> → 10.86, as printed by
    /// MaxEntScan <c>score5.pl</c>).
    /// </para>
    /// <para>
    /// <b>Acceptor:</b> the Yeo &amp; Burge (2004) MaxEntScan <c>score3ss</c> maximum-entropy score in bits,
    /// computed by <see cref="ScoreAcceptorMaxEnt"/> (the motif must be the 23-nt window: 20 intronic +
    /// 3 exonic nt, AG at 1-based positions 19-20; e.g. <c>TTCCAAACGAACTTTTGTAGGGA</c> → 2.89, as
    /// printed by MaxEntScan <c>score3.pl</c>).
    /// </para>
    /// <para>
    /// Any other <paramref name="type"/> returns 0; an empty motif returns 0.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentException">Donor motif not exactly 9 nt, acceptor motif not exactly 23 nt,
    /// or a motif containing a non-A/C/G/T(/U) character (MaxEntScan cannot score it).</exception>
    public static double CalculateMaxEntScore(string motif, SpliceSiteType type)
    {
        if (string.IsNullOrEmpty(motif))
            return 0;

        // Donor: the real MaxEntScan score5 model (Yeo & Burge 2004; score5.pl / me2x5).
        if (type == SpliceSiteType.Donor)
            return ScoreDonorMaxEnt(motif);

        // Acceptor: the real MaxEntScan score3 model (Yeo & Burge 2004; score3.pl / me2x3acc1-9).
        if (type == SpliceSiteType.Acceptor)
            return ScoreAcceptorMaxEnt(motif);

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

    #region MaxEntScan score3ss (maximum-entropy 3' acceptor model)

    // --- MaxEntScan score3ss (Yeo & Burge 2004) ---
    // Maximum-entropy 3' splice-site (acceptor) score for a 23-nt window:
    // 20 intronic + 3 exonic nucleotides, with the invariant AG at window positions
    // 18..19 (0-based). The score is log2( P_maxent(seq) / P_background(seq) ), and is
    // computed by removing the conserved AG dinucleotide and factorising the remaining
    // 21 ("rest") positions over nine overlapping sub-sequences whose precomputed
    // maximum-entropy probabilities are stored in the embedded table.
    //
    // Sources:
    //   Yeo G, Burge CB (2004). "Maximum entropy modeling of short sequence motifs with
    //   applications to RNA splicing signals." J Comput Biol 11(2-3):377-394.
    //   DOI 10.1089/1066527041410418.
    //   Reference factorisation + tables (retrieved 2026-06-24): the MIT-licensed maxentpy
    //   port (kepbod/maxentpy) score3 / score3_matrix.txt. Provenance + licence are recorded
    //   in src/.../Seqeron.Genomics.Annotation/Data/maxent_score3.LICENSE.md.

    // Length of the MaxEntScan 3' (acceptor) scoring window: 20 intron + 3 exon nt.
    private const int MaxEntAcceptorWindowLength = 23;

    // 0-based positions of the conserved AG dinucleotide within the 23-nt window
    // (window[18]=A, window[19]=G). These two positions are scored by the consensus /
    // background dinucleotide model and removed from the "rest" sub-sequence.
    private const int MaxEntAcceptorAgStart = 18;

    // Background single-nucleotide frequencies P_bgd used by MaxEntScan score3
    // (Yeo & Burge 2004; maxentpy bgd_3). Keyed by A/C/G/T (T==U).
    private static readonly IReadOnlyDictionary<char, double> MaxEntBackground = new Dictionary<char, double>
    {
        { 'A', 0.27 }, { 'C', 0.23 }, { 'G', 0.23 }, { 'T', 0.27 }
    };

    // Consensus probabilities for the first AG position (the 'A'); maxentpy cons1_3.
    private static readonly IReadOnlyDictionary<char, double> MaxEntConsensusAg1 = new Dictionary<char, double>
    {
        { 'A', 0.9903 }, { 'C', 0.0032 }, { 'G', 0.0034 }, { 'T', 0.0030 }
    };

    // Consensus probabilities for the second AG position (the 'G'); maxentpy cons2_3.
    private static readonly IReadOnlyDictionary<char, double> MaxEntConsensusAg2 = new Dictionary<char, double>
    {
        { 'A', 0.0027 }, { 'C', 0.0037 }, { 'G', 0.9905 }, { 'T', 0.0030 }
    };

    // The score3 factorisation over the 21-nt "rest" sequence (window with AG removed).
    // Each tuple is (matrixIndex, restStart, length); the maxent probability of the rest
    // sequence is the product of the lookups for the five "numerator" sub-sequences
    // divided by the product of the four "denominator" sub-sequences (inclusion-exclusion
    // over the overlapping windows). Source: maxentpy score3.
    private static readonly (int MatrixIndex, int RestStart, int Length)[] MaxEntNumeratorFactors =
    {
        (0, 0, 7),
        (1, 7, 7),
        (2, 14, 7),
        (3, 4, 7),
        (4, 11, 7)
    };

    private static readonly (int MatrixIndex, int RestStart, int Length)[] MaxEntDenominatorFactors =
    {
        (5, 4, 3),
        (6, 7, 4),
        (7, 11, 3),
        (8, 14, 4)
    };

    private const string MaxEntScore3ResourceName = "Seqeron.Genomics.Annotation.Data.maxent_score3.txt";

    // Lazily-loaded probability tables: tables[matrixIndex][hash] = probability.
    // Nine matrices (4^k entries each). Loaded once from the embedded resource.
    private static IReadOnlyList<IReadOnlyDictionary<int, double>>? _maxEntScore3Tables;
    private static readonly object MaxEntLoadLock = new();

    private static IReadOnlyList<IReadOnlyDictionary<int, double>> MaxEntScore3Tables
    {
        get
        {
            if (_maxEntScore3Tables is not null)
                return _maxEntScore3Tables;

            lock (MaxEntLoadLock)
            {
                if (_maxEntScore3Tables is not null)
                    return _maxEntScore3Tables;
                _maxEntScore3Tables = LoadMaxEntScore3Tables();
                return _maxEntScore3Tables;
            }
        }
    }

    private static IReadOnlyList<IReadOnlyDictionary<int, double>> LoadMaxEntScore3Tables()
    {
        var assembly = typeof(SpliceSitePredictor).Assembly;
        using Stream? stream = assembly.GetManifestResourceStream(MaxEntScore3ResourceName)
            ?? throw new InvalidOperationException(
                $"Embedded MaxEntScan score3 table '{MaxEntScore3ResourceName}' was not found.");

        var tables = new Dictionary<int, double>[9];
        for (int i = 0; i < tables.Length; i++)
            tables[i] = new Dictionary<int, double>();

        using var reader = new StreamReader(stream);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;

            // Each record: "<matrixIndex>\t<hash>\t<probability>" (whitespace-separated).
            string[] parts = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length != 3)
                throw new InvalidDataException($"Malformed MaxEntScan score3 record: '{line}'.");

            int matrixIndex = int.Parse(parts[0], CultureInfo.InvariantCulture);
            int hash = int.Parse(parts[1], CultureInfo.InvariantCulture);
            double probability = double.Parse(parts[2], CultureInfo.InvariantCulture);
            tables[matrixIndex][hash] = probability;
        }

        return tables;
    }

    // Quaternary (base-4) encoding of an A/C/G/T sub-sequence: A=0, C=1, G=2, T=3,
    // most-significant nucleotide first. Matches maxentpy hashseq. T and U are equivalent.
    private static int HashMaxEntSubsequence(ReadOnlySpan<char> seq)
    {
        int hash = 0;
        for (int i = 0; i < seq.Length; i++)
        {
            int code = seq[i] switch
            {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' or 'U' => 3,
                _ => throw new ArgumentException(
                    $"MaxEntScan score3 requires an A/C/G/T(/U) sequence; found '{seq[i]}'.")
            };
            hash = (hash << 2) | code; // hash * 4 + code
        }

        return hash;
    }

    /// <summary>
    /// Computes the Yeo &amp; Burge (2004) MaxEntScan <c>score3ss</c> maximum-entropy score
    /// (in bits) for a 23-nt 3' splice-site (acceptor) window: 20 intronic + 3 exonic
    /// nucleotides, with the conserved <c>AG</c> at window positions 18-19 (0-based).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Also used by <see cref="FindAcceptorSites(string, double, bool, AcceptorScoringModel, ShapiroSenapathyCategory)"/>
    /// with <see cref="AcceptorScoringModel.MaxEntScan"/> (default there: Shapiro &amp; Senapathy 1987).
    /// The returned value is the log<sub>2</sub> ratio of the maximum-entropy
    /// probability of the window to its background probability — higher means a stronger 3'
    /// splice site (the canonical documented example scores 2.89 bits).
    /// </para>
    /// <para>
    /// The score removes the conserved <c>AG</c> dinucleotide (scored by a consensus /
    /// background model) and factorises the remaining 21 positions over nine overlapping
    /// sub-sequences, multiplying the five numerator and dividing by the four denominator
    /// maximum-entropy probabilities (inclusion-exclusion), then takes log<sub>2</sub> of the
    /// product with the AG term. Probability tables are embedded (see
    /// <c>Data/maxent_score3.LICENSE.md</c> for provenance and licence).
    /// </para>
    /// </remarks>
    /// <param name="window">The 23-nt acceptor window (DNA or RNA; case-insensitive). Must be
    /// exactly 23 nucleotides over the A/C/G/T(/U) alphabet.</param>
    /// <returns>The MaxEntScan score3ss in bits.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="window"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentException"><paramref name="window"/> is not exactly 23 nt, or
    /// contains a non-A/C/G/T(/U) character.</exception>
    public static double ScoreAcceptorMaxEnt(string window)
    {
        ArgumentNullException.ThrowIfNull(window);
        if (window.Length != MaxEntAcceptorWindowLength)
            throw new ArgumentException(
                $"MaxEntScan score3 requires a window of exactly {MaxEntAcceptorWindowLength} nt; " +
                $"got {window.Length}.", nameof(window));

        // Uppercase; keep T/U as-is (HashMaxEntSubsequence and the AG model treat them alike).
        string upper = window.ToUpperInvariant();

        // --- Consensus AG dinucleotide term: P_cons(A)*P_cons(G) / (P_bgd(A)*P_bgd(G)) ---
        char ag1 = NormalizeNucleotide(upper[MaxEntAcceptorAgStart]);
        char ag2 = NormalizeNucleotide(upper[MaxEntAcceptorAgStart + 1]);

        double consScore =
            MaxEntConsensusAg1[ag1] * MaxEntConsensusAg2[ag2] /
            (MaxEntBackground[ag1] * MaxEntBackground[ag2]);

        // --- "Rest" sequence: the 21 positions with the AG removed ---
        // rest = window[0..18) + window[20..23)  -> 18 + 3 = 21 nt.
        Span<char> rest = stackalloc char[MaxEntAcceptorWindowLength - 2];
        upper.AsSpan(0, MaxEntAcceptorAgStart).CopyTo(rest);
        upper.AsSpan(MaxEntAcceptorAgStart + 2).CopyTo(rest[MaxEntAcceptorAgStart..]);

        var tables = MaxEntScore3Tables;

        double restScore = 1.0;
        foreach (var (matrixIndex, restStart, length) in MaxEntNumeratorFactors)
        {
            int hash = HashMaxEntSubsequence(rest.Slice(restStart, length));
            restScore *= tables[matrixIndex][hash];
        }

        foreach (var (matrixIndex, restStart, length) in MaxEntDenominatorFactors)
        {
            int hash = HashMaxEntSubsequence(rest.Slice(restStart, length));
            restScore /= tables[matrixIndex][hash];
        }

        // score3.pl / maxentpy compute log2 as ln(x)/ln(2); use the same expression so the
        // result is the reference double (Math.Log2 differs by up to ~1e-14 on some inputs).
        return Math.Log(consScore * restScore) / Math.Log(2.0);
    }

    private static char NormalizeNucleotide(char c, string model = "score3") => c switch
    {
        'A' or 'C' or 'G' or 'T' => c,
        'U' => 'T',
        _ => throw new ArgumentException(
            $"MaxEntScan {model} requires an A/C/G/T(/U) sequence; found '{c}'.")
    };

    #endregion

    #region MaxEntScan score5ss (maximum-entropy 5' donor model)

    // --- MaxEntScan score5ss (Yeo & Burge 2004) ---
    // Maximum-entropy 5' splice-site (donor) score for a 9-nt window:
    // 3 exonic + 6 intronic nucleotides, with the invariant GT at window positions
    // 3..4 (0-based). The score is log2( P_maxent(seq) / P_background(seq) ), and is
    // computed by removing the conserved GT dinucleotide (scored by a consensus /
    // background model) and looking up the maximum-entropy probability of the remaining
    // 7 ("rest") positions in a single embedded table keyed by the 7-nt string.
    //
    // Unlike score3 (nine overlapping sub-matrices), score5 is single-matrix: the 7-mer
    // rest sequence has its full maximum-entropy probability stored directly.
    //
    // Sources:
    //   Yeo G, Burge CB (2004). "Maximum entropy modeling of short sequence motifs with
    //   applications to RNA splicing signals." J Comput Biol 11(2-3):377-394.
    //   DOI 10.1089/1066527041410418.
    //   Reference factorisation + table (retrieved 2026-06-25): the MIT-licensed maxentpy
    //   port (kepbod/maxentpy) score5 / score5_matrix.txt. Provenance + licence are recorded
    //   in src/.../Seqeron.Genomics.Annotation/Data/maxent_score5.LICENSE.md.

    // Length of the MaxEntScan 5' (donor) scoring window: 3 exon + 6 intron nt.
    private const int MaxEntDonorWindowLength = 9;

    // 0-based positions of the conserved GT dinucleotide within the 9-nt window
    // (window[3]=G, window[4]=T). These two positions are scored by the consensus /
    // background dinucleotide model and removed from the "rest" sub-sequence.
    private const int MaxEntDonorGtStart = 3;

    // Consensus probabilities for the first GT position (the 'G'); maxentpy cons1_5.
    private static readonly IReadOnlyDictionary<char, double> MaxEntDonorConsensusGt1 = new Dictionary<char, double>
    {
        { 'A', 0.004 }, { 'C', 0.0032 }, { 'G', 0.9896 }, { 'T', 0.0032 }
    };

    // Consensus probabilities for the second GT position (the 'T'); maxentpy cons2_5.
    private static readonly IReadOnlyDictionary<char, double> MaxEntDonorConsensusGt2 = new Dictionary<char, double>
    {
        { 'A', 0.0034 }, { 'C', 0.0039 }, { 'G', 0.0042 }, { 'T', 0.9884 }
    };

    private const string MaxEntScore5ResourceName = "Seqeron.Genomics.Annotation.Data.maxent_score5.txt";

    // Lazily-loaded probability table: tables5[rest7mer] = probability. 4^7 = 16384 entries.
    private static IReadOnlyDictionary<string, double>? _maxEntScore5Table;
    private static readonly object MaxEntScore5LoadLock = new();

    private static IReadOnlyDictionary<string, double> MaxEntScore5Table
    {
        get
        {
            if (_maxEntScore5Table is not null)
                return _maxEntScore5Table;

            lock (MaxEntScore5LoadLock)
            {
                if (_maxEntScore5Table is not null)
                    return _maxEntScore5Table;
                _maxEntScore5Table = LoadMaxEntScore5Table();
                return _maxEntScore5Table;
            }
        }
    }

    private static IReadOnlyDictionary<string, double> LoadMaxEntScore5Table()
    {
        var assembly = typeof(SpliceSitePredictor).Assembly;
        using Stream? stream = assembly.GetManifestResourceStream(MaxEntScore5ResourceName)
            ?? throw new InvalidOperationException(
                $"Embedded MaxEntScan score5 table '{MaxEntScore5ResourceName}' was not found.");

        var table = new Dictionary<string, double>(StringComparer.Ordinal);

        using var reader = new StreamReader(stream);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;

            // Each record: "<7-nt rest sequence>\t<probability>" (whitespace-separated).
            string[] parts = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length != 2)
                throw new InvalidDataException($"Malformed MaxEntScan score5 record: '{line}'.");

            table[parts[0]] = double.Parse(parts[1], CultureInfo.InvariantCulture);
        }

        return table;
    }

    /// <summary>
    /// Computes the Yeo &amp; Burge (2004) MaxEntScan <c>score5ss</c> maximum-entropy score
    /// (in bits) for a 9-nt 5' splice-site (donor) window: 3 exonic + 6 intronic
    /// nucleotides, with the conserved <c>GT</c> at window positions 3-4 (0-based).
    /// </summary>
    /// <remarks>
    /// <para>
    /// This is the window scorer behind <see cref="DonorScoringModel.MaxEntScan"/> in
    /// <c>FindDonorSites</c> (whose default is the Shapiro &amp; Senapathy 1987 score). The
    /// returned value is the log<sub>2</sub> ratio of the maximum-entropy probability of the
    /// window to its background probability — higher means a stronger 5' splice site (the
    /// canonical documented example <c>cagGTAAGT</c> scores 10.86 bits).
    /// </para>
    /// <para>
    /// The score removes the conserved <c>GT</c> dinucleotide (scored by a consensus /
    /// background model) and looks up the maximum-entropy probability of the remaining 7
    /// positions (<c>window[0:3] + window[5:9]</c>) directly in a single embedded table,
    /// then takes log<sub>2</sub> of the product with the GT term. The probability table is
    /// embedded (see <c>Data/maxent_score5.LICENSE.md</c> for provenance and licence).
    /// </para>
    /// </remarks>
    /// <param name="window">The 9-nt donor window (DNA or RNA; case-insensitive). Must be
    /// exactly 9 nucleotides over the A/C/G/T(/U) alphabet.</param>
    /// <returns>The MaxEntScan score5ss in bits.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="window"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentException"><paramref name="window"/> is not exactly 9 nt, or
    /// contains a non-A/C/G/T(/U) character.</exception>
    public static double ScoreDonorMaxEnt(string window)
    {
        ArgumentNullException.ThrowIfNull(window);
        if (window.Length != MaxEntDonorWindowLength)
            throw new ArgumentException(
                $"MaxEntScan score5 requires a window of exactly {MaxEntDonorWindowLength} nt; " +
                $"got {window.Length}.", nameof(window));

        // Uppercase, then normalise U->T so the rest 7-mer key matches the table's A/C/G/T keys.
        string upper = window.ToUpperInvariant();

        // --- Consensus GT dinucleotide term: P_cons(G)*P_cons(T) / (P_bgd(G)*P_bgd(T)) ---
        char gt1 = NormalizeNucleotide(upper[MaxEntDonorGtStart], "score5");
        char gt2 = NormalizeNucleotide(upper[MaxEntDonorGtStart + 1], "score5");

        double consScore =
            MaxEntDonorConsensusGt1[gt1] * MaxEntDonorConsensusGt2[gt2] /
            (MaxEntBackground[gt1] * MaxEntBackground[gt2]);

        // --- "Rest" sequence: the 7 positions with the GT removed ---
        // rest = window[0..3) + window[5..9)  -> 3 + 4 = 7 nt; normalise U->T for the key.
        Span<char> rest = stackalloc char[MaxEntDonorWindowLength - 2];
        for (int i = 0; i < MaxEntDonorGtStart; i++)
            rest[i] = NormalizeNucleotide(upper[i], "score5");
        for (int i = MaxEntDonorGtStart + 2; i < MaxEntDonorWindowLength; i++)
            rest[i - 2] = NormalizeNucleotide(upper[i], "score5");

        double restScore = MaxEntScore5Table[new string(rest)];

        // score5.pl / maxentpy compute log2 as ln(x)/ln(2); use the same expression so the
        // result is the reference double (Math.Log2 differs by up to ~1e-14 on some inputs).
        return Math.Log(consScore * restScore) / Math.Log(2.0);
    }

    #endregion
}
