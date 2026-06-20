using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides algorithms for RNA secondary structure prediction and analysis.
/// Includes stem-loop detection, hairpin finding, and free energy calculations.
/// </summary>
public static class RnaSecondaryStructure
{
    #region Records and Types

    /// <summary>
    /// Represents a base pair in RNA secondary structure.
    /// </summary>
    public readonly record struct BasePair(int Position1, int Position2, char Base1, char Base2, BasePairType Type);

    /// <summary>
    /// Types of base pairs.
    /// </summary>
    public enum BasePairType
    {
        WatsonCrick,  // A-U, G-C
        Wobble,       // G-U
        NonCanonical  // Other pairs
    }

    /// <summary>
    /// Represents a stem structure (double-stranded region).
    /// </summary>
    public readonly record struct Stem(
        int Start5Prime,
        int End5Prime,
        int Start3Prime,
        int End3Prime,
        int Length,
        IReadOnlyList<BasePair> BasePairs,
        double FreeEnergy);

    /// <summary>
    /// Represents a loop structure.
    /// </summary>
    public readonly record struct Loop(
        LoopType Type,
        int Start,
        int End,
        int Size,
        string Sequence);

    /// <summary>
    /// Types of loops in RNA structure.
    /// </summary>
    public enum LoopType
    {
        Hairpin,      // Terminal loop
        Internal,     // Internal loop
        Bulge,        // Bulge loop (asymmetric)
        MultiLoop,    // Multi-branch loop
        External      // External (unpaired ends)
    }

    /// <summary>
    /// Represents a stem-loop (hairpin) structure.
    /// </summary>
    public readonly record struct StemLoop(
        int Start,
        int End,
        Stem Stem,
        Loop Loop,
        double TotalFreeEnergy,
        string DotBracketNotation);

    /// <summary>
    /// Represents a pseudoknot structure.
    /// </summary>
    public readonly record struct Pseudoknot(
        int Start1, int End1,
        int Start2, int End2,
        IReadOnlyList<BasePair> CrossingPairs);

    /// <summary>
    /// Complete secondary structure prediction result.
    /// </summary>
    public readonly record struct SecondaryStructure(
        string Sequence,
        string DotBracket,
        IReadOnlyList<BasePair> BasePairs,
        IReadOnlyList<StemLoop> StemLoops,
        IReadOnlyList<Pseudoknot> Pseudoknots,
        double MinimumFreeEnergy);

    /// <summary>
    /// Result of the McCaskill equilibrium partition-function computation over the
    /// full ensemble of (pseudoknot-free) secondary structures of a sequence.
    /// </summary>
    /// <param name="PartitionFunction">
    /// Z = Σ_S exp(−E(S)/RT) over all admissible structures S (McCaskill 1990, Z = Q₁ₙ).
    /// </param>
    /// <param name="BasePairProbabilities">
    /// P[i,j] = (Σ_{S ∋ (i,j)} exp(−E(S)/RT)) / Z, the equilibrium base-pair probability
    /// for every pair (i,j) that can form, indexed by 0-based (i,j) with i &lt; j.
    /// </param>
    public readonly record struct PartitionFunctionResult(
        double PartitionFunction,
        IReadOnlyDictionary<(int I, int J), double> BasePairProbabilities);

    #endregion

    #region Free Energy Parameters

    // Turner 2004 nearest-neighbor parameters (kcal/mol at 37°C)
    // Source: NNDB — https://rna.urmc.rochester.edu/NNDB/turner04/
    //
    // Watson-Crick stacking: rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html
    // GU wobble stacking:    rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html
    // Format: 5'-XY-3' / 3'-X'Y'-5' where X pairs with X', Y pairs with Y'
    private static readonly Dictionary<string, double> StackingEnergies = new()
    {
        // Watson-Crick stacking (10 unique + 6 symmetric = 16)
        { "AA/UU", -0.93 }, { "UU/AA", -0.93 },
        { "AU/UA", -1.10 },
        { "UA/AU", -1.33 },
        { "CU/GA", -2.08 }, { "AG/UC", -2.08 },
        { "CA/GU", -2.11 }, { "UG/AC", -2.11 },
        { "GU/CA", -2.24 }, { "AC/UG", -2.24 },
        { "GA/CU", -2.35 }, { "UC/AG", -2.35 },
        { "CG/GC", -2.36 },
        { "GG/CC", -3.26 }, { "CC/GG", -3.26 },
        { "GC/CG", -3.42 },
        // GU wobble stacking (11 unique + 9 symmetric = 20)
        { "AG/UU", -0.55 }, { "UU/GA", -0.55 },
        { "UG/AU", -1.00 }, { "UA/GU", -1.00 },
        { "GA/UU", -1.27 }, { "UU/AG", -1.27 },
        { "AU/UG", -1.36 }, { "GU/UA", -1.36 },
        { "CG/GU", -1.41 }, { "UG/GC", -1.41 },
        { "GG/CU", -1.53 }, { "UC/GG", -1.53 },
        { "CU/GG", -2.11 }, { "GG/UC", -2.11 },
        { "GU/CG", -2.51 }, { "GC/UG", -2.51 },
        { "GG/UU", -0.50 }, { "UU/GG", -0.50 },
        { "UG/GU", +0.30 },
        { "GU/UG", +1.29 },
    };

    // Hairpin loop initiation energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/loop.txt
    // Sizes 3-9 are experimentally determined; 10-30 are extrapolated via
    // ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9) where R=1.987 cal/(mol·K), T=310.15K
    private static readonly Dictionary<int, double> HairpinLoopEnergies = new()
    {
        {  3, 5.4 }, {  4, 5.6 }, {  5, 5.7 }, {  6, 5.4 }, {  7, 6.0 },
        {  8, 5.5 }, {  9, 6.4 }, { 10, 6.5 }, { 11, 6.6 }, { 12, 6.7 },
        { 13, 6.8 }, { 14, 6.9 }, { 15, 6.9 }, { 16, 7.0 }, { 17, 7.1 },
        { 18, 7.1 }, { 19, 7.2 }, { 20, 7.2 }, { 21, 7.3 }, { 22, 7.3 },
        { 23, 7.4 }, { 24, 7.4 }, { 25, 7.5 }, { 26, 7.5 }, { 27, 7.5 },
        { 28, 7.6 }, { 29, 7.6 }, { 30, 7.7 }
    };

    // Special hairpin loops — total free energies that replace normal model calculation.
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/hairpin-special-parameters.html
    // Key = closing5' + loopSequence + closing3' (includes the closing base pair).
    // These values come from experimental data and supersede initiation + mismatch model.
    private static readonly Dictionary<string, double> SpecialHairpinLoops = new()
    {
        // Triloops (closing + 3nt loop + closing = 5-mer key)
        { "CAACG", 6.8 },
        { "GUUAC", 6.9 },
        // Tetraloops (closing + 4nt loop + closing = 6-mer key)
        { "CCUCGG", 2.5 },
        { "CUCCGG", 2.7 },
        { "CUACGG", 2.8 },
        { "CUGCGG", 2.8 },
        { "CCAAGG", 3.3 },
        { "CCCAGG", 3.4 },
        { "CCGAGG", 3.5 },
        { "CUUAGG", 3.5 },
        { "CCGCGG", 3.6 },
        { "CUAAGG", 3.6 },
        { "CCUAGG", 3.7 },
        { "CCACGG", 3.7 },
        { "CUCAGG", 3.7 },
        { "CUUCGG", 3.7 },
        { "CUUUGG", 3.7 },
        { "CAACGG", 5.5 },
        // Hexaloops (closing + 6nt loop + closing = 8-mer key)
        { "ACAGUGUU", 1.8 },
        { "ACAGUACU", 2.8 },
        { "ACAGUGCU", 2.9 },
        { "ACAGUGAU", 3.6 },
    };

    // Terminal mismatch stacking energies (kcal/mol at 37°C).
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/tm-parameters.html
    // Key = closing5' + firstMismatch(5'side) + lastMismatch(3'side) + closing3'.
    // These are the sequence-dependent stacking of the first non-canonical pair on the closing pair.
    // Applied to hairpin loops ≥4nt, internal loops (other sizes), and exterior loops.
    private static readonly Dictionary<string, double> TerminalMismatchEnergies = new()
    {
        // Closing pair A-U (5'AX / 3'UY)
        { "AAAU", -0.8 }, { "AACU", -1.0 }, { "AAGU", -0.8 }, { "AAUU", -1.0 },
        { "ACAU", -0.6 }, { "ACCU", -0.7 }, { "ACGU", -0.6 }, { "ACUU", -0.7 },
        { "AGAU", -0.8 }, { "AGCU", -1.0 }, { "AGGU", -0.8 }, { "AGUU", -1.0 },
        { "AUAU", -0.6 }, { "AUCU", -0.8 }, { "AUGU", -0.6 }, { "AUUU", -0.8 },
        // Closing pair C-G (5'CX / 3'GY)
        { "CAAG", -1.5 }, { "CACG", -1.5 }, { "CAGG", -1.4 }, { "CAUG", -1.5 },
        { "CCAG", -1.0 }, { "CCCG", -1.1 }, { "CCGG", -1.0 }, { "CCUG", -0.8 },
        { "CGAG", -1.4 }, { "CGCG", -1.5 }, { "CGGG", -1.6 }, { "CGUG", -1.5 },
        { "CUAG", -1.0 }, { "CUCG", -1.4 }, { "CUGG", -1.0 }, { "CUUG", -1.2 },
        // Closing pair G-C (5'GX / 3'CY)
        { "GAAC", -1.1 }, { "GACC", -1.5 }, { "GAGC", -1.3 }, { "GAUC", -1.5 },
        { "GCAC", -1.1 }, { "GCCC", -0.7 }, { "GCGC", -1.1 }, { "GCUC", -0.5 },
        { "GGAC", -1.6 }, { "GGCC", -1.5 }, { "GGGC", -1.4 }, { "GGUC", -1.5 },
        { "GUAC", -1.1 }, { "GUCC", -1.0 }, { "GUGC", -1.1 }, { "GUUC", -0.7 },
        // Closing pair G-U (5'GX / 3'UY)
        { "GAAU", -0.3 }, { "GACU", -1.0 }, { "GAGU", -0.8 }, { "GAUU", -1.0 },
        { "GCAU", -0.6 }, { "GCCU", -0.7 }, { "GCGU", -0.6 }, { "GCUU", -0.7 },
        { "GGAU", -0.6 }, { "GGCU", -1.0 }, { "GGGU", -0.8 }, { "GGUU", -1.0 },
        { "GUAU", -0.6 }, { "GUCU", -0.8 }, { "GUGU", -0.6 }, { "GUUU", -0.6 },
        // Closing pair U-A (5'UX / 3'AY)
        { "UAAA", -1.0 }, { "UACA", -0.8 }, { "UAGA", -1.1 }, { "UAUA", -0.8 },
        { "UCAA", -0.7 }, { "UCCA", -0.6 }, { "UCGA", -0.7 }, { "UCUA", -0.5 },
        { "UGAA", -1.1 }, { "UGCA", -0.8 }, { "UGGA", -1.2 }, { "UGUA", -0.8 },
        { "UUAA", -0.7 }, { "UUCA", -0.6 }, { "UUGA", -0.7 }, { "UUUA", -0.5 },
        // Closing pair U-G (5'UX / 3'GY)
        { "UAAG", -1.0 }, { "UACG", -0.8 }, { "UAGG", -1.1 }, { "UAUG", -0.8 },
        { "UCAG", -0.7 }, { "UCCG", -0.6 }, { "UCGG", -0.7 }, { "UCUG", -0.5 },
        { "UGAG", -0.5 }, { "UGCG", -0.8 }, { "UGGG", -0.8 }, { "UGUG", -0.8 },
        { "UUAG", -0.7 }, { "UUCG", -0.6 }, { "UUGG", -0.7 }, { "UUUG", -0.5 },
    };

    // Terminal AU/GU penalty (kcal/mol at 37°C)
    // Source: NNDB — "Per AU end" on wc-parameters.html, "Per GU end" on gu-parameters.html
    // Applied at each end of a helix that terminates with an AU/UA or GU/UG base pair.
    private const double TerminalAU_GU_Penalty = 0.45;

    // Hairpin loop bonus/penalty parameters (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/hairpin-mismatch-parameters.html
    // These are ADDITIVE on top of the terminal mismatch table.
    private const double UU_GA_MismatchBonus = -0.9;
    private const double GG_MismatchBonus = -0.8;
    private const double SpecialGU_ClosureBonus = -2.2;
    private const double AllC_PenaltyPerNt = 0.3;
    private const double AllC_PenaltyConstant_Gt3 = 1.6;
    private const double AllC_Penalty_3nt = 1.5;

    // Special GGUC/CUGG 3-stack context (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html (note b)
    // GU followed by UG is generally unfavorable (+1.29), but when preceded by GC and
    // followed by CG (5'GGUC/3'CUGG), the total for 3 stacks = -4.12 replaces individual sum.
    private const double SpecialGGUC_CUGG_3Stack = -4.12;

    // Internal loop initiation energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/loop.txt (INTERNAL column)
    // For n>6 (experimentally determined up to 6): init(6) + 1.08·ln(n/6)
    // For n>30: extrapolation from size 6.
    private static readonly Dictionary<int, double> InternalLoopInitiation = new()
    {
        {  4, 1.1 }, {  5, 2.0 }, {  6, 2.0 }, {  7, 2.1 }, {  8, 2.3 },
        {  9, 2.4 }, { 10, 2.5 }, { 11, 2.6 }, { 12, 2.7 }, { 13, 2.8 },
        { 14, 2.9 }, { 15, 2.9 }, { 16, 3.0 }, { 17, 3.1 }, { 18, 3.1 },
        { 19, 3.2 }, { 20, 3.3 }, { 21, 3.3 }, { 22, 3.4 }, { 23, 3.4 },
        { 24, 3.5 }, { 25, 3.5 }, { 26, 3.5 }, { 27, 3.6 }, { 28, 3.6 },
        { 29, 3.7 }, { 30, 3.7 }
    };

    // Internal loop parameters (kcal/mol)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/internal-parameters.html
    private const double InternalLoop_AsymmetryPenalty = 0.6;
    private const double InternalLoop_AU_GU_Closure = 0.7;

    // Internal loop mismatch parameters for 2×3 loops (orientation-dependent)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/internal-parameters.html
    // In 2×3 context, R=purine(A,G), Y=pyrimidine(C,U). Orientation of closing pair matters.
    private const double IL_2x3_RA_YG = 0.0;
    private const double IL_2x3_YA_RG = -0.5;
    private const double IL_2x3_RG_YA = -1.2;
    private const double IL_2x3_YG_RA = -1.1;
    private const double IL_2x3_GG = -0.8;
    private const double IL_2x3_UU = -0.4;

    // Internal loop mismatch parameters for other sizes (not 1×1, 1×2, 2×2, 2×3, 1×(n-1))
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/internal-parameters.html
    private const double IL_Other_AG = -0.8;
    private const double IL_Other_GA = -1.0;
    private const double IL_Other_GG = -1.2;
    private const double IL_Other_UU = -0.7;

    // 1×1 internal loop free energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/int11.txt
    // Indexed by [pair1_index, pair2_index, mismatch_X_index, mismatch_Y_index]
    // Pair indices: 0=AU, 1=CG, 2=GC, 3=UA, 4=GU, 5=UG
    // Base indices: 0=A, 1=C, 2=G, 3=U
    // NOTE: Values already include AU/GU closure penalties.
    private static readonly double[,,,] Int11Energies =
    {
        {   // pair1=AU
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.5 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.8 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.8 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.2 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.2 } },
        },
        {   // pair1=CG
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.2 } },
            { { 0.9, -0.4, 0.5, 0.5 }, { 0.3, 0.5, 0.5, 0.6 }, { -0.1, 0.5, -2.2, 0.5 }, { 0.5, 0.0, 0.5, -0.1 } },
            { { 0.9, 0.5, 0.5, 0.5 }, { 0.5, 0.5, 0.5, 0.5 }, { 0.5, 0.5, -1.4, 0.5 }, { 0.5, 0.5, 0.5, 0.4 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.8 } },
            { { 2.2, 1.3, 1.2, 1.2 }, { 1.2, 1.7, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.1 } },
            { { 0.6, 0.5, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { -0.2, 1.2, -1.4, 1.2 }, { 1.2, 1.0, 1.2, 1.1 } },
        },
        {   // pair1=GC
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.2 } },
            { { 0.8, 0.5, 0.5, 0.5 }, { 0.5, 0.5, 0.5, 0.5 }, { 0.5, 0.5, -2.3, 0.5 }, { 0.5, 0.5, 0.5, -0.6 } },
            { { 0.9, 0.3, -0.1, 0.5 }, { -0.4, 0.5, 0.5, 0.0 }, { 0.5, 0.5, -2.2, 0.5 }, { 0.5, 0.6, 0.5, -0.1 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.8 } },
            { { 1.6, 1.2, 1.0, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.7 } },
            { { 1.9, 1.2, 1.5, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.5 } },
        },
        {   // pair1=UA
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.7 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.2 } },
            { { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.2 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.5 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.9 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
        },
        {   // pair1=GU
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
            { { 1.9, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.5, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.5 } },
            { { 0.6, 1.2, -0.2, 1.2 }, { 0.5, 1.2, 1.2, 1.0 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.1 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.2 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.2 } },
        },
        {   // pair1=UG
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.9 } },
            { { 1.6, 1.2, 1.2, 1.2 }, { 1.2, 1.2, 1.2, 1.2 }, { 1.0, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 0.7 } },
            { { 2.2, 1.2, 1.2, 1.2 }, { 1.3, 1.7, 1.2, 1.2 }, { 1.2, 1.2, -1.4, 1.2 }, { 1.2, 1.2, 1.2, 1.1 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.9 } },
            { { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, 1.9, 1.9 }, { 1.9, 1.9, -0.7, 1.9 }, { 1.9, 1.9, 1.9, 1.6 } },
        },
    };

    // Bulge loop initiation energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/loop.txt (BULGE column)
    // For n=1: special formula with stacking continuation.
    // For n>1: just initiation. For n>6: init(6) + 1.75·R·T·ln(n/6).
    private static readonly Dictionary<int, double> BulgeLoopInitiation = new()
    {
        {  1, 3.8 }, {  2, 2.8 }, {  3, 3.2 }, {  4, 3.6 }, {  5, 4.0 },
        {  6, 4.4 }, {  7, 4.6 }, {  8, 4.7 }, {  9, 4.8 }, { 10, 4.9 },
        { 11, 5.0 }, { 12, 5.1 }, { 13, 5.2 }, { 14, 5.3 }, { 15, 5.4 },
        { 16, 5.4 }, { 17, 5.5 }, { 18, 5.5 }, { 19, 5.6 }, { 20, 5.7 },
        { 21, 5.7 }, { 22, 5.8 }, { 23, 5.8 }, { 24, 5.8 }, { 25, 5.9 },
        { 26, 5.9 }, { 27, 6.0 }, { 28, 6.0 }, { 29, 6.0 }, { 30, 6.1 }
    };

    // Special C bulge bonus: bulged C adjacent to at least one paired C.
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/bulge.html
    private const double BulgeSpecialC_Bonus = -0.9;

    // 3' dangling end energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/de-parameters.html
    // Key = closing5' + danglingBase + closing3'. Dangling base is 3' to the closing pair.
    private static readonly Dictionary<string, double> DanglingEnd3 = new()
    {
        // Closing A-U
        { "AAU", -0.8 }, { "ACU", -0.5 }, { "AGU", -0.8 }, { "AUU", -0.6 },
        // Closing C-G
        { "CAG", -1.7 }, { "CCG", -0.8 }, { "CGG", -1.7 }, { "CUG", -1.2 },
        // Closing G-C
        { "GAC", -1.1 }, { "GCC", -0.4 }, { "GGC", -1.3 }, { "GUC", -0.6 },
        // Closing G-U
        { "GAU", -0.8 }, { "GCU", -0.5 }, { "GGU", -0.8 }, { "GUU", -0.6 },
        // Closing U-A
        { "UAA", -0.7 }, { "UCA", -0.1 }, { "UGA", -0.7 }, { "UUA", -0.1 },
        // Closing U-G
        { "UAG", -0.7 }, { "UCG", -0.1 }, { "UGG", -0.7 }, { "UUG", -0.1 },
    };

    // 5' dangling end energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/de-parameters.html
    // Key = closing5' + danglingBase + closing3'. Dangling base is 5' to the closing pair.
    private static readonly Dictionary<string, double> DanglingEnd5 = new()
    {
        // Closing A-U
        { "AAU", -0.3 }, { "ACU", -0.1 }, { "AGU", -0.2 }, { "AUU", -0.2 },
        // Closing C-G
        { "CAG", -0.2 }, { "CCG", -0.3 }, { "CGG",  0.0 }, { "CUG",  0.0 },
        // Closing G-C
        { "GAC", -0.5 }, { "GCC", -0.3 }, { "GGC", -0.2 }, { "GUC", -0.1 },
        // Closing G-U
        { "GAU", -0.3 }, { "GCU", -0.1 }, { "GGU", -0.2 }, { "GUU", -0.2 },
        // Closing U-A
        { "UAA", -0.3 }, { "UCA", -0.3 }, { "UGA", -0.4 }, { "UUA", -0.2 },
        // Closing U-G
        { "UAG", -0.3 }, { "UCG", -0.3 }, { "UGG", -0.4 }, { "UUG", -0.2 },
    };

    // Multibranch loop parameters (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/mb-parameters.html
    // ΔG° = a + b×avg_asymmetry + c×num_helices + strain
    private const double MultibranchLoop_Offset = 9.25;     // a
    private const double MultibranchLoop_Asymmetry = 0.91;  // b
    private const double MultibranchLoop_Helix = -0.63;     // c
    private const double MultibranchLoop_Strain = 3.14;     // 3-way with <2 unpaired

    // Coaxial stacking parameters (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/coax.html
    // Flush coaxial: use WC/GU stacking as if no backbone break.
    // Mismatch-mediated: terminal mismatch + base constant + WC/GU bonus.
    private const double CoaxialMismatch_Base = -2.1;
    private const double CoaxialMismatch_WC_Bonus = -0.4;
    private const double CoaxialMismatch_GU_Bonus = -0.2;

    // Fast pair-type lookup: 0 = no pair, 1 = Watson-Crick, 2 = Wobble.
    // Indexed by [base1 * 128 + base2]. 16 KB — fits L1 cache.
    private static readonly byte[] PairLookup = BuildPairLookup();
    private static readonly double[] PairEnergyByCode = new[] { 0.0, -2.0, -1.0 };

    private static byte[] BuildPairLookup()
    {
        var t = new byte[128 * 128];
        t['A' * 128 + 'U'] = 1; t['U' * 128 + 'A'] = 1;
        t['G' * 128 + 'C'] = 1; t['C' * 128 + 'G'] = 1;
        t['G' * 128 + 'U'] = 2; t['U' * 128 + 'G'] = 2;
        return t;
    }

    #endregion

    #region Base Pairing

    /// <summary>
    /// Determines if two bases can form a pair.
    /// </summary>
    public static bool CanPair(char base1, char base2)
    {
        int b1 = char.ToUpperInvariant(base1);
        int b2 = char.ToUpperInvariant(base2);
        return (b1 | b2) < 128 && PairLookup[b1 * 128 + b2] != 0;
    }

    /// <summary>
    /// Gets the type of base pair, or null if bases cannot pair.
    /// </summary>
    public static BasePairType? GetBasePairType(char base1, char base2)
    {
        int b1 = char.ToUpperInvariant(base1);
        int b2 = char.ToUpperInvariant(base2);
        if ((b1 | b2) >= 128) return null;
        return PairLookup[b1 * 128 + b2] switch
        {
            1 => BasePairType.WatsonCrick,
            2 => BasePairType.Wobble,
            _ => null
        };
    }

    /// <summary>
    /// Gets the complement of a base (RNA).
    /// </summary>
    public static char GetComplement(char base_) => SequenceExtensions.GetRnaComplementBase(base_);

    #endregion

    #region Stem-Loop Finding

    /// <summary>
    /// Finds all potential stem-loop structures in an RNA sequence.
    /// </summary>
    public static IEnumerable<StemLoop> FindStemLoops(
        string rnaSequence,
        int minStemLength = 3,
        int minLoopSize = 3,
        int maxLoopSize = 10,
        bool allowWobble = true)
    {
        // NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops
        // with fewer than 3 nucleotides."
        // Wikipedia Stem-loop: "loops fewer than three bases long are sterically
        // impossible and thus do not form."
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minStemLength * 2 + minLoopSize)
            yield break;

        string upper = rnaSequence.ToUpperInvariant();

        // Scan for potential hairpin loops
        for (int loopStart = minStemLength; loopStart <= upper.Length - minStemLength - minLoopSize; loopStart++)
        {
            for (int loopSize = minLoopSize; loopSize <= Math.Min(maxLoopSize, upper.Length - loopStart - minStemLength); loopSize++)
            {
                int loopEnd = loopStart + loopSize - 1;

                // Try to extend stem on both sides
                var stemLoop = TryBuildStemLoop(upper, loopStart, loopEnd, minStemLength, allowWobble);
                if (stemLoop != null)
                {
                    yield return stemLoop.Value;
                }
            }
        }
    }

    private static StemLoop? TryBuildStemLoop(string sequence, int loopStart, int loopEnd, int minStemLength, bool allowWobble)
    {
        var basePairs = new List<BasePair>();
        int stemLength = 0;

        int left = loopStart - 1;
        int right = loopEnd + 1;

        // Extend stem
        while (left >= 0 && right < sequence.Length)
        {
            var pairType = GetBasePairType(sequence[left], sequence[right]);

            if (pairType == null)
                break;

            if (pairType == BasePairType.Wobble && !allowWobble)
                break;

            basePairs.Add(new BasePair(left, right, sequence[left], sequence[right], pairType.Value));
            stemLength++;
            left--;
            right++;
        }

        if (stemLength < minStemLength)
            return null;

        // Build result
        int stemStart = left + 1;
        int stemEnd5 = loopStart - 1;
        int stemStart3 = loopEnd + 1;
        int stemEnd = right - 1;

        basePairs.Reverse(); // Order from 5' to 3'

        var stem = new Stem(
            Start5Prime: stemStart,
            End5Prime: stemEnd5,
            Start3Prime: stemStart3,
            End3Prime: stemEnd,
            Length: stemLength,
            BasePairs: basePairs,
            FreeEnergy: CalculateStemEnergy(sequence, basePairs));

        string loopSeq = sequence.Substring(loopStart, loopEnd - loopStart + 1);
        var loop = new Loop(
            Type: LoopType.Hairpin,
            Start: loopStart,
            End: loopEnd,
            Size: loopSeq.Length,
            Sequence: loopSeq);

        double loopEnergy = CalculateHairpinLoopEnergy(
            loopSeq, sequence[stemEnd5], sequence[stemStart3],
            IsSpecialGUClosure(sequence, basePairs));
        double totalEnergy = stem.FreeEnergy + loopEnergy;

        string dotBracket = GenerateDotBracket(sequence.Length, basePairs, stemStart, stemEnd);

        return new StemLoop(
            Start: stemStart,
            End: stemEnd,
            Stem: stem,
            Loop: loop,
            TotalFreeEnergy: totalEnergy,
            DotBracketNotation: dotBracket);
    }

    #endregion

    #region Energy Calculations

    /// <summary>
    /// Calculates the free energy of a stem region using Turner 2004 nearest-neighbor
    /// stacking parameters from NNDB, including terminal AU/GU penalties.
    /// </summary>
    public static double CalculateStemEnergy(string sequence, IReadOnlyList<BasePair> basePairs)
    {
        if (basePairs.Count == 0)
            return 0;

        double energy = 0;

        // Sum stacking energies (Turner 2004 nearest-neighbor)
        for (int i = 0; i < basePairs.Count - 1; i++)
        {
            // Special GGUC/CUGG 3-stack context (NNDB note b):
            // When 4 consecutive pairs are G-C, G-U, U-G, C-G, replace
            // 3 individual stacking terms with -4.12 total.
            if (i + 3 < basePairs.Count &&
                basePairs[i].Base1 == 'G' && basePairs[i].Base2 == 'C' &&
                basePairs[i + 1].Base1 == 'G' && basePairs[i + 1].Base2 == 'U' &&
                basePairs[i + 2].Base1 == 'U' && basePairs[i + 2].Base2 == 'G' &&
                basePairs[i + 3].Base1 == 'C' && basePairs[i + 3].Base2 == 'G')
            {
                energy += SpecialGGUC_CUGG_3Stack;
                i += 2; // Skip next 2 iterations (3 stacks consumed)
                continue;
            }

            var pair1 = basePairs[i];
            var pair2 = basePairs[i + 1];

            string stack = $"{pair1.Base1}{pair2.Base1}/{pair1.Base2}{pair2.Base2}";
            if (StackingEnergies.TryGetValue(stack, out double stackEnergy))
                energy += stackEnergy;
            // Unknown stacking pairs contribute 0 (no data available)
        }

        // Terminal AU/GU penalty: +0.45 per helix end that terminates with AU/UA or GU/UG
        // Source: NNDB — "Per AU end" (wc-parameters.html), "Per GU end" (gu-parameters.html)
        if (IsTerminalAUorGU(basePairs[0]))
            energy += TerminalAU_GU_Penalty;
        if (IsTerminalAUorGU(basePairs[^1]))
            energy += TerminalAU_GU_Penalty;

        return Math.Round(energy, 2);
    }

    private static bool IsTerminalAUorGU(BasePair pair)
    {
        return pair.Type == BasePairType.Wobble ||
               (pair.Base1 == 'A' && pair.Base2 == 'U') ||
               (pair.Base1 == 'U' && pair.Base2 == 'A');
    }

    /// <summary>
    /// Detects whether the innermost (closing) pair qualifies for the special GU closure bonus.
    /// Requires: closing pair is G-U (not U-G) AND the two preceding stem pairs have G on the 5' side.
    /// </summary>
    private static bool IsSpecialGUClosure(string sequence, IReadOnlyList<BasePair> basePairs)
    {
        if (basePairs.Count < 3)
            return false;

        var closing = basePairs[^1];
        if (closing.Base1 != 'G' || closing.Base2 != 'U')
            return false;

        return basePairs[^2].Base1 == 'G' && basePairs[^3].Base1 == 'G';
    }

    /// <summary>
    /// Returns the terminal mismatch stacking energy for a closing pair and adjacent mismatched bases.
    /// </summary>
    public static double GetTerminalMismatchEnergy(char closingBase5, char closingBase3, char mismatch5, char mismatch3)
    {
        string key = $"{char.ToUpperInvariant(closingBase5)}{char.ToUpperInvariant(mismatch5)}" +
                     $"{char.ToUpperInvariant(mismatch3)}{char.ToUpperInvariant(closingBase3)}";
        return TerminalMismatchEnergies.TryGetValue(key, out double e) ? e : 0.0;
    }

    /// <summary>
    /// Returns the dangling end energy contribution.
    /// </summary>
    public static double GetDanglingEndEnergy(char closingBase5, char closingBase3, char danglingBase, bool is3Prime)
    {
        string key = $"{char.ToUpperInvariant(closingBase5)}{char.ToUpperInvariant(danglingBase)}" +
                     $"{char.ToUpperInvariant(closingBase3)}";
        var table = is3Prime ? DanglingEnd3 : DanglingEnd5;
        return table.TryGetValue(key, out double e) ? e : 0.0;
    }

    /// <summary>
    /// Calculates the free energy of a hairpin loop using Turner 2004 parameters from NNDB.
    /// For special loops (tri/tetra/hexaloops), the total energy from experimental data
    /// replaces the model calculation. Otherwise uses:
    /// ΔG° = initiation(n) + terminal_mismatch + UU/GA_bonus + GG_bonus + special_GU_closure + all_C_penalty
    /// </summary>
    /// <param name="specialGUClosure">True when closing pair is G-U preceded by two Gs on 5' side.</param>
    public static double CalculateHairpinLoopEnergy(
        string loopSequence, char closingBase5, char closingBase3,
        bool specialGUClosure = false)
    {
        string loopUpper = loopSequence.ToUpperInvariant();
        int size = loopUpper.Length;
        char c5 = char.ToUpperInvariant(closingBase5);
        char c3 = char.ToUpperInvariant(closingBase3);

        // Check for special hairpin loop (NNDB experimental total energy).
        // Key = closing5' + loop + closing3'.
        string specialKey = $"{c5}{loopUpper}{c3}";
        if (SpecialHairpinLoops.TryGetValue(specialKey, out double specialEnergy))
            return Math.Round(specialEnergy, 2);

        // Loop initiation energy
        double energy;
        if (HairpinLoopEnergies.TryGetValue(size, out double loopEnergy))
        {
            energy = loopEnergy;
        }
        else if (size > 30)
        {
            // Jacobson-Stockmayer extrapolation:
            // ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9)
            energy = 6.4 + 1.75 * 1.987 * 310.15 / 1000.0 * Math.Log((double)size / 9.0);
        }
        else
        {
            // Loops < 3 nt are sterically impossible (NNDB Turner 2004, Wikipedia Stem-loop).
            // Return the prohibitive sentinel IMMEDIATELY: the sequence-dependent terms below
            // (terminal mismatch, first-mismatch bonuses, all-C penalty) are defined only for
            // valid loops and must not be layered onto the sentinel. Returning here keeps the
            // sentinel EXACTLY 100.0 (INV-02) and avoids the all-C check's vacuous-truth on an
            // empty loop string (All(c => c == 'C') is true for "" → would add the +1.6 penalty).
            return 100.0;
        }

        // For loops ≥4nt: terminal mismatch stacking + sequence-dependent bonuses + penalties
        if (size >= 4)
        {
            char firstMismatch = loopUpper[0];
            char lastMismatch = loopUpper[^1];

            // Terminal mismatch stacking energy (96-entry table from NNDB)
            string tmKey = $"{c5}{firstMismatch}{lastMismatch}{c3}";
            if (TerminalMismatchEnergies.TryGetValue(tmKey, out double tmEnergy))
                energy += tmEnergy;

            // UU or GA first mismatch bonus — additive on top of terminal mismatch
            if ((firstMismatch == 'U' && lastMismatch == 'U') ||
                (firstMismatch == 'G' && lastMismatch == 'A'))
            {
                energy += UU_GA_MismatchBonus;
            }

            // GG first mismatch bonus — additive on top of terminal mismatch
            if (firstMismatch == 'G' && lastMismatch == 'G')
            {
                energy += GG_MismatchBonus;
            }

            // Special GU closure: -2.2 when closing pair is G-U preceded by two Gs
            if (specialGUClosure && c5 == 'G' && c3 == 'U')
            {
                energy += SpecialGU_ClosureBonus;
            }
        }

        // All-C loop penalty
        if (loopUpper.All(c => c == 'C'))
        {
            if (size == 3)
                energy += AllC_Penalty_3nt;
            else
                energy += AllC_PenaltyPerNt * size + AllC_PenaltyConstant_Gt3;
        }

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the free energy of an internal loop using Turner 2004 parameters.
    /// ΔG° = initiation(n₁+n₂) + |n₁−n₂|·asymmetry + mismatch_1 + mismatch_2 + AU/GU_closure
    /// </summary>
    /// <param name="n1">Number of unpaired bases on the 5' side.</param>
    /// <param name="n2">Number of unpaired bases on the 3' side.</param>
    /// <param name="closingBase5_1">5' base of closing pair 1 (outer).</param>
    /// <param name="closingBase3_1">3' base of closing pair 1 (outer).</param>
    /// <param name="closingBase5_2">5' base of closing pair 2 (inner).</param>
    /// <param name="closingBase3_2">3' base of closing pair 2 (inner).</param>
    /// <param name="mismatch5_1">First unpaired base adjacent to closing pair 1.</param>
    /// <param name="mismatch3_1">Last unpaired base adjacent to closing pair 1.</param>
    /// <param name="mismatch5_2">First unpaired base adjacent to closing pair 2.</param>
    /// <param name="mismatch3_2">Last unpaired base adjacent to closing pair 2.</param>
    public static double CalculateInternalLoopEnergy(
        int n1, int n2,
        char closingBase5_1, char closingBase3_1,
        char closingBase5_2, char closingBase3_2,
        char mismatch5_1, char mismatch3_1,
        char mismatch5_2, char mismatch3_2)
    {
        char c5_1 = char.ToUpperInvariant(closingBase5_1);
        char c3_1 = char.ToUpperInvariant(closingBase3_1);
        char c5_2 = char.ToUpperInvariant(closingBase5_2);
        char c3_2 = char.ToUpperInvariant(closingBase3_2);

        // 1×1 internal loops: use NNDB int11 lookup table (includes AU/GU penalties)
        if (n1 == 1 && n2 == 1)
        {
            int pi1 = Int11PairIndex(c5_1, c3_1);
            int pi2 = Int11PairIndex(c5_2, c3_2);
            int xi = Int11BaseIndex(char.ToUpperInvariant(mismatch5_1));
            int yi = Int11BaseIndex(char.ToUpperInvariant(mismatch3_1));
            if (pi1 >= 0 && pi2 >= 0 && xi >= 0 && yi >= 0)
                return Int11Energies[pi1, pi2, xi, yi];
        }

        int totalSize = n1 + n2;

        // Initiation energy
        double energy;
        if (InternalLoopInitiation.TryGetValue(totalSize, out double initEnergy))
        {
            energy = initEnergy;
        }
        else if (totalSize > 30)
        {
            energy = 2.0 + 1.08 * Math.Log((double)totalSize / 6.0);
        }
        else
        {
            energy = 1.1; // Minimum (size 4)
        }

        // Asymmetry penalty
        energy += InternalLoop_AsymmetryPenalty * Math.Abs(n1 - n2);

        // AU/GU closure penalty on each side
        if (IsAUorGU(c5_1, c3_1))
            energy += InternalLoop_AU_GU_Closure;
        if (IsAUorGU(c5_2, c3_2))
            energy += InternalLoop_AU_GU_Closure;

        // Mismatch contributions (depend on loop geometry)
        char mm5_1 = char.ToUpperInvariant(mismatch5_1);
        char mm3_1 = char.ToUpperInvariant(mismatch3_1);
        char mm5_2 = char.ToUpperInvariant(mismatch5_2);
        char mm3_2 = char.ToUpperInvariant(mismatch3_2);

        if ((n1 == 2 && n2 == 3) || (n1 == 3 && n2 == 2))
        {
            // 2×3 loops: orientation-dependent mismatch
            energy += GetInternalLoop2x3MismatchEnergy(mm5_1, mm3_1);
            energy += GetInternalLoop2x3MismatchEnergy(mm5_2, mm3_2);
        }
        else if (n1 >= 2 && n2 >= 2)
        {
            // Other internal loops: sequence-dependent mismatch
            energy += GetInternalLoopOtherMismatchEnergy(mm5_1, mm3_1);
            energy += GetInternalLoopOtherMismatchEnergy(mm5_2, mm3_2);
        }

        return Math.Round(energy, 2);
    }

    private static bool IsAUorGU(char b5, char b3)
    {
        return (b5 == 'A' && b3 == 'U') || (b5 == 'U' && b3 == 'A') ||
               (b5 == 'G' && b3 == 'U') || (b5 == 'U' && b3 == 'G');
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int Int11PairIndex(char c5, char c3)
    {
        return (c5, c3) switch
        {
            ('A', 'U') => 0,
            ('C', 'G') => 1,
            ('G', 'C') => 2,
            ('U', 'A') => 3,
            ('G', 'U') => 4,
            ('U', 'G') => 5,
            _ => -1
        };
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int Int11BaseIndex(char b)
    {
        return b switch { 'A' => 0, 'C' => 1, 'G' => 2, 'U' => 3, _ => -1 };
    }

    private static double GetInternalLoop2x3MismatchEnergy(char mm5, char mm3)
    {
        bool isPurine5 = mm5 is 'A' or 'G';
        bool isPurine3 = mm3 is 'A' or 'G';

        if (mm5 == 'G' && mm3 == 'G') return IL_2x3_GG;
        if (mm5 == 'U' && mm3 == 'U') return IL_2x3_UU;
        if (isPurine5 && mm3 == 'A') return isPurine3 ? IL_2x3_RA_YG : IL_2x3_RA_YG; // RA
        if (!isPurine5 && mm3 == 'A') return isPurine3 ? IL_2x3_YA_RG : IL_2x3_YA_RG; // YA
        if (isPurine5 && mm3 == 'G') return IL_2x3_RG_YA;
        if (!isPurine5 && mm3 == 'G') return IL_2x3_YG_RA;

        return 0.0;
    }

    private static double GetInternalLoopOtherMismatchEnergy(char mm5, char mm3)
    {
        if (mm5 == 'A' && mm3 == 'G') return IL_Other_AG;
        if (mm5 == 'G' && mm3 == 'A') return IL_Other_GA;
        if (mm5 == 'G' && mm3 == 'G') return IL_Other_GG;
        if (mm5 == 'U' && mm3 == 'U') return IL_Other_UU;
        return 0.0;
    }

    /// <summary>
    /// Calculates the free energy of a bulge loop using Turner 2004 parameters.
    /// n=1: initiation + stacking(as if no bulge) + special_C_bonus − RT·ln(numStates).
    /// n>1: initiation + terminal AU/GU penalties on both sides.
    /// Source: rna.urmc.rochester.edu/NNDB/turner04/bulge.html
    /// </summary>
    /// <param name="bulgeSize">Number of unpaired bases in the bulge.</param>
    /// <param name="bulgedBase">The single bulged nucleotide (used for special C check when n=1).</param>
    /// <param name="pair5_base1">5' base of the pair on the 5' side of the bulge.</param>
    /// <param name="pair5_base2">3' base of the pair on the 5' side.</param>
    /// <param name="pair3_base1">5' base of the pair on the 3' side of the bulge.</param>
    /// <param name="pair3_base2">3' base of the pair on the 3' side.</param>
    /// <param name="numStates">Number of equivalent states for degeneracy (n=1 only). Default 1 = no degeneracy.</param>
    public static double CalculateBulgeLoopEnergy(
        int bulgeSize, char bulgedBase,
        char pair5_base1, char pair5_base2,
        char pair3_base1, char pair3_base2,
        int numStates = 1)
    {
        // Initiation energy
        double energy;
        if (BulgeLoopInitiation.TryGetValue(bulgeSize, out double initEnergy))
        {
            energy = initEnergy;
        }
        else if (bulgeSize > 30)
        {
            // Extrapolation: init(6) + 1.75·R·T·ln(n/6)
            energy = 4.4 + 1.75 * 1.987 * 310.15 / 1000.0 * Math.Log((double)bulgeSize / 6.0);
        }
        else
        {
            energy = 3.8; // Minimum
        }

        char p5b1 = char.ToUpperInvariant(pair5_base1);
        char p5b2 = char.ToUpperInvariant(pair5_base2);
        char p3b1 = char.ToUpperInvariant(pair3_base1);
        char p3b2 = char.ToUpperInvariant(pair3_base2);

        if (bulgeSize == 1)
        {
            // For n=1: stacking as if the bulge didn't exist
            string stack = $"{p5b1}{p3b1}/{p5b2}{p3b2}";
            if (StackingEnergies.TryGetValue(stack, out double stackEnergy))
                energy += stackEnergy;

            // Special C bulge bonus: bulged C adjacent to at least one paired C
            char bulge = char.ToUpperInvariant(bulgedBase);
            if (bulge == 'C' && (p5b1 == 'C' || p5b2 == 'C' || p3b1 == 'C' || p3b2 == 'C'))
                energy += BulgeSpecialC_Bonus;

            // Degeneracy: −RT·ln(numStates) — entropy from equivalent bulge positions
            // Source: NNDB bulge.html, Example 1: 3 C states → −0.616·ln(3) = −0.677
            if (numStates > 1)
            {
                const double RT = 1.987 * 310.15 / 1000.0; // 0.6163 kcal/mol
                energy -= RT * Math.Log(numStates);
            }
        }
        else
        {
            // For n>1: terminal AU/GU penalty on each side
            if (IsAUorGU(p5b1, p5b2))
                energy += TerminalAU_GU_Penalty;
            if (IsAUorGU(p3b1, p3b2))
                energy += TerminalAU_GU_Penalty;
        }

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the free energy of a multibranch loop using Turner 2004 parameters.
    /// ΔG° = a + b·avg_asymmetry + c·num_helices + optimal_stacking + strain
    /// </summary>
    /// <param name="numHelices">Number of helical branches in the multibranch loop.</param>
    /// <param name="numUnpaired">Total number of unpaired nucleotides in the loop.</param>
    /// <param name="hasStrain">True for 3-way junctions with fewer than 2 unpaired bases.</param>
    /// <param name="stackingEnergy">Pre-computed optimal stacking/dangling end energy for the junction.</param>
    public static double CalculateMultibranchLoopEnergy(
        int numHelices, int numUnpaired, bool hasStrain = false, double stackingEnergy = 0.0)
    {
        double asymmetry = numHelices > 0 ? (double)numUnpaired / numHelices : 0.0;
        double energy = MultibranchLoop_Offset
                      + MultibranchLoop_Asymmetry * asymmetry
                      + MultibranchLoop_Helix * numHelices
                      + stackingEnergy;

        if (hasStrain)
            energy += MultibranchLoop_Strain;

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates flush coaxial stacking energy — uses the WC/GU stacking table
    /// as if the two helices were continuous with no backbone break.
    /// </summary>
    public static double CalculateFlushCoaxialStacking(
        char base5_1, char base3_1, char base5_2, char base3_2)
    {
        string stack = $"{char.ToUpperInvariant(base5_1)}{char.ToUpperInvariant(base5_2)}/" +
                       $"{char.ToUpperInvariant(base3_1)}{char.ToUpperInvariant(base3_2)}";
        return StackingEnergies.TryGetValue(stack, out double e) ? Math.Round(e, 2) : 0.0;
    }

    /// <summary>
    /// Calculates mismatch-mediated coaxial stacking energy.
    /// ΔG° = terminal_mismatch + base_constant + WC/GU_bonus
    /// </summary>
    public static double CalculateMismatchCoaxialStacking(
        char closingBase5, char closingBase3, char mismatch5, char mismatch3)
    {
        double tmEnergy = GetTerminalMismatchEnergy(closingBase5, closingBase3, mismatch5, mismatch3);
        double energy = tmEnergy + CoaxialMismatch_Base;

        char c5 = char.ToUpperInvariant(closingBase5);
        char c3 = char.ToUpperInvariant(closingBase3);
        if ((c5 == 'G' && c3 == 'U') || (c5 == 'U' && c3 == 'G'))
            energy += CoaxialMismatch_GU_Bonus;
        else
            energy += CoaxialMismatch_WC_Bonus;

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the minimum free energy (MFE) of an RNA sequence using Zuker-style DP
    /// with Turner 2004 nearest-neighbor parameters from NNDB.
    /// Returns physical kcal/mol values.
    ///
    /// Recurrences:
    ///   W(j) = min{ W(j-1), min over 0≤i≤j { W(i-1) + V(i,j) + AU_penalty(i,j) + dangle(i,j) } }
    ///   V(i,j) = min{ Hairpin(i,j), Stack(i,j)+V(i+1,j-1), InternalOrBulge(i,j), Multi(i,j) }
    ///   WM(i,j) = min{ V(k,j)+b+dangle, WM(i,k-1)+V(k,j)+b, WM(i,j-1)+c, WM(i+1,j)+c }
    ///
    /// Complexity: O(n³) with MAXLOOP=30 for internal/bulge loops.
    /// </summary>
    public static double CalculateMinimumFreeEnergy(string rnaSequence, int minLoopSize = 3)
    {
        // NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops
        // with fewer than 3 nucleotides."
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minLoopSize + 2)
            return 0;

        // Accept DNA input by reading thymine as uracil: T and U are the same base for folding
        // (Watson–Crick A–U/A–T and stacking parameters are identical), so the energy model and the
        // A/C/G/U-indexed tables below treat the two interchangeably.
        string seq = rnaSequence.ToUpperInvariant().Replace('T', 'U');
        int n = seq.Length;
        const double INF = 1e18;
        const int MAXLOOP = 30;

        // Multibranch loop parameters (NNDB Turner 2004)
        const double ML_offset = 9.25;   // a — initiation
        const double ML_helix = -0.63;   // c — per helix
        const double ML_unpaired = 0.0;  // free base cost (simplified to 0)

        var pool = ArrayPool<double>.Shared;
        double[] vBuf = pool.Rent(n * n);  // V(i,j)
        double[] wmBuf = pool.Rent(n * n); // WM(i,j)
        double[] w = pool.Rent(n);          // W(j)

        try
        {
            // Initialize to +INF (no valid structure)
            Array.Fill(vBuf, INF, 0, n * n);
            Array.Fill(wmBuf, INF, 0, n * n);
            Array.Fill(w, 0.0, 0, n);

            // Fill diagonal and short spans
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < Math.Min(i + minLoopSize + 1, n); j++)
                {
                    // Too short to form any pair
                    vBuf[i * n + j] = INF;
                    wmBuf[i * n + j] = INF;
                }
            }

            // Main DP: fill by increasing span length
            for (int span = minLoopSize + 2; span <= n; span++)
            {
                for (int i = 0; i <= n - span; i++)
                {
                    int j = i + span - 1;
                    int ij = i * n + j;

                    // ===== V(i,j): i and j must pair =====
                    byte pairIJ = PairType(seq[i], seq[j]);
                    if (pairIJ != 0)
                    {
                        double vBest = INF;

                        // --- Option 1: Hairpin loop ---
                        int loopLen = j - i - 1;
                        if (loopLen >= minLoopSize)
                        {
                            string loopSeq = seq.Substring(i + 1, loopLen);

                            // Detect special GU closure: closing G-U preceded by two Gs on 5' side
                            // Source: NNDB hairpin-mismatch-parameters.html
                            bool specialGU = seq[i] == 'G' && seq[j] == 'U' &&
                                             i >= 2 && seq[i - 1] == 'G' && seq[i - 2] == 'G';

                            double hEnergy = CalculateHairpinLoopEnergy(loopSeq, seq[i], seq[j], specialGU);

                            // NNDB: terminal AU/GU penalty at inner helix end facing hairpin loop.
                            // Source: rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html
                            if (IsAUorGU(seq[i], seq[j]))
                                hEnergy += TerminalAU_GU_Penalty;

                            if (hEnergy < vBest) vBest = hEnergy;
                        }

                        // --- Option 2: Stacking (i,j) closes over (i+1,j-1) ---
                        if (j - i > 2 && PairType(seq[i + 1], seq[j - 1]) != 0)
                        {
                            double vInner = vBuf[(i + 1) * n + (j - 1)];
                            if (vInner < INF)
                            {
                                string stackKey = $"{seq[i]}{seq[i + 1]}/{seq[j]}{seq[j - 1]}";
                                if (StackingEnergies.TryGetValue(stackKey, out double stackE))
                                {
                                    double sEnergy = stackE + vInner;
                                    if (sEnergy < vBest) vBest = sEnergy;
                                }
                            }
                        }

                        // --- Option 2b: Special GGUC/CUGG 3-stack (NNDB note b) ---
                        // 5'GGUC/3'CUGG: total of 3 stacks = -4.12 replaces individual sum.
                        if (j - i >= 7 &&
                            seq[i] == 'G' && seq[j] == 'C' &&
                            seq[i + 1] == 'G' && seq[j - 1] == 'U' &&
                            PairType(seq[i + 1], seq[j - 1]) != 0 &&
                            seq[i + 2] == 'U' && seq[j - 2] == 'G' &&
                            PairType(seq[i + 2], seq[j - 2]) != 0 &&
                            seq[i + 3] == 'C' && seq[j - 3] == 'G' &&
                            PairType(seq[i + 3], seq[j - 3]) != 0)
                        {
                            double vInner3 = vBuf[(i + 3) * n + (j - 3)];
                            if (vInner3 < INF)
                            {
                                double sEnergy = SpecialGGUC_CUGG_3Stack + vInner3;
                                if (sEnergy < vBest) vBest = sEnergy;
                            }
                        }

                        // --- Option 3: Internal loop or bulge ---
                        // Enumerate (i', j') where i < i' < j' < j, (i',j') paired
                        // with loop size (i'-i-1) + (j-j'-1) ≤ MAXLOOP
                        for (int ip = i + 1; ip < j - 1 && ip - i - 1 <= MAXLOOP; ip++)
                        {
                            int maxN2 = MAXLOOP - (ip - i - 1);
                            int jpMin = Math.Max(ip + 1, j - maxN2 - 1);

                            for (int jp = j - 1; jp >= jpMin; jp--)
                            {
                                if (ip == i + 1 && jp == j - 1)
                                    continue; // Already handled by stacking

                                byte pairIpJp = PairType(seq[ip], seq[jp]);
                                if (pairIpJp == 0) continue;

                                double vInner = vBuf[ip * n + jp];
                                if (vInner >= INF) continue;

                                int n1 = ip - i - 1;
                                int n2 = j - jp - 1;

                                double loopE;
                                if (n1 == 0 || n2 == 0)
                                {
                                    // Bulge loop
                                    int bulgeSize = n1 + n2;
                                    char bulgedBase = n1 > 0 ? seq[i + 1] : seq[j - 1];

                                    // Degeneracy: count identical adjacent bases for n=1 bulges
                                    int numStates = 1;
                                    if (bulgeSize == 1)
                                    {
                                        if (n1 == 1)
                                        {
                                            // Bulge on 5' side at i+1; adjacent bases: seq[i], seq[ip]=seq[i+2]
                                            if (seq[i] == seq[i + 1]) numStates++;
                                            if (seq[ip] == seq[i + 1]) numStates++;
                                        }
                                        else
                                        {
                                            // Bulge on 3' side at j-1; adjacent bases: seq[j], seq[jp]=seq[j-2]
                                            if (seq[j] == seq[j - 1]) numStates++;
                                            if (seq[jp] == seq[j - 1]) numStates++;
                                        }
                                    }

                                    loopE = CalculateBulgeLoopEnergy(
                                        bulgeSize, bulgedBase,
                                        seq[i], seq[j], seq[ip], seq[jp],
                                        numStates);
                                }
                                else
                                {
                                    // Internal loop
                                    loopE = CalculateInternalLoopEnergy(
                                        n1, n2,
                                        seq[i], seq[j], seq[ip], seq[jp],
                                        seq[i + 1], seq[j - 1], seq[ip - 1], seq[jp + 1]);
                                }

                                double total = loopE + vInner;
                                if (total < vBest) vBest = total;
                            }
                        }

                        // --- Option 4: Multiloop ---
                        // V(i,j) = ML_offset + ML_helix + AU_penalty(i,j) + WM(i+1, j-1)
                        // We need at least 2 helices in WM, but approximate with WM decomposition
                        {
                            double wmInner = wmBuf[(i + 1) * n + (j - 1)];
                            if (wmInner < INF)
                            {
                                double mEnergy = ML_offset + ML_helix + wmInner;
                                if (IsAUorGU(seq[i], seq[j]))
                                    mEnergy += TerminalAU_GU_Penalty;
                                if (mEnergy < vBest) vBest = mEnergy;
                            }
                        }

                        vBuf[ij] = vBest;
                    }

                    // ===== WM(i,j): multiloop region =====
                    {
                        double wmBest = INF;

                        // Option A: V(i,j) starts a helix here
                        if (vBuf[ij] < INF)
                        {
                            double e = vBuf[ij] + ML_helix;
                            if (IsAUorGU(seq[i], seq[j]))
                                e += TerminalAU_GU_Penalty;
                            if (e < wmBest) wmBest = e;
                        }

                        // Option A2: i is 5' dangle, helix at (i+1, j)
                        if (i + 1 <= j)
                        {
                            double vij2 = vBuf[(i + 1) * n + j];
                            if (vij2 < INF)
                            {
                                double e = vij2 + ML_helix;
                                if (IsAUorGU(seq[i + 1], seq[j]))
                                    e += TerminalAU_GU_Penalty;
                                e += GetDanglingEndEnergy(seq[i + 1], seq[j], seq[i], false);
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        // Option A3: j is 3' dangle, helix at (i, j-1)
                        if (j - 1 >= i)
                        {
                            double vij3 = vBuf[i * n + (j - 1)];
                            if (vij3 < INF)
                            {
                                double e = vij3 + ML_helix;
                                if (IsAUorGU(seq[i], seq[j - 1]))
                                    e += TerminalAU_GU_Penalty;
                                e += GetDanglingEndEnergy(seq[i], seq[j - 1], seq[j], true);
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        // Option A4: both i and j are dangles, helix at (i+1, j-1)
                        if (i + 1 < j - 1)
                        {
                            double vij4 = vBuf[(i + 1) * n + (j - 1)];
                            if (vij4 < INF)
                            {
                                double e = vij4 + ML_helix;
                                if (IsAUorGU(seq[i + 1], seq[j - 1]))
                                    e += TerminalAU_GU_Penalty;
                                e += GetDanglingEndEnergy(seq[i + 1], seq[j - 1], seq[i], false);
                                e += GetDanglingEndEnergy(seq[i + 1], seq[j - 1], seq[j], true);
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        // Option B: i is unpaired in multiloop
                        if (i + 1 <= j)
                        {
                            double e = wmBuf[(i + 1) * n + j];
                            if (e < INF)
                            {
                                e += ML_unpaired;
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        // Option C: j is unpaired in multiloop
                        if (j - 1 >= i)
                        {
                            double e = wmBuf[i * n + (j - 1)];
                            if (e < INF)
                            {
                                e += ML_unpaired;
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        // Option D: split — WM(i,k) + WM(k+1,j), ensures ≥2 helices
                        for (int k = i; k < j; k++)
                        {
                            double left = wmBuf[i * n + k];
                            double right = wmBuf[(k + 1) * n + j];
                            if (left < INF && right < INF)
                            {
                                double e = left + right;
                                if (e < wmBest) wmBest = e;
                            }
                        }

                        wmBuf[ij] = wmBest;
                    }
                }

                // ===== W(j): overall MFE up to position j =====
                {
                    int j = span - 1; // only fill W for increasing j
                    if (j < minLoopSize + 1)
                    {
                        w[j] = 0; // too short
                        continue;
                    }

                    double wBest = w[j - 1]; // j unpaired

                    for (int i = 0; i <= j; i++)
                    {
                        double vij = vBuf[i * n + j];
                        if (vij >= INF) continue;

                        double auPenalty = IsAUorGU(seq[i], seq[j]) ? TerminalAU_GU_Penalty : 0;
                        double dangle = 0;
                        if (i > 0)
                            dangle += GetDanglingEndEnergy(seq[i], seq[j], seq[i - 1], false);
                        if (j < n - 1)
                            dangle += GetDanglingEndEnergy(seq[i], seq[j], seq[j + 1], true);
                        double left = i > 0 ? w[i - 1] : 0;
                        double total = left + vij + auPenalty + dangle;
                        if (total < wBest) wBest = total;
                    }

                    w[j] = wBest;
                }
            }

            // Final pass: W[j] = min { W[j-1], min_i { W[i-1] + V(i,j) + AU + dangle } }
            for (int j = 0; j < n; j++)
            {
                if (j < minLoopSize + 1)
                {
                    w[j] = 0;
                    continue;
                }

                double wBest = w[j - 1];
                for (int i = 0; i <= j; i++)
                {
                    double vij = vBuf[i * n + j];
                    if (vij >= INF) continue;

                    double auPenalty = IsAUorGU(seq[i], seq[j]) ? TerminalAU_GU_Penalty : 0;
                    double dangle = 0;
                    if (i > 0)
                        dangle += GetDanglingEndEnergy(seq[i], seq[j], seq[i - 1], false);
                    if (j < n - 1)
                        dangle += GetDanglingEndEnergy(seq[i], seq[j], seq[j + 1], true);
                    double left = i > 0 ? w[i - 1] : 0;
                    double total = left + vij + auPenalty + dangle;
                    if (total < wBest) wBest = total;
                }

                w[j] = wBest;
            }

            double result = w[n - 1];
            return Math.Round(result, 2);
        }
        finally
        {
            pool.Return(vBuf);
            pool.Return(wmBuf);
            pool.Return(w);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte PairType(char b1, char b2)
    {
        return (b1, b2) switch
        {
            ('A', 'U') or ('U', 'A') or ('G', 'C') or ('C', 'G') => 1,
            ('G', 'U') or ('U', 'G') => 2,
            _ => 0,
        };
    }

    /// <summary>
    /// Original O(L\u00b3) MFE — retained as benchmark baseline.
    /// </summary>
    internal static double CalculateMinimumFreeEnergyClassic(string rnaSequence, int minLoopSize = 3)
    {
        // NNDB Turner 2004: hairpin loops < 3 nt are prohibited.
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minLoopSize + 2)
            return 0;

        // Accept DNA input by reading thymine as uracil: T and U are the same base for folding
        // (Watson–Crick A–U/A–T and stacking parameters are identical), so the energy model and the
        // A/C/G/U-indexed tables below treat the two interchangeably.
        string seq = rnaSequence.ToUpperInvariant().Replace('T', 'U');
        int n = seq.Length;

        var dp = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                dp[i, j] = 0;

        for (int length = minLoopSize + 2; length <= n; length++)
        {
            for (int i = 0; i <= n - length; i++)
            {
                int j = i + length - 1;
                dp[i, j] = dp[i, j - 1];

                for (int k = i; k < j - minLoopSize; k++)
                {
                    if (CanPair(seq[k], seq[j]))
                    {
                        double pairEnergy = GetPairEnergy(seq[k], seq[j]);
                        double left = k > i ? dp[i, k - 1] : 0;
                        double enclosed = k + 1 < j ? dp[k + 1, j - 1] : 0;
                        double total = pairEnergy + left + enclosed;
                        dp[i, j] = Math.Min(dp[i, j], total);
                    }
                }
            }
        }

        return Math.Round(dp[0, n - 1], 2);
    }

    private static double GetPairEnergy(char base1, char base2)
    {
        var pairType = GetBasePairType(base1, base2);
        return pairType switch
        {
            BasePairType.WatsonCrick => -2.0,
            BasePairType.Wobble => -1.0,
            _ => 0
        };
    }

    #endregion

    #region Structure Prediction

    /// <summary>
    /// Predicts the secondary structure of an RNA sequence.
    /// </summary>
    public static SecondaryStructure PredictStructure(
        string rnaSequence,
        int minStemLength = 3,
        int minLoopSize = 3,
        int maxLoopSize = 10)
    {
        if (string.IsNullOrEmpty(rnaSequence))
        {
            return new SecondaryStructure(
                "", "", new List<BasePair>(), new List<StemLoop>(),
                new List<Pseudoknot>(), 0);
        }

        string seq = rnaSequence.ToUpperInvariant();

        // Find stem-loops
        var stemLoops = FindStemLoops(seq, minStemLength, minLoopSize, maxLoopSize)
            .OrderBy(sl => sl.TotalFreeEnergy)
            .ToList();

        // Select non-overlapping structures greedily by energy
        var selectedStemLoops = SelectNonOverlapping(stemLoops);

        // Collect all base pairs
        var allBasePairs = selectedStemLoops
            .SelectMany(sl => sl.Stem.BasePairs)
            .OrderBy(bp => bp.Position1)
            .ToList();

        // Generate dot-bracket notation
        string dotBracket = GenerateFullDotBracket(seq.Length, allBasePairs);

        // Calculate total MFE
        double mfe = selectedStemLoops.Sum(sl => sl.TotalFreeEnergy);

        // Detect pseudoknots
        var pseudoknots = DetectPseudoknots(allBasePairs).ToList();

        return new SecondaryStructure(
            Sequence: seq,
            DotBracket: dotBracket,
            BasePairs: allBasePairs,
            StemLoops: selectedStemLoops,
            Pseudoknots: pseudoknots,
            MinimumFreeEnergy: mfe);
    }

    private static List<StemLoop> SelectNonOverlapping(List<StemLoop> stemLoops)
    {
        var selected = new List<StemLoop>();
        var used = new HashSet<int>();

        foreach (var sl in stemLoops)
        {
            bool overlaps = false;
            for (int pos = sl.Start; pos <= sl.End; pos++)
            {
                if (used.Contains(pos))
                {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps)
            {
                selected.Add(sl);
                for (int pos = sl.Start; pos <= sl.End; pos++)
                {
                    used.Add(pos);
                }
            }
        }

        return selected;
    }

    #endregion

    #region Pseudoknot Detection

    /// <summary>
    /// Detects pseudoknots in a set of base pairs.
    /// Two base pairs (i,j) and (k,l) — written with open &lt; close positions — form a
    /// pseudoknot when they CROSS, i.e. i &lt; k &lt; j &lt; l. Nested pairs (i &lt; k &lt; l &lt; j) and
    /// disjoint pairs (j &lt; k) are not pseudoknots. Each crossing pair-of-pairs is reported
    /// as one <see cref="Pseudoknot"/>.
    /// </summary>
    /// <remarks>
    /// Crossing condition per Antczak et al. (2018) Bioinformatics 34(8):1304–1312 and
    /// biotite.structure.pseudoknots: a pseudoknot is a non-nested ("crossing") arrangement
    /// of base pairs. Endpoints are normalized to (open &lt; close) before the comparison.
    /// </remarks>
    /// <param name="basePairs">The base pairs to scan. A <c>null</c> set yields no pseudoknots.</param>
    /// <returns>One <see cref="Pseudoknot"/> per crossing pair-of-pairs, in deterministic input order.</returns>
    public static IEnumerable<Pseudoknot> DetectPseudoknots(IReadOnlyList<BasePair> basePairs)
    {
        // A pseudoknot requires two base pairs that cross; fewer than two cannot.
        if (basePairs is null || basePairs.Count < 2)
            yield break;

        for (int a = 0; a < basePairs.Count; a++)
        {
            for (int b = a + 1; b < basePairs.Count; b++)
            {
                var bp1 = basePairs[a];
                var bp2 = basePairs[b];

                // Normalize each pair to (open < close) so storage order does not matter.
                int i = Math.Min(bp1.Position1, bp1.Position2);
                int j = Math.Max(bp1.Position1, bp1.Position2);
                int k = Math.Min(bp2.Position1, bp2.Position2);
                int l = Math.Max(bp2.Position1, bp2.Position2);

                // Order the two pairs by their opening position so the crossing test is
                // direction-independent: the pair that opens first plays the role of (i,j).
                if (k < i)
                {
                    (i, j, k, l) = (k, l, i, j);
                }

                // Crossing (pseudoknot) condition: i < k < j < l.
                if (i < k && k < j && j < l)
                {
                    yield return new Pseudoknot(
                        Start1: i,
                        End1: j,
                        Start2: k,
                        End2: l,
                        CrossingPairs: new List<BasePair> { bp1, bp2 });
                }
            }
        }
    }

    #endregion

    #region Dot-Bracket Notation

    /// <summary>
    /// Generates dot-bracket notation for a structure.
    /// </summary>
    private static string GenerateDotBracket(int length, IReadOnlyList<BasePair> basePairs, int start, int end)
    {
        var notation = new char[end - start + 1];
        for (int i = 0; i < notation.Length; i++)
            notation[i] = '.';

        foreach (var bp in basePairs)
        {
            if (bp.Position1 >= start && bp.Position1 <= end)
                notation[bp.Position1 - start] = '(';
            if (bp.Position2 >= start && bp.Position2 <= end)
                notation[bp.Position2 - start] = ')';
        }

        return new string(notation);
    }

    private static string GenerateFullDotBracket(int length, IReadOnlyList<BasePair> basePairs)
    {
        var notation = new char[length];
        for (int i = 0; i < length; i++)
            notation[i] = '.';

        foreach (var bp in basePairs)
        {
            int left = Math.Min(bp.Position1, bp.Position2);
            int right = Math.Max(bp.Position1, bp.Position2);
            notation[left] = '(';
            notation[right] = ')';
        }

        return new string(notation);
    }

    // Dot-bracket / extended WUSS bracket pairs. Each opening symbol pairs ONLY with
    // its matching closing symbol; the four bracket families are independent pairing
    // systems so crossing (pseudoknotted) helices can be annotated.
    // Sources: ViennaRNA RNA Structure Notations — nested base pairs annotated by matching
    //   pairs of <> () {} []; Infernal/WUSS — "any matching nested pair of (), [], {}
    //   indicates a base pair ... so long as the left and right partners match up".
    private static readonly IReadOnlyDictionary<char, char> ClosingToOpening = new Dictionary<char, char>
    {
        [')'] = '(',
        [']'] = '[',
        ['}'] = '{',
        ['>'] = '<',
    };

    private static readonly IReadOnlyDictionary<char, char> OpeningToClosing = new Dictionary<char, char>
    {
        ['('] = ')',
        ['['] = ']',
        ['{'] = '}',
        ['<'] = '>',
    };

    /// <summary>
    /// Parses dot-bracket / extended WUSS notation to extract base pairs as
    /// (5' position, 3' position) index tuples (zero-based).
    /// </summary>
    /// <remarks>
    /// Each bracket family (<c>()</c>, <c>[]</c>, <c>{}</c>, <c>&lt;&gt;</c>) is matched on
    /// its own stack, and uppercase/lowercase letter pairs (<c>Aa</c>, <c>Bb</c>, ...) are
    /// matched as independent pseudoknot systems, per ViennaRNA / WUSS conventions.
    /// An opening symbol is paired with the nearest unmatched opening of the SAME family
    /// (correct nesting). Dots and any other symbols denote unpaired residues and are skipped.
    /// Closing symbols with no open partner are ignored (parse is best-effort; use
    /// <see cref="ValidateDotBracket"/> to test well-formedness first).
    /// </remarks>
    public static IEnumerable<(int Position1, int Position2)> ParseDotBracket(string dotBracket)
    {
        if (string.IsNullOrEmpty(dotBracket))
            yield break;

        // One independent stack per opening symbol (brackets + uppercase letters).
        var stacks = new Dictionary<char, Stack<int>>();

        for (int i = 0; i < dotBracket.Length; i++)
        {
            char c = dotBracket[i];

            if (OpeningToClosing.ContainsKey(c) || (char.IsLetter(c) && char.IsUpper(c)))
            {
                if (!stacks.TryGetValue(c, out var stack))
                {
                    stack = new Stack<int>();
                    stacks[c] = stack;
                }
                stack.Push(i);
            }
            else if (ClosingToOpening.TryGetValue(c, out char opener))
            {
                if (stacks.TryGetValue(opener, out var stack) && stack.Count > 0)
                    yield return (stack.Pop(), i);
            }
            else if (char.IsLetter(c) && char.IsLower(c))
            {
                char up = char.ToUpperInvariant(c);
                if (stacks.TryGetValue(up, out var stack) && stack.Count > 0)
                    yield return (stack.Pop(), i);
            }
            // '.', ',', '-', ':', '_', '~' and any other symbol => unpaired, skip.
        }
    }

    /// <summary>
    /// Validates dot-bracket / extended WUSS notation: every closing symbol must match an
    /// earlier unmatched opening symbol of the SAME family, and no opening symbol may be left
    /// unclosed. Each bracket family and each letter (case) pair is checked independently, so a
    /// mismatched pairing such as <c>(]</c> is rejected even though the total count is balanced.
    /// </summary>
    /// <remarks>
    /// Sources: ViennaRNA — base pairs are matching pairs of <c>()</c> (and extended <c>[]{}&lt;&gt;</c>,
    /// letters) and the structure must be balanced; WUSS — left and right partners must match up.
    /// </remarks>
    public static bool ValidateDotBracket(string dotBracket)
    {
        if (string.IsNullOrEmpty(dotBracket))
            return true;

        var stacks = new Dictionary<char, Stack<int>>();

        foreach (char c in dotBracket)
        {
            if (OpeningToClosing.ContainsKey(c) || (char.IsLetter(c) && char.IsUpper(c)))
            {
                if (!stacks.TryGetValue(c, out var stack))
                {
                    stack = new Stack<int>();
                    stacks[c] = stack;
                }
                stack.Push(0);
            }
            else if (ClosingToOpening.TryGetValue(c, out char opener))
            {
                if (!stacks.TryGetValue(opener, out var stack) || stack.Count == 0)
                    return false;
                stack.Pop();
            }
            else if (char.IsLetter(c) && char.IsLower(c))
            {
                char up = char.ToUpperInvariant(c);
                if (!stacks.TryGetValue(up, out var stack) || stack.Count == 0)
                    return false;
                stack.Pop();
            }
            // unpaired symbol => ignore.
        }

        foreach (var stack in stacks.Values)
            if (stack.Count != 0)
                return false;
        return true;
    }

    #endregion

    #region Utility Methods

    // Default parameters for inverted-repeat detection.
    // Minimum arm length: an arm shorter than this is not reported (IUPACpal "minimum length").
    // Source: Alamro et al. (2021) IUPACpal, BMC Bioinformatics 22:51. https://doi.org/10.1186/s12859-021-03983-2
    private const int DefaultMinArmLength = 4;
    // Minimum loop length: a stable RNA hairpin loop needs at least 3 unpaired nucleotides.
    // Source: IUPACpal gap |G| >= 0; EMBOSS einverted loop. 3-nt minimum is the repository RNA convention.
    private const int DefaultMinLoopLength = 3;
    // Maximum loop length searched (bounds the gap G; IUPACpal "maximum gap size").
    private const int DefaultMaxLoopLength = 100;

    /// <summary>
    /// Finds inverted repeats (potential RNA stem regions): a left arm <c>W</c> followed by a
    /// loop and then a right arm equal to the <b>reverse complement</b> of <c>W</c>, i.e. the
    /// pattern <c>W G W̄ᴿ</c> with a perfect (zero-mismatch) antiparallel stem. The two arms form
    /// the stem of a potential hairpin and the intervening region is the loop.
    /// </summary>
    /// <param name="sequence">Nucleotide sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="minLength">Minimum arm length; arms shorter than this are not reported (default 4).</param>
    /// <param name="minSpacing">Minimum loop length between the two arms (default 3).</param>
    /// <param name="maxSpacing">Maximum loop length between the two arms (default 100).</param>
    /// <returns>
    /// Maximal, non-overlapping inverted repeats as tuples: left arm <c>[Start1..End1]</c>,
    /// right arm <c>[Start2..End2]</c> (both 0-based inclusive), and the common arm <c>Length</c>.
    /// For every reported repeat and every k in [0,Length), the antiparallel pairing
    /// <c>complement(seq[Start2+Length-1-k]) == seq[Start1+k]</c> holds (strict Watson-Crick/IUPAC).
    /// </returns>
    /// <remarks>
    /// Sources: an inverted repeat is "a single stranded sequence ... followed downstream by its
    /// reverse complement", formally <c>WGW̄ᴿ</c> with |G| >= 0 (Alamro et al. 2021, IUPACpal,
    /// https://doi.org/10.1186/s12859-021-03983-2; Ussery et al. 2008 via Wikipedia). This is the
    /// k = 0 (perfect) case; the scored mismatch/gap variant (EMBOSS einverted) is out of scope.
    /// Search reuse: the repository SuffixTree was evaluated and not used — matching here is
    /// antiparallel reverse-complement extension over a bounded loop window, not exact substring
    /// occurrence enumeration, which is the established approach (EMBOSS einverted uses local DP).
    /// </remarks>
    public static IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)> FindInvertedRepeats(
        string sequence,
        int minLength = DefaultMinArmLength,
        int minSpacing = DefaultMinLoopLength,
        int maxSpacing = DefaultMaxLoopLength)
    {
        if (string.IsNullOrEmpty(sequence) || minLength < 1 || minSpacing < 0 || maxSpacing < minSpacing)
            yield break;

        string upper = sequence.ToUpperInvariant();
        int n = upper.Length;
        // Smallest sequence that can hold two arms of minLength plus the minimum loop.
        if (n < 2 * minLength + minSpacing)
            yield break;

        // Scan loop boundaries: p = last index of the left arm, q = first index of the right arm,
        // with loop length (q - p - 1) in [minSpacing, maxSpacing]. The innermost base pair is
        // (seq[p], seq[q]); the stem extends OUTWARD: seq[p-k] pairs antiparallel with seq[q+k].
        // This realizes the W G W̄ᴿ pattern where the right arm is the reverse complement of the left.
        // Collect the best (maximal) stem per loop boundary, then report non-overlapping repeats
        // preferring the longest stems first (einverted-style maximal local alignment).
        var candidates = new List<(int Start1, int End1, int Start2, int End2, int Length)>();
        for (int p = 0; p < n; p++)
        {
            int firstQ = p + 1 + minSpacing;
            int lastQ = Math.Min(p + 1 + maxSpacing, n - 1);
            int bestLen = 0;
            int bestRightEnd = -1;
            for (int q = firstQ; q <= lastQ; q++)
            {
                // Extend the stem outward while complementary; arms may not cross into the loop
                // (left stays >= 0, right stays <= n-1, and they may not overlap each other).
                int maxExtend = Math.Min(p + 1, n - q);
                int len = 0;
                while (len < maxExtend &&
                       GetComplement(upper[q + len]) == upper[p - len])
                {
                    len++;
                }

                if (len < minLength)
                    continue;

                // Prefer the longest arm; on ties prefer the smaller loop (q closer to p).
                if (len > bestLen)
                {
                    bestLen = len;
                    bestRightEnd = q + len - 1;
                }
            }

            if (bestLen >= minLength)
                candidates.Add((p - bestLen + 1, p, bestRightEnd - bestLen + 1, bestRightEnd, bestLen));
        }

        // Greedy non-overlapping selection: longest stems first; ties broken by left position then loop.
        candidates.Sort((a, b) =>
        {
            int c = b.Length.CompareTo(a.Length);
            if (c != 0) return c;
            c = a.Start1.CompareTo(b.Start1);
            if (c != 0) return c;
            return a.Start2.CompareTo(b.Start2);
        });

        var occupied = new bool[n];
        foreach (var cand in candidates)
        {
            bool overlaps = false;
            for (int idx = cand.Start1; idx <= cand.End1 && !overlaps; idx++)
                if (occupied[idx]) overlaps = true;
            for (int idx = cand.Start2; idx <= cand.End2 && !overlaps; idx++)
                if (occupied[idx]) overlaps = true;
            if (overlaps)
                continue;

            for (int idx = cand.Start1; idx <= cand.End1; idx++) occupied[idx] = true;
            for (int idx = cand.Start2; idx <= cand.End2; idx++) occupied[idx] = true;
            yield return cand;
        }
    }

    // McCaskill partition function — physical constants.
    // R = 1.987 cal/(mol·K) (molar gas constant); RT expressed in kcal/mol matches
    // the kcal/mol energy unit used by the simplified model below.
    // ViennaRNA documents the same Boltzmann constant k ≈ 1.987e-3 kcal/(mol·K)
    // and p(s) = exp(−βE(s))/Z with β = 1/kT.  Source: ViennaRNA pf_fold reference;
    // McCaskill JS (1990) Biopolymers 29:1105-1119 (PMID 1695107).
    private const double GasConstant_CalPerMolK = 1.987;
    // Default folding temperature 37 °C = 310.15 K (ViennaRNA / NNDB default).
    private const double DefaultTemperatureKelvin = 310.15;
    // Minimum number of unpaired bases enclosed by a hairpin loop. A pair (i,j) is
    // sterically admissible only when j − i > MinHairpinLoop. The classic value is 3
    // (Nussinov 1980 / standard nearest-neighbour models forbid loops shorter than 3 nt).
    private const int MinHairpinLoop = 3;
    // Per-base-pair free-energy increment for the SIMPLIFIED energy model (kcal/mol).
    // Each base pair contributes a fixed term E_bp (Freiburg RNA teaching model:
    // "each base pair of a structure contributes a fixed energy term E_bp").
    private const double SimplifiedBasePairEnergy = -1.0;

    /// <summary>
    /// Computes the McCaskill equilibrium partition function Z and the equilibrium
    /// base-pair probabilities P[i,j] over the full ensemble of pseudoknot-free
    /// secondary structures of the given RNA sequence.
    /// </summary>
    /// <remarks>
    /// Implements the McCaskill (1990) inside recursion Q/Q^b with Z = Q₁ₙ in O(n³) time,
    /// O(n²) space, followed by the external base-pair-probability formula
    /// P[i,j] = Q(1,i−1)·Q^b(i,j)·Q(j+1,n)/Z. The energy model is the simplified
    /// fixed-per-pair model (each pair contributes <c>basePairEnergy</c>); this is the
    /// model documented by the Freiburg RNA McCaskill teaching tool. The well-defined
    /// partition-function mathematics (Z ≥ 1, probabilities in [0,1], normalisation,
    /// symmetry) hold independently of the energy parameter.
    /// </remarks>
    /// <param name="sequence">RNA sequence (A/C/G/U, case-insensitive; T is treated as U).</param>
    /// <param name="basePairEnergy">Fixed free energy E_bp (kcal/mol) per base pair.</param>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not positive.</exception>
    public static PartitionFunctionResult CalculatePartitionFunction(
        string sequence,
        double basePairEnergy = SimplifiedBasePairEnergy,
        double temperature = DefaultTemperatureKelvin)
    {
        if (sequence is null)
            throw new ArgumentNullException(nameof(sequence));
        if (temperature <= 0)
            throw new ArgumentOutOfRangeException(nameof(temperature), "Temperature must be positive (Kelvin).");

        int n = sequence.Length;
        var empty = new Dictionary<(int, int), double>();
        if (n == 0)
            return new PartitionFunctionResult(1.0, empty); // empty sequence: only the empty structure, Z = 1.

        string s = sequence.ToUpperInvariant();

        // RT in kcal/mol, β = 1/RT; Boltzmann weight of one base pair = exp(−β·E_bp).
        double rt = GasConstant_CalPerMolK * temperature / 1000.0;
        double pairWeight = Math.Exp(-basePairEnergy / rt);

        // Inside matrices, 0-based inclusive sub-sequences [i..j].
        // Q[i,j]  = partition function of [i..j]; empty interval (i > j) is defined as 1.
        // Qb[i,j] = partition function of [i..j] restricted to structures with (i,j) paired.
        var q = new double[n, n];
        var qb = new double[n, n];

        double GetQ(int i, int j) => i > j ? 1.0 : q[i, j];

        // Fill by increasing sub-sequence length so dependencies are already computed.
        for (int length = 1; length <= n; length++)
        {
            for (int i = 0; i + length - 1 < n; i++)
            {
                int j = i + length - 1;

                // Q^b recursion (simplified model): a closing pair contributes pairWeight
                // times the partition function of the enclosed interval.
                if (j - i > MinHairpinLoop && CanPair(s[i], s[j]))
                    qb[i, j] = pairWeight * GetQ(i + 1, j - 1);
                else
                    qb[i, j] = 0.0;

                // Q recursion: Q[i,j] = Q[i,j-1] + Σ_k Q[i,k-1]·Q^b[k,j]
                // (McCaskill / Will 18.417: structures where j is unpaired, plus structures
                // where j pairs with some k, splitting the interval unambiguously).
                double val = GetQ(i, j - 1);
                for (int k = i; k <= j; k++)
                {
                    if (j - k > MinHairpinLoop && CanPair(s[k], s[j]))
                        val += GetQ(i, k - 1) * qb[k, j];
                }
                q[i, j] = val;
            }
        }

        double z = GetQ(0, n - 1); // Z = Q(1,n) in 1-based notation.

        // McCaskill (1990) base-pair probabilities via the outside recursion.
        // p[i,j] = Q^b(i,j)·O(i,j)/Z, where O(i,j) is the OUTSIDE partition function: the sum
        // of Boltzmann weights of every structure of the parts of the sequence that lie outside
        // [i..j] (and of every pair that encloses (i,j)). It obeys
        //     O(i,j) = Q(0,i-1)·Q(j+1,n-1)                                  (external: (i,j) unenclosed)
        //            + Σ_{k<i, l>j, CanPair(k,l), l-k>m}  w·Q(k+1,i-1)·Q(j+1,l-1)·O(k,l)
        // where w = exp(−β·E_bp) is the per-pair Boltzmann weight: an enclosing pair (k,l)
        // contributes its own weight w, the two gap regions (k+1..i-1) and (j+1..l-1) fold
        // freely (Q), and the structure outside (k,l) is O(k,l) (computed recursively).
        // Using only the external term (the first line) is WRONG for any pair that can be
        // nested inside another pair — it omits the enclosing contributions and under-reports
        // those pairs. Sources: McCaskill JS (1990) Biopolymers 29:1105-1119; Will S., MIT
        // 18.417 McCaskill base-pair-probability slides (mccaskill2.pdf); ViennaRNA pf_fold.
        var probabilities = new Dictionary<(int I, int J), double>();
        var outside = new double[n, n];
        // Probabilities P[i,j] = Qᵇ·O/Z are well-defined only when Z is a FINITE positive
        // number. At sub-Kelvin temperatures RT → 0 and the per-pair Boltzmann weight
        // exp(−βE_bp) overflows IEEE double to +∞, making Z = +∞; computing Qᵇ·O/Z would then
        // yield ∞/∞ = NaN and leak it into the probability map. Require finiteness so an
        // overflowed ensemble returns an empty (non-NaN) probability map rather than NaN mass.
        if (double.IsFinite(z) && z > 0)
        {
            // Collect admissible pairs and process them outermost-first (decreasing span)
            // so that every enclosing pair's outside value O(k,l) is already known.
            var pairs = new List<(int I, int J)>();
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (qb[i, j] > 0)
                        pairs.Add((i, j));
            pairs.Sort((a, b) => (b.J - b.I).CompareTo(a.J - a.I));

            foreach (var (i, j) in pairs)
            {
                // External term: (i,j) is not enclosed by any other pair.
                double o = GetQ(0, i - 1) * GetQ(j + 1, n - 1);

                // Enclosing terms: every pair (k,l) with k<i, j<l that can directly enclose (i,j).
                for (int k = 0; k < i; k++)
                {
                    for (int l = j + 1; l < n; l++)
                    {
                        if (qb[k, l] > 0)
                        {
                            o += pairWeight * GetQ(k + 1, i - 1) * GetQ(j + 1, l - 1) * outside[k, l];
                        }
                    }
                }

                outside[i, j] = o;
                probabilities[(i, j)] = qb[i, j] * o / z;
            }
        }

        return new PartitionFunctionResult(z, probabilities);
    }

    /// <summary>
    /// Calculates the Boltzmann probability of a single structure within the ensemble:
    /// p(S) = exp(−E_S/RT) / exp(−E_ensemble/RT), where E_ensemble = −RT·ln Z is the
    /// ensemble free energy. Equivalent to the McCaskill normalisation p(S) = exp(−E_S/RT)/Z.
    /// </summary>
    /// <remarks>
    /// Source: McCaskill JS (1990) Biopolymers 29:1105-1119, p(P|S) = Z⁻¹·exp(−βE(P));
    /// ViennaRNA pf_fold reference (β = 1/kT, k ≈ 1.987e-3 kcal/(mol·K)).
    /// </remarks>
    /// <param name="temperature">Absolute temperature in Kelvin (default 310.15 K = 37 °C).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not positive.</exception>
    public static double CalculateStructureProbability(double structureEnergy, double ensembleEnergy, double temperature = DefaultTemperatureKelvin)
    {
        // RT appears in the denominator (β = 1/RT); a non-positive Kelvin temperature is
        // physically meaningless and would make RT = 0 (→ exp(±∞) → NaN) or RT < 0 (an
        // inverted Boltzmann weight). Reject it rather than leak a NaN/Inf probability —
        // mirrors CalculatePartitionFunction's temperature guard for caller safety.
        if (temperature <= 0)
            throw new ArgumentOutOfRangeException(nameof(temperature), "Temperature must be positive (Kelvin).");

        double rt = GasConstant_CalPerMolK * temperature / 1000.0; // kcal/mol

        double boltzmann = Math.Exp(-structureEnergy / rt);
        double partition = Math.Exp(-ensembleEnergy / rt);

        return partition > 0 ? boltzmann / partition : 0;
    }

    /// <summary>
    /// Generates a random RNA sequence with the given target GC content.
    /// Uses a non-deterministic seed; for reproducible output use the seeded overload.
    /// </summary>
    public static string GenerateRandomRna(int length, double gcContent = 0.5)
        => GenerateRandomRna(length, new Random(), gcContent);

    /// <summary>
    /// Generates a random RNA sequence with the given target GC content using a
    /// caller-supplied <see cref="Random"/> instance (deterministic when seeded).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="random"/> is null.</exception>
    public static string GenerateRandomRna(int length, Random random, double gcContent = 0.5)
    {
        if (random is null)
            throw new ArgumentNullException(nameof(random));
        var sb = new StringBuilder(length);

        for (int i = 0; i < length; i++)
        {
            double r = random.NextDouble();
            if (r < gcContent / 2)
                sb.Append('G');
            else if (r < gcContent)
                sb.Append('C');
            else if (r < gcContent + (1 - gcContent) / 2)
                sb.Append('A');
            else
                sb.Append('U');
        }

        return sb.ToString();
    }

    #endregion
}
