using System.Buffers;
using System.Runtime.CompilerServices;
using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides algorithms for RNA secondary structure prediction and analysis.
/// Includes stem-loop detection, hairpin finding, and free energy calculations.
/// </summary>
public static partial class RnaSecondaryStructure
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
    /// One crossing (pseudoknotted) pair-of-pairs reported by <see cref="DetectPseudoknots"/>: the pairs
    /// (Start1, End1) and (Start2, End2), normalised to open &lt; close, with Start1 &lt; Start2 &lt; End1 &lt; End2.
    /// Grouping into pseudoknot layers: <see cref="AssignPseudoknotOrders"/>.
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
    /// The minimum-free-energy (MFE) secondary structure recovered by Zuker–Stiegler (1981)
    /// dynamic-programming traceback over the same V/W/WM matrices used to compute the scalar MFE.
    /// </summary>
    /// <param name="Sequence">The folded sequence (upper-cased; T read as U), as scored by the DP.</param>
    /// <param name="DotBracket">
    /// Dot-bracket notation of the optimal structure: <c>(</c>/<c>)</c> for the two ends of every
    /// base pair, <c>.</c> for unpaired positions. Always balanced and pseudoknot-free.
    /// </param>
    /// <param name="BasePairs">
    /// The optimal base pairs as 0-based (5' &lt; 3') index tuples, sorted by 5' position.
    /// </param>
    /// <param name="FreeEnergy">
    /// ΔG° (kcal/mol) of <see cref="DotBracket"/>. By construction this equals the scalar value
    /// returned by <see cref="CalculateMinimumFreeEnergy"/> for the same input and parameters.
    /// </param>
    public readonly record struct MfeStructure(
        string Sequence,
        string DotBracket,
        IReadOnlyList<(int Position1, int Position2)> BasePairs,
        double FreeEnergy);

    #endregion

    #region Free Energy Parameters

    // Turner 2004 nearest-neighbor parameters (kcal/mol at 37°C).
    // The shared tables live in the authoritative single source Seqeron.Genomics.Core.
    // Turner2004Parameters (+ Turner2004InteriorLoopTables / Turner2004Enthalpies) so this folder and
    // MiRnaAnalyzer can never drift apart. Every loop term is evaluated by the one array-backed
    // TurnerModel (at 37 °C or rescaled to another temperature); the dictionaries below are only
    // used by the closed-form NNDB helpers (stem, coaxial stacking) that work on characters.
    private static readonly IReadOnlyDictionary<string, double> StackingEnergies =
        Turner2004Parameters.StackingEnergies;

    private static readonly IReadOnlyDictionary<string, double> TerminalMismatchEnergies =
        Turner2004Parameters.TerminalMismatchEnergies;

    // Terminal AU/GU penalty (kcal/mol at 37°C) — see Turner2004Parameters.
    private const double TerminalAU_GU_Penalty = Turner2004Parameters.TerminalAuGuPenalty;

    // Special GGUC/CUGG 3-stack context (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/gu-parameters.html (note b)
    // GU followed by UG is generally unfavorable (+1.29), but when preceded by GC and
    // followed by CG (5'GGUC/3'CUGG), the total for 3 stacks = -4.12 replaces individual sum.
    private const double SpecialGGUC_CUGG_3Stack = TurnerModel.SpecialGgucCugg3Stack;

    // Multibranch loop initiation, NNDB Turner 2004 (efn2) model (kcal/mol at 37°C).
    // Source: NNDB rna_2004_multibranch_loops.html §9.1 + parameter table turner_2004_multi_branch:
    // ΔG°init = a + b×[average asymmetry] + c×[number of branching helices] + ΔG°strain,
    // average asymmetry = min[2.0, mean |unpaired 5' side − unpaired 3' side| over the helices];
    // strain only for three-way loops with fewer than two unpaired nucleotides.
    private const double MultibranchLoop_Offset = 9.25;     // a
    private const double MultibranchLoop_Asymmetry = 0.91;  // b
    private const double MultibranchLoop_Helix = -0.63;     // c
    private const double MultibranchLoop_Strain = 3.14;     // 3-way with <2 unpaired
    private const double MultibranchLoop_MaxAsymmetry = 2.0;

    // Coaxial stacking parameters (kcal/mol at 37°C)
    // Source: NNDB rna_2004_coaxial_stacking.html §8.2–8.3 + table turner_2004_coax:
    // flush: the helical stack as if there were no backbone break; mismatch-mediated: terminal mismatch
    // (continuous backbone) + −2.1 (discontinuous backbone) + −0.4 / −0.2 when the "mismatch" could form a
    // Watson–Crick / G·U pair.
    private const double CoaxialMismatch_Base = -2.1;
    private const double CoaxialMismatch_WC_Bonus = -0.4;
    private const double CoaxialMismatch_GU_Bonus = -0.2;

    // Fast pair-type lookup: 0 = no pair, 1 = Watson-Crick, 2 = Wobble.
    // Indexed by [base1 * 128 + base2] over UPPERCASE ASCII. 16 KB — fits L1 cache.
    // Canonical pair set of the ViennaRNA default model (md.pair / BP_ENCODING_DEFAULT in
    // src/ViennaRNA/model.c): CG, GC, GU, UG, AU, UA. DNA thymine is the same base as uracil for
    // pairing (ViennaRNA vrna_nucleotide_encode: "make T and U equivalent"), so A·T/T·A are
    // Watson-Crick and G·T/T·G wobble; every other character (IUPAC codes, N, gaps) does not pair.
    private static readonly byte[] PairLookup = BuildPairLookup();

    private static byte[] BuildPairLookup()
    {
        var t = new byte[128 * 128];
        foreach (char u in "UT")
        {
            t['A' * 128 + u] = 1; t[u * 128 + 'A'] = 1;
            t['G' * 128 + u] = 2; t[u * 128 + 'G'] = 2;
        }
        t['G' * 128 + 'C'] = 1; t['C' * 128 + 'G'] = 1;
        return t;
    }

    #endregion

    #region Base Pairing

    /// <summary>
    /// Determines if two bases can form a canonical pair: Watson-Crick A·U / G·C or the G·U wobble
    /// (the ViennaRNA default pair set). Case-insensitive; DNA T is treated as U (A·T, G·T pair).
    /// Any other character (IUPAC ambiguity codes, N, gaps, non-ASCII) never pairs.
    /// </summary>
    public static bool CanPair(char base1, char base2)
    {
        int b1 = char.ToUpperInvariant(base1);
        int b2 = char.ToUpperInvariant(base2);
        return (b1 | b2) < 128 && PairLookup[b1 * 128 + b2] != 0;
    }

    /// <summary>
    /// Gets the type of base pair (<see cref="BasePairType.WatsonCrick"/> for A·U / G·C,
    /// <see cref="BasePairType.Wobble"/> for G·U), or null if the bases cannot pair.
    /// Same pairing rules as <see cref="CanPair"/> (case-insensitive, T treated as U).
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
    /// Enumerates hairpin stem-loops: for every hairpin loop [loopStart, loopEnd] of
    /// <paramref name="minLoopSize"/>..<paramref name="maxLoopSize"/> nucleotides closed by a canonical pair,
    /// the stem that closes it (see <paramref name="maxInternalLoopSize"/>). Overlapping candidates are all
    /// returned, ordered by loop start, then loop size. Case-insensitive; DNA T is read as U (reported bases
    /// and loop sequences use U); any other character never pairs (ViennaRNA default pair set).
    /// </summary>
    /// <remarks>
    /// <para><b>Perfect stems</b> (<paramref name="maxInternalLoopSize"/> = 0, the default): the stem is one
    /// helix — consecutively stacked canonical pairs (NNDB "helix") — grown outward from the loop-closing
    /// pair while the next outer pair can form; candidates whose helix is shorter than
    /// <paramref name="minStemLength"/> are dropped.</para>
    /// <para><b>Stems with bulges and internal loops</b> (<paramref name="maxInternalLoopSize"/> &gt; 0): the
    /// stem is the minimum-free-energy <i>unbranched</i> chain of pairs closing the hairpin loop, i.e. the
    /// Zuker–Stiegler (1981) V(i,j) recursion restricted to degree-2 loops — stacks, bulges and internal loops
    /// with at most <paramref name="maxInternalLoopSize"/> unpaired nucleotides (ViennaRNA MAXLOOP semantics,
    /// default 30 there) — scored with the same Turner 2004 loop functions as the MFE folder, over every outer
    /// closing pair; the optimum is then dropped if it has fewer than <paramref name="minStemLength"/> pairs.</para>
    /// <para><see cref="StemLoop.TotalFreeEnergy"/> is the free energy of the stem-loop as an isolated
    /// structure (its own sub-sequence <c>[Start..End]</c> and <see cref="StemLoop.DotBracketNotation"/>), the
    /// counterpart of ViennaRNA <c>eval_structure</c> on that sub-sequence: stacks, bulge/internal loops,
    /// hairpin loop and the terminal AU/GU penalties of both helix ends (no dangles: the motif has no
    /// neighbours). For a perfect stem it equals <see cref="Stem.FreeEnergy"/>
    /// (<see cref="CalculateStemEnergy"/>, NNDB helix incl. both AU/GU ends) + <see cref="CalculateHairpinLoopEnergy"/>,
    /// except for the four A-U-closed special hexaloops, whose measured loop value already contains the
    /// closure penalty (ViennaRNA <c>vrna_E_hairpin</c> returns <c>Hexaloop_E</c>). For a stem with
    /// bulges/internal loops <see cref="Stem.FreeEnergy"/> = total − hairpin-loop term.</para>
    /// </remarks>
    /// <param name="rnaSequence">RNA/DNA sequence; null or empty yields nothing.</param>
    /// <param name="minStemLength">Minimum number of base pairs in a reported stem (default 3).</param>
    /// <param name="minLoopSize">Minimum hairpin-loop size; values below 3 are raised to 3 (NNDB: hairpin loops
    /// with fewer than 3 nucleotides are prohibited).</param>
    /// <param name="maxLoopSize">Maximum hairpin-loop size (default 10).</param>
    /// <param name="allowWobble">Whether G·U pairs may be part of a stem (default true).</param>
    /// <param name="maxInternalLoopSize">0 (default) = perfect stems only; &gt; 0 = stems may contain bulges and
    /// internal loops with up to this many unpaired nucleotides each (n1 + n2).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxInternalLoopSize"/> is negative.</exception>
    public static IEnumerable<StemLoop> FindStemLoops(
        string rnaSequence,
        int minStemLength = 3,
        int minLoopSize = 3,
        int maxLoopSize = 10,
        bool allowWobble = true,
        int maxInternalLoopSize = 0)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxInternalLoopSize);
        return FindStemLoopsIterator(rnaSequence, minStemLength, minLoopSize, maxLoopSize, allowWobble, maxInternalLoopSize);
    }

    private static IEnumerable<StemLoop> FindStemLoopsIterator(
        string rnaSequence, int minStemLength, int minLoopSize, int maxLoopSize, bool allowWobble, int maxInternalLoopSize)
    {
        // NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops
        // with fewer than 3 nucleotides."
        // Wikipedia Stem-loop: "loops fewer than three bases long are sterically
        // impossible and thus do not form."
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minStemLength * 2 + minLoopSize)
            yield break;

        // Read DNA thymine as uracil (as the MFE engines do): the Turner stacking / loop tables are
        // A/C/G/U-keyed, so T-containing pairs must be scored as their U equivalents.
        string upper = rnaSequence.ToUpperInvariant().Replace('T', 'U');

        int[]? chainPredecessor = maxInternalLoopSize > 0
            ? FillUnbranchedStemDp(upper, minLoopSize, allowWobble, maxInternalLoopSize)
            : null;

        // Scan for potential hairpin loops
        for (int loopStart = minStemLength; loopStart <= upper.Length - minStemLength - minLoopSize; loopStart++)
        {
            for (int loopSize = minLoopSize; loopSize <= Math.Min(maxLoopSize, upper.Length - loopStart - minStemLength); loopSize++)
            {
                int loopEnd = loopStart + loopSize - 1;

                var stemLoop = chainPredecessor == null
                    ? TryBuildStemLoop(upper, loopStart, loopEnd, minStemLength, allowWobble)
                    : TryBuildUnbranchedStemLoop(upper, chainPredecessor, loopStart, loopEnd, minStemLength, allowWobble);
                if (stemLoop != null)
                {
                    yield return stemLoop.Value;
                }
            }
        }
    }

    private static bool IsStemPair(string seq, int i, int j, bool allowWobble)
    {
        byte t = PairType(seq[i], seq[j]);
        return t == 1 || (t == 2 && allowWobble);
    }

    private static StemLoop? TryBuildStemLoop(string sequence, int loopStart, int loopEnd, int minStemLength, bool allowWobble)
    {
        var basePairs = new List<BasePair>();

        int left = loopStart - 1;
        int right = loopEnd + 1;

        // Extend the helix outward while the next pair can form.
        while (left >= 0 && right < sequence.Length && IsStemPair(sequence, left, right, allowWobble))
        {
            basePairs.Add(new BasePair(left, right, sequence[left], sequence[right],
                GetBasePairType(sequence[left], sequence[right])!.Value));
            left--;
            right++;
        }

        if (basePairs.Count < minStemLength)
            return null;

        basePairs.Reverse(); // Order from 5' to 3'
        return BuildStemLoop(sequence, basePairs, loopStart, loopEnd, perfectHelix: true);
    }

    /// <summary>
    /// Outside DP for unbranched stems: G(k,l) = min{ AU/GU end penalty of (k,l) as the outermost pair,
    /// min over enclosing (i,j) with n1 + n2 ≤ maxInternal of G(i,j) + InteriorTerm(i,j,k,l),
    /// G(k−3,l+3) + GGUC/CUGG 3-stack }. The optimal stem closing a hairpin at (p,q) is the argmin chain
    /// of G(p,q) (the hairpin term is fixed by (p,q)). Returns the predecessor of each pair state:
    /// −1 = outermost pair, ≥ 0 = enclosing pair index i·n + j, ≤ −2 = GGUC step from index −2 − value;
    /// int.MinValue = pair not allowed.
    /// </summary>
    private static int[] FillUnbranchedStemDp(string seq, int minLoopSize, bool allowWobble, int maxInternal)
    {
        int n = seq.Length;
        var model = TurnerModel.At37;
        byte[] sc = TurnerModel.Encode(seq);
        var g = new double[n * n];
        var pred = new int[n * n];
        Array.Fill(pred, int.MinValue);
        const double tieEps = 1e-9;

        for (int span = n - 1; span > minLoopSize; span--)
        {
            for (int k = 0; k + span < n; k++)
            {
                int l = k + span;
                if (!IsStemPair(seq, k, l, allowWobble)) continue;

                double best = IsAUorGU(seq[k], seq[l]) ? TerminalAU_GU_Penalty : 0.0;
                int arg = -1;
                for (int i = k - 1; i >= 0 && k - i - 1 <= maxInternal; i--)
                {
                    int n1 = k - i - 1;
                    for (int j = l + 1; j < n && n1 + j - l - 1 <= maxInternal; j++)
                    {
                        int o = i * n + j;
                        if (pred[o] == int.MinValue) continue;
                        double e = g[o] + model.InteriorTerm(sc, i, j, k, l);
                        if (e < best - tieEps) { best = e; arg = o; }
                    }
                }

                if (k >= 3 && l + 3 < n && pred[(k - 3) * n + l + 3] != int.MinValue && allowWobble &&
                    TurnerModel.IsSpecialGgucCuggStack(sc, k - 3, l + 3))
                {
                    int o = (k - 3) * n + l + 3;
                    double e = g[o] + SpecialGGUC_CUGG_3Stack;
                    if (e < best - tieEps) { best = e; arg = -2 - o; }
                }

                g[k * n + l] = best;
                pred[k * n + l] = arg;
            }
        }
        return pred;
    }

    private static StemLoop? TryBuildUnbranchedStemLoop(
        string sequence, int[] pred, int loopStart, int loopEnd, int minStemLength, bool allowWobble)
    {
        int n = sequence.Length;
        int p = loopStart - 1, q = loopEnd + 1;
        if (p < 0 || q >= n || !IsStemPair(sequence, p, q, allowWobble) || pred[p * n + q] == int.MinValue)
            return null;

        // Trace the optimal chain outward (inner → outer).
        var chain = new List<(int I, int J)> { (p, q) };
        int cur = p * n + q;
        bool perfectHelix = true;
        while (pred[cur] != -1)
        {
            int v = pred[cur];
            int ci = cur / n, cj = cur % n;
            if (v <= -2)
            {
                chain.Add((ci - 1, cj + 1));
                chain.Add((ci - 2, cj + 2));
                chain.Add((ci - 3, cj + 3));
                cur = -2 - v;
            }
            else
            {
                if (v / n != ci - 1 || v % n != cj + 1) perfectHelix = false;
                chain.Add((v / n, v % n));
                cur = v;
            }
        }

        if (chain.Count < minStemLength)
            return null;

        chain.Reverse(); // Order from 5' to 3'
        var basePairs = chain
            .Select(c => new BasePair(c.I, c.J, sequence[c.I], sequence[c.J], GetBasePairType(sequence[c.I], sequence[c.J])!.Value))
            .ToList();
        return BuildStemLoop(sequence, basePairs, loopStart, loopEnd, perfectHelix);
    }

    // Builds the result record; basePairs are ordered outer (5') → inner (closing pair of the loop).
    private static StemLoop BuildStemLoop(string sequence, List<BasePair> basePairs, int loopStart, int loopEnd, bool perfectHelix)
    {
        int stemStart = basePairs[0].Position1;
        int stemEnd = basePairs[0].Position2;
        int stemEnd5 = loopStart - 1;
        int stemStart3 = loopEnd + 1;

        string loopSeq = sequence.Substring(loopStart, loopEnd - loopStart + 1);
        var loop = new Loop(
            Type: LoopType.Hairpin,
            Start: loopStart,
            End: loopEnd,
            Size: loopSeq.Length,
            Sequence: loopSeq);

        string dotBracket = GenerateDotBracket(sequence.Length, basePairs, stemStart, stemEnd);

        // Energy of the isolated motif (ViennaRNA eval_structure counterpart, shared loop functions).
        double totalEnergy = EvaluateStructureEnergy(sequence.Substring(stemStart, stemEnd - stemStart + 1), dotBracket);

        double stemEnergy;
        if (perfectHelix)
        {
            stemEnergy = CalculateStemEnergy(sequence, basePairs);
        }
        else
        {
            // Special GU closure in the isolated motif: G·U closing pair preceded by two Gs inside the motif.
            bool specialGU = sequence[stemEnd5] == 'G' && sequence[stemStart3] == 'U' &&
                             stemEnd5 - 2 >= stemStart && sequence[stemEnd5 - 1] == 'G' && sequence[stemEnd5 - 2] == 'G';
            stemEnergy = Math.Round(totalEnergy -
                CalculateHairpinLoopEnergy(loopSeq, sequence[stemEnd5], sequence[stemStart3], specialGU), 2);
        }

        var stem = new Stem(
            Start5Prime: stemStart,
            End5Prime: stemEnd5,
            Start3Prime: stemStart3,
            End3Prime: stemEnd,
            Length: basePairs.Count,
            BasePairs: basePairs,
            FreeEnergy: stemEnergy);

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
    /// Calculates the free energy of a helix (consecutively stacked pairs (i,j), (i+1,j−1), …, listed 5' → 3',
    /// i.e. outermost first) with the Turner 2004 nearest-neighbour model (NNDB "Watson-Crick-Franklin helix"):
    /// the sum of the stacking terms (special 5'GGUC/3'CUGG three-stack −4.12 kcal/mol, NNDB GU note b) plus
    /// the terminal AU/GU penalty (+0.45) for each helix end closed by A·U or G·U. Bases are read
    /// case-insensitively with DNA T as U (as <see cref="CanPair"/>); the AU/GU class is taken from the bases,
    /// not from <see cref="BasePair.Type"/>. A stack with an unknown pair contributes 0. The
    /// <paramref name="sequence"/> argument is not read (the bases come from the pairs).
    /// </summary>
    public static double CalculateStemEnergy(string sequence, IReadOnlyList<BasePair> basePairs)
    {
        if (basePairs.Count == 0)
            return 0;

        var b5 = new char[basePairs.Count];
        var b3 = new char[basePairs.Count];
        for (int k = 0; k < basePairs.Count; k++)
        {
            b5[k] = NormalizeBase(basePairs[k].Base1);
            b3[k] = NormalizeBase(basePairs[k].Base2);
        }

        double energy = 0;

        // Sum stacking energies (Turner 2004 nearest-neighbor).
        // S127: the GGUC/CUGG special context deliberately advances `i` inside the body to consume
        // the extra stacks it already scored; converting this well-tested energy loop is not worth the risk.
#pragma warning disable S127
        for (int i = 0; i < basePairs.Count - 1; i++)
        {
            // Special GGUC/CUGG 3-stack context (NNDB note b):
            // When 4 consecutive pairs are G-C, G-U, U-G, C-G, replace
            // 3 individual stacking terms with -4.12 total.
            if (i + 3 < basePairs.Count &&
                b5[i] == 'G' && b3[i] == 'C' &&
                b5[i + 1] == 'G' && b3[i + 1] == 'U' &&
                b5[i + 2] == 'U' && b3[i + 2] == 'G' &&
                b5[i + 3] == 'C' && b3[i + 3] == 'G')
            {
                energy += SpecialGGUC_CUGG_3Stack;
                i += 2; // Skip next 2 iterations (3 stacks consumed)
                continue;
            }

            string stack = $"{b5[i]}{b5[i + 1]}/{b3[i]}{b3[i + 1]}";
            if (StackingEnergies.TryGetValue(stack, out double stackEnergy))
                energy += stackEnergy;
            // Unknown stacking pairs contribute 0 (no data available)
        }
#pragma warning restore S127

        // Terminal AU/GU penalty: +0.45 per helix end that terminates with AU/UA or GU/UG
        // Source: NNDB — "Per AU end" (wc-parameters.html), "Per GU end" (gu-parameters.html)
        if (IsAUorGU(b5[0], b3[0]))
            energy += TerminalAU_GU_Penalty;
        if (IsAUorGU(b5[^1], b3[^1]))
            energy += TerminalAU_GU_Penalty;

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Returns the Turner 2004 terminal-mismatch stacking energy (NNDB terminal mismatches, kcal/mol) of the
    /// helix end <c>5'-closingBase5 mismatch5-3' / 3'-closingBase3 mismatch3-5'</c>: <paramref name="mismatch5"/>
    /// is the unpaired base 3' of <paramref name="closingBase5"/>, <paramref name="mismatch3"/> the unpaired base
    /// 5' of <paramref name="closingBase3"/>. This is the term of exterior and multibranch helix ends (ViennaRNA
    /// 2.7.2 <c>mismatchExt37</c> = <c>mismatchM37</c> <c>[pair(closingBase3, closingBase5)][mismatch3][mismatch5]</c>).
    /// </summary>
    /// <remarks>
    /// Bases are case-insensitive and DNA T is read as U. An unpaired base other than A/C/G/U (N, IUPAC codes,
    /// gaps) is ViennaRNA's unknown base 'N', whose <c>params/default.c</c> entries are the maximum (least
    /// stabilising) value over A/C/G/U; the same maximum is taken here. A closing pair that is not canonical
    /// (Watson–Crick or G·U, see <see cref="CanPair"/>) is outside the nearest-neighbour model and returns 0.
    /// </remarks>
    public static double GetTerminalMismatchEnergy(char closingBase5, char closingBase3, char mismatch5, char mismatch3)
    {
        return TurnerModel.At37.TerminalMismatch(closingBase5, closingBase3, mismatch5, mismatch3);
    }

    /// <summary>
    /// Returns the Turner 2004 dangling-end energy (NNDB dangling ends, kcal/mol) of a single unpaired base
    /// stacked on the helix end <c>5'-closingBase5-3' / 3'-closingBase3-5'</c>: a 3' dangle
    /// (<paramref name="is3Prime"/> = true) sits 3' of <paramref name="closingBase5"/>, a 5' dangle 5' of
    /// <paramref name="closingBase3"/> (ViennaRNA 2.7.2 <c>dangle3_37</c> / <c>dangle5_37</c>
    /// <c>[pair(closingBase3, closingBase5)][danglingBase]</c>).
    /// </summary>
    /// <remarks>
    /// Same base conventions as <see cref="GetTerminalMismatchEnergy"/>: case-insensitive, T read as U, an unknown
    /// dangling base takes the maximum over A/C/G/U (ViennaRNA's 'N' column), a non-canonical closing pair returns 0.
    /// </remarks>
    public static double GetDanglingEndEnergy(char closingBase5, char closingBase3, char danglingBase, bool is3Prime)
    {
        return TurnerModel.At37.Dangle(closingBase5, closingBase3, danglingBase, is3Prime);
    }

    private static bool IsKnownBase(char b) => b is 'A' or 'C' or 'G' or 'U';

    /// <summary>
    /// Evaluates <paramref name="term"/> for two (normalised) single-base positions; a position holding an unknown
    /// base (anything but A/C/G/U) is ViennaRNA's 'N' and is replaced by every A/C/G/U, keeping the maximum —
    /// the rule by which ViennaRNA 2.7.2 <c>params/default.c</c> fills every 'N' row/column of its Turner 2004
    /// mismatch and dangle tables.
    /// </summary>
    private static double MaxOverUnknownBases(char x, char y, Func<char, char, double> term)
    {
        if (IsKnownBase(x) && IsKnownBase(y))
            return term(x, y);

        string xs = IsKnownBase(x) ? x.ToString() : "ACGU";
        string ys = IsKnownBase(y) ? y.ToString() : "ACGU";
        double max = double.NegativeInfinity;
        foreach (char a in xs)
            foreach (char b in ys)
                max = Math.Max(max, term(a, b));
        return max;
    }

    /// <summary>
    /// Calculates the free energy of a hairpin loop using Turner 2004 parameters from NNDB.
    /// For special loops (tri/tetra/hexaloops), the total energy from experimental data
    /// replaces the model calculation. Otherwise uses:
    /// ΔG° = initiation(n) + terminal_mismatch + UU/GA_bonus + GG_bonus + special_GU_closure + all_C_penalty
    /// (n ≥ 4); ΔG° = initiation(3) + C3 penalty (n = 3, no first-mismatch term). Initiation n = 3..30 is the NNDB
    /// <c>turner_2004_loop_dg.txt</c> table (= ViennaRNA 2.7.2 <c>hairpin37</c>); n &gt; 30 is init(9) + 1.75·RT·ln(n/9).
    /// The closing pair's terminal AU/GU penalty is NOT included (NNDB counts it in the helix; see
    /// <see cref="CalculateStemEnergy"/>); a special loop's tabulated value is the loop total. Case-insensitive,
    /// T read as U; an unknown first-mismatch base takes the maximum over A/C/G/U (ViennaRNA 'N').
    /// Loops shorter than 3 nt return the prohibitive sentinel 100.0.
    /// </summary>
    /// <param name="specialGUClosure">True when the closing pair is G(5')-U(3') and the two nucleotides 5' of that G
    /// are both G (NNDB / Mathews et al. 2004); ignored for other closing pairs and for 3-nt loops.</param>
    public static double CalculateHairpinLoopEnergy(
        string loopSequence, char closingBase5, char closingBase3,
        bool specialGUClosure = false)
    {
        // Case-insensitive; DNA T is the same base as U (ViennaRNA vrna_nucleotide_encode).
        var h = new byte[loopSequence.Length + 2];
        h[0] = TurnerModel.CodeOf(closingBase5);
        for (int k = 0; k < loopSequence.Length; k++) h[k + 1] = TurnerModel.CodeOf(loopSequence[k]);
        h[^1] = TurnerModel.CodeOf(closingBase3);
        return TurnerModel.At37.HairpinLoop(h, specialGUClosure, out _);
    }

    /// <summary>
    /// Calculates the free energy of the degree-2 loop closed by the outer pair (i,j) and the inner pair
    /// (p,q) with the Turner 2004 rules, exactly as ViennaRNA 2.7.2 <c>vrna_E_internal</c>
    /// (<c>eval/eval_internal.c</c>):
    /// <list type="bullet">
    /// <item>1×1, 1×2 (2×1) and 2×2 loops: tabulated int11 / int21 / int22 values (NNDB, all closure terms included);</item>
    /// <item>1×n (n ≥ 3): initiation(n+1) + min(3.0, 0.6·|n1−n2|) + 1×n mismatch terms (0, AU/GU closure +0.7);</item>
    /// <item>2×3: initiation(5) + 0.6 + 2×3 first-mismatch terms;</item>
    /// <item>other: initiation(n1+n2) + min(3.0, 0.6·|n1−n2|) + generic first-mismatch terms;</item>
    /// <item>n1 = n2 = 0: the stacking energy; exactly one side empty: <see cref="CalculateBulgeLoopEnergy"/>.</item>
    /// </list>
    /// Each mismatch term is applied to both closing pairs read from inside the loop:
    /// <c>[pair(i,j)][s(i+1)][s(j−1)] + [pair(q,p)][s(q+1)][s(p−1)]</c>.
    /// </summary>
    /// <remarks>
    /// Sources: NNDB Turner 2004 internal loops (turner_2004_int11_dg.txt, int21_dg.txt, int22_dg.txt,
    /// internal_loop_dg.txt); ViennaRNA 2.7.2 <c>vrna_E_internal</c> and default parameter set. Bases are
    /// case-insensitive, T is read as U, and any other character is ViennaRNA's unknown base 'N'
    /// (table index 0). A non-canonical closing pair uses ViennaRNA's non-standard pair type.
    /// </remarks>
    /// <param name="n1">Number of unpaired bases on the 5' side (between i and p).</param>
    /// <param name="n2">Number of unpaired bases on the 3' side (between q and j).</param>
    /// <param name="closingBase5_1">5' base of closing pair 1 (outer, position i).</param>
    /// <param name="closingBase3_1">3' base of closing pair 1 (outer, position j).</param>
    /// <param name="closingBase5_2">5' base of closing pair 2 (inner, position p).</param>
    /// <param name="closingBase3_2">3' base of closing pair 2 (inner, position q).</param>
    /// <param name="mismatch5_1">Unpaired base adjacent to the outer pair on the 5' side (i+1).</param>
    /// <param name="mismatch3_1">Unpaired base adjacent to the outer pair on the 3' side (j−1).</param>
    /// <param name="mismatch5_2">Unpaired base adjacent to the inner pair on the 5' side (p−1).</param>
    /// <param name="mismatch3_2">Unpaired base adjacent to the inner pair on the 3' side (q+1).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="n1"/> or <paramref name="n2"/> is negative.</exception>
    public static double CalculateInternalLoopEnergy(
        int n1, int n2,
        char closingBase5_1, char closingBase3_1,
        char closingBase5_2, char closingBase3_2,
        char mismatch5_1, char mismatch3_1,
        char mismatch5_2, char mismatch3_2)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n1);
        ArgumentOutOfRangeException.ThrowIfNegative(n2);
        return TurnerModel.At37.InternalLoop(n1, n2,
            TurnerModel.CodeOf(closingBase5_1), TurnerModel.CodeOf(closingBase3_1),
            TurnerModel.CodeOf(closingBase5_2), TurnerModel.CodeOf(closingBase3_2),
            TurnerModel.CodeOf(mismatch5_1), TurnerModel.CodeOf(mismatch3_1),
            TurnerModel.CodeOf(mismatch5_2), TurnerModel.CodeOf(mismatch3_2));
    }

    // Upper-cases a nucleotide and reads DNA thymine as uracil (the shared tables are A/C/G/U keyed).
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static char NormalizeBase(char b)
    {
        char u = char.ToUpperInvariant(b);
        return u == 'T' ? 'U' : u;
    }

    private static bool IsAUorGU(char b5, char b3) => Turner2004Parameters.IsTerminalPenaltyPair(b5, b3);

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
    /// <param name="bulgeOn5PrimeStrand">
    /// Strand of the bulge: <c>true</c> = between <paramref name="pair5_base1"/> and
    /// <paramref name="pair3_base1"/>, <c>false</c> = between <paramref name="pair3_base2"/> and
    /// <paramref name="pair5_base2"/>. NNDB's special C bulge (−0.9) needs a paired C ADJACENT to the bulged C,
    /// i.e. one of the two closing bases on the bulged strand. When <c>null</c> (strand not given) the rule cannot
    /// be decided from the arguments; the bonus is then applied only when it holds on BOTH strands, i.e. the
    /// result is the larger (least stabilising) of the two readings — the convention ViennaRNA 2.7.2 uses for
    /// every unknown in its Turner 2004 tables (maximum over the possibilities). Pass the strand for the exact
    /// NNDB value.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="bulgeSize"/> is less than 1.</exception>
    /// <remarks>Bases are case-insensitive and DNA T is read as U.</remarks>
    public static double CalculateBulgeLoopEnergy(
        int bulgeSize, char bulgedBase,
        char pair5_base1, char pair5_base2,
        char pair3_base1, char pair3_base2,
        int numStates = 1, bool? bulgeOn5PrimeStrand = null)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(bulgeSize, 1);
        return TurnerModel.At37.Bulge(bulgeSize, TurnerModel.CodeOf(bulgedBase),
            TurnerModel.CodeOf(pair5_base1), TurnerModel.CodeOf(pair5_base2),
            TurnerModel.CodeOf(pair3_base1), TurnerModel.CodeOf(pair3_base2),
            numStates, bulgeOn5PrimeStrand);
    }

    /// <summary>
    /// Initiation free energy of a multibranch loop in the linear Turner 2004 model used for dynamic programming —
    /// the model of every folding engine in this class (<see cref="CalculateMinimumFreeEnergy"/>,
    /// <see cref="CalculateMfeStructure"/>, <see cref="CalculatePartitionFunction"/>):
    /// ΔG°init = a + b·[unpaired nucleotides] + c·[branching helices, the closing helix included],
    /// a = 9.3, b = 0, c = −0.9 kcal/mol (Mathews et al. 2004 PNAS 101:7287; ViennaRNA 2.7.2
    /// <c>params/default.c</c> <c>ML_closing37 = 930</c>, <c>ML_BASE37 = 0</c>, <c>ML_intern37 = −90</c>),
    /// plus <paramref name="stackingEnergy"/> (the dangles / terminal mismatches of the branches, see
    /// <see cref="GetTerminalMismatchEnergy"/>). The NNDB efn2 model with the asymmetry and strain terms is the
    /// overload taking the unpaired-nucleotide counts between the helices.
    /// </summary>
    /// <param name="numHelices">Number of helices in the loop, the closing helix included (≥ 3).</param>
    /// <param name="numUnpaired">Total number of unpaired nucleotides in the loop (≥ 0; b = 0 in this model).</param>
    /// <param name="hasStrain">
    /// Must be <c>false</c>: the three-way strain term exists only in the NNDB efn2 model; use
    /// <see cref="CalculateMultibranchLoopEnergy(IReadOnlyList{int}, double)"/>, which derives it from the loop.
    /// </param>
    /// <param name="stackingEnergy">Stacking (dangles / terminal mismatches / coaxial stacks) of the loop, added as given.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="numHelices"/> &lt; 3 or <paramref name="numUnpaired"/> &lt; 0.</exception>
    /// <exception cref="ArgumentException"><paramref name="hasStrain"/> is <c>true</c>.</exception>
    public static double CalculateMultibranchLoopEnergy(
        int numHelices, int numUnpaired, bool hasStrain = false, double stackingEnergy = 0.0)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(numHelices, 3);
        ArgumentOutOfRangeException.ThrowIfNegative(numUnpaired);
        if (hasStrain)
            throw new ArgumentException(
                "The strain term belongs to the NNDB efn2 multibranch model; call the overload with the unpaired " +
                "counts between the helices, which applies it to three-way loops with fewer than two unpaired nucleotides.",
                nameof(hasStrain));

        double energy = Dp_ML_closing + Dp_ML_base * numUnpaired + Dp_ML_intern * numHelices + stackingEnergy;
        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Free energy of a multibranch loop in the NNDB Turner 2004 (efn2) model (NNDB "Multibranch Loops" §9.1;
    /// Mathews &amp; Turner 2002, Biochemistry 41:869):
    /// ΔG° = a + b·[average asymmetry] + c·[branching helices] + ΔG°strain + ΔG°stacking with a = 9.25, b = 0.91,
    /// c = −0.63, strain = +3.14 kcal/mol; average asymmetry = min[2.0, mean over the helices of
    /// |unpaired nucleotides on the 5' side − unpaired nucleotides on the 3' side|]; the strain applies to
    /// three-way loops with fewer than two unpaired nucleotides.
    /// </summary>
    /// <param name="unpairedBetweenHelices">
    /// The loop read 5'→3' starting after the closing pair: element k is the number of unpaired nucleotides
    /// between helix k and helix k+1 (cyclically; helix 0 is the closing helix, the last element is the gap
    /// before the 3' base of the closing pair). Its length is the number of helices (≥ 3); helix k has
    /// element k−1 on its 5' side and element k on its 3' side.
    /// </param>
    /// <param name="stackingEnergy">
    /// Optimal configuration of dangling ends, terminal mismatches and coaxial stacks (NNDB §9.1), added as given
    /// (see <see cref="GetDanglingEndEnergy"/>, <see cref="GetTerminalMismatchEnergy"/>,
    /// <see cref="CalculateFlushCoaxialStacking"/>, <see cref="CalculateMismatchCoaxialStacking"/>).
    /// </param>
    /// <exception cref="ArgumentNullException"><paramref name="unpairedBetweenHelices"/> is null.</exception>
    /// <exception cref="ArgumentException">Fewer than 3 helices, or a negative count.</exception>
    public static double CalculateMultibranchLoopEnergy(IReadOnlyList<int> unpairedBetweenHelices, double stackingEnergy = 0.0)
    {
        ArgumentNullException.ThrowIfNull(unpairedBetweenHelices);
        int h = unpairedBetweenHelices.Count;
        if (h < 3)
            throw new ArgumentException("A multibranch loop has at least 3 helices (the closing helix included).",
                nameof(unpairedBetweenHelices));

        int total = 0;
        double asymmetrySum = 0;
        for (int k = 0; k < h; k++)
        {
            int after = unpairedBetweenHelices[k];
            if (after < 0)
                throw new ArgumentException("Unpaired counts must be non-negative.", nameof(unpairedBetweenHelices));
            int before = unpairedBetweenHelices[(k + h - 1) % h];
            asymmetrySum += Math.Abs(before - after);
            total += after;
        }

        double averageAsymmetry = Math.Min(MultibranchLoop_MaxAsymmetry, asymmetrySum / h);
        double energy = MultibranchLoop_Offset
                      + MultibranchLoop_Asymmetry * averageAsymmetry
                      + MultibranchLoop_Helix * h
                      + stackingEnergy;
        if (h == 3 && total < 2)
            energy += MultibranchLoop_Strain;

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Flush coaxial stacking energy of two helices with no intervening unpaired nucleotide (NNDB Turner 2004
    /// "Coaxial Stacking" §8.2): the Watson–Crick / G·U helical stack as though there were no backbone break.
    /// The continuous strand runs 5'→3' through <paramref name="base5_1"/> then <paramref name="base5_2"/>;
    /// (<paramref name="base5_1"/>, <paramref name="base3_1"/>) and (<paramref name="base5_2"/>,
    /// <paramref name="base3_2"/>) are the two terminal pairs — e.g. an A·U pair followed by a C·G pair is
    /// <c>('A','U','C','G')</c> = −2.24 (NNDB example). Case-insensitive, T read as U; a non-canonical pair
    /// returns 0 (outside the model). The terminal AU/GU penalties of the helices are not included (NNDB note).
    /// </summary>
    public static double CalculateFlushCoaxialStacking(
        char base5_1, char base3_1, char base5_2, char base3_2)
    {
        string stack = $"{NormalizeBase(base5_1)}{NormalizeBase(base5_2)}/{NormalizeBase(base3_1)}{NormalizeBase(base3_2)}";
        return StackingEnergies.TryGetValue(stack, out double e) ? Math.Round(e, 2) : 0.0;
    }

    /// <summary>
    /// Mismatch-mediated coaxial stacking energy (NNDB Turner 2004 "Coaxial Stacking" §8.3 and parameter table):
    /// ΔG° = terminal mismatch of the helix end with the continuous backbone
    /// (<see cref="GetTerminalMismatchEnergy"/>) − 2.1 (stack on the second helix, backbone not continuous)
    /// − 0.4 when the "mismatch" (<paramref name="mismatch5"/>, <paramref name="mismatch3"/>) could form a
    /// Watson–Crick pair, or − 0.2 when it could form a G·U pair; no bonus otherwise. NNDB example: A·U followed
    /// by a GG mismatch = −0.8 − 2.1 = −2.9 kcal/mol. Case-insensitive, T read as U; an unknown mismatch base takes
    /// the maximum over A/C/G/U (ViennaRNA 'N' convention); a non-canonical closing pair returns 0.
    /// </summary>
    public static double CalculateMismatchCoaxialStacking(
        char closingBase5, char closingBase3, char mismatch5, char mismatch3)
    {
        char c5 = NormalizeBase(closingBase5);
        char c3 = NormalizeBase(closingBase3);
        if (!CanPair(c5, c3))
            return 0.0;

        double energy = MaxOverUnknownBases(NormalizeBase(mismatch5), NormalizeBase(mismatch3), (x, y) =>
            TerminalMismatchEnergies[$"{c5}{x}{y}{c3}"] + CoaxialMismatch_Base + GetBasePairType(x, y) switch
            {
                BasePairType.WatsonCrick => CoaxialMismatch_WC_Bonus,
                BasePairType.Wobble => CoaxialMismatch_GU_Bonus,
                _ => 0.0,
            });
        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the minimum free energy (MFE) of an RNA sequence using Zuker-style DP
    /// with Turner 2004 nearest-neighbor parameters (NNDB), the loop decomposition and d2 dangle
    /// model of ViennaRNA 2.7.2 (<c>RNAfold</c> default). Returns physical kcal/mol values.
    ///
    /// Recurrences (see <see cref="FillDp"/>):
    ///   W(j)    = min{ W(j-1), min_i W(i-1) + V(i,j) + Ext(i,j) }            (Ext: AU/GU + terminal mismatch / dangle)
    ///   V(i,j)  = min{ Hairpin(i,j), min_(p,q) Interior(i,j,p,q) + V(p,q), Mclose(i,j) + min_k WM(i+1,k)+WM(k+1,j-1) }
    ///   WM(i,j) = min{ V(i,j)+Mstem(i,j), WM(i+1,j), WM(i,j-1), min_k WM(i,k)+WM(k+1,j) }
    /// Interior loops use the Turner 2004 int11/int21/int22 tables and the 1×n / 2×3 / generic
    /// mismatch rules (<see cref="CalculateInternalLoopEnergy"/>, ViennaRNA <c>vrna_E_internal</c>).
    ///
    /// Complexity: O(n³) with MAXLOOP=30 for internal/bulge loops.
    /// </summary>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <param name="temperature">
    /// Absolute temperature in Kelvin (default 310.15 K = 37 °C). Every Turner 2004 parameter is rescaled with its
    /// enthalpy, ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15 K (ViennaRNA <c>md.temperature</c>); see <see cref="TurnerModel"/>.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="temperature"/> is not a positive finite number.</exception>
    public static double CalculateMinimumFreeEnergy(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        var model = ModelFor(temperature);

        // NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops
        // with fewer than 3 nucleotides."
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minLoopSize + 2)
            return 0;

        // Accept DNA input by reading thymine as uracil: T and U are the same base for folding
        // (Watson–Crick A–U/A–T and stacking parameters are identical), so the energy model and the
        // A/C/G/U-indexed tables treat the two interchangeably.
        byte[] s = TurnerModel.Encode(rnaSequence);
        int n = s.Length;

        var pool = ArrayPool<double>.Shared;
        double[] vBuf = pool.Rent(n * n);  // V(i,j)
        double[] wmBuf = pool.Rent(n * n); // WM(i,j)
        double[] w = pool.Rent(n);          // W(j)

        try
        {
            FillDp(model, s, n, minLoopSize, vBuf, wmBuf, w);
            return RoundEnergy(model, w[n - 1]);
        }
        finally
        {
            pool.Return(vBuf);
            pool.Return(wmBuf);
            pool.Return(w);
        }
    }

    // Validates a Kelvin temperature and returns the matching parameter set (37 °C = NNDB tables).
    private static TurnerModel ModelFor(double temperature)
    {
        if (!(temperature > 0) || double.IsInfinity(temperature))
            throw new ArgumentOutOfRangeException(nameof(temperature), "Temperature must be a positive finite number (Kelvin).");
        return TurnerModel.ForTemperature(temperature);
    }

    // Energies are reported to the 0.01 kcal/mol of the NNDB tables at 37 °C; rescaled parameters are not tabulated.
    private static double RoundEnergy(TurnerModel model, double e) => model.IsReference ? Math.Round(e, 2) : e;

    // Default folding temperature 37 °C = 310.15 K (ViennaRNA / NNDB default; the temperature the parameters were measured at).
    private const double DefaultTemperatureKelvin = TurnerModel.ReferenceTemperatureKelvin;
    // Minimum number of unpaired bases enclosed by a hairpin loop (NNDB Turner 2004: hairpin loops < 3 nt are prohibited).
    private const int MinHairpinLoop = 3;

    // DP constants shared by the matrix fill, the scalar MFE, the traceback and the partition function.
    private const double Dp_INF = TurnerModel.Inf;
    private const int Dp_MAXLOOP = 30;

    // Multibranch loop, linear DP model ΔG = a + b·[unpaired] + c·[branching helices incl. the closing one]
    // (the NNDB efn2 asymmetry/strain form of CalculateMultibranchLoopEnergy cannot be decomposed by an
    // O(n³) recursion). Source: ViennaRNA 2.7.2 params/default.c ML_closing37 = 930, ML_intern37 = −90,
    // ML_BASE37 = 0 (Turner 2004 set, as used by vrna_E_multibranch_stem / E_MLstem).
    private const double Dp_ML_closing = TurnerModel.MlClosing37;   // a
    private const double Dp_ML_intern = TurnerModel.MlIntern37;   // c (per branch, the closing pair included)
    private const double Dp_ML_base = TurnerModel.MlBase37;      // b (per unpaired nucleotide)

    /// <summary>
    /// Fills the Zuker–Stiegler (1981) MFE dynamic-programming matrices in place, with the Turner 2004
    /// loop decomposition and the d2 dangle model of ViennaRNA 2.7.2 (<c>RNAfold -d2</c>, the default):
    /// <code>
    /// V(i,j)  = min{ H(i,j), min_(p,q) I(i,j,p,q) + V(p,q), GGUC3(i,j) + V(i+3,j−3),
    ///                Mclose(i,j) + min_k WM(i+1,k) + WM(k+1,j−1) }
    /// WM(i,j) = min{ V(i,j) + Mstem(i,j), WM(i+1,j) + b, WM(i,j−1) + b, min_k WM(i,k) + WM(k+1,j) }
    /// W(j)    = min{ W(j−1), min_i W(i−1) + V(i,j) + Ext(i,j) }
    /// </code>
    /// The multiloop split demands ≥ 2 inner branches, so a multiloop is never scored with a single
    /// branch (that would be an interior loop). Caller supplies n×n buffers for <paramref name="vBuf"/>
    /// and <paramref name="wmBuf"/> (row-major, index i*n+j) and an n-length buffer for <paramref name="w"/>.
    /// The traceback (<see cref="TracebackMfe"/>) and the scalar MFE reuse the SAME fill.
    /// </summary>
    private static void FillDp(TurnerModel m, byte[] s, int n, int minLoopSize, double[] vBuf, double[] wmBuf, double[] w,
        KnotLayer? knots = null)
    {
        const double INF = Dp_INF;
        const int MAXLOOP = Dp_MAXLOOP;
        double mlBase = m.MlBase;

        Array.Fill(vBuf, INF, 0, n * n);
        Array.Fill(wmBuf, INF, 0, n * n);
        Array.Fill(w, 0.0, 0, n);

        for (int span = minLoopSize + 2; span <= n; span++)
        {
            for (int i = 0; i <= n - span; i++)
            {
                int j = i + span - 1;
                int ij = i * n + j;

                // ===== V(i,j): i and j pair =====
                if (TurnerModel.CanPair(s[i], s[j]))
                {
                    double vBest = m.HairpinTerm(s, i, j);

                    // Stack, bulge and interior loops with inner pair (ip,jp), n1+n2 ≤ MAXLOOP.
                    for (int ip = i + 1; ip < j - 1 && ip - i - 1 <= MAXLOOP; ip++)
                    {
                        int jpMin = Math.Max(ip + 1, j - (MAXLOOP - (ip - i - 1)) - 1);
                        for (int jp = j - 1; jp >= jpMin; jp--)
                        {
                            if (!TurnerModel.CanPair(s[ip], s[jp])) continue;
                            double vInner = vBuf[ip * n + jp];
                            if (vInner >= INF) continue;
                            double total = m.InteriorTerm(s, i, j, ip, jp) + vInner;
                            if (total < vBest) vBest = total;
                        }
                    }

                    // Special GGUC/CUGG 3-stack (NNDB GU pairs, note b).
                    if (TurnerModel.IsSpecialGgucCuggStack(s, i, j))
                    {
                        double vInner3 = vBuf[(i + 3) * n + (j - 3)];
                        if (vInner3 < INF && SpecialGGUC_CUGG_3Stack + vInner3 < vBest)
                            vBest = SpecialGGUC_CUGG_3Stack + vInner3;
                    }

                    // Multiloop closed by (i,j): at least two branches inside.
                    double split = INF;
                    for (int k = i + 1; k < j - 1; k++)
                    {
                        double left = wmBuf[(i + 1) * n + k];
                        if (left >= INF) continue;
                        double right = wmBuf[(k + 1) * n + (j - 1)];
                        if (right >= INF) continue;
                        if (left + right < split) split = left + right;
                    }
                    if (split < INF)
                    {
                        double mEnergy = m.MultiloopClosing(s, i, j) + split;
                        if (mEnergy < vBest) vBest = mEnergy;
                    }

                    vBuf[ij] = vBest;
                }

                // Pseudoknot layer (pknotsRG grammar): K(i,j) and the knot-loop segment C(i,j) of this cell.
                knots?.FillCell(i, j, vBuf);

                // ===== WM(i,j): part of a multiloop with ≥ 1 branch (only meaningful for 0 < i, j < n−1) =====
                if (i > 0 && j < n - 1)
                {
                    double wmBest = INF;
                    if (vBuf[ij] < INF)
                        wmBest = vBuf[ij] + m.MultiloopStem(s, i, j);
                    if (knots is { NestedKnots: true } && knots.K[ij] < INF && knots.K[ij] + Pseudoknot_ComponentPenalty < wmBest)
                        wmBest = knots.K[ij] + Pseudoknot_ComponentPenalty; // pknotsRG pkml: a knot as a multiloop branch

                    double e = wmBuf[(i + 1) * n + j];
                    if (e < INF && e + mlBase < wmBest) wmBest = e + mlBase;
                    e = wmBuf[i * n + (j - 1)];
                    if (e < INF && e + mlBase < wmBest) wmBest = e + mlBase;

                    for (int k = i; k < j; k++)
                    {
                        double left = wmBuf[i * n + k];
                        if (left >= INF) continue;
                        double right = wmBuf[(k + 1) * n + j];
                        if (right >= INF) continue;
                        if (left + right < wmBest) wmBest = left + right;
                    }

                    wmBuf[ij] = wmBest;
                }
            }
        }

        // ===== W(j): exterior loop of the prefix [0..j] =====
        for (int j = 0; j < n; j++)
        {
            double wBest = j > 0 ? w[j - 1] : 0.0;
            for (int i = 0; i + minLoopSize + 1 <= j; i++)
            {
                double vij = vBuf[i * n + j];
                if (vij >= INF) continue;
                double total = (i > 0 ? w[i - 1] : 0.0) + vij + m.ExteriorStem(s, i, j);
                if (total < wBest) wBest = total;
            }
            if (knots is { NestedKnots: true })
                wBest = Math.Min(wBest, knots.BestExteriorKnot(j, w));
            w[j] = wBest;
            knots?.FillExterior(j, vBuf, w);
        }
    }

    /// <summary>
    /// Computes the MFE-OPTIMAL secondary structure of an RNA sequence by Zuker–Stiegler (1981)
    /// dynamic-programming traceback over the same V/W/WM matrices used by
    /// <see cref="CalculateMinimumFreeEnergy"/>. Unlike <see cref="PredictStructure"/> (an assembly of
    /// isolated hairpin stem-loops), this returns the globally optimal pseudoknot-free structure for the
    /// Turner 2004 energy model, and its <see cref="MfeStructure.FreeEnergy"/> equals the scalar
    /// MFE for the same input and parameters.
    /// </summary>
    /// <remarks>
    /// Recurrences and traceback per Zuker M, Stiegler P (1981) Nucleic Acids Res. 9(1):133-148
    /// (W/V matrices; F/C/M/M¹ decomposition as taught in MIT 6.047 Lecture 08, Fig. 13). The
    /// traceback re-derives, at each cell, which recurrence option attained the stored minimum and
    /// recurses into the corresponding sub-problem(s), recording a base pair whenever a V(i,j) cell
    /// is entered — exactly mirroring <see cref="FillDp"/>, so the two are mutually consistent.
    /// </remarks>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    /// <returns>The optimal structure; empty (no pairs, all dots, ΔG = 0) for null/empty/too-short input.</returns>
    public static MfeStructure CalculateMfeStructure(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        var model = ModelFor(temperature);

        // NNDB Turner 2004: hairpin loops < 3 nt are prohibited.
        if (minLoopSize < 3) minLoopSize = 3;

        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minLoopSize + 2)
        {
            string s0 = rnaSequence is null ? "" : rnaSequence.ToUpperInvariant().Replace('T', 'U');
            return new MfeStructure(s0, new string('.', s0.Length),
                new List<(int, int)>(), 0.0);
        }

        // Accept DNA input by reading thymine as uracil (see CalculateMinimumFreeEnergy).
        string seq = rnaSequence.ToUpperInvariant().Replace('T', 'U');
        byte[] s = TurnerModel.Encode(seq);
        int n = seq.Length;

        // Non-pooled arrays: they must survive past the fill so the traceback can read them.
        double[] vBuf = new double[n * n];
        double[] wmBuf = new double[n * n];
        double[] w = new double[n];

        FillDp(model, s, n, minLoopSize, vBuf, wmBuf, w);

        var pairs = TracebackMfe(model, s, n, minLoopSize, vBuf, wmBuf, w);
        pairs.Sort((a, b) => a.Position1.CompareTo(b.Position1));

        string dotBracket = GenerateFullDotBracket(n, pairs);
        double energy = RoundEnergy(model, w[n - 1]);

        return new MfeStructure(seq, dotBracket, pairs, energy);
    }

    /// <summary>
    /// Predicts the MFE-optimal secondary structure via DP traceback and returns it in the same
    /// <see cref="SecondaryStructure"/> shape as <see cref="PredictStructure"/>. This is the
    /// global-optimum counterpart of the stem-loop assembly <see cref="PredictStructure"/>: the dot-bracket and base
    /// pairs are the DP optimum (Zuker–Stiegler 1981) and <see cref="SecondaryStructure.MinimumFreeEnergy"/>
    /// equals <see cref="CalculateMinimumFreeEnergy"/>.
    /// </summary>
    /// <param name="rnaSequence">RNA (or DNA; T read as U) sequence; case-insensitive.</param>
    /// <param name="minLoopSize">Minimum hairpin loop size (NNDB minimum 3; smaller is clamped to 3).</param>
    public static SecondaryStructure PredictStructureMfe(
        string rnaSequence, int minLoopSize = 3, double temperature = DefaultTemperatureKelvin)
    {
        var mfe = CalculateMfeStructure(rnaSequence, minLoopSize, temperature);

        var basePairs = new List<BasePair>(mfe.BasePairs.Count);
        foreach (var (p1, p2) in mfe.BasePairs)
        {
            var type = GetBasePairType(mfe.Sequence[p1], mfe.Sequence[p2]) ?? BasePairType.NonCanonical;
            basePairs.Add(new BasePair(p1, p2, mfe.Sequence[p1], mfe.Sequence[p2], type));
        }

        var pseudoknots = DetectPseudoknots(basePairs).ToList(); // optimal DP is pseudoknot-free → empty

        return new SecondaryStructure(
            Sequence: mfe.Sequence,
            DotBracket: mfe.DotBracket,
            BasePairs: basePairs,
            StemLoops: new List<StemLoop>(),
            Pseudoknots: pseudoknots,
            MinimumFreeEnergy: mfe.FreeEnergy);
    }

    // Comparison tolerance for "this option achieved the stored DP optimum". The DP sums a bounded
    // number of table values (each ≤ 2 decimals), so the optimum is reproduced to well within this.
    private const double Dp_TraceEps = 1e-6;

    /// <summary>
    /// Reconstructs the optimal base-pair set from the filled V/W/WM matrices by Zuker–Stiegler
    /// traceback. At each sub-problem it recomputes the candidate values of the SAME recurrence
    /// options as <see cref="FillDp"/> and follows the one that attains the stored optimum.
    /// </summary>
    private static List<(int Position1, int Position2)> TracebackMfe(
        TurnerModel m, byte[] s, int n, int minLoopSize, double[] vBuf, double[] wmBuf, double[] w,
        KnotLayer? knots = null)
    {
        const double INF = Dp_INF;
        const int MAXLOOP = Dp_MAXLOOP;
        const double EPS = Dp_TraceEps;

        var pairs = new List<(int, int)>();

        // Work items. 'W' → prefix [0..j]; 'V' → pair (i,j), reconstruct what it encloses;
        // 'M' → multiloop region [i..j] with ≥ 1 branch.
        var stack = new Stack<(char Kind, int I, int J)>();
        stack.Push((knots is { NestedKnots: false } ? 'X' : 'W', 0, n - 1)); // 'X': prefix with at most one knot

        while (stack.Count > 0)
        {
            var (kind, i, j) = stack.Pop();

            if (kind is 'X' or 'K' or 'C')
            {
                knots!.Trace(kind, i, j, stack, pairs, vBuf, w);
                continue;
            }

            if (kind == 'W')
            {
                if (j < 0) continue;
                double target = w[j];
                if (j > 0 && Math.Abs(target - w[j - 1]) <= EPS)
                {
                    stack.Push(('W', 0, j - 1));
                    continue;
                }
                bool matched = false;
                for (int k = 0; k + minLoopSize + 1 <= j; k++)
                {
                    double vkj = vBuf[k * n + j];
                    if (vkj >= INF) continue;
                    double left = k > 0 ? w[k - 1] : 0.0;
                    if (Math.Abs(target - (left + vkj + m.ExteriorStem(s, k, j))) <= EPS)
                    {
                        if (k > 0) stack.Push(('W', 0, k - 1));
                        stack.Push(('V', k, j));
                        matched = true;
                        break;
                    }
                }
                if (!matched && knots is { NestedKnots: true })
                {
                    for (int k = 0; k + 10 <= j; k++)
                    {
                        double kk = knots.K[k * n + j];
                        if (kk >= INF) continue;
                        if (Math.Abs(target - ((k > 0 ? w[k - 1] : 0.0) + kk)) <= EPS)
                        {
                            if (k > 0) stack.Push(('W', 0, k - 1));
                            stack.Push(('K', k, j));
                            break;
                        }
                    }
                }
                continue; // nothing matched only for an open chain (target 0)
            }

            if (kind == 'V')
            {
                pairs.Add((i, j));
                double target = vBuf[i * n + j];

                if (Math.Abs(target - m.HairpinTerm(s, i, j)) <= EPS) continue;

                bool found = false;
                for (int ip = i + 1; ip < j - 1 && ip - i - 1 <= MAXLOOP && !found; ip++)
                {
                    int jpMin = Math.Max(ip + 1, j - (MAXLOOP - (ip - i - 1)) - 1);
                    for (int jp = j - 1; jp >= jpMin; jp--)
                    {
                        if (!TurnerModel.CanPair(s[ip], s[jp])) continue;
                        double vInner = vBuf[ip * n + jp];
                        if (vInner >= INF) continue;
                        if (Math.Abs(target - (m.InteriorTerm(s, i, j, ip, jp) + vInner)) <= EPS)
                        {
                            stack.Push(('V', ip, jp));
                            found = true;
                            break;
                        }
                    }
                }
                if (found) continue;

                if (TurnerModel.IsSpecialGgucCuggStack(s, i, j))
                {
                    double vInner3 = vBuf[(i + 3) * n + (j - 3)];
                    if (vInner3 < INF && Math.Abs(target - (SpecialGGUC_CUGG_3Stack + vInner3)) <= EPS)
                    {
                        // The 3-stack records the two intermediate pairs too (i+1,j-1),(i+2,j-2).
                        pairs.Add((i + 1, j - 1));
                        pairs.Add((i + 2, j - 2));
                        stack.Push(('V', i + 3, j - 3));
                        continue;
                    }
                }

                double close = m.MultiloopClosing(s, i, j);
                for (int k = i + 1; k < j - 1; k++)
                {
                    double left = wmBuf[(i + 1) * n + k];
                    double right = wmBuf[(k + 1) * n + (j - 1)];
                    if (left >= INF || right >= INF) continue;
                    if (Math.Abs(target - (close + left + right)) <= EPS)
                    {
                        stack.Push(('M', i + 1, k));
                        stack.Push(('M', k + 1, j - 1));
                        break;
                    }
                }
                continue;
            }

            // kind == 'M'
            {
                double target = wmBuf[i * n + j];
                if (target >= INF) continue;

                double vij = vBuf[i * n + j];
                if (vij < INF && Math.Abs(target - (vij + m.MultiloopStem(s, i, j))) <= EPS)
                {
                    stack.Push(('V', i, j));
                    continue;
                }
                if (knots is { NestedKnots: true } && knots.K[i * n + j] < INF
                    && Math.Abs(target - (knots.K[i * n + j] + Pseudoknot_ComponentPenalty)) <= EPS)
                {
                    stack.Push(('K', i, j));
                    continue;
                }
                double e = wmBuf[(i + 1) * n + j];
                if (e < INF && Math.Abs(target - (e + m.MlBase)) <= EPS)
                {
                    stack.Push(('M', i + 1, j));
                    continue;
                }
                e = wmBuf[i * n + (j - 1)];
                if (e < INF && Math.Abs(target - (e + m.MlBase)) <= EPS)
                {
                    stack.Push(('M', i, j - 1));
                    continue;
                }
                for (int k = i; k < j; k++)
                {
                    double left = wmBuf[i * n + k];
                    double right = wmBuf[(k + 1) * n + j];
                    if (left < INF && right < INF && Math.Abs(target - (left + right)) <= EPS)
                    {
                        stack.Push(('M', i, k));
                        stack.Push(('M', k + 1, j));
                        break;
                    }
                }
            }
        }

        return pairs;
    }

    /// <summary>
    /// Free energy ΔG°37 (kcal/mol) of a given pseudoknot-free secondary structure — the Turner 2004 loop-by-loop
    /// sum (exterior d2 stem terms, hairpins, stacks, bulges, interior loops, GGUC/CUGG three-stack, linear d2
    /// multiloops) of exactly the model minimised by <see cref="CalculateMinimumFreeEnergy"/>; the counterpart of
    /// ViennaRNA 2.7.2 <c>vrna_eval_structure</c> (<c>RNAeval</c>, <c>-d2</c>). For every sequence,
    /// <c>CalculateFreeEnergy(s, CalculateMfeStructure(s).DotBracket)</c> equals the MFE.
    /// </summary>
    /// <param name="rnaSequence">RNA (or DNA; T read as U), case-insensitive.</param>
    /// <param name="dotBracket">Structure of the same length using only '(', ')' and '.'.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">
    /// The structure has another character, is unbalanced, has the wrong length, pairs bases that cannot pair
    /// (see <see cref="CanPair"/>), or contains a hairpin loop of fewer than 3 unpaired nucleotides (prohibited by
    /// the Turner 2004 rules; ViennaRNA evaluates such a structure to INF).
    /// </exception>
    public static double CalculateFreeEnergy(
        string rnaSequence, string dotBracket, double temperature = DefaultTemperatureKelvin)
    {
        ArgumentNullException.ThrowIfNull(rnaSequence);
        ArgumentNullException.ThrowIfNull(dotBracket);
        foreach (char c in dotBracket)
            if (c is not ('(' or ')' or '.'))
                throw new ArgumentException("Only '(', ')' and '.' are allowed (pseudoknot-free structure).", nameof(dotBracket));
        return EvaluateStructureEnergy(rnaSequence, dotBracket, ModelFor(temperature));
    }

    internal static double EvaluateStructureEnergy(string rnaSequence, string dotBracket, TurnerModel? model = null)
    {
        ArgumentNullException.ThrowIfNull(rnaSequence);
        ArgumentNullException.ThrowIfNull(dotBracket);
        model ??= TurnerModel.At37;
        byte[] s = TurnerModel.Encode(rnaSequence);
        int n = s.Length;
        if (dotBracket.Length != n || !ValidateDotBracket(dotBracket))
            throw new ArgumentException("Dot-bracket must be balanced and as long as the sequence.", nameof(dotBracket));

        var partner = new int[n];
        Array.Fill(partner, -1);
        foreach (var (a, b) in ParseDotBracket(dotBracket))
        {
            if (!TurnerModel.CanPair(s[a], s[b]))
                throw new ArgumentException($"Non-canonical pair ({a},{b}).", nameof(dotBracket));
            // NNDB Turner 2004: "The nearest neighbor rules prohibit hairpin loops with fewer than 3 nucleotides"
            // (ViennaRNA 2.7.2 vrna_E_hairpin: hairpin[0..2] = INF; eval_structure reports 100000). A pair spanning
            // fewer than 3 positions cannot enclose another pair, so it always closes such a hairpin.
            if (b - a - 1 < MinHairpinLoop)
                throw new ArgumentException(
                    $"Pair ({a},{b}) closes a hairpin loop of {b - a - 1} nucleotides; Turner 2004 prohibits hairpin loops of fewer than 3.",
                    nameof(dotBracket));
            partner[a] = b;
            partner[b] = a;
        }

        double energy = 0.0;
        int k = 0;
        while (k < n)
        {
            if (partner[k] > k)
            {
                energy += model.ExteriorStem(s, k, partner[k]) + LoopEnergyClosedBy(model, s, partner, k, partner[k]);
                k = partner[k];
            }
            k++;
        }
        return RoundEnergy(model, energy);
    }

    // Energy of the loop closed by (i,j) plus everything nested inside it.
    private static double LoopEnergyClosedBy(TurnerModel m, byte[] s, int[] partner, int i, int j)
    {
        if (TurnerModel.IsSpecialGgucCuggStack(s, i, j) &&
            partner[i + 1] == j - 1 && partner[i + 2] == j - 2 && partner[i + 3] == j - 3)
            return SpecialGGUC_CUGG_3Stack + LoopEnergyClosedBy(m, s, partner, i + 3, j - 3);

        var branches = new List<(int P, int Q)>();
        int k = i + 1;
        while (k < j)
        {
            if (partner[k] > k)
            {
                branches.Add((k, partner[k]));
                k = partner[k];
            }
            k++;
        }

        if (branches.Count == 0)
            return m.HairpinTerm(s, i, j);

        if (branches.Count == 1)
        {
            var (p, q) = branches[0];
            return m.InteriorTerm(s, i, j, p, q) + LoopEnergyClosedBy(m, s, partner, p, q);
        }

        double e = m.MultiloopClosing(s, i, j);
        foreach (var (p, q) in branches)
            e += m.MultiloopStem(s, p, q) + LoopEnergyClosedBy(m, s, partner, p, q);
        return e;
    }

    // Hot-path pair classification for already-normalised (uppercase) sequences: reads the same
    // canonical PairLookup table as CanPair/GetBasePairType (0 = none, 1 = WC, 2 = wobble).
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static byte PairType(char b1, char b2)
        => (b1 | b2) < 128 ? PairLookup[b1 * 128 + b2] : (byte)0;

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
    /// Assembles a secondary structure from hairpin stem-loops: the minimum-free-energy set of mutually
    /// non-overlapping candidates of <see cref="FindStemLoops"/> (stems of ≥ <paramref name="minStemLength"/> stacked
    /// pairs closing a hairpin of <paramref name="minLoopSize"/>..<paramref name="maxLoopSize"/> nt) placed in the
    /// exterior loop. Case-insensitive; DNA T is read as U (the returned sequence is the normalised A/C/G/U string).
    /// </summary>
    /// <remarks>
    /// The selection is the Zuker–Stiegler (1981) exterior-loop recursion restricted to the candidate set:
    /// <c>W(j) = min{ W(j−1), min over candidates (i,j) of W(i−1) + E(i,j) }</c>, where <c>E(i,j)</c> is the free energy of
    /// the whole sequence folded into that single stem-loop (<see cref="CalculateFreeEnergy"/>: stem, hairpin and the
    /// d2 exterior terminal-mismatch/dangle term, which under d2 depends on the sequence only, so the terms of disjoint
    /// stem-loops add up). It is therefore the exact optimum over every combination of candidates — no candidate whose
    /// energy is ≥ 0 is ever chosen and <see cref="SecondaryStructure.MinimumFreeEnergy"/> is the free energy of the
    /// returned structure (= <see cref="CalculateFreeEnergy"/> on its dot-bracket, ≤ 0). The search space has no
    /// multiloops, interior loops, bulges or nested hairpins, so the result is in general LESS stable than the global
    /// optimum: for the minimum-free-energy fold use <see cref="PredictStructureMfe"/> / <see cref="CalculateMfeStructure"/>
    /// (ViennaRNA RNAfold counterpart).
    /// </remarks>
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

        string seq = rnaSequence.ToUpperInvariant().Replace('T', 'U');

        var candidates = FindStemLoops(seq, minStemLength, minLoopSize, maxLoopSize).ToList();
        var selectedStemLoops = SelectStemLoopsByExteriorDp(seq, candidates);

        var allBasePairs = selectedStemLoops
            .SelectMany(sl => sl.Stem.BasePairs)
            .OrderBy(bp => bp.Position1)
            .ToList();

        string dotBracket = GenerateFullDotBracket(seq.Length, allBasePairs);

        // Free energy of the assembled structure with the canonical evaluator (RNAeval -d2 counterpart).
        double mfe = EvaluateStructureEnergy(seq, dotBracket);

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

    /// <summary>
    /// Zuker–Stiegler exterior-loop recursion over a fixed candidate set of stem-loops:
    /// W(j) = min{ W(j−1), min_(c ends at j) W(c.Start−1) + E(c) }, E(c) = free energy of the full sequence folded
    /// into stem-loop c alone. Ties keep the unpaired option (W(j−1)), then the first candidate in input order.
    /// </summary>
    private static List<StemLoop> SelectStemLoopsByExteriorDp(string seq, List<StemLoop> candidates)
    {
        int n = seq.Length;
        var byEnd = new List<(int Index, double Energy)>?[n];
        var db = new char[n];
        for (int c = 0; c < candidates.Count; c++)
        {
            var sl = candidates[c];
            Array.Fill(db, '.');
            foreach (var bp in sl.Stem.BasePairs)
            {
                db[bp.Position1] = '(';
                db[bp.Position2] = ')';
            }
            double e = EvaluateStructureEnergy(seq, new string(db));
            (byEnd[sl.End] ??= new List<(int, double)>()).Add((c, e));
        }

        const double eps = 1e-9;
        var w = new double[n + 1];      // w[j + 1] = W(j); w[0] = 0 (empty prefix)
        var choice = new int[n + 1];    // candidate index ending at j, or −1 = j unpaired
        choice[0] = -1;
        for (int j = 0; j < n; j++)
        {
            w[j + 1] = w[j];
            choice[j + 1] = -1;
            if (byEnd[j] is null) continue;
            foreach (var (index, e) in byEnd[j]!)
            {
                double v = w[candidates[index].Start] + e;
                if (v < w[j + 1] - eps)
                {
                    w[j + 1] = v;
                    choice[j + 1] = index;
                }
            }
        }

        var selected = new List<StemLoop>();
        for (int j = n; j > 0;)
        {
            if (choice[j] < 0) { j--; continue; }
            var sl = candidates[choice[j]];
            selected.Add(sl);
            j = sl.Start;
        }
        selected.Reverse();
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
    /// This is the pairwise conflict relation (a helix of a pairs crossing a helix of b pairs gives a·b
    /// results); the layer decomposition is <see cref="AssignPseudoknotOrders"/> and the bracket rendering
    /// <see cref="ToDotBracket"/>. Positions are not validated and the scan never throws: pairs that share an
    /// endpoint (not a valid secondary structure) never satisfy the strict inequalities. O(p²) for p pairs,
    /// which is also the worst-case output size.
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

    // Window variant used by stem-loop detection: renders only [start, end] (inclusive) of the
    // structure, positions re-based to the window start (pairs lying outside the window are dropped).
    private static string GenerateDotBracket(int length, IReadOnlyList<BasePair> basePairs, int start, int end)
        => ToDotBracket(end - start + 1,
            basePairs.Where(bp => Math.Min(bp.Position1, bp.Position2) >= start && Math.Max(bp.Position1, bp.Position2) <= end)
                     .Select(bp => (bp.Position1 - start, bp.Position2 - start)));

    private static string GenerateFullDotBracket(int length, IReadOnlyList<BasePair> basePairs)
        => ToDotBracket(length, basePairs.Select(bp => (bp.Position1, bp.Position2)));

    // Tuple overload used by the MFE traceback and the pseudoknot predictors.
    private static string GenerateFullDotBracket(int length, IReadOnlyList<(int Position1, int Position2)> basePairs)
        => ToDotBracket(length, basePairs);

    // Bracket families in the order ViennaRNA 2.7.2 vrna_db_from_ptable (structures/structure_dotbracket.c)
    // hands them out: "([{<" then the ASCII letters A..Z / a..z — 30 families, the same set ParseDotBracket reads.
    // Antczak et al. (2018) DBL notation uses the same sequence for pseudoknot orders 0, 1, 2, ….
    private const string BracketFamilyOpen = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private const string BracketFamilyClose = ")]}>abcdefghijklmnopqrstuvwxyz";

    /// <summary>
    /// Renders a base-pair set (pseudoknots allowed) as dot-bracket(-letter) notation, assigning bracket
    /// families exactly as ViennaRNA 2.7.2 <c>vrna_db_from_ptable</c>: repeated 5'→3' passes; in each pass a
    /// pair receives the current family when it nests inside the pairs already opened in that pass, otherwise
    /// it is deferred to the next pass; families are used in the order <c>() [] {} &lt;&gt; Aa … Zz</c>.
    /// The result is the inverse of <see cref="ParseDotBracket"/>: parsing it returns the same pair set.
    /// </summary>
    /// <param name="length">Length of the structure (number of positions).</param>
    /// <param name="basePairs">Base pairs as 0-based position tuples, in any order and either orientation.</param>
    /// <returns>A string of <paramref name="length"/> characters.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="length"/> is negative.</exception>
    /// <exception cref="ArgumentNullException"><paramref name="basePairs"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// A pair has equal ends or an end outside <c>[0, length)</c>; a position occurs in more than one pair (not a
    /// secondary structure: ViennaRNA's pair table cannot hold it); or more than 30 families would be needed
    /// (ViennaRNA logs a warning and silently drops the remaining pairs; here the input is rejected so the
    /// rendering is never lossy).
    /// </exception>
    public static string ToDotBracket(int length, IEnumerable<(int Position1, int Position2)> basePairs)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentNullException.ThrowIfNull(basePairs);

        int[] partner = BuildPartnerTable(length, basePairs, nameof(basePairs));
        var notation = new char[length];
        Array.Fill(notation, '.');

        // vrna_db_from_ptable: one stack of closing positions per pass; a pair that would cross the
        // innermost open pair of this pass is deferred ("recheck") to the next bracket family.
        var stack = new Stack<int>();
        int family = 0;
        bool recheck = true;
        while (recheck)
        {
            recheck = false;
            for (int i = 0; i < length; i++)
            {
                int p = partner[i];
                if (p > i)
                {
                    if (stack.Count > 0 && p > stack.Peek())
                    {
                        recheck = true;
                        continue;
                    }
                    if (family >= BracketFamilyOpen.Length)
                        throw new ArgumentException(
                            $"The structure needs more than {BracketFamilyOpen.Length} bracket families (ViennaRNA's limit).",
                            nameof(basePairs));
                    stack.Push(p);
                    notation[i] = BracketFamilyOpen[family];
                    notation[p] = BracketFamilyClose[family];
                }
                else if (p >= 0 && stack.Count > 0 && i == stack.Peek())
                {
                    partner[i] = partner[p] = -1; // rendered: remove from the pair table
                    stack.Pop();
                }
            }
            stack.Clear();
            family++;
        }

        return new string(notation);
    }

    /// <summary>
    /// Assigns a pseudoknot order to every base pair: order 0 = a maximum-cardinality nested (pseudoknot-free)
    /// subset of the pairs; the removed pairs are processed again to get order 1, and so on until every pair has
    /// an order. This is the procedure of biotite <c>structure.pseudoknots</c> (Smit et al. 2008 optimal
    /// pseudoknot removal applied iteratively to the removed pairs), with the pseudoknot order notion of
    /// Antczak et al. (2018): pairs of equal order never cross.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each layer is an exact maximum nested subset (Nussinov-type interval DP restricted to the given pairs,
    /// O(m²) time and memory for m = 2 × remaining pairs; ViennaRNA <c>vrna_pt_pk_remove</c> keeps a subset of the
    /// same size). When several maximum subsets exist (biotite returns all of them) the one that keeps the pair
    /// opening furthest 5' is chosen at the first difference (lexicographically greatest by opening position),
    /// so an H-type knot <c>([)]</c> of two equal stems keeps the 5' stem at order 0.
    /// </para>
    /// <para>
    /// The number of layers is an upper bound on Antczak's minimum decomposition count: that minimum is the
    /// chromatic number of the crossing (circle) graph, whose computation is NP-complete (Garey, Johnson, Miller
    /// &amp; Papadimitriou 1980), and maximising the order-0 layer can cost one extra layer
    /// (e.g. pairs (0,4),(1,12),(2,6),(3,10),(8,13),(9,14): 4 layers here, 3 colours suffice).
    /// </para>
    /// </remarks>
    /// <param name="basePairs">Base pairs (either orientation, any order).</param>
    /// <returns>The order of each input pair, index-aligned with <paramref name="basePairs"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="basePairs"/> is null.</exception>
    /// <exception cref="ArgumentException">A pair has equal or negative ends, or a position occurs in more than one pair.</exception>
    public static IReadOnlyList<int> AssignPseudoknotOrders(IReadOnlyList<BasePair> basePairs)
    {
        ArgumentNullException.ThrowIfNull(basePairs);
        var tuples = new (int, int)[basePairs.Count];
        int maxPos = -1;
        for (int t = 0; t < basePairs.Count; t++)
        {
            tuples[t] = (basePairs[t].Position1, basePairs[t].Position2);
            maxPos = Math.Max(maxPos, Math.Max(basePairs[t].Position1, basePairs[t].Position2));
        }
        BuildPartnerTable(maxPos + 1, tuples, nameof(basePairs)); // validation only

        var orders = new int[basePairs.Count];
        var remaining = Enumerable.Range(0, basePairs.Count).ToList();
        int order = 0;
        while (remaining.Count > 0)
        {
            var keep = MaximumNestedSubset(remaining.Select(t => tuples[t]).ToList());
            var next = new List<int>();
            for (int r = 0; r < remaining.Count; r++)
            {
                if (keep[r]) orders[remaining[r]] = order;
                else next.Add(remaining[r]);
            }
            remaining = next;
            order++;
        }
        return orders;
    }

    /// <summary>
    /// Renders a base-pair set in Antczak et al. (2018) dot-bracket-letter notation: each pair uses the bracket
    /// family of its <see cref="AssignPseudoknotOrders"/> order (0 <c>()</c>, 1 <c>[]</c>, 2 <c>{}</c>,
    /// 3 <c>&lt;&gt;</c>, 4… <c>Aa</c>, <c>Bb</c>, …).
    /// </summary>
    /// <exception cref="ArgumentException">As <see cref="ToDotBracket"/>, or more than 30 orders are needed.</exception>
    public static string ToDotBracketByPseudoknotOrder(int length, IReadOnlyList<BasePair> basePairs)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(length);
        ArgumentNullException.ThrowIfNull(basePairs);
        BuildPartnerTable(length, basePairs.Select(bp => (bp.Position1, bp.Position2)), nameof(basePairs));
        var orders = AssignPseudoknotOrders(basePairs);
        var notation = new char[length];
        Array.Fill(notation, '.');
        for (int t = 0; t < basePairs.Count; t++)
        {
            if (orders[t] >= BracketFamilyOpen.Length)
                throw new ArgumentException(
                    $"The structure needs more than {BracketFamilyOpen.Length} pseudoknot orders.", nameof(basePairs));
            int i = Math.Min(basePairs[t].Position1, basePairs[t].Position2);
            int j = Math.Max(basePairs[t].Position1, basePairs[t].Position2);
            notation[i] = BracketFamilyOpen[orders[t]];
            notation[j] = BracketFamilyClose[orders[t]];
        }
        return new string(notation);
    }

    // Pair table of a secondary structure (partner[i] = j, −1 = unpaired). Rejects pairs that cannot be a
    // secondary structure: equal ends, ends outside [0, length), or a position used by two pairs.
    private static int[] BuildPartnerTable(int length, IEnumerable<(int Position1, int Position2)> basePairs, string paramName)
    {
        var partner = new int[length];
        Array.Fill(partner, -1);
        foreach (var (p1, p2) in basePairs)
        {
            if (p1 == p2)
                throw new ArgumentException($"Base pair ({p1},{p2}) pairs a position with itself.", paramName);
            if (p1 < 0 || p2 < 0 || p1 >= length || p2 >= length)
                throw new ArgumentException($"Base pair ({p1},{p2}) lies outside [0, {length}).", paramName);
            if (partner[p1] != -1 || partner[p2] != -1)
                throw new ArgumentException(
                    $"Base pair ({p1},{p2}) reuses a position that is already paired; a secondary structure pairs each base at most once.",
                    paramName);
            partner[p1] = p2;
            partner[p2] = p1;
        }
        return partner;
    }

    // Maximum-cardinality nested subset of a valid pair set (no shared positions). Interval DP over the sorted
    // endpoints (Smit et al. 2008 / Nussinov recursion restricted to the given pairs):
    //   N(a,b) = max( N(a+1,b), 1 + N(a+1,c−1) + N(c+1,b) )  where c is the partner of a, a < c ≤ b.
    // Traceback prefers keeping the pair at a on ties (lexicographically greatest subset by opening position).
    // Returns keep flags index-aligned with pairs.
    private static bool[] MaximumNestedSubset(IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        int p = pairs.Count;
        var keep = new bool[p];
        if (p == 0) return keep;

        int m = 2 * p;
        var ends = new (int Pos, int Pair)[m];
        for (int t = 0; t < p; t++)
        {
            ends[2 * t] = (pairs[t].Position1, t);
            ends[2 * t + 1] = (pairs[t].Position2, t);
        }
        Array.Sort(ends, (x, y) => x.Pos.CompareTo(y.Pos));
        var mate = new int[m];  // endpoint index of the partner
        var first = new int[p];
        Array.Fill(first, -1);
        for (int e = 0; e < m; e++)
        {
            int t = ends[e].Pair;
            if (first[t] < 0) first[t] = e;
            else { mate[e] = first[t]; mate[first[t]] = e; }
        }

        // n[a, b] for 0 ≤ a ≤ b < m, row-major; empty intervals (b < a) read as 0.
        var n = new int[m * m];
        int N(int a, int b) => b < a ? 0 : n[a * m + b];
        for (int a = m - 1; a >= 0; a--)
        {
            int c = mate[a];
            for (int b = a; b < m; b++)
            {
                int best = N(a + 1, b);
                if (c > a && c <= b)
                    best = Math.Max(best, 1 + N(a + 1, c - 1) + N(c + 1, b));
                n[a * m + b] = best;
            }
        }

        var work = new Stack<(int A, int B)>();
        work.Push((0, m - 1));
        while (work.Count > 0)
        {
            var (a, b) = work.Pop();
            while (a <= b)
            {
                int c = mate[a];
                if (c > a && c <= b && 1 + N(a + 1, c - 1) + N(c + 1, b) == N(a, b))
                {
                    keep[ends[a].Pair] = true;
                    work.Push((c + 1, b));
                    b = c - 1;
                }
                a++;
            }
        }
        return keep;
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

    // Pseudoknot letter pairs are the ASCII alphabet only: uppercase 'A'..'Z' opens (5'),
    // the matching lowercase 'a'..'z' closes (3'). ViennaRNA 2.7.2 vrna_ptable_from_string
    // (structures/structure_pairtable.c, VRNA_BRACKETS_ALPHA) loops i = 65..90 pairing
    // (char)i with (char)(i + 32); every other character — including non-ASCII letters such
    // as 'É'/'é' — is unpaired.
    private static bool IsLetterOpener(char c) => c is >= 'A' and <= 'Z';
    private static bool IsLetterCloser(char c) => c is >= 'a' and <= 'z';

    /// <summary>
    /// Single canonical dot-bracket scanner shared by <see cref="ParseDotBracket"/> and
    /// <see cref="ValidateDotBracket"/>. Matches every family on its own LIFO stack (ViennaRNA
    /// <c>extract_pairs</c> per bracket type), appends matched pairs to <paramref name="pairs"/>
    /// (when non-null) in closing order, and returns true iff the string is well-formed
    /// (no closer without a same-family opener, no opener left unclosed). Scanning is
    /// best-effort: unmatched closers are skipped and scanning continues.
    /// </summary>
    private static bool ScanDotBracket(string? dotBracket, List<(int Position1, int Position2)>? pairs)
    {
        if (string.IsNullOrEmpty(dotBracket))
            return true;

        // One independent stack per opening symbol (brackets + uppercase letters).
        var stacks = new Dictionary<char, Stack<int>>();
        bool wellFormed = true;

        for (int i = 0; i < dotBracket.Length; i++)
        {
            char c = dotBracket[i];
            char opener;

            if (OpeningToClosing.ContainsKey(c) || IsLetterOpener(c))
            {
                if (!stacks.TryGetValue(c, out var open))
                {
                    open = new Stack<int>();
                    stacks[c] = open;
                }
                open.Push(i);
                continue;
            }

            if (ClosingToOpening.TryGetValue(c, out char bracketOpener))
                opener = bracketOpener;
            else if (IsLetterCloser(c))
                opener = (char)(c - ('a' - 'A'));
            else
                continue; // '.', ',', '-', ':', '_', '~' and any other symbol => unpaired.

            if (stacks.TryGetValue(opener, out var stack) && stack.Count > 0)
            {
                int j = stack.Pop(); // pop unconditionally (not inside a null-conditional call)
                pairs?.Add((j, i));
            }
            else
                wellFormed = false; // closer without a same-family opener
        }

        foreach (var stack in stacks.Values)
            if (stack.Count != 0)
                wellFormed = false; // opener left unclosed

        return wellFormed;
    }

    /// <summary>
    /// Parses dot-bracket / extended WUSS notation to extract base pairs as
    /// (5' position, 3' position) index tuples (zero-based).
    /// </summary>
    /// <remarks>
    /// Each bracket family (<c>()</c>, <c>[]</c>, <c>{}</c>, <c>&lt;&gt;</c>) is matched on
    /// its own stack, and uppercase/lowercase ASCII letter pairs (<c>Aa</c> … <c>Zz</c>) are
    /// matched as independent pseudoknot systems, per ViennaRNA / WUSS conventions
    /// (<c>vrna_ptable_from_string</c> with <c>VRNA_BRACKETS_ANY</c>).
    /// An opening symbol is paired with the nearest unmatched opening of the SAME family
    /// (correct nesting). Dots and any other symbols (including non-ASCII letters) denote
    /// unpaired residues and are skipped. Pairs are returned in the order their closing
    /// symbol is read. Closing symbols with no open partner are ignored (parse is best-effort,
    /// whereas ViennaRNA returns no table; use <see cref="ValidateDotBracket"/> to test
    /// well-formedness first).
    /// </remarks>
    public static IEnumerable<(int Position1, int Position2)> ParseDotBracket(string dotBracket)
    {
        var pairs = new List<(int Position1, int Position2)>();
        ScanDotBracket(dotBracket, pairs);
        return pairs;
    }

    /// <summary>
    /// Validates dot-bracket / extended WUSS notation: every closing symbol must match an
    /// earlier unmatched opening symbol of the SAME family, and no opening symbol may be left
    /// unclosed. Each bracket family and each ASCII letter (case) pair is checked independently,
    /// so a mismatched pairing such as <c>(]</c> is rejected even though the total count is balanced.
    /// Equivalent to ViennaRNA <c>vrna_ptable_from_string(s, VRNA_BRACKETS_ANY) != NULL</c> for
    /// non-empty input; null/empty is accepted as the empty structure.
    /// </summary>
    /// <remarks>
    /// Sources: ViennaRNA — base pairs are matching pairs of <c>()</c> (and extended <c>[]{}&lt;&gt;</c>,
    /// letters A–Z/a–z) and the structure must be balanced; WUSS — left and right partners must match up.
    /// </remarks>
    public static bool ValidateDotBracket(string dotBracket) => ScanDotBracket(dotBracket, null);

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
    /// Finds inverted repeats (potential RNA stem regions): a left arm <c>W</c>, a loop <c>G</c> and a right arm that
    /// pairs antiparallel with <c>W</c> (<c>W G W̄ᴿ</c>). The model, the maximality rule and the search engine are the
    /// canonical ones of <see cref="RepeatFinder.FindInvertedRepeats(string, int, int, int, int, int, bool)"/>
    /// (EMBOSS <c>palindrome</c> stems); only the pairing rule is RNA-specific.
    /// </summary>
    /// <param name="sequence">Nucleotide sequence (DNA or RNA; case-insensitive; T is read as U).</param>
    /// <param name="minLength">Minimum arm length in base pairs (default 4; values &lt; 1 yield no results).</param>
    /// <param name="minSpacing">Minimum loop length (default 3; negative yields no results).</param>
    /// <param name="maxSpacing">Maximum loop length (default 100; &lt; <paramref name="minSpacing"/> yields no results).</param>
    /// <param name="maxMismatches">
    /// Maximum number of non-pairing positions inside a stem (default 0 = perfect stems; must be ≥ 0). EMBOSS
    /// <c>palindrome -nummismatches</c> / IUPACpal "k mismatches" semantics: mismatches are interior only (never at an end
    /// of the stem). Gapped (bulge) alignments are not reported here: use
    /// <see cref="RepeatFinder.FindInvertedRepeatsScored(string, int, int, int, int, int)"/> (EMBOSS <c>einverted</c>).
    /// </param>
    /// <param name="allowWobble">
    /// When true (default) a pair is any pair of <see cref="CanPair(char, char)"/> — A·U, G·C and G·U wobble (ViennaRNA
    /// <c>pair[][]</c>); when false only Watson–Crick A·U (A·T) and G·C pairs are stems. Ambiguity codes (N, R, Y, …) and
    /// any other character never pair in either mode.
    /// </param>
    /// <returns>
    /// The maximal stems (a stem is dropped only when it lies inside another stem in both arms), as tuples: left arm
    /// <c>[Start1..End1]</c>, right arm <c>[Start2..End2]</c> (0-based inclusive) and the common arm <c>Length</c>;
    /// ordered by <c>Start1</c>, then <c>Start2</c>. Different stems may share positions (alternative pairings are
    /// not suppressed). With <c>maxMismatches = 0</c> every pair is <c>CanPair(seq[Start1+k], seq[End2-k])</c>.
    /// </returns>
    /// <remarks>
    /// Sources: inverted repeat = <c>WGW̄ᴿ</c> (Alamro et al. 2021, IUPACpal, doi:10.1186/s12859-021-03983-2; Ussery et al.
    /// 2008); maximal stems and nesting filter = EMBOSS <c>palindrome</c> (Rice et al. 2000); RNA pairing =
    /// ViennaRNA <c>pair[][]</c> as in <see cref="CanPair(char, char)"/>. No separate implementation: the search is
    /// delegated to <c>RepeatFinder</c> (DUP_MAP RNA entry; the former unsourced greedy "longest first, non-overlapping"
    /// selection and the strict-IUPAC self-pairing of N/R/Y… are removed).
    /// </remarks>
    public static IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)> FindInvertedRepeats(
        string sequence,
        int minLength = DefaultMinArmLength,
        int minSpacing = DefaultMinLoopLength,
        int maxSpacing = DefaultMaxLoopLength,
        int maxMismatches = 0,
        bool allowWobble = true)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        if (string.IsNullOrEmpty(sequence) || minLength < 1 || minSpacing < 0 || maxSpacing < minSpacing)
            return Array.Empty<(int, int, int, int, int)>();

        // T is read as U (CanPair); the shared engine pairs ACGT, so U is folded onto T (positions unchanged).
        string seq = sequence.ToUpperInvariant().Replace('U', 'T');
        return RepeatFinder.FindMaximalStemsUpper(seq, minLength, maxSpacing, minSpacing, maxMismatches, allowWobble)
            .Select(r => (r.LeftArmStart, r.LeftArmStart + r.ArmLength - 1,
                          r.RightArmStart, r.RightArmStart + r.ArmLength - 1, r.ArmLength))
            .ToList();
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
        ArgumentNullException.ThrowIfNull(random);
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
