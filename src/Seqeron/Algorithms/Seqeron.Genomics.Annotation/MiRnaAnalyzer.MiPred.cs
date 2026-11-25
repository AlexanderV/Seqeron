using System.Buffers.Binary;
using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Published pre-miRNA classifier: MiPred (Jiang P, Wu H, Wang W, Ma W, Sun X, Lu Z. 2007, "MiPred: classification of real and
/// pseudo microRNA precursors using random forest prediction model with combined features", Nucleic Acids Research
/// 35(Web Server issue):W339-W344, doi:10.1093/nar/gkm368) with the authors' own trained random forest
/// (see <c>MiRnaAnalyzer.MiPredModel.cs</c>) and the feature definitions of the package script
/// <c>microRNAcheck_parallel.pl</c> ported line by line. OPT-IN: the in-repo logistic model behind
/// <see cref="ClassifyPreMiRna"/> and <see cref="ClassifyPreMiRnaTripletSvm"/> are unchanged.
/// </summary>
public static partial class MiRnaAnalyzer
{
    /// <summary>Number of MiPred features: 32 triplet-element frequencies, the MFE and the randomisation p-value.</summary>
    public const int MiPredFeatureCount = 34;

    /// <summary>Number of dinucleotide shuffles of the package script (<c>$shuffle_times = 999</c>).</summary>
    public const int MiPredDefaultShuffleCount = 999;

    /// <summary>Seed of the shuffle RNG used when the caller supplies none (the p-value feature is statistical; fixed seed ⇒ reproducible).</summary>
    public const int MiPredDefaultSeed = 1;

    /// <summary>Package pre-judgement bounds (<c>microRNAcheck_parallel.pl</c>): MFE ≤ −15 kcal/mol, ≥ 18 pairs, no multiloop, length 51..137.</summary>
    private const double MiPredMaxFreeEnergy = -15.0;
    private const int MiPredMinBasePairs = 18;
    private const int MiPredMinLength = 51;
    private const int MiPredMaxLength = 137;

    // Package regex m/\(\.+\)[\.,\(,\)]*\(\.+\)/ : two terminal loops, i.e. a second hairpin (the package's character class also lists ',' which never occurs in a structure).
    private static readonly Regex MiPredSecondHairpin = new(@"\(\.+\)[.,()]*\(\.+\)", RegexOptions.CultureInvariant);

    // Package list2 = ("(((","((.","(.(","(..",".((",".(.","..(","...") and list1 = (A,U,G,C): feature index = nucleotide * 8 + triplet.
    private static readonly string[] MiPredTriplets = { "(((", "((.", "(.(", "(..", ".((", ".(.", "..(", "..." };

    /// <summary>Outcome of the package's hairpin pre-judgement (<see cref="AssessMiPredHairpin"/>).</summary>
    /// <param name="BasePairs">Number of '(' of the structure, or 0 when the structure contains a second hairpin (package leaves the count undefined).</param>
    /// <param name="HasMultiloop">True iff a '(' follows the first ')' (package multiloop test).</param>
    /// <param name="IsHairpinLike">MFE ≤ −15, ≥ 18 pairs, no multiloop and 51 ≤ length ≤ 137.</param>
    public readonly record struct MiPredHairpinCheck(int BasePairs, bool HasMultiloop, bool IsHairpinLike);

    /// <summary>
    /// Result of <see cref="ClassifyPreMiRnaMiPred"/>.
    /// </summary>
    /// <param name="Sequence">Sequence after the package's clean-up (upper case, T→U, characters other than A/C/G/T/U dropped).</param>
    /// <param name="DotBracket">MFE structure of the canonical fold.</param>
    /// <param name="FreeEnergy">MFE ΔG° (kcal/mol).</param>
    /// <param name="Check">Hairpin pre-judgement; when false the package does not classify (not a pre-miRNA-like hairpin).</param>
    /// <param name="PValue">Shuffle p-value (1 + #{shuffle ΔG° ≤ ΔG°})/(N + 1); null when not hairpin-like (not computed, as in the package).</param>
    /// <param name="Features">The 34 features V1..V34 (null when not hairpin-like).</param>
    /// <param name="RealVotes">Trees voting class P (real pre-miRNA); 0 when not classified.</param>
    /// <param name="TreeCount">Number of trees of the forest (1000).</param>
    /// <param name="RealProbability">RealVotes / TreeCount (the package's <c>predict(type="prob")[,2]</c>); null when not classified.</param>
    /// <param name="IsPreMiRna">True iff hairpin-like and <paramref name="RealProbability"/> &gt; 0.5 (package: "It is a real microRNA precursor").</param>
    public readonly record struct MiPredResult(
        string Sequence,
        string DotBracket,
        double FreeEnergy,
        MiPredHairpinCheck Check,
        double? PValue,
        double[]? Features,
        int RealVotes,
        int TreeCount,
        double? RealProbability,
        bool IsPreMiRna);

    /// <summary>
    /// Classifies a candidate hairpin with the published MiPred model (Jiang et al. 2007): the sequence is
    /// folded (canonical MFE, <see cref="RnaSecondaryStructure.CalculateMfeStructure"/>; the package called
    /// RNAfold 1.x), pre-judged as hairpin-like (<see cref="AssessMiPredHairpin"/>) and, if so, coded as the 34
    /// features of <c>microRNAcheck_parallel.pl</c> (<see cref="ComputeMiPredFeatures"/>: 32 triplet-element
    /// frequencies, MFE, p-value of <paramref name="shuffleCount"/> dinucleotide shuffles) and voted on by the
    /// authors' 1000-tree random forest (class P probability &gt; 0.5 ⇒ real pre-miRNA).
    /// <para>
    /// Differences from the package, all documented: the MFE and structure come from this library's Turner-2004
    /// engine (package: RNAfold 1.x, Turner 1999 parameters); the shuffles use <see cref="DinucleotideShuffle"/>
    /// (uShuffle-exact Eulerian shuffle; the package used randfold's triplet-swap shuffle with Perl's
    /// <c>rand</c>), so the p-value feature is statistical and differs between RNGs by its sampling error
    /// (±0.01 at p ≈ 0.1 with 999 shuffles); the default RNG is seeded with <see cref="MiPredDefaultSeed"/>.
    /// </para>
    /// </summary>
    /// <param name="sequence">Candidate (RNA or DNA, case-insensitive; characters other than A/C/G/T/U are dropped like the package's FASTA clean-up).</param>
    /// <param name="random">RNG for the shuffles (default: <c>new Random(MiPredDefaultSeed)</c>).</param>
    /// <param name="shuffleCount">Number of shuffles N (package: 999; the forest was trained with that value).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The result; <c>null</c> when nothing remains of the sequence after the clean-up.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="shuffleCount"/> is &lt; 1 or <paramref name="minLoopSize"/> is negative.</exception>
    public static MiPredResult? ClassifyPreMiRnaMiPred(
        string sequence,
        Random? random = null,
        int shuffleCount = MiPredDefaultShuffleCount,
        int minLoopSize = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(shuffleCount, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopSize);

        string seq = NormaliseMiPredQuery(sequence);
        if (seq.Length == 0)
            return null;

        var mfe = RnaSecondaryStructure.CalculateMfeStructure(seq, minLoopSize);
        var check = AssessMiPredHairpin(mfe.DotBracket, mfe.FreeEnergy);
        int trees = MiPredForestData.Value.TreeCount;
        if (!check.IsHairpinLike)
            return new MiPredResult(seq, mfe.DotBracket, mfe.FreeEnergy, check, null, null, 0, trees, null, false);

        var randfold = CalculateRandfoldPValue(seq, shuffleCount, random ?? new Random(MiPredDefaultSeed), minLoopSize)!.Value;
        double[] features = ComputeMiPredFeatures(seq, mfe.DotBracket, mfe.FreeEnergy, randfold.PValue);
        var (votes, total) = MiPredVotes(features);
        double p = (double)votes / total;
        return new MiPredResult(seq, mfe.DotBracket, mfe.FreeEnergy, check, randfold.PValue, features, votes, total, p, p > 0.5);
    }

    /// <summary>
    /// The package's hairpin pre-judgement (<c>microRNAcheck_parallel.pl</c>, "Check"): hairpin-like iff
    /// MFE ≤ −15 kcal/mol, the structure has ≥ 18 base pairs (count of '(', only evaluated when the structure does
    /// not contain a second terminal loop), no multiloop (no '(' after the first ')') and 51 ≤ length ≤ 137.
    /// </summary>
    /// <param name="dotBracket">MFE structure in <c>.()</c> notation (its length is the sequence length).</param>
    /// <param name="freeEnergy">MFE (kcal/mol).</param>
    public static MiPredHairpinCheck AssessMiPredHairpin(string dotBracket, double freeEnergy)
    {
        ArgumentNullException.ThrowIfNull(dotBracket);
        int pairs = MiPredSecondHairpin.IsMatch(dotBracket) ? 0 : CountChar(dotBracket, '(');
        int firstClose = dotBracket.IndexOf(')');
        bool multiloop = firstClose >= 0 && dotBracket.IndexOf('(', firstClose) >= 0;
        bool like = freeEnergy <= MiPredMaxFreeEnergy && pairs >= MiPredMinBasePairs && !multiloop
                    && dotBracket.Length <= MiPredMaxLength && dotBracket.Length >= MiPredMinLength;
        return new MiPredHairpinCheck(pairs, multiloop, like);
    }

    /// <summary>
    /// The 34 MiPred features V1..V34 (package order): V1..V32 = frequencies of the 32 triplet elements (centre
    /// nucleotide A,U,G,C × local structure "(((", "((.", "(.(", "(..", ".((", ".(.", "..(", "..."; closing brackets are not
    /// counted, positions 2..n−1 only) divided by the sequence LENGTH n; V33 = MFE (kcal/mol); V34 = the shuffle p-value.
    /// </summary>
    /// <param name="sequence">Sequence (RNA or DNA, case-insensitive).</param>
    /// <param name="dotBracket">Structure of the same length (<c>.()</c>).</param>
    /// <param name="freeEnergy">MFE (kcal/mol) → V33.</param>
    /// <param name="pValue">Shuffle p-value → V34.</param>
    public static double[] ComputeMiPredFeatures(string sequence, string dotBracket, double freeEnergy, double pValue)
    {
        double[] f = ComputeMiPredTripletFrequencies(sequence, dotBracket);
        var all = new double[MiPredFeatureCount];
        Array.Copy(f, all, f.Length);
        all[32] = freeEnergy;
        all[33] = pValue;
        return all;
    }

    /// <summary>The 32 triplet-element frequencies V1..V32 of <see cref="ComputeMiPredFeatures"/> (count / sequence length).</summary>
    /// <exception cref="ArgumentException">Lengths differ, the structure has characters other than <c>.()</c> or the sequence has characters other than A/C/G/T/U.</exception>
    public static double[] ComputeMiPredTripletFrequencies(string sequence, string dotBracket)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(dotBracket);
        string seq = NormaliseMiPredQuery(sequence);
        if (seq.Length != sequence.Length)
            throw new ArgumentException("sequence may contain only A, C, G, T/U.", nameof(sequence));
        if (dotBracket.Length != seq.Length)
            throw new ArgumentException("dotBracket must have the same length as the sequence.", nameof(dotBracket));
        foreach (char c in dotBracket)
            if (c != '.' && c != '(' && c != ')')
                throw new ArgumentException("dotBracket may contain only '.', '(' and ')'.", nameof(dotBracket));

        var count = new double[32];
        for (int k = 1; k < seq.Length - 1; k++)
        {
            int nt = seq[k] switch { 'A' => 0, 'U' => 1, 'G' => 2, _ => 3 };
            int j = Array.IndexOf(MiPredTriplets, dotBracket.Substring(k - 1, 3));
            if (j >= 0)
                count[nt * 8 + j]++;
        }
        for (int i = 0; i < count.Length; i++)
            count[i] /= seq.Length;
        return count;
    }

    /// <summary>
    /// Votes of the bundled 1000-tree random forest for a 34-feature vector: the number of trees whose terminal node
    /// is class P (real pre-miRNA) and the tree count. Equals <c>randomForest::predict(model, x, type="vote",
    /// norm.votes=FALSE)[,"P"]</c>; the probability of the package is votes / trees.
    /// </summary>
    public static (int RealVotes, int TreeCount) MiPredVotes(IReadOnlyList<double> features)
    {
        ArgumentNullException.ThrowIfNull(features);
        if (features.Count != MiPredFeatureCount)
            throw new ArgumentException($"Expected {MiPredFeatureCount} features.", nameof(features));
        var x = new double[MiPredFeatureCount];
        for (int i = 0; i < x.Length; i++) x[i] = features[i];
        return (MiPredForestData.Value.CountRealVotes(x), MiPredForestData.Value.TreeCount);
    }

    /// <summary>Port of the package's FASTA clean-up: characters other than A/C/G/T/U (either case) are removed, upper case, T→U.</summary>
    private static string NormaliseMiPredQuery(string sequence)
    {
        var sb = new StringBuilder(sequence.Length);
        foreach (char c in sequence)
        {
            char u = char.ToUpperInvariant(c);
            if (u == 'T') u = 'U';
            if (u is 'A' or 'C' or 'G' or 'U')
                sb.Append(u);
        }
        return sb.ToString();
    }

    /// <summary>The bundled MiPred forest (decoded and structurally validated on first use).</summary>
    internal static MiPredForest BundledMiPredForest => MiPredForestData.Value;

    private static readonly Lazy<MiPredForest> MiPredForestData =
        new(() => MiPredForest.Decode(
            Convert.FromBase64String(string.Concat(MiPredForestBase64)),
            Convert.FromBase64String(string.Concat(MiPredTreeNodeCountsBase64)),
            MiPredFeatureCount), LazyThreadSafetyMode.ExecutionAndPublication);

    /// <summary>
    /// Decoded randomForest classification forest (two classes, numeric variables): per node a code byte
    /// (1..nvar interior = split variable, 101 terminal class N, 102 terminal class P), split point and daughters.
    /// </summary>
    internal sealed class MiPredForest
    {
        private const byte TerminalN = 101;
        private const byte TerminalP = 102;

        private readonly byte[] _code;
        private readonly double[] _split;
        private readonly int[] _left;
        private readonly int[] _right;
        private readonly int[] _treeStart;

        public int TreeCount => _treeStart.Length;

        /// <summary>Total number of nodes of all trees.</summary>
        internal int NodeCount => _code.Length;

        /// <summary>Total number of terminal nodes of all trees.</summary>
        internal int TerminalCount => _code.Count(c => c == TerminalN || c == TerminalP);

        /// <summary>Number of nodes of tree <paramref name="tree"/> (randomForest <c>ndbigtree</c>).</summary>
        internal int TreeNodeCount(int tree) => (tree + 1 < _treeStart.Length ? _treeStart[tree + 1] : _code.Length) - _treeStart[tree];

        /// <summary>Number of interior nodes on the descent paths of <paramref name="x"/> whose feature equals the node's split point (diagnostic for the "&lt;=" rule).</summary>
        internal int CountSplitPointHits(double[] x)
        {
            int hits = 0;
            for (int t = 0; t < _treeStart.Length; t++)
            {
                int k = _treeStart[t];
                while (_code[k] != TerminalN && _code[k] != TerminalP)
                {
                    if (x[_code[k] - 1] == _split[k]) hits++;
                    k = x[_code[k] - 1] <= _split[k] ? _left[k] : _right[k];
                }
            }
            return hits;
        }

        private MiPredForest(byte[] code, double[] split, int[] left, int[] right, int[] treeStart)
        {
            _code = code; _split = split; _left = left; _right = right; _treeStart = treeStart;
        }

        /// <summary>Decodes the packed node table (see <c>MiRnaAnalyzer.MiPredModel.cs</c> header) and validates every tree.</summary>
        internal static MiPredForest Decode(byte[] nodes, byte[] nodeCounts, int variableCount)
        {
            int trees = nodeCounts.Length / 2;
            var counts = new int[trees];
            int total = 0;
            for (int t = 0; t < trees; t++)
            {
                counts[t] = BinaryPrimitives.ReadUInt16LittleEndian(nodeCounts.AsSpan(2 * t, 2));
                total += counts[t];
            }

            var code = new byte[total];
            var split = new double[total];
            var left = new int[total];
            var right = new int[total];
            var start = new int[trees];
            int pos = 0, g = 0;
            for (int t = 0; t < trees; t++)
            {
                start[t] = g;
                for (int k = 0; k < counts[t]; k++, g++)
                {
                    byte c = nodes[pos++];
                    code[g] = c;
                    if (c == TerminalN || c == TerminalP)
                        continue;
                    if (c < 1 || c > variableCount)
                        throw new InvalidDataException($"Invalid node code {c} in tree {t}.");
                    split[g] = BinaryPrimitives.ReadDoubleLittleEndian(nodes.AsSpan(pos, 8));
                    int l = BinaryPrimitives.ReadUInt16LittleEndian(nodes.AsSpan(pos + 8, 2));
                    int r = BinaryPrimitives.ReadUInt16LittleEndian(nodes.AsSpan(pos + 10, 2));
                    pos += 12;
                    if (l >= counts[t] || r >= counts[t] || start[t] + l <= g || start[t] + r <= g)
                        throw new InvalidDataException($"Daughter out of range or not after its parent in tree {t}.");
                    left[g] = start[t] + l;
                    right[g] = start[t] + r;
                }
            }
            if (pos != nodes.Length)
                throw new InvalidDataException("Trailing bytes in the forest table.");
            var forest = new MiPredForest(code, split, left, right, start);
            forest.ValidateTrees(counts);
            return forest;
        }

        /// <summary>Every tree must be a binary tree: each node is reached exactly once from the root and terminals = interior + 1.</summary>
        private void ValidateTrees(int[] counts)
        {
            for (int t = 0; t < _treeStart.Length; t++)
            {
                var seen = new bool[counts[t]];
                var stack = new Stack<int>();
                stack.Push(_treeStart[t]);
                int visited = 0, terminals = 0;
                while (stack.Count > 0)
                {
                    int k = stack.Pop();
                    if (seen[k - _treeStart[t]])
                        throw new InvalidDataException($"Node reached twice in tree {t}.");
                    seen[k - _treeStart[t]] = true;
                    visited++;
                    if (_code[k] == TerminalN || _code[k] == TerminalP) { terminals++; continue; }
                    stack.Push(_left[k]);
                    stack.Push(_right[k]);
                }
                if (counts[t] == 0 || visited != counts[t] || terminals * 2 - 1 != visited)
                    throw new InvalidDataException($"Tree {t} is not a consistent binary tree.");
            }
        }

        /// <summary>Number of trees voting class P: descend from node 0, left iff x[var] ≤ split point (randomForest's predictClassTree).</summary>
        internal int CountRealVotes(double[] x)
        {
            int votes = 0;
            for (int t = 0; t < _treeStart.Length; t++)
            {
                int k = _treeStart[t];
                while (_code[k] != TerminalN && _code[k] != TerminalP)
                    k = x[_code[k] - 1] <= _split[k] ? _left[k] : _right[k];
                if (_code[k] == TerminalP)
                    votes++;
            }
            return votes;
        }
    }
}
