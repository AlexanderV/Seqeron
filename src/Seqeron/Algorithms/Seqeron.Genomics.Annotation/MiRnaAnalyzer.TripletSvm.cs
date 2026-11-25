using System.Text;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Published pre-miRNA classifier: Triplet-SVM (Xue C, Li F, He T, Liu G-P, Li Y, Zhang X. 2005,
/// "Classification of real and pseudo microRNA precursors using local structure-sequence features and
/// support vector machine", BMC Bioinformatics 6:310, doi:10.1186/1471-2105-6-310) with the authors'
/// own trained libsvm model (see <c>MiRnaAnalyzer.TripletSvmModel.cs</c>). OPT-IN: the in-repo
/// logistic model behind <see cref="ClassifyPreMiRna"/> is unchanged.
/// </summary>
public static partial class MiRnaAnalyzer
{
    /// <summary>Number of triplet-element features (4 centre nucleotides × 8 paired/unpaired neighbour states).</summary>
    public const int TripletSvmFeatureCount = 32;

    /// <summary>Default minimum number of paired bases per stem arm of a retrieved stem-loop (Xue's <c>2_get_stemloop.pl</c> argument, the package example uses 22).</summary>
    public const int TripletSvmDefaultMinBasePairs = 22;

    /// <summary>
    /// One stem-loop candidate retrieved from a query and scored by the Triplet-SVM model.
    /// </summary>
    /// <param name="Start">1-based start of the stem-loop in the query (the package's <c>&gt;n_start_end</c> title).</param>
    /// <param name="End">1-based end of the stem-loop in the query.</param>
    /// <param name="Sequence">Stem-loop sequence (RNA, upper case).</param>
    /// <param name="DotBracket">Stem-loop structure.</param>
    /// <param name="Features">The 32 triplet-element frequencies (sum 1; index = nucleotide(A,G,C,U) × 8 + neighbour bits).</param>
    /// <param name="DecisionValue">libsvm decision value Σ y·α·exp(−16‖x−sv‖²) − ρ.</param>
    /// <param name="IsPreMiRna">True iff <paramref name="DecisionValue"/> &gt; 0 (libsvm label +1 = real pre-miRNA).</param>
    public readonly record struct TripletSvmResult(
        int Start,
        int End,
        string Sequence,
        string DotBracket,
        double[] Features,
        double DecisionValue,
        bool IsPreMiRna);

    /// <summary>
    /// Classifies a query with the published Triplet-SVM model (Xue et al. 2005): retrieves the stem-loops
    /// with at least <paramref name="minBasePairs"/> paired bases per arm (<c>2_get_stemloop.pl</c>), codes each
    /// single stem-loop as 32 triplet-element frequencies (<c>3_step_triplet_coding_for_queries.pl</c>) and
    /// evaluates the authors' RBF-SVM (<c>trainset_hsa163_cds168_unite.txt.model</c>, libsvm). The three
    /// package steps are ported line by line and cross-checked against the Perl scripts and libsvm
    /// <c>svm-predict</c>; the only difference from the package is the folder: the structure comes from the
    /// canonical MFE fold of this library (<see cref="RnaSecondaryStructure.CalculateMfeStructure"/>) unless the
    /// caller supplies one (the package used RNAfold).
    /// </summary>
    /// <param name="sequence">Query (RNA or DNA; T read as U). A query with other characters is dropped like the package's step 1 (empty result).</param>
    /// <param name="dotBracket">Optional structure in <c>.()</c> notation (same length); folded with the canonical MFE engine when null.</param>
    /// <param name="minBasePairs">Minimum paired bases per arm of a retrieved stem-loop (default 22).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size of the fold (only used when folding; values below the NNDB
    /// minimum of 3 are raised to 3 by the folder; must be non-negative).</param>
    /// <returns>One result per retrieved stem-loop that has a single arm pair (the package codes only those); empty otherwise.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minBasePairs"/> is &lt; 1 (a stem arm needs at
    /// least one paired base; 0 or a negative value would admit every unpaired piece as a "stem-loop") or
    /// <paramref name="minLoopSize"/> is negative.</exception>
    public static IReadOnlyList<TripletSvmResult> ClassifyPreMiRnaTripletSvm(
        string sequence,
        string? dotBracket = null,
        int minBasePairs = TripletSvmDefaultMinBasePairs,
        int minLoopSize = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(minBasePairs, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopSize);
        var results = new List<TripletSvmResult>();
        string? seq = ToAcguRna(sequence, dropOthers: false);
        if (seq is null || seq.Length == 0)
            return results;

        string structure;
        if (dotBracket is null)
        {
            structure = RnaSecondaryStructure.CalculateMfeStructure(seq, minLoopSize).DotBracket;
        }
        else
        {
            if (dotBracket.Length != seq.Length)
                throw new ArgumentException("dotBracket must have the same length as the sequence.", nameof(dotBracket));
            foreach (char c in dotBracket)
                if (c != '.' && c != '(' && c != ')')
                    throw new ArgumentException("dotBracket may contain only '.', '(' and ')'.", nameof(dotBracket));
            structure = dotBracket;
        }

        foreach (var (start, end, candSeq, candStruct) in RetrieveTripletStemLoops(seq, structure, minBasePairs))
        {
            double[]? features = ComputeTripletFeatures(candSeq, candStruct);
            if (features is null)
                continue;
            double dv = TripletSvmDecisionValue(features);
            results.Add(new TripletSvmResult(start, end, candSeq, candStruct, features, dv, dv > 0));
        }
        return results;
    }

    /// <summary>
    /// The 32 triplet-element features of a single stem-loop (package step 3, <c>reverse_flag = 0</c>,
    /// <c>method = 1</c>): for every nucleotide of the 5′ arm and of the 3′ arm (loop and flanking unpaired
    /// ends excluded) the centre nucleotide (A=0, G=1, C=2, U=3) and the paired(1)/unpaired(0) state of
    /// its two arm neighbours (arm ends padded with unpaired) give the index nucleotide × 8 + bits; counts
    /// are normalised to sum 1.
    /// </summary>
    /// <returns>The features, or <c>null</c> when the structure is not a single arm pair or has non-ACGU characters.</returns>
    public static double[]? ComputeTripletFeatures(string sequence, string dotBracket)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(dotBracket);
        string? seq = ToAcguRna(sequence, dropOthers: false);
        if (seq is null || seq.Length != dotBracket.Length)
            return null;

        var pieces = DivideTripletStructure(dotBracket);
        // get_main_struct: only a structure with exactly one (open-group, close-group) pair is coded.
        int pairs = (pieces.Count + 1) / 2;
        int same = 0;
        for (int i = 0; i < pieces.Count; i += 2)
            if (pieces[i].AsSpan().Count('(') == (i + 1 < pieces.Count ? pieces[i + 1] : "").AsSpan().Count(')'))
                same++;
        if (same != pairs || same != 1 || pieces.Count != 2)
            return null;

        string leftArm = TripletBody(pieces[0]);
        string rightArm = TripletBody(pieces[1]);
        int leftPos = TripletFind(leftArm, dotBracket, 0);
        int rightPos = TripletFind(rightArm, dotBracket, 0);
        string leftSeq = Substr(seq, leftPos, leftArm.Length);
        string rightSeq = Substr(seq, rightPos, rightArm.Length);

        var table = new double[TripletSvmFeatureCount];
        CountTriplets(leftSeq, leftArm, table);
        CountTriplets(rightSeq, rightArm, table);
        double sum = 0;
        foreach (double v in table) sum += v;
        if (sum == 0)
            return null;
        for (int i = 0; i < table.Length; i++)
            table[i] /= sum;
        return table;
    }

    /// <summary>
    /// libsvm decision value of the bundled Triplet-SVM model for a 32-element feature vector:
    /// Σ<sub>i</sub> coef<sub>i</sub>·exp(−γ‖x − sv<sub>i</sub>‖²) − ρ with γ = 16, ρ = 2.74269 (positive = real pre-miRNA).
    /// </summary>
    public static double TripletSvmDecisionValue(IReadOnlyList<double> features)
    {
        ArgumentNullException.ThrowIfNull(features);
        if (features.Count != TripletSvmFeatureCount)
            throw new ArgumentException($"Expected {TripletSvmFeatureCount} features.", nameof(features));
        double sum = 0;
        for (int s = 0; s < TripletSvmSupportVectorCount; s++)
        {
            double d2 = 0;
            int off = s * TripletSvmFeatureCount;
            for (int k = 0; k < TripletSvmFeatureCount; k++)
            {
                double d = features[k] - TripletSvmSupportVectors[off + k];
                d2 += d * d;
            }
            sum += TripletSvmCoefficients[s] * Math.Exp(-TripletSvmGamma * d2);
        }
        return sum - TripletSvmRho;
    }

    // ── Port of 2_get_stemloop.pl ─────────────────────────────────────────────────────────────

    private static IEnumerable<(int Start, int End, string Seq, string Struct)> RetrieveTripletStemLoops(
        string seq, string symbol, int minBasePairs)
    {
        var pieces = DivideTripletStructure(symbol);
        if (pieces.Count < 2)
            yield break;
        var cand = GetTripletCandidateStructs(pieces, minBasePairs);
        int rightPos = 0;
        for (int i = 0; i + 1 < cand.Count; i += 2)
        {
            int leftPos = TripletFind(cand[i], symbol, rightPos);
            rightPos = TripletFind(cand[i + 1], symbol, leftPos);
            int len = rightPos - leftPos + cand[i + 1].Length;
            yield return (1 + leftPos, leftPos + len, Substr(seq, leftPos, len), Substr(symbol, leftPos, len));
        }
    }

    private static List<string> GetTripletCandidateStructs(List<string> pieces, int minBp)
    {
        var left = new List<int>();
        var right = new List<int>();
        for (int i = 0; i < pieces.Count; i += 2)
        {
            left.Add(pieces[i].AsSpan().Count('('));
            right.Add((i + 1 < pieces.Count ? pieces[i + 1] : "").AsSpan().Count(')'));
        }

        var ret = new List<string>();
        for (int i = 0, j = 0; i < left.Count; i++, j += 2)
        {
            if (left[i] < minBp || right[i] < minBp)
                continue;
            string a = pieces[j];
            string b = j + 1 < pieces.Count ? pieces[j + 1] : "";
            if (left[i] == right[i])
            {
                ret.Add(a);
                ret.Add(b);
            }
            else if (left[i] > right[i])
            {
                // trim the surplus 5' pairs: keep from the right[i]-th '(' counted from the end
                int numBp = 0, cutPos = 0;
                for (int k = a.Length - 1; k >= 0; k--)
                {
                    if (numBp != right[i])
                    {
                        if (a[k] == '(') numBp++;
                    }
                    else
                    {
                        cutPos = k + 1;
                        break;
                    }
                }
                ret.Add(a.Substring(cutPos));
                ret.Add(b);
            }
            else
            {
                // trim the surplus 3' pairs: package loop runs k < len-1 and leaves the length undefined (0) when not reached
                int numBp = 0, cutLength = 0;
                for (int k = 0; k < b.Length - 1; k++)
                {
                    if (numBp != left[i])
                    {
                        if (b[k] == ')') numBp++;
                    }
                    else
                    {
                        cutLength = k;
                        break;
                    }
                }
                ret.Add(a);
                ret.Add(b.Substring(0, cutLength));
            }
        }

        for (int i = 0; i < ret.Count; i++)
            ret[i] = TripletBody(ret[i]);
        return ret;
    }

    // ── helpers shared by the ports ───────────────────────────────────────────────────────────

    /// <summary><c>divide_structure</c>: split at every change between '(' and ')' groups (dots stay with the preceding piece).</summary>
    private static List<string> DivideTripletStructure(string symbol)
    {
        var starts = new List<int> { 0 };
        char flag = '(';
        for (int i = 0; i < symbol.Length; i++)
        {
            char c = symbol[i];
            if (c != '.' && c != flag)
            {
                flag = c;
                starts.Add(i);
            }
        }
        var pieces = new List<string>();
        for (int i = 0; i < starts.Count - 1; i++)
            pieces.Add(symbol.Substring(starts[i], starts[i + 1] - starts[i]));
        pieces.Add(symbol.Substring(starts[^1]));
        return pieces;
    }

    /// <summary><c>get_match_all_body</c>: strip leading and trailing dots.</summary>
    private static string TripletBody(string s) => s.Trim('.');

    private static int TripletFind(string piece, string symbol, int from)
    {
        for (int i = from; i < symbol.Length; i++)
            if (string.CompareOrdinal(symbol, i, piece, 0, piece.Length) == 0 && i + piece.Length <= symbol.Length)
                return i;
        return 0; // Perl: undefined position in numeric context
    }

    private static string Substr(string s, int pos, int len)
    {
        if (pos < 0 || pos > s.Length) return "";
        return s.Substring(pos, Math.Min(len, s.Length - pos));
    }

    private static void CountTriplets(string armSeq, string arm, double[] table)
    {
        for (int i = 0; i < arm.Length; i++)
        {
            int nt = armSeq[i] switch { 'A' => 0, 'G' => 1, 'C' => 2, _ => 3 };
            string near3;
            if (i == 0)
                near3 = "." + Substr(arm, 0, 2);
            else if (i == arm.Length - 1)
                near3 = Substr(arm, i - 1, 2) + ".";
            else
                near3 = arm.Substring(i - 1, 3);
            int bits = 0;
            foreach (char c in near3)
                bits = (bits << 1) | (c == '.' ? 0 : 1);
            table[nt * 8 + bits]++;
        }
    }
}
