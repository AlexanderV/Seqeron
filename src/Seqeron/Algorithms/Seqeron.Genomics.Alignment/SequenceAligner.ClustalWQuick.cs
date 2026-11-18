using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region ClustalW fast (k-tuple) pairwise stage (-QUICKTREE)

    // ClustalW 2.1 "fast" pairwise stage, FastPairwiseAlign::pairwiseAlign / pairAlign / makeNPtrs / putFrag / fragRelPos /
    // desQuickSort (src/pairwise/FastPairwiseAlign.cpp of the clustalw 2.1+lgpl-7 Debian source package, opened via
    // sources.debian.org; Wilbur & Lipman 1983 PNAS 80:726; Larkin et al. 2007 Bioinformatics 23:2947). Selected by
    // `clustalw -QUICKTREE` (and the interactive "fast" pairwise menu). It is NOT an estimate of the full-DP identity: it
    // scores the best chain of exactly matching k-mers ("fragments") on the TOPDIAGS diagonals with the most k-mer hits (and
    // the WINDOW diagonals around them), and ClustalW converts that fragment score to a percent: 100 * score / MIN(len1,
    // len2), distance = (100 - percent) / 100 in double precision (the full mode computes it in single precision).
    //
    // Everything that decides the number is ported step by step, because the result depends on details that no description
    // gives: the quicksort (not stable; it decides which diagonals are "top" among equal counts and leaves the last
    // diagonal, index l1 + l2 - 1, out of the sort range), the fragment list ordering of putFrag, the A-loop chaining rule
    // (gap penalty PAIRGAP, same-diagonal extension), and the "(Partial alignment)" early exit when more than
    // 2 * maxAlnLength fragments are created (maxAlnLength = 2 * length of the longest sequence of the whole input, which
    // is why the public method takes the input's longest length).
    //
    // DNA parameters (UserParameters.cpp): KTUP 2, WINDOW 4, SIGNIF (-TOPDIAGS) 4, WINDOWGAP (-PAIRGAP) 5; the command line
    // sets WINDOWGAP = KTUP + 4 whenever -KTUPLE is given (CommandLineParser.cpp), so `-KTUPLE=2` alone gives pair gap 6.
    // Residues are mapped A=0, C=1, G=2, T=U=3; every other IUB code (N, R, Y, ...) has no k-mer code (ziAA2DNA = -1) and
    // breaks the k-mers that contain it.

    private const int ClustalWQuickMaxKtuple = 4;

    // ziAA2DNA of FastPairwiseAlign.cpp over the ClustalW residue codes "ABCDEFGHIKLMNOPQRSTUVWXYZ-".
    private static int ClustalWQuickDnaCode(int clustalWCode) => clustalWCode switch
    {
        0 => 0,    // A
        2 => 1,    // C
        6 => 2,    // G
        18 or 19 => 3,   // T, U
        _ => -1,
    };

    /// <summary>
    /// ClustalW 2.1 fast pairwise stage for DNA (<c>-QUICKTREE</c>; k-tuple diagonal scoring, Wilbur &amp; Lipman 1983) with the percent
    /// score and guide-tree distance ClustalW derives from it.
    /// </summary>
    /// <remarks>
    /// Exact port of <c>FastPairwiseAlign::pairAlign</c> and the score conversion of <c>pairwiseAlign</c>; verified against the
    /// installed <c>clustalw -QUICKTREE -ALIGN</c> (<c>Sequences (i:j) Aligned. Score:</c>, printed as <c>%lg</c>) on the oracle in
    /// <c>clustalw_quick_oracle.tsv</c>. The score is the best k-tuple fragment chain, not an alignment identity.
    /// </remarks>
    /// <param name="sequence1">First sequence (A, C, G, T/U and the IUB ambiguity codes; case-insensitive).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="options">k-tuple parameters (default: ClustalW's DNA defaults).</param>
    /// <param name="maxSequenceLength">
    /// Length of the longest sequence of the whole input set. ClustalW bounds the fragment list by 2 · 2 · this value and
    /// stops with "(Partial alignment)" beyond it; <c>null</c> = the longer of the two sequences (a two-sequence input).
    /// </param>
    public static ClustalWQuickPairwiseResult ClustalWQuickPairwiseAlign(
        string sequence1, string sequence2, ClustalWQuickOptions? options = null, int? maxSequenceLength = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        options ??= new ClustalWQuickOptions();
        ValidateClustalWQuickOptions(options);

        var codes1 = EncodeClustalWDna(sequence1, nameof(sequence1));
        var codes2 = EncodeClustalWDna(sequence2, nameof(sequence2));
        int longest = Math.Max(codes1.Length, codes2.Length);
        if (maxSequenceLength is { } given && given < longest)
            throw new ArgumentOutOfRangeException(nameof(maxSequenceLength), given, "Must be at least the length of the longer sequence.");

        // An empty sequence cannot be read by ClustalW; mirror the full mode (distance 1, not aligned).
        if (codes1.Length == 0 || codes2.Length == 0)
            return new ClustalWQuickPairwiseResult(0, 0.0, 1.0, false, false);

        int maxAlnLength = 2 * (maxSequenceLength ?? longest);
        var (fragmentScore, partial) = ClustalWQuickFragmentScore(codes1, codes2, options, maxAlnLength);

        // FastPairwiseAlign::pairwiseAlign: percent => score / MIN(len1, len2) * 100; distance = (100 - score) / 100.
        double calc = 0.0;
        if (fragmentScore != 0)
        {
            calc = fragmentScore;
            if (options.Percent)
                calc = calc / Math.Min(codes1.Length, codes2.Length) * 100.0;
        }

        double distance = (100.0 - calc) / 100.0;
        return new ClustalWQuickPairwiseResult(fragmentScore, calc, distance, calc > 0.1, partial);
    }

    /// <summary><see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/> for <see cref="DnaSequence"/> inputs.</summary>
    public static ClustalWQuickPairwiseResult ClustalWQuickPairwiseAlign(
        DnaSequence sequence1, DnaSequence sequence2, ClustalWQuickOptions? options = null, int? maxSequenceLength = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWQuickPairwiseAlign(sequence1.Sequence, sequence2.Sequence, options, maxSequenceLength);
    }

    private static void ValidateClustalWQuickOptions(ClustalWQuickOptions options)
    {
        if (options.KTuple < 1 || options.KTuple > ClustalWQuickMaxKtuple)
            throw new ArgumentOutOfRangeException(nameof(options), options.KTuple, $"ClustalW's DNA k-tuple size KTuple is 1..{ClustalWQuickMaxKtuple}.");
        if (options.Window < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.Window, "Window must be positive.");
        if (options.TopDiagonals < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.TopDiagonals, "TopDiagonals must be positive.");
        if (options.PairGap < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.PairGap, "PairGap must be positive.");
    }

    /// <summary>
    /// Distance of the fast mode as stored in ClustalW's double-precision distance matrix; <paramref name="maxSequenceLength"/> is
    /// the longest sequence of the input set (see <see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>).
    /// </summary>
    private static double ClustalWQuickPairwiseDistance(string a, string b, ClustalWQuickOptions options, int maxSequenceLength)
        => ClustalWQuickPairwiseAlign(a, b, options, maxSequenceLength).Distance;

    // Port of FastPairwiseAlign::pairAlign + makeNPtrs; all arrays 1-based as in the C++ source.
    private static (int FragmentScore, bool Partial) ClustalWQuickFragmentScore(
        int[] seq1, int[] seq2, ClustalWQuickOptions options, int maxAlnLength)
    {
        int l1 = seq1.Length, l2 = seq2.Length;
        int ktup = options.KTuple;
        int window = options.Window;
        int signif = options.TopDiagonals;
        int windowGap = options.PairGap;

        // makeNPtrs: for every k-mer code, the position list (most recent position first) of each sequence.
        int limit = 1 << (2 * ktup);   // 4^ktup
        var zza = new int[l1 + 2];
        var zzb = new int[l2 + 2];
        var zzc = new int[limit + 2];
        var zzd = new int[limit + 2];
        BuildClustalWQuickPointers(seq1, ktup, zza, zzc);
        BuildClustalWQuickPointers(seq2, ktup, zzb, zzd);

        int num = Math.Max(2 * maxAlnLength + 1, l1 + l2 + 2);
        var displ = new int[num + 1];
        var slopes = new int[num + 1];
        var diagIndex = new int[num + 1];
        var accum = new int[5][];
        for (int r = 0; r < 5; r++)
            accum[r] = new int[num + 1];

        int tl1 = l1 + l2 - 1;
        for (int i = 1; i <= tl1; i++)
        {
            slopes[i] = displ[i] = 0;
            diagIndex[i] = i;
        }

        // increment the diagonal score for each k-tuple match
        for (int i = 1; i <= limit; i++)
        {
            int vn1 = zzc[i];
            while (vn1 != 0)
            {
                int vn2 = zzd[i];
                while (vn2 != 0)
                {
                    displ[vn1 - vn2 + l2]++;
                    vn2 = zzb[vn2];
                }

                vn1 = zza[vn1];
            }
        }

        // choose the top SIGNIF diagonals (the quicksort leaves index tl1 out of its range, as in ClustalW)
        ClustalWDesQuickSort(displ, diagIndex, tl1);

        int jTop = tl1 - signif + 1;
        if (jTop < 1)
            jTop = 1;

        // flag all diagonals within WINDOW of a top diagonal
        for (int i = tl1; i >= jTop; i--)
        {
            if (displ[i] > 0)
            {
                int pos = diagIndex[i];
                int lo = Math.Max(1, pos - window);
                int hi = Math.Min(tl1, pos + window);
                for (; lo <= hi; lo++)
                    slopes[lo] = 1;
            }
        }

        for (int i = 1; i <= tl1; i++)
            displ[i] = 0;

        int currFrag = 0, maxSoFar = 0;
        int pot2 = 1;   // pot[j] = 4^(j-1)

        for (int i = 1; i <= l1 - ktup + 1; i++)
        {
            int encrypt = 0;
            bool skip = false;
            pot2 = 1;
            for (int j = 1; j <= ktup; j++)
            {
                int residue = ClustalWQuickDnaCode(seq1[i + j - 2]);   // seq1 is 0-based here: (*seq)[i + j - 1] 1-based
                if (residue < 0)
                {
                    skip = true;
                    break;
                }

                encrypt += residue * pot2;
                pot2 *= 4;
            }

            if (skip)
                continue;
            encrypt++;

            int vn2 = zzd[encrypt];
            while (vn2 != 0)
            {
                int osptr = i - vn2 + l2;
                if (slopes[osptr] != 1)
                {
                    vn2 = zzb[vn2];
                    continue;
                }

                int flen = 0;
                int fs = ktup;
                int next = maxSoFar;

                // A-loop
                while (true)
                {
                    if (next == 0)
                    {
                        currFrag++;
                        if (currFrag >= 2 * maxAlnLength)
                            return (accum[0][maxSoFar], true);   // "(Partial alignment)"; maxSoFar is what it was

                        displ[osptr] = currFrag;
                        PutFrag(fs, i, vn2, flen);
                        break;
                    }

                    int tv1 = accum[1][next];
                    int tv2 = accum[2][next];
                    if (FragRelPos(i, vn2, tv1, tv2))
                    {
                        if (i - vn2 == accum[1][next] - accum[2][next])
                        {
                            fs = i > accum[1][next] + (ktup - 1)
                                ? accum[0][next] + ktup
                                : accum[0][next] + (i - accum[1][next]);
                            flen = next;
                            next = 0;
                            continue;
                        }

                        int subt1;
                        if (displ[osptr] == 0)
                        {
                            subt1 = ktup;
                        }
                        else
                        {
                            subt1 = i > accum[1][displ[osptr]] + (ktup - 1)
                                ? accum[0][displ[osptr]] + ktup
                                : accum[0][displ[osptr]] + (i - accum[1][displ[osptr]]);
                        }

                        int subt2 = accum[0][next] - windowGap + ktup;
                        if (subt2 > subt1)
                        {
                            flen = next;
                            fs = subt2;
                        }
                        else
                        {
                            flen = displ[osptr];
                            fs = subt1;
                        }

                        next = 0;
                        continue;
                    }

                    next = accum[4][next];
                }

                vn2 = zzb[vn2];
            }
        }

        return (maxSoFar == 0 ? 0 : accum[0][maxSoFar], false);

        // putFrag: insert the new fragment into the list ordered by decreasing score (accum[4] = link to the next lower one).
        void PutFrag(int fs, int v1, int v2, int flen)
        {
            accum[0][currFrag] = fs;
            accum[1][currFrag] = v1;
            accum[2][currFrag] = v2;
            accum[3][currFrag] = flen;

            if (maxSoFar == 0)
            {
                maxSoFar = 1;
                accum[4][currFrag] = 0;
                return;
            }

            if (fs >= accum[0][maxSoFar])
            {
                accum[4][currFrag] = maxSoFar;
                maxSoFar = currFrag;
                return;
            }

            int nxt = maxSoFar;
            int end;
            while (true)
            {
                end = nxt;
                nxt = accum[4][nxt];
                if (fs >= accum[0][nxt])
                    break;
            }

            accum[4][currFrag] = nxt;
            accum[4][end] = currFrag;
        }

        bool FragRelPos(int a1, int b1, int a2, int b2)
        {
            if (a1 - b1 == a2 - b2)
                return a2 < a1;
            return a2 + ktup - 1 < a1 && b2 + ktup - 1 < b1;
        }
    }

    // makeNPtrs: tptr[i] = previous position with the same k-mer code (0 = none); pl[code] = latest position (1-based).
    private static void BuildClustalWQuickPointers(int[] seq, int ktup, int[] tptr, int[] pl)
    {
        int length = seq.Length;
        for (int i = 1; i <= length - ktup + 1; i++)
        {
            int code = 0, pot = 1;
            bool skip = false;
            for (int j = 1; j <= ktup; j++)
            {
                int residue = ClustalWQuickDnaCode(seq[i + j - 2]);
                if (residue < 0)
                {
                    skip = true;
                    break;
                }

                code += residue * pot;
                pot *= 4;
            }

            if (skip)
                continue;
            code++;
            if (pl[code] != 0)
                tptr[i] = pl[code];
            pl[code] = i;
        }
    }

    // desQuickSort: Kernighan & Plauger quicksort, ascending in array1 with array2 carried along; sorts indices 1..arraySize-1.
    private static void ClustalWDesQuickSort(int[] array1, int[] array2, int arraySize)
    {
        var lst = new int[64];
        var ust = new int[64];
        lst[1] = 1;
        ust[1] = arraySize - 1;
        int p = 1;

        while (p > 0)
        {
            if (lst[p] >= ust[p])
            {
                p--;
                continue;
            }

            int i = lst[p] - 1;
            int j = ust[p];
            int pivlin = array1[j];
            while (i < j)
            {
                i++;
                while (array1[i] < pivlin)
                    i++;

                for (j--; j > i; j--)
                {
                    if (array1[j] <= pivlin)
                        break;
                }

                if (i < j)
                {
                    (array1[i], array1[j]) = (array1[j], array1[i]);
                    (array2[i], array2[j]) = (array2[j], array2[i]);
                }
            }

            j = ust[p];
            (array1[i], array1[j]) = (array1[j], array1[i]);
            (array2[i], array2[j]) = (array2[j], array2[i]);

            if (i - lst[p] < ust[p] - i)
            {
                lst[p + 1] = lst[p];
                ust[p + 1] = i - 1;
                lst[p] = i + 1;
            }
            else
            {
                lst[p + 1] = i + 1;
                ust[p + 1] = ust[p];
                ust[p] = i - 1;
            }

            p++;
        }
    }

    #endregion
}
