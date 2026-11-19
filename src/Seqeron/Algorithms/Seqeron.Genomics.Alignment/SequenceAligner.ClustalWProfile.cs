using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    // ClustalW 2.1 profile alignment commands (audit A2-14, B09 F103): `-PROFILE1=a -PROFILE2=b -PROFILE` (Clustal::profileAlign: two alignments merged by profile-profile alignment) and
    // `-PROFILE1=a -PROFILE2=b -SEQUENCES` (Clustal::sequencesAlignToProfile, MSA::seqsAlignToProfile: the sequences of b are added one by one to the alignment a), on the machinery of the
    // multiple-alignment stage (MyersMillerProfileAlign::profileAlign, MSA::calcPairwiseForProfileAlign / doProfileAlign, TreeInterface::getWeightsForProfileAlign*) of clustalw 2.1+lgpl-7.

    /// <summary>
    /// ClustalW 2.1's profile alignment of two existing alignments (<c>clustalw -TYPE=DNA -PROFILE1=a -PROFILE2=b -PROFILE</c>): the rows of both alignments are kept (gaps '-' allowed, a profile may be
    /// a single sequence), the guide trees of the two profiles weight their sequences, an initial alignment with equal weights measures the identities between the profiles, and the profiles are
    /// aligned by Myers–Miller profile-profile alignment with position-specific gap penalties (<c>Clustal::profileAlign</c> + <c>MSA::doProfileAlign</c>). Rows are returned: profile 1, then profile 2.
    /// </summary>
    /// <remarks>
    /// Verified against the installed <c>clustalw</c> 2.1 on the sets of <c>clustalw_profile_oracle.tsv</c> (profiles produced by <c>clustalw -ALIGN</c> itself). Reproduced defects of the
    /// program: the neighbour-joining tree of a profile 2 that does not start at sequence 1 is built from the identities of the FIRST sequences of the alignment (<c>NJTree::generateTree</c> reads the
    /// entries 1..k of the matrix), and the identities inside the profiles are truncated to integers. <see cref="ClustalWOptions.Iteration"/> = <see cref="ClustalWIteration.Tree"/> is refused: the
    /// program then reads <c>prof1Weight[1]</c> of a one-element vector (undefined behaviour).
    /// </remarks>
    /// <param name="profile1">Rows of the first alignment (at least one).</param>
    /// <param name="profile2">Rows of the second alignment (at least one).</param>
    /// <param name="options">ClustalW parameters (default: ClustalW's DNA defaults).</param>
    /// <param name="scoring">Scoring used for the reported total score only (default <see cref="SimpleDna"/>, affine).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult AlignClustalWProfiles(
        IEnumerable<string> profile1, IEnumerable<string> profile2, ClustalWOptions? options = null, ScoringMatrix? scoring = null)
    {
        options ??= new ClustalWOptions();
        var (p1, p2) = ProfileInputs(profile1, profile2);
        ValidateClustalWOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = p1.Concat(p2).Select(r => EncodeClustalWDna(r, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options));
        return BuildClustalWResult(state.RunProfileAlignment(p1.Count), scoring ?? SimpleDna);
    }

    /// <summary>The PROTEIN counterpart of <see cref="AlignClustalWProfiles(IEnumerable{string}, IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/> (<c>clustalw -TYPE=PROTEIN -PROFILE1=a -PROFILE2=b -PROFILE</c>).</summary>
    public static MultipleAlignmentResult AlignClustalWProfiles(
        IEnumerable<string> profile1, IEnumerable<string> profile2, ClustalWProteinOptions options, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        var (p1, p2) = ProfileInputs(profile1, profile2);
        ValidateClustalWProteinMultipleOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = p1.Concat(p2).Select(r => EncodeClustalWProtein(r, "profiles", allowGaps: true)).ToArray();
        string hyd = new string(options.HydrophilicResidues.Take(ClustalWMaxHydrophilicResidues).TakeWhile(char.IsAsciiLetter).Select(char.ToUpperInvariant).ToArray());
        var state = new ClustalWMsaState(codes, CwConfig.FromProtein(options) with { HydrophilicResidues = hyd.Length > 0 ? hyd : "GPSNDQEKR" });
        return BuildClustalWProteinResult(state.RunProfileAlignment(p1.Count), scoring);
    }

    /// <summary>
    /// ClustalW 2.1's sequence-to-profile alignment (<c>clustalw -TYPE=DNA -PROFILE1=a -PROFILE2=b -SEQUENCES</c>, <c>Clustal::sequencesAlignToProfile</c> / <c>MSA::seqsAlignToProfile</c>): the sequences are
    /// added to the alignment <paramref name="profile"/> one at a time, the one with the highest identity to the aligned set first. Pairwise distances: identities inside the profile (<c>countid</c>),
    /// the full or fast pairwise stage (<see cref="ClustalWOptions.PairwiseMode"/>) between every other pair of rows (computed on the rows with their gaps), then the neighbour-joining tree of all rows
    /// weights them; each sequence is aligned with weights = tree weight x identity to it. Rows are returned: profile, then sequences (input order).
    /// </summary>
    /// <remarks>Verified against the installed <c>clustalw</c> 2.1 on the sets of <c>clustalw_profile_oracle.tsv</c>. <see cref="ClustalWOptions.Iteration"/> = <see cref="ClustalWIteration.Tree"/> is refused (see <see cref="AlignClustalWProfiles(IEnumerable{string}, IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/>).</remarks>
    /// <param name="profile">Rows of the existing alignment (at least one).</param>
    /// <param name="sequences">The sequences to add (at least one; gaps '-' are kept as in the program).</param>
    /// <param name="options">ClustalW parameters (default: ClustalW's DNA defaults).</param>
    /// <param name="scoring">Scoring used for the reported total score only (default <see cref="SimpleDna"/>, affine).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult AlignClustalWSequencesToProfile(
        IEnumerable<string> profile, IEnumerable<string> sequences, ClustalWOptions? options = null, ScoringMatrix? scoring = null)
    {
        options ??= new ClustalWOptions();
        var (p1, p2) = ProfileInputs(profile, sequences);
        ValidateClustalWOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = p1.Concat(p2).Select(r => EncodeClustalWDna(r, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options));
        int longest = codes.Max(c => c.Length);
        Func<int[], int[], double> pairwise;
        if (options.PairwiseMode == ClustalWPairwiseMode.Quick)
        {
            var quick = options.Quick ?? new ClustalWQuickOptions();
            pairwise = (a, b) => ClustalWQuickPairwiseAlignCore(a, b, protein: false, quick, longest).Distance;
        }
        else
        {
            var setup = ClustalWDnaPairwiseSetup(new ClustalWDnaPairwiseOptions
            {
                Matrix = options.PairwiseDnaMatrix,
                UserMatrix = options.PairwiseUserDnaMatrix,
                UserMatrixUse = options.UserDnaMatrixUse,
                TransitionWeight = options.TransitionWeight,
            });
            float go = (float)options.PairwiseGapOpen, ge = (float)options.PairwiseGapExtend;
            int gapOpenInt = (int)(2f * go * (float)ClustalWIntScale * setup.GapOpenScale);
            int gapExtendInt = (int)(ge * (float)ClustalWIntScale * setup.GapExtendScale);
            pairwise = (a, b) => ClustalWGappedPairwiseDistance(a, b, setup.Matrix, gapOpenInt, gapExtendInt);
        }

        return BuildClustalWResult(state.RunSequencesToProfile(p1.Count, pairwise), scoring ?? SimpleDna);
    }

    /// <summary>The PROTEIN counterpart of <see cref="AlignClustalWSequencesToProfile(IEnumerable{string}, IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/> (<c>clustalw -TYPE=PROTEIN ... -SEQUENCES</c>).</summary>
    public static MultipleAlignmentResult AlignClustalWSequencesToProfile(
        IEnumerable<string> profile, IEnumerable<string> sequences, ClustalWProteinOptions options, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        var (p1, p2) = ProfileInputs(profile, sequences);
        ValidateClustalWProteinMultipleOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = p1.Concat(p2).Select(r => EncodeClustalWProtein(r, "profiles", allowGaps: true)).ToArray();
        string hyd = new string(options.HydrophilicResidues.Take(ClustalWMaxHydrophilicResidues).TakeWhile(char.IsAsciiLetter).Select(char.ToUpperInvariant).ToArray());
        var state = new ClustalWMsaState(codes, CwConfig.FromProtein(options) with { HydrophilicResidues = hyd.Length > 0 ? hyd : "GPSNDQEKR" });
        int longest = codes.Max(c => c.Length);
        Func<int[], int[], double> pairwise;
        if (options.PairwiseMode == ClustalWPairwiseMode.Quick)
        {
            pairwise = (a, b) => ClustalWQuickPairwiseAlignCore(a, b, protein: true, options.Pairwise.Quick, longest).Distance;
        }
        else
        {
            var (matrix, intScale, matAvg) = options.Pairwise.UserMatrix is { } user
                ? BuildClustalWUserPairwiseMatrix(user)
                : ClustalWProteinMatrices[(int)options.Pairwise.Matrix].Value;
            float pwGapOpen = (float)options.Pairwise.GapOpen;
            float pwGapExtend = (float)options.Pairwise.GapExtend;
            int gapExtendInt = (int)(pwGapExtend * (float)intScale);
            pairwise = (a, b) =>
            {
                // FullPairwiseAlign: the gap opening depends on the ROW lengths (gaps included)
                double logMin = Math.Log((double)Math.Min(a.Length, b.Length));
                int gapOpenInt = matAvg <= 0
                    ? 2 * (int)(((double)pwGapOpen + logMin) * intScale)
                    : (int)(2 * matAvg * ((double)pwGapOpen + logMin) * 1.0);
                return ClustalWGappedPairwiseDistance(a, b, matrix, gapOpenInt, gapExtendInt);
            };
        }

        return BuildClustalWProteinResult(state.RunSequencesToProfile(p1.Count, pairwise), scoring);
    }

    private static (List<string> P1, List<string> P2) ProfileInputs(IEnumerable<string> profile1, IEnumerable<string> profile2)
    {
        var p1 = ToListNoNull(profile1, nameof(profile1));
        var p2 = ToListNoNull(profile2, nameof(profile2));
        if (p1.Count == 0 || p2.Count == 0)
            throw new ArgumentException("A profile alignment needs at least one sequence in each profile (ClustalW: 'No sequences in memory').", p1.Count == 0 ? nameof(profile1) : nameof(profile2));
        return (p1, p2);
    }

    private static void ValidateClustalWProfileOptions(ClustalWIteration iteration)
    {
        if (iteration == ClustalWIteration.Tree)
            throw new NotSupportedException("ClustalW 2.1 -ITERATION=TREE on a profile alignment reads prof1Weight[1] of a one-element vector (undefined behaviour of the program); use iteration 'none'.");
    }

    private static MultipleAlignmentResult BuildClustalWProteinResult(string[] rows, MatrixScoring? scoring)
    {
        var list = rows.ToList();
        string consensus = BuildConsensus(list, list.Max(r => r.Length), anySymbol: true);
        var sp = scoring ?? MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
        string[] scored = scoring is null ? rows.Select(r => r.Replace('O', 'X').Replace('U', 'X')).ToArray() : rows;
        double exact = rows.Length == 1 ? 0.0 : MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
        return new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
    }

    private sealed partial class ClustalWMsaState
    {
        // Clustal::profileAlign (no guide trees from files, no secondary-structure masks): distances inside the profiles (countid truncated to int), the weights / initial alignment of
        // getWeightsForProfileAlign, MSA::doProfileAlign
        public string[] RunProfileAlignment(int profile1NumSeqs)
        {
            int n = _n;
            var aln = new CwAlignment
            {
                Rows = [.. _seq.Select(r => (int[])r.Clone())],
                Ids = [.. Enumerable.Range(1, n).Select(i => (long)i)],
                Weights = Enumerable.Repeat(100, n + 1).ToArray(),
                Profile1NumSeqs = profile1NumSeqs,
            };

            var distMat = new double[n + 1, n + 1];
            if (profile1NumSeqs >= 2)
                FillProfileDistances(aln, distMat, 1, profile1NumSeqs);
            if (n - profile1NumSeqs >= 2)
                FillProfileDistances(aln, distMat, profile1NumSeqs + 1, n);

            var (prof1Weight, prof2Weight) = GetWeightsForProfileAlign(aln, distMat, duplicateNames: false);
            DoProfileAlign(aln, distMat, prof1Weight, prof2Weight);
            return aln.Rows.Select(CodesToRow).ToArray();
        }

        // Clustal::sequencesAlignToProfile + MSA::seqsAlignToProfile (no guide tree from a file, no structure masks)
        public string[] RunSequencesToProfile(int profile1NumSeqs, Func<int[], int[], double> pairwise)
        {
            int n = _n;
            int p1 = profile1NumSeqs;
            var aln = new CwAlignment { Rows = [.. _seq.Select(r => (int[])r.Clone())], Weights = Enumerable.Repeat(100, n + 1).ToArray() };
            var dist = new double[n + 1, n + 1];
            int newSeq = p1 + 1;
            for (int i = 1; i <= newSeq; i++)
            {
                for (int j = i + 1; j <= newSeq; j++)
                {
                    float dscore = aln.CountId(i, j);
                    dist[i, j] = dist[j, i] = (100.0 - (double)dscore) / 100.0;
                }
            }

            // pairwiseAlign(.., 0, numSeqs, newSeq - 2, numSeqs): every pair (si, sj) with sj in the new sequences
            for (int si = 0; si < n; si++)
            {
                for (int sj = Math.Max(si + 1, p1); sj < n; sj++)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    dist[si + 1, sj + 1] = dist[sj + 1, si + 1] = pairwise(_seq[si], _seq[sj]);
                }
            }

            // getWeightsFromDistMat: the NJ (or UPGMA) tree of all rows, tree weights, similarities from the tree
            var dist0 = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    dist0[i, j] = dist[i + 1, j + 1];
            }

            var sim = new double[n + 1, n + 1];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    sim[i + 1, j + 1] = dist0[i, j];
            }

            var seqWeight = new int[n + 1];
            BuildWeightsAndSteps(dist0, sim, seqWeight);
            _seqWeight = seqWeight;
            _alnWeight = Enumerable.Repeat(100, n + 1).ToArray();
            _distMat = sim;
            SeqsAlignToProfile(sim, p1 - 1);
            return _seq.Select(CodesToRow).ToArray();
        }

        // MSA::seqsAlignToProfile: the first iStart + 1 sequences are the aligned profile; the others are added closest first with weight = tree weight x identity (no -NOWEIGHTS guard)
        private void SeqsAlignToProfile(double[,] distMat, int iStart)
        {
            int n = _n;
            var treeWeight = new int[n];
            for (int i = 0; i < n; i++)
                treeWeight[i] = _seqWeight[i];

            var maxid = new int[n + 1];
            for (int i = 1; i <= n; i++)
            {
                maxid[i] = -1;
                for (int j = 1; j <= n; j++)
                {
                    if (maxid[i] < distMat[i, j])
                        maxid[i] = (int)distMat[i, j];
                }
            }

            var aligned = new int[n + 1];
            int ix = 0;
            for (int i = 1; i <= iStart + 1; i++)
            {
                aligned[i] = 1;
                ix++;
            }

            var group = new int[n + 1];
            while (ix < n)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (ix > 0)
                {
                    for (int i = 1; i <= n; i++)
                    {
                        if (aligned[i] == 0)
                        {
                            maxid[i] = -1;
                            for (int j = 1; j <= n; j++)
                            {
                                if (maxid[i] < distMat[i, j] && aligned[j] != 0)
                                    maxid[i] = (int)distMat[i, j];
                            }
                        }
                    }
                }

                int max = -1, iseq = 0;
                for (int i = 1; i <= n; i++)
                {
                    if (aligned[i] == 0 && maxid[i] > max)
                    {
                        max = maxid[i];
                        iseq = i;
                    }
                }

                for (int j = 1; j <= n; j++)
                {
                    if (aligned[j] != 0)
                        group[j] = 1;
                    else if (iseq == j)
                        group[j] = 2;
                }

                aligned[iseq] = 1;
                for (int j = 0; j < n; j++)
                    _seqWeight[j] = (int)(treeWeight[j] * distMat[j + 1, iseq]);
                int sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (group[j + 1] == 1)
                        sum += _seqWeight[j];
                }

                if (sum == 0)
                {
                    for (int j = 0; j < n; j++)
                        _seqWeight[j] = 1;
                    sum = n;
                }

                for (int j = 0; j < n; j++)
                {
                    _seqWeight[j] = _seqWeight[j] * CwIntScaleFactor / sum;
                    if (_seqWeight[j] < 1)
                        _seqWeight[j] = 1;
                }

                Array.Copy(_seqWeight, _alnWeight, n + 1);
                ProfileAlign(group, aligned);
                ix++;
            }
        }

        // dscore = static_cast<int>(alignmentObj.countid(i, j)); distMat(i, j) = (100.0 - dscore) / 100.0
        private static void FillProfileDistances(CwAlignment aln, double[,] distMat, int first, int last)
        {
            for (int i = first; i <= last; i++)
            {
                for (int j = i + 1; j <= last; j++)
                {
                    int dscore = (int)aln.CountId(i, j);
                    distMat[i, j] = distMat[j, i] = (100.0 - dscore) / 100.0;
                }
            }
        }
    }
}
