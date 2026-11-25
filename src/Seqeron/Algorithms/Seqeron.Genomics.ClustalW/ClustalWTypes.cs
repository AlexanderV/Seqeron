// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.ClustalW;

using System;
using Seqeron.Genomics.Infrastructure;

/// <summary>
/// Source of the pairwise distances that build the guide tree of the progressive, iterative and consistency
/// multiple aligners.
/// </summary>
public enum MsaDistanceMode
{
    /// <summary>
    /// Default (historical behaviour): d = 1 − identities / min(len1, len2) of the repository's own global
    /// alignment under the caller's <see cref="ScoringMatrix"/> (linear-gap Needleman–Wunsch; the affine Gotoh alignment when
    /// <see cref="MsaOptions.GapModel"/> is <see cref="MsaGapModel.Affine"/>).
    /// </summary>
    Linear = 0,

    /// <summary>
    /// ClustalW 2.1 "slow" pairwise stage (Thompson, Higgins &amp; Gibson 1994): affine Smith–Waterman–Gotoh/Myers–Miller
    /// alignment with ClustalW's own DNA matrix (IUB, ×100) and pairwise gap penalties (open 15, extend 6.66), d = 1 −
    /// identities / min(len1, len2) of that alignment (<c>SequenceAligner.ClustalWPairwiseAlign</c>). The caller's scoring matrix is
    /// then used only for the profile alignment, not for the distances.
    /// </summary>
    ClustalWAffine = 1,

    /// <summary>
    /// ClustalW 2.1 "fast" pairwise stage (<c>-QUICKTREE</c>, <c>FastPairwiseAlign.cpp</c>; Wilbur &amp; Lipman 1983 k-tuple diagonal
    /// scoring): d = (100 − 100·fragmentScore / min(len1, len2)) / 100 in double precision, with the k-tuple / window / top-diagonals /
    /// gap parameters of <see cref="MsaOptions.QuickPairwise"/> (ClustalW's DNA defaults 2 / 4 / 4 / 5). The fragment score is not an
    /// alignment identity, so the distances (and therefore the guide tree) differ from <see cref="ClustalWAffine"/>
    /// (<c>SequenceAligner.ClustalWQuickPairwiseAlign</c>). The caller's scoring matrix is used only for the profile alignment.
    /// </summary>
    ClustalWQuick = 2,
}

/// <summary>Tree-building method of the guide tree of the progressive, iterative and consistency multiple aligners.</summary>
public enum MsaGuideTree
{
    /// <summary>
    /// Default (historical behaviour): UPGMA (Sokal &amp; Michener 1958; the Feng–Doolittle 1987 / CLUSTAL-1 choice) from
    /// <c>PhylogeneticAnalyzer</c>; ties between distances equal within 1e-9 go to the lowest (i, j) index pair.
    /// </summary>
    Upgma = 0,

    /// <summary>
    /// Neighbour-joining (Saitou &amp; Nei 1987) as ClustalW 2.1 uses it by default (<c>-CLUSTERING=NJ</c>; also the default
    /// tree of T-Coffee / Clustal Omega-style pipelines): the unrooted NJ tree is built with ClustalW's own neighbour-joining
    /// arithmetic (<c>NJTree::generateTree</c>; same topology as Saitou–Nei, exact ties resolved as clustalw does), re-rooted
    /// with ClustalW's <c>Tree::reRoot</c> rule and its merges are taken in ClustalW's order (left subtree first). With
    /// <see cref="MsaDistanceMode.ClustalWAffine"/> the distances are ClustalW's single-precision values, so the guide tree
    /// equals the one <c>clustalw -ALIGN</c> writes to its <c>.dnd</c> file.
    /// </summary>
    NeighborJoining = 1,
}

/// <summary>
/// Parameters of ClustalW's fast (k-tuple) pairwise stage for DNA: <c>-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE</c>
/// (<c>FastPairwiseAlign.cpp</c>, ClustalW 2.1); defaults are the DNA defaults of <c>UserParameters.cpp</c>.
/// </summary>
public sealed record ClustalWQuickOptions
{
    /// <summary>Word size (<c>-KTUPLE</c>), 1–4 for DNA; default 2.</summary>
    public int KTuple { get; init; } = 2;

    /// <summary>Number of diagonals scanned on each side of a top diagonal (<c>-WINDOW</c>); default 4.</summary>
    public int Window { get; init; } = 4;

    /// <summary>Number of best diagonals (<c>-TOPDIAGS</c>); default 4.</summary>
    public int TopDiagonals { get; init; } = 4;

    /// <summary>
    /// Gap penalty of the fragment chain (<c>-PAIRGAP</c>); default 5 (ClustalW's DNA default). Note: ClustalW's command line sets
    /// the pair gap to <see cref="KTuple"/> + 4 whenever <c>-KTUPLE</c> is given (and a later <c>-PAIRGAP</c> overrides it), so
    /// <c>-KTUPLE=2</c> alone means <c>PairGap = 6</c>, not 5; callers that want the command-line behaviour pass that value.
    /// </summary>
    public int PairGap { get; init; } = 5;

    /// <summary><c>-SCORE=PERCENT</c> (default): score = 100·fragmentScore / min(len1, len2); <c>false</c> = <c>-SCORE=ABSOLUTE</c> (the raw fragment score).</summary>
    public bool Percent { get; init; } = true;
}

/// <summary>ClustalW's guide-tree clustering algorithm for the multiple alignment (<c>-CLUSTERING=NJ | UPGMA</c>).</summary>
public enum ClustalWClustering
{
    /// <summary>Neighbour joining (default).</summary>
    NeighborJoining = 0,

    /// <summary>
    /// UPGMA (<c>-CLUSTERING=UPGMA</c>, <c>UPGMAAlgorithm</c>): ClustalW's own rooted average-linkage tree; its node heights give the sequence weights and its joins the order of the profile steps.
    /// With two sequences the tree step stores no alignment step, but the second pass of the program's progressive stage aligns the pair, as for neighbour joining.
    /// Verified against the program (built at -O0; the installed binary crashes on UPGMA) on 600 sets, F101.
    /// </summary>
    Upgma = 1,
}

/// <summary>Residue type of ClustalW's pairwise stages (<c>-TYPE=DNA | PROTEIN</c>).</summary>
public enum ClustalWSequenceType
{
    /// <summary>Nucleotides (IUB matrix, gap 15 / 6.66, k-tuple 2).</summary>
    Dna = 0,

    /// <summary>Amino acids (GONNET / BLOSUM30 / PAM350 / identity matrix, gap 10 / 0.1, k-tuple 1).</summary>
    Protein = 1,
}

/// <summary>The pairwise substitution matrix of ClustalW's PROTEIN full pairwise stage (<c>-PWMATRIX</c>).</summary>
public enum ClustalWPairwiseMatrix
{
    /// <summary><c>GONNET</c> (default): <c>gon250mt</c> × intScale 10 (<see cref="ClustalWMatrices.ClustalWGonnet250"/>).</summary>
    Gonnet = 0,

    /// <summary><c>BLOSUM</c>: <c>blosum30mt</c> × 100 (<see cref="ClustalWMatrices.ClustalWBlosum30"/>).</summary>
    Blosum = 1,

    /// <summary><c>PAM</c>: <c>pam350mt</c> × 100 (<see cref="ClustalWMatrices.ClustalWPam350"/>).</summary>
    Pam = 2,

    /// <summary><c>ID</c>: identity matrix (<see cref="ClustalWMatrices.ClustalWIdentity"/>) × 100.</summary>
    Identity = 3,
}

/// <summary>
/// Parameters of ClustalW 2.1's PROTEIN pairwise stages (<c>clustalw -TYPE=PROTEIN</c>): the full dynamic-programming stage (<c>-PWMATRIX</c>,
/// <c>-PWGAPOPEN</c>, <c>-PWGAPEXT</c>) and the fast k-tuple stage (<c>-QUICKTREE</c>: <c>-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE</c>); defaults are
/// the PROTEIN defaults of <c>UserParameters.cpp</c>.
/// </summary>
public sealed record ClustalWProteinPairwiseOptions
{
    /// <summary>Pairwise matrix (<c>-PWMATRIX</c>, default <see cref="ClustalWPairwiseMatrix.Gonnet"/>).</summary>
    public ClustalWPairwiseMatrix Matrix { get; init; } = ClustalWPairwiseMatrix.Gonnet;

    /// <summary>Pairwise gap-opening penalty (<c>-PWGAPOPEN</c>, PROTEIN default 10.0; a positive penalty, stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 10.0;

    /// <summary>Pairwise gap-extension penalty (<c>-PWGAPEXT</c>, PROTEIN default 0.1).</summary>
    public double GapExtend { get; init; } = 0.1;

    /// <summary>
    /// Parameters of the fast stage; PROTEIN defaults <c>KTuple</c> 1 (1-2), <c>Window</c> 5, <c>TopDiagonals</c> 5, <c>PairGap</c> 3 (ClustalW's command line sets the pair
    /// gap to <c>ktuple + 3</c> when <c>-KTUPLE</c> is given).
    /// </summary>
    public ClustalWQuickOptions Quick { get; init; } = new() { KTuple = 1, Window = 5, TopDiagonals = 5, PairGap = 3 };

    /// <summary>
    /// <c>-PWMATRIX=file</c>: a user-defined matrix (<see cref="ClustalWUserMatrix.Parse"/> / <see cref="ClustalWUserMatrix.FromSubstitutionMatrix"/>) in place of
    /// <see cref="Matrix"/> (read as ClustalW does: x100 integers, never shifted to non-negative scores; the gap opening is scaled by the matrix average score).
    /// </summary>
    public ClustalWUserMatrix? UserMatrix { get; init; }
}

/// <summary>
/// Result of ClustalW's fast k-tuple pairwise stage (<c>SequenceAligner.ClustalWQuickPairwiseAlign</c>, port of
/// <c>FastPairwiseAlign::pairAlign</c>, ClustalW 2.1).
/// </summary>
/// <param name="FragmentScore">Score of the best k-tuple fragment chain (<c>accum[0][maxSoFar]</c>; 0 when no fragment is found).</param>
/// <param name="Score">The number ClustalW prints as <c>Score:</c> (<c>%lg</c>): the percent 100·FragmentScore / min(len1, len2), or the raw fragment score with <c>-SCORE=ABSOLUTE</c>.</param>
/// <param name="Distance">(100 − Score) / 100 (double precision): ClustalW's guide-tree distance (1 when not aligned; with the absolute score it can be negative, as in ClustalW).</param>
/// <param name="Aligned"><c>Score &gt; 0.1</c>: ClustalW prints <c>Aligned. Score: …</c>, otherwise <c>Not Aligned</c>.</param>
/// <param name="PartialAlignment">The fragment list exceeded ClustalW's bound (2 · 2 · longest sequence) and the scan stopped early (<c>(Partial alignment)</c>).</param>
public sealed record ClustalWQuickPairwiseResult(
    int FragmentScore,
    double Score,
    double Distance,
    bool Aligned,
    bool PartialAlignment);

/// <summary>ClustalW's iterative refinement of the multiple alignment (<c>-ITERATION=</c>; <c>Iteration.cpp</c>).</summary>
public enum ClustalWIteration
{
    /// <summary><c>none</c> (default): no iteration.</summary>
    None = 0,

    /// <summary><c>tree</c>: after every profile-profile step of the progressive alignment, each sequence of the step is removed and re-aligned to the others (<c>Iteration::iterationOnTreeNode</c>).</summary>
    Tree = 1,

    /// <summary><c>alignment</c>: the same "remove first" re-alignment once on the finished alignment (<c>Iteration::removeFirstIterate</c>).</summary>
    Alignment = 2,
}

/// <summary>Pairwise stage that produces the distance matrix of the ClustalW multiple alignment (<see cref="ClustalWOptions.PairwiseMode"/>).</summary>
public enum ClustalWPairwiseMode
{
    /// <summary>ClustalW's default "slow" full dynamic-programming stage (<c>FullPairwiseAlign</c>; <c>SequenceAligner.ClustalWPairwiseAlign</c>).</summary>
    Full = 0,

    /// <summary>ClustalW's fast k-tuple stage (<c>-QUICKTREE</c>, <c>FastPairwiseAlign</c>; <c>SequenceAligner.ClustalWQuickPairwiseAlign</c>).</summary>
    Quick = 1,
}

/// <summary>
/// Parameters of <c>SequenceAligner.MultipleAlignClustalW</c>, the port of ClustalW 2.1's multiple-alignment stage for DNA
/// (<c>clustalw -TYPE=DNA -ALIGN</c>). Defaults are ClustalW's DNA defaults.
/// </summary>
public sealed record ClustalWOptions
{
    /// <summary>
    /// Cancels the alignment (the same convention as <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>): the dynamic programs poll it at their
    /// row / iteration boundaries and the call throws <see cref="OperationCanceledException"/> without a partial result. Default: never cancelled.
    /// It does not take part in the semantics of the options (two option sets that differ only in the token align identically).
    /// </summary>
    public CancellationToken CancellationToken { get; init; }

    /// <summary>Gap-opening penalty of the profile alignment (<c>-GAPOPEN</c>, DNA default 15.0; stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 15.0;

    /// <summary>Gap-extension penalty of the profile alignment (<c>-GAPEXT</c>, DNA default 6.66).</summary>
    public double GapExtend { get; init; } = 6.66;

    /// <summary>Pairwise stage that produces the distances of the guide tree (default <see cref="ClustalWPairwiseMode.Full"/>).</summary>
    public ClustalWPairwiseMode PairwiseMode { get; init; } = ClustalWPairwiseMode.Full;

    /// <summary>Pairwise gap-opening penalty of the full pairwise stage (<c>-PWGAPOPEN</c>, default 15).</summary>
    public double PairwiseGapOpen { get; init; } = 15.0;

    /// <summary>Pairwise gap-extension penalty of the full pairwise stage (<c>-PWGAPEXT</c>, default 6.66).</summary>
    public double PairwiseGapExtend { get; init; } = 6.66;

    /// <summary>Parameters of the fast pairwise stage (<see cref="ClustalWPairwiseMode.Quick"/> only).</summary>
    public ClustalWQuickOptions Quick { get; init; } = new();

    /// <summary>
    /// Weight of transitions (A–G, C–T/U) in the profile substitution matrix (<c>-TRANSWEIGHT</c>, default 0.5; 0 = transitions
    /// score as mismatches, 1 = as matches).
    /// </summary>
    public float TransitionWeight { get; init; } = 0.5f;

    /// <summary>
    /// <c>-MAXDIV</c> (default 30): a sequence whose highest percent identity to any other sequence (as derived from the guide
    /// tree) is at most this value is delayed and added after the rest has been aligned (delayed divergent sequences).
    /// </summary>
    public int DivergenceCutoff { get; init; } = 30;

    /// <summary><c>-CLUSTERING=NJ | UPGMA</c> of the multiple alignment (default neighbour joining); see <see cref="ClustalWClustering"/>.</summary>
    public ClustalWClustering Clustering { get; init; } = ClustalWClustering.NeighborJoining;

    /// <summary>
    /// <c>-ITERATION=TREE | ALIGNMENT | NONE</c> (default none): iterative refinement of the alignment (<see cref="ClustalWIteration"/>). Note: the program only recognises the lower-case
    /// words <c>tree</c> / <c>alignment</c> / <c>none</c> on its command line (upper case prints "Unknown option for iteration. Setting to NONE").
    /// </summary>
    public ClustalWIteration Iteration { get; init; } = ClustalWIteration.None;

    /// <summary><c>-NUMITER=n</c> (default 3; the program refuses n &lt;= 0): the maximum number of iterations of <see cref="Iteration"/>. The program keeps the best re-alignment of the first iteration, so a larger value changes nothing once an iteration finds no improvement.</summary>
    public int NumIterations { get; init; } = 3;

    /// <summary><c>-NOWEIGHTS</c>: no sequence weights from the guide tree (all sequences weigh the same).</summary>
    public bool NoWeights { get; init; }

    /// <summary><c>-NEGATIVE</c>: use the substitution matrix as it is (negative mismatch scores) instead of shifting it to non-negative scores.</summary>
    public bool NegativeMatrix { get; init; }

    /// <summary>
    /// <c>UserParameters::useEndGaps</c> (the menu item "Toggle End Gap Separation"; the command-line flag <c>-ENDGAPS</c> only ever writes <c>false</c>): gaps at the ends
    /// of the sequences count as gaps in the position-specific penalties of <c>ProfileBase::calcGapCoeff</c>. Default <c>false</c>.
    /// </summary>
    public bool UseEndGaps { get; init; }

    /// <summary>
    /// <c>UserParameters::endGapPenalties</c> (no command-line switch; the program's own code sets it only for profile alignments): the end gaps of a profile alignment are
    /// penalised (gap opening / extension at the profile ends instead of free). Default <c>false</c>. With it a profile alignment of the first pass can score below zero,
    /// where ClustalW 2.1 gives up silently: see <see cref="ClustalWAbortedException"/>.
    /// </summary>
    public bool EndGapPenalties { get; init; }

    /// <summary>
    /// <c>-DNAMATRIX=IUB | CLUSTALW</c>: the profile-stage matrix: <see cref="ClustalWDnaMatrix.Iub"/> (default) or <see cref="ClustalWDnaMatrix.Clustalw"/> (ClustalW 1.6's
    /// identity matrix, gap penalties scaled by 0.66). Ignored when <see cref="UserDnaMatrix"/> is set.
    /// </summary>
    public ClustalWDnaMatrix DnaMatrix { get; init; } = ClustalWDnaMatrix.Iub;

    /// <summary><c>-PWDNAMATRIX=IUB | CLUSTALW</c>: the matrix of the full pairwise stage (see <see cref="ClustalWDnaPairwiseOptions.Matrix"/>); ignored when <see cref="PairwiseUserDnaMatrix"/> is set.</summary>
    public ClustalWDnaMatrix PairwiseDnaMatrix { get; init; } = ClustalWDnaMatrix.Iub;

    /// <summary>
    /// <c>-DNAMATRIX=file</c> (command line) / "User defined" of the DNA weight-matrix menu: a user-defined nucleotide matrix (<see cref="ClustalWUserMatrix.Parse"/>) for the profile
    /// stage. How it is used is decided by <see cref="UserDnaMatrixUse"/>.
    /// </summary>
    public ClustalWUserMatrix? UserDnaMatrix { get; init; }

    /// <summary>The user-defined nucleotide matrix of the full pairwise stage (<c>-PWDNAMATRIX=file</c> / the pairwise-parameters menu), see <see cref="UserDnaMatrixUse"/>.</summary>
    public ClustalWUserMatrix? PairwiseUserDnaMatrix { get; init; }

    /// <summary>What ClustalW does with <see cref="UserDnaMatrix"/> / <see cref="PairwiseUserDnaMatrix"/> (default <see cref="ClustalWUserDnaMatrixUse.Applied"/>).</summary>
    public ClustalWUserDnaMatrixUse UserDnaMatrixUse { get; init; } = ClustalWUserDnaMatrixUse.Applied;
}

/// <summary>
/// Parameters of ClustalW 2.1's full pairwise stage for DNA (<c>clustalw -TYPE=DNA</c>): <c>-PWGAPOPEN</c>, <c>-PWGAPEXT</c>, <c>-TRANSWEIGHT</c> and <c>-PWDNAMATRIX</c>.
/// Defaults are ClustalW's DNA defaults.
/// </summary>
public sealed record ClustalWDnaPairwiseOptions
{
    /// <summary>Pairwise gap-opening penalty (<c>-PWGAPOPEN</c>, DNA default 15.0; stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 15.0;

    /// <summary>Pairwise gap-extension penalty (<c>-PWGAPEXT</c>, DNA default 6.66).</summary>
    public double GapExtend { get; init; } = 6.66;

    /// <summary><c>-TRANSWEIGHT</c> (default 0.5): the score of the pairs written to the residue codes (0,4), (2,11), (2,12) is <c>int(TransitionWeight * matrix[0][0])</c>.</summary>
    public float TransitionWeight { get; init; } = 0.5f;

    /// <summary><c>-PWDNAMATRIX=IUB | CLUSTALW</c> (default IUB); CLUSTALW scales the gap opening by 0.6667 and the extension by 0.751. Ignored when <see cref="UserMatrix"/> is applied.</summary>
    public ClustalWDnaMatrix Matrix { get; init; } = ClustalWDnaMatrix.Iub;

    /// <summary><c>-PWDNAMATRIX=file</c> (see <see cref="UserMatrixUse"/>): a user-defined nucleotide matrix (x100, never shifted; <see cref="ClustalWUserMatrix.Parse"/>).</summary>
    public ClustalWUserMatrix? UserMatrix { get; init; }

    /// <summary>What ClustalW does with <see cref="UserMatrix"/> (default <see cref="ClustalWUserDnaMatrixUse.Applied"/>).</summary>
    public ClustalWUserDnaMatrixUse UserMatrixUse { get; init; } = ClustalWUserDnaMatrixUse.Applied;
}

/// <summary>The built-in nucleotide matrices of ClustalW 2.1 (<c>-DNAMATRIX</c>, <c>-PWDNAMATRIX</c>).</summary>
public enum ClustalWDnaMatrix
{
    /// <summary><c>IUB</c> (default): <c>swgapdnamt</c> - +10 for a match (IUB ambiguity codes: any overlap), -9 for a mismatch (x100 in the profile stage).</summary>
    Iub = 0,

    /// <summary><c>CLUSTALW</c> ("CLUSTALW(1.6)"): <c>clustalvdnamt</c> - +10 on the diagonal of A, C, G, T and U, 0 everywhere else; the profile stage scales the gap
    /// penalties by 0.66, the pairwise stage by 0.6667 (opening) / 0.751 (extension).</summary>
    Clustalw = 1,
}

/// <summary>
/// How ClustalW 2.1 treats a user-defined nucleotide matrix. The two behaviours differ because of a defect of the program: on the command line
/// (<c>-DNAMATRIX=file</c>, <c>-PWDNAMATRIX=file</c>) <c>CommandLineParser::interpretParameters</c> reads and validates the file but then registers its name with
/// <c>setCurrentNameAndNum(name, 3, Protein, ...)</c> - the PROTEIN slot - so the DNA stages keep the IUB matrix; the interactive menu ("DNA weight matrix",
/// "User defined") registers it for DNA and the matrix is used.
/// </summary>
public enum ClustalWUserDnaMatrixUse
{
    /// <summary>The matrix is used (the interactive-menu behaviour: <c>userDNAMat</c> / <c>DNAXref</c>, <c>pwUserDNAMat</c> / <c>pwDNAXref</c>).</summary>
    Applied = 0,

    /// <summary><c>clustalw -TYPE=DNA -DNAMATRIX=file -PWDNAMATRIX=file</c> exactly: the file is read and validated (a malformed file ends the program) but the stages keep the IUB
    /// matrix (or the built-in one chosen with <see cref="ClustalWOptions.DnaMatrix"/>).</summary>
    CommandLine = 1,
}

/// <summary>Raised when ClustalW 2.1 would abandon the multiple alignment silently.</summary>
/// <remarks>
/// <c>MSA::multiSeqAlign</c> returns -1 when the score of a profile alignment of its FIRST pass (the guide-tree steps among the sequences that are not delayed) is negative
/// (<c>score &lt; 0</c>; the score is the Myers-Miller optimum divided by 100, truncated), and <c>Clustal::align</c> then returns without writing any alignment (exit
/// status 0, no output file). A profile alignment cannot score below zero while the end gaps are free (the alignment of two profiles consisting only of end gaps costs nothing);
/// with <see cref="ClustalWOptions.EndGapPenalties"/> every gap is charged and a poor pair of profiles can fall below zero. This library reports the situation with
/// this exception instead of returning an alignment ClustalW itself never produces.
/// </remarks>
public sealed class ClustalWAbortedException : InvalidOperationException
{
    /// <summary>Creates the exception.</summary>
    /// <param name="group">1-based number of the profile alignment ("Group n") whose score was negative.</param>
    /// <param name="score">The negative score (optimum / 100, truncated).</param>
    public ClustalWAbortedException(int group, int score)
        : base($"ClustalW 2.1 abandons the alignment: the profile alignment of group {group} scores {score} (< 0) with end-gap penalties; the program writes no output.")
    {
        Group = group;
        Score = score;
    }

    /// <summary>1-based number of the first-pass profile alignment ("Group n") with the negative score.</summary>
    public int Group { get; }

    /// <summary>The negative score that ended the run.</summary>
    public int Score { get; }
}

/// <summary>The substitution-matrix series of ClustalW's PROTEIN profile alignment (<c>-MATRIX</c>; one table of the series is chosen per profile pair by percent identity).</summary>
public enum ClustalWProfileMatrix
{
    /// <summary><c>GONNET</c> (default): GONNET80 / 250 / 120 / 350 / 160 (<see cref="ClustalWMatrices.ClustalWGonnet80"/> ...), intScale 10.</summary>
    Gonnet = 0,

    /// <summary><c>BLOSUM</c>: BLOSUM80 / 62 / 45 / 30 (<see cref="ClustalWMatrices.ClustalWBlosum80"/> ...), 40 with <c>-NEGATIVE</c>.</summary>
    Blosum = 1,

    /// <summary><c>PAM</c>: PAM20 / 60 / 120 / 350 (<see cref="ClustalWMatrices.ClustalWPam20"/> ...).</summary>
    Pam = 2,

    /// <summary><c>ID</c>: the identity matrix (<see cref="ClustalWMatrices.ClustalWIdentity"/>).</summary>
    Identity = 3,
}

/// <summary>
/// Parameters of <c>SequenceAligner.MultipleAlignClustalW</c> for PROTEIN (the port of ClustalW 2.1's multiple-alignment stage for amino acids,
/// <c>clustalw -TYPE=PROTEIN -ALIGN</c>). Defaults are ClustalW's PROTEIN defaults.
/// </summary>
public sealed record ClustalWProteinOptions
{
    /// <summary>
    /// Cancels the alignment (the same convention as <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>): the dynamic programs poll it at their
    /// row / iteration boundaries and the call throws <see cref="OperationCanceledException"/> without a partial result. Default: never cancelled.
    /// It does not take part in the semantics of the options (two option sets that differ only in the token align identically).
    /// </summary>
    public CancellationToken CancellationToken { get; init; }

    /// <summary>Matrix series of the profile alignment (<c>-MATRIX</c>, default <see cref="ClustalWProfileMatrix.Gonnet"/>).</summary>
    public ClustalWProfileMatrix Matrix { get; init; } = ClustalWProfileMatrix.Gonnet;

    /// <summary>Gap-opening penalty of the profile alignment (<c>-GAPOPEN</c>, PROTEIN default 10.0; stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 10.0;

    /// <summary>Gap-extension penalty of the profile alignment (<c>-GAPEXT</c>, PROTEIN default 0.2).</summary>
    public double GapExtend { get; init; } = 0.2;

    /// <summary>Pairwise stage that produces the distances of the guide tree (default <see cref="ClustalWPairwiseMode.Full"/>; <c>-QUICKTREE</c> = <see cref="ClustalWPairwiseMode.Quick"/>).</summary>
    public ClustalWPairwiseMode PairwiseMode { get; init; } = ClustalWPairwiseMode.Full;

    /// <summary>Pairwise stage parameters: <c>-PWMATRIX</c>, <c>-PWGAPOPEN</c>, <c>-PWGAPEXT</c> and the fast-stage <c>-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP</c>.</summary>
    public ClustalWProteinPairwiseOptions Pairwise { get; init; } = new();

    /// <summary><c>-MAXDIV</c> (default 30): see <see cref="ClustalWOptions.DivergenceCutoff"/>.</summary>
    public int DivergenceCutoff { get; init; } = 30;

    /// <summary><c>-CLUSTERING=NJ | UPGMA</c> of the multiple alignment (default neighbour joining); see <see cref="ClustalWClustering"/>.</summary>
    public ClustalWClustering Clustering { get; init; } = ClustalWClustering.NeighborJoining;

    /// <summary>
    /// <c>-ITERATION=TREE | ALIGNMENT | NONE</c> (default none): iterative refinement of the alignment (<see cref="ClustalWIteration"/>). Note: the program only recognises the lower-case
    /// words <c>tree</c> / <c>alignment</c> / <c>none</c> on its command line (upper case prints "Unknown option for iteration. Setting to NONE").
    /// </summary>
    public ClustalWIteration Iteration { get; init; } = ClustalWIteration.None;

    /// <summary><c>-NUMITER=n</c> (default 3; the program refuses n &lt;= 0): the maximum number of iterations of <see cref="Iteration"/>. The program keeps the best re-alignment of the first iteration, so a larger value changes nothing once an iteration finds no improvement.</summary>
    public int NumIterations { get; init; } = 3;

    /// <summary><c>-NOWEIGHTS</c>: no sequence weights from the guide tree.</summary>
    public bool NoWeights { get; init; }

    /// <summary><c>-NEGATIVE</c>: use the matrix as it is (negative scores allowed) and the single matrix BLOSUM40 / PAM120 / GONNET250 instead of the series; the gap-opening penalty is then not scaled by the matrix.</summary>
    public bool NegativeMatrix { get; init; }

    /// <summary><c>-NOPGAP</c>: switch off the residue-specific gap penalties (Pascarella &amp; Argos 1992; ClustalW's default is on).</summary>
    public bool NoPreferredGapPenalties { get; init; }

    /// <summary><c>-NOHGAP</c>: switch off the hydrophilic gap penalties (runs of more than three <see cref="HydrophilicResidues"/> halve the opening penalty; default on).</summary>
    public bool NoHydrophilicGapPenalties { get; init; }

    /// <summary>
    /// <c>-NOVGAP</c> (sic: ClustalW's command line option of this name switches the variable-region penalties ON; ClustalW's default is off): the opening penalty
    /// of a column is scaled by the local divergence of the first two sequences of a two-sequence profile with more than 60 % identity.
    /// </summary>
    public bool VariableRegionGapPenalties { get; init; }

    /// <summary><c>-GAPDIST</c> (default 4): columns within this distance of an existing gap get a higher opening penalty (2 to 4 times at the gap, decreasing with the distance); 0 = off.</summary>
    public int GapSeparationDistance { get; init; } = 4;

    /// <summary>
    /// <c>-HGAPRESIDUES</c> (default <c>GPSNDQEKR</c>): the residues that make up a hydrophilic stretch (letters, case-insensitive; ClustalW keeps at most the
    /// first 9 letters, <c>MAXHYDRESIDUES</c>).
    /// </summary>
    public string HydrophilicResidues { get; init; } = "GPSNDQEKR";

    /// <summary>
    /// The <c>gdist</c> values of the profile objects for which ClustalW never assigns it: with <see cref="VariableRegionGapPenalties"/> and a two-sequence profile of
    /// more than 60 % identity <c>ProfileBase::calcGapCoeff</c> uses the member <c>gdist</c> without initialising it (undefined behaviour of the C++: the member holds
    /// whatever the recycled heap chunk contained, in practice 0 or the <see cref="GapSeparationDistance"/> of an earlier profile, rarely a stale residue code).
    /// <c>null</c> (default) emulates the common outcome (0 for the first profile pair, afterwards <see cref="GapSeparationDistance"/>); a list gives the value of the
    /// 1st, 2nd, ... such profile of the run (later ones use the emulation), which reproduces a specific <c>clustalw</c> run.
    /// </summary>
    public IReadOnlyList<int>? UninitialisedGapDistances { get; init; }

    /// <summary>
    /// <c>-MATRIX=file</c> with a single matrix: the profile stage uses this matrix for every pair of profiles (intScale 100, gap-penalty scale 0.75), in place of
    /// <see cref="Matrix"/>. Ignored when <see cref="UserMatrixSeries"/> is set.
    /// </summary>
    public ClustalWUserMatrix? UserMatrix { get; init; }

    /// <summary>
    /// <c>-MATRIX=file</c> with a <c>CLUSTAL_SERIES</c> file (<see cref="ClustalWMatrixFile.Parse"/>): the matrix of the first entry whose percent-identity range
    /// contains the identity of the profile pair (the first entry if none does), with the gap-penalty scale <c>0.5 + (pcid - llimit) / ((ulimit - llimit) * 2)</c>.
    /// </summary>
    public IReadOnlyList<ClustalWMatrixSeriesEntry>? UserMatrixSeries { get; init; }

    /// <summary><c>UserParameters::useEndGaps</c> (menu item "Toggle End Gap Separation"; <c>-ENDGAPS</c> on the command line only writes <c>false</c>): see <see cref="ClustalWOptions.UseEndGaps"/>.</summary>
    public bool UseEndGaps { get; init; }

    /// <summary><c>UserParameters::endGapPenalties</c>: the end gaps of a profile alignment are penalised; see <see cref="ClustalWOptions.EndGapPenalties"/> and <see cref="ClustalWAbortedException"/>.</summary>
    public bool EndGapPenalties { get; init; }
}

/// <summary>
/// Result of ClustalW's "slow" pairwise alignment (<c>SequenceAligner.ClustalWPairwiseAlign</c>, port of
/// <c>FullPairwiseAlign::pairwiseAlign</c>, ClustalW 2.1): a local affine alignment whose identities define ClustalW's percent
/// identity and distance.
/// </summary>
/// <param name="Score">Alignment score in ClustalW's integer units (matrix ×100, gap penalties ×100).</param>
/// <param name="Identities">Aligned residue pairs with identical residue codes in the local alignment.</param>
/// <param name="PercentIdentity">100 · identities / min(len1, len2) (0 when a sequence is empty).</param>
/// <param name="Distance">1 − identities / min(len1, len2) (1 when a sequence is empty): ClustalW's guide-tree distance.</param>
/// <param name="ClustalWScore">The integer "Score:" ClustalW prints per pair: (int)(100·identities / min(len1, len2)) in single precision.</param>
/// <param name="AlignedSequence1">Aligned region of sequence 1 (gap '-').</param>
/// <param name="AlignedSequence2">Aligned region of sequence 2.</param>
/// <param name="Start1">0-based first aligned residue of sequence 1 (-1 when the alignment is empty).</param>
/// <param name="End1">0-based last aligned residue of sequence 1 (-1 when empty).</param>
/// <param name="Start2">0-based first aligned residue of sequence 2 (-1 when empty).</param>
/// <param name="End2">0-based last aligned residue of sequence 2 (-1 when empty).</param>
public sealed record ClustalWPairwiseResult(
    int Score,
    int Identities,
    double PercentIdentity,
    double Distance,
    int ClustalWScore,
    string AlignedSequence1,
    string AlignedSequence2,
    int Start1,
    int End1,
    int Start2,
    int End2);

/// <summary>The order in which ClustalW lists the rows of an alignment (<c>-OUTORDER</c>).</summary>
public enum ClustalWOutputOrder
{
    /// <summary><c>ALIGNED</c> (the program's default): the order in which the sequences entered the progressive alignment (<c>MSA::multiSeqAlign</c> / <c>seqsAlignToProfile</c>; the profile commands list profile 1 then profile 2).</summary>
    Aligned = 0,

    /// <summary><c>INPUT</c>: the order of the input file.</summary>
    Input = 1,
}

/// <summary>Why ClustalW refuses a guide-tree file (the program prints "ERROR: ..." and ends the run without an alignment).</summary>
public enum ClustalWTreeFileProblem
{
    /// <summary>The first character of the tree is not '(' ("Wrong format in tree file").</summary>
    WrongFormat = 0,

    /// <summary>The tree has a different number of leaves than the alignment has sequences ("tree not compatible with alignment").</summary>
    NotCompatible = 1,

    /// <summary>A sequence name is not a leaf of the tree ("tree not compatible with alignment: name not found"; names are compared case-insensitively, blanks as '_').</summary>
    NameNotFound = 2,

    /// <summary>The tree is unrooted (three branches at the top) and has no branch lengths ("input tree is unrooted and has no distances. Cannot align sequences").</summary>
    UnrootedWithoutDistances = 3,

    /// <summary>
    /// The tree drives the program into undefined behaviour that always ends in a hang or a crash on the oracle sets: the text ends inside a sequence name (the reader then reads past the end of the file for ever),
    /// the tree has more nodes or leaves than the reader's arrays hold (heap overflow), or two branches share a node (nested trichotomy: the tree is freed twice, or has a cycle).
    /// This library refuses these trees instead.
    /// </summary>
    Undefined = 4,
}

/// <summary>Raised when ClustalW 2.1 refuses a guide tree given with <c>-USETREE</c> / <c>-USETREE1</c> / <c>-USETREE2</c> (<c>Tree::readTree</c> returns 0).</summary>
public sealed class ClustalWTreeFileException : InvalidOperationException
{
    /// <summary>Creates the exception.</summary>
    public ClustalWTreeFileException(ClustalWTreeFileProblem problem, string message)
        : base(message)
    {
        Problem = problem;
    }

    /// <summary>The program's reason for the refusal.</summary>
    public ClustalWTreeFileProblem Problem { get; }
}

/// <summary>A named input sequence of a ClustalW run (names are the first word of the FASTA header: no blanks, at most 150 characters, unique).</summary>
/// <param name="Name">Sequence name.</param>
/// <param name="Residues">Residues (gaps '-' only for profiles).</param>
public sealed record ClustalWSequence(string Name, string Residues);

/// <summary>
/// The command-level parameters of a ClustalW run that are independent of the alignment parameters: the guide-tree files and the row order
/// (<c>-NEWTREE</c>, <c>-USETREE</c>, <c>-NEWTREE1/2</c>, <c>-USETREE1/2</c>, <c>-OUTORDER</c>).
/// </summary>
public sealed record ClustalWRunOptions
{
    /// <summary><c>-OUTORDER=ALIGNED | INPUT</c> (default ALIGNED, the program's default).</summary>
    public ClustalWOutputOrder OutputOrder { get; init; } = ClustalWOutputOrder.Aligned;

    /// <summary>
    /// <c>-USETREE=file</c>: the text of a guide tree in Newick / PHYLIP format; the alignment follows it instead of a computed tree (<c>Clustal::doAlignUseOldTree</c>; for <c>-SEQUENCES</c> the tree of all
    /// sequences, <c>getWeightsFromGuideTree</c>). The program reads this tree with its neighbour-joining reader even for <c>-CLUSTERING=UPGMA</c>. Note that <c>-ALIGN -USETREE=x</c> on the command line
    /// does not use the file (it overwrites it with a new tree); this property is the option of <c>-USETREE</c> alone.
    /// </summary>
    public string? UseTree { get; init; }

    /// <summary><c>-USETREE1=file</c> of <c>-PROFILE</c>: the guide tree of profile 1 (ignored for a profile of one sequence, as in the program).</summary>
    public string? UseTree1 { get; init; }

    /// <summary><c>-USETREE2=file</c> of <c>-PROFILE</c>: the guide tree of profile 2.</summary>
    public string? UseTree2 { get; init; }

    /// <summary>The secondary-structure information of profile 1 (<c>!SS_</c> / <c>!GM_</c> lines of a CLUSTAL file, <c>"SS_</c> entries of a GDE file); null = none.</summary>
    public ClustalWProfileStructure? ProfileStructure1 { get; init; }

    /// <summary>The secondary-structure information of profile 2 (or of the sequences of <c>-SEQUENCES</c>).</summary>
    public ClustalWProfileStructure? ProfileStructure2 { get; init; }

    /// <summary>The alignment output parameters (<c>-OUTPUT</c> ... <c>-RANGE</c>, <c>-SECSTROUT</c>) used by <c>ClustalWAlignmentWriter</c> and the command line.</summary>
    public ClustalWOutputOptions Output { get; init; } = new();

    /// <summary>The secondary-structure gap-penalty parameters (<c>-HELIXGAP</c> ... <c>-NOSECSTR2</c>).</summary>
    public ClustalWStructureOptions Structure { get; init; } = new();

    /// <summary>
    /// <c>-NEWTREE</c> without <c>-ALIGN</c> / <c>-NEWTREE1</c> / <c>-NEWTREE2</c>: only the guide tree(s) are produced and no alignment is made (<c>Clustal::doGuideTreeOnly</c>; with the profile commands
    /// the program writes both profile trees and stops). When false the guide tree text of an alignment run is returned as well (<c>-ALIGN -NEWTREE=file</c>).
    /// </summary>
    public bool GuideTreeOnly { get; init; }
}

/// <summary>
/// The secondary-structure information of one profile as ClustalW reads it from the profile file (<c>ClustalFileParser::getSecStructure</c>): the structure line (<c>!SS_name</c>: helix 'A'/'a'/'$', strand 'B'/'b'/'%', anything
/// else loop) or the gap-penalty-mask line (<c>!GM_name</c>: digits 1-9); the program uses the FIRST of the two it finds. Characters are taken left to right over the alignment columns; blanks are skipped.
/// </summary>
/// <param name="Name">The structure name (printed at the beginning of the line in the alignment file).</param>
/// <param name="SecondaryStructure">The structure line, or null when the profile gives a mask only.</param>
/// <param name="GapPenaltyMask">The gap-penalty-mask line (used when <paramref name="SecondaryStructure"/> is null).</param>
public sealed record ClustalWProfileStructure(string Name, string? SecondaryStructure, string? GapPenaltyMask);

/// <summary>The structure information of one profile after the alignment (the masks have got the gaps of the alignment).</summary>
/// <param name="Name">Structure name.</param>
/// <param name="HasSecondaryStructure">True when the profile gave a structure line (the gap-penalty mask was computed from it), false when it gave a mask.</param>
/// <param name="SecondaryStructure">The aligned structure line ('-' for the gaps of the alignment); empty for a mask-only profile.</param>
/// <param name="GapPenaltyMask">The aligned gap-penalty mask.</param>
/// <param name="Used">False after <c>-NOSECSTR1</c> / <c>-NOSECSTR2</c> (the mask is not written).</param>
public sealed record ClustalWAlignmentStructure(string Name, bool HasSecondaryStructure, string SecondaryStructure, string GapPenaltyMask, bool Used);

/// <summary>The alignment ClustalW holds after a run: names, rows in input order, and what its output routines need.</summary>
/// <param name="Names">Sequence names (input order).</param>
/// <param name="Rows">Aligned rows (input order; profile commands: profile 1 then profile 2).</param>
/// <param name="OutputIndex">The (0-based) sequence indices in the order the program lists them (<see cref="ClustalWOutputOrder"/>).</param>
/// <param name="IsProtein">Residue type of the run (<c>-TYPE</c>).</param>
/// <param name="SequenceWeights">The sequence weights of the alignment object at the end of the run (<c>Alignment::getSeqWeight</c>; per mille, summing to about 1000).</param>
public sealed record ClustalWAlignment(
    IReadOnlyList<string> Names,
    IReadOnlyList<string> Rows,
    IReadOnlyList<int> OutputIndex,
    bool IsProtein,
    IReadOnlyList<int> SequenceWeights)
{
    /// <summary>The structure of profile 1 (profile commands with a structure line or mask in the profile file), else null.</summary>
    public ClustalWAlignmentStructure? Structure1 { get; init; }

    /// <summary>The structure of profile 2.</summary>
    public ClustalWAlignmentStructure? Structure2 { get; init; }

    /// <summary>Profile 1 has this many sequences (profile commands; 0 otherwise).</summary>
    public int Profile1Count { get; init; }

    /// <summary>Names in the order of <see cref="OutputIndex"/>.</summary>
    public IReadOnlyList<string> GetOutputNames() => OutputIndex.Select(i => Names[i]).ToList();

    /// <summary>Rows in the order of <see cref="OutputIndex"/>.</summary>
    public IReadOnlyList<string> GetOutputRows() => OutputIndex.Select(i => Rows[i]).ToList();
}

/// <summary>The outcome of a ClustalW command: the alignment (null when only guide trees were requested) and the guide-tree texts the program writes.</summary>
/// <param name="Alignment">The alignment, null for <see cref="ClustalWRunOptions.GuideTreeOnly"/>.</param>
/// <param name="GuideTree">The guide-tree file <c>-NEWTREE</c> would hold (null when a user tree was used or the run has two profiles).</param>
/// <param name="GuideTree1">Guide-tree file of profile 1 (<c>-NEWTREE1</c>; null when a user tree was given or the profile has one sequence).</param>
/// <param name="GuideTree2">Guide-tree file of profile 2 (<c>-NEWTREE2</c>).</param>
public sealed record ClustalWRunResult(ClustalWAlignment? Alignment, string? GuideTree, string? GuideTree1, string? GuideTree2);

/// <summary>The alignment file format of ClustalW (<c>-OUTPUT</c>).</summary>
public enum ClustalWOutputFormat
{
    /// <summary><c>CLUSTAL</c> (default): <c>.aln</c>, blocks of 60 columns with the conservation line.</summary>
    Clustal = 0,

    /// <summary><c>GCG</c>: MSF file (<c>.msf</c>) with checksums and weights.</summary>
    Gcg = 1,

    /// <summary><c>GDE</c>: <c>.gde</c>, '#' (DNA) or '%' (protein) headers; the structure and gap-penalty-mask lines of the profiles are written as <c>"SS_</c> / <c>"GM_</c> entries.</summary>
    Gde = 2,

    /// <summary><c>PHYLIP</c>: interleaved, names cut to 10 characters, blocks of 50 columns.</summary>
    Phylip = 3,

    /// <summary><c>PIR</c>: NBRF / PIR (<c>&gt;DL;</c> / <c>&gt;P1;</c>).</summary>
    Pir = 4,

    /// <summary><c>NEXUS</c>: interleaved DATA block.</summary>
    Nexus = 5,

    /// <summary><c>FASTA</c>: lines of 60 characters.</summary>
    Fasta = 6,
}

/// <summary>What ClustalW writes of the secondary-structure information of the profiles (<c>-SECSTROUT</c>).</summary>
public enum ClustalWStructureOutput
{
    /// <summary><c>STRUCTURE</c> (default): the secondary-structure lines (<c>!SS_</c>).</summary>
    Structure = 0,

    /// <summary><c>MASK</c>: the gap-penalty-mask lines (<c>!GM_</c>).</summary>
    Mask = 1,

    /// <summary><c>BOTH</c>.</summary>
    Both = 2,

    /// <summary><c>NONE</c>.</summary>
    None = 3,
}

/// <summary>The parameters of the alignment output of ClustalW (<c>-OUTPUT</c>, <c>-CASE</c>, <c>-SEQNOS</c>, <c>-SEQNO_RANGE</c>, <c>-RANGE</c>, <c>-SECSTROUT</c>).</summary>
public sealed record ClustalWOutputOptions
{
    /// <summary><c>-OUTPUT</c> (default CLUSTAL).</summary>
    public ClustalWOutputFormat Format { get; init; } = ClustalWOutputFormat.Clustal;

    /// <summary><c>-CASE=LOWER | UPPER</c> (GDE output only; default LOWER: the program's default <c>lowercase = true</c>).</summary>
    public bool Lowercase { get; init; } = true;

    /// <summary><c>-SEQNOS=ON | OFF</c> (CLUSTAL output only; default OFF): the residue count at the end of every line.</summary>
    public bool SequenceNumbers { get; init; }

    /// <summary><c>-SEQNO_RANGE=ON | OFF</c> (all formats; default OFF): <c>name/start-end</c> instead of the name.</summary>
    public bool SequenceNumberRange { get; init; }

    /// <summary><c>-RANGE=m,n</c>: the first and the last column written (null = the whole alignment).</summary>
    public (int From, int To)? Range { get; init; }

    /// <summary><c>-SECSTROUT=STRUCTURE | MASK | BOTH | NONE</c> (default STRUCTURE).</summary>
    public ClustalWStructureOutput Structure { get; init; } = ClustalWStructureOutput.Structure;
}

/// <summary>The gap-penalty parameters of the secondary-structure masks (<c>-HELIXGAP</c> ... <c>-STRANDENDOUT</c>, <c>-NOSECSTR1/2</c>).</summary>
public sealed record ClustalWStructureOptions
{
    /// <summary><c>-HELIXGAP=n</c> (1..9, default 4): gap penalty of helix core residues.</summary>
    public int HelixPenalty { get; init; } = 4;

    /// <summary><c>-STRANDGAP=n</c> (1..9, default 4): gap penalty of strand core residues.</summary>
    public int StrandPenalty { get; init; } = 4;

    /// <summary><c>-LOOPGAP=n</c> (1..9, default 1): gap penalty of loop regions.</summary>
    public int LoopPenalty { get; init; } = 1;

    /// <summary><c>-TERMINALGAP=n</c> (1..9, default 2): gap penalty at the ends of helices and strands.</summary>
    public int TerminalPenalty { get; init; } = 2;

    /// <summary><c>-HELIXENDIN=n</c> (0..3, default 3): residues inside a helix treated as terminal.</summary>
    public int HelixEndInside { get; init; } = 3;

    /// <summary><c>-HELIXENDOUT=n</c> (0..3, default 0): residues outside a helix treated as terminal.</summary>
    public int HelixEndOutside { get; init; }

    /// <summary><c>-STRANDENDIN=n</c> (0..3, default 1): residues inside a strand treated as terminal.</summary>
    public int StrandEndInside { get; init; } = 1;

    /// <summary><c>-STRANDENDOUT=n</c> (0..3, default 1): residues outside a strand treated as terminal.</summary>
    public int StrandEndOutside { get; init; } = 1;

    /// <summary>Set by <c>-NOSECSTR1</c> to false: the secondary-structure mask of profile 1 is not used.</summary>
    public bool UseStructure1 { get; init; } = true;

    /// <summary>Set by <c>-NOSECSTR2</c> to false: the secondary-structure mask of profile 2 is not used.</summary>
    public bool UseStructure2 { get; init; } = true;
}

/// <summary>Options of the guide-tree multiple aligners (<c>MultipleAlignProgressive</c>, <c>MultipleAlignIterative</c>): the options of <see cref="MsaOptions"/> plus the distance stage and the guide tree.</summary>
public record ProgressiveMsaOptions : MsaOptions
{
    /// <summary>Pairwise distance used for the guide tree of the progressive family (default <see cref="MsaDistanceMode.Linear"/>).</summary>
    public MsaDistanceMode Distance { get; init; } = MsaDistanceMode.Linear;

    /// <summary>
    /// <see cref="MsaDistanceMode.Linear"/> only: the work limit of the all-pairs distance stage (B09 WP78, R14-16). The stage aligns every pair of sequences (len1 · len2 cells per pair, O(k² L²) in all); when the sum of
    /// len1 · len2 over all pairs is above this many cells the distance of a pair is the k-tuple estimate instead (the share of the k-mers of the shorter sequence that the longer one has, MUSCLE's k-mer
    /// similarity, converted to the identity scale: d = 1 − F<sup>1/k</sup>; k = 8 for nucleotides, 4 for a protein matrix), which costs O(len1 + len2) per pair and gives the same ranking of the pairs on
    /// related sequences. Pairs with a sequence shorter than 2k keep the alignment distance. Default 400 million cells (about 1.3 s of the alignment stage; see docs/algorithms/Alignment/Multiple_Sequence_Alignment.md for
    /// the measurement); 0 or less = never (always the alignments, the historical behaviour for any size). Below the limit the distances are exactly those of the alignments.
    /// </summary>
    public long FastDistanceCells { get; init; } = 400_000_000;

    /// <summary>Guide-tree method of the progressive family (default <see cref="MsaGuideTree.Upgma"/>).</summary>
    public MsaGuideTree GuideTree { get; init; } = MsaGuideTree.Upgma;

    /// <summary>
    /// Parameters of the fast pairwise stage, used only with <see cref="MsaDistanceMode.ClustalWQuick"/> (default: ClustalW's DNA defaults).
    /// </summary>
    public ClustalWQuickOptions QuickPairwise { get; init; } = new();

    /// <summary>
    /// ClustalW PROTEIN pairwise parameters (matrix, <c>-PWGAPOPEN</c> / <c>-PWGAPEXT</c>, fast-stage options), used by
    /// <see cref="MsaDistanceMode.ClustalWAffine"/> / <see cref="MsaDistanceMode.ClustalWQuick"/> when the sequences are scored with a protein
    /// <see cref="MatrixScoring"/> (the DNA stage is used for nucleotide matrices and scalar scoring).
    /// </summary>
    public ClustalWProteinPairwiseOptions ClustalWProtein { get; init; } = new();
}
