using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The kind of sequences for MAFFT (<c>--nuc</c> / <c>--amino</c>; none = the program's guess from the letter frequencies).</summary>
public enum MafftSequenceType
{
    /// <summary>No option: nucleotide when more than 75 % of the letters are A, C, G, T, U or N (<c>io.c</c> <c>getnumlen</c>), otherwise amino acids.</summary>
    Auto = 0,
    /// <summary><c>--nuc</c> (<c>-D</c>).</summary>
    Nucleotide = 1,
    /// <summary><c>--amino</c> (<c>-P</c>).</summary>
    Amino = 2,
}

/// <summary>The scoring matrix of the amino-acid alignment (<c>--bl n</c>: BLOSUM, <c>--jtt n</c>: the JTT PAM matrix, <c>--tm n</c>: the transmembrane PAM matrix; for nucleotides <c>--jtt n</c> / <c>--kimura n</c> set the PAM distance of the Kimura matrix).</summary>
public enum MafftScoringModel
{
    /// <summary><c>--bl n</c> (default 62).</summary>
    Blosum = 0,
    /// <summary><c>--jtt n</c> / <c>--kimura n</c>.</summary>
    Jtt = 1,
    /// <summary><c>--tm n</c>.</summary>
    Tm = 2,
}

/// <summary>The memory option of the script: the default, <c>--memsave</c> (<c>-M -B</c>: the memory-saving dynamic program of the programs) or <c>--nomemsave</c> (<c>-N</c>).</summary>
public enum MafftMemoryMode
{
    /// <summary>No option.</summary>
    Default = 0,
    /// <summary><c>--memsave</c>.</summary>
    MemSave = 1,
    /// <summary><c>--nomemsave</c>.</summary>
    NoMemSave = 2,
}

/// <summary>The pairwise alignment information that the progressive alignment incorporates (<c>--globalpair</c> G-INS, <c>--localpair</c> L-INS, <c>--genafpair</c> E-INS); <see cref="None"/> = the 6-mer distances of FFT-NS.</summary>
public enum MafftPairLibrary
{
    /// <summary>No pair library (the default: FFT-NS / NW-NS).</summary>
    None = 0,
    /// <summary><c>--globalpair</c>: the global alignment (Needleman-Wunsch, <c>Galign11</c>) of every pair.</summary>
    Global = 1,
    /// <summary><c>--localpair</c>: the local alignment (Smith-Waterman, <c>Lalign11</c>) of every pair.</summary>
    Local = 2,
    /// <summary><c>--genafpair</c> / <c>--localgenafpair</c>: the local alignment with generalised affine gap costs (<c>genalign11</c>) of every pair.</summary>
    Genaf = 3,
    /// <summary><c>--fastapair</c>: the optimised score and the local alignments of FASTA 3.4 (<c>fasta34 -z3 -m10</c>, ported) of every sequence with the ones that follow it (<c>dndfast7</c>): the scores give the guide tree, the alignments the pair library (F-INS-1).</summary>
    Fasta = 4,
    /// <summary><c>--fastswpair</c>: as <see cref="Fasta"/> with the Smith-Waterman alignment of <c>fasta34 -A</c> (H-INS-1).</summary>
    FastaSmithWaterman = 5,
    /// <summary><c>--hybridpair</c> (no fragment option): <c>pairlocalalign -Y</c> (the local alignments of the pairs: <c>hat3</c>), <c>disttbfast -E 1 -T -y</c> (the 6-mer distances: <c>hat2</c>) and <c>tbfast</c> that reads both files (Hybrid-INS-1).</summary>
    Hybrid = 6,
    /// <summary><c>--lastpair</c> (no fragment option): <c>pairlocalalign -R</c> (the alignments of the program <c>lastal</c> of LAST for every pair: <c>hat2</c>, <c>hat3</c>; see <see cref="MafftOptions.LastPrograms"/>) and <c>tbfast</c> that loads both files (Last-INS-1).</summary>
    Last = 7,
}

/// <summary>How the distances between the new and the old sequences are made when fragments are added (<c>--addfragments</c> and the other options of the program <c>addsingle</c>; <c>--6merpair</c> ... <c>--hybridpair</c>).</summary>
public enum MafftFragmentPairs
{
    /// <summary><c>--6merpair</c> / <c>--10merpair</c> (<c>addsingle -W n</c>): the k-mer distances of the program itself.</summary>
    KmerPair = 0,
    /// <summary><c>--6mermultipair</c> / <c>--10mermultipair</c> (<c>addsingle -d -W n</c>): as <see cref="KmerPair"/>, the old sequences are also compared by the score of their alignment.</summary>
    KmerMultiPair = 1,
    /// <summary>No pair option, or <c>--multipair</c>: <c>dndpre -M 2</c> (old sequences) and <c>pairlocalalign -Y</c> (old against new) write the distance files.</summary>
    MultiPair = 2,
    /// <summary><c>--hybridpair</c>: <c>pairlocalalign -Y</c>, <c>disttbfast -E 1 -T -y</c>.</summary>
    HybridPair = 3,
    /// <summary><c>--localpair</c> with fragments: <c>pairlocalalign</c> of every pair (the local alignments are kept in <c>hat3</c>).</summary>
    LocalPair = 4,
    /// <summary><c>--lastpair</c> with fragments: <c>pairlocalalign -R</c> (LAST), then <c>addsingle</c>.</summary>
    LastPair = 5,
    /// <summary><c>--lastmultipair</c> with fragments: <c>dndpre -M 2 ... -Q</c> (which stops: the option <c>-Q</c> is not known to <c>dndpre</c>), then <c>pairlocalalign -r</c> and <c>addsingle -d</c>.</summary>
    LastMultiPair = 6,
}

/// <summary>The tree-splitting methods of MAFFT for very many sequences (<c>splittbfast</c>): the sequences are classified by their distance to a few representatives, recursively.</summary>
public enum MafftPartTree
{
    /// <summary>None (the progressive / iterative methods).</summary>
    None = 0,
    /// <summary><c>--parttree</c>: the 6-mer distance.</summary>
    KmerDistance = 1,
    /// <summary><c>--dpparttree</c>: the score of the global alignment of every sequence with the representatives (dynamic programming).</summary>
    DynamicProgramming = 2,
    /// <summary><c>--fastaparttree</c>: the optimised score of FASTA 3.4 (<c>fasta34</c>, ported).</summary>
    Fasta = 3,
}

/// <summary>
/// The guide tree of the progressive stage when it is not the UPGMA tree of the full distance matrix (<c>mafft --treein</c>, <c>--pileup</c>, <c>--randomchain</c>, <c>--large</c> / <c>--memsavetree</c>, <c>--memsavetreex</c>, <c>--stepadd</c>, <c>--youngestlinkage</c>): the script writes the method to the file <c>_guidetree</c>,
/// the program <c>disttbfast</c> (or <c>tbfast</c> with a pair library) reads it with <c>-U</c>.
/// </summary>
public enum MafftGuideTree
{
    /// <summary>The UPGMA tree of the full matrix of distances (the default).</summary>
    Default = 0,
    /// <summary><c>--treein file</c>: the tree of <see cref="MafftOptions.GuideTreeText"/> (one node per line: two sequence numbers and two branch lengths, the output of <c>newick2mafft.rb</c>).</summary>
    UserTree = 1,
    /// <summary><c>--pileup</c>: a chain, the sequences are added in the order of the input (6-mer distance methods only).</summary>
    Pileup = 2,
    /// <summary><c>--randomchain</c>: a chain in a random order (<c>rand()</c> of glibc seeded with <see cref="MafftOptions.RandomSeed"/>).</summary>
    RandomChain = 3,
    /// <summary><c>--large</c> / <c>--memsavetree</c>: every sequence is added above its nearest earlier sequence (<c>compacttree</c>; no distance matrix in memory).</summary>
    MemSaveTree = 4,
    /// <summary><c>--memsavetreex</c>: <c>compacttree</c> with a limit of RAM for the matrix (the script's 20GB: the UPGMA of the full matrix as long as it fits).</summary>
    MemSaveTreeX = 5,
    /// <summary><c>--stepadd</c>: the same tree as <see cref="MemSaveTree"/> (the message differs).</summary>
    StepAdd = 6,
    /// <summary><c>--youngestlinkage</c>: the UPGMA over the nearest neighbours only.</summary>
    YoungestLinkage = 7,
}

/// <summary>
/// Options of the MAFFT 7.505 command <c>mafft [--anysymbol] [--nuc | --amino] [--op x] [--ep y] [--bl n | --jtt n | --tm n] [--retree n] [--fft | --nofft] [--6merpair | --10merpair] ... input</c>: the plain progressive
/// methods (FFT-NS-1, FFT-NS-2, NW-NS-1, NW-NS-2) for any number of sequences. The options that need the later stages of the port (<c>--maxiterate</c>, <c>--localpair</c>, <c>--globalpair</c>,
/// <c>--genafpair</c>) are not options of this record.
/// </summary>
public sealed record MafftOptions
{
    /// <summary><c>--nuc</c> / <c>--amino</c>.</summary>
    public MafftSequenceType SequenceType { get; init; } = MafftSequenceType.Auto;

    /// <summary><c>--op</c>: gap opening penalty (default 1.53; the program passes <c>-f -op</c>, scaled by 0.6, or by 1.8 for nucleotides).</summary>
    public double GapOpen { get; init; } = 1.53;

    /// <summary><c>--ep</c>: offset value (default 0.0; <c>-h -ep</c>: subtracted from the scaled substitution scores, which makes the matrix more positive when negative values are given, i.e. favours gaps less).</summary>
    public double Offset { get; init; }

    /// <summary><c>--bl</c>: BLOSUM matrix for amino acids (30, 45, 50, 62 (default), 80, 90, 100 or 0).</summary>
    public int Blosum { get; init; } = 62;

    /// <summary><c>--anysymbol</c> (default true here: letters outside the alphabet are replaced by X / n for the alignment and restored in the output, as T-Coffee's <c>mafft_pair</c> calls the program). False is the plain command: letters are upper-cased (amino acids) or lower-cased (nucleotides) and unusual letters make the program stop.</summary>
    public bool AnySymbol { get; init; } = true;

    /// <summary><c>--retree n</c> (accepted; the tree of two sequences is trivial).</summary>
    public int Retree { get; init; } = 2;

    /// <summary><c>--fft</c> (default) / <c>--nofft</c> (accepted; verified equal on the committed pairs, see the Evidence).</summary>
    public bool Fft { get; init; } = true;

    /// <summary><c>--opdist x</c>: the gap opening penalty of the distance calculation (<c>-V</c>; default = <see cref="GapOpen"/>).</summary>
    public double? GapOpenDistance { get; init; }

    /// <summary><c>--exp x</c>: gap extension penalty (the script passes <c>-g -x</c>; default 0).</summary>
    public double GapExtension { get; init; }

    /// <summary><c>--termgappenalty</c>: penalise gaps at the ends of the alignment (the script then leaves out <c>-O</c>).</summary>
    public bool TerminalGapPenalty { get; init; }

    /// <summary><c>--unweight</c> (<c>-u</c>): no weighting of the sequences by the guide tree.</summary>
    public bool Unweight { get; init; }

    /// <summary><c>--bl</c> / <c>--jtt</c> / <c>--tm</c>: the model of the scoring matrix of amino acids.</summary>
    public MafftScoringModel ScoringModel { get; init; }

    /// <summary>The PAM distance of <c>--jtt n</c> / <c>--tm n</c> / <c>--kimura n</c> (default 200).</summary>
    public int Pam { get; init; } = 200;

    /// <summary><c>--kappa n</c>: the transition / transversion ratio of the Kimura matrix of nucleotides (default 2).</summary>
    public int? Kappa { get; init; }

    /// <summary><c>--fmodel</c>: the frequencies of the letters of the data (not the model's) in the scoring matrix.</summary>
    public bool FrequencyModel { get; init; }

    /// <summary><c>--nwildcard</c>: <c>n</c> scores like a wildcard in nucleotide sequences.</summary>
    public bool WildcardN { get; init; }

    /// <summary><c>--fft</c> given explicitly (<c>-G</c>: the FFT is used for every pair of groups, not only those with more residues than sequences).</summary>
    public bool ForceFft { get; init; }

    /// <summary><c>--6merpair</c> (6, default) / <c>--10merpair</c> (10, nucleotides only): the length of the k-mers of the distance.</summary>
    public int TupleSize { get; init; } = 6;

    /// <summary><c>--mixedlinkage x</c> (<c>--averagelinkage</c> = 1, <c>--minimumlinkage</c> = 0): the share of the average linkage of the UPGMA guide tree (default 0.1).</summary>
    public double? TreeMixing { get; init; }

    /// <summary><c>--leavegappyregion</c> / <c>--legacygappenalty</c>: the gap penalty of MAFFT before 7.110 (not scaled by the gap frequency of the columns).</summary>
    public bool LegacyGapPenalty { get; init; }

    /// <summary><c>--quiet</c>: the progress messages are not collected.</summary>
    public bool Quiet { get; init; }

    /// <summary><c>--thread n</c> (default 0 = one thread; -1 = the number of processors). The alignment does not depend on it (verified: the 220 committed sets of <c>--thread 4</c> equal the single-thread runs); the program's progress text
    /// with several threads is written in the order the threads run, which is not reproducible, so the port writes the text of the single-thread run with the thread counts.</summary>
    public int Threads { get; init; }

    /// <summary><c>--threadit n</c>: the number of threads of the iterative refinement (default: the script's rule, <c>min(--thread, 8)</c>). With more than zero threads the installed program is NON-DETERMINISTIC (the branches are refined in the order the threads finish: ten runs of the same input gave ten different alignments, B09.md F284), so such a run is refused; <c>--thread 4 --threadit 0</c> is deterministic and equals the single-thread run.</summary>
    public int? ThreadsIterate { get; init; }

    /// <summary><c>--maxiterate n</c>: the number of cycles of iterative refinement (FFT-NS-i: <c>dndpre</c> distances of the first alignment and <c>dvtditr</c>, the tree-dependent refinement of every branch of the guide tree); default 0 = the progressive method. The script lets at most 16 iterations run.</summary>
    public int MaxIterate { get; init; }

    /// <summary><c>--globalpair</c> / <c>--localpair</c> / <c>--genafpair</c>: every pair of sequences is aligned first (all against all) and the blocks of these alignments constrain the progressive alignment (G-INS-1, L-INS-1, E-INS-1; with <see cref="MaxIterate"/> the refinement of G-INS-i, L-INS-i, E-INS-i).</summary>
    public MafftPairLibrary PairLibrary { get; init; }

    /// <summary><c>--weighti x</c>: the weight of the pair library against the substitution scores (default 2.7). 0 switches the constraints off (<c>tbfast -c 0</c>, F317).</summary>
    public double WeightI { get; init; } = 2.7;

    /// <summary><c>--gop x</c> / <c>--gep x</c> / <c>--gexp x</c>: gap opening penalty, offset and gap extension of the global pair alignments (as the script passes them, the opening and extension are negative; defaults -2.00, 0.10, -0.10).</summary>
    public double GlobalPairOpen { get; init; } = -2.00;
    /// <summary>See <see cref="GlobalPairOpen"/>: <c>--gep</c> (default 0.10).</summary>
    public double GlobalPairOffset { get; init; } = 0.10;
    /// <summary>See <see cref="GlobalPairOpen"/>: <c>--gexp</c> (default -0.10).</summary>
    public double GlobalPairExtension { get; init; } = -0.10;

    /// <summary><c>--lop x</c> / <c>--lep x</c> / <c>--lexp x</c>: gap opening penalty, offset and gap extension of the local pair alignments (defaults -2.00, 0.100, -0.100).</summary>
    public double LocalPairOpen { get; init; } = -2.00;
    /// <summary>See <see cref="LocalPairOpen"/>: <c>--lep</c> (default 0.100).</summary>
    public double LocalPairOffset { get; init; } = 0.100;
    /// <summary>See <see cref="LocalPairOpen"/>: <c>--lexp</c> (default -0.100).</summary>
    public double LocalPairExtension { get; init; } = -0.100;

    /// <summary><c>--LOP x</c> / <c>--LEXP x</c>: the second gap opening penalty and the extension of the generalised affine gap costs of <c>--genafpair</c> (defaults -6.00 and -0.000).</summary>
    public double GenafOpen { get; init; } = -6.00;
    /// <summary>See <see cref="GenafOpen"/>: <c>--LEXP</c> (default -0.000).</summary>
    public double GenafExtension { get; init; } = -0.000;

    /// <summary>
    /// <c>--auto</c>: the script chooses the strategy from the number of sequences and the longest sequence (without gaps): fewer than 100 sequences of fewer than 3 000 letters = L-INS-i (<c>--localpair --maxiterate 1000</c>); fewer than 200 of fewer than 1 000 =
    /// L-INS-2 iterations (<c>--localpair --maxiterate 2</c>); fewer than 500 of fewer than 10 000 = FFT-NS-i with 2 iterations and two trees; fewer than 20 000 = FFT-NS-2. The pair library, the iterations and the number of trees of the options are then replaced;
    /// 20 000 to 199 999 sequences use the memory-saving tree (<c>--large</c>, <see cref="MafftGuideTree.MemSaveTree"/>), 200 000 or more <c>--dpparttree</c> (shorter than 3 000 letters) or <c>--parttree</c>.
    /// </summary>
    public bool Auto { get; init; }

    /// <summary><c>--parttree</c> / <c>--dpparttree</c> / <c>--fastaparttree</c> (given after a pair library option, which it then replaces; <see cref="Auto"/> replaces it): <c>splittbfast</c> classifies the sequences by their distance to <see cref="PartSize"/> representatives, follows the UPGMA tree of the representatives recursively and aligns the groups in the order of the tree
    /// (<see cref="Retree"/> passes; every pass after the first takes its distances from the alignment of the previous one). <c>--maxiterate</c> and a pair library cannot be combined with it (the program stops).</summary>
    public MafftPartTree PartTree { get; init; }

    /// <summary><c>--partsize n</c>: the number of sequences that represent a group in <see cref="PartTree"/> (default 50).</summary>
    public int PartSize { get; init; } = 50;

    /// <summary><c>--groupsize n</c>: groups of <c>n</c> or more sequences are not aligned as a whole (with <see cref="PartTree"/>; the alignment is then not complete, the rows have different lengths; default -1 = everything is aligned; 1 = only the order of the sequences and the tree).</summary>
    public int GroupSize { get; init; } = -1;

    /// <summary><c>--treein</c> / <c>--pileup</c> / <c>--randomchain</c> / <c>--large</c> ...: the method that builds the guide tree (default: the UPGMA tree of the matrix of distances).</summary>
    public MafftGuideTree GuideTree { get; init; }

    /// <summary>The tree of <c>--treein file</c> (<see cref="MafftGuideTree.UserTree"/>): the lines of the file.</summary>
    public string? GuideTreeText { get; init; }

    /// <summary><c>--lhlimit n</c> (with <see cref="MafftGuideTree.MemSaveTree"/> / <see cref="MafftGuideTree.StepAdd"/> and a pair library: <c>nodepair</c>): only the <c>n</c> longest sequences (shuffled within each half) get local homologies with the others; default: all.</summary>
    public int LhLimit { get; init; } = int.MaxValue;

    /// <summary><c>--treeout</c>: the guide tree is returned in <see cref="MafftRunResult.TreeFile"/> (the file <c>input.tree</c> of the program).</summary>
    public bool TreeOut { get; init; }

    /// <summary><c>--nodeout</c>: as <see cref="TreeOut"/> with the numbers of the nodes (progressive methods only).</summary>
    public bool NodeOut { get; init; }

    /// <summary><c>--distout</c>: the matrix of distances is returned in <see cref="MafftRunResult.DistanceFile"/> (the file <c>input.hat2</c>).</summary>
    public bool DistOut { get; init; }

    /// <summary>The names of the sequences (the lines after '&gt;'): they appear in the files <see cref="MafftRunResult.TreeFile"/> and <see cref="MafftRunResult.DistanceFile"/> (<c>--treeout</c>, <c>--nodeout</c>, <c>--distout</c>); default <c>seq1</c>, <c>seq2</c> ...</summary>
    public IReadOnlyList<string>? Names { get; init; }

    /// <summary><c>--minimumweight x</c>: the smallest weight of a sequence in the progressive and iterative steps (<c>-W</c>; default 0.00001).</summary>
    public double MinimumWeight { get; init; } = 0.00001;

    /// <summary><c>--skipiterate x</c> (<c>dvtditr -E</c>): a group is not refined when ... (the program's "fix" threshold; default 0.0).</summary>
    public double FixThreshold { get; init; }

    /// <summary><c>--shiftpenalty x</c> (<c>-Q</c> of the programs; <c>--allowshift</c> sets 2.0): the penalty of a shift of the alignment (a jump to a better cell of an earlier part of the matrix) in units of the gap opening penalty; a value below 10
    /// switches the dynamic programs of the global pair method to the "trywarp" mode (default 100: no shift).</summary>
    public double ShiftPenalty { get; init; } = 100.0;

    /// <summary><c>--unalignlevel x</c> (<c>-s</c> of tbfast, disttbfast and dvtditr, <c>-u</c> of the pair programs; <c>--allowshift</c> sets 0.8): the level at which dissimilar sequences are left unaligned (scoring matrices by the distance of a join); default 0.0.</summary>
    public double UnalignLevel { get; init; }

    /// <summary><c>--randomseed n</c> (<c>dvtditr -t</c>; default 0).</summary>
    public int RandomSeed { get; init; }

    /// <summary><c>--bunkatsunashi</c> (<c>dvtditr -B</c>): the alignment is not divided at anchors before the iterative refinement.</summary>
    public bool NoSegmentation { get; init; }

    /// <summary><c>--sp</c> (<c>-@</c> of the programs, <c>alg = 'd'</c>): the dynamic program <c>D__align</c> / <c>D__align_variousdist</c> (<c>Dalignmm.c</c>), which weights the gap penalty of every cell by the probability that the gap is new
    /// (the distribution of the lengths of the gaps of both groups); no FFT, no terminal gap penalty, no memory-saving mode; the iterative refinement (<c>-@ -B -Z -z 1000</c>) runs on one thread.</summary>
    public bool SumOfPairs { get; init; }

    /// <summary>The iterative refinement of <c>--sp</c> (<c>$algoptit</c> = <c>-@ -B -Z -z 1000</c>): the script sets it only when <c>--maxiterate</c> was positive at the place where it handles <c>--sp</c> (an iteration that <c>--seed</c> adds later runs with the default dynamic program).</summary>
    public bool SumOfPairsRefinement { get; init; }

    /// <summary><c>--oneiteration</c> (<c>disttbfast -r</c>).</summary>
    public bool OneIteration { get; init; }

    /// <summary><c>--noscore</c> (<c>dvtditr -Z</c>): no score is calculated during the iterative refinement.</summary>
    public bool NoScore { get; init; }

    /// <summary><c>--bestfirst</c> / <c>--adhoc0</c> / <c>--adhoc1</c> / <c>--adhoc2</c> (<c>dvtditr -p</c>): BESTFIRST, BAATARI0, BAATARI1 or BAATARI2 (default); they differ only with several threads in the refinement.</summary>
    public string ParallelizationStrategy { get; init; } = "BAATARI2";

    /// <summary><c>--scoreout</c> (<c>-S -B</c>): the program reports the score of the alignment.</summary>
    public bool ScoreOut { get; init; }

    /// <summary><c>--outnum</c> (<c>-n</c>): the sequences are numbered in the files the programs exchange.</summary>
    public bool OutNumber { get; init; }

    /// <summary><c>--sparsepickup n</c> (<c>disttbfast -q</c>).</summary>
    public int SparsePickup { get; init; }

    /// <summary><c>--threadtb n</c>: the number of threads of the progressive step (default: the script's rule).</summary>
    public int? ThreadsProgressive { get; init; }

    /// <summary><c>--memsave</c> (<c>-M -B</c>) / <c>--nomemsave</c> (<c>-N</c>): the memory-saving dynamic program / never.</summary>
    public MafftMemoryMode MemoryMode { get; init; }

    /// <summary>
    /// <c>--addfragments</c> (1), <c>--addfull</c> (-1), <c>--addlong</c> (-2), <c>--addtotop</c> (-3): the last <see cref="AddCount"/> sequences are added one by one by the program <c>addsingle</c> (0 = <c>--add</c> / <c>--addprofile</c> / no addition).
    /// </summary>
    public int Fragment { get; init; }

    /// <summary>How the distances for <see cref="Fragment"/> are made.</summary>
    public MafftFragmentPairs FragmentPairs { get; init; }

    /// <summary>
    /// The programs <c>lastdb</c> / <c>lastal</c> of LAST that <c>--lastpair</c> / <c>--lastmultipair</c> start (the GPL-licensed port lives in a separate assembly; this assembly only knows the interface). <c>null</c>: LAST is not installed and the script
    /// stops with its "Install LAST" message.
    /// </summary>
    public ILastPrograms? LastPrograms { get; init; }

    /// <summary>The directory of the LAST programs (the path that <c>lastal</c> writes in front of its messages).</summary>
    public string LastDirectory { get; init; } = "/opt/last/bin";

    /// <summary>The library directory of MAFFT that <c>pairlocalalign</c> prints as <c>whereispairalign</c> (the Debian installation).</summary>
    public string LibraryDirectory { get; init; } = "/usr/lib/mafft/lib/mafft";

    /// <summary><c>--last_e</c>: the minimum score of the alignments of <c>lastal -e</c> (default 5000).</summary>
    public int LastE { get; init; } = 5000;

    /// <summary><c>--last_m</c>: the multiplier of <c>lastal -m</c> (default 3).</summary>
    public int LastM { get; init; } = 3;

    /// <summary><c>--last_subopt</c>: keep the suboptimal alignments of a pair too (<c>pairlocalalign -S</c>).</summary>
    public bool LastSubOptimal { get; init; }

    /// <summary><c>--last_once</c>: one call of <c>lastal</c> for all the sequences (<c>pairlocalalign -U</c>).</summary>
    public bool LastOnce { get; init; }

    /// <summary><c>--adjustdirection</c> (1) / <c>--adjustdirectionaccurately</c> (2): <c>makedirectionlist</c> decides for every (new) nucleotide sequence whether it is closer to the others reverse complemented, <c>setdirection</c> reverse complements it and puts _R_ in front of its title.</summary>
    public int AdjustDirection { get; init; }

    /// <summary><c>--anchors file</c>: the text of the file (one anchor per line: "i j starti endi startj endj score", counted from 1); the script gives <c>disttbfast -l</c>.</summary>
    public string? ExternalAnchors { get; init; }

    /// <summary><c>--skipanchorsremoterthan n</c> (<c>disttbfast -x n</c>): an anchor farther than n columns from the nearest stronger one is given up (default 1000).</summary>
    public int MaxAnchorSeparation { get; init; } = 1000;

    /// <summary><c>-Z</c> of <c>dndpre</c> / <c>pairlocalalign</c> (<c>$usenaivepairscore</c> of the script, <c>--addlong</c>): the distances come from the score of the alignment counted without gap penalties.</summary>
    public bool NaivePairScore { get; init; }

    /// <summary><c>--codonpos file</c>: the text of the file (one codon position per column of the reference sequence); the script gives <c>addsingle -R</c>.</summary>
    public string? CodonPositions { get; init; }

    /// <summary><c>--codonscore file</c>: the text of the file (the scores of the codon pairs); the script gives <c>addsingle -S</c>.</summary>
    public string? CodonScores { get; init; }

    /// <summary><c>--add</c> / <c>--addprofile</c> (<c>disttbfast -K -I n</c> / <c>-I n</c>): the last <see cref="AddCount"/> of the sequences are added to the alignment that the others form (they are aligned; the new ones are not).</summary>
    public int AddCount { get; init; }

    /// <summary><c>--keeplength</c> (<c>-Y</c>): the length of the alignment of the original sequences is kept: the residues of the added sequences that fall into new columns are deleted.</summary>
    public bool KeepLength { get; init; }

    /// <summary><c>--mapout</c> (1) / <c>--compactmapout</c> (2): the file that lists the deleted residues (<see cref="MafftRunResult.MapFile"/>).</summary>
    public int MapOut { get; init; }

    /// <summary><c>--smoothing</c> (<c>-p</c>).</summary>
    public bool Smoothing { get; init; }

    /// <summary><c>--addprofile</c>: the added sequences are an alignment themselves (<c>-I n</c> without <c>-K</c>).</summary>
    public bool AddProfile { get; init; }

    /// <summary>The names are the complete titles the programs hold after <c>replaceu</c> and <c>--outnum</c> (the script layer computes them).</summary>
    internal bool NamesAreProgramTitles { get; init; }

    /// <summary><c>--oldgenafpair</c>: the parameters of E-INS-i before version 7.243 (the offset and the extension of the pair stage as given, the alignment score of the pair).</summary>
    public bool OldGenafParameters { get; init; }

    /// <summary><c>--reorder</c> (accepted; the two sequences keep their input order).</summary>
    public bool Reorder { get; init; }

    /// <summary>
    /// <c>--aamatrix file</c> (the script writes <c>-b -1</c> and copies the file to <c>_aamtx</c>): the text of that file, the lower half of a 20 x 20 matrix with the letters in the first line, optionally a line "frequency" or "norescale";
    /// used when <see cref="Blosum"/> is -1.
    /// </summary>
    public string? AaMatrix { get; init; }

    /// <summary>
    /// <c>--text</c> (the script writes <c>-P -b -2 -a</c>, no <c>--fft</c> and reads the output with <c>f2cl -E</c>): the sequences are texts: every byte except gap, blank and line end is a letter of an alphabet of 256 characters (the case counts),
    /// the matrix scores +1 for a match and -1 for a mismatch (<see cref="AaMatrix"/>, <c>--textmatrix</c>: the lines "0x41 0x42 3.5" of the file <c>_aamtx</c> replace entries), and the frequencies are those of the data.
    /// </summary>
    public bool TextMode { get; init; }

    /// <summary>
    /// <c>--seed file</c> / <c>--seedtable file</c> (<c>hat3.seed</c> of the script): the blocks of the structural alignments of the seeds, one line per block "i j overlapaa score start1 end1 start2 end2 k" (the sequence numbers are 0-based and count the
    /// seed sequences first; <c>multi2hat3s</c> writes them for every seed alignment). With it the progressive stage of <c>disttbfast</c> is followed by <see cref="SeedCycle"/> - 1 passes of <c>tbfast -J</c> (the distances from the alignment, the
    /// blocks of kind 'k' weighted by the structure weights) and the iterative refinement uses the blocks too (<c>dvtditr -l</c>).
    /// </summary>
    public string? SeedTable { get; init; }

    /// <summary>With <see cref="SeedTable"/>: the number of tree builds of the script (at least 2; <c>disttbfast -E 1</c> is the first, <c>tbfast -J</c> the others).</summary>
    public int SeedCycle { get; init; } = 2;

    /// <summary>The standard error of the runs of <c>replaceu</c> that the script made before the first program (one per seed file and one for the input, with <c>--seed --anysymbol</c>), replacing the line of the single run.</summary>
    internal string? ReplaceuLogOverride { get; init; }

    /// <summary>
    /// <c>--merge table</c> (the script writes <c>-H 0</c> for every program and copies the file to <c>_subalignmentstable</c>): the text of that file, one line per subalignment with the numbers (1-based, in the input order) of its sequences,
    /// a negative number to keep the columns that only gaps have; <c>#</c> starts a comment. The sequences of one subalignment must be an alignment; the programs keep their columns and the subalignment is a group of the guide tree.
    /// </summary>
    public string? MergeTable { get; init; }

    /// <summary><c>--core</c>: only the core of the alignment is written (<c>setcore</c>; it needs the file hat2 of the run: the iterative methods, the pair libraries with <c>--maxiterate</c> or <c>--distout</c>).</summary>
    public bool CoreOut { get; init; }

    /// <summary><c>--coreext</c>: the next <see cref="CoreWindow"/> residues of every sequence on both sides of the core are written too.</summary>
    public bool CoreExtend { get; init; }

    /// <summary><c>--corethr x</c> (default 0.5; the text of the script's argument).</summary>
    public string CoreThreshold { get; init; } = "0.5";

    /// <summary><c>--corewin n</c> (default 100; the text of the script's argument).</summary>
    public string CoreWindow { get; init; } = "100";

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
}

/// <summary>The result of a run of <c>mafft</c>: the alignment and the text that the program writes to standard error (progress, strategy; empty with <c>--quiet</c>).</summary>
/// <param name="Alignment">Rows in the order of the output (the input order, or the order of the guide tree with <c>--reorder</c>).</param>
/// <param name="StandardError">The messages of the binaries and of the script, byte for byte.</param>
public sealed record MafftRunResult(MultipleAlignmentResult Alignment, string StandardError)
{
    /// <summary>The index in the input of every row of the output (the identity unless <c>--reorder</c> or a tree-splitting method changed the order).</summary>
    public int[] InputIndex { get; init; } = Array.Empty<int>();

    /// <summary>The file <c>_deletemap</c> of <see cref="MafftOptions.MapOut"/> (<c>name.map</c>), or null.</summary>
    public string? MapFile { get; init; }

    /// <summary>With <c>--compactmapout</c> the program titles of the added sequences are changed ("Nins:...|title"): the new title of every row of the output, null where it is unchanged.</summary>
    public string?[]? ChangedNames { get; init; }

    /// <summary>The file <c>notused</c> of <see cref="MafftOptions.SparsePickup"/> (the sequences that were not picked up), or null.</summary>
    public string? NotUsedFile { get; init; }

    /// <summary>The guide tree file of <see cref="MafftOptions.TreeOut"/> / <see cref="MafftOptions.NodeOut"/> (<c>input.tree</c>), or null.</summary>
    public string? TreeFile { get; init; }

    /// <summary>--core: the names of the rows as <c>setcore</c> writes them (with the leading '=' and the prefix of <c>replaceu</c>), or null.</summary>
    public string[]? CoreNames { get; init; }

    /// <summary>--core: the text of <c>setcore</c> on the terminal. The script does not redirect it to the progress file, so it is not part of <see cref="MultipleAlignmentResult"/>'s standard error here, and it shows with <c>--quiet</c> too; it follows the text of the programs.</summary>
    public string TerminalText { get; init; } = "";

    /// <summary>The distance matrix file of <see cref="MafftOptions.DistOut"/> (<c>input.hat2</c>), or null.</summary>
    public string? DistanceFile { get; init; }

    /// <summary>The text that <c>disttbfast</c> writes to standard output before the alignment (the scores of the test of the guide tree file <c>auto</c>, <c>varpairscore</c>): it stays in front of the FASTA output, the other formats do not show it.</summary>
    public string OutputPrefix { get; init; } = "";

    /// <summary>--merge: the file <c>GuideTree</c> that the script copies to the current directory when disttbfast stopped because a subalignment is no group of the guide tree.</summary>
    internal string? GuideTreeCopy { get; init; }

    /// <summary>With <see cref="MafftOptions.GroupSize"/> 2 or more the program writes <c>Group-n</c> / <c>Group-para</c> in front of the name of every sequence: that label per row of the output; otherwise null.</summary>
    public string[]? GroupLabels { get; init; }
}

public static partial class SequenceAligner
{
    #region MAFFT 7.505 (stage 1: the alignment of a pair of sequences)

    // Port of MAFFT 7.505 (Katoh & Standley 2013, MBE 30:772; Debian mafft_7.505.orig.tar.gz, core/*.c) for the command `mafft --anysymbol [--nuc|--amino] [--op --ep --bl] in > out` on TWO sequences:
    //   replaceu.c (alphabet guess by getnumlen_casepreserve, unusual letters -> X / n), disttbfast.c (getnumlen guess, constants(), treebase: the pair is aligned by G__align11 with outgap = 0 because the
    //   script passes -O), constants.c (BLOSUM / Kimura-PAM scoring matrices, penalties scaled to ints), blosum.c, Galign11.c (G__align11 + Atracking), restoreu.c (original letters back).
    // The Falign (FFT) step of the default FFT-NS-2 is not needed for pairs in the tested range (it falls back to the same dynamic program); see the Evidence for the measured limits.
    // Everything that MAFFT keeps in double is double here.

    private const string MafftAminoAlphabet = "ARNDCQEGHILKMFPSTWYVBZX.-J";
    private const string MafftNucAlphabet = "agctuAGCTUnNbdhkmnrsvwyx-O";

    /// <summary>
    /// MAFFT 7.505's alignment of two sequences (<c>mafft --anysymbol in.fa</c>): the program's scoring matrix (BLOSUM scaled by 600 / average, or the 200-PAM Kimura matrix for nucleotides), gap opening
    /// penalty, free end gaps (<c>-O</c>) and the dynamic program <c>G__align11</c>. More than two sequences need the progressive stage (<c>A__align</c>, the 6-mer distance and the guide tree) and are not ported yet.
    /// </summary>
    /// <param name="sequences">Exactly two sequences of letters ('-' and '.' are accepted).</param>
    /// <param name="options">Parameters; default = <c>mafft --anysymbol</c>.</param>
    /// <exception cref="ArgumentException">An empty sequence among two, a letter outside the alphabet with <c>AnySymbol = false</c>, an unusable parameter.</exception>
    /// <exception cref="NotSupportedException">Three or more sequences (only the alignment of two sequences is available).</exception>
    public static MultipleAlignmentResult MultipleAlignMafft(IEnumerable<string> sequences, MafftOptions? options = null) => RunMafft(sequences, options).Alignment;

    /// <summary>
    /// <see cref="MultipleAlignMafft"/> with the messages of the program: the same alignment, and the standard error of <c>mafft</c> (the banner lines of the binaries <c>disttbfast</c> ... and the strategy text of the script,
    /// empty with <see cref="MafftOptions.Quiet"/>).
    /// </summary>
    public static MafftRunResult RunMafft(IEnumerable<string> sequences, MafftOptions? options = null)
    {
        var o = options ?? new MafftOptions();
        var list = new List<string>();
        foreach (var s in sequences ?? throw new ArgumentNullException(nameof(sequences)))
            list.Add(s ?? throw new ArgumentException("A sequence is null.", nameof(sequences)));
        // Empty input: the program accepts it ("Warning: Only 0 sequence found.", exit status 0, empty output: installed mafft 7.505), so does the port (Kalign, POA and ProbCons refuse it, as their
        // programs do: see EmptyInputMessage). One sequence: the program writes it back ("Only 1 sequence found" is a warning), here too.
        if (list.Count == 0)
            return new MafftRunResult(MultipleAlignmentResult.Empty, "");
        if (list.Count == 1)
            return new MafftRunResult(MafftEngine.SingleSequence(list[0], o), "");
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        return MafftEngine.AlignMany(list, o);
    }

    internal static class MafftEngine
    {
        /// <summary>The 26 x 26 integer scoring matrix of <c>constants()</c> (<c>n_dis</c>) in the order of the program's alphabet (<c>ARNDCQEGHILKMFPSTWYVBZX.-J</c> for amino acids, <c>agctuAGCTUnNbdhkmnrsvwyx-O</c> for nucleotides) and the integer gap opening penalty / extension.</summary>
        internal static (int[,] Matrix, int Penalty, int PenaltyEx) ScoringMatrix(bool nucleotide, MafftOptions o)
        {
            var g = MafftGlobals.Parse(MafftScript.DisttbfastArguments(o, 3));
            g.Dorp = nucleotide ? 'd' : 'p';
            g.Constants(Array.Empty<string>());
            var m = new int[26, 26];
            for (int i = 0; i < 26; i++)
                for (int j = 0; j < 26; j++)
                    m[i, j] = g.NDis[i][j];
            return (m, g.Penalty, g.PenaltyEx);
        }

        // -------------------------------------------------------------------- front end (replaceu.c, io.c)

        /// <summary>countATGC over the sequences (<c>getnumlen*</c>): the letters A T G C U N (any case) among all letters; more than 75 % = nucleotides (replaceu guesses on the original letters, disttbfast again on the replaced ones).</summary>
        internal static bool GuessNucleotide(IReadOnlyList<string> seqs, MafftOptions o)
        {
            if (o.SequenceType == MafftSequenceType.Nucleotide)
                return true;
            if (o.SequenceType == MafftSequenceType.Amino)
                return false;
            int atgc = 0, total = 0;
            foreach (var s in seqs)
            {
                if (total >= 1000000)
                    break;
                int n = 0, a = 0;
                foreach (char ch in s)
                {
                    if (!IsAsciiLetter(ch))
                        continue;
                    n++;
                    char c = char.ToLowerInvariant(ch);
                    if (c is 'a' or 't' or 'g' or 'c' or 'u' or 'n')
                        a++;
                }
                atgc += a;
                total += n;
            }
            return (double)atgc / total > 0.75; // NaN (no letters) is false: amino acids
        }

        /// <summary>The residues of a sequence as <c>countlen</c> and the plain mode of the programs read it (letters, '.', '-', '*' as '-').</summary>
        internal static string PlainResidues(string raw) => LoadResidues(raw, false);

        /// <summary>The residues as <c>replaceu</c> (<paramref name="anySymbol"/>) or the plain programs read them.</summary>
        internal static string PlainOrAnyResidues(string raw, bool anySymbol) => LoadResidues(raw, anySymbol);

        /// <summary>The standard error of <c>replaceu</c> (<c>--anysymbol</c>): the input file and "n x longest - shortest d|p".</summary>
        internal static string ReplaceuLog(IReadOnlyList<string> originals, bool nucleotide) =>
            "inputfile = orig\n" + originals.Count + " x " + originals.Max(x => x.Length) + " - " + originals.Min(x => x.Length) + " " + (nucleotide ? 'd' : 'p') + "\n";

        private static bool IsAsciiLetter(char c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

        /// <summary>replaceu.c: letters outside the usual alphabet are replaced (X for amino acids, n for nucleotides), the rest is upper- / lower-cased.</summary>
        internal static string ReplaceUnusual(string s, bool nucleotide)
        {
            const string usualP = "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";
            const string usualN = "ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (char ch in s)
            {
                if (ch >= 256)
                    throw new ArgumentException("MAFFT works on bytes; the sequence contains a character above U+00FF.", nameof(s));
                if (nucleotide)
                    sb.Append(usualN.Contains(ch) ? char.ToLowerInvariant(ch) : 'n');
                else
                    sb.Append(usualP.Contains(ch) ? char.ToUpperInvariant(ch) : 'X');
            }
            return sb.ToString();
        }

        /// <summary>One sequence: the program writes it back without aligning (installed mafft 7.505, "Warning: Only 1 sequence found."): the residues of the input, gaps removed; the original letters with <c>--anysymbol</c>, otherwise nucleotides in lower and amino acids in upper case.</summary>
        internal static MultipleAlignmentResult SingleSequence(string sequence, MafftOptions o)
        {
            string residues = NoGaps(LoadResidues(sequence, o.AnySymbol));
            if (!o.AnySymbol)
                residues = GuessNucleotide(new[] { residues }, o) ? residues.ToLowerInvariant() : residues.ToUpperInvariant();
            return new MultipleAlignmentResult(new[] { residues }, BuildConsensus(new List<string> { residues }, residues.Length, anySymbol: true), 0) { ExactTotalScore = 0 };
        }

        /// <summary>Three or more sequences: <c>replaceu</c>, <c>disttbfast</c> (<see cref="MafftRun.Disttbfast"/>), <c>restoreu</c>.</summary>
        internal static MafftRunResult AlignMany(IReadOnlyList<string> list, MafftOptions o, bool epilogue = true)
        {
            using var aaMatrix = MafftGlobals.UseAaMatrix(o.AaMatrix);
            using var mergeTable = MafftGlobals.UseSubalignmentsTable(o.MergeTable);
            return AlignManyCore(list, o, epilogue);
        }

        private static MafftRunResult AlignManyCore(IReadOnlyList<string> list, MafftOptions o, bool epilogue)
        {
            if (o.TextMode)
            {
                if (o.AnySymbol)
                    throw new ArgumentException("--text and --anysymbol cannot be combined (MAFFT: 'The combination of --text and ( --anysymbol or --preservecase ) is impossible.').", nameof(o));
                o = o with { ScoringModel = MafftScoringModel.Blosum, Blosum = -2, SequenceType = MafftSequenceType.Amino, FrequencyModel = true };
            }
            if (o.Blosum is not (30 or 45 or 50 or 62 or 80 or 90 or 100 or 0) && !(o.Blosum == -1 && o.AaMatrix != null) && !(o.Blosum == -2 && o.TextMode))
                throw new ArgumentException($"--bl must be 30, 45, 50, 62, 80, 90, 100 or 0, got {o.Blosum}.", nameof(o));
            if (!double.IsFinite(o.GapOpen) || !double.IsFinite(o.Offset) || !double.IsFinite(o.GapExtension))
                throw new ArgumentException("--op, --ep and --exp must be finite.", nameof(o));
            if (o.NodeOut && o.MaxIterate > 0 && o.Retree != 0)
                throw new ArgumentException("The --nodeout option supports only the progressive method (--maxiterate 0) for now (MAFFT stops).");
            if (o.Retree == 0 && o.PairLibrary != MafftPairLibrary.None)
                o = o with { WeightI = 0.0, MaxIterate = 0 }; // the script: with --retree 0 "weighti" is 0.0 and no iteration runs
            string[] originals;
            if (o.TextMode)
            {
                try
                {
                    originals = list.Select(MafftRun.TextResidues).ToArray();
                }
                catch (MafftProgramFailureException e)
                {
                    throw new MafftProgramFailureException(TextReadPrefix(o, list.Count) + e.Text, e.ExitCode);
                }
            }
            else
            {
                originals = list.Select(x => LoadResidues(x, o.AnySymbol)).ToArray();
            }
            if (o.Auto)
                o = ApplyAuto(o, list.Count, originals.Max(x => NoGaps(x).Length));
            bool nuc1 = GuessNucleotide(originals, o);
            var log = new System.Text.StringBuilder();
            if (o.AnySymbol)
                log.Append(o.ReplaceuLogOverride ?? ReplaceuLog(originals, nuc1));
            string[] replaced = originals;
            if (!o.TextMode)
                replaced = o.AnySymbol
                    ? originals.Select(s => ReplaceUnusual(s, nuc1)).ToArray()
                    : originals.Select(s => nuc1 ? s.ToLowerInvariant() : s.ToUpperInvariant()).ToArray();
            bool nuc = GuessNucleotide(replaced, o);
            if (!o.TextMode)
                replaced = nuc ? replaced.Select(s => s.ToLowerInvariant()).ToArray() : replaced.Select(s => s.ToUpperInvariant()).ToArray();
            var header = o.Names ?? Enumerable.Range(0, list.Count).Select(i => "seq" + (i + 1).ToString(System.Globalization.CultureInfo.InvariantCulture)).ToArray();
            if (header.Count != list.Count)
                throw new ArgumentException("Names must have one entry per sequence.", nameof(o));
            string TitlePrefix(int i) => o.AnySymbol && !o.NamesAreProgramTitles ? "_os_" + i.ToString(System.Globalization.CultureInfo.InvariantCulture) + "_oe_" : "";
            var programNames = header.Select((h, i) => "=" + TitlePrefix(i) + h).ToArray();
            if (o.PartTree != MafftPartTree.None)
                return AlignPartTree(list.Count, originals, replaced, nuc, o, log, epilogue, programNames);
            string?[]? directionNames = null;
            if (o.AdjustDirection > 0)
                (replaced, programNames, originals, directionNames) = AdjustDirections(o, replaced, programNames, originals, header, log);
            MafftRun.Result result;
            string guideTreeText = MafftScript.GuideTreeFile(o);
            string outputPrefix = "";
            string? treeFile = null, distFile = null;
            double[][]? pairHat2 = null;
            MafftLocalHom[][]? pairHat3 = null;
            bool[]? seedKozo = null;
            if (o.PairLibrary is MafftPairLibrary.Fasta or MafftPairLibrary.FastaSmithWaterman)
            {
                // the script: dndfast7 [-A] < infile (hat2, hat3), then tbfast without a pair stage: it loads hat3 and hat2 (hat3.seed in front of hat3)
                var d7 = MafftRun.Dndfast7(programNames, replaced, o.PairLibrary == MafftPairLibrary.FastaSmithWaterman);
                log.Append(d7.Err);
                if (d7.Segfault)
                    throw new MafftProgramFailureException(log + MafftProgram.ProgramPath + ": line 2816: <pid> Segmentation fault      \"$prefix/dndfast7\" $swopt < infile > /dev/null 2>> \"$progressfile\"\n", 1);
                if (d7.ExitCode != 0)
                    throw new MafftProgramFailureException(log.ToString(), d7.ExitCode);
                string hat3Text = (o.SeedTable ?? "") + d7.Hat3;
                var (_, fastaTbArgs) = MafftScript.TbfastArguments(o, list.Count);
                var tg = MafftGlobals.ParseTbfast(fastaTbArgs);
                if (tg.Dorp == '\0')
                    tg.Dorp = nuc ? 'd' : 'p';
                var trun = new MafftRun(tg) { GuideTreeText = MafftScript.GuideTreeFile(o) };
                result = GuardRun(trun, log, () => trun.Tbfast(replaced, null, programNames, null, hat3Text, d7.Hat2));
                log.Append(trun.Log);
                treeFile = trun.TreeFile;
                distFile = trun.DistFile ?? d7.Hat2;
                if (o.WeightI >= 0.0001 && result.Rows.Length > 0)
                {
                    seedKozo = new bool[list.Count];
                    pairHat3 = MafftRun.ReadLocalHomTable2Half(hat3Text, list.Count, seedKozo);
                }

                pairHat2 = MafftRun.ReadHat2Full(d7.Hat2, list.Count);
            }
            else if (o.PairLibrary == MafftPairLibrary.Last)
            {
                // --lastpair without a fragment option (mafft.tmpl 2575 - 2588): pairlocalalign -R (the alignments of LAST for every pair: hat2, hat3), hat3.seed in front of hat3, tbfast that loads both files
                if (o.LastPrograms is null)
                    throw new NotSupportedException("--lastpair needs the programs lastdb and lastal of LAST: supply an implementation of ILastPrograms in MafftOptions.LastPrograms or MafftInstallation.Last.");
                var pg = MafftGlobals.ParsePairlocalalign(MafftScript.PairlocalalignLastArguments(o));
                pg.Alg = 'R';
                if (pg.Dorp == '\0')
                    pg.Dorp = nuc ? 'd' : 'p';
                var prun = new MafftRun(pg);
                var pairs = GuardRun(prun, log, () => prun.PairlocalalignLastProgram(replaced, programNames, o.LastPrograms, o.LastDirectory, o.LibraryDirectory, o.LastE, o.LastM, o.LastSubOptimal));
                log.Append(prun.Log);
                string hat3Text = (o.SeedTable ?? "") + pairs.Hat3;
                var (_, lastTbArgs) = MafftScript.TbfastArguments(o, list.Count);
                var tg = MafftGlobals.ParseTbfast(lastTbArgs);
                if (tg.Dorp == '\0')
                    tg.Dorp = nuc ? 'd' : 'p';
                var trun = new MafftRun(tg) { GuideTreeText = MafftScript.GuideTreeFile(o) };
                result = GuardRun(trun, log, () => trun.Tbfast(replaced, null, programNames, null, hat3Text, pairs.Hat2));
                log.Append(trun.Log);
                treeFile = trun.TreeFile;
                distFile = trun.DistFile ?? pairs.Hat2;
                if (o.WeightI >= 0.0001 && result.Rows.Length > 0)
                {
                    seedKozo = new bool[list.Count];
                    pairHat3 = MafftRun.ReadLocalHomTable2Half(hat3Text, list.Count, seedKozo);
                }

                pairHat2 = MafftRun.ReadHat2Full(pairs.Hat2, list.Count);
            }
            else if (o.PairLibrary == MafftPairLibrary.Hybrid)
            {
                // --hybridpair without a fragment option (mafft.tmpl 2615 - 2624): pairlocalalign -Y (the blocks of the old sequences against the new ones: hat3; nothing without --add), hat3.seed in front of hat3,
                // disttbfast -E 1 -T -y (the 6-mer distances: hat2) and tbfast, which loads both files
                var pg = MafftGlobals.ParsePairlocalalign(MafftScript.PairlocalalignFragmentArguments(o));
                if (pg.Dorp == '\0')
                    pg.Dorp = nuc ? 'd' : 'p';
                var prun = new MafftRun(pg);
                var pairs = GuardRun(prun, log, () => prun.PairlocalalignProgram(replaced, programNames));
                log.Append(prun.Log);
                string hat3Text = (o.SeedTable ?? "") + pairs.Hat3;
                var hg = MafftGlobals.Parse(MafftScript.HybridDisttbfastArguments(o));
                if (hg.Dorp == '\0')
                    hg.Dorp = nuc ? 'd' : 'p';
                var hrun = new MafftRun(hg) { GuideTreeText = guideTreeText, DistFileRewrittenLater = true };
                GuardRun(hrun, log, () => hrun.Disttbfast(replaced, programNames));
                log.Append(hrun.Log);
                string hat2Text = hrun.DistFile ?? pairs.Hat2;
                var (_, hybridTbArgs) = MafftScript.TbfastArguments(o, list.Count);
                var tg = MafftGlobals.ParseTbfast(hybridTbArgs);
                if (tg.Dorp == '\0')
                    tg.Dorp = nuc ? 'd' : 'p';
                var trun = new MafftRun(tg) { GuideTreeText = MafftScript.GuideTreeFile(o) };
                result = GuardRun(trun, log, () => trun.Tbfast(replaced, null, programNames, null, hat3Text, hat2Text));
                log.Append(trun.Log);
                treeFile = trun.TreeFile;
                distFile = trun.DistFile ?? hat2Text;
                if (o.WeightI >= 0.0001 && result.Rows.Length > 0)
                {
                    seedKozo = new bool[list.Count];
                    pairHat3 = MafftRun.ReadLocalHomTable2Half(hat3Text, list.Count, seedKozo);
                }

                pairHat2 = MafftRun.ReadHat2Full(hat2Text, list.Count);
            }
            else if (o.PairLibrary != MafftPairLibrary.None)
            {
                var (pairArgs, tbArgs) = MafftScript.TbfastArguments(o, list.Count);
                NodePairData? nodePair = null;
                string tbGuideTree = MafftScript.GuideTreeFile(o);
                if (o.GuideTree is MafftGuideTree.MemSaveTree or MafftGuideTree.StepAdd)
                {
                    // --large / --memsavetree / --stepadd with a pair library: nodepair builds the tree and the pairwise alignments of every node, tbfast -U reads them
                    MafftGlobals ng;
                    try
                    {
                        ng = MafftGlobals.ParseNodepair(MafftScript.NodepairArguments(o));
                    }
                    catch (MafftProgramFailureException e)
                    {
                        throw new MafftProgramFailureException(log.ToString() + e.Text, e.ExitCode);
                    }
                    if (ng.Dorp == '\0')
                        ng.Dorp = nuc ? 'd' : 'p';
                    var nrun = new MafftRun(ng) { GuideTreeText = o.GuideTreeText != null ? MafftScript.NormalizeScriptFile(o.GuideTreeText) : "" };
                    nodePair = GuardRun(nrun, log, () => nrun.NodePair(replaced, programNames));
                    log.Append(nrun.Log);
                    tbGuideTree = "nodepair\n";
                    treeFile = nodePair.TreeFile;
                }
                var tg = MafftGlobals.ParseTbfast(tbArgs);
                if (tg.Dorp == '\0')
                    tg.Dorp = nuc ? 'd' : 'p';
                var trun = new MafftRun(tg) { GuideTreeText = tbGuideTree };
                result = GuardRun(trun, log, () => trun.Tbfast(replaced, MafftGlobals.ParsePairlocalalign(pairArgs), programNames, nodePair, o.SeedTable));
                log.Append(trun.Log);
                pairHat2 = trun.Hat2Table;
                pairHat3 = trun.Hat3Table;
                seedKozo = trun.Kozo;
                treeFile = trun.TreeFile ?? treeFile;
                distFile = trun.DistFile;
            }
            else if (o.Fragment != 0)
            {
                // --addfragments / --addfull / --addlong / --addtotop: the programs of the chain write the files, the program addsingle adds the new sequences one by one
                MafftRun.AddSingleFiles? files = null;
                if (o.CodonPositions != null || o.CodonScores != null)
                {
                    string? codonPos = o.CodonPositions == null ? null : MafftScript.NormalizeScriptFile(o.CodonPositions);
                    string? codonScore = o.CodonScores == null ? null : MafftScript.NormalizeScriptFile(o.CodonScores);
                    files = new MafftRun.AddSingleFiles(CodonPos: codonPos, CodonScore: codonScore);
                }
                if (o.FragmentPairs is MafftFragmentPairs.LastPair or MafftFragmentPairs.LastMultiPair)
                {
                    if (o.LastPrograms is null)
                        throw new NotSupportedException("--lastpair / --lastmultipair need the programs lastdb and lastal of LAST: supply an implementation of ILastPrograms in MafftInstallation.Last.");
                    if (o.FragmentPairs == MafftFragmentPairs.LastMultiPair)
                    {
                        // dndpre $model -M 2 $addarg ... -Q $spfactor ...: the program does not know -Q ("options: Check source file !"), the script stops
                        try
                        {
                            MafftGlobals.ParseDndpre(MafftScript.DndpreLastArguments(o));
                        }
                        catch (MafftProgramFailureException e)
                        {
                            throw new MafftProgramFailureException(log + e.Text, e.ExitCode);
                        }

                        throw new NotSupportedException("dndpre accepted -Q: the chain of --lastmultipair (pairlocalalign -r, addsingle -d) is not ported.");
                    }

                    var lg = MafftGlobals.ParsePairlocalalign(MafftScript.PairlocalalignLastArguments(o));
                    lg.Alg = 'R';
                    lg.Nadd = o.AddCount;
                    if (lg.Dorp == '\0')
                        lg.Dorp = nuc ? 'd' : 'p';
                    var lrun = new MafftRun(lg);
                    var lastPairs = GuardRun(lrun, log, () => lrun.PairlocalalignLastProgram(replaced, programNames, o.LastPrograms, o.LastDirectory, o.LibraryDirectory, o.LastE, o.LastM, o.LastSubOptimal));
                    log.Append(lrun.Log);
                    files = (files ?? new MafftRun.AddSingleFiles()) with { Hat2 = lastPairs.Hat2, Hat3 = lastPairs.Hat3 };
                }
                else if (o.FragmentPairs is MafftFragmentPairs.MultiPair or MafftFragmentPairs.HybridPair or MafftFragmentPairs.LocalPair)
                {
                    string? hat2i = null;
                    if (o.FragmentPairs == MafftFragmentPairs.MultiPair)
                    {
                        var dg = MafftGlobals.ParseDndpre(MafftScript.DndpreFragmentArguments(o));
                        var drun = new MafftRun(dg);
                        GuardRun(drun, log, () => drun.Dndpre(replaced, programNames));
                        log.Append(drun.Log);
                        hat2i = drun.DistFile;
                    }

                    var pg = MafftGlobals.ParsePairlocalalign(MafftScript.PairlocalalignFragmentArguments(o));
                    if (pg.Dorp == '\0')
                        pg.Dorp = nuc ? 'd' : 'p';
                    var prun = new MafftRun(pg);
                    var pairs = GuardRun(prun, log, () => prun.PairlocalalignProgram(replaced, programNames));
                    log.Append(prun.Log);
                    string hat2 = pairs.Hat2;
                    if (o.FragmentPairs == MafftFragmentPairs.HybridPair)
                    {
                        var hg = MafftGlobals.Parse(MafftScript.HybridDisttbfastArguments(o));
                        if (hg.Dorp == '\0')
                            hg.Dorp = nuc ? 'd' : 'p';
                        var hrun = new MafftRun(hg) { GuideTreeText = guideTreeText, DistFileRewrittenLater = true };
                        GuardRun(hrun, log, () => hrun.Disttbfast(replaced, programNames));
                        log.Append(hrun.Log);
                        hat2 = hrun.DistFile ?? hat2;
                    }

                    files = o.FragmentPairs == MafftFragmentPairs.MultiPair
                        ? (files ?? new MafftRun.AddSingleFiles()) with { Hat2I = hat2i, Hat2N = hat2, Hat3 = pairs.Hat3 }
                        : (files ?? new MafftRun.AddSingleFiles()) with { Hat2 = hat2, Hat3 = pairs.Hat3 };
                }

                MafftGlobals ag;
                try
                {
                    ag = MafftGlobals.ParseAddsingle(MafftScript.AddsingleArguments(o));
                }
                catch (MafftProgramFailureException e)
                {
                    throw new MafftProgramFailureException(log + e.Text, e.ExitCode); // --sp: addsingle does not know -@ (the stop follows the output of the programs before it)
                }
                if (ag.Dorp == '\0')
                    ag.Dorp = nuc ? 'd' : 'p';
                var arun = new MafftRun(ag) { GuideTreeText = guideTreeText };
                result = GuardRun(arun, log, () => arun.AddSingle(replaced, programNames, files));
                log.Append(arun.Log);
                treeFile = arun.TreeFile;
                distFile = arun.DistFile;
            }
            else
            {
                var g = MafftGlobals.Parse(MafftScript.DisttbfastArguments(o, list.Count));
                if (g.Dorp == '\0')
                    g.Dorp = nuc ? 'd' : 'p';
                g.ExternalAnchorsText = o.ExternalAnchors == null ? null : MafftScript.NormalizeScriptFile(o.ExternalAnchors);
                var run = new MafftRun(g) { GuideTreeText = guideTreeText, DistFileRewrittenLater = MafftScript.Iterations(o) > 0 };
                result = GuardRun(run, log, () => run.Disttbfast(replaced, programNames));
                guideTreeText = run.GuideTreeText;
                outputPrefix = run.StdoutPrefix;
                treeFile = run.TreeFile;
                distFile = run.DistFile;
                log.Append(run.Log);
                if (o.SeedTable != null)
                {
                    // --seed / --seedtable: the later passes of the script are tbfast -J on the alignment of the pass before (distances from the alignment, the blocks of the seeds as constraints)
                    for (int pass = 1; pass < o.SeedCycle; pass++)
                    {
                        var tg = MafftGlobals.ParseTbfast(MafftScript.TbfastJArguments(o));
                        if (tg.Dorp == '\0')
                            tg.Dorp = nuc ? 'd' : 'p';
                        var trun = new MafftRun(tg);
                        var previousRows = result.Rows;
                        result = GuardRun(trun, log, () => trun.Tbfast(previousRows, null, programNames, null, o.SeedTable));
                        log.Append(trun.Log);
                        treeFile = trun.TreeFile ?? treeFile;
                        distFile = trun.DistFile ?? distFile;
                    }
                    seedKozo = new bool[list.Count];
                    try
                    {
                        if (o.WeightI >= 0.0001 && result.Rows.Length > 0) // fastathreshold < 0.0001: constraint = 0, dvtditr does not read hat3 (and it does not run when the disttbfast stopped)
                            pairHat3 = MafftRun.ReadLocalHomTable2Half(o.SeedTable, list.Count, seedKozo);
                    }
                    catch (MafftProgramFailureException e)
                    {
                        throw new MafftProgramFailureException(log + e.Text, e.ExitCode);
                    }
                }
            }
            if (result.Rows.Length == 0 && list.Count > 0)
            {
                // the program stopped after its message (the original sequences of --add are not aligned)
                if (MafftScript.Iterations(o) > 0 && o.PairLibrary == MafftPairLibrary.None)
                {
                    // the script goes on after the stop of disttbfast (exit status 0, empty file pre): dndpre stops on it
                    throw new MafftProgramFailureException(log + "At least 2 sequences should be input!\nOnly 0 sequence found.\n");
                }
                string quit = "";
                if (!o.Quiet)
                    quit = epilogue ? log + MafftScript.Epilogue(o, list.Count) : log.ToString();
                return new MafftRunResult(MultipleAlignmentResult.Empty, quit) { InputIndex = Array.Empty<int>(), GuideTreeCopy = result.GuideTreeCopy, TreeFile = o.TreeOut || o.NodeOut ? treeFile : null, DistanceFile = o.DistOut ? distFile : null };
            }
            int[]? picked = result.Picked;
            string? notUsed = result.NotUsed;
            var deletedRuns = result.Deleted;
            string? mapFile = result.DeleteMap;
            var changedNames = new string?[list.Count];
            if (directionNames != null)
                for (int k = 0; k < directionNames.Length; k++)
                    changedNames[k] = directionNames[k];
            if (result.NewNames != null)
                for (int k = 0; k < result.NewNames.Length; k++)
                    changedNames[list.Count - result.NewNames.Length + k] = result.NewNames[k];
            string[] activeNames = picked == null ? programNames : picked.Select(i => programNames[i]).ToArray();
            if (MafftScript.Iterations(o) > 0)
            {
                int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
                int threadsIt = o.ThreadsIterate ?? Math.Min(threads, 8);
                var previousOrder = result.Order;
                if (threadsIt > 0)
                    throw new NotSupportedException("The iterative refinement of MAFFT with threads (--thread n --maxiterate m without --threadit 0) is non-deterministic in the program (the branches are refined in the order the threads finish: repeated runs of one input give different alignments); use Threads = 0 or ThreadsIterate = 0, which equal the single-thread program.");
                double[][] hat2;
                if (pairHat2 != null)
                {
                    hat2 = pairHat2; // the distances of the pair library (tbfast wrote hat2 itself)
                }
                else
                {
                    var run2 = new MafftRun(MafftGlobals.ParseDndpre(MafftScript.DndpreArguments(o)));
                    hat2 = GuardRun(run2, log, () => run2.Dndpre(result.Rows, activeNames));
                    log.Append(run2.Log);
                    distFile = run2.DistFile; // dndpre rewrites the file hat2
                }
                var run3 = new MafftRun(MafftGlobals.ParseDvtditr(MafftScript.DvtditrArguments(o))) { GuideTreeText = guideTreeText };
                result = GuardRun(run3, log, () => run3.Dvtditr(result.Rows, hat2, pairHat3, activeNames, seedKozo));
                if (result.Order.Length == 0)
                    result = result with { Order = previousOrder }; // --skipiterate with a large value: dvtditr leaves the file 'order' as it was
                log.Append(run3.Log);
                treeFile = run3.TreeFile; // dvtditr writes the file infile.tree (fixed_musclesupg_double_treeout) whatever the options say
            }
            string[] rows = o.AnySymbol && !o.CoreOut
                ? result.Rows.Select((r, i) => Restore(DeleteRuns(originals[picked?[i] ?? i], deletedRuns, i), r)).ToArray()
                : result.Rows;
            string[]? coreNames = null;
            string coreText = "";
            if (o.CoreOut)
            {
                // setcore replaces restoreu: the alignment keeps the symbols of replaceu and the names keep '='
                var crun = new MafftRun(MafftGlobals.ParseSetcore(MafftScript.SetcoreArguments(o)));
                (rows, coreNames) = crun.SetCore(result.Rows, activeNames, distFile);
                coreText = crun.Log.ToString();
            }
            if (picked != null && o.Reorder)
                throw new NotSupportedException("--sparsepickup with --reorder: the program stops before it writes the file 'order'.");
            if (o.Reorder && rows.Length == result.Order.Length)
            {
                var previous = rows;
                var previousNames = coreNames;
                rows = result.Order.Select(i => previous[i]).ToArray();
                if (previousNames != null)
                    coreNames = result.Order.Select(i => previousNames[i]).ToArray();
            }
            string stderr = "";
            if (!o.Quiet)
                stderr = epilogue ? log + MafftScript.Epilogue(o, list.Count) : log.ToString();
            int[] inputIndex = picked ?? (o.Reorder ? result.Order.ToArray() : Enumerable.Range(0, rows.Length).ToArray());
            if (coreNames != null && rows.Length != inputIndex.Length)
                inputIndex = Enumerable.Range(0, rows.Length).ToArray();
            return new MafftRunResult(ToResult(rows, nuc), stderr)
            {
                CoreNames = coreNames,
                TerminalText = coreText,
                InputIndex = inputIndex,
                OutputPrefix = outputPrefix,
                TreeFile = o.TreeOut || o.NodeOut ? treeFile : null,
                DistanceFile = o.DistOut ? distFile : null,
                NotUsedFile = notUsed,
                MapFile = mapFile,
                ChangedNames = changedNames.Any(x => x != null) ? changedNames : null,
            };
        }

        /// <summary>
        /// The block of the script for <c>--adjustdirection</c> (<c>--adjustdirectionaccurately</c>): <c>makedirectionlist -F -C n -m -I nadd -i infile -t 0.00 -r 5000 -o a [-d] &gt; _direction</c> (its exit status is not looked at), then <c>setdirection [-H 0] -d _direction -i infile</c>
        /// (the sequences and the titles that the programs read) and, with <c>--anysymbol</c>, <c>setdirection ... -i orig -r</c> (the sequences that <c>restoreu</c> puts back; their titles stay).
        /// </summary>
        private static (string[] Replaced, string[] ProgramNames, string[] Originals, string?[] Renamed) AdjustDirections(MafftOptions o, string[] replaced, string[] programNames, string[] originals, IReadOnlyList<string> header, System.Text.StringBuilder log)
        {
            var dargs = new List<string> { "-F", "-C", (o.Threads < 0 ? Environment.ProcessorCount : o.Threads).ToString(System.Globalization.CultureInfo.InvariantCulture), "-m", "-I", o.AddCount.ToString(System.Globalization.CultureInfo.InvariantCulture), "-i", "infile", "-t", "0.00" };
            dargs.AddRange(o.AdjustDirection == 1 ? new[] { "-r", "5000", "-o", "a" } : new[] { "-r", "100", "-o", "a", "-d" });
            var dg = MafftGlobals.ParseMakedirectionlist(dargs);
            var drun = new MafftRun(dg);
            string[] lines;
            try
            {
                lines = drun.MakeDirectionList(replaced, programNames).Lines;
                log.Append(drun.Log);
            }
            catch (MafftProgramFailureException e)
            {
                // the program stops: the file _direction is empty
                log.Append(drun.Log).Append(e.Text);
                lines = Array.Empty<string>();
            }

            bool merge = o.MergeTable != null;
            var srun = new MafftRun(new MafftGlobals { Subalignment = merge });
            (string[] rows, string[] names) = GuardRun(srun, log, () => srun.SetDirectionProgram(replaced, programNames, lines, false, merge));
            log.Append(srun.Log);
            if (o.AnySymbol)
            {
                var origNames = header.Select(h => "=" + h).ToArray();
                var orun = new MafftRun(new MafftGlobals { Subalignment = merge });
                (originals, _) = GuardRun(orun, log, () => orun.SetDirectionProgram(originals, origNames, lines, true, merge));
                log.Append(orun.Log);
            }

            var renamed = new string?[names.Length];
            for (int i = 0; i < names.Length; i++)
            {
                string cut = names[i].Length > MafftProgramIo.NameLimit + 1 ? names[i][..(MafftProgramIo.NameLimit + 1)] : names[i];
                names[i] = cut;
                if (!ReferenceEquals(cut, programNames[i]) && cut != programNames[i])
                    renamed[i] = cut;
            }

            return (rows, names, originals, renamed);
        }

        /// <summary>
        /// <c>mafft --parttree | --dpparttree | --fastaparttree</c>: <c>splittbfast</c> once, and <c>-Z</c> passes (the distances of the alignment of the previous pass) as the script runs them (<c>cycle</c> passes in all).
        /// </summary>
        private static MafftRunResult AlignPartTree(int nseq, string[] originals, string[] replaced, bool nuc, MafftOptions o, System.Text.StringBuilder log, bool epilogue, string[] partNames)
        {
            if (o.MaxIterate > 0)
                throw new ArgumentException("--parttree / --dpparttree / --fastaparttree cannot be combined with --maxiterate (MAFFT: 'Impossible' for more than 1 iteration; with 1 iteration dvtditr stops, no distance matrix 'hat2' exists).");
            if (o.DistOut)
                throw new ArgumentException("--distout with --parttree / --dpparttree / --fastaparttree: the program has no distance matrix to copy and its 'cp' fails with a message that names a random temporary directory.");
            var (cycle, groupsize) = MafftScript.PartTreeCycle(o, nseq);
            int passes = Math.Max(cycle, 1);
            string[] rows = replaced;
            int[] index = Enumerable.Range(0, nseq).ToArray();
            string? tree = null;
            string[]? labels = null;
            for (int pass = 0; pass < passes; pass++)
            {
                MafftGlobals g;
                try
                {
                    g = MafftGlobals.ParseSplittbfast(MafftScript.SplittbfastArguments(o, nseq, pass > 0, cycle, groupsize));
                }
                catch (MafftProgramFailureException e)
                {
                    throw new MafftProgramFailureException(log + e.Text, e.ExitCode); // --sp: splittbfast does not know -@
                }
                var run = new MafftRun(g);
                var res = GuardRun(run, log, () => run.Splittbfast(rows, partNames));
                log.Append(run.Log);
                rows = res.Rows;
                index = res.Order.Select(i => index[i]).ToArray();
                tree = res.Tree ?? tree;
                labels = res.GroupLabels;
            }
            string[]? coreNames = null;
            string coreText = "";
            if (o.AnySymbol && !o.CoreOut)
                rows = rows.Select((r, k) => Restore(originals[index[k]], r)).ToArray();
            if (o.CoreOut)
            {
                // setcore: the program splittbfast writes no file hat2 ("Make hat2." and an empty output)
                if (groupsize > 1 && groupsize <= nseq)
                    throw new NotSupportedException("--core with --groupsize: the sequences are not aligned with each other and setcore reads them as an alignment (undefined).");
                var crun = new MafftRun(MafftGlobals.ParseSetcore(MafftScript.SetcoreArguments(o)));
                (rows, coreNames) = crun.SetCore(rows, index.Select(k => partNames[k]).ToArray(), null);
                coreText = crun.Log.ToString();
                if (rows.Length == 0)
                    index = Array.Empty<int>();
            }
            string stderr = "";
            if (!o.Quiet)
                stderr = epilogue ? log + MafftScript.Epilogue(o, nseq) : log.ToString();
            bool labelled = groupsize > 1 && groupsize <= nseq;
            return new MafftRunResult(ToResult(rows, nuc), stderr)
            {
                CoreNames = coreNames,
                TerminalText = coreText,
                InputIndex = index,
                TreeFile = o.TreeOut || o.NodeOut ? tree : null,
                GroupLabels = labelled ? labels : null,
            };
        }

        /// <summary>The standard error that the first program of the script has written when it stops while it reads a text (<c>--text</c>): the lines of its <c>arguments()</c> (and, for <c>tbfast</c>, the two lines before the input is read).</summary>
        private static string TextReadPrefix(MafftOptions o, int nseq)
        {
            if (o.PairLibrary == MafftPairLibrary.None)
                return MafftGlobals.Parse(MafftScript.DisttbfastArguments(o, nseq)).Log.ToString();
            var (_, tb) = MafftScript.TbfastArguments(o, nseq);
            return MafftGlobals.ParseTbfast(tb).Log + "treein = " + (o.GuideTree != MafftGuideTree.Default ? "1" : "0") + "\ncompacttree = 0\n";
        }

        /// <summary>The rule of <c>--auto</c> in <c>mafft.tmpl</c> (<c>nlen</c> and <c>nseq</c> of <c>countlen</c>).</summary>
        private static MafftOptions ApplyAuto(MafftOptions o, int nseq, int nlen)
        {
            if (nlen < 3000 && nseq < 100)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.Local, MaxIterate = 1000, Retree = 1 };
            if (nlen < 1000 && nseq < 200)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.Local, MaxIterate = 2, Retree = 1 };
            if (nlen < 10000 && nseq < 500)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.None, MaxIterate = 2, Retree = 2 };
            if (nseq < 20000)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.None, MaxIterate = 0, Retree = 2 };
            if (nseq < 100000)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.None, MaxIterate = 0, Retree = 2, GuideTree = MafftGuideTree.MemSaveTree };
            if (nseq < 200000)
                return o with { Auto = false, PartTree = MafftPartTree.None, PairLibrary = MafftPairLibrary.None, MaxIterate = 0, Retree = 1, GuideTree = MafftGuideTree.MemSaveTree };
            return o with { Auto = false, PairLibrary = MafftPairLibrary.None, MaxIterate = 0, Retree = 1, PartTree = nlen < 3000 ? MafftPartTree.DynamicProgramming : MafftPartTree.KmerDistance };
        }

        internal static MultipleAlignmentResult ToResult(string[] rows, bool nucleotide)
        {
            if (rows.Length == 0)
                return MultipleAlignmentResult.Empty;
            if (rows.Select(r => r.Length).Distinct().Count() > 1)
                return new MultipleAlignmentResult(rows, "", 0); // --groupsize: the groups are not aligned with each other
            var aligned = rows.ToList();
            string consensus = BuildConsensus(aligned, aligned.Max(r => r.Length), anySymbol: true);
            var sp = nucleotide ? MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5) : MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
            string[] scored = rows.Select(r => r.ToUpperInvariant()).ToArray();
            double exact;
            try
            {
                exact = MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
            }
            catch (ArgumentException)
            {
                exact = double.NaN;
            }
            return new MultipleAlignmentResult(rows, consensus, double.IsNaN(exact) ? 0 : (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
        }

        /// <summary>
        /// The sequence as the program holds it after reading. <c>--anysymbol</c> (replaceu reads the raw lines, <c>charfilter</c>): every character except blanks, '\r' and '\n' stays (digits, '*', '#' ...; '=', '&lt;', '&gt;'
        /// stop the program); the gap '-' stays too (its positions count in the lengths and are removed by <c>gappick0</c> before the alignment). Plain: <c>onlyAlpha_*</c> and <c>kake2hiku</c> keep letters, '.', '-' and '*' (read as '-').
        /// </summary>
        private static string LoadResidues(string s, bool anySymbol)
        {
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (char ch in s)
            {
                if (ch >= 256)
                    throw new ArgumentException("MAFFT works on bytes; the sequence contains a character above U+00FF.", nameof(s));
                if (anySymbol)
                {
                    if (ch is '=' or '<' or '>')
                        throw new ArgumentException("Characters '= < >' can be used only in the title lines in the --anysymbol mode.", nameof(s));
                    if (ch is not (' ' or '\r' or '\n'))
                        sb.Append(ch);
                }
                else if (IsAsciiLetter(ch) || ch is '.' or '-')
                {
                    sb.Append(ch);
                }
                else if (ch == '*')
                {
                    sb.Append('-'); // kake2hiku: a star is read as a gap
                }
            }
            return sb.ToString();
        }

        public static bool IsNucleotide(IReadOnlyList<string> seqs, MafftOptions o) => GuessNucleotide(seqs.Select(x => LoadResidues(x, o.AnySymbol)).ToList(), o);

        /// <summary><c>restoreu -Y</c>: the residues that <c>--keeplength</c> deleted are taken out of the original letters (positions count the residues, '-' excluded).</summary>
        private static string DeleteRuns(string original, List<(int Row, int Pos, int Len)>? runs, int row)
        {
            if (runs == null || !runs.Any(r => r.Row == row))
                return original;
            var letters = NoGaps(original).ToCharArray();
            var gone = new bool[letters.Length];
            foreach (var r in runs.Where(r => r.Row == row))
                for (int k = 0; k < r.Len; k++)
                    gone[r.Pos + k] = true;
            return new string(letters.Where((_, i) => !gone[i]).ToArray());
        }

        private static string NoGaps(string s) => s.Replace("-", "", StringComparison.Ordinal);

        private static string Restore(string original, string alignedReplaced)
        {
            original = NoGaps(original);
            var sb = new System.Text.StringBuilder(alignedReplaced.Length);
            int k = 0;
            foreach (char c in alignedReplaced)
                sb.Append(c == '-' ? '-' : original[k++]);
            return sb.ToString();
        }

    }

    #endregion
}
