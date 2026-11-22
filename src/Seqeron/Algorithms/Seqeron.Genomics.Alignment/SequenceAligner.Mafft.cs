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

    /// <summary><c>--weighti x</c>: the weight of the pair library against the substitution scores (default 2.7). 0 switches the constraints off, which this port does not support.</summary>
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

    /// <summary><c>--oldgenafpair</c>: the parameters of E-INS-i before version 7.243 (the offset and the extension of the pair stage as given, the alignment score of the pair).</summary>
    public bool OldGenafParameters { get; init; }

    /// <summary><c>--reorder</c> (accepted; the two sequences keep their input order).</summary>
    public bool Reorder { get; init; }

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
}

/// <summary>The result of a run of <c>mafft</c>: the alignment and the text that the program writes to standard error (progress, strategy; empty with <c>--quiet</c>).</summary>
/// <param name="Alignment">Rows in the order of the output (the input order, or the order of the guide tree with <c>--reorder</c>).</param>
/// <param name="StandardError">The messages of the binaries and of the script, byte for byte.</param>
public sealed record MafftRunResult(MultipleAlignmentResult Alignment, string StandardError);

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

        private static bool IsAsciiLetter(char c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

        /// <summary>replaceu.c: letters outside the usual alphabet are replaced (X for amino acids, n for nucleotides), the rest is upper- / lower-cased.</summary>
        private static string ReplaceUnusual(string s, bool nucleotide)
        {
            const string usualP = "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";
            const string usualN = "ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (char ch in s)
            {
                if (ch >= 128)
                    throw new ArgumentException("MAFFT works on ASCII letters; the sequence contains a character above U+007F.", nameof(s));
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
        internal static MafftRunResult AlignMany(IReadOnlyList<string> list, MafftOptions o)
        {
            if (o.Blosum is not (30 or 45 or 50 or 62 or 80 or 90 or 100 or 0))
                throw new ArgumentException($"--bl must be 30, 45, 50, 62, 80, 90, 100 or 0, got {o.Blosum}.", nameof(o));
            if (!double.IsFinite(o.GapOpen) || !double.IsFinite(o.Offset) || !double.IsFinite(o.GapExtension))
                throw new ArgumentException("--op, --ep and --exp must be finite.", nameof(o));
            if (o.Retree == 0)
                throw new NotSupportedException("--retree 0 (the guide tree only) of MAFFT is not available in this port.");
            if (o.PairLibrary != MafftPairLibrary.None && !(o.WeightI >= 0.0001))
                throw new NotSupportedException("--weighti 0 (the pair library without constraints) of MAFFT is not available in this port.");
            var originals = list.Select(x => LoadResidues(x, o.AnySymbol)).ToArray();
            foreach (var s in originals)
                if (NoGaps(s).Length == 0)
                    throw new ArgumentException("A sequence without residues cannot be aligned.", nameof(list));
            bool nuc1 = GuessNucleotide(originals, o);
            var log = new System.Text.StringBuilder();
            if (o.AnySymbol)
                log.Append("inputfile = orig\n").Append(originals.Length).Append(" x ").Append(originals.Max(x => x.Length)).Append(" - ").Append(originals.Min(x => x.Length)).Append(' ').Append(nuc1 ? 'd' : 'p').Append('\n');
            string[] replaced = o.AnySymbol
                ? originals.Select(s => ReplaceUnusual(s, nuc1)).ToArray()
                : originals.Select(s => nuc1 ? s.ToLowerInvariant() : s.ToUpperInvariant()).ToArray();
            bool nuc = GuessNucleotide(replaced, o);
            replaced = nuc ? replaced.Select(s => s.ToLowerInvariant()).ToArray() : replaced.Select(s => s.ToUpperInvariant()).ToArray();
            MafftRun.Result result;
            double[][]? pairHat2 = null;
            MafftLocalHom[][]? pairHat3 = null;
            if (o.PairLibrary != MafftPairLibrary.None)
            {
                var (pairArgs, tbArgs) = MafftScript.TbfastArguments(o, list.Count);
                var tg = MafftGlobals.ParseTbfast(tbArgs);
                if (tg.Dorp == '\0')
                    tg.Dorp = nuc ? 'd' : 'p';
                var trun = new MafftRun(tg);
                result = trun.Tbfast(replaced, MafftGlobals.ParsePairlocalalign(pairArgs));
                log.Append(trun.Log);
                pairHat2 = trun.Hat2Table;
                pairHat3 = trun.Hat3Table;
            }
            else
            {
                var g = MafftGlobals.Parse(MafftScript.DisttbfastArguments(o, list.Count));
                if (g.Dorp == '\0')
                    g.Dorp = nuc ? 'd' : 'p';
                var run = new MafftRun(g);
                result = run.Disttbfast(replaced);
                log.Append(run.Log);
            }
            if (MafftScript.Iterations(o) > 0)
            {
                int threads = o.Threads < 0 ? Environment.ProcessorCount : o.Threads;
                int threadsIt = o.ThreadsIterate ?? Math.Min(threads, 8);
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
                    hat2 = run2.Dndpre(result.Rows);
                    log.Append(run2.Log);
                }
                var run3 = new MafftRun(MafftGlobals.ParseDvtditr(MafftScript.DvtditrArguments(o)));
                result = run3.Dvtditr(result.Rows, hat2, pairHat3);
                log.Append(run3.Log);
            }
            string[] rows = o.AnySymbol
                ? result.Rows.Select((r, i) => Restore(originals[i], r)).ToArray()
                : result.Rows;
            if (o.Reorder)
                rows = result.Order.Select(i => rows[i]).ToArray();
            string stderr = o.Quiet ? "" : log.ToString() + MafftScript.Epilogue(o, list.Count);
            return new MafftRunResult(ToResult(rows, nuc), stderr);
        }

        internal static MultipleAlignmentResult ToResult(string[] rows, bool nucleotide)
        {
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
                if (ch >= 128)
                    throw new ArgumentException("MAFFT works on ASCII letters; the sequence contains a character above U+007F.", nameof(s));
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

        internal static bool IsNucleotide(IReadOnlyList<string> seqs, MafftOptions o) => GuessNucleotide(seqs.Select(x => LoadResidues(x, o.AnySymbol)).ToList(), o);

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
