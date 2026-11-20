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

/// <summary>
/// Options of the MAFFT 7.505 command <c>mafft --anysymbol [--nuc | --amino] [--op x] [--ep y] [--bl n] [--retree n] [--fft | --nofft] input</c> (stage 1 of the port: the alignment of TWO sequences,
/// B09.md A2-19). The options of the command line that have no effect on the alignment of a pair are accepted and recorded (<see cref="Threads"/>, <see cref="Retree"/>, <see cref="Fft"/>,
/// <see cref="Reorder"/>); the options that need the later stages (<c>--maxiterate</c>, <c>--localpair</c>, <c>--globalpair</c>, <c>--genafpair</c>, <c>--jtt</c>, <c>--tm</c>, three and more sequences) are refused with
/// <see cref="NotSupportedException"/> naming the registered sub-item.
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

    /// <summary><c>--thread n</c> (accepted, no effect).</summary>
    public int Threads { get; init; } = 1;

    /// <summary><c>--reorder</c> (accepted; the two sequences keep their input order).</summary>
    public bool Reorder { get; init; }

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
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
    public static MultipleAlignmentResult MultipleAlignMafft(IEnumerable<string> sequences, MafftOptions? options = null)
    {
        var o = options ?? new MafftOptions();
        var list = new List<string>();
        foreach (var s in sequences ?? throw new ArgumentNullException(nameof(sequences)))
            list.Add(s ?? throw new ArgumentException("A sequence is null.", nameof(sequences)));
        // Empty input: the program accepts it ("Warning: Only 0 sequence found.", exit status 0, empty output: installed mafft 7.505), so does the port (Kalign, POA and ProbCons refuse it, as their
        // programs do: see EmptyInputMessage). One sequence: the program writes it back ("Only 1 sequence found" is a warning), here too.
        if (list.Count == 0)
            return MultipleAlignmentResult.Empty;
        if (list.Count == 1)
            return MafftEngine.SingleSequence(list[0], o);
        if (list.Count > 2)
            throw new NotSupportedException("MAFFT's progressive stage for three or more sequences (6-mer distances, guide tree, A__align) is not implemented in this library; only the alignment of two sequences is available.");
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var rows = MafftEngine.AlignPair(list[0], list[1], o);
        var aligned = rows.ToList();
        string consensus = BuildConsensus(aligned, aligned.Max(r => r.Length), anySymbol: true);
        bool nucleotide = MafftEngine.IsNucleotide(list, o);
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

    internal static class MafftEngine
    {
        /// <summary>The 26 x 26 integer scoring matrix of <c>constants()</c> (<c>n_dis</c>) in the order of the program's alphabet (<c>ARNDCQEGHILKMFPSTWYVBZX.-J</c> for amino acids, <c>agctuAGCTUnNbdhkmnrsvwyx-O</c> for nucleotides) and the integer gap opening penalty / extension.</summary>
        internal static (int[,] Matrix, int Penalty, int PenaltyEx) ScoringMatrix(bool nucleotide, MafftOptions o)
        {
            var sc = Scoring.Create(nucleotide, o);
            string alphabet = nucleotide ? MafftNucAlphabet : MafftAminoAlphabet;
            var m = new int[26, 26];
            for (int i = 0; i < 26; i++)
                for (int j = 0; j < 26; j++)
                    m[i, j] = (int)sc.Matrix[(byte)alphabet[i]][(byte)alphabet[j]];
            return (m, (int)sc.Penalty, (int)sc.PenaltyEx);
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
            string residues = LoadResidues(sequence, o.AnySymbol);
            if (!o.AnySymbol)
                residues = GuessNucleotide(new[] { residues }, o) ? residues.ToLowerInvariant() : residues.ToUpperInvariant();
            return new MultipleAlignmentResult(new[] { residues }, BuildConsensus(new List<string> { residues }, residues.Length, anySymbol: true), 0) { ExactTotalScore = 0 };
        }

        /// <summary>The alignment of a pair; rows carry the original letters (restoreu.c).</summary>
        internal static string[] AlignPair(string a, string b, MafftOptions o)
        {
            if (o.Blosum is not (30 or 45 or 50 or 62 or 80 or 90 or 100 or 0))
                throw new ArgumentException($"--bl must be 30, 45, 50, 62, 80, 90, 100 or 0, got {o.Blosum}.", nameof(o));
            if (!double.IsFinite(o.GapOpen) || !double.IsFinite(o.Offset))
                throw new ArgumentException("--op and --ep must be finite.", nameof(o));
            var originals = new[] { LoadResidues(a, o.AnySymbol), LoadResidues(b, o.AnySymbol) };
            if (originals[0].Length == 0 || originals[1].Length == 0)
                throw new ArgumentException("A sequence without residues cannot be aligned.", nameof(a));
            // replaceu: guess on the original letters (unless forced), then replace
            bool nuc1 = GuessNucleotide(originals, o);
            string[] replaced;
            if (o.AnySymbol)
                replaced = originals.Select(s => ReplaceUnusual(s, nuc1)).ToArray();
            else
                replaced = originals.Select(s => nuc1 ? s.ToLowerInvariant() : s.ToUpperInvariant()).ToArray();
            // disttbfast: the guess again, on the replaced letters
            bool nuc = GuessNucleotide(replaced, o);
            if (nuc)
                replaced = replaced.Select(s => s.ToLowerInvariant()).ToArray();
            else
                replaced = replaced.Select(s => s.ToUpperInvariant()).ToArray();
            var scoring = Scoring.Create(nuc, o);
            foreach (var s in replaced)
                foreach (char ch in s)
                    if (!scoring.Known.Contains(ch))
                        throw new ArgumentException($"Alphabet '{ch}' is unknown (MAFFT stops: try --anysymbol).", nameof(a));

            var (m1, m2) = G__align11(scoring, replaced[0], replaced[1], headgp: 0, tailgp: 0);
            // --anysymbol: restoreu puts the original letters back; the plain command writes the aligned letters as they are (nucleotides lower case, amino acids upper case)
            return o.AnySymbol ? new[] { Restore(originals[0], m1), Restore(originals[1], m2) } : new[] { m1, m2 };
        }

        /// <summary>
        /// The residues of an input sequence. <c>--anysymbol</c> (replaceu reads the raw lines): every character except blanks, '\r', '\n' and the gap '-' is a residue (digits, '*', '#' ...; '=', '&lt;', '&gt;'
        /// stop the program). Plain: <c>onlyAlpha_*</c> keeps letters, '.' and '*' (read as a gap), so only letters and '.' remain.
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
                    if (ch is not (' ' or '\r' or '\n' or '-'))
                        sb.Append(ch);
                }
                else if (IsAsciiLetter(ch) || ch == '.')
                {
                    sb.Append(ch);
                }
            }
            return sb.ToString();
        }

        internal static bool IsNucleotide(IReadOnlyList<string> seqs, MafftOptions o) => GuessNucleotide(seqs.Select(x => LoadResidues(x, o.AnySymbol)).ToList(), o);

        private static string Restore(string original, string alignedReplaced)
        {
            var sb = new System.Text.StringBuilder(alignedReplaced.Length);
            int k = 0;
            foreach (char c in alignedReplaced)
                sb.Append(c == '-' ? '-' : original[k++]);
            return sb.ToString();
        }

        // -------------------------------------------------------------------- constants.c

        private sealed class Scoring
        {
            public readonly double[][] Matrix = new double[256][];
            public readonly HashSet<char> Known = new();
            public double Penalty, PenaltyEx;

            private Scoring()
            {
                for (int i = 0; i < 256; i++)
                    Matrix[i] = new double[256];
            }

            private static int ToInt(double v) => (int)v; // C cast: truncation toward zero

            private static double Shishagonyuu(double v)
            {
                if (v > 0.0)
                    return (int)(v + 0.5);
                if (v == 0.0)
                    return 0;
                if (v < 0.0)
                    return (int)(v - 0.5);
                return 0;
            }

            /// <summary>constants.c, protein branch with scoremtx == 1 (BLOSUM): the matrix scaled to 600 / (frequency-weighted mean of the diagonal), shifted by the offset, rounded to integers.</summary>
            private static int[,] AminoMatrix(int blosum, int offset)
            {
                var tri = MafftTables.BlosumLowerTriangle[blosum];
                var freq = MafftTables.BlosumFrequencies;
                var d = new double[20, 20];
                int count = 0;
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j <= i; j++)
                        d[i, j] = d[j, i] = tri[count++];
                double average = 0.0;
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        average += d[i, j] * freq[i] * freq[j];
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        d[i, j] -= average;
                average = 0.0;
                for (int i = 0; i < 20; i++)
                    average += d[i, i] * freq[i];
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        d[i, j] *= 600.0 / average;
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        d[i, j] -= offset;
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        d[i, j] = Shishagonyuu(d[i, j]);
                var n = new int[26, 26];
                for (int i = 0; i < 20; i++)
                    for (int j = 0; j < 20; j++)
                        n[i, j] = (int)d[i, j];
                return n;
            }

            /// <summary>constants.c, nucleotide branch: the Kimura two-parameter PAM matrix of 200 steps (kappa 2, equal frequencies), log10 odds * 1000, centred, scaled to 600 / mean diagonal, shifted by the offset; ambiguity codes by <c>ambiguousscore</c>.</summary>
            private static int[,] NucleotideMatrix(string amino, int offset)
            {
                const int kimuraR = 2, pamN = 200;
                var freq = new[] { 0.25, 0.25, 0.25, 0.25 };
                var r = new double[4, 4];
                r[0, 0] = 0.0; r[0, 1] = kimuraR; r[0, 2] = 1.0; r[0, 3] = 1.0;
                r[1, 0] = kimuraR; r[1, 1] = 0.0; r[1, 2] = 1.0; r[1, 3] = 1.0;
                r[2, 0] = 1.0; r[2, 1] = 1.0; r[2, 2] = 0.0; r[2, 3] = kimuraR;
                r[3, 0] = 1.0; r[3, 1] = 1.0; r[3, 2] = kimuraR; r[3, 3] = 0.0;
                var mut = new double[4];
                double total = 0.0;
                for (int i = 0; i < 4; i++)
                {
                    double tmp = 0.0;
                    for (int j = 0; j < 4; j++)
                        tmp += r[i, j] * freq[j];
                    mut[i] = tmp;
                    total += tmp * freq[i];
                }
                var pam1 = new double[4, 4];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pam1[i, j] = i != j ? 0.01 / total * r[i, j] * freq[j] : 1.0 - 0.01 / total * mut[i];
                var pamx = new double[11, 11];
                for (int i = 0; i < 4; i++)
                    pamx[i, i] = 1.0;
                var row = new double[4];
                for (int x = 0; x < pamN; x++)
                {
                    for (int i = 0; i < 4; i++)
                    {
                        for (int k = 0; k < 4; k++)
                            row[k] = pamx[i, k];
                        for (int j = 0; j < 4; j++)
                        {
                            double s = 0.0;
                            for (int k = 0; k < 4; k++)
                                s += row[k] * pam1[k, j];
                            pamx[i, j] = s;
                        }
                    }
                }
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pamx[i, j] /= freq[j];
                for (int i = 0; i < 4; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        if (pamx[i, j] == 0.0)
                            pamx[i, j] = 0.00001;
                        pamx[i, j] = Math.Log10(pamx[i, j]) * 1000.0;
                    }
                }
                double average = 0.0;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        average += pamx[i, j] * freq[i] * freq[j];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pamx[i, j] -= average;
                average = 0.0;
                for (int i = 0; i < 4; i++)
                    average += pamx[i, i] * 1.0 / 4.0;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pamx[i, j] *= 600.0 / average;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pamx[i, j] -= offset;
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        pamx[i, j] = Shishagonyuu(pamx[i, j]);
                for (int i = 0; i < 5; i++)
                {
                    pamx[4, i] = pamx[3, i];
                    pamx[i, 4] = pamx[i, 3];
                }
                for (int i = 5; i < 10; i++)
                    for (int j = 5; j < 10; j++)
                        pamx[i, j] = pamx[i - 5, j - 5];
                var n = new int[26, 26];
                for (int i = 0; i < 10; i++)
                    for (int j = 0; j < 10; j++)
                        n[i, j] = (int)Shishagonyuu(pamx[i, j]);
                AmbiguousScore(amino, n);
                return n;
            }

            private static void AmbiguousScore(string amino, int[,] n)
            {
                int Idx(char c) => amino.LastIndexOf(c); // amino_n[c] = the last index of the letter
                int a = Idx('a'), g = Idx('g'), c = Idx('c'), t = Idx('t');
                int r = Idx('r'), y = Idx('y'), k = Idx('k'), m = Idx('m'), s = Idx('s'), w = Idx('w'), b = Idx('b'), d = Idx('d'), h = Idx('h'), v = Idx('v');
                for (int i = 0; i < 26; i++)
                {
                    n[i, r] = (int)Shishagonyuu((double)1 / 2 * (n[a, i] + n[g, i]));
                    n[i, y] = (int)Shishagonyuu((double)1 / 2 * (n[c, i] + n[t, i]));
                    n[i, k] = (int)Shishagonyuu((double)1 / 2 * (n[g, i] + n[t, i]));
                    n[i, m] = (int)Shishagonyuu((double)1 / 2 * (n[a, i] + n[c, i]));
                    n[i, s] = (int)Shishagonyuu((double)1 / 2 * (n[g, i] + n[c, i]));
                    n[i, w] = (int)Shishagonyuu((double)1 / 2 * (n[a, i] + n[t, i]));
                    n[i, b] = (int)Shishagonyuu((double)1 / 3 * (n[c, i] + n[g, i] + n[t, i]));
                    n[i, d] = (int)Shishagonyuu((double)1 / 3 * (n[a, i] + n[g, i] + n[t, i]));
                    n[i, h] = (int)Shishagonyuu((double)1 / 3 * (n[a, i] + n[c, i] + n[t, i]));
                    n[i, v] = (int)Shishagonyuu((double)1 / 3 * (n[a, i] + n[c, i] + n[g, i]));
                    n[r, i] = n[i, r];
                    n[y, i] = n[i, y];
                    n[k, i] = n[i, k];
                    n[m, i] = n[i, m];
                    n[s, i] = n[i, s];
                    n[w, i] = n[i, w];
                    n[b, i] = n[i, b];
                    n[d, i] = n[i, d];
                    n[h, i] = n[i, h];
                    n[v, i] = n[i, v];
                }
                n[r, r] = (int)Shishagonyuu((double)1 / 2 * (n[a, a] + n[g, g]));
                n[y, y] = (int)Shishagonyuu((double)1 / 2 * (n[c, c] + n[t, t]));
                n[k, k] = (int)Shishagonyuu((double)1 / 2 * (n[g, g] + n[t, t]));
                n[m, m] = (int)Shishagonyuu((double)1 / 2 * (n[a, a] + n[c, c]));
                n[s, s] = (int)Shishagonyuu((double)1 / 2 * (n[g, g] + n[c, c]));
                n[w, w] = (int)Shishagonyuu((double)1 / 2 * (n[a, a] + n[t, t]));
                n[b, b] = (int)Shishagonyuu((double)1 / 3 * (n[c, c] + n[g, g] + n[t, t]));
                n[d, d] = (int)Shishagonyuu((double)1 / 3 * (n[a, a] + n[g, g] + n[t, t]));
                n[h, h] = (int)Shishagonyuu((double)1 / 3 * (n[a, a] + n[c, c] + n[t, t]));
                n[v, v] = (int)Shishagonyuu((double)1 / 3 * (n[a, a] + n[c, c] + n[g, g]));
            }

            internal static Scoring Create(bool nuc, MafftOptions o)
            {
                var sc = new Scoring();
                // the script passes -f -op, -h -ep (and -g 0): arguments(): (int)(atof(x) * 1000 - 0.5)
                int ppenalty = ToInt(-o.GapOpen * 1000 - 0.5);
                int ppenaltyEx = ToInt(0.0 * 1000 - 0.5);
                int poffset = ToInt(-o.Offset * 1000 - 0.5);
                string amino;
                int[,] nDis;
                if (nuc)
                {
                    amino = MafftNucAlphabet;
                    sc.Penalty = ToInt(3 * 600.0 / 1000.0 * ppenalty + 0.5);
                    sc.PenaltyEx = ToInt(3 * 600.0 / 1000.0 * ppenaltyEx + 0.5);
                    int offset = ToInt(1 * 600.0 / 1000.0 * poffset + 0.5);
                    nDis = NucleotideMatrix(amino, offset);
                }
                else
                {
                    amino = MafftAminoAlphabet;
                    sc.Penalty = ToInt(600.0 / 1000.0 * ppenalty + 0.5);
                    sc.PenaltyEx = ToInt(600.0 / 1000.0 * ppenaltyEx + 0.5);
                    int offset = ToInt(600.0 / 1000.0 * poffset + 0.5);
                    nDis = AminoMatrix(o.Blosum, offset);
                }
                for (int i = 0; i < 26; i++)
                {
                    sc.Known.Add(amino[i]);
                    for (int j = 0; j < 26; j++)
                        sc.Matrix[(byte)amino[i]][(byte)amino[j]] = (double)nDis[i, j] * 1.0; // consweight_multi = 1.0
                }
                return sc;
            }
        }

        // -------------------------------------------------------------------- Galign11.c (trywarp = 0: -Q 100 >= 10)

        private static (string, string) G__align11(Scoring sc, string seq1, string seq2, int headgp, int tailgp)
        {
            int lgth1 = seq1.Length, lgth2 = seq2.Length;
            double fpenalty = sc.Penalty, fpenaltyEx = sc.PenaltyEx;
            const double termGapFac = 0.0, termGapFacEx = 0.0;
            RequireCells(((long)lgth1 + 1) * ((long)lgth2 + 1), MafftMaxPairCells, nameof(MafftMaxPairCells), "MAFFT", "keeps an int trace-back matrix ijp[len1 + 1][len2 + 1] (G__align11)",
                string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{lgth1:N0} x {lgth2:N0} residues = {((long)lgth1 + 1) * ((long)lgth2 + 1):N0} cells"), "The program has no limit and is killed when the allocation (4 bytes per cell plus the row arrays) fails", "use the progressive / FFT variants of MAFFT (not ported) or Kalign for sequences of this length");
            var m = new double[lgth2 + 2];
            var mp = new int[lgth2 + 2];
            var w1 = new double[lgth2 + 2];
            var w2 = new double[lgth2 + 2];
            var initverticalw = new double[lgth1 + 2];
            var lastverticalw = new double[lgth1 + 2];
            var ijp = new int[lgth1 + 1][];
            for (int i = 0; i <= lgth1; i++)
                ijp[i] = new int[lgth2 + 1];

            double[] currentw = w1, previousw = w2;
            for (int j = 0; j < lgth1; j++)
                initverticalw[j] = sc.Matrix[(byte)seq2[0]][(byte)seq1[j]];
            for (int j = 0; j < lgth2; j++)
                currentw[j] = sc.Matrix[(byte)seq1[0]][(byte)seq2[j]];
            if (headgp == 1)
            {
                for (int i = 1; i < lgth1 + 1; i++)
                    initverticalw[i] += fpenalty;
                for (int j = 1; j < lgth2 + 1; j++)
                    currentw[j] += fpenalty;
            }
            else
            {
                for (int i = 1; i < lgth1 + 1; i++)
                {
                    initverticalw[i] += fpenalty * termGapFac;
                    initverticalw[i] += fpenaltyEx * i * termGapFacEx;
                }
                for (int j = 1; j < lgth2 + 1; j++)
                {
                    currentw[j] += fpenalty * termGapFac;
                    currentw[j] += fpenaltyEx * j * termGapFacEx;
                }
            }
            for (int j = 1; j < lgth2 + 1; ++j)
            {
                m[j] = currentw[j - 1];
                mp[j] = 0;
            }
            lastverticalw[0] = lgth2 == 0 ? 0.0 : currentw[lgth2 - 1];
            int lasti = tailgp != 0 ? lgth1 + 1 : lgth1;
            int lastj = lgth2 + 1;
            double mi;
            int mpi;
            for (int i = 1; i < lasti; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                (previousw, currentw) = (currentw, previousw);
                previousw[0] = initverticalw[i - 1];
                for (int j = 0; j < lgth2; j++)
                    currentw[j] = sc.Matrix[(byte)seq1[i]][(byte)seq2[j]];
                currentw[0] = initverticalw[i];
                mi = previousw[0];
                mpi = 0;
                double fpenaltyExI = i < lgth1 ? fpenaltyEx : 0.0;
                var ijpRow = ijp[i];
                for (int j = 1; j < lastj; j++)
                {
                    double wm = previousw[j - 1];
                    ijpRow[j] = 0;
                    double g;
                    if ((g = mi + fpenalty) > wm)
                    {
                        wm = g;
                        ijpRow[j] = -(j - mpi);
                    }
                    if ((g = previousw[j - 1]) >= mi)
                    {
                        mi = g;
                        mpi = j - 1;
                    }
                    mi += fpenaltyExI;
                    if ((g = m[j] + fpenalty) > wm)
                    {
                        wm = g;
                        ijpRow[j] = +(i - mp[j]);
                    }
                    if ((g = previousw[j - 1]) >= m[j])
                    {
                        m[j] = g;
                        mp[j] = i - 1;
                    }
                    if (j < lgth2)
                        m[j] += fpenaltyEx;
                    currentw[j] += wm;
                }
                lastverticalw[i] = currentw[lgth2 - 1];
            }
            return Atracking(currentw, lastverticalw, seq1, seq2, ijp, tailgp);
        }

        private static (string, string) Atracking(double[] lasthorizontalw, double[] lastverticalw, string seq1, string seq2, int[][] ijp, int tailgp)
        {
            int lgth1 = seq1.Length, lgth2 = seq2.Length;
            for (int i = 0; i < lgth1 + 1; i++)
                ijp[i][0] = i + 1;
            for (int j = 0; j < lgth2 + 1; j++)
                ijp[0][j] = -(j + 1);
            if (tailgp != 1)
            {
                double wm = lasthorizontalw[lgth2 - 1] - 1.0;
                double g;
                for (int j = lgth2 - 2; j >= 0; j--)
                {
                    if ((g = lasthorizontalw[j] + 0.0) > wm)
                    {
                        wm = g;
                        ijp[lgth1][lgth2] = -(lgth2 - j);
                    }
                }
                for (int i = lgth1 - 2; i >= 0; i--)
                {
                    if ((g = lastverticalw[i] + 0.0) > wm)
                    {
                        wm = g;
                        ijp[lgth1][lgth2] = +(lgth1 - i);
                    }
                }
                if (lasthorizontalw[lgth2 - 1] > wm)
                    ijp[lgth1][lgth2] = 0;
            }
            var a = new char[lgth1 + lgth2 + 1];
            var b = new char[lgth1 + lgth2 + 1];
            int pa = lgth1 + lgth2, pb = lgth1 + lgth2;
            int iin = lgth1, jin = lgth2;
            int limk = lgth1 + lgth2 + 1;
            int k = 0;
            while (k < limk)
            {
                int ifi, jfi;
                if (ijp[iin][jin] < 0)
                {
                    ifi = iin - 1;
                    jfi = jin + ijp[iin][jin];
                }
                else if (ijp[iin][jin] > 0)
                {
                    ifi = iin - ijp[iin][jin];
                    jfi = jin - 1;
                }
                else
                {
                    ifi = iin - 1;
                    jfi = jin - 1;
                }
                int l = iin - ifi;
                while (--l > 0)
                {
                    a[--pa] = seq1[ifi + l];
                    b[--pb] = '-';
                    k++;
                }
                l = jin - jfi;
                while (--l > 0)
                {
                    a[--pa] = '-';
                    b[--pb] = seq2[jfi + l];
                    k++;
                }
                if (iin <= 0 || jin <= 0)
                    break;
                a[--pa] = seq1[ifi];
                b[--pb] = seq2[jfi];
                k++;
                iin = ifi;
                jin = jfi;
                k++;
            }
            return (new string(a, pa, lgth1 + lgth2 - pa), new string(b, pb, lgth1 + lgth2 - pb));
        }
    }

    #endregion
}
