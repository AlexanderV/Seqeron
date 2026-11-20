using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Alignment;

/// <summary>What a run of <see cref="PoaProgram.Run"/> produced.</summary>
/// <param name="ExitCode">0 = success, 1 = the program stops with an error, 255 = called without arguments (the usage screen, <c>exit(-1)</c>).</param>
/// <param name="Stdout">Standard output: the <c>-v</c> gap tables, the <c>-printmatrix</c> table, the permutation line of <c>-preserve_seqorder</c>.</param>
/// <param name="Stderr">The error messages and the <c>...Read N sequences...</c> / <c>...Wrote N sequences...</c> lines (the alignment progress lines of the program, which print uninitialised <c>%ld</c> arguments, are not reproduced).</param>
/// <param name="Files">The files written: the <c>-pir</c>, <c>-clustal</c> and <c>-po</c> outputs.</param>
public sealed record PoaProgramResult(int ExitCode, string Stdout, string Stderr, IReadOnlyDictionary<string, string> Files);

/// <summary>
/// POA 2.0 as a command line (<c>poa -read_fasta in.fa [-do_global] [-do_progressive] [-read_pairscores F] [-fuse_all] [-hb] [-hbmin X] [-pir F] [-clustal F] [-po F]
/// [-preserve_seqorder] [-best] [-tolower] [-toupper] [-printmatrix LETTERS] [-v] MATRIXFILE</c>): the argument loop of <c>main.c</c> (every argument is tested against every flag, unknown flags
/// are ignored, the first argument that does not start with a dash is the matrix file), the FASTA reader of <c>fasta_format.c</c> and the writers of <c>lpo_format.c</c>. Input and output files are
/// text given as a dictionary.
/// </summary>
/// <remarks>
/// Not ported (accepted by the program, unsupported here): <c>-read_msa</c>, <c>-read_msa2</c>, <c>-read_msa_list</c> (MSA / PO file input) and <c>-subset</c> / <c>-remove</c> (they raise
/// <see cref="NotSupportedException"/>).
/// </remarks>
public static class PoaProgram
{
    private const string Usage =
        "\nUsage: poa [OPTIONS] MATRIXFILE\nAlign a set of sequences or alignments using the scores in MATRIXFILE.\nExample: poa -read_fasta multidom.seq -clustal m.aln blosum80.mat\n\n" +
        "INPUT:\n  -read_fasta FILE       Read in FASTA sequence file.\n  -read_msa FILE         Read in MSA alignment file.\n  -read_msa2 FILE        Read in second MSA file. \n" +
        "  -subset FILE           Filter MSA to include list of seqs in file.\n  -subset2 FILE          Filter second MSA to include list of seqs in file.\n" +
        "  -remove FILE           Filter MSA to exclude list of seqs in file.\n  -remove2 FILE          Filter second MSA to exclude list of seqs in file.\n" +
        "  -read_msa_list FILE    Read an MSA from each filename listed in file.\n  -tolower               Force FASTA/MSA sequences to lowercase\n" +
        "                           (nucleotides in our matrix files)\n  -toupper               Force FASTA/MSA sequences to UPPERCASE\n" +
        "                           (amino acids in our matrix files)\n\nALIGNMENT:\n  -do_global             Do global alignment.\n" +
        "  -do_progressive        Perform progressive alignment using a guide tree\n                           built by neighbor joining from a set of\n                           sequence-sequence similarity scores.\n" +
        "  -read_pairscores FILE  Read tab-delimited file of similarity scores.\n                           (If not provided, scores are constructed\n                           using pairwise sequence alignment.)\n" +
        "  -fuse_all              Fuse identical letters on align rings.\n\nANALYSIS:\n  -hb                    Perform heaviest bundling to generate consensi.\n" +
        "  -hbmin VALUE           Include in heaviest bundle sequences with\n                           percent ID (as a fraction) >= value.\n\nOUTPUT:\n" +
        "  -pir FILE              Write out MSA in PIR format.\n  -clustal FILE          Write out MSA in CLUSTAL format.\n  -po FILE               Write out MSA in PO format.\n" +
        "  -preserve_seqorder     Write out MSA with sequences in their input order.\n  -printmatrix LETTERS   Print score matrix to stdout.\n" +
        "  -best                  Restrict MSA output to heaviest bundles (PIR only).\n  -v                     Run in verbose mode (e.g. output gap penalties).\n\n" +
        "  NOTE:  One of the -read_fasta, -read_msa, or -read_msa_list arguments\n         must be used, since a sequence or alignment file is required.\n\n" +
        "For more information, see http://www.bioinformatics.ucla.edu/poa.\n\n";

    /// <summary>The usage screen (what the program writes to standard error when called without arguments).</summary>
    public static string UsageText => Usage;

    /// <summary><c>read_fasta</c>: the records of a FASTA text (names are the first word of the header, titles the rest, <c>untitled</c> when there is none; <c>#</c> lines are comments; lines starting with <c>*</c> are ignored).</summary>
    public static IReadOnlyList<PoaSequence> ReadFasta(string text, out string? comment)
    {
        ArgumentNullException.ThrowIfNull(text);
        comment = null;
        var seqs = new List<PoaSequence>();
        string name = string.Empty, title = string.Empty;
        var tmp = new StringBuilder();
        int pos = 0;
        string? cmt = null;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            string line = nl < 0 ? text[pos..] : text[pos..nl];
            pos = nl < 0 ? text.Length : nl + 1;
            if (line.Length > 0 && line[0] == '#')
            {
                cmt = line[1..];
            }
            else if (line.Length > 0 && line[0] == '>')
            {
                if (name.Length > 0 && tmp.Length > 0) seqs.Add(new PoaSequence(name, title, tmp.ToString()));
                name = string.Empty;
                var m = Regex.Match(line[1..], @"^\s*(\S+)\s*(.*)$", RegexOptions.Singleline);
                if (m.Success)
                {
                    name = m.Groups[1].Value;
                    title = m.Groups[2].Value.Length > 0 ? m.Groups[2].Value : "untitled";
                }
                else
                {
                    title = "untitled";
                }

                tmp.Clear();
            }
            else if (line.Length > 0 && line[0] == '*')
            {
                // ignored
            }
            else if (name.Length > 0)
            {
                tmp.Append(line);
            }

            if (pos < text.Length && text[pos] == '#' && seqs.Count > 0) break; // unigene cluster terminator
        }

        if (name.Length > 0 && tmp.Length > 0) seqs.Add(new PoaSequence(name, title, tmp.ToString()));
        comment = cmt;
        return seqs;
    }

    /// <summary>Runs the program on text "files". See <see cref="PoaProgramResult"/>.</summary>
    /// <param name="arguments">The command-line arguments (without the program name).</param>
    /// <param name="files">Input files (name -> text): the FASTA file, the matrix file, the pair-score file.</param>
    /// <exception cref="NotSupportedException">An MSA input option (unsupported; use FASTA input).</exception>
    public static PoaProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        files ??= new Dictionary<string, string>();
        var err = new StringBuilder();
        var outFiles = new Dictionary<string, string>();
        if (arguments.Count < 1)
            return new PoaProgramResult(255, string.Empty, Usage, outFiles);

        string? fastaFile = null, matrixFile = null, pairFile = null, pirOut = null, clustalOut = null, poOut = null, printMatrix = null, hbmin = null;
        string? msa = null, msa2 = null, msaList = null, subset = null, subset2 = null, remove = null, remove2 = null;
        bool lower = false, upper = false, verbose = false, best = false, hb = false, fuseAll = false, global = false, progressive = false, preserve = false;
        for (int i = 0; i < arguments.Count; i++)
        {
            string a = arguments[i];
            if (a == "-tolower") lower = true;
            if (a == "-toupper") upper = true;
            if (a == "-v") verbose = true;
            if (a == "-best") best = true;
            if (a == "-hb") hb = true;
            if (a == "-preserve_seqorder") preserve = true;
            if (a == "-fuse_all") fuseAll = true;
            if (a == "-do_global") global = true;
            if (a == "-do_progressive") progressive = true;
            if (TakeValue(arguments, ref i, "-printmatrix", ref printMatrix) || TakeValue(arguments, ref i, "-read_msa", ref msa) || TakeValue(arguments, ref i, "-read_msa2", ref msa2)
                || TakeValue(arguments, ref i, "-read_msa_list", ref msaList) || TakeValue(arguments, ref i, "-pir", ref pirOut) || TakeValue(arguments, ref i, "-clustal", ref clustalOut)
                || TakeValue(arguments, ref i, "-po", ref poOut) || TakeValue(arguments, ref i, "-hbmin", ref hbmin) || TakeValue(arguments, ref i, "-read_pairscores", ref pairFile)
                || TakeValue(arguments, ref i, "-subset", ref subset) || TakeValue(arguments, ref i, "-subset2", ref subset2) || TakeValue(arguments, ref i, "-remove", ref remove)
                || TakeValue(arguments, ref i, "-remove2", ref remove2) || TakeValue(arguments, ref i, "-read_fasta", ref fastaFile))
                continue;
            if (a.Length > 0 && a[0] != '-' && matrixFile == null) matrixFile = a;
        }

        if (msaList != null && (msa != null || msa2 != null))
        {
            err.Append("Error: The -read_po_list and -read_po flags cannot be used at the same time.\nExiting.\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        if (((subset != null || remove != null) && msa == null) || ((subset2 != null || remove2 != null) && msa2 == null))
        {
            err.Append("Error: Each -subset/-remove flag must have a corresponding -read_po flag.\nExiting.\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        if ((subset != null && remove != null) || (subset2 != null && remove2 != null))
        {
            err.Append("Error: The -subset and -remove flags cannot be used at the same time.\nExiting.\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        if (msa != null || msa2 != null || msaList != null)
            throw new NotSupportedException("-read_msa / -read_msa2 / -read_msa_list (MSA and PO file input) are accepted by the program but are not supported by this library; use -read_fasta.");

        float hbThreshold = 0.9f;
        if (hbmin != null) hbThreshold = (float)Atof(hbmin);
        PoaScoreMatrix matrix;
        try
        {
            if (matrixFile == null || !files.TryGetValue(matrixFile, out var mtext)) throw new ArgumentException("no matrix");
            matrix = PoaScoreMatrix.Parse(mtext);
        }
        catch (ArgumentException)
        {
            err.Append("Error reading matrix file ").Append(matrixFile ?? "because none specified").Append(".\nExiting\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        IReadOnlyList<PoaSequence> seqs;
        string? comment;
        if (fastaFile != null)
        {
            if (!files.TryGetValue(fastaFile, out var ftext))
            {
                err.Append("Couldn't open sequence file ").Append(fastaFile).Append(".\nExiting\n");
                return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
            }

            seqs = ReadFasta(ftext, out comment);
            if (seqs.Count == 0)
            {
                err.Append("Error reading sequence file ").Append(fastaFile).Append(".\nExiting\n");
                return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
            }

            err.Append("...Read ").Append(seqs.Count.ToString(CultureInfo.InvariantCulture)).Append(" sequences from sequence file ").Append(fastaFile).Append("...\n");
        }
        else
        {
            err.Append("No input sequences were provided; use one of the -read_ flags.\nExiting.\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        string? pairText = null;
        if (pairFile != null && !files.TryGetValue(pairFile, out pairText))
        {
            err.Append("Error reading pair score file ").Append(pairFile).Append(".\nExiting\n");
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        var options = new PoaOptions
        {
            Global = global,
            Progressive = progressive,
            PairScores = pairText,
            FuseAll = fuseAll,
            HeaviestBundle = hb,
            HeaviestBundleMinimum = hbThreshold,
            Best = best,
            PreserveSequenceOrder = preserve,
            Case = CaseOf(lower, upper),
            Verbose = verbose,
            PrintMatrix = printMatrix,
            Comment = comment,
        };
        PoaRunResult run;
        try
        {
            run = SequenceAligner.PoaRun(seqs, matrix, options);
        }
        catch (ArgumentException e)
        {
            err.Append(e.Message).Append('\n');
            return new PoaProgramResult(1, string.Empty, err.ToString(), outFiles);
        }

        int nsource = run.Names.Length - (hb ? run.Order.Count(o => o < 0) : 0);
        if (poOut != null)
        {
            outFiles[poOut] = run.Po;
            err.Append("...Wrote ").Append(nsource.ToString(CultureInfo.InvariantCulture)).Append(" sequences to PO file ").Append(poOut).Append("...\n");
        }

        if (pirOut != null)
        {
            outFiles[pirOut] = run.Pir;
            err.Append("...Wrote ").Append(nsource.ToString(CultureInfo.InvariantCulture)).Append(" sequences to FASTA-PIR file ").Append(pirOut).Append("...\n");
        }

        if (clustalOut != null)
        {
            outFiles[clustalOut] = run.Clustal;
            err.Append("...Wrote ").Append(nsource.ToString(CultureInfo.InvariantCulture)).Append(" sequences to CLUSTAL file ").Append(clustalOut).Append("...\n");
        }

        return new PoaProgramResult(0, run.Stdout, err.ToString(), outFiles);
    }

    private static PoaCase CaseOf(bool lower, bool upper)
    {
        if (upper) return PoaCase.Upper;
        return lower ? PoaCase.Lower : PoaCase.Keep;
    }

    private static bool TakeValue(IReadOnlyList<string> args, ref int i, string flag, ref string? target)
    {
        if (args[i] != flag) return false;
        target = i + 1 < args.Count ? args[++i] : null;
        return true;
    }

    private static double Atof(string s)
    {
        var m = Regex.Match(s, @"^\s*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?");
        return m.Success ? double.Parse(m.Value.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture) : 0.0;
    }
}
