using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The meaning of a Kalign 3.3.5 command line (<see cref="KalignProgram.Parse"/>; <c>src/run_kalign.c</c> <c>main</c>): the alignment parameters and the options that only concern input and output.
/// </summary>
/// <param name="Options">Alignment parameters (<c>--type</c>, <c>--gpo</c>, <c>--gpe</c>, <c>--tgpe</c>, <c>-n</c>).</param>
/// <param name="InputFiles">Input file names: <c>-i / --input / --infile / --in</c> first, then the non-option arguments (the program also reads standard input when it is not a terminal).</param>
/// <param name="OutFile"><c>-o / --output / --outfile / --out</c>; null = standard output.</param>
/// <param name="Format"><c>-f / --format</c> (<c>fasta</c> default; <c>msf</c>, <c>clu</c> are accepted by the program).</param>
/// <param name="Quiet"><c>-q / --quiet</c> (no banner / log lines).</param>
/// <param name="Help"><c>-h / --help</c>.</param>
/// <param name="Version"><c>-v / -V / --version</c>.</param>
/// <param name="ShowWarranty"><c>--showw</c>.</param>
/// <param name="ParamSet"><c>--set</c> (stored by the program, never used).</param>
public sealed record KalignCommand(KalignOptions Options, IReadOnlyList<string> InputFiles, string? OutFile, string? Format, bool Quiet, bool Help, bool Version, bool ShowWarranty, int? ParamSet);

/// <summary>What a run of <see cref="KalignProgram.Run"/> produced.</summary>
/// <param name="ExitCode">0 = success, 1 = the program stopped with an error.</param>
/// <param name="Stdout">What the program writes to standard output (the alignment when no output file is named, the version / help text).</param>
/// <param name="Files">Files the program wrote (<c>-o</c>): name -> text.</param>
/// <param name="Messages">What the program prints to the error stream (errors, warnings).</param>
public sealed record KalignProgramResult(int ExitCode, string Stdout, IReadOnlyDictionary<string, string> Files, IReadOnlyList<string> Messages);

/// <summary>
/// Kalign 3.3.5 as a command line (<c>kalign -i in.fa -o out.fa --type dna --gpo 5 --gpe 2 --tgpe 1 -n 4 -f fasta</c>): the argument parser of <c>run_kalign.c</c> (<c>getopt_long_only</c>: long
/// options with one or two dashes and unambiguous abbreviations, <c>--name=value</c> / <c>--name value</c>, short options <c>-i -o -f -n -h -q -v -V</c>), the type word of <c>set_aln_type</c>
/// (<c>strstr</c> in the order <c>rna</c>, <c>dna</c>, <c>internal</c>, <c>protein</c>, <c>divergent</c>: <c>--type internal</c> contains "rna" and therefore means RNA), the format check of
/// <c>parameters.c</c> and the FASTA reader / writer of <c>msa_io.c</c> (60 columns per line). Input and output are text, the "files" are given as a dictionary.
/// </summary>
/// <remarks>
/// Accepted and documented but without effect on the alignment: <c>-n / --nthreads</c> beyond the one-thread / many-threads distinction of <see cref="KalignOptions.Threads"/>, <c>--set</c>, <c>-q</c>.
/// <c>--format msf</c> and <c>--format clu</c> (output) and MSF / Clustal inputs are accepted by the program but not yet written / read here: <see cref="NotSupportedException"/> (registered, B09.md A2-20b).
/// </remarks>
public static class KalignProgram
{
    private const string Version = "kalign 3.3.5";

    private static readonly (string Name, bool HasArg, char Val)[] LongOptions =
    {
        ("showw", false, '\x05'), ("set", true, '\x01'), ("format", true, 'f'), ("type", true, '\x0d'), ("gpo", true, '\x06'), ("gpe", true, '\x07'), ("tgpe", true, '\x08'),
        ("nthreads", true, 'n'), ("input", true, 'i'), ("infile", true, 'i'), ("in", true, 'i'), ("output", true, 'o'), ("outfile", true, 'o'), ("out", true, 'o'),
        ("help", false, 'h'), ("version", false, 'v'), ("quiet", false, 'q'),
    };

    /// <summary>C <c>atof</c>: the longest numeric prefix (decimal, optional exponent), 0 when there is none.</summary>
    private static double Atof(string s)
    {
        var m = Regex.Match(s, @"^\s*[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?");
        return m.Success ? double.Parse(m.Value.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture) : 0.0;
    }

    /// <summary>C <c>atoi</c>.</summary>
    private static int Atoi(string s)
    {
        var m = Regex.Match(s, @"^\s*[+-]?\d+");
        return m.Success && long.TryParse(m.Value.Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long v) ? (int)Math.Clamp(v, int.MinValue, int.MaxValue) : 0;
    }

    /// <summary>The alignment type of the <c>--type</c> word (<c>set_aln_type</c>).</summary>
    /// <exception cref="ArgumentException">The word contains none of the five type names (the program stops: "In %s not recognized.").</exception>
    public static KalignType ParseType(string word)
    {
        ArgumentNullException.ThrowIfNull(word);
        if (word.Contains("rna", StringComparison.Ordinal)) return KalignType.Rna;
        if (word.Contains("dna", StringComparison.Ordinal)) return KalignType.Dna;
        if (word.Contains("internal", StringComparison.Ordinal)) return KalignType.Rna; // sic: "internal" contains "rna", so the program's DNA_INTERNAL branch is unreachable
        if (word.Contains("protein", StringComparison.Ordinal)) return KalignType.Protein;
        if (word.Contains("divergent", StringComparison.Ordinal)) return KalignType.ProteinDivergent;
        throw new ArgumentException($"In {word} not recognized.", nameof(word));
    }

    /// <summary>Parses the arguments like the program's <c>getopt_long_only</c> loop and checks what <c>main</c> checks.</summary>
    /// <exception cref="ArgumentException">An unknown / ambiguous option, a missing option argument, a type word or format that is not recognised, fewer than one thread.</exception>
    public static KalignCommand Parse(IReadOnlyList<string> arguments)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        var o = new KalignOptions();
        var inputs = new List<string>();
        var rest = new List<string>();
        string? inFile = null, outFile = null, format = null, typeWord = null;
        bool quiet = false, help = false, version = false, showw = false;
        int? paramSet = null;
        int a = 0;
        while (a < arguments.Count)
        {
            string arg = arguments[a++];
            if (arg == "--")
            {
                rest.AddRange(arguments.Skip(a));
                break;
            }
            if (arg.Length < 2 || arg[0] != '-')
            {
                rest.Add(arg);
                continue;
            }
            string body = arg.TrimStart('-');
            string? inlineValue = null;
            int eq = body.IndexOf('=');
            if (eq >= 0)
            {
                inlineValue = body[(eq + 1)..];
                body = body[..eq];
            }
            char val;
            bool hasArg;
            var exact = LongOptions.Where(l => l.Name == body).ToList();
            var candidates = exact.Count > 0 ? exact : LongOptions.Where(l => l.Name.StartsWith(body, StringComparison.Ordinal)).ToList();
            if (candidates.Count > 0)
            {
                if (candidates.Select(c => (c.HasArg, c.Val)).Distinct().Count() > 1)
                    throw new ArgumentException($"option '{arg}' is ambiguous");
                (_, hasArg, val) = candidates[0];
            }
            else if (arg[0] == '-' && arg[1] != '-' && "iofnhqvV".Contains(body[0]))
            {
                val = body[0];
                hasArg = "iofn".Contains(val);
                if (hasArg && body.Length > 1 && inlineValue == null)
                    inlineValue = body[1..]; // -n4
                else if (!hasArg && body.Length > 1)
                    throw new ArgumentException($"invalid option -- '{body[1]}'");
            }
            else
                throw new ArgumentException($"unrecognized option '{arg}'");
            string? value = null;
            if (hasArg)
            {
                if (inlineValue != null)
                    value = inlineValue;
                else if (a < arguments.Count)
                    value = arguments[a++];
                else
                    throw new ArgumentException($"option '{arg}' requires an argument");
            }
            switch (val)
            {
                case '\x0d': typeWord = value; break;
                case '\x05': showw = true; break;
                case '\x01': paramSet = Atoi(value!); break;
                case 'f': format = value; break;
                case 'n': o = o with { Threads = Atoi(value!) }; break;
                case 'q': quiet = true; break;
                case '\x06': o = o with { GapOpen = (float)Atof(value!) }; break;
                case '\x07': o = o with { GapExtend = (float)Atof(value!) }; break;
                case '\x08': o = o with { TerminalGapExtend = (float)Atof(value!) }; break;
                case 'h': help = true; break;
                case 'v':
                case 'V': version = true; break;
                case 'i': inFile = value; break;
                case 'o': outFile = value; break;
                default: throw new InvalidOperationException();
            }
        }
        if (inFile != null)
            inputs.Add(inFile);
        inputs.AddRange(rest);
        if (!version && !showw && !help)
        {
            if (o.Threads < 1)
                throw new ArgumentException("Number of threads has to be >= 1.");
            CheckFormat(format);
            if (typeWord != null)
                o = o with { Type = ParseType(typeWord) };
        }
        return new KalignCommand(o, inputs, outFile, format, quiet, help, version, showw, paramSet);
    }

    private static void CheckFormat(string? format)
    {
        if (format == null)
            return;
        if (!(format.Contains("msf", StringComparison.Ordinal) || format.Contains("clu", StringComparison.Ordinal) || format.Contains("fasta", StringComparison.Ordinal) || format.Contains("fa", StringComparison.Ordinal)))
            throw new ArgumentException($"Format {format} not recognized.");
    }

    /// <summary>
    /// Runs the program: reads the input files (name -> FASTA text; a missing name is an error; <paramref name="stdin"/> is read when no file is named), aligns with <see cref="SequenceAligner.KalignRun"/>
    /// and writes the alignment as FASTA (60 residues per line, <c>write_msa_fasta</c>) to standard output or to the <c>-o</c> file.
    /// </summary>
    public static KalignProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null, string? stdin = null)
    {
        var messages = new List<string>();
        KalignCommand cmd;
        try
        {
            cmd = Parse(arguments);
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, string.Empty, new Dictionary<string, string>(), messages);
        }
        if (cmd.Version)
            return new KalignProgramResult(0, Version + "\n", new Dictionary<string, string>(), messages);
        if (cmd.ShowWarranty)
            return new KalignProgramResult(0, WarrantyText, new Dictionary<string, string>(), messages);
        if (cmd.Help)
            return new KalignProgramResult(0, HelpText, new Dictionary<string, string>(), messages);
        files ??= new Dictionary<string, string>();

        var names = new List<string>();
        var seqs = new List<string>();
        var sources = new List<string>();
        foreach (var f in cmd.InputFiles)
        {
            if (!files.TryGetValue(f, out var text))
            {
                messages.Add($"File: {f} does not exist.");
                return new KalignProgramResult(1, string.Empty, new Dictionary<string, string>(), messages);
            }
            sources.Add(text);
        }
        if (cmd.InputFiles.Count == 0 && stdin != null)
            sources.Add(stdin);
        if (sources.Count == 0)
        {
            messages.Add("No input files");
            return new KalignProgramResult(0, HelpText, new Dictionary<string, string>(), messages);
        }
        foreach (var text in sources)
            ReadFasta(text, names, seqs);

        if (cmd.Format != null && !cmd.Format.Contains("fasta", StringComparison.Ordinal) && !cmd.Format.Contains("fa", StringComparison.Ordinal))
            throw new NotSupportedException($"Kalign output format '{cmd.Format}' is accepted by the program but not written here (B09.md A2-20b); use --format fasta.");

        KalignRunResult result;
        try
        {
            result = SequenceAligner.KalignRun(seqs, cmd.Options);
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, string.Empty, new Dictionary<string, string>(), messages);
        }
        var sb = new StringBuilder();
        for (int i = 0; i < result.KeptIndices.Length; i++)
        {
            sb.Append('>').Append(names[result.KeptIndices[i]]).Append('\n');
            string row = result.Alignment.AlignedSequences[i];
            for (int p = 0; p < row.Length; p += 60)
                sb.Append(row, p, Math.Min(60, row.Length - p)).Append('\n');
        }
        if (cmd.OutFile != null)
            return new KalignProgramResult(0, string.Empty, new Dictionary<string, string> { [cmd.OutFile] = sb.ToString() }, messages);
        return new KalignProgramResult(0, sb.ToString(), new Dictionary<string, string>(), messages);
    }

    /// <summary>read_file_stdin / read_fasta: a '>' line starts a sequence (its whole text is the name), other lines add residues; a line ends at its first control character.</summary>
    private static void ReadFasta(string text, List<string> names, List<string> seqs)
    {
        var current = new List<StringBuilder>();
        foreach (var raw in text.Split('\n'))
        {
            int cut = raw.Length;
            for (int i = 0; i < raw.Length; i++)
            {
                if (raw[i] < 32 || raw[i] == 127)
                {
                    cut = i;
                    break;
                }
            }
            string line = raw[..cut];
            if (line.Length == 0)
                continue;
            if (line[0] == '>')
            {
                names.Add(line[1..]);
                current.Add(new StringBuilder());
            }
            else
            {
                if (current.Count == 0)
                    throw new ArgumentException("Encountered a sequence before encountering it's name");
                current[^1].Append(line);
            }
        }
        seqs.AddRange(current.Select(sb => sb.ToString()));
    }

    private const string WarrantyText =
        "Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n\n15. Disclaimer of Warranty.\n" +
        "THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\nAPPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n" +
        "HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\nOF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n" +
        "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\nPURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n" +
        "IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\nALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n" +
        "A complete copy of the GPL can be found in the \"COPYING\" file.\n";

    private const string HelpText =
        "\nUsage: kalign  -i <seq file> -o <out aln> \n\nOptions:\n\n   --format           : Output format. [Fasta]\n   --type             : Alignment type (rna, dna, internal). [rna]\n" +
        "                        Options: protein, divergent (protein) \n                                 rna, dna, internal (nuc). \n   --gpo              : Gap open penalty. []\n" +
        "   --gpe              : Gap extension penalty. []\n   --tgpe             : Terminal gap extension penalty. []\n   -n/--nthreads      : Number of threads. [4]\n" +
        "   --version (-V/-v)  : Prints version. [NA]\n";
}
