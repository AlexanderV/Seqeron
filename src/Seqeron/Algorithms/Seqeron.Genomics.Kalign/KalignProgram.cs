// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from Kalign 3.3.5 (https://github.com/TimoLassmann/kalign), Copyright (C) 2006, 2019, 2020 Timo Lassmann; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Kalign, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Kalign;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>
/// The meaning of a Kalign 3.3.5 command line (<see cref="KalignProgram.Parse"/>; <c>src/run_kalign.c</c> <c>main</c>): the alignment parameters and the options that only concern input and output.
/// </summary>
/// <param name="Options">Alignment parameters (<c>--type</c>, <c>--gpo</c>, <c>--gpe</c>, <c>--tgpe</c>, <c>-n</c>).</param>
/// <param name="InputFiles">Input file names: <c>-i / --input / --infile / --in</c> first, then the non-option arguments (the program also reads standard input when it is not a terminal).</param>
/// <param name="OutFile"><c>-o / --output / --outfile / --out</c>; null = standard output.</param>
/// <param name="Format"><c>-f / --format</c> (<c>fasta</c> default, <c>msf</c>, <c>clu</c>).</param>
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
/// Accepted and documented but without effect on the alignment: <c>-n / --nthreads</c> beyond the one-thread / many-threads distinction of <see cref="KalignOptions.Threads"/>, <c>--set</c>.
/// Formats: input FASTA, MSF and Clustal (detected, the readers of <c>msa_io.c</c> with their quirks: a name ends at the first white space, a Clustal line without white space is read as residues), output
/// FASTA / MSF / Clustal. Not reproduced: letters outside the alphabet of the detected type (the program leaves the internal code of such a residue uninitialised; this library refuses), the text of error
/// messages with the program's source locations, the exact CPU times (<see cref="KalignRuntime"/>).
/// </remarks>
public static class KalignProgram
{
    /// <summary>Key of <see cref="Exception.Data"/> on the exception of <c>Number of threads has to be &gt;= 1.</c>: the program prints its help and the log line before it fails (the value is the <c>-q</c> flag).</summary>
    private const string UsageKey = "kalign.usage.quiet";


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
        SequenceAligner.RequireNoNullElements(arguments, nameof(arguments), "Argument");
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
                throw new ArgumentException("Number of threads has to be >= 1.") { Data = { [UsageKey] = quiet } };
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
    /// Runs the program: reads the input files (name -> text in FASTA, MSF or Clustal format, detected from the first 100 lines as <c>detect_alignment_format</c> does; a missing name is an error;
    /// <paramref name="stdin"/> is read first when it is given, as the program reads a standard input that is not a terminal), aligns with <see cref="KalignAligner.KalignRun"/> and writes the alignment in the
    /// <c>--format</c> (<c>fasta</c> 60 columns, <c>msf</c> with the GCG checksums and the date line, <c>clu</c>) to standard output or to the <c>-o</c> file. Without <c>-q</c> the banner, the log lines
    /// (<c>Detected DNA sequences.</c>, <c>Read n sequences from f.</c>, the stopwatch lines, the warnings about unaligned / aligned input) and, with <c>-n 1</c>, the <c>Aligned a and b</c> lines precede the alignment
    /// on standard output exactly as the program prints them (the time of the day and the CPU time come from <paramref name="runtime"/>).
    /// </summary>
    public static KalignProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string>? files = null, string? stdin = null, KalignRuntime? runtime = null)
    {
        SequenceAligner.RequireProgramCall(arguments, files);
        runtime ??= new KalignRuntime();
        var messages = new List<string>();
        var stdout = new StringBuilder();
        KalignCommand cmd;
        try
        {
            cmd = Parse(arguments);
        }
        catch (ArgumentException e) when (e.Data.Contains(UsageKey))
        {
            // the program prints its banner (without -q), the help text and the log line, then fails
            bool quiet = (bool)e.Data[UsageKey]!;
            if (!quiet)
                stdout.Append(BannerText);
            stdout.Append(HelpText);
            new KalignLog(runtime, stdout, quiet).Log(e.Message);
            return new KalignProgramResult(1, stdout.ToString(), new Dictionary<string, string>(), messages);
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, string.Empty, new Dictionary<string, string>(), messages);
        }
        if (cmd.Version)
            return new KalignProgramResult(0, Version + "\n", new Dictionary<string, string>(), messages);
        if (!cmd.Quiet)
            stdout.Append(BannerText);
        if (cmd.ShowWarranty)
            return new KalignProgramResult(0, stdout + WarrantyText, new Dictionary<string, string>(), messages);
        if (cmd.Help)
            return new KalignProgramResult(0, stdout + HelpText, new Dictionary<string, string>(), messages);
        files ??= new Dictionary<string, string>();
        var log = new KalignLog(runtime, stdout, cmd.Quiet);

        // the inputs in the program's order: standard input (when there is one), -i, the other arguments
        var sources = new List<(string? Name, string? Text)>();
        if (stdin != null)
            sources.Add((null, stdin));
        foreach (var f in cmd.InputFiles)
            sources.Add((f, files.TryGetValue(f, out var text) ? text : null));
        if (sources.Count == 0)
        {
            stdout.Append(HelpText);
            log.Log("No input files");
            return new KalignProgramResult(0, stdout.ToString(), new Dictionary<string, string>(), messages);
        }

        var all = new List<KalignInputSequence>();
        bool haveMsa = false;
        bool? firstIsDna = null;
        try
        {
            foreach (var (name, text) in sources)
            {
                if (text == null)
                    throw new ArgumentException($"File: {name} does not exist.");
                log.Start();
                var lines = KalignMsaIo.ReadLines(text);
                // an empty input and an input of unknown format both set the alignment read so far to NULL (the program loses it)
                if (KalignMsaIo.IsEmptyInput(lines))
                {
                    all.Clear();
                    haveMsa = false;
                    continue;
                }
                var format = KalignMsaIo.Detect(lines);
                if (format == KalignMsaIo.Format.Fail)
                {
                    log.Warning(name != null ? $"Could not detect input in file: {name}" : "Could not detect input in standard input", "msa_io.c", name != null ? 129 : 131);
                    all.Clear();
                    haveMsa = false;
                    continue;
                }
                var read = format switch
                {
                    KalignMsaIo.Format.Fasta => KalignMsaIo.ReadFasta(lines),
                    KalignMsaIo.Format.Msf => KalignMsaIo.ReadMsf(lines),
                    _ => KalignMsaIo.ReadClustal(lines),
                };
                if (read.Count == 0)
                    throw new ArgumentException("The input holds no sequence (the program crashes with a segmentation fault).");
                // detect_alphabet / detect_aligned of the sequences of this file, then (for the second and later files) of the merged set
                bool? dna = DetectAlphabet(LetterFrequencies(read), log);
                DetectAligned(read, log);
                log.Stop();
                if (!cmd.Quiet)
                {
                    log.Log(name != null ? $"Read {read.Count} sequences from {name}." : $"Read {read.Count} sequences from standard input.");
                    log.Timing();
                }
                if (haveMsa)
                {
                    if (firstIsDna != null && firstIsDna != dna)
                        throw new ArgumentException("Input alignments have different alphabets");
                    all.AddRange(read);
                    firstIsDna = DetectAlphabet(LetterFrequencies(all), log);
                    DetectAligned(all, log);
                }
                else
                {
                    all.AddRange(read);
                    firstIsDna = dna;
                    haveMsa = true;
                }
                if (all.Count < 2)
                    throw new ArgumentException(all.Count == 0 ? "No sequences were found in the input files or standard input." : "Only 1 sequence was found in the input files or standard input");
            }
            if (!haveMsa)
                throw new ArgumentException("No alignment");
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, stdout.ToString(), new Dictionary<string, string>(), messages);
        }

        KalignRunResult result;
        int[] residueCounts;
        try
        {
            result = KalignAligner.KalignRunCore(all.Select(a => a.Raw), cmd.Options, log, all.Select(a => a.Name).ToList(), out residueCounts);
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, stdout.ToString(), new Dictionary<string, string>(), messages);
        }
        var names = result.KeptIndices.Select(i => all[i].Name).ToList();
        var rows = result.Alignment.AlignedSequences;
        string outFormat = cmd.Format ?? "fasta";
        string written;
        try
        {
            if (outFormat.Contains("msf", StringComparison.Ordinal))
                written = KalignMsaIo.WriteMsf(names, rows, residueCounts, cmd.OutFile, runtime.Now());
            else if (outFormat.Contains("clu", StringComparison.Ordinal))
                written = KalignMsaIo.WriteClustal(names, rows);
            else
                written = KalignMsaIo.WriteFasta(names, rows);
        }
        catch (ArgumentException e)
        {
            messages.Add(e.Message);
            return new KalignProgramResult(1, stdout.ToString(), new Dictionary<string, string>(), messages);
        }
        if (cmd.OutFile != null)
            return new KalignProgramResult(0, stdout.ToString(), new Dictionary<string, string> { [cmd.OutFile] = written }, messages);
        stdout.Append(written);
        return new KalignProgramResult(0, stdout.ToString(), new Dictionary<string, string>(), messages);
    }

    /// <summary>The <c>letter_freq</c> of the sequences (every character after the name enters, letters and punctuation alike).</summary>
    private static int[] LetterFrequencies(IReadOnlyList<KalignInputSequence> seqs)
    {
        var freq = new int[128];
        foreach (var s in seqs)
            foreach (char c in s.Raw)
                if (c < 128)
                    freq[c]++;
        return freq;
    }

    /// <summary><c>detect_alphabet</c> (<c>msa_op.c</c>): true = DNA, false = protein, null = the letter counts fit both equally (the program warns and leaves the alphabet undefined).</summary>
    private static bool? DetectAlphabet(int[] freq, KalignLog log)
    {
        var dna = new double[128];
        var protein = new double[128];
        for (int i = 0; i < 128; i++)
        {
            dna[i] = Math.Log(0.0001 * 1.0 / 116.0);
            protein[i] = Math.Log(0.0001 * 1.0 / 88.0);
        }
        foreach (char c in "acgtunACGTUN")
            dna[c] = Math.Log(0.9999 * 1.0 / 12.0);
        foreach (char c in "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY")
            protein[c] = Math.Log(0.9999 * 1.0 / 40.0);
        double d = 0.0, p = 0.0;
        for (int i = 0; i < 128; i++)
        {
            if (freq[i] != 0)
            {
                d += dna[i] * freq[i];
                p += protein[i] * freq[i];
            }
        }
        if (d == p)
        {
            log.Warning("Could not determine whether we have a DNA or Protein alignment", "msa_op.c", 124);
            return null;
        }
        if (!log.Quiet)
            log.Log(d > p ? "Detected DNA sequences." : "Detected protein sequences.");
        return d > p;
    }

    /// <summary><c>detect_aligned</c>: the warnings about input that has gaps but rows of different total length, or equal lengths without gaps.</summary>
    private static void DetectAligned(IReadOnlyList<KalignInputSequence> seqs, KalignLog log)
    {
        int min = int.MaxValue, max = 0;
        long gaps = 0;
        var punctuation = new int[128];
        foreach (var s in seqs)
        {
            int letters = 0, punct = 0;
            foreach (char c in s.Raw)
            {
                if (c is >= 'a' and <= 'z' or >= 'A' and <= 'Z')
                    letters++;
                else if (c > ' ' && c < 127 && !(c >= '0' && c <= '9'))
                {
                    punct++;
                    punctuation[c]++;
                }
            }
            gaps += punct;
            min = Math.Min(min, letters + punct);
            max = Math.Max(max, letters + punct);
        }
        if (log.Quiet)
            return;
        if (gaps > 0)
        {
            if (min != max)
            {
                log.Warning("--------------------------------------------", "msa_op.c", 413);
                log.Warning("The input sequences contain gap characters: ", "msa_op.c", 414);
                for (int i = 0; i < 128; i++)
                    if (punctuation[i] > 0)
                        log.Warning($"\"{(char)i}\" : {punctuation[i],4} found                            ", "msa_op.c", 418);
                log.Warning("BUT the presumably aligned sequences do not ", "msa_op.c", 422);
                log.Warning("have the same length.                       ", "msa_op.c", 423);
                log.Warning("                                            ", "msa_op.c", 424);
                log.Warning("Kalign will remove the gap characters and   ", "msa_op.c", 425);
                log.Warning("align the sequences.                        ", "msa_op.c", 426);
                log.Warning("--------------------------------------------", "msa_op.c", 427);
            }
        }
        else if (min == max)
        {
            log.Warning("--------------------------------------------", "msa_op.c", 434);
            log.Warning("All input sequences have the same length.   ", "msa_op.c", 435);
            log.Warning("BUT there are no gap characters.            ", "msa_op.c", 436);
            log.Warning("                                            ", "msa_op.c", 437);
            log.Warning("Unable to determine whether the sequences   ", "msa_op.c", 438);
            log.Warning("are already aligned.                        ", "msa_op.c", 439);
            log.Warning("Kalign will align the sequences.            ", "msa_op.c", 440);
            log.Warning("--------------------------------------------", "msa_op.c", 441);
        }
    }

    private const string BannerText =
        "\nKalign (3.3.5)\n\nCopyright (C) 2006,2019,2020,2021 Timo Lassmann\n\nThis program comes with ABSOLUTELY NO WARRANTY; for details type:\n`kalign -showw'.\n" +
        "This is free software, and you are welcome to redistribute it\nunder certain conditions; consult the COPYING file for details.\n\nPlease cite:\n  Lassmann, Timo.\n" +
        "  \"Kalign 3: multiple sequence alignment of large data sets.\"\n  Bioinformatics (2019) \n  https://doi.org/10.1093/bioinformatics/btz795\n\n";

    private const string WarrantyText =
        "Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n\n15. Disclaimer of Warranty.\n" +
        "THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\nAPPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n" +
        "HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\nOF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n" +
        "THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\nPURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n" +
        "IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\nALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\n" +
        "A complete copy of the GPL can be found in the \"COPYING\" file.\n";

    private const string HelpText =
        "\nUsage: kalign-avx2  -i <seq file> -o <out aln> \n\nOptions:\n\n   --format           : Output format. [Fasta]\n   --type             : Alignment type (rna, dna, internal). [rna]\n" +
        "                        Options: protein, divergent (protein) \n                                 rna, dna, internal (nuc). \n   --gpo              : Gap open penalty. []\n" +
        "   --gpe              : Gap extension penalty. []\n   --tgpe             : Terminal gap extension penalty. []\n   -n/--nthreads      : Number of threads. [4]\n" +
        "   --version (-V/-v)  : Prints version. [NA]\n\nExamples:\n\nPassing sequences via stdin:\n\n   cat input.fa | kalign -f fasta > out.afa\n\n" +
        "Combining multiple input files:\n\n   kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa\n\n";
}
