using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-26c: SPOA 4.0.8 (port of the program, verified against the installed binary); delegates to SpoaProgram.Run.
public partial class AlignmentTools
{
    [McpServerTool(Name = "spoa_command", Title = "SPOA 4.0.8 - command line (FASTA / FASTQ files)", ReadOnly = true)]
    [Description("SPOA 4.0.8 (Vaser et al. 2017, Genome Research 27:737; partial order alignment with vectorized dynamic programming) as a function of its command line and the files it names (SpoaProgram.Run): spoa [-m MATCH] [-n MISMATCH] [-g GAP_OPEN] [-e GAP_EXTEND] [-q GAP_OPEN2] [-c GAP_EXTEND2] [-l|--algorithm 0 local / 1 global / 2 semi-global] [-r|--result 0 consensus / 1 MSA / 2 both / 3 GFA / 4 consensus + GFA (repeatable)] [-d|--dot FILE] [-s|--strand-ambiguous] [--version] [-h|--help] SEQUENCES.fa|.fasta|.fna|.faa|.fq|.fastq (optionally .gz). The argument loop is glibc's getopt_long (permutation, unique abbreviations of the long options, the diagnostics), the numeric options are C `atoi` values stored in 8-bit integers, the gap model is linear (g >= e), affine (g <= q or e >= c) or convex; the file parsers are bioparser's (names end at the first white space, a FASTQ file gives Phred weights, an invalid file or an empty input aborts with status 134 and the `terminate called` text of the C++ runtime). The alignment engine is the SIMD engine of the program emulated lane by lane (int16 / int32 vectors chosen by the worst-case score, padding lanes, shifted prefix maxima, the traceback of the program, lanes = 16 / 8 as on AVX2; the result does not depend on the vector width), so the ties are the program's. Returns the exit status, standard output (consensus / MSA / GFA text, --help, --version), standard error and the written DOT file. VERIFICATION: exit status, standard output, standard error and the DOT file equal the installed Debian spoa 4.0.8-1 on 441 / 441 committed command lines (scripts/alignment/generate_spoa_oracle.py).")]
    public static PoaCommandResultDto SpoaCommand(
        [Description("The command-line arguments of spoa, one per element (for example [\"-l\", \"1\", \"-r\", \"2\", \"in.fa\"]).")] string[] arguments,
        [Description("The text of the input files the arguments name: file name -> text (FASTA / FASTQ; bytes are read as Latin-1).")] Dictionary<string, string>? files = null,
        [Description("Binary input files (for example gzip-compressed FASTA): file name -> base64.")] Dictionary<string, string>? binaryFiles = null)
    {
        if (arguments == null)
            throw new ArgumentException("arguments must not be null.", nameof(arguments));
        var bytes = new Dictionary<string, byte[]>();
        foreach (var kv in files ?? [])
            bytes[kv.Key] = Encoding.Latin1.GetBytes(kv.Value);
        foreach (var kv in binaryFiles ?? [])
            bytes[kv.Key] = Convert.FromBase64String(kv.Value);
        var result = global::Seqeron.Genomics.Alignment.SpoaProgram.Run(arguments, p => bytes.TryGetValue(p, out var b) ? b : null);
        return new PoaCommandResultDto(result.ExitCode, result.StdOut, result.StdErr, new Dictionary<string, string>(result.Files));
    }
}
