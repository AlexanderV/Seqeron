using System.IO.Compression;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The outcome of one run of the <c>spoa</c> command line (exit status, standard output / error, the files the program wrote).</summary>
public sealed record SpoaRunResult(int ExitCode, string StdOut, string StdErr, IReadOnlyDictionary<string, string> Files);

/// <summary>
/// SPOA 4.0.8 as a command line (<c>spoa [options ...] &lt;sequences&gt;</c>, <c>src/main.cpp</c>): the <c>getopt_long</c> loop of glibc (permutation, abbreviations, the diagnostics), the FASTA / FASTQ parsers of bioparser
/// (64 kB chunks, gzip input, names shortened at the first white space, the invalid-format aborts), the alignment engine (<see cref="SpoaEngine"/>) and the output modes consensus / MSA / GFA / DOT of <see cref="SpoaGraph"/>.
/// An uncaught C++ exception ends the program with the status 134 (SIGABRT) and the <c>terminate called</c> text on standard error.
/// </summary>
public static class SpoaProgram
{
    private const string Version = "4.0.8";

    private const string HelpText =
        "usage: spoa [options ...] <sequences>\n"
        + "\n"
        + "  # default output is stdout\n"
        + "  <sequences>\n"
        + "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        + "\n"
        + "  options:\n"
        + "    -m <int>\n"
        + "      default: 5\n"
        + "      score for matching bases\n"
        + "    -n <int>\n"
        + "      default: -4\n"
        + "      score for mismatching bases\n"
        + "    -g <int>\n"
        + "      default: -8\n"
        + "      gap opening penalty (must be non-positive)\n"
        + "    -e <int>\n"
        + "      default: -6\n"
        + "      gap extension penalty (must be non-positive)\n"
        + "    -q <int>\n"
        + "      default: -10\n"
        + "      gap opening penalty of the second affine function\n"
        + "      (must be non-positive)\n"
        + "    -c <int>\n"
        + "      default: -4\n"
        + "      gap extension penalty of the second affine function\n"
        + "      (must be non-positive)\n"
        + "    -l, --algorithm <int>\n"
        + "      default: 0\n"
        + "      alignment mode:\n"
        + "        0 - local (Smith-Waterman)\n"
        + "        1 - global (Needleman-Wunsch)\n"
        + "        2 - semi-global\n"
        + "    -r, --result <int> (option can be used multiple times)\n"
        + "      default: 0\n"
        + "      result mode:\n"
        + "        0 - consensus (FASTA)\n"
        + "        1 - multiple sequence alignment (FASTA)\n"
        + "        2 - 0 & 1 (FASTA)\n"
        + "        3 - partial order graph (GFA)\n"
        + "        4 - 0 & 3 (GFA)\n"
        + "    -d, --dot <file>\n"
        + "      output file for the partial order graph in DOT format\n"
        + "    -s, --strand-ambiguous\n"
        + "      for each sequence pick the strand with the better alignment\n"
        + "    --version\n"
        + "      prints the version number\n"
        + "    -h, --help\n"
        + "      prints the usage\n"
        + "\n"
        + "  gap mode:\n"
        + "    linear if g >= e\n"
        + "    affine if g <= q or e >= c\n"
        + "    convex otherwise (default)\n";

    /// <summary>An uncaught C++ exception of the program (<c>std::terminate</c>).</summary>
    public sealed class SpoaAbortException(string type, string what) : Exception(what)
    {
        /// <summary>The C++ exception type shown by the terminate handler.</summary>
        public string Type { get; } = type;
    }

    private sealed record Seq(string Name, string Data, string Quality);

    private static readonly (string Name, bool HasArg, char Code)[] LongOptions =
    [
        ("algorithm", true, 'l'),
        ("result", true, 'r'),
        ("dot", true, 'd'),
        ("strand-ambiguous", false, 's'),
        ("version", false, 'v'),
        ("help", false, 'h'),
    ];

    /// <summary>C's <c>atoi</c>: leading white space, a sign and digits (saturating at the <c>long</c> range, then truncated to <c>int</c> as glibc does).</summary>
    private static int Atoi(string s)
    {
        int i = 0;
        while (i < s.Length && (s[i] == ' ' || (s[i] >= '\t' && s[i] <= '\r')))
            i++;
        bool neg = false;
        if (i < s.Length && (s[i] == '+' || s[i] == '-'))
        {
            neg = s[i] == '-';
            i++;
        }

        decimal v = 0;
        while (i < s.Length && s[i] >= '0' && s[i] <= '9')
        {
            v = (v * 10) + (s[i] - '0');
            if (v > long.MaxValue)
                v = (decimal)long.MaxValue + 1;
            i++;
        }

        bool saturated = v >= (decimal)long.MaxValue + 1;
        long l;
        if (saturated)
            l = neg ? long.MinValue : long.MaxValue;
        else
            l = neg ? -(long)v : (long)v;
        return unchecked((int)l);
    }

    /// <summary>Runs the program. <paramref name="args"/> are the arguments after the program name; <paramref name="readFile"/> returns the bytes of a path or null; <paramref name="registerBits"/> = 256 (AVX2, default) or 128 (SSE4.1): the vector width of the emulated SIMD engine.</summary>
    public static SpoaRunResult Run(IReadOnlyList<string> args, Func<string, byte[]?> readFile, int registerBits = 256)
    {
        ArgumentNullException.ThrowIfNull(args);
        ArgumentNullException.ThrowIfNull(readFile);
        var stdout = new StringBuilder();
        var stderr = new StringBuilder();
        var files = new Dictionary<string, string>();
        int code;
        try
        {
            code = Main(args, readFile, stdout, stderr, files, registerBits);
        }
        catch (SpoaAbortException a)
        {
            stderr.Append("terminate called after throwing an instance of '").Append(a.Type).Append("'\n  what():  ").Append(a.Message).Append('\n');
            code = 134;
        }

        return new SpoaRunResult(code, stdout.ToString(), stderr.ToString(), files);
    }

    private static int Main(IReadOnlyList<string> args, Func<string, byte[]?> readFile, StringBuilder stdout, StringBuilder stderr, Dictionary<string, string> files, int registerBits)
    {
        int m = 5, n = -4, g = -8, e = -6, q = -10, c = -4;
        int algorithm = 0;
        var results = new List<int> { 0 };
        string dotPath = "";
        bool isStrandAmbiguous = false;

        // getopt_long: permute the arguments (options first, then the operands in order)
        var operands = new List<string>();
        int pos = 0;
        var argv = args.ToList();
        while (pos < argv.Count)
        {
            string arg = argv[pos];
            if (arg == "--")
            {
                operands.AddRange(argv.Skip(pos + 1));
                break;
            }

            if (arg.Length < 2 || arg[0] != '-')
            {
                operands.Add(arg);
                pos++;
                continue;
            }

            if (arg.StartsWith("--", StringComparison.Ordinal))
            {
                string body = arg[2..];
                string? value = null;
                int eq = body.IndexOf('=');
                string name = body;
                if (eq >= 0)
                {
                    name = body[..eq];
                    value = body[(eq + 1)..];
                }

                var matches = LongOptions.Where(o => o.Name == name).ToList();
                if (matches.Count == 0)
                    matches = [.. LongOptions.Where(o => o.Name.StartsWith(name, StringComparison.Ordinal))];
                if (matches.Count == 0)
                {
                    stderr.Append("spoa: unrecognized option '--").Append(name).Append("'\n");
                    return 1;
                }

                if (matches.Count > 1 && matches.Select(o => o.Code).Distinct().Count() > 1)
                {
                    stderr.Append("spoa: option '--").Append(name).Append("' is ambiguous; possibilities:").Append(string.Concat(matches.Select(o => " '--" + o.Name + "'"))).Append('\n');
                    return 1;
                }

                var opt = matches[0];
                pos++;
                if (opt.HasArg)
                {
                    if (value is null)
                    {
                        if (pos >= argv.Count)
                        {
                            stderr.Append("spoa: option '--").Append(opt.Name).Append("' requires an argument\n");
                            return 1;
                        }

                        value = argv[pos++];
                    }
                }
                else if (value is not null)
                {
                    stderr.Append("spoa: option '--").Append(opt.Name).Append("' doesn't allow an argument\n");
                    return 1;
                }

                int rc = Apply(opt.Code, value, stdout, ref m, ref n, ref g, ref e, ref q, ref c, ref algorithm, results, ref dotPath, ref isStrandAmbiguous);
                if (rc >= 0)
                    return rc;
                continue;
            }

            // a cluster of short options
            pos++;
            int k = 1;
            while (k < arg.Length)
            {
                char ch = arg[k++];
                bool takesArg = "mngeqclrd".Contains(ch);
                bool known = "mngeqclrdsh".Contains(ch) && ch != ':';
                if (!known)
                {
                    stderr.Append("spoa: invalid option -- '").Append(ch).Append("'\n");
                    return 1;
                }

                string? value = null;
                if (takesArg)
                {
                    if (k < arg.Length)
                    {
                        value = arg[k..];
                        k = arg.Length;
                    }
                    else if (pos < argv.Count)
                    {
                        value = argv[pos++];
                    }
                    else
                    {
                        stderr.Append("spoa: option requires an argument -- '").Append(ch).Append("'\n");
                        return 1;
                    }
                }

                int rc = Apply(ch, value, stdout, ref m, ref n, ref g, ref e, ref q, ref c, ref algorithm, results, ref dotPath, ref isStrandAmbiguous);
                if (rc >= 0)
                    return rc;
            }
        }

        if (results.Count > 1)
            results.RemoveAt(0);
        if (operands.Count == 0)
        {
            stderr.Append("[spoa::] error: missing input file!\n");
            stdout.Append(HelpText);
            return 1;
        }

        string path = operands[0];
        bool fasta = new[] { ".fasta", ".fasta.gz", ".fna", ".fna.gz", ".faa", ".faa.gz", ".fa", ".fa.gz" }.Any(path.EndsWith);
        bool fastq = new[] { ".fastq", ".fastq.gz", ".fq", ".fq.gz" }.Any(path.EndsWith);
        if (!fasta && !fastq)
        {
            stderr.Append("[spoa::CreateParser] error: file ").Append(path).Append(" has unsupported format extension (valid extensions: .fasta, .fasta.gz, .fna, .fna.gz, .faa, .faa.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)\n");
            return 1;
        }

        var bytes = readFile(path);
        if (bytes is null)
        {
            stderr.Append("[bioparser::Parser::Create] error: unable to open file ").Append(path).Append('\n');
            return 1;
        }

        SpoaEngine engine;
        try
        {
            engine = SpoaEngine.Create((byte)algorithm, (sbyte)m, (sbyte)n, (sbyte)g, (sbyte)e, (sbyte)q, (sbyte)c, registerBits);
        }
        catch (ArgumentException ex)
        {
            stderr.Append(ex.Message).Append('\n');
            return 1;
        }

        var sequences = Parse(Decompress(bytes), fasta);
        uint maxLen = 0;
        foreach (var s in sequences)
            maxLen = Math.Max(maxLen, (uint)s.Data.Length);
        if (SpoaEngine.PreallocFails(maxLen))
            throw new SpoaAbortException("std::length_error", "vector::_M_fill_insert");
        var graph = new SpoaGraph();
        var isReversed = new List<bool>();
        for (int idx = 0; idx < sequences.Count; idx++)
        {
            var it = sequences[idx];
            string data = it.Data, quality = it.Quality;
            try
            {
                var alignment = engine.Align(data, graph, out int score);
                if (isStrandAmbiguous)
                {
                    string rc = ReverseAndComplement(data);
                    string rq = new(quality.Reverse().ToArray());
                    var alignmentRev = engine.Align(rc, graph, out int scoreRev);
                    if (score >= scoreRev)
                    {
                        data = ReverseAndComplement(rc);
                        isReversed.Add(false);
                    }
                    else
                    {
                        data = rc;
                        quality = rq;
                        alignment = alignmentRev;
                        isReversed.Add(true);
                    }
                }

                if (quality.Length == 0)
                    graph.AddAlignment(alignment, data);
                else
                    graph.AddAlignment(alignment, data, quality);
            }
            catch (ArgumentException ex)
            {
                stderr.Append(ex.Message).Append('\n');
                return 1;
            }
        }

        foreach (int it in results)
        {
            switch (it)
            {
                case 0:
                    var consensus = graph.GenerateConsensus();
                    stdout.Append(">Consensus LN:i:").Append(consensus.Length).Append('\n').Append(consensus).Append('\n');
                    break;
                case 1:
                case 2:
                    var msa = graph.GenerateMultipleSequenceAlignment(it == 2);
                    for (int i = 0; i < msa.Count; i++)
                    {
                        string name = i < sequences.Count ? sequences[i].Name : "Consensus";
                        stdout.Append('>').Append(name).Append('\n').Append(msa[i]).Append('\n');
                    }

                    break;
                case 3:
                case 4:
                    graph.GenerateConsensus();
                    PrintGfa(graph, sequences.Select(s => s.Name).ToList(), isReversed, it == 4, stdout, stderr);
                    break;
            }
        }

        if (dotPath.Length > 0)
            files[dotPath] = graph.PrintDot();
        return 0;
    }

    /// <summary>One option of the loop; returns -1 to continue or the exit status.</summary>
    private static int Apply(char opt, string? optarg, StringBuilder stdout, ref int m, ref int n, ref int g, ref int e, ref int q, ref int c, ref int algorithm, List<int> results, ref string dotPath, ref bool strand)
    {
        switch (opt)
        {
            case 'm': m = (sbyte)Atoi(optarg!); break;
            case 'n': n = (sbyte)Atoi(optarg!); break;
            case 'g': g = (sbyte)Atoi(optarg!); break;
            case 'e': e = (sbyte)Atoi(optarg!); break;
            case 'q': q = (sbyte)Atoi(optarg!); break;
            case 'c': c = (sbyte)Atoi(optarg!); break;
            case 'l': algorithm = (byte)Atoi(optarg!); break;
            case 'r': results.Add((byte)Atoi(optarg!)); break;
            case 'd': dotPath = optarg!; break;
            case 's': strand = true; break;
            case 'v':
                stdout.Append(Version).Append('\n');
                return 0;
            case 'h':
                stdout.Append(HelpText);
                return 0;
        }

        return -1;
    }

    private static byte[] Decompress(byte[] bytes)
    {
        if (bytes.Length >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b)
        {
            using var input = new MemoryStream(bytes);
            using var gz = new GZipStream(input, CompressionMode.Decompress);
            using var output = new MemoryStream();
            gz.CopyTo(output);
            return output.ToArray();
        }

        return bytes;
    }

    private static string ReverseAndComplement(string data)
    {
        var chars = data.ToCharArray();
        for (int i = 0; i < chars.Length; i++)
        {
            switch (char.ToUpperInvariant(chars[i]))
            {
                case 'A': chars[i] = 'T'; break;
                case 'C': chars[i] = 'G'; break;
                case 'G': chars[i] = 'C'; break;
                case 'T':
                case 'U': chars[i] = 'A'; break;
                case 'R': chars[i] = 'Y'; break;
                case 'Y': chars[i] = 'R'; break;
                case 'K': chars[i] = 'M'; break;
                case 'M': chars[i] = 'K'; break;
                case 'B': chars[i] = 'V'; break;
                case 'D': chars[i] = 'H'; break;
                case 'H': chars[i] = 'D'; break;
                case 'V': chars[i] = 'B'; break;
            }
        }

        Array.Reverse(chars);
        return new string(chars);
    }

    private static void PrintGfa(SpoaGraph graph, List<string> headers, List<bool> isReversed, bool includeConsensus, StringBuilder stdout, StringBuilder stderr)
    {
        if (headers.Count < graph.Sequences.Count)
        {
            stderr.Append("[spoa::PrintGfa] error: missing header(s)\n");
            return;
        }

        if (isReversed.Count > 0 && isReversed.Count < graph.Sequences.Count)
        {
            stderr.Append("[spoa::PringGfa] error: missing reversion flag(s)\n");
            return;
        }

        var isConsensusNode = new bool[graph.Nodes.Count];
        foreach (var it in graph.Consensus)
            isConsensusNode[it.Id] = true;
        stdout.Append("H\tVN:Z:1.0\n");
        foreach (var it in graph.Nodes)
        {
            stdout.Append("S\t").Append(it.Id + 1).Append('\t').Append(graph.Decoder(it.Code));
            if (isConsensusNode[it.Id])
                stdout.Append("\tic:Z:true");
            stdout.Append('\n');
            foreach (var jt in it.OutEdges)
            {
                stdout.Append("L\t").Append(it.Id + 1).Append("\t+\t").Append(jt.Head.Id + 1).Append("\t+\tOM\tew:f:").Append(jt.Weight);
                if (isConsensusNode[it.Id] && isConsensusNode[jt.Head.Id])
                    stdout.Append("\tic:Z:true");
                stdout.Append('\n');
            }
        }

        for (int i = 0; i < graph.Sequences.Count; i++)
        {
            stdout.Append("P\t").Append(headers[i]).Append('\t');
            var path = new List<uint>();
            var curr = graph.Sequences[i];
            while (curr is not null)
            {
                path.Add(curr.Id + 1);
                curr = curr.Successor((uint)i);
            }

            bool ir = isReversed.Count > 0 && isReversed[i];
            if (ir)
                path.Reverse();
            for (int j = 0; j < path.Count; j++)
            {
                if (j != 0)
                    stdout.Append(',');
                stdout.Append(path[j]).Append(ir ? '-' : '+');
            }

            stdout.Append("\t*\n");
        }

        if (includeConsensus)
        {
            stdout.Append("P\tConsensus\t");
            for (int i = 0; i < graph.Consensus.Count; i++)
            {
                if (i != 0)
                    stdout.Append(',');
                stdout.Append(graph.Consensus[i].Id + 1).Append('+');
            }

            stdout.Append("\t*\n");
        }
    }

    // ------------------------------------------------------------------ bioparser

    private sealed class Chunker
    {
        private readonly byte[] _input;
        private int _inPos;
        public readonly byte[] Buffer = new byte[65536];
        public uint BufferPtr;
        public uint BufferBytes;
        public byte[] Storage = new byte[4194304];
        public uint StoragePtr;

        public Chunker(byte[] input) => _input = input;

        /// <summary><c>Parser::Read</c>: true when the buffer was not filled (end of input).</summary>
        public bool Read()
        {
            BufferPtr = 0;
            int count = Math.Min(Buffer.Length, _input.Length - _inPos);
            Array.Copy(_input, _inPos, Buffer, 0, count);
            _inPos += count;
            BufferBytes = (uint)count;
            return count < Buffer.Length;
        }

        public void Store(uint count, bool strip)
        {
            if (BufferPtr + count > Buffer.Length)
                throw new SpoaAbortException("std::invalid_argument", "[bioparser::Parser::Store] error: buffer overflow");
            if (StoragePtr + count > Storage.Length)
                Array.Resize(ref Storage, 2 * Storage.Length);
            Array.Copy(Buffer, BufferPtr, Storage, StoragePtr, count);
            StoragePtr += strip ? RightStrip(Storage, StoragePtr, count) : count;
            BufferPtr += count + 1;
        }

        public void Clear() => StoragePtr = 0;

        private static bool IsSpace(byte b) => b is (byte)' ' or (byte)'\t' or (byte)'\n' or 0x0B or 0x0C or (byte)'\r';

        public static uint RightStrip(byte[] data, uint start, uint len)
        {
            while (len > 0 && IsSpace(data[start + len - 1]))
                len--;
            return len;
        }

        public uint Shorten(uint len)
        {
            for (uint i = 0; i < len; i++)
            {
                if (IsSpace(Storage[i]))
                    return i;
            }

            return len;
        }

        public string Text(uint start, uint len) => Encoding.Latin1.GetString(Storage, (int)start, (int)len);
    }

    private static List<Seq> Parse(byte[] input, bool fasta) => fasta ? ParseFasta(input) : ParseFastq(input);

    private static List<Seq> ParseFasta(byte[] input)
    {
        var p = new Chunker(input);
        var dst = new List<Seq>();
        uint dataPtr = 0;

        void Create()
        {
            if (dataPtr == 0)
                throw new SpoaAbortException("std::invalid_argument", "[bioparser::FastaParser] error: invalid file format");
            uint nameLen = p.Shorten(dataPtr);
            uint dataLen = p.StoragePtr - dataPtr;
            if (nameLen == 0 || p.Storage[0] != (byte)'>' || dataLen == 0)
                throw new SpoaAbortException("std::invalid_argument", "[bioparser::FastaParser] error: invalid file format");
            dst.Add(new Seq(p.Text(1, nameLen - 1), p.Text(dataPtr, dataLen), ""));
            dataPtr = 0;
            p.Clear();
        }

        bool isEof = false, isName = true;
        while (true)
        {
            uint bufferPtr = p.BufferPtr;
            for (; bufferPtr < p.BufferBytes; bufferPtr++)
            {
                byte c = p.Buffer[bufferPtr];
                if (c == (byte)'\n')
                {
                    p.Store(bufferPtr - p.BufferPtr, !isName);
                    if (isName)
                    {
                        dataPtr = p.StoragePtr;
                        isName = false;
                    }
                }
                else if (!isName && c == (byte)'>')
                {
                    isName = true;
                    Create();
                }
            }

            if (p.BufferPtr < bufferPtr)
                p.Store(bufferPtr - p.BufferPtr, !isName);
            if (isEof)
                break;
            isEof = p.Read();
        }

        if (p.StoragePtr != 0)
            Create();
        return dst;
    }

    private static List<Seq> ParseFastq(byte[] input)
    {
        var p = new Chunker(input);
        var dst = new List<Seq>();
        uint dataPtr = 0, commentPtr = 0, qualityPtr = 0;

        void Create()
        {
            if (dataPtr == 0 || commentPtr == 0 || qualityPtr == 0)
                throw new SpoaAbortException("std::invalid_argument", "[bioparser::FastqParser] error: invalid file format");
            uint nameLen = p.Shorten(dataPtr);
            uint dataLen = commentPtr - dataPtr;
            uint qualityLen = p.StoragePtr - qualityPtr;
            if (nameLen == 0 || p.Storage[0] != (byte)'@' || dataLen == 0 || qualityLen == 0 || dataLen != qualityLen)
                throw new SpoaAbortException("std::invalid_argument", "[bioparser::FastqParser] error: invalid file format");
            dst.Add(new Seq(p.Text(1, nameLen - 1), p.Text(dataPtr, dataLen), p.Text(qualityPtr, qualityLen)));
            dataPtr = 0;
            commentPtr = 0;
            qualityPtr = 0;
            p.Clear();
        }

        bool isEof = false, isName = true, isData = false, isComment = false, isQuality = false;
        while (true)
        {
            uint bufferPtr = p.BufferPtr;
            for (; bufferPtr < p.BufferBytes; bufferPtr++)
            {
                byte c = p.Buffer[bufferPtr];
                if (c == (byte)'\n')
                {
                    p.Store(bufferPtr - p.BufferPtr, !isName);
                    if (isName)
                    {
                        isName = false;
                        isData = true;
                        dataPtr = p.StoragePtr;
                    }
                    else if (isComment)
                    {
                        isComment = false;
                        isQuality = true;
                        qualityPtr = p.StoragePtr;
                    }
                    else if (isQuality && p.StoragePtr - qualityPtr == commentPtr - dataPtr)
                    {
                        isQuality = false;
                        isName = true;
                        Create();
                    }
                }
                else if (isData && c == (byte)'+')
                {
                    isData = false;
                    isComment = true;
                    commentPtr = p.StoragePtr;
                }
            }

            if (p.BufferPtr < bufferPtr)
                p.Store(bufferPtr - p.BufferPtr, !isName);
            if (isEof)
                break;
            isEof = p.Read();
        }

        if (p.StoragePtr != 0)
            Create();
        return dst;
    }
}
