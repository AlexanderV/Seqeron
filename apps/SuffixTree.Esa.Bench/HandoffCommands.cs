using System.Diagnostics;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;
using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.9f, the hand-off package (<c>analysis/s1b-handoff</c>): the three commands that the production chain of the owner's machine needs besides <c>sesa resume build</c>. <c>sesa handoff env</c>: the CRC-32C path
/// that this CPU takes (the SSE 4.2 instruction on x86-64, the CRC extension on arm64) checked against the managed table path, and the vector widths. <c>sesa handoff open</c>: the container opened at the level <c>Header</c>,
/// <c>Sections</c> or <c>Structure</c> with the seconds of each level and the anonymous resident memory before, at the peak and after (a container of 2^31 - 1 symbols or more has no reader: reported as <c>not_applicable</c>,
/// exit code 4). <c>sesa handoff queries</c>: a fixed deterministic query set answered from the container (the engine, and a search over its suffix array) and compared with an oracle that reads only the input text.
/// </summary>
internal static class HandoffCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    internal static void Run(string[] args)
    {
        if (args.Length == 0 || args[0] is not ("env" or "open" or "queries"))
        {
            throw new ArgumentException("usage: sesa handoff env [--label L] | open --container <file> [--level header|sections|structure] [--access checked|trusted] [--sha false] [--label L] | queries --container <file> --input <text file> [--seed 20261122] [--lengths 12,20,32] [--present 300] [--absent 200] [--cap 1000] [--engine auto|true|false] [--label L]");
        }
        var o = new Dictionary<string, string>();
        for (int i = 1; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            if (!args[i].StartsWith("--", StringComparison.Ordinal)) throw new ArgumentException($"unknown option '{args[i]}'");
            o[args[i][2..]] = value;
        }
        string label = o.GetValueOrDefault("label", args[0]);
        switch (args[0])
        {
            case "env": Env(label); break;
            case "open": Environment.ExitCode = Open(o, label); break;
            default: Environment.ExitCode = Queries(o, label); break;
        }
    }

    private static string Need(Dictionary<string, string> o, string key) => o.TryGetValue(key, out string? v) ? v : throw new ArgumentException($"--{key} is required");

    private static void Row(string metric, string value, string unit, string label, string notes = "") => Emit.Row("handoff", metric, value, unit, string.IsNullOrEmpty(notes) ? label : $"{label} {notes}");

    private static void Env(string label)
    {
        Row("env_process_arch", RuntimeInformation.ProcessArchitecture.ToString(), "", label);
        Row("env_os_arch", RuntimeInformation.OSArchitecture.ToString(), "", label);
        Row("env_os", RuntimeInformation.OSDescription.Replace('\t', ' '), "", label);
        Row("env_little_endian", BitConverter.IsLittleEndian ? "yes" : "no", "bool", label);
        Row("env_crc32c_hardware", Crc32C.HardwareAccelerated ? "yes" : "no", "bool", label, "Crc32C.HardwareAccelerated: the section checksums run on an instruction");
        string instruction = "none (managed slice-by-8 table)";
        if (Sse42.X64.IsSupported) instruction = "Sse42.X64.Crc32";
        else if (Crc32.Arm64.IsSupported) instruction = "Crc32.Arm64.ComputeCrc32C";
        Row("env_crc32c_instruction", instruction, "", label);
        var vector = new byte[(1 << 20) + 13];
        new Random(20261122).NextBytes(vector);
        bool known = Crc32C.Compute("123456789"u8) == 0xE3069283u;
        Row("env_crc32c_check_value", known ? "ok" : "WRONG", "", label, "CRC-32C of 123456789 is E3069283");
        bool same = true;
        foreach (int length in new[] { 0, 1, 7, 8, 9, 31, 32, 33, 255, 4096, 65537, vector.Length })
        {
            uint table = Crc32C.AppendTable(0, vector.AsSpan(0, length));
            uint chosen = Crc32C.Append(0, vector.AsSpan(0, length));
            if (table != chosen || (Crc32C.HardwareAccelerated && Crc32C.AppendHardware(0, vector.AsSpan(0, length)) != table)) same = false;
        }
        Row("env_crc32c_table_equals_hardware", same ? "yes" : "NO", "bool", label, "twelve lengths up to 1 MiB: the table path, the chosen path and the hardware path give one value");
        Row("env_vector128_hardware", System.Runtime.Intrinsics.Vector128.IsHardwareAccelerated ? "yes" : "no", "bool", label);
        Row("env_vector256_hardware", System.Runtime.Intrinsics.Vector256.IsHardwareAccelerated ? "yes" : "no", "bool", label);
        Row("env_cores", Environment.ProcessorCount.ToString(Inv), "n", label, "Environment.ProcessorCount");
    }

    private static int Open(Dictionary<string, string> o, string label)
    {
        string path = Need(o, "container");
        EsaValidationLevel level = o.GetValueOrDefault("level", "header").ToLowerInvariant() switch
        {
            "header" => EsaValidationLevel.Header,
            "sections" => EsaValidationLevel.Sections,
            "structure" => EsaValidationLevel.Structure,
            var other => throw new ArgumentException($"--level must be header, sections or structure, not '{other}'"),
        };
        EsaAccessMode access = o.GetValueOrDefault("access", "checked").ToLowerInvariant() == "trusted" ? EsaAccessMode.Trusted : EsaAccessMode.Checked;
        bool sha = bool.Parse(o.GetValueOrDefault("sha", "false"));
        long n = -1;
        try
        {
            using SesaFile header = SesaFile.Open(path, new EsaOpenOptions());
            n = header.Header.TextLength;
        }
        catch (Exception e) when (e is InvalidDataException or NotSupportedException)
        {
            Row("open_result", "corrupt", "", label, e.Message.Replace('\t', ' ').Replace('\n', ' '));
            return 3;
        }
        Row("open_n", n.ToString(Inv), "symbols", label);
        Row("open_level", level.ToString(), "", label, $"access {access}, text SHA-256 {(sha ? "compared" : "not compared")}");
        var reader = new RssReader();
        long before = reader.Read().RssAnonKb;
        using var sampler = new RssSampler(reader, 20);
        sampler.Start();
        long t0 = Stopwatch.GetTimestamp();
        string result = "ok";
        string message = "";
        SesaOpenReport? report = null;
        try
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = level, Access = access, VerifyTextSha256 = sha });
            report = file.Report;
        }
        catch (NotSupportedException e)
        {
            result = "not_applicable";
            message = e.Message;
        }
        catch (InvalidDataException e)
        {
            result = "corrupt";
            message = e.Message;
        }
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        RssStats stats = sampler.Stop();
        long after = reader.Read().RssAnonKb;
        if (report is { } r)
        {
            Row("open_header_s", r.HeaderSeconds.ToString("R", Inv), "s", label);
            Row("open_sections_s", r.SectionsSeconds.ToString("R", Inv), "s", label, $"{r.SectionBytesChecksummed} bytes checksummed, hardware CRC {(r.HardwareCrc ? "yes" : "no")}");
            Row("open_structure_s", r.StructureSeconds.ToString("R", Inv), "s", label, $"bitmap {r.BitmapBytes} bytes, buffers {r.BufferBytes} bytes");
        }
        Row("open_seconds", seconds.ToString("R", Inv), "s", label, "wall time of the open at the level, including the mapping");
        if (stats.Available)
        {
            Row("open_rssanon_before_mb", (before * 1024.0 / 1e6).ToString("R", Inv), "MB", label);
            Row("open_rssanon_peak_mb", (stats.MaxRssAnonKb * 1024.0 / 1e6).ToString("R", Inv), "MB", label, $"{stats.Samples} samples at 20 ms; NOT the gate of the build (320 MB), the structure pass holds a bitmap of n / 8 bytes");
            Row("open_rssanon_after_mb", (after * 1024.0 / 1e6).ToString("R", Inv), "MB", label);
        }
        Row("open_result", result, "", label, message.Replace('\t', ' ').Replace('\n', ' '));
        return result switch { "ok" => 0, "not_applicable" => 4, _ => 3 };
    }

    private static int Queries(Dictionary<string, string> o, string label)
    {
        var options = new HandoffQueryOptions
        {
            Container = Need(o, "container"),
            Input = Need(o, "input"),
            Seed = ulong.Parse(o.GetValueOrDefault("seed", "20261122"), Inv),
            Lengths = [.. o.GetValueOrDefault("lengths", "12,20,32").Split(',', StringSplitOptions.RemoveEmptyEntries).Select(x => int.Parse(x, Inv))],
            Present = int.Parse(o.GetValueOrDefault("present", "300"), Inv),
            Absent = int.Parse(o.GetValueOrDefault("absent", "200"), Inv),
            Cap = int.Parse(o.GetValueOrDefault("cap", "1000"), Inv),
            Engine = o.GetValueOrDefault("engine", "auto"),
        };
        if (options.Engine is not ("auto" or "true" or "false")) throw new ArgumentException("--engine must be auto, true or false");
        return HandoffQueryRun.Execute(options, (metric, value, notes) => Row(metric, value, "", label, notes));
    }
}
