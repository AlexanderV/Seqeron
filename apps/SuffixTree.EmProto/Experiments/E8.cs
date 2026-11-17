using System.Diagnostics;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Verify;

namespace SuffixTree.EmProto.Experiments;

/// <summary>E8 of ADR-004 section 5: the hand-off package. Commands <c>e8-gen</c> (the streamed genome-like text) and <c>e8-verify</c> (the verifier of the outputs of <c>e6</c>).</summary>
internal static class E8
{
    /// <summary><c>e8-gen --n &lt;symbols&gt; --out &lt;file&gt; [--seed s] [--sums file] [--buffer-kb 1024] [--force]</c>: the genome-like text, streamed; the file appears under its name only when it is complete.</summary>
    internal static int Gen(Options o, Report r)
    {
        long n = E8Parsing.Count(o.Require("n"));
        string outPath = o.Require("out");
        string sumsPath = o.Get("sums") ?? outPath + ".sums";
        ulong seed = (ulong)o.Long("seed", (long)GenomeStream.DefaultSeed);
        int bufferBytes = o.Int("buffer-kb", 1024) * 1024;
        string c = o.Str("case", "e8");
        if (bufferBytes <= 0) throw new ArgumentException("--buffer-kb must be positive");
        if (File.Exists(outPath) && !o.Flag("force")) throw new ArgumentException($"'{outPath}' exists (use --force to overwrite)");
        string part = outPath + ".part";
        GenomeStreamResult res;
        try
        {
            using (var fs = new FileStream(part, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 16, FileOptions.SequentialScan))
            {
                res = GenomeStream.Generate(n, seed, block => fs.Write(block), bufferBytes);
                fs.Flush(true);
            }
            File.Move(part, outPath, true);
        }
        catch
        {
            if (File.Exists(part)) File.Delete(part);
            throw;
        }
        File.WriteAllText(sumsPath, GenomeStream.SumsText(res));
        r.Row("E8", c, "gen_n", res.N, "symbols", $"genome-like family, seed {res.Seed}, written to {Path.GetFileName(outPath)}");
        r.Row("E8", c, "gen_sha256", res.Sha256, "", "SHA-256 of the text file (also in " + Path.GetFileName(sumsPath) + ")");
        r.Row("E8", c, "gen_counts", string.Join(',', res.Counts), "", "counts of " + GenomeStream.Letters);
        r.Row("E8", c, "gen_seconds", res.Seconds, "s", "generation including the file write and the SHA-256");
        r.Row("E8", c, "gen_ns_per_symbol", res.Seconds * 1e9 / res.N, "ns", "gen_seconds / n");
        return 0;
    }

    /// <summary><c>e8-env</c>: what the .NET process sees (architecture, cores, memory the GC may use, page size, which vector instructions the program will use); rows <c>env_*</c>, for the preflight and the rows file.</summary>
    internal static int Env(Options o, Report r)
    {
        string c = o.Str("case", "env");
        r.Row("E8", c, "env_framework", RuntimeInformation.FrameworkDescription, "", "RuntimeInformation.FrameworkDescription");
        r.Row("E8", c, "env_os", RuntimeInformation.OSDescription, "", "RuntimeInformation.OSDescription");
        r.Row("E8", c, "env_os_arch", RuntimeInformation.OSArchitecture.ToString(), "", "architecture of the operating system (x64 / Arm64)");
        r.Row("E8", c, "env_process_arch", RuntimeInformation.ProcessArchitecture.ToString(), "", "architecture of this process (a difference to the OS architecture means emulation by the runtime)");
        r.Row("E8", c, "env_runtime_id", RuntimeInformation.RuntimeIdentifier, "", "runtime identifier of the framework");
        r.Row("E8", c, "env_cores", Environment.ProcessorCount, "n", "Environment.ProcessorCount (respects the container's CPU limit)");
        r.Row("E8", c, "env_gc_available_memory", GC.GetGCMemoryInfo().TotalAvailableMemoryBytes, "bytes", "memory the GC believes it may use (the container's memory limit when there is one)");
        r.Row("E8", c, "env_gc_hard_limit", Environment.GetEnvironmentVariable("DOTNET_GCHeapHardLimit") ?? "not set", "", "DOTNET_GCHeapHardLimit of this process");
        r.Row("E8", c, "env_page_size", Environment.SystemPageSize, "bytes", "system page size (4096 on x86-64 and on Docker Desktop's arm64 VM)");
        r.Row("E8", c, "env_little_endian", BitConverter.IsLittleEndian ? "yes" : "NO", "", "the output files are little endian by construction of the writers");
        r.Row("E8", c, "env_avx2", Avx2.IsSupported ? "yes" : "no", "", "System.Runtime.Intrinsics.X86.Avx2.IsSupported");
        r.Row("E8", c, "env_advsimd", AdvSimd.IsSupported ? "yes" : "no", "", "System.Runtime.Intrinsics.Arm.AdvSimd.IsSupported");
        r.Row("E8", c, "env_lce_path", Avx2.IsSupported ? "avx2 256-bit vectors" : "scalar 64-bit words", "", "how PackedText compares suffixes here (the scalar path gives the same results, tested; it is slower on tie-rich text)");
        return 0;
    }

    /// <summary>
    /// <c>e8-verify --prefix &lt;e6 prefix&gt; [--input text] [--threads 4] [--lce-cap 65536] [--stack-cap-m 64] [--scratch dir] [--no-hash] [--no-child] [--sums file] [--quiet]</c>: one row per check, the facts, the hashes;
    /// exit code 0 only when every listed check passed.
    /// </summary>
    internal static int Verify(Options o, Report r)
    {
        string prefix = o.Require("prefix");
        string c = o.Str("case", "e8");
        var options = new VerifyOptions
        {
            Prefix = prefix,
            InputPath = o.Get("input"),
            ScratchDir = o.Get("scratch"),
            Threads = o.Int("threads", 4),
            LceCap = o.Long("lce-cap", 65536),
            StackCapEntries = o.Long("stack-cap-m", 64) << 20,
            Hash = !o.Flag("no-hash"),
            Child = !o.Flag("no-child"),
            BlockRanks = o.Int("block-ranks", 1 << 20),
            Probes = o.Int("probes", 8),
            IoBufferBytes = o.Int("io-kb", 4096) * 1024,
        };
        long t0 = Stopwatch.GetTimestamp();
        Action<string>? progress = o.Flag("quiet") ? null : msg => Console.Error.WriteLine($"# verify {(Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency,8:F1} s  {msg}");
        VerifyResult res = OutputVerifier.Run(options, progress);
        foreach (CheckResult check in res.Checks)
        {
            string status = check.Status switch { CheckStatus.Pass => "pass", CheckStatus.Fail => "FAIL", _ => "skipped" };
            r.Row("E8", c, "verify_" + check.Name, status, "", check.Detail);
            r.Row("E8", c, $"verify_{check.Name}_s", check.Seconds, "s", "wall time of the check (the checks of the rank pass share one pass and show its time)");
        }
        foreach ((string key, string value) in res.Facts) r.Row("E8", c, "verify_fact_" + key, value, "", "number found by the verifier");
        foreach ((string ext, string sha) in res.Sha256) r.Row("E8", c, "sha256_" + ext.Replace('.', '_'), sha, "", $"SHA-256 of {Path.GetFileName(prefix)}.{ext}");
        if (o.Get("sums") is { } sums && res.Sha256.Count > 0) File.WriteAllLines(sums, new[] { "pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats" }.Where(res.Sha256.ContainsKey).Select(e => $"{e}={res.Sha256[e]}"));
        double wall = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        r.Row("E8", c, "verify_wall_s", wall, "s", $"whole verification, {options.Threads} threads, LCE cap {options.LceCap}");
        r.Row("E8", c, "verify_all", res.AllPassed ? "pass" : "FAIL", "", $"{res.Checks.Count(x => x.Status == CheckStatus.Pass)} of {res.Checks.Count} checks passed; --no-child / --no-hash omit checks (a run for the gate G1 uses neither)");
        return res.AllPassed ? 0 : 1;
    }
}

/// <summary>Parsing of counts with decimal suffixes (<c>k</c>, <c>m</c>, <c>g</c>), the unit of the E8 protocol (1 Gbp = 10^9 symbols).</summary>
internal static class E8Parsing
{
    internal static long Count(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        string t = text.Trim();
        decimal unit = 1;
        if (t.Length > 0 && char.IsLetter(t[^1]))
        {
            unit = char.ToLowerInvariant(t[^1]) switch
            {
                'k' => 1_000m,
                'm' => 1_000_000m,
                'g' => 1_000_000_000m,
                _ => throw new ArgumentException($"unknown unit in '{text}' (k, m or g)"),
            };
            t = t[..^1];
        }
        if (!decimal.TryParse(t, NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out decimal value)) throw new ArgumentException($"'{text}' is not a count");
        decimal total = value * unit;
        if (total < 1 || total != decimal.Truncate(total) || total > long.MaxValue) throw new ArgumentException($"'{text}' is not a positive whole number of symbols");
        return (long)total;
    }
}
