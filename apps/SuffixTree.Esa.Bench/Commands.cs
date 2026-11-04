using System.Diagnostics;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

internal static class Commands
{
    private const string TimerNote = "Stopwatch around TextPreparation.Prepare / SaisReference.Build / LcpReference.Build only; input already in memory";

    internal static void Build(Options o)
    {
        string text = Pipeline.LoadText(o.Input);
        int n = text.Length;
        Emit.Row("build", "input_chars", n, "chars");
        Emit.Row("build", "bounds_checked", EsaBuild.BoundsChecked ? 1 : 0, "bool", "1 = the Esa assembly was compiled with ESA_BOUNDS (checked accessors): such a build is never timed for acceptance");
        Emit.Row("build", "peak_rss_hwm_after_load", Process.StatusBytes("VmHWM"), "bytes", "VmHWM after loading the text as a managed string (2 B/char), before the build");
        long allocated0 = GC.GetTotalAllocatedBytes(precise: true);
        int gen0 = GC.CollectionCount(0), gen1 = GC.CollectionCount(1), gen2 = GC.CollectionCount(2);
        TimeSpan pause0 = GC.GetTotalPauseDuration();

        BuildResult r = Pipeline.Build(text, saisKernel: o.SaisKernel, lcpKernel: o.LcpKernel);

        long allocated = GC.GetTotalAllocatedBytes(precise: true) - allocated0;
        long hwm = Process.StatusBytes("VmHWM");
        Emit.Row("build", "build_prepare_s", r.PrepareSeconds, "s", "dense order-preserving remap of the UTF-16 text to byte symbols; " + TimerNote);
        Emit.Row("build", "build_sais_s", r.SaisSeconds, "s", $"SA-IS kernel {o.SaisKernel}");
        Emit.Row("build", "build_sa_s", r.SaSeconds, "s", "prepare + SA-IS: the SA from the string; libsais row = SA only from bytes");
        Emit.Row("build", "build_lcp_s", r.LcpSeconds, "s", $"LCP kernel {o.LcpKernel} (Phi/PLCP into the compact LCP)");
        Emit.Row("build", "build_lcp_jit_s", r.LcpJitSeconds, "s", "JIT compilation time of this thread inside the LCP step (part of build_lcp_s)");
        Emit.Row("build", "build_lcp_jit_methods", r.LcpJitMethods, "count", "methods compiled by this thread inside the LCP step");
        Emit.Row("build", "build_total_s", r.TotalSeconds, "s", "SA + LCP: the S1a gate metric (<= 3.3 s at 10 Mbp)");
        Emit.Row("build", "throughput", n / r.TotalSeconds, "chars/s", "SA + LCP");
        Emit.Row("build", "peak_rss_hwm_after_build", hwm, "bytes", "VmHWM of the whole process after the build (runtime, text, SA, LCP, transient arrays)");
        Emit.Row("build", "peak_rss_bytes_per_char", hwm / (double)n, "bytes/char", "gate: <= 24 B/char (PROGRESS-S1.md)");
        Emit.Row("build", "gc_allocated_bytes", allocated, "bytes", "GC.GetTotalAllocatedBytes delta during the build (every transient array counts)");
        Emit.Row("build", "gc_allocated_bytes_per_char", allocated / (double)n, "bytes/char");
        Emit.Row("build", "gc_collections_gen0", GC.CollectionCount(0) - gen0, "count");
        Emit.Row("build", "gc_collections_gen1", GC.CollectionCount(1) - gen1, "count");
        Emit.Row("build", "gc_collections_gen2", GC.CollectionCount(2) - gen2, "count");
        Emit.Row("build", "gc_pause_ms", (GC.GetTotalPauseDuration() - pause0).TotalMilliseconds, "ms", "total GC pause during the build");

        long saBytes = 4L * (n + 1);
        long lcpBytes = r.Lcp.PayloadBytes;
        long textBytes = ((long)n * (r.Prepared.Width == SymbolWidth.Byte ? 1 : 2)) + (2L * r.Prepared.AlphabetSize);
        Emit.Row("build", "index_bytes_sa", saBytes, "bytes", "int32 x (n + 1)");
        Emit.Row("build", "index_bytes_lcp", lcpBytes, "bytes", $"compact LCP kind {r.Lcp.Kind}, {r.Lcp.OverflowCount} escapes (8 B per escape)");
        Emit.Row("build", "index_bytes_text", textBytes, "bytes", $"remapped symbols ({r.Prepared.Width}) + alphabet table of {r.Prepared.AlphabetSize} code units");
        Emit.Row("build", "index_bytes", saBytes + lcpBytes + textBytes, "bytes", "SA + LCP + text (Core profile without child table and navigation pack)");
        Emit.Row("build", "index_bytes_per_char", (saBytes + lcpBytes + textBytes) / (double)n, "bytes/char", "gate: <= 16 B/char for the Full profile (S1b)");
        Emit.Row("build", "lcp_kind", r.Lcp.Kind.ToString(), "kind");
        Emit.Row("build", "lcp_escapes", r.Lcp.OverflowCount, "count");
        Emit.Row("build", "alphabet_size", r.Prepared.AlphabetSize, "count");

        // outside the timers: the result must be the libsais result
        string saHash = Hashing.Sa(r.Sa), lcpHash = Hashing.Lcp(r.Lcp, n);
        Emit.Row("build", "sa_sha256", saHash, "hex", "SHA-256 of SA[1..n] as little-endian int32 (= libsais SA)");
        Emit.Row("build", "lcp_sha256", lcpHash, "hex", "SHA-256 of LCP[1..n] as little-endian int32 (= libsais LCP)");
        if (o.Verify is not null)
        {
            var expected = Hashing.ExpectedVector(o.Verify, o.Name);
            if (expected is null) Emit.Row("build", "vector_match", "NA", "bool", $"{o.Name} is not listed in {Path.GetFileName(o.Verify)}");
            else
            {
                bool match = expected.Value.Sa == saHash && expected.Value.Lcp == lcpHash;
                Emit.Row("build", "vector_match", match ? 1 : 0, "bool", "SA and LCP hashes equal the libsais vectors of analysis/sa_vectors.tsv");
                if (!match) Environment.ExitCode = 3;
            }
        }
    }

    internal static void Steady(Options o)
    {
        int builds = o.Builds > 0 ? o.Builds : 5;
        if (builds < 4) throw new ArgumentException("--builds must be at least 4 (the last three are reported)");
        string text = Pipeline.LoadText(o.Input);
        for (int k = 1; k <= builds; k++)
        {
            BuildResult r = Pipeline.Build(text, saisKernel: o.SaisKernel, lcpKernel: o.LcpKernel);
            string note = $"build {k} of {builds} in one process";
            if (k <= 2) note += " (JIT warm-up, not reported as steady state)";
            Emit.Row(k > builds - 3 ? "steady" : "warmup", "build_sa_s", r.SaSeconds, "s", note);
            Emit.Row(k > builds - 3 ? "steady" : "warmup", "build_lcp_s", r.LcpSeconds, "s", note);
            Emit.Row(k > builds - 3 ? "steady" : "warmup", "build_total_s", r.TotalSeconds, "s", note);
        }
        Emit.Row("steady", "peak_rss_hwm_after_builds", Process.StatusBytes("VmHWM"), "bytes", $"VmHWM after {builds} builds");
    }

    internal static void Profile(Options o)
    {
        int builds = o.Builds > 0 ? o.Builds : 3;
        if (builds < 2) throw new ArgumentException("--builds must be at least 2 (first = cold, last = warm)");
        string text = Pipeline.LoadText(o.Input);
        for (int k = 1; k <= builds; k++)
        {
            var sais = new SaisTrace();
            var lcpTrace = new LcpTrace();
            BuildResult r = Pipeline.Build(text, sais, lcpTrace, o.SaisKernel, o.LcpKernel);
            if (k != 1 && k != builds) continue;
            string op = k == 1 ? "profile_cold" : "profile_warm";
            double tick = 1.0 / Stopwatch.Frequency;
            string note = $"instrumented build {k} of {builds}; exclusive phase times, the recursion is accounted by its own levels";
            Emit.Row(op, "prepare_s", r.PrepareSeconds, "s", note);
            foreach (SaisPhase phase in Enum.GetValues<SaisPhase>())
            {
                Emit.Row(op, $"sais_{phase}_s", sais.TotalTicks(phase) * tick, "s", note + "; all levels");
                Emit.Row(op, $"sais_{phase}_level0_s", sais.Ticks(0, phase) * tick, "s", note + "; level 0 only");
            }
            Emit.Row(op, "lcp_phi_s", lcpTrace.PhiTicks * tick, "s", note);
            Emit.Row(op, "lcp_plcp_s", lcpTrace.PlcpTicks * tick, "s", note);
            Emit.Row(op, "lcp_count_s", lcpTrace.CountTicks * tick, "s", note);
            Emit.Row(op, "lcp_output_s", lcpTrace.OutputTicks * tick, "s", note);
            Emit.Row(op, "sais_s", r.SaisSeconds, "s", note + "; whole SaisReference.Build");
            Emit.Row(op, "lcp_s", r.LcpSeconds, "s", note + "; whole LcpReference.Build");
            Emit.Row(op, "total_s", r.TotalSeconds, "s", note);
            Emit.Row(op, "sais_levels", sais.Depth, "count", note);
            for (int level = 0; level < Math.Min(sais.Depth, 4); level++)
                Emit.Row(op, $"sais_level{level}_length", sais.Levels[level].Length, "count", $"{note}; text length of recursion level {level} (m = {sais.Levels[level].LmsCount})");
        }
    }
}
