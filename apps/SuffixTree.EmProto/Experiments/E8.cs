using System.Globalization;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;

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

    internal static int Verify(Options o, Report r) => throw new NotImplementedException();
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
