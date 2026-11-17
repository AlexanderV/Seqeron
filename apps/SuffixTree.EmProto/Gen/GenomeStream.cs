using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Gen;

/// <summary>What a streamed generation produced: the length, the seed, the SHA-256 of the whole output and the count of every letter of <see cref="GenomeStream.Letters"/>.</summary>
internal sealed record GenomeStreamResult(long N, ulong Seed, string Sha256, long[] Counts, double Seconds);

/// <summary>
/// E8 (ADR-004 section 5, row E8): the genome-like text of <see cref="Families.GenomeLike"/> produced as a stream of ASCII bytes with a bounded buffer, for any n (1 to 3x10^9 and beyond), byte for byte equal to
/// <c>Families.ToAscii(Families.GenomeLike(n, seed))</c> wherever that can be built in RAM, and identical for any buffer size.
/// </summary>
internal static class GenomeStream
{
    /// <summary>The seed of the family in <c>FamilyLoader</c> (<c>genome:N</c>), i.e. of the inputs of E3 to E7.</summary>
    internal const ulong DefaultSeed = 20240002;

    /// <summary>The letters of the codes 0..4 (A C G N T), as <see cref="Families.GenomeLike"/>.</summary>
    internal const string Letters = "ACGNT";

    /// <summary>Streams the text to <paramref name="sink"/> in blocks of at most <paramref name="bufferBytes"/> bytes.</summary>
    internal static GenomeStreamResult Generate(long n, ulong seed, Action<ReadOnlySpan<byte>> sink, int bufferBytes = 1 << 20)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(n, 0);
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bufferBytes, 0);
        ArgumentNullException.ThrowIfNull(sink);
        long t0 = Stopwatch.GetTimestamp();
        using var emitter = new Emitter(sink, bufferBytes);
        var rng = new SplitMix64(seed);
        var families = new byte[24][];
        for (int f = 0; f < families.Length; f++)
        {
            families[f] = new byte[300 + rng.NextInt(5700)];
            for (int i = 0; i < families[f].Length; i++) families[f][i] = (byte)Four(rng);
        }
        var monomer = new byte[171];
        for (int i = 0; i < monomer.Length; i++) monomer[i] = (byte)Four(rng);
        var nRuns = new List<(long Start, long Length)>
        {
            (n / 20, n / 200),
            (n / 3, n / 1000),
            (n / 2, n / 2000),
            (n * 3 / 5, n / 1000),
            (n * 4 / 5, n / 2000),
            (n - (n / 2000), n / 2000),
        };
        nRuns.Sort((a, b) => a.Start.CompareTo(b.Start));
        long satMean = Math.Max(2000, n / 400) + (Math.Max(1000, n / 800) / 2);
        double wBackground = 0.55 / 1000, wRepeat = 0.40 / 1500, wSatellite = 0.04 / satMean;
        double total = wBackground + wRepeat + wSatellite;
        var hor = new byte[171 * 12];
        long pos = 0;
        int nextRun = 0;
        while (pos < n)
        {
            if (nextRun < nRuns.Count && pos >= nRuns[nextRun].Start)
            {
                long end = Math.Min(n, pos + nRuns[nextRun].Length);
                emitter.Fill((byte)'N', end - pos);
                pos = end;
                nextRun++;
                continue;
            }
            long limit = nextRun < nRuns.Count ? Math.Min(n, nRuns[nextRun].Start) : n;
            double pick = rng.NextDouble() * total;
            if (pick < wBackground)
            {
                long len = Math.Min(limit - pos, 200 + rng.NextInt(1600));
                for (long i = 0; i < len; i++) emitter.Put(Code((byte)Four(rng)));
                pos += len;
            }
            else if (pick < wBackground + wRepeat)
            {
                long len = Math.Min(limit - pos, 500 + rng.NextInt(2000));
                RepeatCopy(emitter, len, families, rng);
                pos += len;
            }
            else
            {
                long len = Math.Min(limit - pos, Math.Max(2000, n / 400) + rng.NextInt((int)Math.Max(1000, n / 800)));
                for (int m = 0; m < 12; m++)
                {
                    for (int i = 0; i < 171; i++) hor[(m * 171) + i] = rng.NextDouble() < 0.25 ? (byte)Four(rng) : monomer[i];
                }
                int h = 0;
                for (long k = 0; k < len; k++)
                {
                    byte v = hor[h];
                    if (++h == hor.Length) h = 0;
                    if (rng.NextDouble() < 0.01) v = (byte)Four(rng);
                    emitter.Put(Code(v));
                }
                pos += len;
            }
        }
        (string sha, long[] counts) = emitter.Finish();
        if (emitter.Written != n) throw new InvalidOperationException($"the generator wrote {emitter.Written} symbols, expected {n}");
        return new GenomeStreamResult(n, seed, sha, counts, (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency);
    }

    /// <summary>The text of the checksum file (<c>key=value</c> lines).</summary>
    internal static string SumsText(GenomeStreamResult r) => string.Create(CultureInfo.InvariantCulture,
        $"format=e8-gen-v1\ngenerator=genome-like family of ADR-004 (Families.GenomeLike), streamed\nn={r.N}\nseed={r.Seed}\nalphabet={Letters}\ncounts={string.Join(',', r.Counts)}\nsha256={r.Sha256}\n");

    /// <summary><c>NextInt(4)</c> without the division: the same value (<c>x % 4 == x &amp; 3</c> for unsigned x), the same single draw.</summary>
    private static int Four(SplitMix64 rng) => (int)(rng.Next() & 3);

    /// <summary>Final code of a base 0..3 (A C G T) as <see cref="Families"/>: T moves to 4 so that N (3) sorts between G and T; returned as the ASCII letter.</summary>
    private static byte Code(byte b) => (byte)Letters[b == 3 ? 4 : b];

    private static void RepeatCopy(Emitter emitter, long length, byte[][] families, SplitMix64 rng)
    {
        byte[] f = families[rng.NextInt(families.Length)];
        double divergence = 0.10 + (0.05 * rng.NextDouble());
        int at = rng.NextInt(Math.Max(1, f.Length - (int)Math.Min(length, f.Length) + 1));
        long written = 0;
        int src = at;
        while (written < length)
        {
            if (src >= f.Length) src = 0;
            double r = rng.NextDouble();
            if (r < 0.005)
            {
                src++; // deletion in the copy
                continue;
            }
            if (r < 0.01)
            {
                emitter.Put(Code((byte)Four(rng))); // insertion
                written++;
                continue;
            }
            byte v = f[src++];
            if (rng.NextDouble() < divergence) v = (byte)((v + 1 + rng.NextInt(3)) & 3);
            emitter.Put(Code(v));
            written++;
        }
    }

    /// <summary>The output buffer: collects bytes, hands full blocks to the sink, hashes and counts what it hands over.</summary>
    private sealed class Emitter : IDisposable
    {
        private readonly Action<ReadOnlySpan<byte>> _sink;
        private readonly byte[] _buffer;
        private readonly IncrementalHash _hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        private readonly long[] _histogram = new long[256];
        private int _fill;

        internal Emitter(Action<ReadOnlySpan<byte>> sink, int bufferBytes)
        {
            _sink = sink;
            _buffer = new byte[bufferBytes];
        }

        /// <summary>Symbols handed to <see cref="Put"/> and <see cref="Fill"/> so far.</summary>
        internal long Written { get; private set; }

        internal void Put(byte letter)
        {
            _buffer[_fill++] = letter;
            Written++;
            if (_fill == _buffer.Length) Flush();
        }

        internal void Fill(byte letter, long count)
        {
            Written += count;
            while (count > 0)
            {
                int take = (int)Math.Min(count, _buffer.Length - _fill);
                _buffer.AsSpan(_fill, take).Fill(letter);
                _fill += take;
                count -= take;
                if (_fill == _buffer.Length) Flush();
            }
        }

        internal (string Sha256, long[] Counts) Finish()
        {
            Flush();
            string sha = Convert.ToHexStringLower(_hash.GetHashAndReset());
            return (sha, [.. Letters.Select(c => _histogram[c])]);
        }

        public void Dispose() => _hash.Dispose();

        private void Flush()
        {
            if (_fill == 0) return;
            ReadOnlySpan<byte> block = _buffer.AsSpan(0, _fill);
            foreach (byte b in block) _histogram[b]++;
            _hash.AppendData(block);
            _sink(block);
            _fill = 0;
        }
    }
}
