// Statistics of the real sections of a container that decide whether a delta form pays (ADR-005 section 5): the child cells (the owner's candidate "delta from the index with escapes"),
// PSI (increasing inside each first-symbol bucket) and LCP escapes. usage: stats --container <file>
internal static unsafe class Stats
{
    internal static int Run(Dictionary<string, string> o)
    {
        using var c = new Container(Program.Get(o, "container"));
        long n = c.N;
        Console.WriteLine($"n={n} sigma={c.Sigma} posWidth={c.PosWidth}");
        Section child = c.Sections[6];
        uint* cld = (uint*)c.At(6);
        long[] hist = new long[4 * 4]; // [kind][bucket]: 0 = |d| <= 127, 1 = <= 32767, 2 = <= 2^31, 3 = more
        long empty = 0;
        for (long k = 0; k <= n; k++)
        {
            uint cell = cld[k];
            int kind = (int)(cell >> 30);
            long v = cell & ((1u << 30) - 1);
            if (kind == 0) { empty++; continue; }
            long d = Math.Abs(v - k);
            int b = d <= 127 ? 0 : d <= 32767 ? 1 : d <= (1L << 31) ? 2 : 3;
            hist[(kind * 4) + b]++;
        }
        long total = n + 1;
        Console.WriteLine($"child cells: {total}, empty {empty} ({100.0 * empty / total:F2} %), codec {child.Codec}");
        string[] names = { "-", "next", "up", "down" };
        for (int kind = 1; kind <= 3; kind++)
        {
            long sum = hist[(kind * 4)..((kind * 4) + 4)].Sum();
            Console.WriteLine($"  {names[kind],-5} {sum,12} cells: |v-k|<=127 {100.0 * hist[kind * 4] / total:F2} % of all cells, <=32767 {100.0 * hist[(kind * 4) + 1] / total:F2} %, larger {100.0 * (hist[(kind * 4) + 2] + hist[(kind * 4) + 3]) / total:F2} %");
        }
        long inByte = hist[4] + hist[8] + hist[12];
        Console.WriteLine($"  cells whose value is within 127 of the index (a signed byte): {100.0 * inByte / total:F2} % of all cells; the rest {100.0 * (total - empty - inByte) / total:F2} % would need an escape (of which empty cells need none)");
        if (c.Has(12))
        {
            uint* psi = (uint*)c.At(12);
            long le255 = 0, bigger = 0, decreases = 0;
            for (long i = 1; i <= n; i++)
            {
                uint d = psi[i] - psi[i - 1];
                if (psi[i] < psi[i - 1]) { decreases++; continue; }
                if (d <= 255) le255++; else bigger++;
            }
            Console.WriteLine($"psi: increments <= 255 {100.0 * le255 / total:F3} %, larger {100.0 * bigger / total:F3} %, decreases (bucket starts) {decreases} ({100.0 * decreases / total:F6} %)");
        }
        if (c.Has(5)) Console.WriteLine($"lcp overflow entries {c.Sections[5].Count} ({100.0 * c.Sections[5].Count / total:F4} % of ranks), lcp codec {c.Sections[4].Codec}");
        return 0;
    }
}
