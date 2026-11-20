using System.Globalization;
using System.Text;

namespace ReadBench;

/// <summary>
/// The texts of the benchmark definition of work package 1b.9b (BENCHMARK-1b9b.md): deterministic generators, one per text kind, all over SplitMix64, written as one byte per symbol (code units below 256).
/// usage: ReadBench gen &lt;kind&gt; &lt;length&gt; &lt;seed&gt; &lt;output file&gt;   with kind one of sigma5, prot-uniform, prot-real, text90, bytes256 (the DNA text is the S0 input dna_10m, a file that already exists).
/// </summary>
internal static class Generate
{
    internal const string Dna = "ACGT";
    internal const string Sigma5 = "ACGTN";
    internal const string Amino = "ACDEFGHIKLMNPQRSTVWY";

    // composition of the UniProtKB/Swiss-Prot release statistics, in percent, in the order of Amino
    private static readonly double[] AminoPercent =
    [
        8.25, 1.37, 5.45, 6.75, 3.86, 7.07, 2.27, 5.96, 5.84, 9.66, 2.42, 4.06, 4.70, 3.93, 5.53, 6.56, 5.34, 6.87, 1.08, 2.92,
    ];

    internal static int Run(string[] args)
    {
        if (args.Length != 5)
        {
            Console.Error.WriteLine("usage: ReadBench gen <sigma5|prot-uniform|prot-real|text90|bytes256> <length> <seed> <output>");
            return 2;
        }
        string kind = args[1];
        int n = int.Parse(args[2], CultureInfo.InvariantCulture);
        var rng = new Rng(ulong.Parse(args[3], CultureInfo.InvariantCulture));
        byte[] text = kind switch
        {
            "sigma5" => Uniform(n, Sigma5, rng),
            "prot-uniform" => Uniform(n, Amino, rng),
            "prot-real" => ProteinRealistic(n, rng),
            "text90" => Text90(n, rng),
            "bytes256" => Bytes256(n, rng),
            _ => throw new ArgumentException($"unknown kind '{kind}'"),
        };
        File.WriteAllBytes(args[4], text);
        int distinct = text.Distinct().Count();
        Console.WriteLine($"{kind}\t{n}\t{distinct} distinct symbols\t{args[4]}");
        return 0;
    }

    private static byte[] Uniform(int n, string alphabet, Rng rng)
    {
        var text = new byte[n];
        for (int i = 0; i < n; i++) text[i] = (byte)alphabet[rng.Int(alphabet.Length)];
        return text;
    }

    private static byte[] Bytes256(int n, Rng rng)
    {
        var text = new byte[n];
        for (int i = 0; i < n; i += 8)
        {
            ulong v = rng.Next();
            for (int k = 0; k < 8 && i + k < n; k++) text[i + k] = (byte)(v >> (8 * k));
        }
        return text;
    }

    /// <summary>
    /// A proteome-like text over the 20 amino-acid letters: proteins of exponentially distributed length (60 + mean 300), each starting with M; a residue is drawn from the Swiss-Prot composition; a segment is, with
    /// probability 0.35, a copy of an earlier stretch of the text (40-300 residues) with a per-residue substitution rate of 2, 10, 30 or 50 % and indels at 1 % (paralogous domains), with probability 0.04 a low-complexity run
    /// (a residue 5-30 times, or a unit of 2-6 residues repeated 3-15 times), and otherwise 10-60 fresh residues. The proteins are concatenated without a separator.
    /// </summary>
    private static byte[] ProteinRealistic(int n, Rng rng)
    {
        double[] cumulative = new double[AminoPercent.Length];
        double sum = 0;
        for (int i = 0; i < cumulative.Length; i++) cumulative[i] = sum += AminoPercent[i];
        byte Residue()
        {
            double u = rng.Unit() * sum;
            int i = Array.BinarySearch(cumulative, u);
            if (i < 0) i = ~i;
            return (byte)Amino[Math.Min(i, Amino.Length - 1)];
        }
        double[] rates = [0.02, 0.10, 0.30, 0.50];
        var text = new byte[n];
        int length = 0;
        while (length < n)
        {
            int protein = 60 + (int)(-Math.Log(1 - rng.Unit()) * 300);
            int end = (int)Math.Min(n, (long)length + protein);
            text[length++] = (byte)'M';
            while (length < end)
            {
                double kind = rng.Unit();
                if (kind < 0.35 && length > 400)
                {
                    int copy = 40 + rng.Int(261);
                    int from = rng.Int(length - copy);
                    double rate = rates[rng.Int(rates.Length)];
                    for (int k = 0; k < copy && length < end; k++)
                    {
                        double u = rng.Unit();
                        if (u < 0.005) continue; // a deletion
                        if (u < 0.01) text[length++] = Residue(); // an insertion
                        if (length >= end) break;
                        text[length++] = rng.Unit() < rate ? Residue() : text[from + k];
                    }
                }
                else if (kind < 0.39)
                {
                    if (rng.Int(2) == 0)
                    {
                        byte r = Residue();
                        int run = 5 + rng.Int(26);
                        for (int k = 0; k < run && length < end; k++) text[length++] = r;
                    }
                    else
                    {
                        int unit = 2 + rng.Int(5);
                        var u = new byte[unit];
                        for (int k = 0; k < unit; k++) u[k] = Residue();
                        int copies = 3 + rng.Int(13);
                        for (int c = 0; c < copies; c++)
                            for (int k = 0; k < unit && length < end; k++) text[length++] = u[k];
                    }
                }
                else
                {
                    int fresh = 10 + rng.Int(51);
                    for (int k = 0; k < fresh && length < end; k++) text[length++] = Residue();
                }
            }
        }
        return text;
    }

    /// <summary>
    /// A natural-language-like text over about 90 printable ASCII symbols (letters of both cases, digits, space and punctuation; the generator prints the count): a vocabulary of 20,000 words (letter strings of 2-12 letters with English letter frequencies), drawn with Zipf probability 1 / rank, in sentences of 5-20 words
    /// (the first word capitalised, the sentence closed by one of ". ? ! ", words separated by a space or by a punctuation mark and a space), 5 % of the tokens a number of 1-6 digits, and 1 % a token of 1-4 rare marks.
    /// </summary>
    private static byte[] Text90(int n, Rng rng)
    {
        const string Letters = "etaoinshrdlcumwfgypbvkjxqz";
        double[] letterWeights = [12.7, 9.1, 8.2, 7.5, 7.0, 6.7, 6.3, 6.1, 6.0, 4.3, 4.0, 2.8, 2.8, 2.4, 2.4, 2.2, 2.0, 2.0, 1.9, 1.5, 1.0, 0.8, 0.15, 0.15, 0.1, 0.07];
        double[] letterCumulative = new double[letterWeights.Length];
        double lsum = 0;
        for (int i = 0; i < letterWeights.Length; i++) letterCumulative[i] = lsum += letterWeights[i];
        char Letter()
        {
            int i = Array.BinarySearch(letterCumulative, rng.Unit() * lsum);
            if (i < 0) i = ~i;
            return Letters[Math.Min(i, Letters.Length - 1)];
        }
        const int Words = 20000;
        var vocabulary = new string[Words];
        for (int w = 0; w < Words; w++)
        {
            int len = 2 + Math.Min(10, (int)(-Math.Log(1 - rng.Unit()) * 3.0));
            var chars = new char[len];
            for (int k = 0; k < len; k++) chars[k] = Letter();
            vocabulary[w] = new string(chars);
        }
        double[] zipf = new double[Words];
        double zsum = 0;
        for (int w = 0; w < Words; w++) zipf[w] = zsum += 1.0 / (w + 1);
        const string Rare = "!#$%&()*+/:;<=>?@[]_{}\"'-,.0123456789";
        var sb = new StringBuilder(n + 64);
        while (sb.Length < n)
        {
            int words = 5 + rng.Int(16);
            for (int w = 0; w < words; w++)
            {
                double kind = rng.Unit();
                if (kind < 0.05)
                {
                    int digits = 1 + rng.Int(6);
                    for (int d = 0; d < digits; d++) sb.Append((char)('0' + rng.Int(10)));
                }
                else if (kind < 0.06)
                {
                    int marks = 1 + rng.Int(4);
                    for (int d = 0; d < marks; d++) sb.Append(Rare[rng.Int(Rare.Length)]);
                }
                else
                {
                    int i = Array.BinarySearch(zipf, rng.Unit() * zsum);
                    if (i < 0) i = ~i;
                    string word = vocabulary[Math.Min(i, Words - 1)];
                    if (w == 0) sb.Append(char.ToUpperInvariant(word[0])).Append(word, 1, word.Length - 1);
                    else if (rng.Int(40) == 0) sb.Append(char.ToUpperInvariant(word[0])).Append(word, 1, word.Length - 1);
                    else sb.Append(word);
                }
                if (w < words - 1) sb.Append(rng.Int(8) == 0 ? ", " : " ");
            }
            sb.Append(rng.Int(3) switch { 0 => ". ", 1 => "? ", _ => "! " });
        }
        var text = new byte[n];
        for (int i = 0; i < n; i++) text[i] = (byte)sb[i];
        return text;
    }

    /// <summary>SplitMix64 (the generator of the harness).</summary>
    private sealed class Rng(ulong seed)
    {
        private ulong _state = seed;

        internal ulong Next()
        {
            unchecked
            {
                _state += 0x9E3779B97F4A7C15UL;
                ulong z = _state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
                return z ^ (z >> 31);
            }
        }

        internal int Int(int bound) => (int)(Next() % (ulong)bound);

        internal double Unit() => (Next() >> 11) * (1.0 / (1UL << 53));
    }
}
