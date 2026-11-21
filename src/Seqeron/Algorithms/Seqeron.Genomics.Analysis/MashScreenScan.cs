using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The pass of <c>mash screen</c> / <c>mash taxscreen</c> over the mixture (<c>CommandScreen.cpp</c> <c>run</c> + <c>hashSequence</c>, marbl/Mash v2.3): the records of the mixture files are read round-robin (one record of each open file in turn, a file leaves the
/// rotation at its end) and appended to a chunk as <c>'*'</c> + sequence (a record shorter than k is not appended); a chunk is hashed when the next record would make it longer than 2^20 characters and at the end. A window of the chunk is hashed when
/// it holds only characters of the alphabet, or, for an amino-acid query sketch, the chunk is translated in six frames (frame = offset in the CHUNK, so the frames of a record depend on the records before it, and a codon spanning the separator
/// translates to a stop) and the windows without a stop are hashed. A5-10b2 (B09 F214).
/// </summary>
internal static class MashScreenScan
{
    private const int ChunkSize = 1 << 20;

    // genetic code in the order of "ACGT" for the three positions (aaFromCodon of CommandScreen.cpp: only upper-case A, C, G, T; anything else is a stop '*')
    private const string Code = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

    private static int Base(byte b) => b switch { (byte)'A' => 0, (byte)'C' => 1, (byte)'G' => 2, (byte)'T' => 3, _ => -1 };

    private static byte AminoAcid(ReadOnlySpan<byte> codon)
    {
        int a = Base(codon[0]), b = Base(codon[1]), c = Base(codon[2]);
        return a < 0 || b < 0 || c < 0 ? (byte)'*' : (byte)Code[16 * a + 4 * b + c];
    }

    /// <summary>
    /// The chunks the program hashes and the value of its <c>count</c> (the number of loop turns that read a record or the end of a file, with at least one file still open: 0 means "Did not find sequence records in inputs").
    /// </summary>
    internal static List<string> Chunks(IReadOnlyList<MashInput> mixtures, int k, out long count)
    {
        var open = new List<(IReadOnlyList<MashRecord> Records, int Next)>();
        foreach (var m in mixtures)
            open.Add((m.Records, 0));
        var chunks = new List<string>();
        var input = new StringBuilder();
        count = 0;
        int it = 0;
        while (true)
        {
            int l;
            string seq = "";
            if (open.Count == 0)
            {
                l = 0;
            }
            else
            {
                var (records, next) = open[it];
                if (next >= records.Count)
                {
                    l = -1;                                    // end of the file: it leaves the rotation, the turn goes on with the file that followed it
                    open.RemoveAt(it);
                    if (it >= open.Count)
                        it = 0;
                }
                else
                {
                    seq = records[next].Sequence ?? "";
                    l = seq.Length;
                    open[it] = (records, next + 1);
                }
            }

            if (input.Length + (l >= k ? l + 1 : 0) > ChunkSize || open.Count == 0)
            {
                chunks.Add(input.ToString());
                input.Clear();
            }

            if (open.Count == 0)
                break;
            count++;
            if (l >= k)
                input.Append('*').Append(seq);
            it++;
            if (it >= open.Count)
                it = 0;
        }

        return chunks;
    }

    /// <summary>
    /// <c>hashSequence</c> of every chunk: <paramref name="onHash"/> is called with the hash of every window (every window of every frame with <paramref name="trans"/>); the bottom-<paramref name="sketchSize"/> set of all hashes is returned
    /// (the union of the per-chunk heaps of the program).
    /// </summary>
    internal static ulong[] Scan(IReadOnlyList<string> chunks, int k, int sketchSize, bool canonical, uint seed, MashAlphabet alphabet, bool trans, Action<ulong> onHash, CancellationToken cancellationToken)
    {
        var bottom = new SortedSet<ulong>();
        bool use64 = alphabet.Uses64BitHash(k);
        void Sink(ulong h)
        {
            ulong hash = use64 ? h : (uint)h;
            onHash(hash);
            if (bottom.Count >= sketchSize && hash >= bottom.Max)
                return;
            if (!bottom.Add(hash))
                return;
            if (bottom.Count > sketchSize)
                bottom.Remove(bottom.Max);
        }

        foreach (string chunk in chunks)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (!trans)
            {
                if (chunk.Length >= k)
                    KmerAnalyzer.StreamHashesInternal(chunk, k, canonical, seed, Sink, cancellationToken, alphabet);
                continue;
            }

            HashTranslated(chunk, k, seed, alphabet, Sink, cancellationToken);
        }

        return [.. bottom];
    }

    private static void HashTranslated(string chunk, int k, uint seed, MashAlphabet alphabet, Action<ulong> sink, CancellationToken cancellationToken)
    {
        int l = chunk.Length;
        var seq = new byte[l];
        for (int i = 0; i < l; i++)
        {
            char c = chunk[i];
            if (c > 127)
            {
                seq[i] = 0x80;                                  // no defined meaning in the program (signed char table index); never a base
                continue;
            }

            if (!alphabet.PreserveCase && c is >= 'a' and <= 'z')
                c = (char)(c - 32);
            seq[i] = (byte)c;
        }

        var seqRev = new byte[l];
        for (int i = 0; i < l; i++)
            seqRev[i] = MashAlphabet.Complement(seq[l - i - 1]);
        for (int f = 0; f < 6; f++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            int frame = f % 3;
            bool rev = f > 2;
            int lenTrans = (l - frame) / 3;
            var src = rev ? seqRev : seq;
            var trans = new byte[Math.Max(lenTrans, 0)];
            for (int a = 0; a < lenTrans; a++)
                trans[a] = AminoAcid(new ReadOnlySpan<byte>(src, frame + 3 * a, 3));
            int run = 0;
            for (int i = 0; i < trans.Length; i++)
            {
                if (trans[i] == (byte)'*')
                {
                    run = 0;
                    continue;
                }

                if (++run < k)
                    continue;
                run = k;
                sink(KmerAnalyzer.MurmurHash3X64_128(new ReadOnlySpan<byte>(trans, i - k + 1, k), seed).H1);
            }
        }
    }
}
