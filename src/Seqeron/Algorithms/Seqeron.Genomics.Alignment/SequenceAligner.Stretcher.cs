using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// EMBOSS 6.6 stretcher (audit A2-28, B09 F121): the Myers-Miller (1988, CABIOS 4:11) linear-space global alignment with affine gaps, written as emboss/stretcher.c does - the same
// midpoint search with its own tie rules (CC[j] + RR[j] >= midc taking the later j when the sum is larger or a gap-closed forward cell meets an open reverse cell, then the type-2 scan
// from the right), the same edit-script macros (STRETCHERDEL / STRETCHERINS / STRETCHERREP) and the same script-to-rows conversion (stretcher_Calcons), so that the rows are the rows of the program.
public static partial class SequenceAligner
{
    /// <summary>
    /// EMBOSS <c>stretcher</c>: global alignment by the Myers-Miller linear-space algorithm (Myers &amp; Miller 1988) with integer scores. A gap of length k costs
    /// <c>gapOpen + (k - 1) * gapExtend</c> (EMBOSS <c>-gapopen</c> / <c>-gapextend</c>, the same numbers as needle's); the alignment is end to end with end gaps charged like any gap
    /// (stretcher has no end-gap options), the sequences are upper-cased, a residue that is not in the matrix scores 0 against everything (EMBOSS "bad character", code 0).
    /// The rows are the rows the installed stretcher 6.6.0 prints (its midpoint tie rules differ from <see cref="GlobalAlignAffine(string, string, MatrixScoring)"/>; the score is the same).
    /// </summary>
    /// <param name="sequence1">First sequence (<c>-asequence</c>).</param>
    /// <param name="sequence2">Second sequence (<c>-bsequence</c>).</param>
    /// <param name="matrix">Residue-pair scores (<c>-datafile</c>); default EDNAFULL for nucleotide sequences, BLOSUM62 (EBLOSUM62) otherwise (EMBOSS picks by sequence type).</param>
    /// <param name="gapOpen">Cost of a gap of length one (<c>-gapopen</c>, default 16 for nucleotide sequences and 12 for protein, at least 0).</param>
    /// <param name="gapExtend">Cost of each further gap position (<c>-gapextend</c>, default 4 for nucleotide sequences and 2 for protein, at least 0).</param>
    /// <returns>The global alignment: both rows, the integer score (<see cref="AlignmentResult.Score"/>) and the full coordinates.</returns>
    public static AlignmentResult StretcherAlign(string sequence1, string sequence2, SubstitutionMatrix? matrix = null, int? gapOpen = null, int? gapExtend = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        string s1 = EmbossSequenceReader.Clean(sequence1).ToUpperInvariant();
        string s2 = EmbossSequenceReader.Clean(sequence2).ToUpperInvariant();
        bool nucleotide = EmbossSequenceReader.IsNucleotide(s1) && EmbossSequenceReader.IsNucleotide(s2);
        matrix ??= nucleotide ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        int go = gapOpen ?? (nucleotide ? 16 : 12);
        int ge = gapExtend ?? (nucleotide ? 4 : 2);
        if (go < 0)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "The gap opening penalty must be at least 0 (EMBOSS minimum).");
        if (ge < 0)
            throw new ArgumentOutOfRangeException(nameof(gapExtend), "The gap extension penalty must be at least 0 (EMBOSS minimum).");
        RequireStretcherRange(matrix, go, ge, s1.Length, s2.Length);
        var engine = new StretcherEngine(s1, s2, matrix, go - ge, ge);
        var (score, rows1, rows2) = engine.Run();
        return new AlignmentResult(rows1, rows2, score, AlignmentType.Global, 0, 0, s1.Length, s2.Length);
    }

    // static state of stretcher.c made instance state: g = G, hh = H, m = g + hh
    private sealed class StretcherEngine
    {
        private readonly string _a;
        private readonly string _b;
        private readonly int _g;
        private readonly int _hh;
        private readonly int _m;
        private readonly int[][] _subRows;     // _subRows[i] (1-based residue i of A): scores against every B code
        private readonly int[] _bCode;         // 1-based codes of B
        private readonly int[] _cc;
        private readonly int[] _dd;
        private readonly int[] _rr;
        private readonly int[] _ss;
        private readonly int[] _script;
        private int _sapp;
        private int _last;

        public StretcherEngine(string a, string b, SubstitutionMatrix matrix, int g, int hh)
        {
            _a = a;
            _b = b;
            _g = g;
            _hh = hh;
            _m = g + hh;
            string alphabet = matrix.Alphabet;
            int unknown = alphabet.Length;
            var rowOfChar = new Dictionary<char, int[]>();
            _bCode = new int[b.Length + 1];
            for (int j = 1; j <= b.Length; j++)
            {
                int k = alphabet.IndexOf(b[j - 1], StringComparison.Ordinal);
                _bCode[j] = k < 0 ? unknown : k;
            }

            _subRows = new int[a.Length + 1][];
            for (int i = 1; i <= a.Length; i++)
            {
                char c = a[i - 1];
                if (!rowOfChar.TryGetValue(c, out var row))
                {
                    row = new int[alphabet.Length + 1];
                    for (int k = 0; k < alphabet.Length; k++)
                        row[k] = matrix.TryGetScore(c, alphabet[k], out int s) ? s : 0;
                    rowOfChar[c] = row;
                }

                _subRows[i] = row;
            }

            _cc = new int[b.Length + 1];
            _dd = new int[b.Length + 1];
            _rr = new int[b.Length + 1];
            _ss = new int[b.Length + 1];
            _script = new int[a.Length + b.Length + 2];
        }

        private int Gap(int k) => k <= 0 ? 0 : _g + (_hh * k);

        private void Del(int k)
        {
            if (_last < 0)
            {
                _script[_sapp - 1] -= k;
                _last = _script[_sapp - 1];
            }
            else
            {
                _script[_sapp++] = -k;
                _last = -k;
            }
        }

        private void Ins(int k)
        {
            if (_last < 0)
            {
                _script[_sapp - 1] = k;
                _script[_sapp++] = _last;
            }
            else
            {
                _script[_sapp++] = k;
                _last = k;
            }
        }

        private void Rep()
        {
            _script[_sapp++] = 0;
            _last = 0;
        }

        public (int Score, string Row1, string Row2) Run()
        {
            _sapp = 0;
            _last = 0;
            int score = Align(0, 0, _a.Length, _b.Length, -_g, -_g);
            // stretcher_Calcons
            var r1 = new StringBuilder(_a.Length + _b.Length);
            var r2 = new StringBuilder(_a.Length + _b.Length);
            int i0 = 0, i1 = 0, op = 0, rp = 0;
            while (i0 < _a.Length || i1 < _b.Length)
            {
                if (op == 0 && _script[rp] == 0)
                {
                    op = _script[rp++];
                    r1.Append(_a[i0++]);
                    r2.Append(_b[i1++]);
                }
                else
                {
                    if (op == 0)
                        op = _script[rp++];
                    if (op > 0)
                    {
                        r1.Append('-');
                        r2.Append(_b[i1++]);
                        op--;
                    }
                    else
                    {
                        r1.Append(_a[i0++]);
                        r2.Append('-');
                        op++;
                    }
                }
            }

            return (score, r1.ToString(), r2.ToString());
        }

        // stretcher_Align over A[ao+1..ao+M], B[bo+1..bo+N]
        private int Align(int ao, int bo, int bigM, int n, int tb, int te)
        {
            if (n <= 0)
            {
                if (bigM > 0)
                    Del(bigM);
                return -Gap(bigM);
            }

            if (bigM <= 1)
            {
                if (bigM <= 0)
                {
                    Ins(n);
                    return -Gap(n);
                }

                if (tb < te)
                    tb = te;
                int midc0 = (tb - _hh) - Gap(n);
                int midj0 = 0;
                var wa0 = _subRows[ao + 1];
                for (int j = 1; j <= n; j++)
                {
                    int c0 = -Gap(j - 1) + wa0[_bCode[bo + j]] - Gap(n - j);
                    if (c0 > midc0)
                    {
                        midc0 = c0;
                        midj0 = j;
                    }
                }

                if (midj0 == 0)
                {
                    Ins(n);
                    Del(1);
                }
                else
                {
                    if (midj0 > 1)
                        Ins(midj0 - 1);
                    Rep();
                    if (midj0 < n)
                        Ins(n - midj0);
                }

                return midc0;
            }

            int midi = bigM / 2;
            _cc[0] = 0;
            int t = -_g;
            for (int j = 1; j <= n; j++)
            {
                t -= _hh;
                _cc[j] = t;
                _dd[j] = t - _g;
            }

            t = tb;
            int c, d, e, s;
            for (int i = 1; i <= midi; i++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per row of the forward / reverse sweep
                s = _cc[0];
                t -= _hh;
                _cc[0] = c = t;
                e = t - _g;
                var wa = _subRows[ao + i];
                for (int j = 1; j <= n; j++)
                {
                    c -= _m;
                    e -= _hh;
                    if (c > e)
                        e = c;
                    c = _cc[j] - _m;
                    d = _dd[j] - _hh;
                    if (c > d)
                        d = c;
                    c = s + wa[_bCode[bo + j]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    s = _cc[j];
                    _cc[j] = c;
                    _dd[j] = d;
                }
            }

            _dd[0] = _cc[0];
            _rr[n] = 0;
            t = -_g;
            for (int j = n - 1; j >= 0; j--)
            {
                t -= _hh;
                _rr[j] = t;
                _ss[j] = t - _g;
            }

            t = te;
            for (int i = bigM - 1; i >= midi; i--)
            {
                AlignmentCancellation.ThrowIfRequested();
                s = _rr[n];
                t -= _hh;
                _rr[n] = c = t;
                e = t - _g;
                var wa = _subRows[ao + i + 1];
                for (int j = n - 1; j >= 0; j--)
                {
                    c -= _m;
                    e -= _hh;
                    if (c > e)
                        e = c;
                    c = _rr[j] - _m;
                    d = _ss[j] - _hh;
                    if (c > d)
                        d = c;
                    c = s + wa[_bCode[bo + j + 1]];
                    if (e > c)
                        c = e;
                    if (d > c)
                        c = d;
                    s = _rr[j];
                    _rr[j] = c;
                    _ss[j] = d;
                }
            }

            _ss[n] = _rr[n];
            int midc = _cc[0] + _rr[0];
            int midj = 0;
            int type = 1;
            for (int j = 0; j <= n; j++)
            {
                c = _cc[j] + _rr[j];
                if (c >= midc && (c > midc || (_cc[j] != _dd[j] && _rr[j] == _ss[j])))
                {
                    midc = c;
                    midj = j;
                }
            }

            for (int j = n; j >= 0; j--)
            {
                c = _dd[j] + _ss[j] + _g;
                if (c > midc)
                {
                    midc = c;
                    midj = j;
                    type = 2;
                }
            }

            if (type == 1)
            {
                Align(ao, bo, midi, midj, tb, -_g);
                Align(ao + midi, bo + midj, bigM - midi, n - midj, -_g, te);
            }
            else
            {
                Align(ao, bo, midi - 1, midj, tb, 0);
                Del(2);
                Align(ao + midi + 1, bo + midj, bigM - midi - 1, n - midj, 0, te);
            }

            return midc;
        }
    }
}
