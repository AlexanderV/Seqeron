// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;

/// <summary>The consensus and the statistics of an alignment of any number of rows: <c>alignConsStats</c> of EMBOSS 6.6.0 <c>ajax/core/ajalign.c</c> (float weights, plurality 51 %, the matrix EDNAFULL / EBLOSUM62).</summary>
public static class EmbossAlignmentConsensus
{
    private const float PluralityPercent = 51.0f;
    private const float Epsilon = 1.192e-6f;

    /// <summary>The result of <c>alignConsStats</c>: the consensus string (upper case = identical, lower case = similar, <c>#</c> = no plurality of residues, <c>-</c> = gap), the identical / similar / gap column counts.</summary>
    public readonly record struct Stats(string Consensus, int Identical, int Similar, int Gaps);

    private static bool Eq(float a, float b) => (b - Epsilon) < a && a < (b + Epsilon);

    private static bool Zero(float a) => Math.Abs((double)a) <= Epsilon;

    /// <summary>Calculates the consensus of the rows (all of equal length).</summary>
    internal static Stats Calculate(IReadOnlyList<string> rows, IReadOnlyList<float> weights, SubstitutionMatrix matrix)
    {
        string alphabet = matrix.Alphabet;
        int matsize = alphabet.Length + 1;
        int Code(char c)
        {
            int i = alphabet.IndexOf(char.ToUpperInvariant(c));
            return i < 0 ? 0 : i + 1;
        }

        int Score(int m1, int m2) => matrix.GetScore(alphabet[m1 - 1], alphabet[m2 - 1]);
        int nseqs = rows.Count;
        int mlen = rows[0].Length;
        float totweight = 0f;
        foreach (float w in weights)
            totweight += w;
        float fplural = totweight * PluralityPercent / 100f;
        float ident = totweight;
        var identical = new float[matsize];
        var matching = new float[matsize];
        var posScore = new float[nseqs];
        var cons = new System.Text.StringBuilder(mlen);
        int retident = 0, retsim = 0, retgap = 0;
        for (int k = 0; k < mlen; k++)
        {
            char res = '-';
            bool isident = false, issim = false, isgap = false;
            Array.Clear(identical);
            Array.Clear(matching);
            Array.Clear(posScore);
            for (int iseq = 0; iseq < nseqs; iseq++)
            {
                int m1 = Code(rows[iseq][k]);
                if (m1 != 0)
                    identical[m1] += weights[iseq];
                for (int jseq = iseq + 1; jseq < nseqs; jseq++)
                {
                    int m2 = Code(rows[jseq][k]);
                    if (m1 != 0 && m2 != 0)
                    {
                        float contri = (Score(m1, m2) * weights[jseq]) + posScore[iseq];
                        float contrj = (Score(m1, m2) * weights[iseq]) + posScore[jseq];
                        posScore[iseq] = contri;
                        posScore[jseq] = contrj;
                    }
                }
            }

            int highindex = -1;
            float max = -float.MaxValue;
            int numres = 0;
            for (int iseq = 0; iseq < nseqs; iseq++)
            {
                if (rows[iseq][k] != ' ' && rows[iseq][k] != '-')
                    numres++;
                if (posScore[iseq] > max)
                {
                    highindex = iseq;
                    max = posScore[iseq];
                }
            }

            for (int iseq = 0; iseq < nseqs; iseq++)
            {
                int m1 = Code(rows[iseq][k]);
                if (!Zero(matching[m1]))
                    continue;
                for (int jseq = 0; jseq < nseqs; jseq++)
                {
                    int m2 = Code(rows[jseq][k]);
                    if (m1 != 0 && m2 != 0 && Score(m1, m2) > 0)
                        matching[m1] += weights[jseq];
                }
            }

            int matchingmaxindex = 0, identicalmaxindex = 0;
            for (int iseq = 0; iseq < nseqs; iseq++)
            {
                int m1 = Code(rows[iseq][k]);
                if (identical[m1] > identical[identicalmaxindex])
                    identicalmaxindex = m1;
            }

            for (int iseq = 0; iseq < nseqs; iseq++)
            {
                int m1 = Code(rows[iseq][k]);
                if (matching[m1] > matching[matchingmaxindex])
                    matchingmaxindex = m1;
                else if (Eq(matching[m1], matching[matchingmaxindex]) && identical[m1] > identical[matchingmaxindex])
                    matchingmaxindex = m1;
                if (rows[iseq][k] is '-' or ' ')
                    isgap = true;
            }

            float himatch = matching[Code(rows[highindex][k])];
            if (identical[identicalmaxindex] >= ident)
                isident = true;
            if (matching[matchingmaxindex] >= fplural)
                issim = true;
            if (himatch >= fplural && rows[highindex][k] != '-')
                res = char.ToUpperInvariant(rows[highindex][k]);
            if (Eq(himatch, weights[highindex]))
                res = numres > 1 ? '#' : '-';
            if (issim && !isident)
                res = char.ToLowerInvariant(res);
            cons.Append(res);
            if (isident)
                retident++;
            if (issim)
                retsim++;
            if (isgap)
                retgap++;
        }

        return new Stats(cons.ToString(), retident, retsim, retgap);
    }

    /// <summary><c>alignSim</c>: upper case -> <paramref name="id"/>, lower case -> <paramref name="sim"/>, <c>#</c> -> <paramref name="mis"/>, anything else -> <paramref name="gap"/>.</summary>
    internal static string Markup(string cons, char id, char sim, char mis, char gap)
    {
        var sb = new System.Text.StringBuilder(cons.Length);
        foreach (char c in cons)
        {
            if (c == '#')
                sb.Append(mis);
            else if (char.IsUpper(c))
                sb.Append(id);
            else if (char.IsLower(c))
                sb.Append(sim);
            else
                sb.Append(gap);
        }

        return sb.ToString();
    }
}
