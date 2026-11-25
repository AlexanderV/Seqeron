// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

using static Seqeron.Genomics.Analysis.MultipleAlignmentStatistics;

/// <summary>Per-sequence statistics of one row of a multiple alignment against a reference row or the consensus
/// (EMBOSS <c>infoalign</c> columns).</summary>
/// <param name="Index">0-based row index.</param>
/// <param name="SequenceLength">Residues of the row: <c>Identities + Similarities + Differences</c> (infoalign <c>SeqLen</c>).</param>
/// <param name="AlignmentLength">Columns from the row's first to its last residue (leading/trailing gap columns excluded; infoalign <c>AlignLen</c>).</param>
/// <param name="GapRuns">Internal gap runs of the row (infoalign <c>Gaps</c>).</param>
/// <param name="GapLength">Internal gap positions of the row (infoalign <c>GapLen</c>).</param>
/// <param name="Identities">Residue = reference residue (case-insensitive; infoalign <c>Ident</c>).</param>
/// <param name="Similarities">Not identical, substitution score &gt; 0 (infoalign <c>Similar</c>).</param>
/// <param name="Differences">Every other residue, including residues facing a reference gap (infoalign <c>Differ</c>).</param>
/// <param name="PercentChange">100·(AlignmentLength − Identities)/AlignmentLength (infoalign <c>% Change</c>; similarities and gap positions count as change).</param>
public sealed record MsaSequenceStatistics(
    int Index,
    int SequenceLength,
    int AlignmentLength,
    int GapRuns,
    int GapLength,
    int Identities,
    int Similarities,
    int Differences,
    double PercentChange);

public static partial class EmbossAlignmentStatistics
{
    /// <summary>
    /// EMBOSS <c>infoalign</c> statistics of every row against a reference row or against the EMBOSS <c>cons</c> consensus.
    /// </summary>
    /// <param name="alignedRows">Aligned rows (at least 2 when the consensus is the reference).</param>
    /// <param name="referenceIndex">0-based reference row, or null for the consensus (infoalign <c>-refseq 0</c>).</param>
    /// <param name="residueType">Matrix choice (EDNAFULL / EBLOSUM62) and consensus type; Auto = type of the first row as EMBOSS.</param>
    /// <remarks>
    /// infoalign_Compare (EMBOSS 6.6.0 <c>emboss/infoalign.c</c>): leading and trailing gaps of the row are ignored; inside, each gap run is
    /// one "Gap" and every '-' a gap position; each residue is Identical (case-insensitive equal to the reference), else Similar when the
    /// substitution score with the reference symbol is &gt; 0, else Different (also a reference gap, or a row longer than the reference);
    /// %Change = 100·(AlignLen − Ident)/AlignLen. A row consisting only of gaps yields zeros (EMBOSS reads past the buffer).
    /// </remarks>
    public static IReadOnlyList<MsaSequenceStatistics> CalculateSequenceStatistics(
        IReadOnlyList<string> alignedRows, int? referenceIndex = null, ConsensusResidueType residueType = ConsensusResidueType.Auto)
    {
        var rows = Normalize(alignedRows);
        bool nucleotide = IsNucleotide(rows, residueType);
        var matrix = nucleotide ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;

        string reference;
        if (referenceIndex is int r)
        {
            if (r < 0 || r >= rows.Length)
                throw new ArgumentOutOfRangeException(nameof(referenceIndex), "Reference sequence number is outside the alignment.");
            reference = rows[r];
        }
        else
        {
            // infoalign: embConsCalc(plurality 50 % of the total weight, setcase 0, identity 0) == cons defaults except setcase.
            reference = EmbossConsensus.GenerateEmbossConsensus(rows, residueType, plurality: null, identity: 0, setcase: 0f);
        }

        var result = new List<MsaSequenceStatistics>(rows.Length);
        for (int i = 0; i < rows.Length; i++)
            result.Add(Compare(i, reference, rows[i], matrix));
        return result;
    }

    private static MsaSequenceStatistics Compare(int index, string reference, string seq, SubstitutionMatrix matrix)
    {
        int begin = 0;
        while (begin < seq.Length && seq[begin] == Gap) begin++;
        int end = seq.Length - 1;
        while (end >= 0 && seq[end] == Gap) end--;
        if (begin > end)
            return new MsaSequenceStatistics(index, 0, 0, 0, 0, 0, 0, 0, 0);

        int gaps = 0, gapLength = 0, ident = 0, sim = 0, diff = 0;
        bool inGap = false;
        for (int i = begin; i <= end; i++)
        {
            char s = seq[i];
            if (s == Gap)
            {
                if (!inGap)
                {
                    inGap = true;
                    gaps++;
                }
                gapLength++;
                continue;
            }
            inGap = false;
            if (i >= reference.Length)
            {
                diff++;
                continue;
            }
            char r = char.ToUpperInvariant(reference[i]);
            if (r == s)
                ident++;
            else if (matrix.IsPositive(r, s))
                sim++;
            else
                diff++;
        }

        int alignLength = end - begin + 1;
        return new MsaSequenceStatistics(index, ident + sim + diff, alignLength, gaps, gapLength, ident, sim, diff,
            (float)(alignLength - ident) * 100f / alignLength);
    }

}
