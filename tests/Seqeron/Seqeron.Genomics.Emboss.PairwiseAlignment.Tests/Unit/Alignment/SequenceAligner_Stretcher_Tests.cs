// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 (audit A2-28, B09 F121): <see cref="EmbossPairwisePrograms.StretcherAlign"/> = EMBOSS 6.6.0 <c>stretcher</c> (Myers-Miller linear-space global alignment with the program's own midpoint tie rules and edit-script
/// conversion) against the installed binary: TestData/Alignment/stretcher_oracle.tsv (scripts/alignment/generate_stretcher_oracle.py), 689 random DNA / IUPAC / protein pairs (1 - 260 residues, integer gap penalties 0 - 20 / 0 - 6, 7 protein
/// and 3 nucleotide matrices): both rows and the score are identical on every pair.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SequenceAligner_Stretcher_Tests
{
    private sealed record Row(string Id, string Matrix, int GapOpen, int GapExtend, string Seq1, string Seq2, string Row1, string Row2, int Score);


    [Test]
    public void StretcherAlign_RowsAndScore_AreIdenticalToEmbossOnAllOracleCases()
    {
        var bad = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = EmbossPairwisePrograms.StretcherAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend);
            if (a.AlignedSequence1 != r.Row1 || a.AlignedSequence2 != r.Row2 || a.Score != r.Score)
                bad.Add($"{r.Id} {r.Matrix} {r.GapOpen}/{r.GapExtend}");
        }

        Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(100));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    [Test]
    public void StretcherAlign_Score_IsTheEmbossAffineGlobalOptimum_OnEveryOneOfTheFirst300Rows()
    {
        // the same optimum as the quadratic-space Gotoh program under the EMBOSS gap convention (gapopen + (k - 1) * gapextend): stretcher's Myers-Miller program has no forbidden gap-state transition
        // (measured, A3-14b: 300 of 300 sample rows are equal - no tolerance)
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => x.GapOpen >= x.GapExtend).Take(300))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var s = EmbossPairwisePrograms.StretcherAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend);
            var g = SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend));
            Assert.That(s.Score, Is.EqualTo(g.Score), r.Id);
        }

        Assert.That(total, Is.EqualTo(300));
    }

    [Test]
    public void StretcherAlign_Defaults_FollowTheSequenceType()
    {
        // EMBOSS: nucleotide EDNAFULL 16 / 4, protein EBLOSUM62 12 / 2
        var dna = EmbossPairwisePrograms.StretcherAlign("ACGTACGTAGCTAGCTAGCATCGATCGATCGAT", "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT");
        Assert.That(dna.Score, Is.EqualTo(140));
        Assert.That(dna.AlignedSequence1, Is.EqualTo("ACGTACGTAGCTAGCTAGCATCGA-TCGATCGAT"));
        var prot = EmbossPairwisePrograms.StretcherAlign("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV", "MKTAYIAKQRQISFVKSHFSRQLERLGLIEVQAPILSRVGDGTQDNLSGAEKAVKVKVKALPDAQFEV");
        Assert.That(prot.Score, Is.EqualTo(303));
    }

    [Test]
    public void StretcherAlign_Contract_NullEmptyAndPenalties()
    {
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.StretcherAlign(null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.StretcherAlign("ACGT", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.StretcherAlign("ACGT", "ACGT", null, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.StretcherAlign("ACGT", "ACGT", null, 10, -2));
        var empty = EmbossPairwisePrograms.StretcherAlign("", "ACG");
        Assert.That(empty.AlignedSequence1, Is.EqualTo("---"));
        Assert.That(empty.AlignedSequence2, Is.EqualTo("ACG"));
        Assert.That(empty.Score, Is.EqualTo(-(16 + (2 * 4))));
        Assert.That(EmbossPairwisePrograms.StretcherAlign("", "").AlignedSequence1, Is.Empty);
    }

    [Test]
    public void StretcherAlign_LowerCaseInput_IsUpperCased()
    {
        var a = EmbossPairwisePrograms.StretcherAlign("acgtacgt", "ACGTACGT");
        Assert.That(a.AlignedSequence1, Is.EqualTo("ACGTACGT"));
        Assert.That(a.Score, Is.EqualTo(40));
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.stretcher_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], int.Parse(f[2], CultureInfo.InvariantCulture), int.Parse(f[3], CultureInfo.InvariantCulture), f[4], f[5], f[6].ToUpperInvariant(), f[7].ToUpperInvariant(), int.Parse(f[8], CultureInfo.InvariantCulture)));
        }

        return list;
    });
}
