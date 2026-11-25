// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from Kalign 3.3.5 (https://github.com/TimoLassmann/kalign), Copyright (C) 2006, 2019, 2020 Timo Lassmann; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Kalign, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Kalign;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-20, B09 F118): Kalign 3.3.5 port - the pieces that the whole-program oracle (<see cref="SequenceAligner_KalignOracle_Tests"/>) does not name:
/// the bit-parallel distance against an independent dynamic program (Myers 1999 / the semi-global edit distance), the program's quirks (terminal-gap option, internal type, dropped sequences,
/// alphabet detection), contracts.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_Kalign_Tests
{
    /// <summary>Independent reference: semi-global edit distance of the pattern's first 255 symbols against the best window of the text (free text prefix and suffix, unit costs).</summary>
    private static int SemiGlobalEditDistance(byte[] text, byte[] pattern)
    {
        int m = Math.Min(255, pattern.Length);
        var prev = new int[m + 1];
        var cur = new int[m + 1];
        for (int i = 0; i <= m; i++)
            prev[i] = i;
        int best = prev[m];
        foreach (byte t in text)
        {
            cur[0] = 0;
            for (int i = 1; i <= m; i++)
                cur[i] = Math.Min(Math.Min(prev[i] + 1, cur[i - 1] + 1), prev[i - 1] + (pattern[i - 1] == t ? 0 : 1));
            best = Math.Min(best, cur[m]);
            (prev, cur) = (cur, prev);
        }
        return best;
    }

    [TestCase(1, 40, 60)]
    [TestCase(2, 5, 300)]
    [TestCase(3, 64, 64)]
    [TestCase(4, 65, 200)]
    [TestCase(5, 128, 400)]
    [TestCase(6, 255, 255)]
    [TestCase(7, 256, 600)]
    [TestCase(8, 300, 700)]
    [TestCase(9, 1, 50)]
    public void Bpm_EqualsSemiGlobalEditDistance_Myers1999(int seed, int patternLength, int textLength)
    {
        var rng = new Random(seed);
        for (int rep = 0; rep < 60; rep++)
        {
            int alpha = rep % 3 == 0 ? 2 : rep % 3 == 1 ? 4 : 13;
            var pattern = Enumerable.Range(0, patternLength).Select(_ => (byte)rng.Next(alpha)).ToArray();
            var text = Enumerable.Range(0, textLength).Select(_ => (byte)rng.Next(alpha)).ToArray();
            if (rep % 4 == 0 && textLength >= patternLength)
                Array.Copy(pattern, 0, text, rng.Next(textLength - patternLength + 1), patternLength); // an exact occurrence: distance 0
            Assert.That(KalignAligner.KalignBpm(text, text.Length, pattern, pattern.Length), Is.EqualTo(SemiGlobalEditDistance(text, pattern)), $"seed {seed} rep {rep}");
        }
    }

    [Test]
    public void Bpm_ExactOccurrence_IsZero_AndAbsentPatternIsLimitedByItsLength()
    {
        byte[] text = [0, 1, 2, 3, 0, 1, 2, 3];
        Assert.Multiple(() =>
        {
            Assert.That(KalignAligner.KalignBpm(text, text.Length, [2, 3, 0], 3), Is.Zero);
            Assert.That(KalignAligner.KalignBpm(text, text.Length, [3, 3, 3], 3), Is.EqualTo(2));
        });
    }

    [Test]
    public void TerminalGapOption_IsAppliedOnlyWhenGapExtendIsSet_AQuirkOf335()
    {
        string[] seqs = ["ACGTACGTAGCTAGCTAGCATCGATCGA", "ACGTACGTTAGCTAGCAGCATCGATCGA", "TTACGTACGTAGCTAGCTAGCATC"];
        var def = KalignAligner.MultipleAlignKalign(seqs).AlignedSequences;
        // --tgpe alone: ignored (the program tests gpe >= 0 before applying tgpe)
        Assert.That(KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { TerminalGapExtend = 3000f }).AlignedSequences, Is.EqualTo(def));
        // --gpe alone sets the terminal extension to the unset value -1: equal to giving the default gpe and tgpe -1 explicitly
        Assert.That(KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { GapExtend = 39.4f }).AlignedSequences,
            Is.EqualTo(KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { GapExtend = 39.4f, TerminalGapExtend = -1f }).AlignedSequences));
    }

    [Test]
    public void DnaInternalType_IsDnaWithTerminalExtension8_LibraryOnly()
    {
        string[] seqs = ["ACGTACGTAGCTAGCTAGCATCGATCGA", "ACGTTCGTTAGCTAGCAGCATCGATCGA", "TTACGTACGTAGCTAGCTAGCATC", "ACGTACGTAGCTAGCTAGCATCGATCGATT"];
        var internalType = KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { Type = KalignType.DnaInternal }).AlignedSequences;
        var explicitDna = KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { Type = KalignType.Dna, GapOpen = 8, GapExtend = 6, TerminalGapExtend = 8 }).AlignedSequences;
        Assert.That(internalType, Is.EqualTo(explicitDna));
    }

    [Test]
    public void SequencesWithoutResidues_AreDropped_AndKeptIndicesSaySo()
    {
        var res = KalignAligner.KalignRun(["ACGTACGTAC", "----", "ACGTCGTAC", "...", "ACGTACGAC"]);
        Assert.Multiple(() =>
        {
            Assert.That(res.KeptIndices, Is.EqualTo(new[] { 0, 2, 4 }));
            Assert.That(res.Alignment.AlignedSequences, Has.Length.EqualTo(3));
            Assert.That(res.Alphabet, Is.EqualTo("dna"));
        });
    }

    [Test]
    public void Alphabet_IsDetectedFromTheLetterCounts_NotFromTheTypeOption()
    {
        // letters A C G T only: the DNA letter model wins, even for what is meant as a peptide
        Assert.That(KalignAligner.KalignRun(["ACGTACGTAC", "ACGTCGTAC"]).Alphabet, Is.EqualTo("dna"));
        Assert.That(KalignAligner.KalignRun(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ"]).Alphabet, Is.EqualTo("protein"));
        Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["ACGTACGTAC", "ACGTCGTAC"], new KalignOptions { Type = KalignType.Protein }), "DNA with --type protein");
        Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ"], new KalignOptions { Type = KalignType.Rna }), "protein with --type rna");
    }

    [Test]
    public void Rows_KeepTheInputLettersAndCase_AndCharactersThatAreNotLettersAreGaps()
    {
        string[] seqs = ["acgtACGTac-gt", "ACGTa.cgtAC", "ACGTTCGTAC"];
        var rows = KalignAligner.MultipleAlignKalign(seqs).AlignedSequences;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "acgtACGTacgt", "ACGTacgtAC", "ACGTTCGTAC" }));
            Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
        });
    }

    [Test]
    public void Contracts_RefuseWhatTheProgramRefuses()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["ACGT"]), "only one sequence");
            Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["ACGT", "----"]), "one sequence left");
            Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["ACGT", "ACGTé"]), "non-ASCII");
            Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["MKVLAAGIVGLLLAQJ", "MKVLAGIVALLAQ"]), "J is not in the protein alphabet");
            Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(["ACGT", "ACGA"], new KalignOptions { Threads = 0 }));
            Assert.Throws<ArgumentNullException>(() => KalignAligner.MultipleAlignKalign(null!));
            Assert.Throws<OperationCanceledException>(() => KalignAligner.MultipleAlignKalign(["ACGT", "ACGA"], new KalignOptions { CancellationToken = new CancellationToken(true) }));
        });
    }

    [Test]
    public void OneThreadAndManyThreads_GiveTheProgramsTwoVariants()
    {
        // found with the installed program (see the oracle): -n 1 and -n 4 differ on this set because only the OpenMP variant swaps the shorter sequence first
        string[] seqs = ["ACGTACGTTAGCTAGCTAGCATCGATCGATGCATGCATCG", "ACGTACGTAGCTAGCTAGCATCGATCGATGCATGCCATCG", "ACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGCATCG"];
        var one = KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { Threads = 1 }).AlignedSequences;
        Assert.Multiple(() =>
        {
            Assert.That(one, Is.EqualTo(new[]
            {
                "ACGTACGTTAGCTAGCTAGCATCGATCGATGCATG-CATCG",
                "ACGTACG-TAGCTAGCTAGCATCGATCGATGCATGCCATCG",
                "ACGTTCGTTAGCTAGCTTGCATCGATCGATGCATG-CATCG",
            }));
            Assert.That(KalignAligner.MultipleAlignKalign(seqs, new KalignOptions { Threads = 4, Type = KalignType.Dna }).AlignedSequences[0], Is.EqualTo("ACGTACGTTAGCTAGCTAGCATCGATCGATGCATGC-ATCG"));
        });
    }
}
