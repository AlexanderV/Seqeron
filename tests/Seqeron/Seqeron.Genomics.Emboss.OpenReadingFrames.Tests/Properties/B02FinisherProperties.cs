// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.CodonOptimizer;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviours added by the B02 finisher
/// (docs/Validation/review-2026-09/B02.md, F31–F61 and the WP24/WP25 sweep): ASCII-only case folding
/// (F31, SequenceCase), Sec/Pyl acceptance and additive molecular weight (F32), Biopython nucleotide X
/// (F33), explicit not-calculable ENC and the ENC estimator variants (F36/F50/F51/F57), CodonW GC
/// convention (F37), RSCU pseudocount / cubar sub-families (F52/F59), Fop/CBI (F61), getorf modes
/// (F44–F46), transeq conventions and post-processing (F53/F58), getorf header ranges (F56) and
/// <see cref="GeneticCode.GetByTableNumberOrThrow"/>.
///
/// Every property asserts an exact algebraic relation, a documented range, or the equality of two
/// public paths; floating-point comparisons use 1e-9 only where a sum is re-associated.
///
/// Test Units: TRANS-PROT-001, TRANS-CODON-001, TRANS-SIXFRAME-001, CODON-ENC-001, CODON-RSCU-001,
/// CODON-STATS-001, B02-SWEEP.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Translation")]
public class B02FinisherProperties
{

    /// <summary>
    /// F44: with -minsize 0 the getorf -find 0 records of each frame tile that frame's complete codons
    /// stop-to-stop: contiguous, from the frame offset to the last complete codon; each protein has one residue
    /// per non-stop codon, no '*', and a terminated record ends on a codon the code reads as '*'.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode0_PartitionsFramesStopToStop()
    {
        return Prop.ForAll(OrfInput("ACGTACGTACGTN", 90).ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            int len = upper.Length;
            var records = EmbossGetorf.FindOrfs(t.S, OrfFindMode.TranslationStopToStop, t.Code, minSize: 0).ToList();
            foreach (int frame in new[] { 1, 2, 3, -1, -2, -3 })
            {
                string strand = Strand(upper, frame);
                int offset = Math.Abs(frame) - 1;
                var inFrame = records.Where(r => r.Frame == frame).ToList();
                int codons = Math.Max(0, (len - offset) / 3);
                if (codons == 0)
                {
                    if (inFrame.Count != 0) return false;
                    continue;
                }
                int expectedStart = offset;
                foreach (var r in inFrame)
                {
                    if (r.StartPosition != expectedStart || (r.EndPosition - r.StartPosition + 1) % 3 != 0) return false;
                    int n = (r.EndPosition - r.StartPosition + 1) / 3;
                    if (r.Sequence.Length != n - (r.IsTerminated ? 1 : 0) || r.Sequence.Contains('*')) return false;
                    string last = strand.Substring(r.EndPosition - 2, 3);
                    if (r.IsTerminated != (t.Code.Translate(last) == '*')) return false;
                    expectedStart = r.EndPosition + 1;
                }
                // Covered up to the last complete codon (a frame ending on a STOP leaves no trailing record).
                if (expectedStart != offset + 3 * codons) return false;
            }
            return true;
        });
    }

    /// <summary>F44: mode 2 nucleotides are the scanned strand's bases of the record without the STOP; modes 0 and 2 give the same regions.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode2_SequenceIsStrandSubstring()
    {
        return Prop.ForAll(OrfInput("ACGTACGTN", 90).ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            var nuc = EmbossGetorf.FindOrfs(t.S, OrfFindMode.NucleotidesStopToStop, t.Code, minSize: 0).ToList();
            var prot = EmbossGetorf.FindOrfs(t.S, OrfFindMode.TranslationStopToStop, t.Code, minSize: 0).ToList();
            bool sameRegions = nuc.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated))
                .SequenceEqual(prot.Select(r => (r.StartPosition, r.EndPosition, r.Frame, r.IsTerminated)));
            return sameRegions && nuc.All(r =>
            {
                int length = r.EndPosition - r.StartPosition + 1 - (r.IsTerminated ? 3 : 0);
                return r.Sequence == Strand(upper, r.Frame).Substring(r.StartPosition, length);
            });
        });
    }

    /// <summary>
    /// F44: flank windows (modes 4–6) lie inside the sequence, are 3 + 2·flanking long, are centred on the
    /// codon, hold the strand's bases, and the flanked codon is a START (mode 4) or a STOP (modes 5, 6).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfFlankWindows_WithinBounds()
    {
        var gen =
            from input in OrfInput("ACGTACGTACGTN", 120)
            from mode in Gen.Elements(OrfFindMode.FlankingStart, OrfFindMode.FlankingInitialStop, OrfFindMode.FlankingEndingStop)
            from flanking in Gen.Choose(0, 12)
            select (input.Code, input.S, mode, flanking);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string upper = Dna(t.S);
            return EmbossGetorf.FindOrfs(t.S, t.mode, t.Code, minSize: 0, flanking: t.flanking).All(r =>
            {
                string strand = Strand(upper, r.Frame);
                string codon = strand.Substring(r.CodonPosition, 3);
                bool codonKind = t.mode == OrfFindMode.FlankingStart
                    ? t.Code.IsStartCodon(codon)
                    : t.Code.Translate(codon) == '*';
                return r.StartPosition >= 0 && r.EndPosition < upper.Length
                    && r.EndPosition - r.StartPosition + 1 == 3 + 2 * t.flanking
                    && r.CodonPosition == r.StartPosition + t.flanking
                    && r.Sequence == strand.Substring(r.StartPosition, 3 + 2 * t.flanking)
                    && codonKind;
            });
        });
    }

    /// <summary>
    /// F44/F56: getorf -find 1 equals the legacy START→STOP overload (minSize = 3·minLength), and the two
    /// ToGetOrfRange helpers agree record by record.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property GetOrfMode1_EqualsLegacy_AndRangesAgree()
    {
        var gen = from input in OrfInput("ACGTACGTACGTATGATG", 120) from minLength in Gen.Choose(0, 6) select (input.Code, input.S, minLength);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            int len = t.S.Length;
            var legacy = Translator.FindOrfs(t.S, t.Code, t.minLength)
                .Select(o => (o.StartPosition, o.EndPosition, o.Frame, Seq: o.Protein.Sequence,
                    Range: len > 0 ? o.ToGetOrfRange(len) : default))
                .Order().ToList();
            var getorf = EmbossGetorf.FindOrfs(t.S, OrfFindMode.TranslationStartToStop, t.Code, minSize: 3 * t.minLength)
                .Select(o => (o.StartPosition, o.EndPosition, o.Frame, Seq: o.Sequence,
                    Range: len > 0 ? o.ToGetOrfRange(len) : default))
                .Order().ToList();
            return legacy.SequenceEqual(getorf);
        });
    }

    /// <summary>
    /// F46: -circular keeps every linear terminated START→STOP ORF, and every linear terminated STOP→STOP ORF
    /// that does not start the frame (the leading region is not an ORF on a circle).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Circular_PreservesLinearOrfsNotCrossingOrigin()
    {
        var gen = from input in OrfInput("ACGTACGTATG", 90)
                  from mode in Gen.Elements(OrfFindMode.TranslationStartToStop, OrfFindMode.TranslationStopToStop,
                      OrfFindMode.NucleotidesStartToStop, OrfFindMode.NucleotidesStopToStop)
                  select (input.Code, input.S, mode);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var circular = EmbossGetorf.FindOrfs(t.S, t.mode, t.Code, minSize: 0, circular: true).ToHashSet();
            bool stopToStop = t.mode is OrfFindMode.TranslationStopToStop or OrfFindMode.NucleotidesStopToStop;
            return EmbossGetorf.FindOrfs(t.S, t.mode, t.Code, minSize: 0)
                .Where(r => r.IsTerminated && (!stopToStop || r.StartPosition >= 3))
                .All(circular.Contains);
        });
    }



    private static string Strand(string upper, int frame) =>
        frame > 0 ? upper : DnaSequence.GetReverseComplementString(upper);

    private static string Dna(string s) => SequenceCase.ToUpperDna(s);

    private static Gen<(GeneticCode Code, string S)> OrfInput(string alphabet, int maxLen) =>
        from code in OrfCodeGen() from s in StringGen(alphabet, 0, maxLen) select (code, s);

 // ſ ı K é Å а

    private static Gen<string> StringGen(string alphabet, int minLen, int maxLen) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alphabet.ToCharArray()).ArrayOf(n)
        select new string(chars);

    // Tables whose stops are not dual-coding (the getorf / ORF API rejects 27, 28, 31).
    private static Gen<GeneticCode> OrfCodeGen() =>
        Gen.Elements(GeneticCode.SupportedTableNumbers.Where(t => t is not (27 or 28 or 31)).ToArray())
            .Select(GeneticCode.GetByTableNumber);
}
