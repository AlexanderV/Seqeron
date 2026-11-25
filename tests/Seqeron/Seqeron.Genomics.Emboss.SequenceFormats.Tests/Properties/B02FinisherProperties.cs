// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
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
    /// F53: on A/C/G/T, every frame of the transeq -alternative convention is the Biopython frame plus the
    /// dangling-base residue (none for 0 dangling bases, X for one, the b1 b2 N codon for two); forward frames
    /// are identical in both EMBOSS conventions.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property TranseqAlternative_IsBiopythonPlusDangle()
    {
        var gen = from code in CodeGen() from s in StringGen(Acgt, 1, 50) select (code, s);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var bio = EmbossTranslator.TranslateSixFrames(t.s, t.code, SixFrameConvention.Biopython);
            var alt = EmbossTranslator.TranslateSixFrames(t.s, t.code, SixFrameConvention.EmbossTranseqAlternative);
            var locked = EmbossTranslator.TranslateSixFrames(t.s, t.code, SixFrameConvention.EmbossTranseq);
            string rc = DnaSequence.GetReverseComplementString(t.s);
            return FrameKeys.All(k =>
            {
                string strand = k > 0 ? t.s : rc;
                int offset = Math.Abs(k) - 1;
                int dangle = Math.Max(0, strand.Length - offset) % 3;
                string extra = dangle switch
                {
                    0 => "",
                    1 => "X",
                    _ => t.code.Translate(string.Concat(strand.AsSpan(strand.Length - 2), "N")).ToString(),
                };
                return alt[k].Sequence == bio[k].Sequence + extra
                    && (k < 0 || locked[k].Sequence == alt[k].Sequence);
            });
        });
    }

    /// <summary>
    /// F58: TranseqOptions.Default equals the plain convention overload; -trim removes exactly the trailing
    /// X/* run; -clean (after trim) turns every '*' into 'X'; both are idempotent on their own output.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property TranseqTrimClean_ExactAndIdempotent()
    {
        var gen =
            from code in CodeGen()
            from s in StringGen("ACGTACGTACGTN", 1, 50)
            from convention in Gen.Elements(SixFrameConvention.EmbossTranseq, SixFrameConvention.EmbossTranseqAlternative)
            select (code, s, convention);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var plain = EmbossTranslator.TranslateSixFrames(t.s, t.code, t.convention);
            var dflt = EmbossTranslator.TranslateSixFrames(t.s, t.code, t.convention, TranseqOptions.Default);
            var trim = EmbossTranslator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Trim = true });
            var clean = EmbossTranslator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Clean = true });
            var both = EmbossTranslator.TranslateSixFrames(t.s, t.code, t.convention, new TranseqOptions { Trim = true, Clean = true });
            return FrameKeys.All(k =>
            {
                string p = plain[k].Sequence;
                string trimmed = p.TrimEnd('X', '*');
                return dflt[k].Sequence == p
                    && trim[k].Sequence == trimmed
                    && clean[k].Sequence == p.Replace('*', 'X')
                    && both[k].Sequence == trimmed.Replace('*', 'X')
                    && trimmed.TrimEnd('X', '*') == trimmed
                    && both[k].Sequence.Replace('*', 'X') == both[k].Sequence;
            });
        });
    }

    private static readonly int[] FrameKeys = { 1, 2, 3, -1, -2, -3 };

    private const string Acgt = "ACGT";

 // ſ ı K é Å а

    private static Gen<string> StringGen(string alphabet, int minLen, int maxLen) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alphabet.ToCharArray()).ArrayOf(n)
        select new string(chars);

    private static Gen<GeneticCode> CodeGen() =>
        Gen.Elements(GeneticCode.SupportedTableNumbers.ToArray()).Select(GeneticCode.GetByTableNumber);
}
