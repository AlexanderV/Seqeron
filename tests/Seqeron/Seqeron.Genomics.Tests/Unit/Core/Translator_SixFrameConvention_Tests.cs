// TRANS-SIXFRAME-001 (audit A2-10) — EMBOSS transeq six-frame conventions:
//   TranslateSixFrames(DnaSequence | string, GeneticCode?, SixFrameConvention).
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md; TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md (M39–M44)
// References:
//   EMBOSS 6.6.0 transeq.c (frame -1: ajTrnSeqOrig(seq, 2-i) "uses frame 1 codons"; -alternative: ajTrnSeqOrig(seq, -i-1)
//   "uses codons starting at end"), acd/transeq.acd (-alternative help), ajax/core/ajtranslate.c (ajTrnSeqFrameC,
//   ajTrnSeqRevC, ajTrnSeqAltRevC, ajTrnSeqDangleS). Locked TSV: TestData/Transeq/transeq_six_frames.tsv =
//   `transeq -frame 6 -table T -alternative N|Y` binary output, 80 sequences × tables 1/2/11 (480 rows).

using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class Translator_SixFrameConvention_Tests
{
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Transeq.transeq_six_frames.tsv";
    private static readonly int[] FrameOrder = { 1, 2, 3, -1, -2, -3 };

    private static IEnumerable<string[]> LoadRows()
    {
        var asm = Assembly.GetExecutingAssembly();
        using var stream = asm.GetManifestResourceStream(Resource)
            ?? throw new InvalidOperationException($"Embedded oracle '{Resource}' not found.");
        using var reader = new StreamReader(stream);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            yield return line.Split('\t');
        }
    }

    // M39 — EMBOSS 6.6.0 transeq binary lock: 80 sequences × tables 1/2/11 × -alternative N/Y, all six frames.
    [Test]
    public void TranslateSixFrames_EmbossConventions_MatchTranseqBinary()
    {
        var rows = LoadRows().ToList();
        Assert.That(rows, Has.Count.EqualTo(480), "80 sequences × 3 tables × 2 conventions");

        Assert.Multiple(() =>
        {
            foreach (var r in rows)
            {
                var code = GeneticCode.GetByTableNumber(int.Parse(r[0]));
                var convention = r[1] == "Y"
                    ? SixFrameConvention.EmbossTranseqAlternative
                    : SixFrameConvention.EmbossTranseq;
                var frames = Translator.TranslateSixFrames(r[2], code, convention);
                for (int i = 0; i < 6; i++)
                    Assert.That(frames[FrameOrder[i]].Sequence, Is.EqualTo(r[3 + i]),
                        $"table {r[0]}, -alternative {r[1]}, {r[2]}, frame {FrameOrder[i]}");
            }
        });
    }

    // M40 — hand-traced ATGAAACCCGGGTT (14 nt), transeq -frame 6 (Evidence):
    // forward frames translate the dangling bases (TT → TTN = F/L → X; T → X);
    // default -1 = revcomp(frame-1 codons ATG AAA CCC GGG) = PGFH, -2 = revcomp(TGA AAC CCG GGT) + X (leading A),
    // -3 = revcomp(GAA ACC CGG GTT) + comp(T)comp(A)N = ATN → X; -alternative: -1 = NPGF + X, -2 = TRVS + X, -3 = PGFH.
    [Test]
    public void TranslateSixFrames_EmbossConventions_HandTracedExample()
    {
        const string seq = "ATGAAACCCGGGTT";
        var locked = Translator.TranslateSixFrames(seq, null, SixFrameConvention.EmbossTranseq);
        var alt = Translator.TranslateSixFrames(new DnaSequence(seq), null, SixFrameConvention.EmbossTranseqAlternative);

        Assert.Multiple(() =>
        {
            Assert.That(locked[1].Sequence, Is.EqualTo("MKPGX"));
            Assert.That(locked[2].Sequence, Is.EqualTo("*NPGX"));
            Assert.That(locked[3].Sequence, Is.EqualTo("ETRV"));
            Assert.That(locked[-1].Sequence, Is.EqualTo("PGFH"), "phase-locked: frame-1 codons reverse-complemented");
            Assert.That(locked[-2].Sequence, Is.EqualTo("TRVSX"));
            Assert.That(locked[-3].Sequence, Is.EqualTo("NPGFX"));
            Assert.That(alt[1].Sequence, Is.EqualTo("MKPGX"), "forward frames are convention-independent");
            Assert.That(alt[-1].Sequence, Is.EqualTo("NPGFX"), "alternative: starts at the last base");
            Assert.That(alt[-2].Sequence, Is.EqualTo("TRVSX"));
            Assert.That(alt[-3].Sequence, Is.EqualTo("PGFH"));
        });
    }

    // M41 — two-base dangling end translated as codon b1 b2 N when the third base does not matter
    // (GGN → G, CTN → L under table 1), transeq ajTrnSeqDangleS.
    [Test]
    public void TranslateSixFrames_Emboss_DecidableDanglingPair_IsTranslated()
    {
        var frames = Translator.TranslateSixFrames("ATGGG", null, SixFrameConvention.EmbossTranseq);
        var biopython = Translator.TranslateSixFrames("ATGGG", null, SixFrameConvention.Biopython);

        Assert.Multiple(() =>
        {
            Assert.That(frames[1].Sequence, Is.EqualTo("MG"), "ATG + GGN → M G");
            Assert.That(biopython[1].Sequence, Is.EqualTo("M"), "Biopython drops the partial codon");
            Assert.That(frames[-1].Sequence, Is.EqualTo("H"), "revcomp(ATG) = CAT → H; trailing GG not translated");
            Assert.That(frames[-3].Sequence, Is.EqualTo("PX"), "revcomp(GGG) = CCC → P; dangle comp(T)comp(A)N = ATN → X");
        });
    }

    // M42 — Biopython convention through the new overloads is identical to the existing API.
    [TestCase("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 1)]
    [TestCase("ATGAAACCCGGGTT", 2)]
    [TestCase("acgunRYatg", 11)]
    public void TranslateSixFrames_BiopythonConvention_EqualsDefault(string seq, int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        var expected = Translator.TranslateSixFrames(seq, code);
        var actual = Translator.TranslateSixFrames(seq, code, SixFrameConvention.Biopython);

        Assert.Multiple(() =>
        {
            foreach (int k in FrameOrder)
                Assert.That(actual[k].Sequence, Is.EqualTo(expected[k].Sequence), $"frame {k}");
        });
    }

    // M43 — -alternative equals Biopython on whole codons: when (len − k + 1) is a multiple of 3 for every frame
    // is impossible, but for len ≡ 0 (mod 3) frame ±1 have no dangling bases, so they coincide.
    [Test]
    public void TranslateSixFrames_Alternative_EqualsBiopythonWithoutDanglingBases()
    {
        const string seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"; // 39 nt
        var alt = Translator.TranslateSixFrames(seq, null, SixFrameConvention.EmbossTranseqAlternative);
        var bio = Translator.TranslateSixFrames(seq, null);

        Assert.Multiple(() =>
        {
            Assert.That(alt[1].Sequence, Is.EqualTo(bio[1].Sequence));
            Assert.That(alt[-1].Sequence, Is.EqualTo(bio[-1].Sequence));
            Assert.That(alt[-2].Sequence, Is.EqualTo(bio[-2].Sequence + "X"), "38 bases: 12 codons + 2 dangling (→ X or aa)");
        });
    }

    // M44 — documented ambiguity difference: the EMBOSS conventions keep this library's IUPAC codon
    // resolution (Biopython tables: RAY → B, SAR → Z), where transeq prints X.
    [Test]
    public void TranslateSixFrames_Emboss_KeepsIupacResolution()
    {
        var frames = Translator.TranslateSixFrames("RAYSARGCN", null, SixFrameConvention.EmbossTranseq);
        Assert.That(frames[1].Sequence, Is.EqualTo("BZA"), "transeq 6.6.0 prints XXA (documented difference)");
    }

    // S9 — edge cases and argument checks.
    [Test]
    public void TranslateSixFrames_Emboss_ShortInputs_AndArguments()
    {
        var empty = Translator.TranslateSixFrames("", null, SixFrameConvention.EmbossTranseq);
        var one = Translator.TranslateSixFrames("A", null, SixFrameConvention.EmbossTranseq);
        var oneAlt = Translator.TranslateSixFrames("A", null, SixFrameConvention.EmbossTranseqAlternative);

        Assert.Multiple(() =>
        {
            Assert.That(empty.Values.All(p => p.Sequence.Length == 0), "empty input → six empty frames");
            Assert.That(FrameOrder.Select(k => one[k].Sequence), Is.EqualTo(new[] { "X", "", "", "", "X", "X" }),
                "transeq on 'A' (phase-locked)");
            Assert.That(FrameOrder.Select(k => oneAlt[k].Sequence), Is.EqualTo(new[] { "X", "", "", "X", "", "" }),
                "transeq -alternative on 'A'");
            Assert.That(() => Translator.TranslateSixFrames((DnaSequence)null!, null, SixFrameConvention.EmbossTranseq),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, (SixFrameConvention)7),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }
}
