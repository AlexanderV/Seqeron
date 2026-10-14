// TRANS-SIXFRAME-001 (audit A3-3) — EMBOSS transeq post-processing qualifiers:
//   TranslateSixFrames(DnaSequence | string, GeneticCode?, SixFrameConvention, TranseqOptions).
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md; TestSpec: tests/TestSpecs/TRANS-SIXFRAME-001.md (M45–M50)
// References:
//   EMBOSS 6.6.0 emboss/transeq.c (main: ajRangeSeqExtract when !wholeseq; ajTrnSeqInit when methionine && !wholeseq;
//   transeq_Trim, then transeq_Clean), emboss/acd/transeq.acd (-regions, -trim, -clean, -methionine default Y),
//   ajax/core/ajrange.c (ajRangeSeqExtract), ajax/core/ajtranslate.c (ajTrnSeqInitC, ajTrnSeqRevinitC "iend > 2",
//   ajTrnSeqAltRevinitC). Locked TSV: TestData/Transeq/transeq_options.tsv = transeq binary output (1320 rows).

using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class Translator_TranseqOptions_Tests
{
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Transeq.transeq_options.tsv";
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

    private static IReadOnlyList<(int Start, int End)>? ParseRegions(string spec) =>
        spec == "-"
            ? null
            : spec.Split(',').Select(p =>
            {
                var a = p.Split('-');
                return (int.Parse(a[0]), int.Parse(a[1]));
            }).ToList();

    private static string[] Frames(IReadOnlyDictionary<int, ProteinSequence> f) =>
        FrameOrder.Select(k => f[k].Sequence).ToArray();

    // M45 — EMBOSS 6.6.0 transeq binary lock: 5 sequences × tables 1/2/11 × -alternative N/Y × -trim × -clean ×
    // -regions (6 specs) × -methionine Y/N, all six frames.
    [Test]
    public void TranslateSixFrames_TranseqOptions_MatchTranseqBinary()
    {
        var rows = LoadRows().ToList();
        Assert.That(rows, Has.Count.EqualTo(1320));

        Assert.Multiple(() =>
        {
            foreach (var r in rows)
            {
                var options = new TranseqOptions
                {
                    Trim = r[2] == "Y",
                    Clean = r[3] == "Y",
                    Methionine = r[4] == "Y",
                    Regions = ParseRegions(r[5]),
                };
                var convention = r[1] == "Y"
                    ? SixFrameConvention.EmbossTranseqAlternative
                    : SixFrameConvention.EmbossTranseq;
                var frames = Translator.TranslateSixFrames(r[6], GeneticCode.GetByTableNumber(int.Parse(r[0])),
                    convention, options);
                for (int i = 0; i < 6; i++)
                    Assert.That(frames[FrameOrder[i]].Sequence, Is.EqualTo(r[7 + i]),
                        $"table {r[0]} alt {r[1]} trim {r[2]} clean {r[3]} met {r[4]} regions {r[5]} {r[6]} frame {FrameOrder[i]}");
            }
        });
    }

    // M46 — hand-checked TTGAAATAGCCCTGA, table 11, frame 1 (transeq binary):
    // whole sequence LK*P*; -regions 1-6,13-15,7-9 → TTGAAA TGA TAG → MK** (initial TTG → M);
    // + -trim -clean → MK (trailing * removed before cleaning).
    [Test]
    public void TranslateSixFrames_Regions_ConcatenatedInOrder_MethionineTrimClean()
    {
        const string seq = "TTGAAATAGCCCTGA";
        var code = GeneticCode.GetByTableNumber(11);
        var regions = new[] { (1, 6), (13, 15), (7, 9) };
        var whole = Translator.TranslateSixFrames(seq, code, SixFrameConvention.EmbossTranseq, TranseqOptions.Default);
        var reg = Translator.TranslateSixFrames(seq, code, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Regions = regions });
        var regNoMet = Translator.TranslateSixFrames(seq, code, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Regions = regions, Methionine = false });
        var trimmed = Translator.TranslateSixFrames(seq, code, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Regions = regions, Trim = true, Clean = true });
        var cleanOnly = Translator.TranslateSixFrames(seq, code, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Clean = true });

        Assert.Multiple(() =>
        {
            Assert.That(whole[1].Sequence, Is.EqualTo("LK*P*"), "-methionine acts only with -regions");
            Assert.That(reg[1].Sequence, Is.EqualTo("MK**"));
            Assert.That(regNoMet[1].Sequence, Is.EqualTo("LK**"));
            Assert.That(trimmed[1].Sequence, Is.EqualTo("MK"));
            Assert.That(cleanOnly[1].Sequence, Is.EqualTo("LKXPX"));
        });
    }

    // M47 — ajTrnSeqRevinitC (phase-locked reverse frames) replaces the first codon only when the frame has at
    // least two whole codons (iend > 2); ajTrnSeqAltRevinitC needs one. transeq -table 11 -regions 1-len:
    // CAA → -1: L (phase-locked) / M (-alternative); GGGCAA → -1: MP under both.
    [Test]
    public void TranslateSixFrames_Methionine_PhaseLockedReverseNeedsTwoCodons()
    {
        var code = GeneticCode.GetByTableNumber(11);
        TranseqOptions Opt(string s) => new() { Regions = new[] { (1, s.Length) } };

        Assert.Multiple(() =>
        {
            Assert.That(Frames(Translator.TranslateSixFrames("CAA", code, SixFrameConvention.EmbossTranseq, Opt("CAA"))),
                Is.EqualTo(new[] { "Q", "X", "X", "L", "X", "X" }));
            Assert.That(Frames(Translator.TranslateSixFrames("CAA", code, SixFrameConvention.EmbossTranseqAlternative, Opt("CAA"))),
                Is.EqualTo(new[] { "Q", "X", "X", "M", "X", "X" }));
            Assert.That(Frames(Translator.TranslateSixFrames("GGGCAA", code, SixFrameConvention.EmbossTranseq, Opt("GGGCAA"))),
                Is.EqualTo(new[] { "GQ", "GX", "AX", "MP", "AX", "CP" }));
            Assert.That(Frames(Translator.TranslateSixFrames(new DnaSequence("GGGCAA"), code,
                    SixFrameConvention.EmbossTranseqAlternative, Opt("GGGCAA"))),
                Is.EqualTo(new[] { "GQ", "GX", "AX", "MP", "CP", "AX" }));
        });
    }

    // M48 — region end past the sequence end is clamped (ajStrAppendSubS); transeq -regions 5-40 on
    // ATGAAACCCGGGTTTTAG (18 nt) = -regions 5-18: NPGFX, TRVLX, PGF*, KTRV, *NPGX, LKPGX.
    [Test]
    public void TranslateSixFrames_RegionEndPastSequenceEnd_IsClamped()
    {
        const string seq = "ATGAAACCCGGGTTTTAG";
        var clamped = Translator.TranslateSixFrames(seq, null, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Regions = new[] { (5, 40) } });
        Assert.That(Frames(clamped), Is.EqualTo(new[] { "NPGFX", "TRVLX", "PGF*", "KTRV", "*NPGX", "LKPGX" }));
    }

    // M49 — defaults reproduce the convention overload; -trim on an all-X/* frame empties it.
    [TestCase("ATGAAACCCGGGTT", 1)]
    [TestCase("acgunRYatgTAG", 2)]
    [TestCase("TTGATTGTGTAA", 11)]
    public void TranslateSixFrames_DefaultOptions_EqualConventionOverload(string seq, int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        foreach (var conv in new[] { SixFrameConvention.EmbossTranseq, SixFrameConvention.EmbossTranseqAlternative })
        {
            Assert.That(Frames(Translator.TranslateSixFrames(seq, code, conv, TranseqOptions.Default)),
                Is.EqualTo(Frames(Translator.TranslateSixFrames(seq, code, conv))), conv.ToString());
            Assert.That(Frames(Translator.TranslateSixFrames(seq, code, conv, new TranseqOptions { Methionine = false })),
                Is.EqualTo(Frames(Translator.TranslateSixFrames(seq, code, conv))), "no regions: methionine inert");
        }

        var trimmedA = Translator.TranslateSixFrames("A", null, SixFrameConvention.EmbossTranseq,
            new TranseqOptions { Trim = true });
        Assert.That(trimmedA.Values.All(p => p.Sequence.Length == 0), "'A' → X,'','','',X,X; all trimmed");
    }

    // M50 — argument checks: null options, Biopython convention, invalid regions.
    [Test]
    public void TranslateSixFrames_TranseqOptions_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, SixFrameConvention.EmbossTranseq, null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => Translator.TranslateSixFrames((DnaSequence)null!, null, SixFrameConvention.EmbossTranseq,
                TranseqOptions.Default), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, SixFrameConvention.Biopython,
                TranseqOptions.Default), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, SixFrameConvention.EmbossTranseq,
                new TranseqOptions { Regions = new[] { (0, 3) } }), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, SixFrameConvention.EmbossTranseq,
                new TranseqOptions { Regions = new[] { (3, 2) } }), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, SixFrameConvention.EmbossTranseq,
                new TranseqOptions { Regions = new[] { (5, 6) } }), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
                "start past the end (transeq would silently read the last base)");
            Assert.That(() => Translator.TranslateSixFrames("ACGT", null, (SixFrameConvention)7,
                TranseqOptions.Default), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }
}
