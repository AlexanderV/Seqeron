using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.IO.QualityScoreAnalyzer;
using QualityEncoding = Seqeron.Genomics.IO.QualityScoreAnalyzer.QualityEncoding;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// PARSE-FASTQ-001 — file-level (multi-read) Phred-offset detection,
/// <see cref="Seqeron.Genomics.IO.QualityScoreAnalyzer.DetectEncoding(System.Collections.Generic.IEnumerable{string})"/>.
///
/// The single-string detector cannot resolve a read whose characters all fall in the Phred+33/Phred+64
/// overlap range (ASCII 64-74); the file-level detector resolves it from the rest of the reads, exactly as
/// FastQC / Biopython do. Expected encodings are derived from the canonical ASCII ranges in
/// Cock et al. (2010), Nucleic Acids Research 38(6):1767-1771 (https://doi.org/10.1093/nar/gkp1137):
/// Sanger/Phred+33 = ASCII 33-126; Illumina 1.3+/Phred+64 = ASCII 64-126. A character &lt; 64 is therefore
/// outside Phred+64 (proven Phred+33); a character &gt; 74 (Phred 41, the Illumina 1.8+ ceiling) with no
/// low character indicates Phred+64; everything else is genuinely ambiguous.
/// </summary>
[TestFixture]
public class QualityScoreAnalyzer_DetectEncodingFile_Tests
{
    // Build a quality string from raw ASCII codes (so the byte ranges under test are explicit).
    private static string Ascii(params int[] codes) => new string(codes.Select(c => (char)c).ToArray());

    // ---- the case the limitation named: a lone overlap-only read is ambiguous; the file resolves it ----

    [Test]
    public void SingleOverlapOnlyRead_IsAmbiguous_DefaultsPhred33()
    {
        // Every character in [64, 73] ('@'..'I'): valid under BOTH encodings -> undeterminable.
        var read = Ascii(64, 70, 73, 66);

        var r = DetectEncoding(new[] { read });

        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Ambiguous));
        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33)); // documented safe default
        Assert.That(r.MinAscii, Is.EqualTo(64));
        Assert.That(r.MaxAscii, Is.EqualTo(73));
    }

    [Test]
    public void OverlapOnlyRead_ResolvedByAnotherReadWithLowChar_Phred33Definitive()
    {
        // read1 is overlap-only (ambiguous alone); read2 carries '#' (ASCII 35), impossible in Phred+64.
        var read1 = Ascii(64, 70, 73);           // overlap range only
        var read2 = Ascii(35, 40, 64, 73);       // '#' = ASCII 35 < 64 -> proves Phred+33

        var r = DetectEncoding(new[] { read1, read2 });

        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Definitive));
        Assert.That(r.MinAscii, Is.EqualTo(35));
        // The lone read alone could not have decided this: the per-read detector merely defaults.
        Assert.That(DetectEncoding(read1), Is.EqualTo(QualityEncoding.Phred33)); // default, not proof
    }

    [Test]
    public void OverlapOnlyRead_ResolvedByAnotherReadWithHighChar_Phred64Inferred()
    {
        // read1 overlap-only; read2 has 'h' (ASCII 104) and no char below 64 -> inferred Phred+64.
        var read1 = Ascii(64, 70, 73);
        var read2 = Ascii(66, 90, 104);          // min 66 (>=64), max 104 (>74)

        var r = DetectEncoding(new[] { read1, read2 });

        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred64));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Inferred));
        Assert.That(r.MinAscii, Is.EqualTo(64));
        Assert.That(r.MaxAscii, Is.EqualTo(104));
    }

    // ---- boundary of the "proven Phred+33" rule: char < 64 ----

    [TestCase(63, EncodingConfidence.Definitive, QualityEncoding.Phred33)] // '?' = 63 < 64 -> proven Phred+33
    [TestCase(64, EncodingConfidence.Ambiguous, QualityEncoding.Phred33)]  // '@' = 64, rest overlap -> ambiguous
    public void LowCharBoundary_At64(int minCode, EncodingConfidence expectedConf, QualityEncoding expectedEnc)
    {
        // minCode plus two overlap-range chars; max stays <= 74 so the only lever is the low char.
        var r = DetectEncoding(new[] { Ascii(minCode, 70, 73) });

        Assert.That(r.Confidence, Is.EqualTo(expectedConf));
        Assert.That(r.Encoding, Is.EqualTo(expectedEnc));
    }

    // ---- boundary of the "inferred Phred+64" rule: char > 74, with no char < 64 ----

    [TestCase(74, EncodingConfidence.Ambiguous, QualityEncoding.Phred33)] // 'J' = 74, NOT > 74 -> still ambiguous
    [TestCase(75, EncodingConfidence.Inferred, QualityEncoding.Phred64)]  // 'K' = 75 > 74 -> inferred Phred+64
    public void HighCharBoundary_At74(int maxCode, EncodingConfidence expectedConf, QualityEncoding expectedEnc)
    {
        // min fixed at 64 (no proof of Phred+33); max is the lever.
        var r = DetectEncoding(new[] { Ascii(64, 70, maxCode) });

        Assert.That(r.Confidence, Is.EqualTo(expectedConf));
        Assert.That(r.Encoding, Is.EqualTo(expectedEnc));
    }

    [Test]
    public void LowCharWins_OverHighChar_Phred33Definitive()
    {
        // Both a sub-64 char (proof of Phred+33) and a >74 char are present: proof wins, never Phred+64.
        var r = DetectEncoding(new[] { Ascii(35, 80, 104) });

        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Definitive));
    }

    // ---- reported span / count ----

    [Test]
    public void ReportsGlobalMinMaxAndCharacterCount()
    {
        var r = DetectEncoding(new[] { Ascii(40, 73), Ascii(35), Ascii(104, 90) });

        Assert.That(r.MinAscii, Is.EqualTo(35));
        Assert.That(r.MaxAscii, Is.EqualTo(104));
        Assert.That(r.CharactersExamined, Is.EqualTo(5));
    }

    // ---- empty / null handling ----

    [Test]
    public void EmptyInput_YieldsAmbiguousPhred33_ZeroCount()
    {
        var r = DetectEncoding(Enumerable.Empty<string>());

        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Ambiguous));
        Assert.That(r.CharactersExamined, Is.EqualTo(0));
        Assert.That(r.MinAscii, Is.EqualTo(0));
        Assert.That(r.MaxAscii, Is.EqualTo(0));
    }

    [Test]
    public void NullAndEmptyStringsAreSkipped()
    {
        // Only the '#'-bearing read contributes; null and "" are ignored, not counted.
        var input = new string?[] { null, "", Ascii(35, 64) };
        var r = DetectEncoding(input!);

        Assert.That(r.CharactersExamined, Is.EqualTo(2));
        Assert.That(r.MinAscii, Is.EqualTo(35));
        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Definitive));
    }

    [Test]
    public void AllEmptyStrings_AreLikeEmptyInput()
    {
        var r = DetectEncoding(new[] { "", "", "" });

        Assert.That(r.CharactersExamined, Is.EqualTo(0));
        Assert.That(r.Confidence, Is.EqualTo(EncodingConfidence.Ambiguous));
        Assert.That(r.Encoding, Is.EqualTo(QualityEncoding.Phred33));
    }

    [Test]
    public void NullEnumerable_Throws()
        => Assert.Throws<ArgumentNullException>(() => DetectEncoding((IEnumerable<string>)null!));

    // ---- single-read parity: file-level encoding == per-read encoding for one read ----

    [TestCase("!I")]      // min 33 < 64
    [TestCase("@K")]      // min 64, max 75 > 74
    [TestCase("@J")]      // min 64, max 74 -> overlap
    [TestCase("@F")]      // min 64, max 70 -> overlap
    [TestCase("?@A")]     // min 63 < 64
    [TestCase("hhhh")]    // min 104 -> max also > 74
    public void SingleReadEncoding_MatchesPerStringDetector(string quality)
    {
        var fileLevel = DetectEncoding(new[] { quality }).Encoding;
        var perString = DetectEncoding(quality);

        Assert.That(fileLevel, Is.EqualTo(perString),
            $"file-level and per-string detection must agree for a single read \"{quality}\"");
    }
}
