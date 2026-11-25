// ASSEMBLY-TRIM-001 — Quality Trimming (BWA / cutadapt running-sum)
// Evidence: docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-TRIM-001.md
// Source: Cutadapt algorithm docs (quality trimming); Li H. BWA bwa_trim_read;
//         Cock et al. (2010) NAR 38(6):1767-1771 (Phred+33).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_QualityTrimReads_Tests
{
    #region QualityTrimReads

    // Phred+33 (Cock et al. 2010): char = Phred + 33.
    // Cutadapt worked example qualities 42,40,26,27,8,7,11,4,2,3 -> "KI;<)(,%#$".
    private const string ExampleQuality = "KI;<)(,%#$";

    // M1 — Cutadapt worked example: subtract cutoff 10, partial sums from end
    //      (70),(38),8,-8,-25,-23,-20,-21,-15,-7; minimum -25 at index 4 ->
    //      read trimmed to the first four bases. (Cutadapt algorithm docs.)
    [Test]
    public void QualityTrimReads_CutadaptWorkedExample_TrimsToFirstFourBases()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "single read survives min-length 1");
        Assert.That(result[0], Is.EqualTo("ACGT"),
            "cutadapt example: minimum partial sum -25 at index 4 keeps the first four bases");
    }

    // M2 — All-high-quality read: every (q - cutoff) >= 0, partial sums never go
    //      below 0, minimum is at the end -> nothing trimmed.
    [Test]
    public void QualityTrimReads_AllHighQuality_ReturnsUnchanged()
    {
        // Arrange: 'I' = Phred 40; cutoff 20.
        var reads = new List<(string, string)> { ("ACGTAC", "IIIIII") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "high-quality read is kept");
        Assert.That(result[0], Is.EqualTo("ACGTAC"),
            "all qualities above cutoff -> no bases removed");
    }

    // M3 — All-low-quality read: every (q - cutoff) < 0, both end passes consume the
    //      read, trimmed length 0 < minLength -> dropped.
    [Test]
    public void QualityTrimReads_AllLowQuality_DropsRead()
    {
        // Arrange: '!' = Phred 0; cutoff 20.
        var reads = new List<(string, string)> { ("ACGTAC", "!!!!!!") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        // Assert
        Assert.That(result, Is.Empty,
            "all qualities below cutoff -> read trimmed to length 0 and dropped");
    }

    // M4 — Threshold below 1 disables trimming (BWA bwa_trim_read: trim_qual < 1 -> 0).
    [Test]
    public void QualityTrimReads_ThresholdZero_ReturnsUnchanged()
    {
        // Arrange: mixed qualities that WOULD trim under a positive cutoff.
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 0, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "read kept");
        Assert.That(result[0], Is.EqualTo("ACGTACGTAC"),
            "cutoff < 1 disables trimming per BWA guard");
    }

    // M5 — Min-length filter: the cutadapt example trims to 4 bases; with minLength 5
    //      the survivor (length 4) is dropped.
    [Test]
    public void QualityTrimReads_TrimmedShorterThanMinLength_DropsRead()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 5);

        // Assert
        Assert.That(result, Is.Empty,
            "trimmed length 4 < minLength 5 -> read dropped");
    }

    // S1 — 5'-end trimming: low-quality prefix, high-quality suffix.
    //      Qualities 3,2,4,11,7,8,27,26,40,42 ("$#%,()<;IK"), cutoff 10.
    //      3' pass: partial sums from end never negative -> no 3' trim.
    //      5' pass: partial sums from start -7,-15,-21,-20,-23,-25,...; minimum -25 at
    //      index 5 -> cut after index 5 -> keep the last four bases.
    [Test]
    public void QualityTrimReads_LowQualityPrefix_TrimsFivePrimeEnd()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", "$#%,()<;IK") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "read survives");
        Assert.That(result[0], Is.EqualTo("GTAC"),
            "5' running-sum minimum at index 5 keeps the last four (high-quality) bases");
    }

    // S2 — Min-length boundary: trimmed length exactly equals minLength -> kept.
    [Test]
    public void QualityTrimReads_TrimmedEqualsMinLength_KeepsRead()
    {
        // Arrange: cutadapt example trims to 4; minLength 4.
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 4);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "trimmed length 4 == minLength 4 -> kept");
        Assert.That(result[0], Is.EqualTo("ACGT"),
            "boundary: length-4 survivor retained at minLength 4");
    }

    // S3 — Multiple reads trimmed independently: one trimmable, one droppable, one clean.
    [Test]
    public void QualityTrimReads_MultipleReads_TrimsEachIndependently()
    {
        // Arrange
        var reads = new List<(string, string)>
        {
            ("ACGTACGTAC", ExampleQuality), // -> "ACGT"
            ("ACGTAC", "!!!!!!"),           // all-low -> dropped
            ("ACGTAC", "IIIIII"),           // all-high -> unchanged
        };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(2), "the all-low read is dropped");
            Assert.That(result[0], Is.EqualTo("ACGT"), "first read trimmed to four bases");
            Assert.That(result[1], Is.EqualTo("ACGTAC"), "high-quality read unchanged, order preserved");
        });
    }

    // C4 — Refinement / s<0 early break: a high-quality base at the 3' end protects the
    //      bases to its 5' side from the 3' pass. Qualities 20,0,8,40 ("5!)I"), cutoff 10.
    //      cutadapt/BWA bwa_trim_read: the 3' pass accumulates (cutoff - q); at index 3 the
    //      sum is 10-40 = -30 < 0 so it breaks immediately (stop stays 4); the 5' pass breaks
    //      at index 0 (10-20 = -10 < 0, start stays 0). Result (0,4): NOTHING is trimmed.
    //      A global-minimum scan WITHOUT the early break would over-trim to the last base.
    //      Source: cutadapt quality_trim_index (qualtrim.pyx); BWA bwa_trim_read `if (s<0) break`.
    [Test]
    public void QualityTrimReads_GoodBaseAt3PrimeEnd_KeepsLowQualityInterior()
    {
        // Arrange: Phred 20,0,8,40 = ASCII "5!)I".
        var reads = new List<(string, string)> { ("ACGT", "5!)I") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "read survives");
        Assert.That(result[0], Is.EqualTo("ACGT"),
            "high-quality 3' base triggers the s<0 break; no bases are trimmed (cutadapt/BWA)");
    }

    // C5 — "A few good-quality bases among the bad ones": qualities 2,2,40,2,2 ("##I##"),
    //      cutoff 10. 3' pass: i4 -> cutoff-q=8 (s=8,stop=4), i3 -> +8 (s=16,stop=3),
    //      i2 -> 10-40=-30 (s=-14<0) break => stop=3. 5' pass: i0 -> 8 (s=8,start=1),
    //      i1 -> 8 (s=16,start=2), i2 -> -30 (s=-14<0) break => start=2. Result (2,3):
    //      only the single isolated high-quality base is kept.
    //      Source: cutadapt algorithm docs ("allowing some good-quality bases among the bad ones").
    [Test]
    public void QualityTrimReads_IsolatedGoodBaseInBadRead_KeepsOnlyGoodBase()
    {
        // Arrange: Phred 2,2,40,2,2 = ASCII "##I##".
        var reads = new List<(string, string)> { ("ACGTA", "##I##") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "the isolated good base survives");
        Assert.That(result[0], Is.EqualTo("G"),
            "both passes break at the high-quality base, keeping only index 2");
    }

    // C6 — start >= stop => empty good-quality segment (read dropped). cutadapt
    //      quality_trim_index final guard `if start >= stop: start, stop = 0, 0`.
    //      Qualities 20,5,5 ("5&&"), cutoff 20. 5' pass: i0 -> 20-20=0 (s=0, start stays 0),
    //      i1 -> 20-5=15 (s=15,start=2), i2 -> 15 (s=30,start=3) => start=3. 3' pass:
    //      i2 -> 15 (s=15,stop=2), i1 -> 15 (s=30,stop=1), i0 -> 0 (s=30) => stop=1.
    //      start(3) >= stop(1) -> (0,0) -> dropped.
    [Test]
    public void QualityTrimReads_WindowsCross_DropsRead()
    {
        // Arrange: Phred 20,5,5 = ASCII "5&&".
        var reads = new List<(string, string)> { ("ACG", "5&&") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        // Assert
        Assert.That(result, Is.Empty,
            "start >= stop after both passes -> empty good-quality segment -> read dropped");
    }

    // F-TRIM-1 — Quality decoding goes through the canonical Core FASTQ decoder (QualitySequence,
    //      Cock et al. 2010 / Biopython QualityIO): a quality string whose length differs from the
    //      sequence, or a character outside ['!', '~'], is malformed FASTQ. Previously a short quality
    //      string threw IndexOutOfRangeException and a long one was silently ignored; '\u001f' decoded to
    //      Phred -2. cutadapt/dnaio reject the length mismatch ("lengths differ").
    [TestCase("ACGT", "II")]
    [TestCase("ACGT", "IIIIII")]
    [TestCase("ACGT", "II\u001fI")]
    [TestCase("ACGT", "II\u007fI")]
    public void QualityTrimReads_MalformedQuality_ThrowsArgumentException(string sequence, string quality)
    {
        var reads = new List<(string, string)> { (sequence, quality) };

        Assert.That(() => SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1),
            NUnit.Framework.Throws.InstanceOf<ArgumentException>());
    }

    // F-TRIM-1 — a null sequence or quality inside the batch is rejected (was NullReferenceException).
    [Test]
    public void QualityTrimReads_NullSequenceOrQuality_ThrowsArgumentNullException()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => SequenceAssembler.QualityTrimReads(
                    new List<(string, string)> { (null!, "II") }, 10, 1),
                NUnit.Framework.Throws.InstanceOf<ArgumentNullException>());
            Assert.That(() => SequenceAssembler.QualityTrimReads(
                    new List<(string, string)> { ("AC", null!) }, 10, 1),
                NUnit.Framework.Throws.InstanceOf<ArgumentNullException>());
        });
    }

    // F-TRIM-2 / F-TRIM-3 — cutadapt quality_trim_index(qualities, cutoff_front, cutoff_back, base)
    //      (cutadapt 5.2, src/cutadapt/qualtrim.pyx), values produced by the installed cutadapt 5.2:
    //      separate 5'/3' cutoffs (`-q 5',3'`; a single `-q C` = (0, C), 3' only, as BWA) and
    //      Phred+64 (`--quality-base 64`). "jhZ[GFJCAB" is the cutadapt worked example in Phred+64.
    [TestCase("KI;<)(,%#$", 10, 10, 33, 0, 4)]
    [TestCase("KI;<)(,%#$", 0, 10, 33, 0, 4)]
    [TestCase("$#%,()<;IK", 10, 0, 33, 6, 10)]
    [TestCase("$#%,()<;IK", 0, 10, 33, 0, 10)]
    [TestCase("##I##", 10, 10, 33, 2, 3)]
    [TestCase("##I##", 0, 10, 33, 0, 3)]
    [TestCase("##I##", 10, 0, 33, 2, 5)]
    [TestCase("5!)I", 10, 10, 33, 0, 4)]
    [TestCase("5&&", 20, 20, 33, 0, 0)]
    [TestCase("jhZ[GFJCAB", 10, 10, 64, 0, 4)]
    [TestCase("jhZ[GFJCAB", 10, 0, 64, 0, 10)]
    [TestCase("", 10, 10, 33, 0, 0)]
    public void QualityTrimIndex_MatchesCutadaptQualityTrimIndex(
        string quality, int cutoffFront, int cutoffBack, int phredOffset, int expectedStart, int expectedStop)
    {
        var sequence = new string('A', quality.Length);
        var decoded = new QualitySequence(sequence, quality, phredOffset);

        var (start, stop) = SequenceAssembler.QualityTrimIndex(decoded.Qualities, cutoffFront, cutoffBack);
        IReadOnlyList<string> reads = SequenceAssembler.QualityTrimReads(
            new List<(string, string)> { (sequence, quality) }, cutoffFront, cutoffBack, minLength: 0, phredOffset);

        Assert.Multiple(() =>
        {
            Assert.That((start, stop), Is.EqualTo((expectedStart, expectedStop)), "cutadapt quality_trim_index");
            Assert.That(reads, Has.Count.EqualTo(1));
            Assert.That(reads[0], Has.Length.EqualTo(expectedStop - expectedStart));
        });
    }

    // F-TRIM-2 — the single-cutoff overload equals cutoffFront = cutoffBack (cutadapt `-q C,C`),
    //      and the five-argument overload with cutoffFront 0 is the cutadapt/BWA default 3'-only trim.
    [Test]
    public void QualityTrimReads_SeparateCutoffs_ThreePrimeOnlyKeepsLowQualityPrefix()
    {
        var reads = new List<(string, string)> { ("ACGTACGTAC", "$#%,()<;IK") };

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.QualityTrimReads(reads, 10, 1),
                Is.EqualTo(SequenceAssembler.QualityTrimReads(reads, 10, 10, 1)));
            Assert.That(SequenceAssembler.QualityTrimReads(reads, cutoffFront: 0, cutoffBack: 10, minLength: 1),
                Is.EqualTo(new[] { "ACGTACGTAC" }), "3'-only (cutadapt -q 10): the low-quality 5' run is kept");
        });
    }

    // F-TRIM-3 — Phred+64: the cutadapt example "jhZ[GFJCAB" decodes to 42,40,26,27,8,7,11,4,2,3
    //      and trims to the first four bases, case of the input sequence preserved; under the default
    //      Phred+33 the same string decodes to 73..33 (all >= 10) and nothing is trimmed.
    [Test]
    public void QualityTrimReads_Phred64_DecodesWithOffset64()
    {
        var reads = new List<(string, string)> { ("acgtACGTac", "jhZ[GFJCAB") };

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.QualityTrimReads(reads, 10, 10, 1, phredOffset: 64),
                Is.EqualTo(new[] { "acgt" }));
            Assert.That(SequenceAssembler.QualityTrimReads(reads, 10, 1),
                Is.EqualTo(new[] { "acgtACGTac" }));
            Assert.That(() => SequenceAssembler.QualityTrimReads(
                    new List<(string, string)> { ("ACGT", "5555") }, 10, 10, 1, phredOffset: 64),
                NUnit.Framework.Throws.InstanceOf<ArgumentException>(), "'5' (ASCII 53) is below the Phred+64 range");
        });
    }

    // C1 — Empty reads list -> empty result.
    [Test]
    public void QualityTrimReads_EmptyList_ReturnsEmpty()
    {
        // Arrange
        var reads = new List<(string, string)>();

        // Act
        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads);

        // Assert
        Assert.That(result, Is.Empty, "no reads -> empty output");
    }

    // C2 — Empty sequence/quality read -> length 0 < minLength (default) -> dropped.
    [Test]
    public void QualityTrimReads_EmptySequence_DropsRead()
    {
        // Arrange
        var reads = new List<(string, string)> { (string.Empty, string.Empty) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Is.Empty, "empty read has length 0 < minLength -> dropped");
    }

    // C3 — Null reads argument -> ArgumentNullException.
    [Test]
    public void QualityTrimReads_NullReads_Throws()
    {
        // Act / Assert
        Assert.That(
            () => SequenceAssembler.QualityTrimReads(null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "null reads collection is rejected");
    }

    #endregion
}
