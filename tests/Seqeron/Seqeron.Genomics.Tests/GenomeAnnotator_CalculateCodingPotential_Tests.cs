// ANNOT-CODING-001 — Coding Potential Calculation (CPAT hexamer usage-bias score)
// Evidence: docs/Evidence/ANNOT-CODING-001-Evidence.md
// TestSpec: tests/TestSpecs/ANNOT-CODING-001.md
// Source: Wang L et al. (2013). CPAT. Nucleic Acids Res 41(6):e74. https://doi.org/10.1093/nar/gkt006
//         Reference impl: cpmodule/FrameKmer.py kmer_ratio (frame 0).
//         https://github.com/WGLab/lncScore/blob/master/tools/cpmodule/FrameKmer.py

namespace Seqeron.Genomics.Tests;

using System;
using System.Collections.Generic;
using Seqeron.Genomics.Annotation;

[TestFixture]
public class GenomeAnnotator_CalculateCodingPotential_Tests
{
    #region CalculateCodingPotential

    // M1 — Two in-frame hexamers. ATGAAA: ln(8/2)=1.3862943611198906; AAACCC: ln(2/4)=-0.6931471805599453.
    //      sum=0.6931471805599453, count=2 → 0.34657359027997264 (FrameKmer.kmer_ratio).
    [Test]
    public void CalculateCodingPotential_TwoInFrameHexamers_ReturnsMeanLogRatio()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 4 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        Assert.That(score, Is.EqualTo(0.34657359027997264).Within(1e-10),
            "Score must be the mean of ln(8/2) and ln(2/4) per FrameKmer.kmer_ratio.");
    }

    // M2 — Single in-both hexamer: ln(4/1)=ln 4=1.3862943611198906; count=1 → score = that value.
    [Test]
    public void CalculateCodingPotential_SingleHexamerBothTables_ReturnsNaturalLogRatio()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 4 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 1 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, noncoding);

        Assert.That(score, Is.EqualTo(1.3862943611198906).Within(1e-10),
            "A single in-both hexamer scores ln(coding/noncoding) = ln(4) (natural log).");
    }

    // M3 — Coding-only hexamer (noncoding==0) contributes +1 (kmer_ratio elif).
    [Test]
    public void CalculateCodingPotential_CodingOnlyHexamer_ContributesPlusOne()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 5 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 0 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, noncoding);

        Assert.That(score, Is.EqualTo(1.0).Within(1e-10),
            "coding>0 & noncoding==0 contributes +1; single hexamer → score 1.");
    }

    // M4 — Noncoding-only hexamer (coding==0) contributes -1 (kmer_ratio elif).
    [Test]
    public void CalculateCodingPotential_NoncodingOnlyHexamer_ContributesMinusOne()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 0 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 5 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAA", coding, noncoding);

        Assert.That(score, Is.EqualTo(-1.0).Within(1e-10),
            "coding==0 & noncoding>0 contributes -1; single hexamer → score -1.");
    }

    // M5 — Sign invariant: coding-biased tables give a positive score (CPAT interpretation).
    [Test]
    public void CalculateCodingPotential_CodingBiasedTables_ReturnsPositiveScore()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 8 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 8 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        // ATGAAA: ln(4)=1.3862943611198906; AAACCC: ln(1)=0; mean=0.6931471805599453.
        Assert.That(score, Is.EqualTo(0.6931471805599453).Within(1e-10),
            "Coding-biased table must yield a positive score (mean of ln4 and ln1).");
    }

    // M6 — Sign invariant: noncoding-biased tables give a negative score.
    [Test]
    public void CalculateCodingPotential_NoncodingBiasedTables_ReturnsNegativeScore()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 4 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 2 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        // ATGAAA: ln(2/8)=ln(0.25)=-1.3862943611198906; AAACCC: ln(4/2)=ln2=0.6931471805599453.
        // mean = -0.34657359027997264.
        Assert.That(score, Is.EqualTo(-0.34657359027997264).Within(1e-10),
            "Noncoding-biased table must yield a negative score.");
    }

    // M7 — Sequence shorter than wordSize → 0 (kmer_ratio guard: len(seq) < word_size).
    [Test]
    public void CalculateCodingPotential_SequenceShorterThanWord_ReturnsZero()
    {
        var coding = new Dictionary<string, double> { ["ATGAA"] = 5 };
        var noncoding = new Dictionary<string, double> { ["ATGAA"] = 5 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAA", coding, noncoding);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "A 5-nt sequence has no 6-mer; kmer_ratio returns 0.");
    }

    // M8 — A hexamer missing from a table is skipped and not counted; only the in-both hexamer is scored.
    [Test]
    public void CalculateCodingPotential_HexamerMissingFromTable_IsSkipped()
    {
        // ATGAAA in both tables; AAACCC only in coding → skipped.
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 4, ["AAACCC"] = 9 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 1 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        // Only ATGAAA counts: ln(4/1)=1.3862943611198906; count=1.
        Assert.That(score, Is.EqualTo(1.3862943611198906).Within(1e-10),
            "AAACCC absent from noncoding table must be skipped, leaving only ATGAAA scored.");
    }

    // M9 — In-frame stepping: an out-of-frame hexamer (offset not a multiple of 3) is never scored.
    [Test]
    public void CalculateCodingPotential_OutOfFrameHexamer_IsNotScored()
    {
        // Sequence GATGAA AAA (len 9). In-frame hexamers (step 3): i=0 GATGAA, i=3 GAAAAA.
        // The out-of-frame hexamer ATGAAA (offset 1) must NOT be scored even if present in tables.
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 100 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 1 };

        double score = GenomeAnnotator.CalculateCodingPotential("GATGAAAAA", coding, noncoding);

        // GATGAA and GAAAAA are absent from both tables → skipped; ATGAAA is out of frame → never seen.
        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "Only frame-0 hexamers (step 3) are scored; the out-of-frame ATGAAA is ignored.");
    }

    // M10 — Null sequence → ArgumentNullException.
    [Test]
    public void CalculateCodingPotential_NullSequence_Throws()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        Assert.Throws<ArgumentNullException>(
            () => GenomeAnnotator.CalculateCodingPotential(null!, t, t),
            "Null sequence must throw ArgumentNullException.");
    }

    // M11 — Null coding table → ArgumentNullException.
    [Test]
    public void CalculateCodingPotential_NullCodingTable_Throws()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        Assert.Throws<ArgumentNullException>(
            () => GenomeAnnotator.CalculateCodingPotential("ATGAAA", null!, t),
            "Null coding table must throw ArgumentNullException.");
    }

    // M12 — Null noncoding table → ArgumentNullException.
    [Test]
    public void CalculateCodingPotential_NullNoncodingTable_Throws()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        Assert.Throws<ArgumentNullException>(
            () => GenomeAnnotator.CalculateCodingPotential("ATGAAA", t, null!),
            "Null noncoding table must throw ArgumentNullException.");
    }

    // M13 — Non-positive wordSize → ArgumentOutOfRangeException.
    [Test]
    public void CalculateCodingPotential_NonPositiveWordSize_Throws()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        Assert.Throws<ArgumentOutOfRangeException>(
            () => GenomeAnnotator.CalculateCodingPotential("ATGAAA", t, t, wordSize: 0),
            "wordSize <= 0 must throw.");
    }

    // M14 — Non-positive stepSize → ArgumentOutOfRangeException.
    [Test]
    public void CalculateCodingPotential_NonPositiveStepSize_Throws()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        Assert.Throws<ArgumentOutOfRangeException>(
            () => GenomeAnnotator.CalculateCodingPotential("ATGAAA", t, t, stepSize: 0),
            "stepSize <= 0 must throw.");
    }

    // M15 — No scorable hexamer (all hexamers absent from tables) → 0 (port choice; ASSUMPTION-1).
    [Test]
    public void CalculateCodingPotential_NoScorableHexamer_ReturnsZero()
    {
        var coding = new Dictionary<string, double> { ["GGGGGG"] = 5 };
        var noncoding = new Dictionary<string, double> { ["GGGGGG"] = 5 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "When no in-frame hexamer is in the tables, the score is 0 (no information).");
    }

    // S1 — Case-insensitivity: lowercase input gives the same score as uppercase.
    [Test]
    public void CalculateCodingPotential_LowercaseSequence_MatchesUppercase()
    {
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 4 };

        double lower = GenomeAnnotator.CalculateCodingPotential("atgaaaccc", coding, noncoding);
        double upper = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        Assert.Multiple(() =>
        {
            Assert.That(lower, Is.EqualTo(0.34657359027997264).Within(1e-10),
                "Lowercase input is upper-cased internally and scores identically.");
            Assert.That(lower, Is.EqualTo(upper).Within(1e-10),
                "Lowercase and uppercase inputs must produce the same score.");
        });
    }

    // S2 — Empty sequence → 0 (length 0 < wordSize).
    [Test]
    public void CalculateCodingPotential_EmptySequence_ReturnsZero()
    {
        var t = new Dictionary<string, double> { ["ATGAAA"] = 1 };

        double score = GenomeAnnotator.CalculateCodingPotential("", t, t);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "Empty sequence has no hexamer; score is 0.");
    }

    // C1 — In-both hexamer with both values 0: contributes 0 but IS counted (reference fall-through).
    [Test]
    public void CalculateCodingPotential_BothTablesZero_HexamerIsCountedAsZero()
    {
        // ATGAAA: coding=0, noncoding=0 → contributes 0, counted. AAACCC: ln(4/1)=1.3862943611198906, counted.
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 0, ["AAACCC"] = 4 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 0, ["AAACCC"] = 1 };

        double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);

        // sum = 0 + 1.3862943611198906; count = 2 → 0.6931471805599453.
        Assert.That(score, Is.EqualTo(0.6931471805599453).Within(1e-10),
            "A both-zero in-both hexamer adds 0 but is counted, so it halves the AAACCC contribution.");
    }

    #endregion
}
