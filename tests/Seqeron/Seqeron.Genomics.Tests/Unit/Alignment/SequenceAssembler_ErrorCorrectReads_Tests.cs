// ASSEMBLY-CORRECT-001 — K-mer spectrum (two-sided) read error correction
// Evidence: docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-CORRECT-001.md
// Source: Liu, Schmidt & Maskell (2013). Musket. Bioinformatics 29(3):308-315.
//         Kelley, Schatz & Salzberg (2010). Quake. Genome Biol 11:R116.

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_ErrorCorrectReads_Tests
{
    #region ErrorCorrectReads

    // Worked example (k=3, cut-off=2): three true reads ACGTACGT plus one read with a
    // single substitution A->T at index 4. Spectrum: ACG=7,CGT=8,GTA=3,TAC=3,GTT=1,TTC=1,
    // TCG=1; trusted (>=2): ACG,CGT,GTA,TAC. Only index 4 of the error read is covered
    // solely by untrusted k-mers; the unique base making GTA,TAC,ACG trusted is A.
    private const string TrueRead = "ACGTACGT";
    private const string ErrorRead = "ACGTTCGT"; // A->T at index 4

    // M1 — Single-substitution read is corrected to the trusted sequence (Musket two-sided
    //      rule: unique alternative base making all covering k-mers trusted).
    [Test]
    public void ErrorCorrectReads_SingleSubstitution_CorrectsToTrustedSequence()
    {
        // Arrange
        var reads = new List<string> { TrueRead, TrueRead, TrueRead, ErrorRead };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(4), "one corrected read per input read");
            Assert.That(result[3], Is.EqualTo(TrueRead),
                "index-4 error corrected: A is the unique base making GTA,TAC,ACG trusted");
            Assert.That(result[0], Is.EqualTo(TrueRead), "true read unchanged");
            Assert.That(result[1], Is.EqualTo(TrueRead), "true read unchanged");
            Assert.That(result[2], Is.EqualTo(TrueRead), "true read unchanged");
        });
    }

    // M2 — Every position of identical reads is covered by a trusted k-mer; the trusted-base
    //      rule forbids modification (Musket/Quake).
    [Test]
    public void ErrorCorrectReads_AllReadsTrusted_LeavesUnchanged()
    {
        // Arrange: 3 identical reads -> all k-mers have multiplicity 3 (>=2 trusted).
        var reads = new List<string> { TrueRead, TrueRead, TrueRead };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result[0], Is.EqualTo(TrueRead), "trusted base never modified");
            Assert.That(result[1], Is.EqualTo(TrueRead), "trusted base never modified");
            Assert.That(result[2], Is.EqualTo(TrueRead), "trusted base never modified");
        });
    }

    // M3 — Ambiguous position (k=1, cut-off=2): 1-mers A=2,C=2 are both trusted; read "T"
    //      has two valid alternatives (A and C) so it is left unchanged (Musket ambiguity).
    [Test]
    public void ErrorCorrectReads_AmbiguousAlternative_LeavesBaseUnchanged()
    {
        // Arrange
        var reads = new List<string> { "A", "A", "C", "C", "G", "T" };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 1, minKmerFrequency: 2);

        // Assert
        Assert.That(result[5], Is.EqualTo("T"),
            "T has two trusted alternatives (A,C) -> ambiguous -> unchanged");
    }

    // M4 — No valid alternative (Quake: no correcting set). The error read carries TWO adjacent
    //      substitutions (pos 2,3) inside the same window, so NO single base substitution makes all
    //      k-mers covering an untrusted position trusted -> the position is left unchanged.
    //      Genome reads "AAAAAAAA" x3 make only the 4-mer AAAA trusted (mult=16>=2). In the error
    //      read "AACCAAAA" positions 2 and 3 are each covered solely by untrusted 4-mers, and for
    //      each the candidate search yields ZERO valid bases (the other error keeps a covering k-mer
    //      weak). Verified independently by hand and by a reference reimplementation:
    //          pos 2 valid candidates = {}; pos 3 valid candidates = {} -> unchanged.
    //      NOTE: this case must NOT use a read that is itself fully trusted (e.g. "TTTTTTTT", where
    //      TTT is trusted) — that would exercise the trusted-base rule (M2), not the no-correction
    //      branch tested here.
    [Test]
    public void ErrorCorrectReads_NoTrustedAlternative_LeavesBaseUnchanged()
    {
        // Arrange
        const string genomeRead = "AAAAAAAA";
        const string twoErrorRead = "AACCAAAA"; // two substitutions at indices 2 and 3
        var reads = new List<string> { genomeRead, genomeRead, genomeRead, twoErrorRead };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 4, minKmerFrequency: 2);

        // Assert
        Assert.That(result[3], Is.EqualTo(twoErrorRead),
            "untrusted positions 2,3 have no single substitution making all covering k-mers trusted -> unchanged");
    }

    // M5 — Substitution-only model: output count and per-read length always preserved
    //      (INV-1, INV-2).
    [Test]
    public void ErrorCorrectReads_Always_PreservesCountAndLength()
    {
        // Arrange
        var reads = new List<string> { TrueRead, TrueRead, TrueRead, ErrorRead };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(reads.Count), "count preserved (INV-1)");
            for (int i = 0; i < reads.Count; i++)
            {
                Assert.That(result[i].Length, Is.EqualTo(reads[i].Length),
                    $"read {i} length preserved (no indels, INV-2)");
            }
        });
    }

    // M6 — Null reads -> ArgumentNullException (documented failure mode).
    [Test]
    public void ErrorCorrectReads_NullReads_Throws()
    {
        Assert.That(() => SequenceAssembler.ErrorCorrectReads(null!),
            NUnit.Framework.Throws.ArgumentNullException, "null reads is invalid input");
    }

    // M7 — kmerSize < 1 -> ArgumentOutOfRangeException (documented failure mode).
    [Test]
    public void ErrorCorrectReads_KmerSizeBelowOne_Throws()
    {
        var reads = new List<string> { TrueRead };
        Assert.That(() => SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 0),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "k must be at least 1");
    }

    // S1 — Error-free reads pass through unchanged (trusted reads not altered).
    [Test]
    public void ErrorCorrectReads_ErrorFreeReads_PassThroughUnchanged()
    {
        // Arrange
        var reads = new List<string> { TrueRead, TrueRead };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result[0], Is.EqualTo(TrueRead), "no error -> unchanged");
            Assert.That(result[1], Is.EqualTo(TrueRead), "no error -> unchanged");
        });
    }

    // S2 — Case-insensitive: a lowercase erroneous read is corrected and upper-cased.
    [Test]
    public void ErrorCorrectReads_LowercaseInput_CorrectsAndUpperCases()
    {
        // Arrange: same spectrum as M1; the error read is lowercase.
        var reads = new List<string> { TrueRead, TrueRead, TrueRead, "acgttcgt" };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.That(result[3], Is.EqualTo(TrueRead),
            "lowercase read is upper-cased and the index-4 error corrected to A");
    }

    // S3 — Determinism: identical inputs produce identical outputs (INV-5).
    [Test]
    public void ErrorCorrectReads_RunTwice_IsDeterministic()
    {
        // Arrange
        var reads = new List<string> { TrueRead, TrueRead, TrueRead, ErrorRead };

        // Act
        IReadOnlyList<string> first =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
        IReadOnlyList<string> second =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.That(second, Is.EqualTo(first), "same inputs -> same output (INV-5)");
    }

    // C1 — Read shorter than k contributes no k-mers and is returned unchanged (upper-cased).
    [Test]
    public void ErrorCorrectReads_ReadShorterThanK_ReturnsUnchanged()
    {
        // Arrange: read length 3 < k=5.
        var reads = new List<string> { "acg" };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 5, minKmerFrequency: 2);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1), "count preserved");
            Assert.That(result[0], Is.EqualTo("ACG"),
                "read shorter than k -> no k-mers cover it -> returned upper-cased unchanged");
        });
    }

    // P1 — Property (INV-2/INV-3, fixed seed): for randomly generated reads with injected
    //      single substitutions, every output read keeps its length, and any position that
    //      was already covered by a trusted k-mer in the input is never changed. O(n*r*k^2)
    //      algorithm -> property-based invariant per Definition of Done.
    [Test]
    public void ErrorCorrectReads_RandomReads_PreservesLengthAndTrustedPositions()
    {
        // Arrange: deterministic generator (fixed seed) -> reproducible.
        const int kmerSize = 7;
        const int cutoff = 2;
        var rng = new Random(20260613);
        const string bases = "ACGT";
        string genome = string.Concat(System.Linq.Enumerable.Range(0, 200)
            .Select(_ => bases[rng.Next(bases.Length)]));

        var reads = new List<string>();
        for (int i = 0; i < 60; i++)
        {
            int start = rng.Next(genome.Length - 40);
            var read = new System.Text.StringBuilder(genome.Substring(start, 40));
            if (rng.NextDouble() < 0.3) // inject a single substitution into some reads
            {
                int p = rng.Next(read.Length);
                char c;
                do { c = bases[rng.Next(bases.Length)]; } while (c == read[p]);
                read[p] = c;
            }
            reads.Add(read.ToString());
        }

        // Build the same k-mer spectrum the algorithm uses, to identify trusted positions.
        var spectrum = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (string r in reads)
        {
            for (int i = 0; i + kmerSize <= r.Length; i++)
            {
                string km = r.Substring(i, kmerSize);
                spectrum[km] = spectrum.GetValueOrDefault(km, 0) + 1;
            }
        }

        bool IsTrustedPos(string r, int pos)
        {
            int first = Math.Max(0, pos - kmerSize + 1);
            int last = Math.Min(pos, r.Length - kmerSize);
            for (int s = first; s <= last; s++)
            {
                if (spectrum.GetValueOrDefault(r.Substring(s, kmerSize), 0) >= cutoff)
                {
                    return true;
                }
            }

            return false;
        }

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(reads, kmerSize, cutoff);

        // Assert
        Assert.Multiple(() =>
        {
            for (int i = 0; i < reads.Count; i++)
            {
                Assert.That(result[i].Length, Is.EqualTo(reads[i].Length),
                    $"read {i} length preserved (INV-2)");
                for (int pos = 0; pos < reads[i].Length; pos++)
                {
                    if (IsTrustedPos(reads[i], pos))
                    {
                        Assert.That(result[i][pos], Is.EqualTo(reads[i][pos]),
                            $"read {i} pos {pos} was trusted in input -> must be unchanged (INV-3)");
                    }
                }
            }
        });
    }

    // C2 — Empty read list -> empty result.
    [Test]
    public void ErrorCorrectReads_EmptyReadList_ReturnsEmpty()
    {
        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(new List<string>(), kmerSize: 3, minKmerFrequency: 2);

        // Assert
        Assert.That(result, Is.Empty, "no reads -> empty result");
    }

    #endregion

    #region Musket 1.1 port (B10 review 2026-11: F31-F35)

    // F31 — Musket two-sided rule probes only the LEFTMOST and RIGHTMOST k-mer covering position i
    //       (paper: "evaluates both the leftmost and the rightmost k-mers that cover position i";
    //       musket-1.1 paracorrect.cpp paraCorrectCoreTwoSides). k=3, cut-off 2, reads TCTCT, TCGCT:
    //       forward spectrum TCT=2 (trusted), CTC/TCG/CGC/GCT=1. In TCGCT no base is trusted; index 2
    //       (G): leftmost TCG@0 -> {T: TCT}, rightmost GCT@2 -> {T: TCT}; intersection {T} -> G->T.
    //       The middle k-mer CGC->CTC stays untrusted (1); the former "all covering k-mers" rule
    //       rejected this correction.
    [Test]
    public void ErrorCorrectReads_TwoSidedRule_ProbesOnlyLeftmostAndRightmostKmers()
    {
        IReadOnlyList<string> result =
            SequenceAssembler.ErrorCorrectReads(new[] { "TCTCT", "TCGCT" }, kmerSize: 3, minKmerFrequency: 2);

        Assert.That(result, Is.EqualTo(new[] { "TCTCT", "TCTCT" }));
    }

    // F32 — Musket counts canonical k-mers (main.cpp: rep = min(km, twin)). Three reverse-strand copies
    //       of GATTACAGGCTTC plus one forward read with C->G at index 6. Forward-only spectrum: the read's
    //       genomic k-mers are never seen on the forward strand, while AAT (from ...GTAATC) is, so the
    //       untrusted GAT at index 0 is "corrected" to AAT (false correction) and the real error stays.
    //       Canonical spectrum: every genomic k-mer is trusted (3), index 6 has the unique fix A.
    [Test]
    public void ErrorCorrectReads_ReverseStrandReads_CanonicalSpectrumCorrects_ForwardSpectrumMiscorrects()
    {
        const string reverseCopy = "GAAGCCTGTAATC"; // reverse complement of GATTACAGGCTTC
        var reads = new[] { reverseCopy, reverseCopy, reverseCopy, "GATTACGGGCTTC" };

        // A9-1 default flip: the 3-argument overload counts canonical k-mers (Musket); the forward-only
        // spectrum is the explicit opt-out.
        IReadOnlyList<string> canonical = SequenceAssembler.ErrorCorrectReads(reads, 3, 2);
        IReadOnlyList<string> forward = SequenceAssembler.ErrorCorrectReads(reads, 3, 2,
            SequenceAssembler.ErrorCorrectionOptions.TwoSidedForwardStrand);
        IReadOnlyList<string> musket = SequenceAssembler.ErrorCorrectReads(reads, 3, 2,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.Multiple(() =>
        {
            Assert.That(forward[3], Is.EqualTo("AATTACGGGCTTC"), "forward-only spectrum (explicit opt-out, documented)");
            Assert.That(canonical[3], Is.EqualTo("GATTACAGGCTTC"), "canonical spectrum (default) corrects index 6");
            Assert.That(canonical, Is.EqualTo(SequenceAssembler.ErrorCorrectReads(reads, 3, 2,
                SequenceAssembler.ErrorCorrectionOptions.TwoSided)), "3-arg overload = TwoSided preset");
            Assert.That(musket[3], Is.EqualTo("GATTACAGGCTTC"), "Musket preset");
            Assert.That(musket.Take(3), Is.All.EqualTo(reverseCopy), "trusted reads unchanged");
        });
    }

    // F33 — Musket stage 2 (aggressive one-sided correction, -maxerr). M4 data: two adjacent errors in one
    //       k-mer are beyond the two-sided stage, but the right extension of the trusted run AAAA@4
    //       leftwards fixes C@3 then C@2 (each the unique trusted alternative A). Two corrections within k
    //       bases need -maxerr >= 2 (musket paraCorrectCore numCorrectionsPerKmer > maxErrorPerKmer).
    [TestCase(0, "AACCAAAA")]
    [TestCase(1, "AACCAAAA")]
    [TestCase(2, "AAAAAAAA")]
    [TestCase(4, "AAAAAAAA")]
    public void ErrorCorrectReads_AggressiveStage_TwoErrorsInOneKmer_NeedMaxErrorsTwo(int maxErrors, string expected)
    {
        var reads = new[] { "AAAAAAAA", "AAAAAAAA", "AAAAAAAA", "AACCAAAA" };

        IReadOnlyList<string> result = SequenceAssembler.ErrorCorrectReads(reads, 4, 2,
            SequenceAssembler.ErrorCorrectionOptions.Musket with { MaxErrorsPerKmer = maxErrors });

        Assert.That(result[3], Is.EqualTo(expected));
    }

    // Musket -maxiter: stage 1 applies ONE substitution per iteration (paraCorrectCoreTwoSides returns at
    // the first unique correction). Read with errors at 3 (G->T) and 14 (C->A), k=4, forward, cut-off 2.
    [TestCase(1, "ACGGTCATTGCAGTACAT")]
    [TestCase(2, "ACGGTCATTGCAGTCCAT")]
    public void ErrorCorrectReads_TwoSidedStage_OneCorrectionPerIteration(int iterations, string expected)
    {
        const string genome = "ACGGTCATTGCAGTCCAT";
        var reads = new[] { genome, genome, genome, "ACGTTCATTGCAGTACAT" };

        IReadOnlyList<string> result = SequenceAssembler.ErrorCorrectReads(reads, 4, 2,
            SequenceAssembler.ErrorCorrectionOptions.TwoSidedForwardStrand with { MaxIterations = iterations });

        Assert.That(result[3], Is.EqualTo(expected));
    }

    // Musket stage 3 (voteErroneousRead + fixErroneousBase, minVote 3): read 0 carries two adjacent errors
    // that stages 1-2 leave; voting fixes them. Oracle: Python transcription of musket-1.1 (Evidence).
    [Test]
    public void ErrorCorrectReads_VotingStage_FixesWhatOneSidedStageLeaves()
    {
        var reads = new[] { "ACAATCTTCA", "ACAATTGTCA", "ACTATTGTCA" };

        IReadOnlyList<string> result = SequenceAssembler.ErrorCorrectReads(reads, 3, 2,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.That(result, Is.All.EqualTo("ACAATTGTCA"));
    }

    // Deviation (b) from musket-1.1: predecessor() look-ahead tests the k-mers that really precede the
    // corrected k-mer. k=3, canonical, cut-off 2, reads CGCGT, ACCGT: left extension of CGT@2 meets CCG@1
    // with two trusted alternatives G (GCG=2) and A (ACG=2); the preceding k-mers AGC / AAC are untrusted,
    // so neither is accepted and ACCGT stays. musket-1.1 calls predecessor(seq, ipos-1, ..), sees an
    // empty prefix, accepts G and then C@0 (-> CGCGT).
    [Test]
    public void ErrorCorrectReads_AggressiveStage_LookAheadUsesTruePrecedingKmers()
    {
        IReadOnlyList<string> result = SequenceAssembler.ErrorCorrectReads(new[] { "CGCGT", "ACCGT" }, 3, 2,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.That(result, Is.EqualTo(new[] { "CGCGT", "ACCGT" }));
    }

    // Oracle table: Python line-by-line transcription of musket-1.1 ParaCorrect (paracorrect.cpp,
    // voting.cpp) with deviations (a)/(b); random 30-bp genomes, 8 reads of 15 bp from both strands with
    // 0-2 substitutions, k=5, cut-off 2 (seed 2026). Musket = canonical, -maxiter 2, -maxerr 4;
    // Default = three-argument overload (two-sided stage to convergence, canonical spectrum; oracle
    // canonical=True, max_iters=10**9, num_errors=0); Forward = TwoSidedForwardStrand (canonical=False).
    [TestCase("Musket", "TACCCGGCATATCTT,CCTCTATCCAGCAAA,GACCTCTACCCAGGA,TGCTGGGTAGAGGTC,CAAAGATATGCTGGG,CTTTCAAAGATATGC,CTTTCAAAGACTTGC,AAAGATATTCTGGGT", "TACCCAGCATATCTT,CCTCTACCCAGCAAA,GACCTCTACCCAGCA,TGCTGGGTAGAGGTC,CAAAGATATGCTGGG,CTTTCAAAGATATGC,CTTTCAAAGATATGC,AAAGATATGCTGGGT")]
    [TestCase("Musket", "CATGCGTTTTAGGTA,AGACATACCTAAAAC,AAACGCAAGTTCCGA,CCTAAAACGCAAGTT,CTAGCGTTTTAGGTA,TGTGTTTTAGGGATG,ACCTAAAACGCAAGT,ATACCTAAAACGCAA", "CTTGCGTTTTAGGTA,AGACATACCTAAAAC,AAACGCAAGTTCCGA,CCTAAAACGCAAGTT,CTTGCGTTTTAGGTA,TGCGTTTTAGGTATG,ACCTAAAACGCAAGT,ATACCTAAAACGCAA")]
    [TestCase("Musket", "CACCTAAGCATACCT,AAGGTATGCATAGGT,CACCTAAGCATACCT,GCATACCTTCTGTTT,GCATACCTTTTGTTT,ACCTAAGCATACCTT,CCAGCTAAGCCTACC,AGCATACCTTCTGTT", "CACCTAAGCATACCT,AAGGTATGCATAGGT,CACCTAAGCATACCT,GCATACCTTCTGTTT,GCATACCTTCTGTTT,ACCTAAGCATACCTT,CCACCTAAGCATACC,AGCATACCTTCTGTT")]
    [TestCase("Musket", "GTTCCCTGCGGCAAC,GTGGTTCCCTGTGCA,CACAGGGAACAAACT,CCGGTGCCAACTTGG,CAAGTTGGCACAGGG,GTTGGCACAAGGAAC,CAGGTAGTTCCCTGT,CCCTCGGCCAACTTG", "GTTCCCTGTGCCAAC,GTGGTTCCCTGTGCC,CACAGGGAACAAACT,CCTGTGCCAACTTGG,CAAGTTGGCACAGGG,GTTGGCACAGGGAAC,CAGGTAGTTCCCTGT,CCCTGTGCCAACTTG")]
    [TestCase("Musket", "GCGTGCACTGAGCCC,GTGCAATGAGACCCG,CGCTTTGGAGTGCAA,TTGCAAGCCAAAGAT,CGGGGCTCATTGCAC,ACGCCAAAGAGGTCC,TTGGCCTGTAATGAG,GGGCCTCTTTGGCGT", "GAGTGCAATGAGCCC,GTGCAATGAGCCCCG,CTCTTTGGAGTGCAA,TTGCAAGCCAAAGAG,CGGGGCTCATTGCAC,ACGCCAAAGAGGTCC,TTGGCCTGCAATGAG,AGGCCTCTTTGGCGT")]
    [TestCase("Musket", "AACGCTCCGATCAAG,TTGACCGGAGCGTTG,AGGTCCTTGCTCGGA,AAGTACCGGGGGTGT,ACGCTCCCATCAAGT,CCCCAGGTACTTGAT,CTCCGATCAAATACT,AGGTACTTGATCGGA", "AACGCTCCGATCAAG,TTGATCGGAGCGTTG,AGGTACTTGATCGGA,AAGTACCGGTGGTGT,ACGCTCCGATCAAGT,CCCCAGGTACTTGAT,CTCCGATCAAGTACT,AGGTACTTGATCGGA")]
    [TestCase("Default", "GCAAGGCCGGCTCCT,GGTCCCCCAAGGCCG,CATGGCCGGCTCCAC,TGGGGGACCGACGGA,CTTGGAGGACCGACG,TCGGTCCCCCAAGGC,CCAAGGCCGGCTCCA,GCTCTGGGGGACCGA", "CCAAGGCCGGCTCCA,GGTCCCCCAAGGCCG,CAAGGCCGGCTCCAC,TGGGGGACCGACGGA,CTTGGAGGACCGACG,TCGGTCCCCCAAGGC,CCAAGGCCGGCTCCA,GCTCTGGGGGACCGA")]
    [TestCase("Default", "ACAAAGCACCTCAAC,TACCAATCACTTCTA,TTGCGTTACCAACCA,GTGCTTGGGAACGCA,ACGTTACCAAGCACC,TGCTTGTTACCGCAA,CACCTCTACGGGTAT,CGTTACCAAGCACCT", "ACCAAGCACCTCAAC,TACCAATCACTTCTA,TTGCGTTACCAAGCA,GTGCTTGGTAACGCA,GCGTTACCAAGCACC,TGCTTGTTACCGCAA,CACCTCTACGGGTAA,CGTTACCAAGCACCT")]
    [TestCase("Forward", "GCAAGGCCGGCTCCT,GGTCCCCCAAGGCCG,CATGGCCGGCTCCAC,TGGGGGACCGACGGA,CTTGGAGGACCGACG,TCGGTCCCCCAAGGC,CCAAGGCCGGCTCCA,GCTCTGGGGGACCGA", "CCAAGGCCGGCTCCA,GGTCCCCCAAGGCCG,CAAGGCCGGCTCCAC,TGGGGGACCGACGGA,CTTGGAGGACCGACG,TCGGTCCCCCAAGGC,CCAAGGCCGGCTCCA,GCTCTGGGGGACCGA")]
    [TestCase("Forward", "ACAAAGCACCTCAAC,TACCAATCACTTCTA,TTGCGTTACCAACCA,GTGCTTGGGAACGCA,ACGTTACCAAGCACC,TGCTTGTTACCGCAA,CACCTCTACGGGTAT,CGTTACCAAGCACCT", "ACCAAGCACCTCAAC,TACCAATCACTTCTA,TTGCGTTACCAAGCA,GTGCTTGGGAACGCA,ACGTTACCAAGCACC,TGCTTGTTACCGCAA,CACCTCTACGGGTAT,CGTTACCAAGCACCT")]
    public void ErrorCorrectReads_MatchesMusketTranscriptionOracle(string mode, string input, string expected)
    {
        string[] reads = input.Split(',');

        IReadOnlyList<string> result = mode switch
        {
            "Musket" => SequenceAssembler.ErrorCorrectReads(reads, 5, 2, SequenceAssembler.ErrorCorrectionOptions.Musket),
            "Forward" => SequenceAssembler.ErrorCorrectReads(reads, 5, 2, SequenceAssembler.ErrorCorrectionOptions.TwoSidedForwardStrand),
            _ => SequenceAssembler.ErrorCorrectReads(reads, 5, 2),
        };

        Assert.That(result, Is.EqualTo(expected.Split(',')));
    }

    // F34 — validation: a cut-off below 1 makes every (even absent) k-mer "trusted"; null read elements.
    [TestCase(0)]
    [TestCase(-1)]
    public void ErrorCorrectReads_MinKmerFrequencyBelowOne_Throws(int cutoff)
    {
        Assert.That(() => SequenceAssembler.ErrorCorrectReads(new[] { "ACGT" }, 2, cutoff),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    [Test]
    public void ErrorCorrectReads_NullReadElement_ThrowsArgumentException()
    {
        Assert.That(() => SequenceAssembler.ErrorCorrectReads(new[] { "ACGT", null! }, 2, 2),
            NUnit.Framework.Throws.ArgumentException.With.Property("ParamName").EqualTo("reads"));
    }

    [Test]
    public void ErrorCorrectReads_InvalidOptions_Throw()
    {
        var reads = new[] { "ACGT" };
        Assert.Multiple(() =>
        {
            Assert.That(() => SequenceAssembler.ErrorCorrectReads(reads, 2, 2, null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => SequenceAssembler.ErrorCorrectReads(reads, 2, 2,
                    SequenceAssembler.ErrorCorrectionOptions.Musket with { MaxIterations = 0 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => SequenceAssembler.ErrorCorrectReads(reads, 2, 2,
                    SequenceAssembler.ErrorCorrectionOptions.Musket with { MaxErrorsPerKmer = -1 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // Musket presets and an explicit default instance agree (Musket 1.1 defaults: both strands,
    // -maxiter 2, -maxerr 4).
    [Test]
    public void ErrorCorrectionOptions_Presets_MatchMusketDefaults()
    {
        var musket = SequenceAssembler.ErrorCorrectionOptions.Musket;
        var twoSided = SequenceAssembler.ErrorCorrectionOptions.TwoSidedForwardStrand;
        Assert.Multiple(() =>
        {
            Assert.That(musket, Is.EqualTo(new SequenceAssembler.ErrorCorrectionOptions()));
            Assert.That((musket.CountReverseComplement, musket.MaxIterations, musket.MaxErrorsPerKmer),
                Is.EqualTo((true, 2, 4)));
            Assert.That((twoSided.CountReverseComplement, twoSided.MaxIterations, twoSided.MaxErrorsPerKmer),
                Is.EqualTo((false, int.MaxValue, 0)));
            var canonicalTwoSided = SequenceAssembler.ErrorCorrectionOptions.TwoSided;
            Assert.That((canonicalTwoSided.CountReverseComplement, canonicalTwoSided.MaxIterations, canonicalTwoSided.MaxErrorsPerKmer),
                Is.EqualTo((true, int.MaxValue, 0)), "3-argument overload preset: canonical spectrum, stage 1 only");
        });
    }

    // A9-2 — default parameters are the sourced Musket ones: k = 21 (musket-1.1 option.cpp when -k is absent) and
    //        cut-off 3 (smallest value of the main.cpp valley estimator: valleyIndex starts at 2, then ++); library
    //        and MCP wrapper use the same values. Reads of 25 bases with 3 identical copies (every 21-mer has
    //        multiplicity 3 >= 3) plus one copy with a substitution at index 12 (corrected back); with k = 15 /
    //        cut-off 2 (old library defaults) the result is the same here, so the defaults are pinned by equality
    //        with the explicit call and by a case where cut-off 2 vs 3 differ.
    [Test]
    public void ErrorCorrectReads_DefaultParameters_AreMusketK21AndCutoff3()
    {
        const string genome = "ACGTTGCAAGGCTTACCGATGCATG"; // 25 bases
        string error = genome.Remove(12, 1).Insert(12, genome[12] == 'A' ? "C" : "A");
        var threeCopies = new[] { genome, genome, genome, error };
        var twoCopies = new[] { genome, genome, error };

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.ErrorCorrectReads(threeCopies),
                Is.EqualTo(SequenceAssembler.ErrorCorrectReads(threeCopies, 21, 3,
                    SequenceAssembler.ErrorCorrectionOptions.TwoSided)));
            Assert.That(SequenceAssembler.ErrorCorrectReads(threeCopies)[3], Is.EqualTo(genome),
                "21-mers seen 3x are trusted; the unique alternative at index 12 restores the genome");
            // Two copies: 21-mers have multiplicity 2 < default cut-off 3 -> nothing trusted -> read unchanged.
            Assert.That(SequenceAssembler.ErrorCorrectReads(twoCopies),
                Is.EqualTo(twoCopies), "cut-off 3: multiplicity-2 k-mers are untrusted");
        });
    }

    // A2-5 — Musket 2-bit encoding (kmer.cpp Kmer::encode: every symbol but A/C/G/T/a/c/g/t is 0 = A): the spectrum and
    // every k-mer lookup see non-ACGT symbols (N, IUPAC codes, U, ...) as A, and a base no correction assigns keeps
    // its (upper-cased) input symbol. Oracle: tools/review/oracles/B10/wp4_oracle.py musket_fold_correct (Musket
    // stages 1-3 transcription + A-encoding); 375 random cases per seed x 7 seeds (N/IUPAC/U/lower-case noise,
    // k 5..9, cut-off 2..3, Musket options) equal the C# output.
    private const string NonAcgtGenome = "ACGTTGCAAGCTTGACCGTAGGATC"; // idx 7 = 'A', idx 15 = 'C'

    private static string[] SixCopiesPlus(params string[] extra) =>
        Enumerable.Repeat(NonAcgtGenome, 6).Concat(extra).ToArray();

    [Test]
    public void ErrorCorrectReads_NAtABaseThatIsA_IsSeenAsA_AndKeepsItsSymbol()
    {
        // N at index 7 where the genome has A: the folded read equals the genome, every k-mer is trusted, nothing is
        // assigned, so the N survives (a bare "N is its own symbol" spectrum would find untrusted k-mers).
        string withN = NonAcgtGenome[..7] + "N" + NonAcgtGenome[8..];

        var corrected = SequenceAssembler.ErrorCorrectReads(SixCopiesPlus(withN), 7, 3,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.That(corrected[6], Is.EqualTo(withN));
        Assert.That(corrected.Take(6), Is.All.EqualTo(NonAcgtGenome));
    }

    [Test]
    public void ErrorCorrectReads_NAtABaseThatIsNotA_IsCorrectedToTheTrustedBase()
    {
        // N at index 15 where the genome has C: the folded read has A there (an ordinary substitution error) and the
        // unique trusted alternative is C.
        string withN = NonAcgtGenome[..15] + "N" + NonAcgtGenome[16..];

        var corrected = SequenceAssembler.ErrorCorrectReads(SixCopiesPlus(withN), 7, 3,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.That(corrected[6], Is.EqualTo(NonAcgtGenome));
    }

    [Test]
    public void ErrorCorrectReads_IupacCodeUracilAndLowerCase_AreFoldedUpperCasedAndCorrected()
    {
        // r at 5 (genome G -> corrected), lower-case prefix (upper-cased), U at 12 (genome T -> corrected), and an IUPAC K at 7
        // (genome A = the folded base -> trusted, keeps its symbol K).
        string iupac = NonAcgtGenome[..5] + "r" + NonAcgtGenome[6..];
        string lower = NonAcgtGenome[..10].ToLowerInvariant() + NonAcgtGenome[10..];
        string uracil = NonAcgtGenome[..12] + "U" + NonAcgtGenome[13..];
        string keepK = NonAcgtGenome[..7] + "K" + NonAcgtGenome[8..];

        var corrected = SequenceAssembler.ErrorCorrectReads(SixCopiesPlus(iupac, lower, uracil, keepK), 7, 3,
            SequenceAssembler.ErrorCorrectionOptions.Musket);

        Assert.Multiple(() =>
        {
            Assert.That(corrected[6], Is.EqualTo(NonAcgtGenome), "r (folded to A) is replaced by the trusted G");
            Assert.That(corrected[7], Is.EqualTo(NonAcgtGenome), "lower case is upper-cased");
            Assert.That(corrected[8], Is.EqualTo(NonAcgtGenome), "U (folded to A) is replaced by the trusted T");
            Assert.That(corrected[9], Is.EqualTo(NonAcgtGenome[..7] + "K" + NonAcgtGenome[8..]), "K folds to the trusted A: untouched");
        });
    }

    #endregion
}
