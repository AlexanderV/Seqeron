namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class ISequenceTests
{
    #region IupacDnaSequence Tests

    [Test]
    public void IupacDnaSequence_StandardBases_Valid()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq.Length, Is.EqualTo(4));
        Assert.That(seq.Sequence, Is.EqualTo("ACGT"));
        Assert.That(seq.IsValid(), Is.True);
        Assert.That(seq.Type, Is.EqualTo(SequenceType.IupacDna));
    }

    [Test]
    public void IupacDnaSequence_AmbiguityCodes_Valid()
    {
        var seq = new IupacDnaSequence("ACGTNRYSWKMBDHV");

        Assert.That(seq.IsValid(), Is.True);
    }

    [Test]
    public void IupacDnaSequence_GetComplement_CorrectlyComplements()
    {
        var seq = new IupacDnaSequence("ACGTRYWSKMBDHVN");
        var comp = seq.GetComplement() as IupacDnaSequence;

        Assert.That(comp, Is.Not.Null);
        // Correct IUPAC complements:
        // A->T, C->G, G->C, T->A, R->Y, Y->R, W->W, S->S, K->M, M->K, B->V, D->H, H->D, V->B, N->N
        Assert.That(comp!.Sequence, Is.EqualTo("TGCAYRWSMKVHDBN"));
    }

    [Test]
    [Description("SEQ-COMP-001: IupacDnaSequence delegates to canonical GetComplementBase — Biopython reverse_complement('ACGTRYWSKMBDHVN') = 'NBDHVKMSWRYACGT'; complement('AC-.G') = 'TG-.C'")]
    public void IupacDnaSequence_Complement_MatchesBiopythonIncludingGaps()
    {
        var rc = new IupacDnaSequence("ACGTRYWSKMBDHVN").GetReverseComplement() as IupacDnaSequence;
        var gaps = new IupacDnaSequence("AC-.G").GetComplement() as IupacDnaSequence;
        var nonIupac = new IupacDnaSequence("AXG").GetComplement() as IupacDnaSequence;

        Assert.Multiple(() =>
        {
            Assert.That(rc!.Sequence, Is.EqualTo("NBDHVKMSWRYACGT"));
            Assert.That(gaps!.Sequence, Is.EqualTo("TG-.C"), "gaps are self-complementary");
            Assert.That(nonIupac!.Sequence, Is.EqualTo("TNC"), "non-alphabet symbol becomes N");
        });
    }

    [Test]
    public void IupacDnaSequence_GetReverseComplement_Works()
    {
        var seq = new IupacDnaSequence("ACGT");
        var rc = seq.GetReverseComplement() as IupacDnaSequence;

        Assert.That(rc, Is.Not.Null);
        Assert.That(rc!.Sequence, Is.EqualTo("ACGT")); // Palindrome
    }

    [Test]
    public void IupacDnaSequence_ExpandCode_ExpandsCorrectly()
    {
        Assert.That(IupacDnaSequence.ExpandCode('A'), Is.EquivalentTo(new[] { 'A' }));
        Assert.That(IupacDnaSequence.ExpandCode('N'), Is.EquivalentTo(new[] { 'A', 'C', 'G', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('R'), Is.EquivalentTo(new[] { 'A', 'G' }));
        Assert.That(IupacDnaSequence.ExpandCode('Y'), Is.EquivalentTo(new[] { 'C', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('W'), Is.EquivalentTo(new[] { 'A', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('S'), Is.EquivalentTo(new[] { 'G', 'C' }));
    }

    [Test]
    public void IupacDnaSequence_GetIupacCode_EncodesCorrectly()
    {
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A' }), Is.EqualTo('A'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A', 'G' }), Is.EqualTo('R'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'C', 'T' }), Is.EqualTo('Y'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A', 'C', 'G', 'T' }), Is.EqualTo('N'));
    }

    [Test]
    public void IupacDnaSequence_CodesMatch_MatchesCorrectly()
    {
        Assert.That(IupacDnaSequence.CodesMatch('A', 'A'), Is.True);
        Assert.That(IupacDnaSequence.CodesMatch('A', 'R'), Is.True);  // A matches R (A/G)
        Assert.That(IupacDnaSequence.CodesMatch('A', 'Y'), Is.False); // A doesn't match Y (C/T)
        Assert.That(IupacDnaSequence.CodesMatch('N', 'A'), Is.True);  // N matches anything
        Assert.That(IupacDnaSequence.CodesMatch('R', 'Y'), Is.False); // No overlap
        Assert.That(IupacDnaSequence.CodesMatch('W', 'M'), Is.True);  // A in common
    }

    [Test]
    public void IupacDnaSequence_MatchesAt_WithWildcards()
    {
        var seq = new IupacDnaSequence("ACGTACGT");

        Assert.That(seq.MatchesAt("ACGT", 0), Is.True);
        Assert.That(seq.MatchesAt("NNNN", 0), Is.True);  // N matches any
        Assert.That(seq.MatchesAt("RCGT", 0), Is.True);  // R=A/G, A matches
        Assert.That(seq.MatchesAt("YCGT", 0), Is.False); // Y=C/T, A doesn't match
    }

    [Test]
    public void IupacDnaSequence_FindPattern_WithWildcards()
    {
        var seq = new IupacDnaSequence("AAAAACGTAAAA");

        var positions = seq.FindPattern("RCGT").ToList();

        Assert.That(positions, Contains.Item(4));
    }

    [Test]
    public void IupacDnaSequence_GetAmbiguityLevel_CalculatesCorrectly()
    {
        var unambiguous = new IupacDnaSequence("ACGT");
        var halfAmbiguous = new IupacDnaSequence("ACNN");
        var allAmbiguous = new IupacDnaSequence("NNNN");

        Assert.That(unambiguous.GetAmbiguityLevel(), Is.EqualTo(1.0));
        Assert.That(halfAmbiguous.GetAmbiguityLevel(), Is.EqualTo(0.5));
        Assert.That(allAmbiguous.GetAmbiguityLevel(), Is.EqualTo(0.0));
    }

    // scikit-bio 0.7.4: DNA(s, lowercase=True).definites().mean() / .degenerates().mean()
    [TestCase("ACNR-G.T", 0.5, 0.25)]
    [TestCase("ACGT", 1.0, 0.0)]
    [TestCase("NNNN", 0.0, 1.0)]
    [TestCase("acgt", 1.0, 0.0)]
    [TestCase("ACGTN", 0.8, 0.2)]
    [TestCase("RYSWKMBDHVN", 0.0, 1.0)]
    [TestCase("----", 0.0, 0.0)]
    [TestCase("A-", 0.5, 0.0)]
    [TestCase("ACGTNNNN", 0.5, 0.5)]
    public void IupacDnaSequence_DefiniteAndDegenerateFractions_MatchScikitBio(string seq, double definite, double degenerate)
    {
        var s = new IupacDnaSequence(seq);
        Assert.Multiple(() =>
        {
            Assert.That(s.GetAmbiguityLevel(), Is.EqualTo(definite).Within(1e-15), "definites().mean()");
            Assert.That(s.GetDegenerateFraction(), Is.EqualTo(degenerate).Within(1e-15), "degenerates().mean()");
        });
    }

    [Test]
    public void IupacDnaSequence_DefiniteAndDegenerateFractions_EmptyAndU_DocumentedConventions()
    {
        // scikit-bio gives NaN for the empty mean; Seqeron documents 1.0 / 0.0 (no degenerate position).
        var empty = new IupacDnaSequence("");
        // scikit-bio DNA rejects U; this container tolerates it and counts it as neither definite nor degenerate.
        var withU = new IupacDnaSequence("ACGU");
        Assert.Multiple(() =>
        {
            Assert.That(empty.GetAmbiguityLevel(), Is.EqualTo(1.0));
            Assert.That(empty.GetDegenerateFraction(), Is.EqualTo(0.0));
            Assert.That(withU.GetAmbiguityLevel(), Is.EqualTo(0.75));
            Assert.That(withU.GetDegenerateFraction(), Is.EqualTo(0.0));
        });
    }

    [Test]
    public void IupacDnaSequence_ExpandAll_GeneratesAllPossibilities()
    {
        var seq = new IupacDnaSequence("AN");

        var expanded = seq.ExpandAll().ToList();

        Assert.That(expanded, Has.Count.EqualTo(4));
        Assert.That(expanded, Contains.Item("AA"));
        Assert.That(expanded, Contains.Item("AC"));
        Assert.That(expanded, Contains.Item("AG"));
        Assert.That(expanded, Contains.Item("AT"));
    }

    [Test]
    public void IupacDnaSequence_ExpandAll_LimitsResults()
    {
        var seq = new IupacDnaSequence("NNNNNN"); // 4^6 = 4096 possibilities

        var expanded = seq.ExpandAll(maxResults: 100).ToList();

        Assert.That(expanded.Count, Is.LessThanOrEqualTo(100));
    }

    [Test]
    public void IupacDnaSequence_Subsequence_ReturnsCorrectType()
    {
        var seq = new IupacDnaSequence("ACGTNRYSWKM");
        var sub = seq.Subsequence(2, 5);

        Assert.That(sub, Is.TypeOf<IupacDnaSequence>());
        Assert.That(sub.Sequence, Is.EqualTo("GTNRY"));
    }

    #endregion

    #region QualitySequence Tests

    [Test]
    public void QualitySequence_Constructor_SetsSequenceAndQuality()
    {
        var qual = new byte[] { 30, 30, 30, 30 };
        var seq = new QualitySequence("ACGT", qual);

        Assert.That(seq.Sequence, Is.EqualTo("ACGT"));
        Assert.That(seq.Qualities, Is.EquivalentTo(qual));
        Assert.That(seq.Type, Is.EqualTo(SequenceType.Quality));
    }

    [Test]
    public void QualitySequence_FromQualityString_ParsesCorrectly()
    {
        // Phred+33: 'I' = 40, '5' = 20
        var seq = new QualitySequence("ACGT", "II55", phredOffset: 33);

        Assert.That(seq.GetQuality(0), Is.EqualTo(40));
        Assert.That(seq.GetQuality(1), Is.EqualTo(40));
        Assert.That(seq.GetQuality(2), Is.EqualTo(20));
        Assert.That(seq.GetQuality(3), Is.EqualTo(20));
    }

    [Test]
    public void QualitySequence_MismatchedLength_ThrowsException()
    {
        Assert.Throws<ArgumentException>(() =>
            new QualitySequence("ACGT", new byte[] { 30, 30 }));
    }

    [Test]
    public void QualitySequence_MeanQuality_CalculatesCorrectly()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });

        Assert.That(seq.MeanQuality, Is.EqualTo(25.0));
    }

    [Test]
    public void QualitySequence_MeanQuality_Empty_ThrowsDocumentedException()
    {
        // Python statistics.mean([]) raises StatisticsError ("mean requires at least one data point").
        var empty = new QualitySequence("", Array.Empty<byte>());
        var ex = Assert.Throws<InvalidOperationException>(() => _ = empty.MeanQuality);
        Assert.That(ex!.Message, Does.Contain("empty"));
        // statistics.mean([0, 93, 41]) = 44.666666666666664
        Assert.That(new QualitySequence("ACG", new byte[] { 0, 93, 41 }).MeanQuality,
            Is.EqualTo(44.666666666666664).Within(1e-12));
    }

    [Test]
    public void QualitySequence_GetQualityString_EncodesCorrectly()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 0, 10, 30, 40 });

        var qualStr = seq.GetQualityString(33);

        Assert.That(qualStr[0], Is.EqualTo('!')); // 0 + 33 = '!'
        Assert.That(qualStr[3], Is.EqualTo('I')); // 40 + 33 = 'I'
    }

    [Test]
    public void QualitySequence_TrimByQuality_TrimsLowQualityEnds()
    {
        var seq = new QualitySequence("AACGTAA", new byte[] { 5, 5, 30, 30, 30, 5, 5 });

        var trimmed = seq.TrimByQuality(minQuality: 20);

        Assert.That(trimmed.Sequence, Is.EqualTo("CGT"));
        Assert.That(trimmed.Length, Is.EqualTo(3));
    }

    [Test]
    public void QualitySequence_TrimByQuality_AllLowQuality_ReturnsEmpty()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 5, 5, 5, 5 });

        var trimmed = seq.TrimByQuality(minQuality: 20);

        Assert.That(trimmed.Length, Is.EqualTo(0));
    }

    [Test]
    public void QualitySequence_MaskLowQuality_MasksWithN()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 30, 5, 30, 5 });

        var masked = seq.MaskLowQuality(minQuality: 20);

        Assert.That(masked.Sequence, Is.EqualTo("ANGN"));
        Assert.That(masked.Sequence[0], Is.EqualTo('A'));
        Assert.That(masked.Sequence[1], Is.EqualTo('N'));
        Assert.That(masked.Sequence[2], Is.EqualTo('G'));
        Assert.That(masked.Sequence[3], Is.EqualTo('N'));
    }

    [Test]
    public void QualitySequence_PhredToErrorProbability_ConvertsCorrectly()
    {
        Assert.That(QualitySequence.PhredToErrorProbability(10), Is.EqualTo(0.1).Within(0.001));
        Assert.That(QualitySequence.PhredToErrorProbability(20), Is.EqualTo(0.01).Within(0.0001));
        Assert.That(QualitySequence.PhredToErrorProbability(30), Is.EqualTo(0.001).Within(0.00001));
    }

    [Test]
    public void QualitySequence_ErrorProbabilityToPhred_ConvertsCorrectly()
    {
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.1), Is.EqualTo(10));
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.01), Is.EqualTo(20));
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.001), Is.EqualTo(30));
    }

    [Test]
    public void QualitySequence_ExpectedErrors_CalculatesCorrectly()
    {
        // Q30 = 0.001 error probability per base
        var seq = new QualitySequence("ACGT", new byte[] { 30, 30, 30, 30 });

        var expected = seq.ExpectedErrors();

        Assert.That(expected, Is.EqualTo(0.004).Within(0.0001));
    }

    [Test]
    public void QualitySequence_GetComplement_PreservesQuality()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });
        var comp = seq.GetComplement() as QualitySequence;

        Assert.That(comp, Is.Not.Null);
        Assert.That(comp!.Sequence, Is.EqualTo("TGCA"));
        Assert.That(comp.Qualities, Is.EquivalentTo(new byte[] { 10, 20, 30, 40 }));
    }

    [Test]
    [Description("SEQ-COMP-001: self-complementary IUPAC codes S/W/N are not 'unknown' — Biopython complement('ASWRN') = 'TSWYN', reverse_complement = 'NYWST'")]
    public void QualitySequence_Complement_SelfComplementaryIupacCodesPreserved()
    {
        var seq = new QualitySequence("aswrn", new byte[] { 10, 20, 30, 40, 50 });
        var comp = seq.GetComplement() as QualitySequence;
        var rc = seq.GetReverseComplement() as QualitySequence;
        var withGap = new QualitySequence("A-G", new byte[] { 1, 2, 3 }).GetComplement() as QualitySequence;

        Assert.Multiple(() =>
        {
            Assert.That(comp!.Sequence, Is.EqualTo("TSWYN"));
            Assert.That(rc!.Sequence, Is.EqualTo("NYWST"));
            Assert.That(rc.Qualities, Is.EqualTo(new byte[] { 50, 40, 30, 20, 10 }));
            Assert.That(withGap!.Sequence, Is.EqualTo("TNC"), "non-IUPAC characters become N");
        });
    }

    [Test]
    public void QualitySequence_GetReverseComplement_ReversesQuality()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });
        var rc = seq.GetReverseComplement() as QualitySequence;

        Assert.That(rc, Is.Not.Null);
        Assert.That(rc!.Sequence, Is.EqualTo("ACGT"));
        Assert.That(rc.Qualities, Is.EquivalentTo(new byte[] { 40, 30, 20, 10 }));
    }

    // B01 finisher dedup sweep (uncovered public method): GetReverse. Reference Biopython 1.88:
    // Seq("ACNR-G.T")[::-1] = "T.G-RNCA"; SeqRecord(Seq("ACGTN"), phred_quality=[10,20,30,40,2])[::-1]
    // → "NTGCA", [2,40,30,20,10] (letter annotations are reversed with the sequence).
    [Test]
    public void IupacAndQualitySequence_GetReverse_MatchesBiopythonSliceReversal()
    {
        var iupac = new IupacDnaSequence("ACNR-G.T").GetReverse();
        var qual = new QualitySequence("ACGTN", new byte[] { 10, 20, 30, 40, 2 }).GetReverse() as QualitySequence;

        Assert.Multiple(() =>
        {
            Assert.That(iupac, Is.InstanceOf<IupacDnaSequence>());
            Assert.That(iupac.Sequence, Is.EqualTo("T.G-RNCA"));
            Assert.That(qual, Is.Not.Null);
            Assert.That(qual!.Sequence, Is.EqualTo("NTGCA"));
            Assert.That(qual.Qualities, Is.EqualTo(new byte[] { 2, 40, 30, 20, 10 }));
        });
    }

    [Test]
    public void QualitySequence_Subsequence_PreservesQuality()
    {
        var seq = new QualitySequence("AACGTAA", new byte[] { 10, 20, 30, 40, 50, 60, 70 });
        var sub = seq.Subsequence(2, 3) as QualitySequence;

        Assert.That(sub, Is.Not.Null);
        Assert.That(sub!.Sequence, Is.EqualTo("CGT"));
        Assert.That(sub.Qualities, Is.EquivalentTo(new byte[] { 30, 40, 50 }));
    }

    #endregion

    #region SequenceBase Tests

    [Test]
    public void SequenceBase_Indexer_ReturnsCorrectChar()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq[0], Is.EqualTo('A'));
        Assert.That(seq[2], Is.EqualTo('G'));
    }

    [Test]
    public void SequenceBase_ToString_ReturnsSequence()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq.ToString(), Is.EqualTo("ACGT"));
    }

    [Test]
    public void SequenceBase_Equals_ComparesSequences()
    {
        var seq1 = new IupacDnaSequence("ACGT");
        var seq2 = new IupacDnaSequence("ACGT");
        var seq3 = new IupacDnaSequence("TGCA");

        Assert.That(seq1.Equals(seq2), Is.True);
        Assert.That(seq1.Equals(seq3), Is.False);
    }

    [Test]
    public void SequenceBase_GetHashCode_SameForEqualSequences()
    {
        var seq1 = new IupacDnaSequence("ACGT");
        var seq2 = new IupacDnaSequence("ACGT");

        Assert.That(seq1.GetHashCode(), Is.EqualTo(seq2.GetHashCode()));
    }

    #endregion

    #region B01-SWEEP — sourced fixes (Biopython 1.88 / scikit-bio 0.7.4 references)

    [TestCase('r', 'a', true)]
    [TestCase('n', 't', true)]
    [TestCase('k', 'Y', true)]
    [TestCase('s', 'w', false)]
    [TestCase('b', 'a', false)]
    [TestCase('d', 'c', false)]
    [TestCase('m', 'K', false)]
    [TestCase('h', 'g', false)]
    [Description("B01-SWEEP: CodesMatch = intersection of NC-IUB base sets (Biopython IUPACData.ambiguous_dna_values), ASCII case-insensitive like ExpandCode")]
    public void IupacDnaSequence_CodesMatch_CaseInsensitive_MatchesBiopythonSetIntersection(char c1, char c2, bool expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacDnaSequence.CodesMatch(c1, c2), Is.EqualTo(expected));
            Assert.That(IupacDnaSequence.CodesMatch(c2, c1), Is.EqualTo(expected), "symmetric");
            Assert.That(IupacDnaSequence.CodesMatch(char.ToUpperInvariant(c1), char.ToUpperInvariant(c2)), Is.EqualTo(expected));
        });
    }

    [Test]
    [Description("A1-6/F19: ExpandCode/MatchesAt fold ASCII only — U+017F 'ſ' is not S (Biopython ambiguous_dna_values has no 'ſ' key, nt_search raises KeyError; scikit-bio DNA('ſ') raises)")]
    public void IupacDnaSequence_ExpandCodeAndMatchesAt_NonAsciiLettersNotFolded()
    {
        var seq = new IupacDnaSequence("ACGTSN\u017f");
        Assert.Multiple(() =>
        {
            // Unknown-character contract: { 'N' } (as for 'X', '1', '-').
            Assert.That(IupacDnaSequence.ExpandCode('\u017f'), Is.EqualTo(new[] { 'N' }));
            Assert.That(IupacDnaSequence.ExpandCode('\u0131'), Is.EqualTo(new[] { 'N' }));
            Assert.That(IupacDnaSequence.ExpandCode('\u212a'), Is.EqualTo(new[] { 'N' }));
            Assert.That(IupacDnaSequence.ExpandCode('X'), Is.EqualTo(new[] { 'N' }));
            Assert.That(IupacDnaSequence.ExpandCode('s'), Is.EquivalentTo(new[] { 'G', 'C' }));
            Assert.That(IupacDnaSequence.ExpandCode('k'), Is.EquivalentTo(new[] { 'G', 'T' }));

            // Pattern 'ſ' is not S: no match at G (2), C (1) or S (4); it matches only the literal 'ſ' (6),
            // a non-code symbol matching itself per the CodesMatch contract.
            Assert.That(seq.MatchesAt("\u017f", 1), Is.False);
            Assert.That(seq.MatchesAt("\u017f", 2), Is.False);
            Assert.That(seq.MatchesAt("\u017f", 4), Is.False);
            Assert.That(seq.FindPattern("\u017f").ToList(), Is.EqualTo(new[] { 6 }));
            Assert.That(seq.FindPattern("s").ToList(), Is.EqualTo(new[] { 1, 2, 4, 5 }));
            Assert.That(seq.FindPattern("\u212a").ToList(), Is.Empty);
            Assert.That(seq.MatchesAt("acgtsn", 0), Is.True);
        });
    }

    [TestCase("\u017f")]
    [TestCase("AC\u017f")]
    [Description("A1-6/F19: DnaSequence/RnaSequence ctor reports the rejected non-ASCII character verbatim (not the ToUpperInvariant fold 'S')")]
    public void DnaRnaSequence_Ctor_NonAsciiLetter_RejectedAndReportedVerbatim(string sequence)
    {
        Assert.Multiple(() =>
        {
            var dna = Assert.Throws<ArgumentException>(() => new DnaSequence(sequence));
            Assert.That(dna!.Message, Does.Contain("'\u017f'"));
            var rna = Assert.Throws<ArgumentException>(() => new RnaSequence(sequence));
            Assert.That(rna!.Message, Does.Contain("'\u017f'"));
            Assert.That(new DnaSequence("acgt").Sequence, Is.EqualTo("ACGT"));
            Assert.That(new RnaSequence("acgu").Sequence, Is.EqualTo("ACGU"));
        });
    }

    [TestCase("AU", 'W')]
    [TestCase("cu", 'Y')]
    [TestCase("gu", 'K')]
    [TestCase("ACU", 'H')]
    [TestCase("acgu", 'N')]
    [TestCase("ag", 'R')]
    [Description("B01-SWEEP: U is T (Biopython ambiguous_rna_values W='AU', Y='CU', K='GU', H='ACU'); consistent with ExpandCode('U') = T")]
    public void IupacDnaSequence_GetIupacCode_RnaUracilTreatedAsThymine(string bases, char expected)
    {
        Assert.That(IupacDnaSequence.GetIupacCode(bases), Is.EqualTo(expected));
    }

    [TestCase(0.2, 7)]
    [TestCase(0.0011, 30)]
    [TestCase(0.5, 3)]
    [TestCase(0.3, 5)]
    [TestCase(0.05, 13)]
    [TestCase(0.9, 0)]
    [TestCase(1.0, 0)]
    [TestCase(1e-9, 90)]
    [TestCase(1e-10, 93)]
    [TestCase(0.0, 93)]
    [Description("B01-SWEEP: Q = round(-10 log10 p), capped at Q93 (Cock et al. 2010; Biopython _get_sanger_quality_str rounds: 0.2 -> 7, 0.0011 -> 30, 1e-10 -> 93)")]
    public void QualitySequence_ErrorProbabilityToPhred_RoundsLikeBiopython(double p, int expected)
    {
        Assert.That((int)QualitySequence.ErrorProbabilityToPhred(p), Is.EqualTo(expected));
    }

    [TestCase(double.NaN)]
    [TestCase(-0.1)]
    [TestCase(1.1)]
    public void QualitySequence_ErrorProbabilityToPhred_OutsideProbabilityDomain_Throws(double p)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QualitySequence.ErrorProbabilityToPhred(p));
    }

    [Test]
    [Description("B01-SWEEP: Biopython FASTQ parsing — 'II!' -> [40,40,0], '~~~' -> [93,93,93]; fastq-illumina 'h@~' -> [40,0,62]")]
    public void QualitySequence_QualityStringCtor_DecodesLikeBiopython()
    {
        Assert.Multiple(() =>
        {
            Assert.That(new QualitySequence("ACG", "II!").Qualities, Is.EqualTo(new byte[] { 40, 40, 0 }));
            Assert.That(new QualitySequence("ACG", "~~~").Qualities, Is.EqualTo(new byte[] { 93, 93, 93 }));
            Assert.That(new QualitySequence("ACG", "h@~", phredOffset: 64).Qualities, Is.EqualTo(new byte[] { 40, 0, 62 }));
        });
    }

    [TestCase("ACG", "II", 33)]
    [TestCase("ACG", "IIII", 33)]
    [TestCase("ACG", "I I", 33)]
    [TestCase("ACG", "II\u007f", 33)]
    [TestCase("ACG", "h?~", 64)]
    [Description("B01-SWEEP: Biopython raises on quality/sequence length mismatch and on characters outside [offset, 126]")]
    public void QualitySequence_QualityStringCtor_InvalidQualityString_Throws(string seq, string qual, int offset)
    {
        qual = System.Text.RegularExpressions.Regex.Unescape(qual);
        Assert.Throws<ArgumentException>(() => new QualitySequence(seq, qual, offset));
    }

    [TestCase(32)]
    [TestCase(127)]
    public void QualitySequence_InvalidPhredOffset_Throws(int offset)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new QualitySequence("A", "I", offset));
        Assert.Throws<ArgumentOutOfRangeException>(() => new QualitySequence("A", new byte[] { 40 }).GetQualityString(offset));
    }

    [Test]
    [Description("B01-SWEEP: Biopython _get_sanger_quality_str([50,40,30,20,10,0]) = 'SI?5+!'; Q>93 capped at '~' (Sanger) and Q>62 at '~' (Illumina 1.3+)")]
    public void QualitySequence_GetQualityString_MatchesBiopythonAndCapsAtTilde()
    {
        Assert.Multiple(() =>
        {
            Assert.That(new QualitySequence("ACGTAN", new byte[] { 50, 40, 30, 20, 10, 0 }).GetQualityString(), Is.EqualTo("SI?5+!"));
            Assert.That(new QualitySequence("AC", new byte[] { 100, 93 }).GetQualityString(33), Is.EqualTo("~~"));
            Assert.That(new QualitySequence("AC", new byte[] { 70, 62 }).GetQualityString(64), Is.EqualTo("~~"));
        });
    }

    #endregion
}
