namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Tests for Translator protein translation.
/// Test Unit: TRANS-PROT-001
/// 
/// Evidence Sources:
///   - Wikipedia: Translation (biology), Reading frame, Open reading frame
///   - NCBI: The Genetic Codes (translation tables)
///   - Lodish H et al. (2007): Molecular Cell Biology
/// 
/// Key Invariants:
///   - Codons are read in triplets from 5' to 3'
///   - Frame parameter offsets reading start position
///   - Stop codons (UAA, UAG, UGA) terminate translation
///   - Six-frame translation covers both strands
/// </summary>
[TestFixture]
public class TranslatorTests
{
    #region Basic Translation

    [Test]
    public void Translate_SingleCodon_ReturnsSingleAminoAcid()
    {
        var dna = new DnaSequence("ATG");
        var protein = Translator.Translate(dna);
        Assert.That(protein.Sequence, Is.EqualTo("M"));
    }

    [Test]
    public void Translate_MultipleCodens_ReturnsProtein()
    {
        // ATG GCT TAA = M A *
        var dna = new DnaSequence("ATGGCTTAA");
        var protein = Translator.Translate(dna);
        Assert.That(protein.Sequence, Is.EqualTo("MA*"));
    }

    [Test]
    public void Translate_ToFirstStop_StopsAtStopCodon()
    {
        // ATG GCT TAA GCT = M A * A
        var dna = new DnaSequence("ATGGCTTAAGCT");
        var protein = Translator.Translate(dna, toFirstStop: true);
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    [Test]
    public void Translate_Frame1_ShiftsReading()
    {
        // A ATG GCT = skip A, then ATG GCT = M A
        var dna = new DnaSequence("AATGGCT");
        var protein = Translator.Translate(dna, frame: 1);
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    [Test]
    public void Translate_Frame2_ShiftsReading()
    {
        // AA ATG GCT = skip AA, then ATG GCT = M A
        var dna = new DnaSequence("AAATGGCT");
        var protein = Translator.Translate(dna, frame: 2);
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    [Test]
    public void Translate_InvalidFrame_ThrowsException()
    {
        var dna = new DnaSequence("ATGGCT");
        Assert.Throws<ArgumentOutOfRangeException>(() => Translator.Translate(dna, frame: 3));
    }

    [Test]
    public void Translate_EmptySequence_ReturnsEmpty()
    {
        var protein = Translator.Translate("");
        Assert.That(protein.Sequence, Is.Empty);
    }

    [Test]
    public void Translate_NullDna_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() => Translator.Translate((DnaSequence)null!));
    }

    [Test]
    public void Translate_NullRna_ThrowsException()
    {
        // Source: Implementation spec - null input handling
        Assert.Throws<ArgumentNullException>(() => Translator.Translate((RnaSequence)null!));
    }

    [Test]
    public void Translate_SequenceShorterThan3_ReturnsEmpty()
    {
        // Less than one complete codon — no amino acid produced
        // Source: Wikipedia - codons are read in triplets
        var protein = Translator.Translate("AT");
        Assert.That(protein.Sequence, Is.Empty);
    }

    #endregion

    #region RNA Translation

    [Test]
    public void Translate_Rna_Works()
    {
        var rna = new RnaSequence("AUGGCUUAA");
        var protein = Translator.Translate(rna);
        Assert.That(protein.Sequence, Is.EqualTo("MA*"));
    }

    [Test]
    public void Translate_RnaToFirstStop_Works()
    {
        var rna = new RnaSequence("AUGGCUUAAGCU");
        var protein = Translator.Translate(rna, toFirstStop: true);
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    #endregion

    #region String Translation

    [Test]
    public void Translate_DnaString_ConvertsTToU()
    {
        // DNA string with T must produce identical result to RNA with U
        // Source: Wikipedia - DNA T corresponds to RNA U in the genetic code
        var dnaResult = Translator.Translate("ATGGCTTAA");
        var rnaResult = Translator.Translate(new RnaSequence("AUGGCUUAA"));
        Assert.That(dnaResult.Sequence, Is.EqualTo("MA*"));
        Assert.That(dnaResult.Sequence, Is.EqualTo(rnaResult.Sequence));
    }

    [Test]
    public void Translate_LowercaseString_Works()
    {
        var protein = Translator.Translate("atggct");
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    #endregion

    #region Alternative Genetic Codes

    [Test]
    public void Translate_VertebrateMitochondrial_UsesDifferentCode()
    {
        // AGA is Arg in standard, but Stop in vertebrate mitochondrial
        var dna = new DnaSequence("ATGAGA");

        var standardProtein = Translator.Translate(dna, GeneticCode.Standard);
        Assert.That(standardProtein.Sequence, Is.EqualTo("MR"));

        var mitoProtein = Translator.Translate(dna, GeneticCode.VertebrateMitochondrial);
        Assert.That(mitoProtein.Sequence, Is.EqualTo("M*"));
    }

    [Test]
    public void Translate_YeastMitochondrial_CUU_IsThreonine()
    {
        // CUU is Leu in standard, but Thr in yeast mitochondrial
        var dna = new DnaSequence("ATGCTT");

        var standardProtein = Translator.Translate(dna, GeneticCode.Standard);
        Assert.That(standardProtein.Sequence, Is.EqualTo("ML"));

        var yeastProtein = Translator.Translate(dna, GeneticCode.YeastMitochondrial);
        Assert.That(yeastProtein.Sequence, Is.EqualTo("MT"));
    }

    #endregion

    #region Six Frame Translation

    [Test]
    public void TranslateSixFrames_ReturnsAllSixFrames()
    {
        var dna = new DnaSequence("ATGGCTAAA");
        var frames = Translator.TranslateSixFrames(dna);

        Assert.That(frames.Count, Is.EqualTo(6));
        Assert.That(frames.ContainsKey(1), Is.True);
        Assert.That(frames.ContainsKey(2), Is.True);
        Assert.That(frames.ContainsKey(3), Is.True);
        Assert.That(frames.ContainsKey(-1), Is.True);
        Assert.That(frames.ContainsKey(-2), Is.True);
        Assert.That(frames.ContainsKey(-3), Is.True);
    }

    [Test]
    public void TranslateSixFrames_Frame1_MatchesDirect()
    {
        var dna = new DnaSequence("ATGGCTAAA");
        var frames = Translator.TranslateSixFrames(dna);
        var direct = Translator.Translate(dna, frame: 0);

        Assert.That(frames[1].Sequence, Is.EqualTo(direct.Sequence));
    }

    [Test]
    public void TranslateSixFrames_NegativeFrames_UseReverseComplement()
    {
        // All 3 negative frames must match direct translation of reverse complement
        // Source: Wikipedia Reading frame - 6 frames from double-stranded DNA
        var dna = new DnaSequence("ATGGCTAAA");
        var revComp = dna.ReverseComplement();
        var frames = Translator.TranslateSixFrames(dna);

        for (int f = 0; f < 3; f++)
        {
            var expected = Translator.Translate(revComp, frame: f);
            Assert.That(frames[-(f + 1)].Sequence, Is.EqualTo(expected.Sequence),
                $"Frame -{f + 1} should match reverse complement frame {f}");
        }
    }

    [Test]
    public void TranslateSixFrames_NullInput_ThrowsException()
    {
        // Source: Implementation spec - null input handling
        Assert.Throws<ArgumentNullException>(() => Translator.TranslateSixFrames(null!));
    }

    [Test]
    public void TranslateSixFrames_EmptySequence_ReturnsEmptyFrames()
    {
        var dna = new DnaSequence("");
        var frames = Translator.TranslateSixFrames(dna);

        Assert.That(frames.Count, Is.EqualTo(6));
        Assert.That(frames.Values.All(p => p.Sequence == ""), Is.True);
    }

    #endregion

    #region ORF Finding (Smoke Tests — canonical tests in GenomeAnnotator_ORF_Tests.cs)

    [Test]
    public void FindOrfs_SimpleOrf_FindsIt()
    {
        // ATG GCT TTC TAA = M A F * → protein "MAF" (3 amino acids)
        // Source: Wikipedia ORF - start codon to stop codon
        var dna = new DnaSequence("ATGGCTTTCTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1));
        Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MAF"));
    }

    [Test]
    public void FindOrfs_NoStartCodon_ReturnsEmpty()
    {
        // No ATG/TTG/CTG (start codons) in any reading frame → no ORF
        // Source: Wikipedia ORF - requires start codon
        var dna = new DnaSequence("GCCGCCGCCTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty);
    }

    [Test]
    public void FindOrfs_RespectMinLength_FindsSmallOrfs()
    {
        // ATG GCT TAA = 2 amino acids (M A), meets minLength=2
        var dna = new DnaSequence("ATGGCTTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 2, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1));
        Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MA"));
    }

    [Test]
    public void FindOrfs_ShortOrf_FilteredByMinLength()
    {
        // ATG GCT TAA = protein "MA" (2 aa) < minLength 5 → filtered out
        // Source: Wikipedia ORF - short ORFs excluded by length threshold
        var dna = new DnaSequence("ATGGCTTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 5, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty);
    }

    [Test]
    public void FindOrfs_ForwardOnly_DoesNotSearchReverseStrand()
    {
        // TTAGCCGCCCAT has no start codons on forward strand (any frame)
        // Reverse complement = ATGGGCGGCTAA which has ATG...TAA ORF
        // With searchBothStrands=false, reverse-strand ORF must not appear
        var dna = new DnaSequence("TTAGCCGCCCAT");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Is.Empty);
    }

    [Test]
    public void FindOrfs_BothStrands_SearchesReverseComplement()
    {
        // Forward: TTA GCC GCC CAT → no ATG → no ORFs
        // Reverse complement: ATG GGC GGC TAA → M G G * → ORF found
        // Source: Wikipedia ORF - six-frame search covers both strands
        var dna = new DnaSequence("TTAGCCGCCCAT");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: true).ToList();

        Assert.That(orfs, Has.Count.EqualTo(1));
        Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MGG"));
        Assert.That(orfs[0].Frame, Is.Negative, "ORF should be on reverse strand");
    }

    [Test]
    public void FindOrfs_OrfResult_HasCorrectPositions()
    {
        // ATG GCT TAA at positions 0-8 in frame 0 (reported as frame 1)
        // Source: Wikipedia ORF - ORF spans from start to stop codon
        var dna = new DnaSequence("ATGGCTTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(orfs, Has.Count.EqualTo(1));
            Assert.That(orfs[0].StartPosition, Is.EqualTo(0));
            Assert.That(orfs[0].EndPosition, Is.EqualTo(8));
            Assert.That(orfs[0].Frame, Is.EqualTo(1));
            Assert.That(orfs[0].NucleotideLength, Is.EqualTo(9));
            Assert.That(orfs[0].AminoAcidLength, Is.EqualTo(2));
        });
    }

    [Test]
    public void FindOrfs_NullDna_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() => Translator.FindOrfs(null!).ToList());
    }

    [Test]
    public void FindOrfs_MultipleOrfs_FindsAll()
    {
        // Two ORFs in same frame: ATG GCT TAA | ATG GCT TAA
        // Source: Wikipedia ORF - multiple ORFs can exist in same sequence
        var dna = new DnaSequence("ATGGCTTAAATGGCTTAA");
        var orfs = Translator.FindOrfs(dna, minLength: 1, searchBothStrands: false).ToList();

        Assert.That(orfs, Has.Count.EqualTo(2));
        Assert.That(orfs[0].Protein.Sequence, Is.EqualTo("MA"));
        Assert.That(orfs[1].Protein.Sequence, Is.EqualTo("MA"));
    }

    #endregion

    #region Real Sequences

    [Test]
    public void Translate_RosalindPROT_SampleDataset_ProducesExpectedProtein()
    {
        // Rosalind PROT sample dataset.
        // Source: https://rosalind.info/problems/prot/
        // Input RNA: AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
        // Codons: AUG GCC AUG GCG CCC AGA ACU GAG AUC AAU AGU ACC CGU AUU AAC GGG UGA
        //         M   A   M   A   P   R   T   E   I   N   S   T   R   I   N   G   *
        // Rosalind excludes the terminating stop codon → MAMAPRTEINSTRING
        var rna = new RnaSequence("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");

        // toFirstStop:true follows the Rosalind convention (stop terminates, excluded).
        var protein = Translator.Translate(rna, toFirstStop: true);
        Assert.That(protein.Sequence, Is.EqualTo("MAMAPRTEINSTRING"));

        // Default (toFirstStop:false) shows the stop codon as '*'.
        var withStop = Translator.Translate(rna);
        Assert.That(withStop.Sequence, Is.EqualTo("MAMAPRTEINSTRING*"));
    }

    [Test]
    public void Translate_AmbiguousIupacCodon_ResolvedAsBiopython()
    {
        // Biopython Seq.translate (ambiguous codon tables): an IUPAC-ambiguous codon is
        // expanded; if every expansion gives the same amino acid that residue is emitted
        // (GCN -> A, since GCA/GCC/GCG/GCU are all Ala); a codon that may be a stop or
        // several unrelated residues (NNN) is emitted as 'X'.
        // Oracle (Biopython 1.88): Seq("AUGNNNGCNUAA").translate() == "MXA*".
        // (Review 2026-09 TRANS-CODON-001: previous expectation "MXX*" contradicted Biopython.)
        // NOTE: the typed DnaSequence/RnaSequence overloads reject IUPAC ambiguity
        // codes at construction, so this path is reachable only via the string overload.
        var protein = Translator.Translate("AUGNNNGCNUAA");
        Assert.That(protein.Sequence, Is.EqualTo("MXA*"));
    }

    [Test]
    public void Translate_InternalStop_DefaultConventionEmitsStarAndContinues()
    {
        // Default (toFirstStop:false) follows Biopython's default: an in-frame
        // stop is emitted as '*' and translation of the whole frame continues.
        // AUG (M) · UAA (*) · GCU (A) · UAG (*) · GGG (G)
        var rna = new RnaSequence("AUGUAAGCUUAGGGG");
        var withStops = Translator.Translate(rna);
        Assert.That(withStops.Sequence, Is.EqualTo("M*A*G"));

        // toFirstStop:true terminates at the first stop and excludes it.
        var toStop = Translator.Translate(rna, toFirstStop: true);
        Assert.That(toStop.Sequence, Is.EqualTo("M"));
    }

    [Test]
    public void Translate_TrailingPartialCodon_IsDropped()
    {
        // A trailing 1–2 nt remainder that cannot form a full codon is ignored
        // (Biopython six_frame_translations / EMBOSS transeq triplet truncation).
        // AUG (M) · GCU (A) · "CC" (2 nt remainder, dropped)
        var rna = new RnaSequence("AUGGCUCC");
        var protein = Translator.Translate(rna);
        Assert.That(protein.Sequence, Is.EqualTo("MA"));
    }

    [Test]
    public void Translate_InsulinBChain_ProducesCorrectProtein()
    {
        // Human insulin B chain coding sequence
        // Source: UniProt P01308, positions 25-54 of preproinsulin
        // DNA from NCBI RefSeq NM_000207.3
        var dna = new DnaSequence("TTCGTGAACCAGCACCTGTGCGGCTCCCACCTGGTGGAAGCTCTGTACCTGGTGTGTGGGGAGCGTGGCTTCTTCTACACACCCAAGACC");
        var protein = Translator.Translate(dna);

        Assert.That(protein.Sequence, Is.EqualTo("FVNQHLCGSHLVEALYLVCGERGFFYTPKT"));
    }

    #endregion

    #region Review 2026-09 — Biopython differential (TRANS-PROT-001)

    /// <summary>
    /// Differential oracle: Biopython 1.88 <c>Bio.Seq.translate(seq[frame:], table=t, to_stop=s)</c>
    /// for two sequences per NCBI table (one ACGT, one IUPAC-ambiguous, some RNA / lower case).
    /// Trailing partial codons are dropped as Biopython does. Values generated with Biopython 1.88.
    /// </summary>
    [TestCase(1, 2, true, "CTGTTTGCGTTC", "VCV")]
    [TestCase(1, 0, true, "AACAADTCGGNGGWTS", "NXSXX")]
    [TestCase(2, 1, true, "CTGGAGCCCGCAGTGCTC", "WSPQC")]
    [TestCase(2, 1, true, "uacaaydauggcdugghnuurckug", "TXMAWXXX")]
    [TestCase(3, 0, true, "CCAGGCGCTCCGTTG", "PGAPL")]
    [TestCase(3, 1, false, "GMGSVCHCGSVRHGCBAHA", "XXRXAX")]
    [TestCase(4, 2, false, "TAGTTGGACGTTCGAAGTTGAGTTTCCT", "VGRSKLSF")]
    [TestCase(4, 2, false, "aadaggwaugcgavukauggu", "XXCXXW")]
    [TestCase(5, 1, true, "CATAGTCTAGTAGTGTATCCCACCCC", "MV")]
    [TestCase(5, 1, false, "CTGTTCACAAGGTMBCABGYKCT", "CSQGXXX")]
    [TestCase(6, 2, true, "CCGACTAGTGAGGCTCCCACTTCAAAA", "DQ")]
    [TestCase(6, 0, true, "kuruhccgndybyhgvgsg", "XXRXXX")]
    [TestCase(9, 1, false, "CGTTCCGATAATAAAGGTCCACCTGTAAGGG", "VPIIKVHL*G")]
    [TestCase(9, 2, false, "ACGCCVTBTRVTCTWC", "AXXS")]
    [TestCase(10, 2, false, "TCTGTTTCTTATCACGGCTGAGATTTTGT", "CFLSRLRFC")]
    [TestCase(10, 1, true, "cucsahubyca", "SXX")]
    [TestCase(11, 2, false, "TGGAGAGCCAGTACGCTAGAGCCTTT", "ESQYARAF")]
    [TestCase(11, 0, true, "TCTCDTDTAGGCATANTDSGTRGHATMYGTAK", "SXXGIXXXIX")]
    [TestCase(12, 0, true, "CTACGCTGTCATAGCTCTCAGGACTCTCAGATG", "LRCHSSQDSQM")]
    [TestCase(12, 0, false, "uuguggawuyubwugrag", "LWXXXX")]
    [TestCase(13, 1, true, "GCCGCTCTCAAATCAAC", "PLSNQ")]
    [TestCase(13, 1, false, "TCWNANVCRYA", "XXX")]
    [TestCase(14, 1, true, "ATAGGATCGCCC", "")]
    [TestCase(14, 1, true, "auamcwccbgcagasry", "YXXQX")]
    [TestCase(15, 1, true, "AGGCCCATTACTTGCAGACAGGGTCTCG", "GPLLADRVS")]
    [TestCase(15, 1, false, "GTTATNSAMRARAKVCCTGGCACATTTAC", "LXXXXLAHL")]
    [TestCase(16, 2, true, "GCTAGCACTGTCAGTAGATAGCC", "LHCQLIA")]
    [TestCase(16, 1, false, "gbynugvuadsaasuacucvnanrhkcgaudb", "XXXZXLXXXX")]
    [TestCase(21, 1, false, "TGTCGGTGCCGC", "VGA")]
    [TestCase(21, 0, false, "CAATGAKKVTTACAGCB", "QWXLQ")]
    [TestCase(22, 2, false, "TTGGTCATAGCAATCAAAAACCGCAGGTCTCT", "GHSNQKPQVS")]
    [TestCase(22, 0, false, "uuvagcacgaangacckuha", "XSTXDX")]
    [TestCase(23, 0, true, "GGGACGCTTCTGATT", "GTLLI")]
    [TestCase(23, 0, false, "CCACCRGCDS", "PPA")]
    [TestCase(24, 0, true, "TGACAACGCGAAATCTGCAATTGG", "WQREICNW")]
    [TestCase(24, 0, true, "akcuamayurgvayacadvuuncwaahcgar", "XXXXXXXXXR")]
    [TestCase(25, 0, true, "AAGATAGTCGCA", "KIVA")]
    [TestCase(25, 0, true, "MBVGGVCGTGBMATCNCNNCYGK", "XGRXIXX")]
    [TestCase(26, 2, false, "GGCGCGTGCTAGA", "RVL")]
    [TestCase(26, 1, true, "cgugguacaugcgnuumv", "VVHAX")]
    [TestCase(27, 1, false, "TGAGCTCTAGCGTAACCTAG", "ELQRNL")]
    [TestCase(27, 0, false, "TNAGGGTDA", "XGX")]
    [TestCase(28, 2, false, "CATTGCCGTATTCAAACTT", "LPYSN")]
    [TestCase(28, 1, false, "waagdauggaacucucavaauc", "KXGTLXI")]
    [TestCase(29, 2, false, "GCTCGGGCTTCG", "SGF")]
    [TestCase(29, 1, true, "GATTCSBCHSGGTVCDCHTGCHTRG", "IXXGXXAX")]
    [TestCase(30, 0, false, "CACTCGGGCTCCCCCCGTCTTACCCTGA", "HSGSPRLTL")]
    [TestCase(30, 1, false, "cagbugruuaa", "XXL")]
    [TestCase(31, 1, false, "GAATTTGACTGT", "NLT")]
    [TestCase(31, 1, false, "HTMGHGATNCGAGTMGTKCKKG", "XXXEXXX")]
    [TestCase(32, 1, false, "CGAATCTACCTTAGAGACCGTAAGTTAACACG", "ESTLETVS*H")]
    [TestCase(32, 1, false, "wvghuucnungvsgh", "XFXX")]
    [TestCase(33, 0, false, "CTTGAATTACCTGTTCC", "LELPV")]
    [TestCase(33, 0, false, "ANCTCCNMYTTABTTYYARDCVTTKT", "XSXLXXXX")]
    public void Translate_String_AllTables_MatchBiopython(int table, int frame, bool toFirstStop, string sequence, string expected)
    {
        var protein = Translator.Translate(sequence, GeneticCode.GetByTableNumber(table), frame, toFirstStop);
        Assert.That(protein.Sequence, Is.EqualTo(expected));
    }

    /// <summary>
    /// Ambiguous codons that resolve to IUPAC ambiguous amino acids (RAY → B = Asx, SAR → Z = Glx,
    /// MTH → J = Xle) must be representable in the translated protein. Before the review fix,
    /// <see cref="ProteinSequence"/> rejected B/Z/J and Translate threw ArgumentException.
    /// Oracle: Biopython 1.88 <c>translate(seq, table)</c> / <c>translate(seq, table, to_stop=True)</c>.
    /// </summary>
    [TestCase("RAYSARMTH", 1, false, "BZJ")]
    [TestCase("ATGRAYSARMTHTAA", 1, false, "MBZJ*")]
    [TestCase("ATGRAYSARMTHTAA", 1, true, "MBZJ")]
    [TestCase("atgraysarmthtaa", 1, false, "MBZJ*")]
    [TestCase("AUGRAYSARMUHUAA", 1, false, "MBZJ*")]
    [TestCase("GCNRAYTARGGN", 1, false, "AB*G")]
    [TestCase("GCNRAYTARGGN", 1, true, "AB")]
    [TestCase("ATGRAYTAAMTH", 1, false, "MB*J")]
    [TestCase("ATGAGRRAY", 2, false, "M*B")]
    [TestCase("ATGAGRRAY", 2, true, "M")]
    [TestCase("WTARAY", 2, false, "XB")]
    [TestCase("CTNRAY", 3, false, "TB")]
    public void Translate_AmbiguousCodonsToAsxGlxXle_MatchBiopython(string sequence, int table, bool toFirstStop, string expected)
    {
        var protein = Translator.Translate(sequence, GeneticCode.GetByTableNumber(table), 0, toFirstStop);
        Assert.That(protein.Sequence, Is.EqualTo(expected));
    }

    /// <summary>
    /// NCBI tables 27, 28 and 31 contain dual-coding stop codons (translated as an amino acid,
    /// but also terminators in context). "Translate to first stop" is undefined for them:
    /// Biopython 1.88 raises ValueError ("You cannot use 'to_stop=True' with this table ...").
    /// Previously Seqeron silently read through the terminator (table 27 "ATGTGAGGG" → "MWG").
    /// </summary>
    [TestCase(27)]
    [TestCase(28)]
    [TestCase(31)]
    public void Translate_ToFirstStop_DualCodingStopTable_ThrowsArgumentException(int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => Translator.Translate("ATGTAAGGGTGA", code, 0, toFirstStop: true));
            Assert.Throws<ArgumentException>(() => Translator.Translate(new DnaSequence("ATGTAAGGGTGA"), code, 0, toFirstStop: true));
            Assert.Throws<ArgumentException>(() => Translator.Translate(new RnaSequence("AUGUAAGGGUGA"), code, 0, toFirstStop: true));
            Assert.Throws<ArgumentException>(() => Translator.Translate("", code, 0, toFirstStop: true));
        });
    }

    /// <summary>
    /// Without toFirstStop the dual-coding codons are translated as their amino acid
    /// (Biopython 1.88: table 27 "ATGTAAGGGTGA" → "MQGW", 28 → "MQGW", 31 → "MEGW").
    /// </summary>
    [TestCase(27, "MQGW")]
    [TestCase(28, "MQGW")]
    [TestCase(31, "MEGW")]
    public void Translate_DualCodingStopTable_WithoutToFirstStop_TranslatesAsAminoAcid(int table, string expected)
    {
        var protein = Translator.Translate("ATGTAAGGGTGA", GeneticCode.GetByTableNumber(table));
        Assert.That(protein.Sequence, Is.EqualTo(expected));
    }

    [Test]
    public void Translate_ToFirstStop_AllTablesWithoutDualCodingStops_DoNotThrow()
    {
        // Biopython 1.88: only tables 27, 28, 31 have a stop codon that is also in forward_table.
        foreach (int table in GeneticCode.SupportedTableNumbers.Where(t => t is not (27 or 28 or 31)))
        {
            var code = GeneticCode.GetByTableNumber(table);
            Assert.DoesNotThrow(() => Translator.Translate("ATGGCC", code, 0, toFirstStop: true), $"table {table}");
        }
    }

    [TestCase(-1)]
    [TestCase(3)]
    public void Translate_EmptyOrNullString_InvalidFrame_ThrowsLikeTypedOverloads(int frame)
    {
        // The string overload validates its arguments exactly like the DnaSequence overload.
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.Translate("", frame: frame));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.Translate((string)null!, frame: frame));
            Assert.Throws<ArgumentOutOfRangeException>(() => Translator.Translate(new DnaSequence(""), frame: frame));
        });
    }

    #endregion
}
