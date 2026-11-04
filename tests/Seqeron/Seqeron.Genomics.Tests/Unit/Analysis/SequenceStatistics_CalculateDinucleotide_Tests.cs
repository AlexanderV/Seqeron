// SEQ-DINUC-001 — Dinucleotide Analysis (frequencies, O/E relative abundance, codon frequencies)
// Evidence: docs/Evidence/SEQ-DINUC-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-DINUC-001.md
// Source: Karlin S. (1998). Pervasive properties of the genomic signature. PMC126251.
//         Karlin S, Burge C (1995). Trends Genet 11(7):283-290.
//         Nakamura Y et al. Kazusa Codon Usage Database (CUTG) readme.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateDinucleotide_Tests
{
    // Expected values are exact rationals derived from the Karlin odds-ratio
    // definition and count/total frequency definitions in the Evidence, computed
    // independently of the implementation.
    private const double Tolerance = 1e-10;

    #region CalculateDinucleotideRatios

    // M1 — rho_XY = f_XY/(f_X f_Y) for ATGCGCGT (A=1,T=2,G=3,C=2; N=8; 7 dinucleotide positions).
    // Evidence: Karlin PMC126251. rho_GC=rho_CG=(2/7)/((3/8)(2/8))=64/21; rho_AT=(1/7)/((1/8)(2/8))=32/7;
    //           rho_TG=rho_GT=(1/7)/((2/8)(3/8))=(1/7)/((3/8)(2/8))=32/21.
    [Test]
    public void CalculateDinucleotideRatios_AtgcgcgtSequence_ReturnsExactOddsRatios()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_GC = (2/7)/((3/8)(2/8)) = 64/21 per Karlin odds ratio (INV-03)");
            Assert.That(ratios["CG"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_CG = (2/7)/((2/8)(3/8)) = 64/21");
            Assert.That(ratios["AT"], Is.EqualTo(32.0 / 7.0).Within(Tolerance),
                "rho_AT = (1/7)/((1/8)(2/8)) = 32/7");
            Assert.That(ratios["TG"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_TG = (1/7)/((2/8)(3/8)) = 32/21");
            Assert.That(ratios["GT"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_GT = (1/7)/((3/8)(2/8)) = 32/21");
        });
    }

    // M7 — no-bias baseline: homopolymer AAAA has f_AA=1 and f_A=1, so rho_AA = 1/(1*1) = 1.0 exactly.
    // Evidence: Karlin PMC126251 (r=1 means no bias, INV-03).
    [Test]
    public void CalculateDinucleotideRatios_Homopolymer_ReturnsOddsRatioOne()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("AAAA");

        Assert.That(ratios["AA"], Is.EqualTo(1.0).Within(Tolerance),
            "rho_AA = f_AA/(f_A*f_A) = 1/(1*1) = 1.0, the no-bias baseline (INV-03)");
    }

    // S4 — robustness when some bases are absent: a present dinucleotide always has both of its
    //      constituent bases present, so the expected product f_X*f_Y is strictly positive and the
    //      ratio is finite (no division by zero) even when other bases (here G, C) never occur.
    //      "ATAT": A=2,T=2,G=0,C=0. Dinucs: AT,TA,AT. f_AT=2/3, f_A=2/4, f_T=2/4 -> rho_AT=(2/3)/(1/4)=8/3.
    //      f_TA=1/3, rho_TA=(1/3)/((2/4)(2/4))=4/3. Independently recomputed with exact rationals.
    // Evidence: Karlin PMC126251 (odds-ratio definition); no infinity arises because expected>0 for
    //           every dinucleotide that actually occurs.
    [Test]
    public void CalculateDinucleotideRatios_AbsentBaseContext_DoesNotProduceInfinity()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATAT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["AT"], Is.EqualTo(8.0 / 3.0).Within(Tolerance),
                "rho_AT = (2/3)/((2/4)(2/4)) = 8/3");
            Assert.That(ratios["TA"], Is.EqualTo(4.0 / 3.0).Within(Tolerance),
                "rho_TA = (1/3)/((2/4)(2/4)) = 4/3");
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "AT", "TA" }),
                "only dinucleotides over present bases appear; absent bases G and C yield no keys");
            Assert.That(ratios.Values, Has.All.Matches<double>(double.IsFinite),
                "no ratio is infinite/NaN even though G and C are absent");
        });
    }

    // S5 — non-alphabet (ambiguous) bases are excluded from ratio inputs. "ATGNCG": the N breaks
    //      the GN and NC contexts, leaving valid dinucleotides AT, TG, CG (each f = 1/3). Single-base
    //      frequencies are over the 5 unambiguous bases (A1,T1,G2,C1; N excluded), N_total = 5.
    //      rho_AT = (1/3)/((1/5)(1/5)) = 25/3; rho_TG = (1/3)/((1/5)(2/5)) = 25/6;
    //      rho_CG = (1/3)/((1/5)(2/5)) = 25/6. Independently recomputed with exact rationals.
    // Evidence: Karlin PMC126251 (normalized base/dinucleotide frequencies over the valid alphabet).
    [Test]
    public void CalculateDinucleotideRatios_AmbiguousBase_ExcludedFromCounts()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGNCG");

        Assert.Multiple(() =>
        {
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "AT", "TG", "CG" }),
                "GN and NC are excluded; only ACGT dinucleotides remain");
            Assert.That(ratios["AT"], Is.EqualTo(25.0 / 3.0).Within(Tolerance),
                "rho_AT = (1/3)/((1/5)(1/5)) = 25/3");
            Assert.That(ratios["TG"], Is.EqualTo(25.0 / 6.0).Within(Tolerance),
                "rho_TG = (1/3)/((1/5)(2/5)) = 25/6");
            Assert.That(ratios["CG"], Is.EqualTo(25.0 / 6.0).Within(Tolerance),
                "rho_CG = (1/3)/((1/5)(2/5)) = 25/6");
        });
    }

    // S1 — guards: null, empty, and length < 2 all return an empty dictionary.
    [Test]
    public void CalculateDinucleotideRatios_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(null!), Is.Empty,
                "null input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(string.Empty), Is.Empty,
                "empty input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios("A"), Is.Empty,
                "length < 2 yields empty dictionary (no dinucleotide positions)");
        });
    }

    // C1 — case-insensitivity: lowercase input produces the same ratios as uppercase (M1 values).
    [Test]
    public void CalculateDinucleotideRatios_LowercaseInput_MatchesUppercase()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("atgcgcgt");

        Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
            "lowercase is normalized via ToUpperInvariant; rho_GC equals the uppercase value 64/21");
    }

    #endregion

    #region CalculateDinucleotideFrequencies

    // M2 — normalized dinucleotide frequencies for ATGCGCGT (7 positions): count/(N-1).
    // Evidence: Karlin PMC126251 (f_XY normalized). GC=CG=2/7; AT=TG=GT=1/7.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_GC = 2/7 (2 of 7 dinucleotide positions)");
            Assert.That(freq["CG"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_CG = 2/7");
            Assert.That(freq["AT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_AT = 1/7");
            Assert.That(freq["TG"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_TG = 1/7");
            Assert.That(freq["GT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_GT = 1/7");
        });
    }

    // M3 — INV-01: dinucleotide frequencies sum to 1.0.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_SumsToOne()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
            "normalized frequencies over all valid dinucleotide positions sum to 1.0 (INV-01)");
    }

    // S2 — guards: null, empty, length < 2 return empty.
    [Test]
    public void CalculateDinucleotideFrequencies_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(null!), Is.Empty,
                "null input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(string.Empty), Is.Empty,
                "empty input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies("A"), Is.Empty,
                "length < 2 yields empty dictionary");
        });
    }

    // S5f — non-alphabet (ambiguous) bases excluded from frequency counts. "ATGNCG": dinucleotides
    //       GN and NC contain N and are dropped, leaving AT, TG, CG over 3 valid positions, each 1/3.
    //       Independently recomputed (count/(valid positions)).
    // Evidence: Karlin PMC126251 (frequency over valid dinucleotide positions); ambiguous units excluded.
    [Test]
    public void CalculateDinucleotideFrequencies_AmbiguousBase_ExcludedFromCounts()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGNCG");

        Assert.Multiple(() =>
        {
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "AT", "TG", "CG" }),
                "GN and NC contain the ambiguous base N and are excluded");
            Assert.That(freq["AT"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_AT = 1/3");
            Assert.That(freq["TG"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_TG = 1/3");
            Assert.That(freq["CG"], Is.EqualTo(1.0 / 3.0).Within(Tolerance), "f_CG = 1/3");
            Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
                "valid-position frequencies still sum to 1.0 (INV-01)");
        });
    }

    // C2 — RNA U handling: U is a valid base; AUGCGC includes the AU dinucleotide.
    [Test]
    public void CalculateDinucleotideFrequencies_RnaInput_IncludesUracilDinucleotides()
    {
        // AUGCGC -> dinucs AU,UG,GC,CG,GC (5). AU=1/5, GC=2/5.
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("AUGCGC");

        Assert.Multiple(() =>
        {
            Assert.That(freq.ContainsKey("AU"), Is.True,
                "U is in the dinucleotide alphabet, so AU is counted (RNA support)");
            Assert.That(freq["AU"], Is.EqualTo(1.0 / 5.0).Within(Tolerance),
                "f_AU = 1/5 (1 of 5 dinucleotide positions)");
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 5.0).Within(Tolerance),
                "f_GC = 2/5");
        });
    }

    #endregion

    #region Reference implementations (seqinr rho, EMBOSS compseq) and Karlin rho* (review 2026-09, F13)

    private static void AssertRatios(IReadOnlyDictionary<string, double> actual,
        IReadOnlyDictionary<string, double> expected, string source)
    {
        Assert.Multiple(() =>
        {
            Assert.That(actual.Keys, Is.EquivalentTo(expected.Keys), $"key set ({source})");
            foreach (var (k, v) in expected)
                Assert.That(actual.GetValueOrDefault(k, double.NaN), Is.EqualTo(v).Within(Tolerance), $"{k} ({source})");
        });
    }

    // EMBOSS 6.6.0 `compseq -word 2 -calcfreq` Obs/Exp (executed) = seqinr rho (cran/seqinr R/rho.R,
    // ported literally) for pure-ACGT input: 32/27 = 1.1851852, 16/9 = 1.7777778, 16/27 = 0.5925926.
    [Test]
    public void CalculateDinucleotideRatios_MatchesEmbossCompseqAndSeqinrRho()
    {
        var expected = new Dictionary<string, double>
        {
            ["AA"] = 32.0 / 27, ["AG"] = 32.0 / 27, ["AT"] = 16.0 / 9, ["CA"] = 32.0 / 27,
            ["CC"] = 16.0 / 9, ["CG"] = 16.0 / 27, ["CT"] = 16.0 / 27, ["GA"] = 16.0 / 27,
            ["GC"] = 16.0 / 9, ["GG"] = 32.0 / 27, ["TA"] = 32.0 / 27, ["TC"] = 16.0 / 27,
            ["TG"] = 16.0 / 27, ["TT"] = 16.0 / 9,
        };
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("GATCCTTAAAGGCGCATTTAGGCCCATG"),
            expected, "compseq/seqinr");
    }

    // seqinr rho/count: N, IUPAC codes and case are dropped from both numerator and denominators
    // (11 in-alphabet dinucleotides, 14 in-alphabet bases). ρ_CG = (3/11)/((3/14)(4/14)) = 49/11.
    [Test]
    public void CalculateDinucleotideRatios_AmbiguityAndLowercase_MatchesSeqinrRho()
    {
        const string seq = "ATGCGCGTNNacgtRAT";
        var expected = new Dictionary<string, double>
        {
            ["AC"] = 196.0 / 99, ["AT"] = 98.0 / 33, ["CG"] = 49.0 / 11,
            ["GC"] = 98.0 / 33, ["GT"] = 49.0 / 22, ["TG"] = 49.0 / 44,
        };
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios(seq), expected, "seqinr rho");

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        Assert.Multiple(() =>
        {
            Assert.That(freq["CG"], Is.EqualTo(3.0 / 11).Within(Tolerance), "seqinr count(freq=TRUE): 3/11");
            Assert.That(freq["AT"], Is.EqualTo(2.0 / 11).Within(Tolerance), "seqinr count(freq=TRUE): 2/11");
        });
    }

    // Karlin rho*: f*_A=f*_T=(f_A+f_T)/2, f*_C=f*_G=(f_C+f_G)/2, f*_GT=(f_GT+f_AC)/2 (Karlin PMC126251; restated
    // in PMC1829274, WebSearch snippet). Reference = seqinr rho on seq + separator + reverse complement (executed).
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_MatchesKarlinRhoStar()
    {
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGT", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AC"] = 128.0 / 105, ["AT"] = 256.0 / 63, ["CA"] = 128.0 / 105, ["CG"] = 512.0 / 175,
                ["GC"] = 512.0 / 175, ["GT"] = 128.0 / 105, ["TG"] = 128.0 / 105,
            }, "rho* ATGCGCGT");

        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("GATCCTTAAAGGCGCATTTAGGCCCATG", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AA"] = 40.0 / 27, ["AG"] = 8.0 / 9, ["AT"] = 16.0 / 9, ["CA"] = 8.0 / 9,
                ["CC"] = 40.0 / 27, ["CG"] = 16.0 / 27, ["CT"] = 8.0 / 9, ["GA"] = 16.0 / 27,
                ["GC"] = 16.0 / 9, ["GG"] = 40.0 / 27, ["TA"] = 32.0 / 27, ["TC"] = 16.0 / 27,
                ["TG"] = 8.0 / 9, ["TT"] = 40.0 / 27,
            }, "rho* 28-mer");

        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGTNNacgtRAT", strandSymmetric: true),
            new Dictionary<string, double>
            {
                ["AC"] = 24.0 / 11, ["AT"] = 32.0 / 11, ["CA"] = 8.0 / 11, ["CG"] = 48.0 / 11,
                ["GC"] = 32.0 / 11, ["GT"] = 24.0 / 11, ["TG"] = 8.0 / 11,
            }, "rho* with N/R/lowercase");
    }

    // rho* is strand-invariant: rho*_XY = rho*_{revcomp(XY)} and rho*(S) = rho*(revcomp(S)).
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_IsReverseComplementInvariant()
    {
        const string seq = "GATCCTTAAAGGCGCATTTAGGCCCATGNNACGTTGCA";
        var forward = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: true);
        var reverse = SequenceStatistics.CalculateDinucleotideRatios(
            DnaSequence.GetReverseComplementString(seq), strandSymmetric: true);

        AssertRatios(reverse, forward, "rho*(revcomp S)");
        Assert.Multiple(() =>
        {
            Assert.That(forward["GT"], Is.EqualTo(forward["AC"]).Within(Tolerance), "rho*_GT = rho*_AC");
            Assert.That(forward["AG"], Is.EqualTo(forward["CT"]).Within(Tolerance), "rho*_AG = rho*_CT");
        });
    }

    // rho* is a double-stranded DNA statistic: U is excluded like N. "AUGC": only GC counted
    // (symmetric count 2 of 2), bases A,G,C -> f*_C=f*_G=2/6 -> rho*_GC = 1/(1/9) = 9.
    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetric_ExcludesUracil()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("AUGC", strandSymmetric: true);
        Assert.Multiple(() =>
        {
            Assert.That(ratios.Keys, Is.EquivalentTo(new[] { "GC" }));
            Assert.That(ratios["GC"], Is.EqualTo(9.0).Within(Tolerance));
        });
    }

    [Test]
    public void CalculateDinucleotideRatios_StrandSymmetricFalse_EqualsSingleArgumentOverload()
    {
        foreach (var seq in new[] { "ATGCGCGT", "AUGCGC", "ATGCGCGTNNacgtRAT", "A", "", "NNNN" })
            AssertRatios(SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: false),
                SequenceStatistics.CalculateDinucleotideRatios(seq), seq);
    }

    // D3: counting delegates to the canonical overlapping k-mer counter (k = 2), filtered to {A,C,G,T,U}.
    [Test]
    public void CalculateDinucleotideFrequencies_EqualsCanonicalKmerCounterFiltered()
    {
        var rng = new Random(20260928);
        const string alphabet = "ACGTUNacgtu-RY";
        for (int n = 0; n < 200; n++)
        {
            var seq = new string(Enumerable.Range(0, rng.Next(2, 60)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            var counts = KmerAnalyzer.CountKmers(seq, 2)
                .Where(kv => kv.Key.All(c => "ACGTU".Contains(c)))
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            int total = counts.Values.Sum();
            var expected = counts.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);
            AssertRatios(SequenceStatistics.CalculateDinucleotideFrequencies(seq), expected, seq);
        }
    }

    #endregion
    #region Karlin representation classes (FIN-B03 WP11, audit row A19, F37)

    // Criterion (Karlin & Burge 1995, quoted verbatim in MBE 19(6):964 — SEQ-DINUC-001 Evidence):
    // "underrepresented if ρ*XY ≤ 0.78 and overrepresented if ρ*XY ≥ 1.23" (bounds inclusive).
    [TestCase(0.0, SequenceStatistics.DinucleotideRepresentation.UnderRepresented)]
    [TestCase(0.78, SequenceStatistics.DinucleotideRepresentation.UnderRepresented)]
    [TestCase(0.78000000000000014, SequenceStatistics.DinucleotideRepresentation.Normal)]
    [TestCase(1.0, SequenceStatistics.DinucleotideRepresentation.Normal)]
    [TestCase(1.2299999999999998, SequenceStatistics.DinucleotideRepresentation.Normal)]
    [TestCase(1.23, SequenceStatistics.DinucleotideRepresentation.OverRepresented)]
    [TestCase(double.PositiveInfinity, SequenceStatistics.DinucleotideRepresentation.OverRepresented)]
    public void ClassifyDinucleotideRatio_KarlinThresholds(double rho, SequenceStatistics.DinucleotideRepresentation expected)
    {
        Assert.That(SequenceStatistics.ClassifyDinucleotideRatio(rho), Is.EqualTo(expected));
    }

    [TestCase(double.NaN)]
    [TestCase(-0.1)]
    public void ClassifyDinucleotideRatio_InvalidRatio_Throws(double rho)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceStatistics.ClassifyDinucleotideRatio(rho));
    }

    // ρ* of ATGCATGCAATTGGCCATGCAAACGTTTAAGG computed independently in Python (seqinr rho of
    // seq + separator + reverse complement): AA 1.4273197929111907, AC 0.5243215565796211,
    // AG 0.26216077828981055, AT 1.6312226204699323, CA 2.0972862263184844, CC 1.011191573403555,
    // CG 0.67412771560237, CT 0.26216077828981055, GC 2.69651086240948, GG 1.011191573403555,
    // GT 0.5243215565796211, TA 0.40780565511748307, TG 2.0972862263184844, TT 1.4273197929111907.
    [Test]
    public void ClassifyDinucleotideRatios_StrandSymmetric_MatchesKarlinClassesOfIndependentRhoStar()
    {
        const string seq = "ATGCATGCAATTGGCCATGCAAACGTTTAAGG";
        var under = SequenceStatistics.DinucleotideRepresentation.UnderRepresented;
        var over = SequenceStatistics.DinucleotideRepresentation.OverRepresented;
        var normal = SequenceStatistics.DinucleotideRepresentation.Normal;
        var expected = new Dictionary<string, SequenceStatistics.DinucleotideRepresentation>
        {
            ["AA"] = over, ["AC"] = under, ["AG"] = under, ["AT"] = over, ["CA"] = over, ["CC"] = normal,
            ["CG"] = under, ["CT"] = under, ["GC"] = over, ["GG"] = normal, ["GT"] = under, ["TA"] = under,
            ["TG"] = over, ["TT"] = over,
        };

        var classes = SequenceStatistics.ClassifyDinucleotideRatios(seq);
        var ratios = SequenceStatistics.CalculateDinucleotideRatios(seq, strandSymmetric: true);

        Assert.Multiple(() =>
        {
            Assert.That(classes, Is.EquivalentTo(expected));
            Assert.That(ratios["CG"], Is.EqualTo(0.67412771560237).Within(1e-12));
            Assert.That(ratios["TA"], Is.EqualTo(0.40780565511748307).Within(1e-12));
        });
    }

    // Single-strand ρ (Python): AAAATTTTCCCCGGGGACGT -> AA/CC/GG/TT 2.526315789473684 (over),
    // CG 1.6842105263157894 (over), AC/AT/GA/GT/TC 0.8421052631578947 (normal).
    [Test]
    public void ClassifyDinucleotideRatios_SingleStrand_KeysAndClassesFollowRatioMap()
    {
        var classes = SequenceStatistics.ClassifyDinucleotideRatios("AAAATTTTCCCCGGGGACGT", strandSymmetric: false);
        var over = SequenceStatistics.DinucleotideRepresentation.OverRepresented;
        var normal = SequenceStatistics.DinucleotideRepresentation.Normal;

        Assert.Multiple(() =>
        {
            Assert.That(classes.Keys, Is.EquivalentTo(new[] { "AA", "AC", "AT", "CC", "CG", "GA", "GG", "GT", "TC", "TT" }));
            foreach (var d in new[] { "AA", "CC", "GG", "TT", "CG" }) Assert.That(classes[d], Is.EqualTo(over), d);
            foreach (var d in new[] { "AC", "AT", "GA", "GT", "TC" }) Assert.That(classes[d], Is.EqualTo(normal), d);
            Assert.That(SequenceStatistics.ClassifyDinucleotideRatios(""), Is.Empty);
            Assert.That(SequenceStatistics.ClassifyDinucleotideRatios("A"), Is.Empty);
        });
    }

    #endregion

    #region Absent dinucleotides, rho = 0 (FIN-B03 R2-WP4, round 2 B7, F51)

    // seqinr 4.2.36 rho() executed 2026-10-14 (R 4, Debian r-cran-seqinr): every word over the occurring
    // bases is reported; an unobserved word with both bases present is 0, a word with an absent base NaN.
    private static void AssertRatios(IReadOnlyDictionary<string, double> actual, IReadOnlyDictionary<string, double> expected)
    {
        Assert.That(actual.Keys, Is.EquivalentTo(expected.Keys));
        Assert.Multiple(() =>
        {
            foreach (var (k, v) in expected)
                Assert.That(actual[k], Is.EqualTo(v).Within(1e-12), k);
        });
    }

    [Test]
    public void DinucleotideRatios_IncludeAbsent_SingleStrand_MatchesSeqinrRhoAll16()
    {
        const double hi = 32.0 / 15; // seqinr 2.1333333333333333
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("AACCTTGGAACCTTGG", strandSymmetric: false, includeAbsent: true),
            new Dictionary<string, double>
            {
                ["AA"] = hi, ["AC"] = hi, ["CC"] = hi, ["CT"] = hi, ["GG"] = hi, ["TG"] = hi, ["TT"] = hi,
                ["GA"] = 16.0 / 15, // seqinr 1.0666666666666667
                ["AG"] = 0, ["AT"] = 0, ["CA"] = 0, ["CG"] = 0, ["GC"] = 0, ["GT"] = 0, ["TA"] = 0, ["TC"] = 0,
            });
    }

    [Test]
    public void DinucleotideRatios_IncludeAbsent_StrandSymmetric_MatchesSeqinrRhoOfSeqPlusReverseComplement()
    {
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("AACCTTGGAACCTTGG", strandSymmetric: true, includeAbsent: true),
            new Dictionary<string, double>
            {
                ["AA"] = 32.0 / 15, ["CC"] = 32.0 / 15, ["GG"] = 32.0 / 15, ["TT"] = 32.0 / 15,
                ["AC"] = 16.0 / 15, ["AG"] = 16.0 / 15, ["CA"] = 16.0 / 15, ["CT"] = 16.0 / 15, ["GT"] = 16.0 / 15, ["TG"] = 16.0 / 15,
                ["GA"] = 8.0 / 15, ["TC"] = 8.0 / 15, // seqinr 0.53333333333333333
                ["AT"] = 0, ["CG"] = 0, ["GC"] = 0, ["TA"] = 0,
            });

        // F37 sequence: seqinr rho* reports GA = TC = 0; the other 14 values are the F37 ones.
        var r = SequenceStatistics.CalculateDinucleotideRatios("ATGCATGCAATTGGCCATGCAAACGTTTAAGG", true, true);
        Assert.That(r, Has.Count.EqualTo(16));
        Assert.Multiple(() =>
        {
            Assert.That(r["GA"], Is.EqualTo(0.0));
            Assert.That(r["TC"], Is.EqualTo(0.0));
            Assert.That(r["GC"], Is.EqualTo(2.69651086240947979).Within(1e-12));
            Assert.That(r["AG"], Is.EqualTo(0.26216077828981055).Within(1e-12));
        });
    }

    // seqinr: rho(s2c("aaaatttt")) → cg etc. NaN (C and G absent): undefined, so not keys; TA = 0 is.
    // seqinr: rho(s2c("aacng")) → aa 2, ac 4, ag/ca/cc/cg/ga/gc/gg 0, words with t NaN.
    [Test]
    public void DinucleotideRatios_IncludeAbsent_UndefinedRatiosExcluded_AsSeqinrNaN()
    {
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("AAAATTTT", false, true),
            new Dictionary<string, double> { ["AA"] = 12.0 / 7, ["AT"] = 4.0 / 7, ["TA"] = 0, ["TT"] = 12.0 / 7 });
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("AAAATTTT", true, true),
            new Dictionary<string, double> { ["AA"] = 12.0 / 7, ["AT"] = 4.0 / 7, ["TA"] = 0, ["TT"] = 12.0 / 7 });
        AssertRatios(SequenceStatistics.CalculateDinucleotideRatios("aacng", false, true),
            new Dictionary<string, double>
            {
                ["AA"] = 2, ["AC"] = 4, ["AG"] = 0, ["CA"] = 0, ["CC"] = 0, ["CG"] = 0, ["GA"] = 0, ["GC"] = 0, ["GG"] = 0,
            });

        // seqinr rho(s2c("an")): all NaN (no counted dinucleotide) → nothing added.
        Assert.That(SequenceStatistics.CalculateDinucleotideRatios("AN", false, true), Is.Empty);
        Assert.That(SequenceStatistics.CalculateDinucleotideRatios("A", true, true), Is.Empty);
        Assert.That(SequenceStatistics.CalculateDinucleotideRatios("", false, true), Is.Empty);
    }

    // Karlin: rho <= 0.78 under-represented; rho = 0 (CpG-free) is the strongest under-representation.
    [Test]
    public void ClassifyDinucleotideRatios_IncludeAbsent_CpGFree_CGUnderRepresented()
    {
        var withAbsent = SequenceStatistics.ClassifyDinucleotideRatios("AACCTTGGAACCTTGG", strandSymmetric: true, includeAbsent: true);
        var observedOnly = SequenceStatistics.ClassifyDinucleotideRatios("AACCTTGGAACCTTGG", strandSymmetric: true);
        Assert.Multiple(() =>
        {
            Assert.That(withAbsent, Has.Count.EqualTo(16));
            foreach (var d in new[] { "AT", "CG", "GC", "TA", "GA", "TC" })
                Assert.That(withAbsent[d], Is.EqualTo(SequenceStatistics.DinucleotideRepresentation.UnderRepresented), d);
            foreach (var d in new[] { "AA", "CC", "GG", "TT" })
                Assert.That(withAbsent[d], Is.EqualTo(SequenceStatistics.DinucleotideRepresentation.OverRepresented), d);
            Assert.That(observedOnly.ContainsKey("CG"), Is.False, "default keeps observed dinucleotides only");
        });

        var f37 = SequenceStatistics.ClassifyDinucleotideRatios("ATGCATGCAATTGGCCATGCAAACGTTTAAGG", true, true);
        Assert.That(f37["GA"], Is.EqualTo(SequenceStatistics.DinucleotideRepresentation.UnderRepresented));
        Assert.That(f37["TC"], Is.EqualTo(SequenceStatistics.DinucleotideRepresentation.UnderRepresented));
    }

    [Test]
    public void DinucleotideRatios_IncludeAbsent_DefaultUnchanged_AndAddsOnlyZerosOverOccurringBases()
    {
        var rng = new Random(51);
        const string alphabet = "ACGTUNacgt-";
        for (int iteration = 0; iteration < 400; iteration++)
        {
            var chars = new char[rng.Next(0, 40)];
            for (int i = 0; i < chars.Length; i++) chars[i] = alphabet[rng.Next(rng.Next(2, alphabet.Length + 1))];
            string seq = new(chars);
            foreach (bool symmetric in new[] { false, true })
            {
                var baseline = SequenceStatistics.CalculateDinucleotideRatios(seq, symmetric);
                Assert.That(SequenceStatistics.CalculateDinucleotideRatios(seq, symmetric, includeAbsent: false), Is.EquivalentTo(baseline), seq);
                Assert.That(SequenceStatistics.ClassifyDinucleotideRatios(seq, symmetric, includeAbsent: false),
                    Is.EquivalentTo(SequenceStatistics.ClassifyDinucleotideRatios(seq, symmetric)), seq);

                var complete = SequenceStatistics.CalculateDinucleotideRatios(seq, symmetric, includeAbsent: true);
                string upper = seq.ToUpperInvariant();
                foreach (var (k, v) in complete)
                {
                    if (baseline.TryGetValue(k, out double observed))
                    {
                        Assert.That(v, Is.EqualTo(observed), seq + " " + k);
                        continue;
                    }

                    Assert.That(v, Is.EqualTo(0.0), seq + " " + k);
                    foreach (char b in k)
                    {
                        bool occurs = symmetric
                            ? (b is 'A' or 'T' ? upper.IndexOfAny(new[] { 'A', 'T' }) >= 0 : upper.IndexOfAny(new[] { 'C', 'G' }) >= 0)
                            : upper.Contains(b);
                        Assert.That(occurs, Is.True, $"{seq} {k}: base {b} must occur");
                    }
                }

                if (baseline.Count > 0)
                {
                    int bases = symmetric
                        ? (upper.IndexOfAny(new[] { 'A', 'T' }) >= 0 ? 2 : 0) + (upper.IndexOfAny(new[] { 'C', 'G' }) >= 0 ? 2 : 0)
                        : "ACGTU".Count(upper.Contains);
                    Assert.That(complete, Has.Count.EqualTo(bases * bases), seq);
                }
                else
                {
                    Assert.That(complete, Is.Empty, seq);
                }
            }
        }
    }

    #endregion
}
