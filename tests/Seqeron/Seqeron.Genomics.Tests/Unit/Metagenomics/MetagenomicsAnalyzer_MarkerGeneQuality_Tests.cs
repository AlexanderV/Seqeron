// META-BIN-001 — CheckM-style marker-gene completeness / contamination (opt-in addendum)
// Evidence: docs/Evidence/META-BIN-001-MarkerQC-Evidence.md
// TestSpec: tests/TestSpecs/META-BIN-001.md
// Source: Parks DH et al. (2015). CheckM. Genome Res 25:1043–1055.
//         Reference impl: Ecogenomics/CheckM checkm/markerSets.py (MarkerSet.genomeCheck).
namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

/// <summary>
/// Tests for the CheckM-style single-copy marker-gene completeness/contamination metrics added to
/// <see cref="MetagenomicsAnalyzer"/> (META-BIN-001 addendum). Completeness/contamination are the
/// per-collocated-marker-set averages of Parks et al. (2015) Eqs. 1–2:
///   Completeness = 100·(1/|M|)·Σ_{s∈M} |s∩G_M|/|s|;
///   Contamination = 100·(1/|M|)·Σ_{s∈M} Σ_{g∈s}(N_g−1)/|s|.
/// </summary>
[TestFixture]
[Category("Metagenomics")]
[Category("META-BIN-001")]
public class MetagenomicsAnalyzer_MarkerGeneQuality_Tests
{
    // E. coli small ribosomal subunit protein uS8 (RpsH), UniProt P0A7W7 — the true positive for
    // Pfam PF00410 (Ribosomal_S8). Retrieved 2026-06-25 from rest.uniprot.org/uniprotkb/P0A7W7.fasta.
    private const string EcoliRpsH_uS8 =
        "MSMQDPIADMLTRIRNGQAANKAAVTMPSSKLKVAIANVLKEEGFIEDFKVEGDTKPELE" +
        "LTLKYFQGKAVVESIQRVSRPGLRIYKRKDELPKVMAGLGIAVVSTSKGVMTDRAARQAG" +
        "LGGEIICYVA";

    // E. coli K-12 chaperone GrpE, UniProt P09372 — the true positive for the bac120 Pfam marker
    // PF01025 (GrpE). Retrieved 2026-06-25 from rest.uniprot.org/uniprotkb/P09372.fasta.
    private const string EcoliGrpE =
        "MSSKEQKTPEGQAPEEIIMDQHEEIEAVEPEASAEQVDPRDEKVANLEAQLAEAQTRERD" +
        "GILRVKAEMENLRRRTELDIEKAHKFALEKFINELLPVIDSLDRALEVADKANPDMSAMV" +
        "EGIELTLKSMLDVVRKFGVEVIAETNVPLDPNVHQAIAMVESDDVAPGNVLGIMQKGYTL" +
        "NGRTIRAAMVTVAKAKA";

    private static MetagenomicsAnalyzer.MarkerSet Set(params string[] ids)
        => new(ids);

    #region EstimateBinQualityFromMarkerCounts — exact CheckM formula

    // M-FORMULA — pins Parks 2015 Eqs. 1–2 on the hand-derived synthetic bin.
    // M = { {A,B}, {C,D,E}, {F} }; counts A=1,B=0,C=2,D=1,E=1,F=1.
    //   s1: present 1/2, multiCopy 0/2 ; s2: present 3/3, multiCopy 1/3 ; s3: 1/1, 0/1.
    //   comp = 0.5+1+1 = 2.5 ⇒ 100·2.5/3 = 250/3 = 83.333…% ;
    //   cont = 0+1/3+0 = 1/3 ⇒ 100·(1/3)/3 = 100/9 = 11.111…%.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_SyntheticBin_MatchesHandDerivedCheckMFormula()
    {
        var sets = new[] { Set("A", "B"), Set("C", "D", "E"), Set("F") };
        var counts = new Dictionary<string, int>
        {
            ["A"] = 1, ["B"] = 0, ["C"] = 2, ["D"] = 1, ["E"] = 1, ["F"] = 1,
        };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(250.0 / 3.0).Within(1e-10),
                "Completeness = 100·(0.5+1+1)/3 = 250/3 % per CheckM Eq.1 on the hand-derived bin.");
            Assert.That(q.Contamination, Is.EqualTo(100.0 / 9.0).Within(1e-10),
                "Contamination = 100·(0+1/3+0)/3 = 100/9 % per CheckM Eq.2 (C_g=N-1 for the duplicated C).");
            Assert.That(q.MarkerSetCount, Is.EqualTo(3), "|M| = 3 non-empty marker sets.");
            Assert.That(q.MarkerCount, Is.EqualTo(6), "Total markers across sets = 2+3+1 = 6.");
            Assert.That(q.MarkersPresent, Is.EqualTo(5),
                "Distinct markers found ≥1 = {A,C,D,E,F}; B is absent.");
        });
    }

    // M-ALLPRESENT — every single-copy marker present exactly once ⇒ perfect genome.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_AllSingleCopyPresent_Is100And0()
    {
        var sets = new[] { Set("A", "B"), Set("C") };
        var counts = new Dictionary<string, int> { ["A"] = 1, ["B"] = 1, ["C"] = 1 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(100.0).Within(1e-10),
                "All markers present single-copy ⇒ each set fully recovered ⇒ 100% complete.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10),
                "No multi-copy markers ⇒ 0% contamination.");
        });
    }

    // M-ALLABSENT — no markers found ⇒ empty genome.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_NoMarkersPresent_Is0And0()
    {
        var sets = new[] { Set("A", "B"), Set("C") };
        var counts = new Dictionary<string, int>(); // nothing found

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(0.0).Within(1e-10),
                "No markers found ⇒ 0% completeness.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10),
                "No markers found ⇒ 0% contamination.");
            Assert.That(q.MarkersPresent, Is.EqualTo(0), "No distinct markers present.");
        });
    }

    // M-TRIPLICATE — C_g = N-1 = 2 extra copies of a marker in a singleton set ⇒ 200% contamination.
    // Set {A} with A found 3 times: present 1/1, multiCopy (3-1)/1 = 2 ⇒ comp 100%, cont 100·2/1 = 200%.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_TriplicatedMarker_ContaminationIsTwoExtraCopies()
    {
        var sets = new[] { Set("A") };
        var counts = new Dictionary<string, int> { ["A"] = 3 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(100.0).Within(1e-10),
                "A multi-copy marker still counts once toward present ⇒ 100% complete.");
            Assert.That(q.Contamination, Is.EqualTo(200.0).Within(1e-10),
                "C_A = N-1 = 2 extra copies; contamination = 100·2/1 = 200% per CheckM Eq.2.");
        });
    }

    // M-EMPTYSET — an empty marker set must be skipped (it would divide by |s|=0) and not counted in |M|.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_EmptyMarkerSet_IsIgnored()
    {
        var sets = new[] { Set("A"), Set( /* empty */ ) };
        var counts = new Dictionary<string, int> { ["A"] = 1 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        Assert.Multiple(() =>
        {
            Assert.That(q.MarkerSetCount, Is.EqualTo(1),
                "The empty set is excluded from |M| (no div-by-zero).");
            Assert.That(q.Completeness, Is.EqualTo(100.0).Within(1e-10),
                "Only the non-empty fully-recovered set contributes ⇒ 100%.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10), "No multi-copy markers ⇒ 0%.");
        });
    }

    // M-NOMARKERSETS — |M| = 0 ⇒ undefined; return 0/0 rather than NaN/div-by-zero.
    [Test]
    public void EstimateBinQualityFromMarkerCounts_NoMarkerSets_Is0And0()
    {
        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
            System.Array.Empty<MetagenomicsAnalyzer.MarkerSet>(),
            new Dictionary<string, int>());

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(0.0).Within(1e-10), "|M|=0 ⇒ 0% (no div-by-zero).");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10), "|M|=0 ⇒ 0%.");
            Assert.That(q.MarkerSetCount, Is.EqualTo(0), "No marker sets used.");
        });
    }

    [Test]
    public void EstimateBinQualityFromMarkerCounts_NullArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
                    null!, new Dictionary<string, int>()),
                NUnit.Framework.Throws.ArgumentNullException, "Null markerSets must throw.");
            Assert.That(() => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
                    System.Array.Empty<MetagenomicsAnalyzer.MarkerSet>(), null!),
                NUnit.Framework.Throws.ArgumentNullException, "Null markerCounts must throw.");
        });
    }

    #endregion

    #region DetectMarkers / EstimateBinQualityFromMarkers — Plan7 HMM integration

    // H-TRUEPOSITIVE — the bundled CC0 PF00410 (Ribosomal_S8) HMM detects E. coli uS8, and NO other
    // bundled ribosomal family does (specificity). Reuses the verified Plan7ProfileHmm engine.
    [Test]
    public void DetectMarkers_EcoliuS8_HitsOnlyPF00410()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms();

        var counts = MetagenomicsAnalyzer.DetectMarkers(new[] { EcoliRpsH_uS8 }, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(counts["PF00410"], Is.EqualTo(1),
                "E. coli uS8 is the true positive for Ribosomal_S8 (PF00410) ⇒ exactly 1 hit.");
            foreach (var marker in hmms.Where(m => m.MarkerId != "PF00410"))
                Assert.That(counts[marker.MarkerId], Is.EqualTo(0),
                    $"uS8 must NOT match a different ribosomal family ({marker.MarkerId}).");
        });
    }

    // H-COMPLETENESS — end-to-end: 9 bundled markers as singleton sets, only PF00410 present
    // ⇒ exactly 1 of 9 sets recovered ⇒ Completeness = 100/9 %, Contamination = 0.
    [Test]
    public void EstimateBinQualityFromMarkers_BinWithOnlyuS8_CompletenessIsOneOfNine()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledRibosomalMarkerSets();

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkers(
            new[] { EcoliRpsH_uS8 }, sets, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(q.MarkerSetCount, Is.EqualTo(9), "Nine bundled singleton marker sets.");
            Assert.That(q.Completeness, Is.EqualTo(100.0 / 9.0).Within(1e-10),
                "Only PF00410 present (1 of 9 singleton sets) ⇒ Completeness = 100/9 %.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10),
                "Single copy of the one detected marker ⇒ 0% contamination.");
            Assert.That(q.MarkersPresent, Is.EqualTo(1), "Exactly one bundled marker detected.");
        });
    }

    // H-DUPLICATE — two copies of uS8 ⇒ PF00410 detected twice ⇒ that singleton set contributes
    // (2-1)/1 contamination ⇒ Contamination = 100·1/9 = 100/9 %.
    [Test]
    public void EstimateBinQualityFromMarkers_TwoCopiesOfuS8_ContaminationReflectsDuplication()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledRibosomalMarkerSets();

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkers(
            new[] { EcoliRpsH_uS8, EcoliRpsH_uS8 }, sets, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(q.Completeness, Is.EqualTo(100.0 / 9.0).Within(1e-10),
                "Duplicate still counts once toward present ⇒ Completeness unchanged at 100/9 %.");
            Assert.That(q.Contamination, Is.EqualTo(100.0 / 9.0).Within(1e-10),
                "PF00410 found N=2 ⇒ C=N-1=1 extra copy in a singleton set ⇒ 100·1/9 %.");
        });
    }

    [Test]
    public void LoadBundledRibosomalMarkerHmms_LoadsNineMarkersWithGatheringThresholds()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms();

        Assert.Multiple(() =>
        {
            Assert.That(hmms, Has.Count.EqualTo(9), "Nine bundled CC0 ribosomal-protein marker HMMs.");
            // PF00410 GA1 = 24 bits (from the embedded HMM 'GA 24 24;' line).
            var s8 = hmms.Single(m => m.MarkerId == "PF00410");
            Assert.That(s8.BitScoreThreshold, Is.EqualTo(24.0).Within(1e-9),
                "PF00410 uses its Pfam GA1 gathering threshold (24 bits).");
            Assert.That(s8.Hmm.Length, Is.EqualTo(125), "PF00410 LENG = 125 match states.");
        });
    }

    [Test]
    public void DetectMarkers_NullArguments_Throw()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms();
        Assert.Multiple(() =>
        {
            Assert.That(() => MetagenomicsAnalyzer.DetectMarkers(null!, hmms),
                NUnit.Framework.Throws.ArgumentNullException, "Null proteins must throw.");
            Assert.That(() => MetagenomicsAnalyzer.DetectMarkers(new[] { EcoliRpsH_uS8 }, null!),
                NUnit.Framework.Throws.ArgumentNullException, "Null markerHmms must throw.");
        });
    }

    #endregion

    #region Caller-supplied loader

    // L-CALLERSUPPLIED — a user with the full CheckM data supplies their own HMM via a reader;
    // it is keyed by accession and usable for detection. Round-trip a bundled profile through the
    // public Parse/loader path to prove the caller-supplied loader works.
    [Test]
    public void LoadMarkerHmms_CallerSuppliedProfile_KeyedByAccessionAndDetects()
    {
        // Obtain the raw PF00410 profile text via the same public Plan7 engine the loader uses,
        // then re-load it through the caller-supplied reader path.
        string hmmText = ReadBundledHmmText("PF00410_Ribosomal_S8.hmm");
        using var reader = new System.IO.StringReader(hmmText);

        var loaded = MetagenomicsAnalyzer.LoadMarkerHmms(new[] { reader });

        Assert.Multiple(() =>
        {
            Assert.That(loaded, Has.Count.EqualTo(1), "One caller-supplied marker loaded.");
            Assert.That(loaded[0].MarkerId, Is.EqualTo("PF00410.25"),
                "Marker is keyed by the profile ACC (PF00410.25).");
            var counts = MetagenomicsAnalyzer.DetectMarkers(new[] { EcoliRpsH_uS8 }, loaded);
            Assert.That(counts["PF00410.25"], Is.EqualTo(1),
                "The caller-supplied S8 profile detects E. coli uS8.");
        });
    }

    [Test]
    public void LoadMarkerHmms_NullReaders_Throws()
        => Assert.That(() => MetagenomicsAnalyzer.LoadMarkerHmms(null!),
            NUnit.Framework.Throws.ArgumentNullException, "Null hmmReaders must throw.");

    private static string ReadBundledHmmText(string resourceFile)
    {
        var asm = typeof(MetagenomicsAnalyzer).Assembly;
        using var stream = asm.GetManifestResourceStream(
            $"Seqeron.Genomics.Metagenomics.Resources.{resourceFile}")!;
        using var sr = new System.IO.StreamReader(stream);
        return sr.ReadToEnd();
    }

    #endregion

    #region Bundled GTDB domain-level universal marker sets (bac120 / ar122 Pfam subsets)

    // B-SETSIZE — the bundled bac120 Pfam subset is exactly the 6 CC0 Pfam markers of GTDB bac120
    // (PF00380, PF00410, PF00466, PF01025, PF02576, PF03726); GA1-gated; one singleton set each.
    [Test]
    public void LoadBundledBacterialMarkerHmms_LoadsSixBac120PfamMarkers()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledBacterialMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledBacterialMarkerSets();

        var expected = new[] { "PF00380", "PF00410", "PF00466", "PF01025", "PF02576", "PF03726" };
        Assert.Multiple(() =>
        {
            Assert.That(hmms.Select(m => m.MarkerId), Is.EquivalentTo(expected),
                "Bundled bac120 Pfam subset = the 6 CC0 Pfam markers of GTDB bac120.");
            Assert.That(sets, Has.Count.EqualTo(6), "One singleton collocated set per bac120 Pfam marker.");
            // PF01025 (GrpE) GA1 = 25.8 bits, LENG = 165 (from the embedded HMM).
            var grpe = hmms.Single(m => m.MarkerId == "PF01025");
            Assert.That(grpe.BitScoreThreshold, Is.EqualTo(25.8).Within(1e-9),
                "PF01025 uses its Pfam GA1 gathering threshold (25.8 bits).");
            Assert.That(grpe.Hmm.Length, Is.EqualTo(165), "PF01025 LENG = 165 match states.");
        });
    }

    // A-SETSIZE — the bundled ar122 Pfam subset is exactly the 35 CC0 Pfam markers of GTDB ar122.
    [Test]
    public void LoadBundledArchaealMarkerHmms_LoadsThirtyFiveAr122PfamMarkers()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledArchaealMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledArchaealMarkerSets();

        var expected = new[]
        {
            "PF00368", "PF00410", "PF00466", "PF00687", "PF00827", "PF00900", "PF01000", "PF01015",
            "PF01090", "PF01092", "PF01157", "PF01191", "PF01194", "PF01198", "PF01200", "PF01269",
            "PF01280", "PF01282", "PF01496", "PF01655", "PF01798", "PF01864", "PF01866", "PF01868",
            "PF01984", "PF01990", "PF02006", "PF02978", "PF03874", "PF04019", "PF04104", "PF04919",
            "PF07541", "PF13656", "PF13685",
        };
        Assert.Multiple(() =>
        {
            Assert.That(hmms.Select(m => m.MarkerId), Is.EquivalentTo(expected),
                "Bundled ar122 Pfam subset = the 35 CC0 Pfam markers of GTDB ar122.");
            Assert.That(sets, Has.Count.EqualTo(35), "One singleton collocated set per ar122 Pfam marker.");
            // PF00410 (Ribosomal_S8) GA1 = 24 bits, LENG = 125 (shared universal family).
            var s8 = hmms.Single(m => m.MarkerId == "PF00410");
            Assert.That(s8.BitScoreThreshold, Is.EqualTo(24.0).Within(1e-9),
                "PF00410 uses its Pfam GA1 gathering threshold (24 bits).");
        });
    }

    // B-TRUEPOSITIVE — the bundled bac120 PF01025 (GrpE) HMM detects E. coli GrpE, and NO other
    // bundled bac120 Pfam family does (specificity). Real-marker detection via the Plan7 engine.
    [Test]
    public void DetectMarkers_EcoliGrpE_HitsOnlyPF01025InBac120()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledBacterialMarkerHmms();

        var counts = MetagenomicsAnalyzer.DetectMarkers(new[] { EcoliGrpE }, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(counts["PF01025"], Is.EqualTo(1),
                "E. coli GrpE is the true positive for the GrpE family (PF01025) ⇒ exactly 1 hit.");
            foreach (var marker in hmms.Where(m => m.MarkerId != "PF01025"))
                Assert.That(counts[marker.MarkerId], Is.EqualTo(0),
                    $"GrpE must NOT match a different bac120 Pfam family ({marker.MarkerId}).");
        });
    }

    // B-COMPLETENESS — end-to-end over the 6 bac120 Pfam singleton sets, only PF01025 present
    // ⇒ exactly 1 of 6 sets recovered ⇒ Completeness = 100/6 %, Contamination = 0.
    [Test]
    public void EstimateBinQualityFromMarkers_BacterialBinWithOnlyGrpE_CompletenessIsOneOfSix()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledBacterialMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledBacterialMarkerSets();

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkers(new[] { EcoliGrpE }, sets, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(q.MarkerSetCount, Is.EqualTo(6), "Six bundled bac120 Pfam singleton sets.");
            Assert.That(q.Completeness, Is.EqualTo(100.0 / 6.0).Within(1e-10),
                "Only PF01025 present (1 of 6 singleton sets) ⇒ Completeness = 100/6 % per CheckM Eq.1.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10),
                "Single copy of the one detected marker ⇒ 0% contamination.");
            Assert.That(q.MarkersPresent, Is.EqualTo(1), "Exactly one bundled bac120 marker detected.");
        });
    }

    // A-TRUEPOSITIVE — E. coli uS8 hits the universal PF00410 in the archaeal set too (PF00410 is a
    // universal S8 family shared by bac120 and ar122). End-to-end completeness = 1 of 35 sets.
    [Test]
    public void EstimateBinQualityFromMarkers_ArchaealBinWithUniversaluS8_CompletenessIsOneOf35()
    {
        var hmms = MetagenomicsAnalyzer.LoadBundledArchaealMarkerHmms();
        var sets = MetagenomicsAnalyzer.BundledArchaealMarkerSets();

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkers(new[] { EcoliRpsH_uS8 }, sets, hmms);

        Assert.Multiple(() =>
        {
            Assert.That(q.MarkerSetCount, Is.EqualTo(35), "Thirty-five bundled ar122 Pfam singleton sets.");
            Assert.That(q.Completeness, Is.EqualTo(100.0 / 35.0).Within(1e-10),
                "Only the universal PF00410 present (1 of 35 singleton sets) ⇒ Completeness = 100/35 %.");
            Assert.That(q.Contamination, Is.EqualTo(0.0).Within(1e-10),
                "Single copy of the one detected marker ⇒ 0% contamination.");
            Assert.That(q.MarkersPresent, Is.EqualTo(1), "Exactly one bundled ar122 marker detected.");
        });
    }

    // B-DEFAULTS-UNCHANGED — the original 9-ribosomal accessor is untouched and disjoint in count
    // from the new domain sets (regression guard: existing API + defaults preserved).
    [Test]
    public void BundledMarkerSets_DomainSetsAreAdditive_RibosomalSetUnchanged()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MetagenomicsAnalyzer.LoadBundledRibosomalMarkerHmms(), Has.Count.EqualTo(9),
                "The original 9-marker ribosomal set is unchanged (additive opt-in expansion).");
            Assert.That(MetagenomicsAnalyzer.BundledRibosomalMarkerSets(), Has.Count.EqualTo(9),
                "The original 9 ribosomal singleton sets are unchanged.");
        });
    }

    #endregion
}
