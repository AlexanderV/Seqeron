namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Metagenomics area — CheckM-style genome-bin completeness/contamination from
/// single-copy marker genes (META-CHECKM-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, malformed and boundary inputs to a unit and asserts that the code NEVER
/// fails in an undisciplined way: no hang, no state corruption, and no *unhandled* runtime
/// exception (DivideByZeroException on an empty/zero-set marker collection,
/// IndexOutOfRangeException, KeyNotFoundException, NullReferenceException, OverflowException). Every
/// input must yield EITHER a well-defined, theory-correct result OR a *documented, intentional*
/// validation exception (ArgumentNullException from the entry-point guards; FormatException from the
/// HMMER3/f profile parser on malformed HMM text). It must also never emit *out-of-contract output*:
/// completeness ∈ [0, 100], contamination ≥ 0, neither NaN; both 0 (not NaN) when there is no
/// information. A passing "no crash" result must additionally be theory-correct — a sentinel that
/// returns 0/0 for everything would be a bug, so positive control cases pin the exact CheckM values.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-CHECKM-001 — CheckM marker-gene completeness / contamination (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 248.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — a marker id that is not in any marker set / not in the count map;
///          a marker count map whose keys disagree with the set ids; a malformed HMM profile text
///          fed to the Plan7 parser (bad version line, truncated body, non-amino alphabet, bad LENG).
///   • BE = Boundary Exploitation — EMPTY marker set (a set with 0 ids → would divide by |s| = 0);
///          NO marker sets at all (|M| = 0 → would divide by |M| = 0); NO genes detected
///          (empty count map / all counts ≤ 0); ALL markers duplicated (every N_g = 2 → maximal
///          single-duplication redundancy); extreme copy counts (int.MaxValue) on a singleton set.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes MC, BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The CheckM contract under test (derived independently from the doc, not read off the code)
/// ───────────────────────────────────────────────────────────────────────────
/// MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts implements the exact CheckM formula
/// (Parks et al. 2015, Genome Res 25:1043, Eqs. 1–2; reference impl MarkerSet.genomeCheck) over a
/// collection M of collocated single-copy marker sets and a per-marker copy-count map. With G_M the
/// markers found ≥ 1× and, for a marker found N ≥ 1 times, C_g = N − 1 (and C_g = 0 if absent):
///   • present_s   = |s ∩ G_M|                  (markers of set s found at least once)
///   • multiCopy_s = Σ_{g ∈ s} C_g              (extra/redundant copies within s)
///   • Completeness  = 100 · (1/|M|) · Σ_{s∈M} present_s/|s|        ∈ [0, 100]
///   • Contamination = 100 · (1/|M|) · Σ_{s∈M} multiCopy_s/|s|      ≥ 0 (unbounded above)
///   — docs/algorithms/Metagenomics/Genome_Binning.md §"CheckM-style marker-gene …";
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (EstimateBinQualityFromMarkerCounts ~1987, DetectMarkers ~2053,
///      EstimateBinQualityFromMarkers ~2090; records MarkerSet ~1849, MarkerHmm ~1874).
///
/// Boundary / malformed-input handling fixed by the doc + method XML, pinned here so it can never
/// silently drift:
///   • EMPTY marker set (ids.Count == 0) → contributes nothing and is SKIPPED before the /|s|
///     division (line ~2008), so it never divides by zero; it is also excluded from |M| (usedSets).
///   • |M| = 0 after skipping empties (or no sets at all) → "no information" → completeness AND
///     contamination are both DEFINED as 0.0 (NOT NaN, NOT a DivideByZeroException) via the
///     `usedSets > 0 ? … : 0.0` guard (lines ~2032–2033). This is the empty-marker-set div-by-zero
///     boundary the row calls out.
///   • NO genes detected (empty / all-≤0 count map) → present_s = 0, multiCopy_s = 0 for every set
///     → completeness 0%, contamination 0%, MarkersPresent 0.
///   • ALL markers duplicated (every present marker has N = 2) → present_s = |s| and
///     multiCopy_s = Σ(2−1) = |s| for every set → completeness 100%, contamination 100%.
///   • A marker id present in the count map but in NO set → ignored (sets define the expected
///     markers); a set id absent from the map → counted as missing (TryGetValue false ⇒ 0). No
///     KeyNotFoundException either way.
///   • Malformed HMM text → Plan7ProfileHmm.Parse throws a *documented* FormatException
///     (bad/missing "HMMER3/" version line, truncated body, non-amino ALPH, invalid LENG). This is
///     the only way a MarkerHmm can be "malformed"; EstimateBinQualityFromMarkers cannot receive a
///     half-built profile because construction fails first. No raw parse crash, no hang.
///   • Argument validation: null markerSets / markerCounts / proteins / markerHmms →
///     ArgumentNullException (documented entry-point guards).
///
/// LimitationPolicy: EstimateBinQualityFromMarkerCounts enforces "META-BIN-001" and so throws under
/// Strict. The test assembly's _LimitationPolicyTestBootstrap (ModuleInitializer) sets DefaultMode =
/// Permissive, so these tests exercise the real estimate rather than the guard (which has its own
/// dedicated tests). No per-test bootstrap is therefore required here.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Metagenomics")]
[Category("META-CHECKM-001")]
public class MetagenomicsCheckMFuzzTests
{
    private static MetagenomicsAnalyzer.MarkerSet Set(params string[] ids)
        => new(ids);

    // Assert the universal output contract: completeness ∈ [0,100], contamination ≥ 0, neither NaN.
    private static void AssertInContract(MetagenomicsAnalyzer.BinMarkerQuality q)
    {
        double.IsNaN(q.Completeness).Should().BeFalse("completeness must never be NaN");
        double.IsNaN(q.Contamination).Should().BeFalse("contamination must never be NaN");
        q.Completeness.Should().BeGreaterThanOrEqualTo(0.0).And.BeLessThanOrEqualTo(100.0,
            "CheckM completeness is a fraction of expected markers present ⇒ [0,100]%");
        q.Contamination.Should().BeGreaterThanOrEqualTo(0.0,
            "CheckM contamination (redundancy) is ≥ 0 (unbounded above with heavy duplication)");
        q.MarkerSetCount.Should().BeGreaterThanOrEqualTo(0);
        q.MarkerCount.Should().BeGreaterThanOrEqualTo(0);
        q.MarkersPresent.Should().BeGreaterThanOrEqualTo(0);
    }

    #region META-CHECKM-001 — empty / zero marker sets (BE: div-by-zero boundary)

    // BE — a marker SET with zero ids must be skipped, never divide by |s| = 0. Mixed with one real
    // singleton set {A} (A present once): the empty set contributes nothing, so |M| = 1, present 1/1.
    [Test]
    public void EmptyMarkerSet_IsSkipped_NoDivideByZero_OnlyNonEmptyContributes()
    {
        var sets = new[] { Set(), Set("A") };
        var counts = new Dictionary<string, int> { ["A"] = 1 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.MarkerSetCount.Should().Be(1, "the empty set is excluded from |M|");
        q.MarkerCount.Should().Be(1, "only the {A} set's single marker is counted");
        q.Completeness.Should().BeApproximately(100.0, 1e-10,
            "the one non-empty set {A} is fully recovered ⇒ 100% over |M| = 1");
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
    }

    // BE — EVERY marker set is empty ⇒ |M| (usedSets) = 0 ⇒ "no information" ⇒ 0/0, NOT NaN, NOT a
    // DivideByZeroException. This is the row's "empty marker set" div-by-zero boundary.
    [Test]
    public void AllMarkerSetsEmpty_GivesZeroZero_NotNaN_NoCrash()
    {
        var sets = new[] { Set(), Set(), Set() };
        var counts = new Dictionary<string, int> { ["A"] = 5 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.MarkerSetCount.Should().Be(0, "all sets empty ⇒ |M| = 0");
        q.Completeness.Should().Be(0.0, "no markers expected ⇒ completeness defined as 0");
        q.Contamination.Should().Be(0.0, "no markers expected ⇒ contamination defined as 0");
        q.MarkersPresent.Should().Be(0);
    }

    // BE — NO marker sets at all (empty M). Same |M| = 0 boundary, reached without any iteration.
    [Test]
    public void NoMarkerSets_GivesZeroZero_NotNaN()
    {
        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
            Array.Empty<MetagenomicsAnalyzer.MarkerSet>(),
            new Dictionary<string, int> { ["A"] = 2 });

        AssertInContract(q);
        q.MarkerSetCount.Should().Be(0);
        q.Completeness.Should().Be(0.0);
        q.Contamination.Should().Be(0.0);
    }

    #endregion

    #region META-CHECKM-001 — no genes detected (BE)

    // BE — empty count map (no genes detected at all) ⇒ present 0, multiCopy 0 ⇒ 0/0.
    [Test]
    public void NoGenesDetected_EmptyCountMap_GivesZeroCompletenessZeroContamination()
    {
        var sets = new[] { Set("A", "B"), Set("C", "D", "E") };
        var counts = new Dictionary<string, int>();

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(0.0, 1e-10);
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
        q.MarkerSetCount.Should().Be(2);
        q.MarkerCount.Should().Be(5, "|{A,B}| + |{C,D,E}| = 2 + 3");
        q.MarkersPresent.Should().Be(0);
    }

    // MC/BE — non-positive copy counts (0, negative) are treated as ABSENT, never as present and
    // never as negative redundancy. A whole map of ≤ 0 values is identical to "no genes detected".
    [Test]
    public void NonPositiveCounts_AreTreatedAsAbsent_NoNegativeContamination()
    {
        var sets = new[] { Set("A", "B", "C") };
        var counts = new Dictionary<string, int> { ["A"] = 0, ["B"] = -3, ["C"] = -1 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(0.0, 1e-10, "all counts ≤ 0 ⇒ nothing present");
        q.Contamination.Should().BeApproximately(0.0, 1e-10,
            "negative/zero counts must not produce negative or spurious redundancy");
        q.MarkersPresent.Should().Be(0);
    }

    #endregion

    #region META-CHECKM-001 — all markers duplicated (BE: maximal single-duplication redundancy)

    // BE — every marker present in exactly N = 2 copies. Hand-derived: for each set s,
    //   present_s = |s| (all present)  ⇒ present_s/|s| = 1
    //   multiCopy_s = Σ_{g∈s}(2−1) = |s| ⇒ multiCopy_s/|s| = 1
    // ⇒ Completeness = 100·(1/|M|)·Σ 1 = 100%,  Contamination = 100·(1/|M|)·Σ 1 = 100%.
    [Test]
    public void AllMarkersDuplicated_TwoCopies_Gives100Complete_100Contaminated()
    {
        var sets = new[] { Set("A", "B"), Set("C", "D", "E"), Set("F") };
        var counts = new Dictionary<string, int>
        {
            ["A"] = 2, ["B"] = 2, ["C"] = 2, ["D"] = 2, ["E"] = 2, ["F"] = 2,
        };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(100.0, 1e-10,
            "every marker present ⇒ each set fully recovered ⇒ 100%");
        q.Contamination.Should().BeApproximately(100.0, 1e-10,
            "every marker in 2 copies ⇒ multiCopy_s/|s| = 1 for every set ⇒ 100%");
        q.MarkersPresent.Should().Be(6, "all six distinct markers found");
    }

    // BE — extreme single-set duplication: a singleton set {A} with A in int.MaxValue copies.
    //   present 1/1 = 1 ⇒ completeness 100%;  multiCopy = (MaxValue−1)/1 ⇒ contamination = 100·(MaxValue−1).
    // Pins that contamination is unbounded-above-but-finite (no overflow to NaN/Inf, no negative wrap).
    [Test]
    public void ExtremeCopyCount_SingletonSet_FiniteUnboundedContamination_NoOverflow()
    {
        var sets = new[] { Set("A") };
        var counts = new Dictionary<string, int> { ["A"] = int.MaxValue };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(100.0, 1e-6);
        double expected = 100.0 * (int.MaxValue - 1);
        q.Contamination.Should().BeApproximately(expected, 1.0,
            "C_g = N − 1 for the duplicated marker; contamination = 100·(MaxValue−1) over |M| = 1");
        double.IsInfinity(q.Contamination).Should().BeFalse();
    }

    #endregion

    #region META-CHECKM-001 — marker / set id mismatch (MC)

    // MC — count map keys that are NOT in any set are ignored (the sets define the expected markers);
    // set ids not in the map are counted as absent. Hand-derived on set {A,B} with A=1, B absent,
    // plus garbage keys X=9, Y=4 that belong to no set:
    //   present 1/2 = 0.5 ⇒ completeness 50%;  multiCopy 0 ⇒ contamination 0%. No KeyNotFoundException.
    [Test]
    public void CountMapWithUnknownMarkerIds_AreIgnored_SetDefinesExpected()
    {
        var sets = new[] { Set("A", "B") };
        var counts = new Dictionary<string, int> { ["A"] = 1, ["X"] = 9, ["Y"] = 4 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(50.0, 1e-10,
            "only A∈{A,B} is present ⇒ 1/2; the stray X,Y are not expected markers");
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
        q.MarkersPresent.Should().Be(1, "only A is a present *expected* marker");
        q.MarkerCount.Should().Be(2, "MarkerCount counts set ids, unaffected by stray map keys");
    }

    // MC — randomized differential fuzz: build a random marker set and a random count map (with
    // overlapping, disjoint and ≤ 0 entries) and recompute the CheckM formula independently. The
    // implementation must agree exactly, stay in contract, and never throw on any of 200 inputs.
    [Test]
    public void RandomMarkerSetsAndCounts_NeverCrash_AlwaysInContract_MatchIndependentFormula()
    {
        var rng = new Random(248_2026); // locally seeded for determinism
        const string alphabet = "ABCDEFGHIJ";

        for (int iter = 0; iter < 200; iter++)
        {
            int setCount = rng.Next(0, 5);
            var sets = new List<MetagenomicsAnalyzer.MarkerSet>(setCount);
            for (int s = 0; s < setCount; s++)
            {
                int size = rng.Next(0, 5); // 0 ⇒ empty set (skipped)
                var ids = new List<string>(size);
                for (int j = 0; j < size; j++)
                    ids.Add(alphabet[rng.Next(alphabet.Length)].ToString());
                sets.Add(new MetagenomicsAnalyzer.MarkerSet(ids));
            }

            var counts = new Dictionary<string, int>(StringComparer.Ordinal);
            int entries = rng.Next(0, 12);
            for (int e = 0; e < entries; e++)
            {
                // Keys from the alphabet plus occasional out-of-set garbage; counts in [-2, 4].
                string key = rng.Next(4) == 0
                    ? "Z" + rng.Next(100) // out-of-set garbage id
                    : alphabet[rng.Next(alphabet.Length)].ToString();
                counts[key] = rng.Next(-2, 5);
            }

            MetagenomicsAnalyzer.BinMarkerQuality q = default;
            var act = () => q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);
            act.Should().NotThrow($"input must never crash (iter {iter})");
            AssertInContract(q);

            // Independent re-derivation of CheckM Eqs. 1–2 (counted distinct present per set).
            double comp = 0, cont = 0;
            int used = 0;
            var distinctPresent = new HashSet<string>(StringComparer.Ordinal);
            foreach (var set in sets)
            {
                if (set.MarkerIds.Count == 0) continue;
                int present = 0, multi = 0;
                foreach (var id in set.MarkerIds)
                {
                    int c = counts.TryGetValue(id, out var v) && v > 0 ? v : 0;
                    if (c >= 1) { present++; distinctPresent.Add(id); }
                    if (c > 1) multi += c - 1;
                }
                comp += (double)present / set.MarkerIds.Count;
                cont += (double)multi / set.MarkerIds.Count;
                used++;
            }
            double expComp = used > 0 ? 100.0 * comp / used : 0.0;
            double expCont = used > 0 ? 100.0 * cont / used : 0.0;

            q.Completeness.Should().BeApproximately(expComp, 1e-9, $"iter {iter}");
            q.Contamination.Should().BeApproximately(expCont, 1e-9, $"iter {iter}");
            q.MarkerSetCount.Should().Be(used, $"iter {iter}");
            q.MarkersPresent.Should().Be(distinctPresent.Count, $"iter {iter}");
        }
    }

    #endregion

    #region META-CHECKM-001 — null-argument guards (documented validation exceptions)

    [Test]
    public void NullMarkerSets_ThrowsArgumentNullException()
    {
        var act = () => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
            null!, new Dictionary<string, int>());
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void NullMarkerCounts_ThrowsArgumentNullException()
    {
        var act = () => MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(
            Array.Empty<MetagenomicsAnalyzer.MarkerSet>(), null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DetectMarkers_NullProteinsOrHmms_ThrowsArgumentNullException()
    {
        var a = () => MetagenomicsAnalyzer.DetectMarkers(null!, Array.Empty<MetagenomicsAnalyzer.MarkerHmm>());
        a.Should().Throw<ArgumentNullException>();
        var b = () => MetagenomicsAnalyzer.DetectMarkers(Array.Empty<string>(), null!);
        b.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region META-CHECKM-001 — malformed HMM (MC: Plan7 profile parser)

    // MC — a MarkerHmm can only be "malformed" by feeding malformed HMMER3/f text to the Plan7
    // parser. Every malformed profile must throw a *documented* FormatException (never a raw crash,
    // never a hang), so a half-built MarkerHmm can never reach DetectMarkers/EstimateBinQualityFromMarkers.
    [TestCase("", TestName = "MalformedHmm_Empty")]
    [TestCase("not a hmm file at all", TestName = "MalformedHmm_NoVersionLine")]
    [TestCase("HMMER3/f [3.4]\n", TestName = "MalformedHmm_VersionOnly_Truncated")]
    [TestCase("HMMER3/f [3.4]\nNAME  X\nLENG  3\nALPH  amino\n", TestName = "MalformedHmm_NoMainModel")]
    [TestCase("HMMER3/f [3.4]\nNAME  X\nLENG  3\nALPH  DNA\nHMM\n", TestName = "MalformedHmm_NonAminoAlphabet")]
    [TestCase("HMMER3/f [3.4]\nNAME  X\nLENG  0\nALPH  amino\nHMM\n", TestName = "MalformedHmm_ZeroLength")]
    [TestCase("HMMER3/f [3.4]\nNAME  X\nLENG  -1\nALPH  amino\nHMM\n", TestName = "MalformedHmm_NegativeLength")]
    public void MalformedHmmText_ThrowsFormatException_NeverRawCrashOrHang(string hmmText)
    {
        var act = () => Plan7ProfileHmm.Parse(hmmText);
        act.Should().Throw<FormatException>(
            "a malformed HMMER3/f profile is a documented validation failure, not a runtime crash");
    }

    // MC — feeding malformed HMM text through the marker loader (LoadMarkerHmms) surfaces the SAME
    // documented FormatException at construction time — the bad profile cannot become a MarkerHmm.
    [Test]
    public void LoadMarkerHmms_FromMalformedReader_ThrowsFormatException()
    {
        var readers = new TextReader[] { new StringReader("this is not a HMMER3 profile") };
        var act = () => MetagenomicsAnalyzer.LoadMarkerHmms(readers);
        act.Should().Throw<FormatException>();
    }

    // BE — DetectMarkers with NO HMMs and/or NO proteins detects nothing: an empty count map, which
    // feeds EstimateBinQualityFromMarkers to a clean 0%/0% (no genes), with no crash.
    [Test]
    public void EstimateBinQualityFromMarkers_NoHmmsNoProteins_GivesZeroQuality()
    {
        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkers(
            proteins: Array.Empty<string>(),
            markerSets: new[] { Set("PF00318"), Set("PF00177") },
            markerHmms: Array.Empty<MetagenomicsAnalyzer.MarkerHmm>());

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(0.0, 1e-10, "no proteins / no HMMs ⇒ nothing detected");
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
        q.MarkersPresent.Should().Be(0);
        q.MarkerSetCount.Should().Be(2);
    }

    #endregion

    #region META-CHECKM-001 — positive control (a 0/0 sentinel must fail this)

    // A bin recovering 2 of 4 expected single-copy markers (each once) over two 2-marker sets:
    //   s1 = {A,B}: A=1,B=0 ⇒ present 1/2 ;  s2 = {C,D}: C=1,D=0 ⇒ present 1/2.
    //   completeness = 100·(0.5+0.5)/2 = 50% ; contamination = 0%.
    [Test]
    public void HalfMarkersPresentSingleCopy_Gives50Complete_0Contaminated()
    {
        var sets = new[] { Set("A", "B"), Set("C", "D") };
        var counts = new Dictionary<string, int> { ["A"] = 1, ["C"] = 1 };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);

        AssertInContract(q);
        q.Completeness.Should().BeApproximately(50.0, 1e-10,
            "2 of 4 expected single-copy markers present ⇒ 50% complete");
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
        q.MarkersPresent.Should().Be(2);
    }

    #endregion
}
