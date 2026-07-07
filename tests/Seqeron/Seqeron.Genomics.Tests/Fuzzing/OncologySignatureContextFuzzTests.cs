using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology SBS-96 trinucleotide-context catalog — ONCO-SIG-001.
/// The unit under test is the deterministic, specification-driven single-base-substitution (SBS)
/// context classifier and catalogue builder
/// <see cref="OncologyAnalyzer.ClassifySbsContext"/> (one SBS → its 96-channel label),
/// <see cref="OncologyAnalyzer.EnumerateSbs96Channels"/> (the 96 canonical labels) and
/// <see cref="OncologyAnalyzer.Build96ContextCatalog"/> (tally SBS variants into the 96-channel
/// spectrum), implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime fault
/// (KeyNotFoundException / IndexOutOfRangeException / ArgumentOutOfRange on a
/// substring / DivideByZero). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a *documented, intentional* outcome — here, an
/// <see cref="ArgumentException"/> for any base that is not A/C/G/T or for a
/// non-mutation (ref == alt), and an <see cref="ArgumentNullException"/> for a
/// null variant enumerable.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SIG-001 — SBS-96 trinucleotide-context catalog (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 96.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content — невалідний контент.
///     Targets (checklist row 96): "non-SNV variant, ambiguous base, no flanking
///     context, empty set". Mapped onto this char-per-base, single-base-only API:
///       — NON-SNV VARIANT: indels / MNVs / doublet substitutions are *out of
///         scope* (the unit is single-base only, §3.3, §6.2). A caller expressing
///         such a variant supplies a non-ACGT placeholder for the gap / extra base
///         (e.g. '-', '*', a length>1 intent), which is rejected with a documented
///         ArgumentException — never silently miscounted into the SBS catalogue and
///         never a crash.
///       — AMBIGUOUS BASE: an IUPAC ambiguity code (N, R, Y, S, W, K, M, B, D, H,
///         V) anywhere in ref / alt / 5' / 3' has no defined trinucleotide context
///         ⇒ documented ArgumentException (§6.1: "Non-ACGT base ⇒ ArgumentException").
///         It must NOT raise a KeyNotFoundException (catalogue indexing) nor leak.
///       — NO FLANKING CONTEXT: a mutation at the very start / end of a contig has
///         a missing 5' or 3' flank. The caller signals "no flank" with a non-ACGT
///         sentinel (NUL '\0', '.', ' ', 'X'); that is guarded by the same A/C/G/T
///         validation ⇒ ArgumentException, with NO IndexOutOfRange / substring
///         crash (classification is a constant-time base computation — no substring
///         is taken, §5.2).
///       — EMPTY SET: Build96ContextCatalog over no variants ⇒ all 96 channels
///         present, every count 0 (§6.1 "Empty variant collection"), no
///         DivideByZero, fixed 96-dimensional shape.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// SBS96_Trinucleotide_Context_Catalog.md (docs/algorithms/Oncology/SBS96_Trinucleotide_Context_Catalog.md):
///   • Six pyrimidine substitutions: C>A, C>G, C>T, T>A, T>C, T>G.            (§2.2, §4.2)
///   • Channel label `5'[REF>ALT]3'`; mutated base centred; 96 = 6×4×4.       (§2.2, §2.4 INV-02)
///   • Pyrimidine-strand folding: a PURINE reference (A/G) is reverse-
///       complemented — REF'=comp(REF), ALT'=comp(ALT), 5''=comp(3'),
///       3''=comp(5'); complement map A↔T, C↔G. So G>T at 5'-T G A-3' ⇒
///       `T[C>A]A`.                                                            (§2.2, §4.1, INV-01, INV-04)
///   • Build96ContextCatalog: all 96 channels present (zeros included);
///       Σ counts = number of classifiable variants (a partition).            (§3.2, INV-03)
///   • Bases upper-cased; non-ACGT ⇒ ArgumentException; ref==alt ⇒
///       ArgumentException; null variants ⇒ ArgumentNullException.            (§3.3, §6.1)
///   • Single-base substitutions ONLY — indels / DBS / MNV are not handled.   (§3.3, §6.2)
///   • Classification is constant-time base arithmetic — NO substring search,
///       so the "no flank" case is a base-validation guard, not a slice.      (§5.2)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySignatureContextFuzzTests
{
    private static readonly char[] AcgtBases = { 'A', 'C', 'G', 'T' };

    // The six documented pyrimidine substitutions (§2.2).
    private static readonly (char Ref, char Alt)[] PyrimidineSubs =
    {
        ('C', 'A'), ('C', 'G'), ('C', 'T'), ('T', 'A'), ('T', 'C'), ('T', 'G'),
    };

    // IUPAC ambiguity codes (everything that is NOT one of the four unambiguous DNA bases).
    // Source: IUPAC nucleotide notation (R,Y,S,W,K,M,B,D,H,V,N).
    private static readonly char[] AmbiguityCodes =
    {
        'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'n', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
    };

    // Sentinels a caller would use to signal a MISSING flank (start/end of contig) or a
    // NON-SNV gap (indel / multi-base) in this char-per-base API.
    private static readonly char[] MissingOrGapSentinels =
    {
        '\0', '.', ' ', '-', '*', 'X', '?', '\t', '0', '1', 'U', 'I',
    };

    private static char Complement(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"'{b}' is not a DNA base"),
    };

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY successfully classified
    // channel label: it is exactly `5'[REF>ALT]3'`, the centre REF is a PYRIMIDINE
    // (C or T, INV-01), (REF,ALT) is one of the six pyrimidine substitutions, and
    // the two flanks are A/C/G/T. This is what stops a fuzz test from rubber-
    // stamping a malformed / non-folded label green.
    private static void AssertWellFormedChannel(string channel)
    {
        channel.Should().HaveLength(7, "a channel is exactly 5'[REF>ALT]3' (7 chars)");
        channel[1].Should().Be('[');
        channel[3].Should().Be('>');
        channel[5].Should().Be(']');

        char five = channel[0];
        char reference = channel[2];
        char alternate = channel[4];
        char three = channel[6];

        AcgtBases.Should().Contain(five, "the 5' flank must be A/C/G/T");
        AcgtBases.Should().Contain(three, "the 3' flank must be A/C/G/T");
        reference.Should().BeOneOf(new[] { 'C', 'T' },
            "the centre reference base is folded to a pyrimidine (INV-01)");
        PyrimidineSubs.Should().Contain((reference, alternate),
            "(REF,ALT) must be one of the six pyrimidine substitutions (§2.2)");
    }

    #region ONCO-SIG-001 — positive sanity (exact documented classes; pyrimidine-strand folding)

    // The headline documented worked example (§7.1): a G>T at 5'-T G A-3' folds onto the pyrimidine
    // strand to EXACTLY T[C>A]A (reverse-complement of TGA is TCA; G>T becomes C>A).
    [Test]
    public void ClassifySbsContext_PurineRefWorkedExample_FoldsTo_TCA()
    {
        ClassifySbsContext('T', 'G', 'T', 'A').Should().Be("T[C>A]A",
            "G>T at 5'-T G A-3' reverse-complements to the pyrimidine-strand channel T[C>A]A (§7.1, INV-04)");
    }

    // A pyrimidine-reference SBS is kept as-is: a C>T at context A_G ⇒ exactly A[C>T]G.
    [Test]
    public void ClassifySbsContext_PyrimidineRef_KeptAsIs_ACG()
    {
        ClassifySbsContext('A', 'C', 'T', 'G').Should().Be("A[C>T]G",
            "a pyrimidine (C) reference is already on the pyrimidine strand — unchanged (§2.2)");
    }

    // Strand-equivalence (INV-04): a purine-ref mutation and its reverse-complement pyrimidine form map to
    // the SAME channel. G>A at 5'-A G C-3' folds (RC of AGC = GCT; G>A ⇒ C>T) to G[C>T]T, which equals the
    // direct pyrimidine classification of C>T at 5'-G C T-3'.
    [Test]
    public void ClassifySbsContext_PurineAndPyrimidineForms_SameChannel()
    {
        string folded = ClassifySbsContext('A', 'G', 'A', 'C');   // purine ref G>A
        string direct = ClassifySbsContext('G', 'C', 'T', 'T');   // its pyrimidine-strand twin C>T

        folded.Should().Be("G[C>T]T");
        direct.Should().Be("G[C>T]T");
        folded.Should().Be(direct, "a substitution and its reverse complement are biologically one channel (INV-04)");
    }

    // Case-insensitivity (§3.3): lower-case bases classify identically to upper-case.
    [Test]
    public void ClassifySbsContext_LowercaseBases_ClassifiedSameAsUppercase()
    {
        ClassifySbsContext('a', 'c', 't', 'g').Should().Be(ClassifySbsContext('A', 'C', 'T', 'G'),
            "bases are upper-cased before classification (§3.3)");
        ClassifySbsContext('t', 'g', 't', 'a').Should().Be("T[C>A]A");
    }

    // EnumerateSbs96Channels is the exact 96-label channel space: 96 distinct labels, all well-formed,
    // and equal to the cross-product of the six pyrimidine subs × 4 × 4 (INV-02).
    [Test]
    public void EnumerateSbs96Channels_Is_The_Exact_96_Channel_Space()
    {
        var channels = EnumerateSbs96Channels();

        channels.Should().HaveCount(96, "6 substitutions × 4 5'-bases × 4 3'-bases = 96 (INV-02)");
        channels.Should().OnlyHaveUniqueItems("the 96 channels are distinct labels");

        var expected = new List<string>();
        foreach (var (r, a) in PyrimidineSubs)
        {
            foreach (char five in AcgtBases)
            {
                foreach (char three in AcgtBases)
                {
                    expected.Add($"{five}[{r}>{a}]{three}");
                }
            }
        }

        channels.Should().Equal(expected, "channels are in substitution-major then 5' then 3' order (§4.2)");
        foreach (string c in channels)
        {
            AssertWellFormedChannel(c);
        }
    }

    // Every classifiable input lands in a channel that EnumerateSbs96Channels declares — closure of the
    // classifier over the 96-channel space (a sanity bridge between the two public methods).
    [Test]
    public void ClassifySbsContext_AllAcgtInputs_LandInEnumeratedChannelSpace()
    {
        var channelSet = EnumerateSbs96Channels().ToHashSet(StringComparer.Ordinal);
        var produced = new HashSet<string>(StringComparer.Ordinal);

        foreach (char five in AcgtBases)
        {
            foreach (char reference in AcgtBases)
            {
                foreach (char alternate in AcgtBases)
                {
                    if (reference == alternate)
                    {
                        continue;
                    }

                    foreach (char three in AcgtBases)
                    {
                        string channel = ClassifySbsContext(five, reference, alternate, three);
                        AssertWellFormedChannel(channel);
                        channelSet.Should().Contain(channel, "every classification lands in the enumerated space");
                        produced.Add(channel);
                    }
                }
            }
        }

        // 4×4×3×4 = 192 ordered (5',ref,alt,3') inputs collapse onto all 96 channels exactly (each channel
        // is reached by exactly two strand-equivalent inputs).
        produced.Should().HaveCount(96, "the 192 ACGT inputs fold onto all 96 channels (each reached twice)");
    }

    #endregion

    #region ONCO-SIG-001 — MC: ambiguous base (IUPAC code ⇒ ArgumentException, never KeyNotFound)

    // A single ambiguous base in ANY of the four positions has no defined trinucleotide context ⇒
    // ArgumentException (§6.1). Critically it must NOT be a KeyNotFoundException (catalogue indexing)
    // nor an IndexOutOfRange — the guard is base validation, before any dictionary lookup.
    [Test]
    public void ClassifySbsContext_AmbiguousBase_InAnyPosition_ThrowsArgumentException_NeverKeyNotFound()
    {
        foreach (char amb in AmbiguityCodes)
        {
            // 5' flank ambiguous
            AssertClassifyRejected(amb, 'C', 'T', 'G', amb, "5'");
            // reference ambiguous
            AssertClassifyRejected('A', amb, 'T', 'G', amb, "ref");
            // alternate ambiguous
            AssertClassifyRejected('A', 'C', amb, 'G', amb, "alt");
            // 3' flank ambiguous
            AssertClassifyRejected('A', 'C', 'T', amb, amb, "3'");
        }
    }

    // The most common ambiguity code, N (unknown base), in every position.
    [Test]
    public void ClassifySbsContext_NBase_EveryPosition_ThrowsArgumentException()
    {
        AssertClassifyRejected('N', 'C', 'T', 'G', 'N', "5'");
        AssertClassifyRejected('A', 'N', 'T', 'G', 'N', "ref");
        AssertClassifyRejected('A', 'C', 'N', 'G', 'N', "alt");
        AssertClassifyRejected('A', 'C', 'T', 'N', 'N', "3'");
    }

    // Build96ContextCatalog must surface the SAME ArgumentException for an ambiguous variant — never a
    // KeyNotFoundException from catalog[channel]++ on an unclassifiable variant, and it must NOT silently
    // miscount the variant. (Validation precedes the dictionary increment.)
    [Test]
    public void Build96ContextCatalog_AmbiguousVariant_ThrowsArgumentException_NotKeyNotFound()
    {
        var variants = new[] { ('A', 'C', 'T', 'G'), ('A', 'N', 'T', 'G') }; // second is ambiguous

        // Throw<ArgumentException> already excludes KeyNotFoundException (not an ArgumentException),
        // so a KeyNotFound from catalog[channel]++ on an unclassifiable variant would fail this assertion.
        FluentActions.Invoking(() => Build96ContextCatalog(variants))
            .Should().Throw<ArgumentException>(
                "an ambiguous base must be rejected by validation before any catalogue indexing (§6.1)");
    }

    // Fuzz: random placement of a random IUPAC code among otherwise-valid ACGT bases — always an
    // ArgumentException, never a KeyNotFound / IndexOutOfRange / unhandled fault.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifySbsContext_RandomAmbiguousFuzz_AlwaysArgumentException_NeverCrash()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            char[] pos =
            {
                AcgtBases[rng.Next(4)], AcgtBases[rng.Next(4)], AcgtBases[rng.Next(4)], AcgtBases[rng.Next(4)],
            };

            // Make ref ≠ alt so the ONLY defect is the injected ambiguity (isolates the cause).
            if (pos[1] == pos[2])
            {
                pos[2] = AcgtBases[(Array.IndexOf(AcgtBases, pos[1]) + 1) % 4];
            }

            int corrupt = rng.Next(4);
            pos[corrupt] = AmbiguityCodes[rng.Next(AmbiguityCodes.Length)];

            // Throw<ArgumentException> excludes KeyNotFoundException / IndexOutOfRangeException (neither is
            // an ArgumentException) — so any such crash would fail this assertion.
            FluentActions.Invoking(() => ClassifySbsContext(pos[0], pos[1], pos[2], pos[3]))
                .Should().Throw<ArgumentException>(
                    $"seed {seed}: ambiguous base at {corrupt} must be rejected, not crash");
        }
    }

    #endregion

    #region ONCO-SIG-001 — MC: no flanking context (missing 5'/3' ⇒ guarded, no IndexOutOfRange/substring crash)

    // A mutation at the very start of a contig has no 5' flank; at the very end, no 3' flank. The caller
    // signals "no flank" with a non-ACGT sentinel ('\0', '.', ' ', 'X', …). The unit guards it with the
    // A/C/G/T validation ⇒ ArgumentException, and — because classification is constant-time base arithmetic
    // with NO substring (§5.2) — NEVER an IndexOutOfRangeException / substring ArgumentOutOfRange.
    [Test]
    public void ClassifySbsContext_MissingFlankSentinel_ThrowsArgumentException_NoIndexOutOfRange()
    {
        foreach (char sentinel in MissingOrGapSentinels)
        {
            // Missing 5' flank (start of contig): a real C>T but no upstream base.
            AssertNoIndexCrash(() => ClassifySbsContext(sentinel, 'C', 'T', 'G'), sentinel, "5' flank");
            // Missing 3' flank (end of contig): a real C>T but no downstream base.
            AssertNoIndexCrash(() => ClassifySbsContext('A', 'C', 'T', sentinel), sentinel, "3' flank");
            // Both flanks missing (single-base contig).
            AssertNoIndexCrash(() => ClassifySbsContext(sentinel, 'C', 'T', sentinel), sentinel, "both flanks");
        }
    }

    // Even when BOTH the substitution and the flanks would be valid except the flank is the NUL sentinel,
    // the documented guard fires — the result is an ArgumentException, never a partial/garbled label.
    [Test]
    public void ClassifySbsContext_NulFlank_PurineRef_StillGuarded_NoFold_NoCrash()
    {
        // A purine-ref G>A whose 5' flank is missing: the fold path reads the 3' neighbour as the new 5'
        // — must still validate the missing flank up front and throw, never IndexOutOfRange in the fold.
        FluentActions.Invoking(() => ClassifySbsContext('\0', 'G', 'A', 'C'))
            .Should().Throw<ArgumentException>("a missing flank is non-ACGT ⇒ guarded before folding (§6.1)");
    }

    // Fuzz: random insertion of a missing-flank sentinel into the 5' and/or 3' position — always a guarded
    // ArgumentException, never an IndexOutOfRange / substring fault, across random valid substitutions.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifySbsContext_RandomMissingFlankFuzz_AlwaysGuarded_NoIndexCrash()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            var (reference, alternate) = PyrimidineSubs[rng.Next(PyrimidineSubs.Length)];
            // Randomly also use the purine-strand twin to exercise the folding code path.
            if (rng.Next(2) == 0)
            {
                reference = Complement(reference);
                alternate = Complement(alternate);
            }

            char five = rng.Next(2) == 0
                ? MissingOrGapSentinels[rng.Next(MissingOrGapSentinels.Length)]
                : AcgtBases[rng.Next(4)];
            char three = rng.Next(2) == 0
                ? MissingOrGapSentinels[rng.Next(MissingOrGapSentinels.Length)]
                : AcgtBases[rng.Next(4)];

            // Ensure at least one flank is a missing sentinel.
            bool fiveMissing = Array.IndexOf(MissingOrGapSentinels, five) >= 0;
            bool threeMissing = Array.IndexOf(MissingOrGapSentinels, three) >= 0;
            if (!fiveMissing && !threeMissing)
            {
                five = MissingOrGapSentinels[rng.Next(MissingOrGapSentinels.Length)];
            }

            char r = reference;
            char a = alternate;
            char f = five;
            char t = three;
            AssertNoIndexCrash(() => ClassifySbsContext(f, r, a, t), f, $"seed {seed}");
        }
    }

    #endregion

    #region ONCO-SIG-001 — MC: non-SNV variant (indel/MNV out of scope ⇒ guarded, never miscounted)

    // Single-base substitutions ONLY (§3.3, §6.2). A non-SNV variant — an indel or multi-base substitution
    // — is expressed in this char-per-base API by a gap / placeholder character ('-', '*', '.', '0', …) in
    // the reference or alternate slot. It must be rejected with ArgumentException — NEVER folded into a
    // false SBS channel, NEVER crash.
    [Test]
    public void ClassifySbsContext_NonSnvGapPlaceholder_InRefOrAlt_ThrowsArgumentException()
    {
        foreach (char gap in MissingOrGapSentinels)
        {
            // Deletion-style: alternate is a gap (REF deleted) — not a substitution.
            AssertNoIndexCrash(() => ClassifySbsContext('A', 'C', gap, 'G'), gap, "alt gap (deletion)");
            // Insertion-style / MNV: reference is a gap / placeholder.
            AssertNoIndexCrash(() => ClassifySbsContext('A', gap, 'T', 'G'), gap, "ref gap (insertion)");
        }
    }

    // ref == alt is "not a substitution" (§3.3, §6.1) — a degenerate non-SNV (no base change). Every
    // ACGT base equal to itself must raise ArgumentException, distinctly from the non-ACGT guard.
    [Test]
    public void ClassifySbsContext_RefEqualsAlt_ThrowsArgumentException_NotAMutation()
    {
        foreach (char b in AcgtBases)
        {
            FluentActions.Invoking(() => ClassifySbsContext('A', b, b, 'G'))
                .Should().Throw<ArgumentException>($"ref == alt ('{b}') is not a substitution (§6.1)");
        }
    }

    // Build96ContextCatalog with a non-SNV (gap) variant must throw — and must NOT have miscounted any
    // VALID variant that preceded it: the partition invariant (INV-03) is never violated by silent
    // acceptance of a malformed variant. (Throwing before incrementing keeps the catalogue honest.)
    [Test]
    public void Build96ContextCatalog_NonSnvVariant_ThrowsArgumentException_NeverMiscounts()
    {
        var variants = new[] { ('A', 'C', 'T', 'G'), ('A', 'C', '-', 'G') }; // second is a deletion-style gap

        FluentActions.Invoking(() => Build96ContextCatalog(variants))
            .Should().Throw<ArgumentException>("a non-SNV (gap) variant is out of scope for SBS-96 (§6.2)");
    }

    // Fuzz: random gap/placeholder injected into ref or alt of an otherwise-valid SBS — always rejected,
    // never miscounted, never crash.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifySbsContext_RandomNonSnvFuzz_AlwaysArgumentException()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            char five = AcgtBases[rng.Next(4)];
            char three = AcgtBases[rng.Next(4)];
            char reference = AcgtBases[rng.Next(4)];
            char alternate = AcgtBases[rng.Next(4)];

            // Corrupt exactly ref OR alt with a gap so the only defect is the non-SNV placeholder.
            char gap = MissingOrGapSentinels[rng.Next(MissingOrGapSentinels.Length)];
            bool corruptRef = rng.Next(2) == 0;
            if (corruptRef)
            {
                reference = gap;
            }
            else
            {
                alternate = gap;
            }

            FluentActions.Invoking(() => ClassifySbsContext(five, reference, alternate, three))
                .Should().Throw<ArgumentException>(
                    $"seed {seed}: non-SNV gap in {(corruptRef ? "ref" : "alt")} must be rejected");
        }
    }

    #endregion

    #region ONCO-SIG-001 — MC: empty set (empty catalogue of the correct shape, no DivideByZero)

    // The empty catalogue: no variants ⇒ all 96 channels present, every count exactly 0 (§6.1).
    // Fixed 96-dimensional shape, no DivideByZero, no missing keys.
    [Test]
    public void Build96ContextCatalog_EmptySet_AllNinetySixChannelsZero()
    {
        var catalog = Build96ContextCatalog(Enumerable.Empty<(char, char, char, char)>());

        catalog.Should().HaveCount(96, "the empty spectrum still has the full fixed 96-channel shape (§5.2)");
        catalog.Keys.Should().BeEquivalentTo(EnumerateSbs96Channels(),
            "every one of the 96 canonical channels is present even when empty");
        catalog.Values.Should().OnlyContain(v => v == 0, "no variants ⇒ a partition of the empty set ⇒ all zero");
        catalog.Values.Sum().Should().Be(0, "Σ counts = number of classifiable variants = 0 (INV-03)");
    }

    // An empty array (distinct from an empty enumerable) yields the same all-zero 96-channel spectrum.
    [Test]
    public void Build96ContextCatalog_EmptyArray_AllZero_SameShape()
    {
        var catalog = Build96ContextCatalog(Array.Empty<(char, char, char, char)>());

        catalog.Should().HaveCount(96);
        catalog.Values.Sum().Should().Be(0);
    }

    // Null variants ⇒ documented ArgumentNullException (§3.3) — distinct from the empty case.
    [Test]
    public void Build96ContextCatalog_NullVariants_ThrowsArgumentNull()
    {
        FluentActions.Invoking(() => Build96ContextCatalog(null!))
            .Should().Throw<ArgumentNullException>("a null variant enumerable is a documented guard (§3.3)");
    }

    // Partition invariant on a well-formed multiset (the positive counterpart to the empty case): Σ counts
    // equals the number of variants, all 96 keys present, and a known mutation lands in its exact channel.
    [Test]
    public void Build96ContextCatalog_WellFormedMultiset_PartitionsExactly()
    {
        var variants = new[]
        {
            ('A', 'C', 'A', 'A'),  // A[C>A]A  (pyrimidine, unchanged)
            ('T', 'G', 'T', 'A'),  // folds to T[C>A]A
            ('A', 'C', 'A', 'A'),  // A[C>A]A again
        };

        var catalog = Build96ContextCatalog(variants);

        catalog.Should().HaveCount(96, "fixed shape regardless of input (§5.2)");
        catalog["A[C>A]A"].Should().Be(2, "two A[C>A]A variants");
        catalog["T[C>A]A"].Should().Be(1, "the purine-ref G>T folds to T[C>A]A (§7.1)");
        catalog.Values.Sum().Should().Be(variants.Length, "Σ counts = number of variants (INV-03)");
        catalog.Values.Where(v => v > 0).Should().HaveCount(2, "only the two reached channels are non-zero");
    }

    // Fuzz: random well-formed multisets always partition exactly (Σ counts = n, fixed 96 shape) and the
    // empty draw (n = 0) is the all-zero spectrum — no DivideByZero across the size range.
    [Test]
    [CancelAfter(20_000)]
    public void Build96ContextCatalog_RandomWellFormedSets_PartitionInvariant_IncludingEmpty()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 40); // includes n = 0 (the empty set)
            var variants = new List<(char, char, char, char)>(n);
            for (int i = 0; i < n; i++)
            {
                char five = AcgtBases[rng.Next(4)];
                char three = AcgtBases[rng.Next(4)];
                char reference = AcgtBases[rng.Next(4)];
                char alternate = AcgtBases[rng.Next(4)];
                if (reference == alternate)
                {
                    alternate = AcgtBases[(Array.IndexOf(AcgtBases, reference) + 1) % 4];
                }

                variants.Add((five, reference, alternate, three));
            }

            IReadOnlyDictionary<string, int> catalog = null!;
            FluentActions.Invoking(() => catalog = Build96ContextCatalog(variants))
                .Should().NotThrow($"seed {seed}: a well-formed multiset of size {n} must not crash");

            catalog.Should().HaveCount(96, $"seed {seed}: fixed 96-channel shape");
            catalog.Values.Sum().Should().Be(n, $"seed {seed}: Σ counts = n = {n} (INV-03)");
            catalog.Values.Should().OnlyContain(v => v >= 0, $"seed {seed}: counts are non-negative");
            foreach (string key in catalog.Keys)
            {
                AssertWellFormedChannel(key);
            }
        }
    }

    #endregion

    #region helpers

    private static void AssertClassifyRejected(char five, char reference, char alternate, char three,
        char offending, string where)
    {
        // Throw<ArgumentException> excludes KeyNotFoundException (not an ArgumentException): the base guard
        // fires before any catalogue indexing, so a KeyNotFound would fail this assertion.
        FluentActions.Invoking(() => ClassifySbsContext(five, reference, alternate, three))
            .Should().Throw<ArgumentException>(
                $"ambiguous/non-ACGT base '{offending}' in {where} ⇒ ArgumentException");
    }

    // Asserts the call throws a *documented* ArgumentException (the base-validation guard) and explicitly
    // NOT an IndexOutOfRangeException / KeyNotFoundException (the crash modes the "no flank" / non-SNV
    // targets are guarding against).
    private static void AssertNoIndexCrash(Action act, char offending, string where)
    {
        Exception? captured = null;
        try
        {
            act();
        }
        catch (Exception ex)
        {
            captured = ex;
        }

        captured.Should().NotBeNull($"a malformed base '{(int)offending}' in {where} must be guarded, not pass silently");
        captured.Should().BeOfType<ArgumentException>(
            $"'{(int)offending}' in {where} ⇒ documented ArgumentException (§6.1), never an undisciplined crash");
        captured.Should().NotBeOfType<IndexOutOfRangeException>(
            "classification takes no substring (§5.2) — a missing flank never causes IndexOutOfRange");
        captured.Should().NotBeOfType<ArgumentOutOfRangeException>(
            "no substring/slice is performed — never an ArgumentOutOfRange from indexing");
        captured.Should().NotBeOfType<KeyNotFoundException>(
            "the base guard precedes catalogue indexing — never KeyNotFound");
    }

    #endregion
}
