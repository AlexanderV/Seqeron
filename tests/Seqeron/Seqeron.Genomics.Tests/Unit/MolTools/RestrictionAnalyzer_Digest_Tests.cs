// Canonical Test File: RESTR-DIGEST-001 - Restriction Digest Simulation
// Methods Under Test: RestrictionAnalyzer.Digest, GetDigestSummary, CreateMap, AreCompatible, FindCompatibleEnzymes
//
// Evidence Sources:
// - Wikipedia: Restriction digest (https://en.wikipedia.org/wiki/Restriction_digest)
// - Wikipedia: Restriction enzyme (https://en.wikipedia.org/wiki/Restriction_enzyme)
// - Wikipedia: Restriction map (https://en.wikipedia.org/wiki/Restriction_map)
// - Addgene: Restriction Digest Protocol (https://www.addgene.org/protocols/restriction-digest/)
// - Roberts RJ (1976) Restriction endonucleases
// - REBASE: The Restriction Enzyme Database
//
// Algorithm Documentation: docs/algorithms/MolTools/Restriction_Digest_Simulation.md
// TestSpec: TestSpecs/RESTR-DIGEST-001.md

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Canonical test class for RESTR-DIGEST-001: Restriction Digest Simulation.
/// Tests Digest, GetDigestSummary, CreateMap, AreCompatible, and FindCompatibleEnzymes methods.
/// For FindSites/GetEnzyme tests, see RestrictionAnalyzer_FindSites_Tests.cs (RESTR-FIND-001).
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Digest_Tests
{
    #region Digest Core Tests

    /// <summary>
    /// Evidence: Wikipedia - Single cut produces two fragments.
    /// Invariant #1: k cut positions → k+1 fragments.
    /// Invariant #2: Fragment sum = original sequence length (Addgene Protocol).
    /// </summary>
    [Test]
    public void Digest_SingleCut_ReturnsTwoFragmentsWithCorrectSum()
    {
        // Arrange: EcoRI (G↓AATTC) at position 3, cuts at position 4
        var sequence = new DnaSequence("AAAGAATTCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2 with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(2), "Single cut → 2 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(12),
                "Fragment sum must equal original sequence length (Addgene)");
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(4), "First fragment: AAAG");
            Assert.That(fragments[1].StartPosition, Is.EqualTo(4));
            Assert.That(fragments[1].Length, Is.EqualTo(8), "Second fragment: AATTCAAA");
        });
    }

    /// <summary>
    /// Evidence: Implementation contract - No cuts returns whole sequence.
    /// Invariant #7: Zero cut sites returns single fragment equal to original.
    /// </summary>
    [Test]
    public void Digest_NoCuts_ReturnsWholeSequenceAsSingleFragment()
    {
        // Arrange: Sequence without any EcoRI sites
        var sequence = new DnaSequence("AAAAAAAAAAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariant #7
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1), "No cuts should return single fragment");
            Assert.That(fragments[0].Length, Is.EqualTo(sequence.Length), "Fragment equals original");
            Assert.That(fragments[0].Sequence, Is.EqualTo(sequence.Sequence), "Content matches original");
            Assert.That(fragments[0].LeftEnzyme, Is.Null, "First fragment has no left enzyme");
            Assert.That(fragments[0].RightEnzyme, Is.Null, "Last fragment has no right enzyme");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Multiple cuts produce expected fragment count.
    /// Invariant #1: k cut positions → k+1 fragments.
    /// Invariant #2: Fragment sum = original sequence length.
    /// </summary>
    [Test]
    public void Digest_MultipleCuts_ReturnsCorrectFragmentCount()
    {
        // Arrange: Two EcoRI sites at positions 0 and 9, cuts at 1 and 10
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2 with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(18));
            Assert.That(fragments[0].Length, Is.EqualTo(1), "Fragment 1: G");
            Assert.That(fragments[1].Length, Is.EqualTo(9), "Fragment 2: AATTCAAAG");
            Assert.That(fragments[2].Length, Is.EqualTo(8), "Fragment 3: AATTCAAA");
        });
    }

    /// <summary>
    /// Evidence: Comprehensive fragment property verification with exact values.
    /// Invariants #3 (sequential numbering), #4 (monotonic positions),
    /// #5 (boundary enzymes), #6 (positive length).
    /// </summary>
    [Test]
    public void Digest_FragmentsHaveCorrectProperties()
    {
        // Arrange: 3 EcoRI sites at positions 0, 9, 18 → cuts at 1, 10, 19 → 4 fragments
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAAGAATTCAAA"); // 27 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: All fragment properties with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(4), "3 cuts → 4 fragments");

            // Fragment 1: [0, 1) = "G"
            Assert.That(fragments[0].FragmentNumber, Is.EqualTo(1));
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(1));
            Assert.That(fragments[0].LeftEnzyme, Is.Null, "First fragment: LeftEnzyme null (Invariant #5)");
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 2: [1, 10) = "AATTCAAAG"
            Assert.That(fragments[1].FragmentNumber, Is.EqualTo(2));
            Assert.That(fragments[1].StartPosition, Is.EqualTo(1));
            Assert.That(fragments[1].Length, Is.EqualTo(9));
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 3: [10, 19) = "AATTCAAAG"
            Assert.That(fragments[2].FragmentNumber, Is.EqualTo(3));
            Assert.That(fragments[2].StartPosition, Is.EqualTo(10));
            Assert.That(fragments[2].Length, Is.EqualTo(9));
            Assert.That(fragments[2].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[2].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 4: [19, 27) = "AATTCAAA"
            Assert.That(fragments[3].FragmentNumber, Is.EqualTo(4));
            Assert.That(fragments[3].StartPosition, Is.EqualTo(19));
            Assert.That(fragments[3].Length, Is.EqualTo(8));
            Assert.That(fragments[3].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[3].RightEnzyme, Is.Null, "Last fragment: RightEnzyme null (Invariant #5)");

            // Invariant #4: Monotonic start positions (absorbed from dedicated test)
            for (int i = 1; i < fragments.Count; i++)
                Assert.That(fragments[i].StartPosition, Is.GreaterThan(fragments[i - 1].StartPosition),
                    $"Fragment {i + 1} start > Fragment {i} start");

            // Invariant #6: All positive length
            Assert.That(fragments.All(f => f.Length > 0), "All fragments have positive length");
        });
    }

    /// <summary>
    /// Evidence: Verify fragment sequence content matches exact known substrings.
    /// </summary>
    [Test]
    public void Digest_FragmentSequenceContent_MatchesExpectedSubstring()
    {
        // Arrange: EcoRI (G↓AATTC) at position 3, cuts at position 4
        var sequence = new DnaSequence("AAAGAATTCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Exact known sequences
        Assert.Multiple(() =>
        {
            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAG"), "Fragment 1: positions [0, 4)");
            Assert.That(fragments[1].Sequence, Is.EqualTo("AATTCAAA"), "Fragment 2: positions [4, 12)");
        });
    }

    /// <summary>
    /// Evidence: Multiple enzymes produce combined cut sites with correct enzyme attribution.
    /// </summary>
    [Test]
    public void Digest_MultipleEnzymes_CutsWithBoth()
    {
        // Arrange: EcoRI at position 3 (cuts at 4), BamHI at position 12 (cuts at 13)
        var sequence = new DnaSequence("AAAGAATTCAAAGGATCCAAA"); // 21 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI", "BamHI").ToList();

        // Assert: Exact fragment properties with enzyme attribution
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(21));
            Assert.That(fragments[0].Length, Is.EqualTo(4), "Fragment 1: AAAG");
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].Length, Is.EqualTo(9), "Fragment 2: AATTCAAAG");
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.EqualTo("BamHI"));
            Assert.That(fragments[2].Length, Is.EqualTo(8), "Fragment 3: GATCCAAA");
            Assert.That(fragments[2].LeftEnzyme, Is.EqualTo("BamHI"));
        });
    }

    #endregion

    #region Digest Edge Cases

    /// <summary>
    /// Evidence: API contract - No enzymes provided must throw.
    /// </summary>
    [Test]
    public void Digest_NoEnzymes_ThrowsArgumentException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.Digest(sequence).ToList());
    }

    /// <summary>
    /// Evidence: API contract - Null sequence must throw.
    /// </summary>
    [Test]
    public void Digest_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.Digest(null!, "EcoRI").ToList());
    }

    /// <summary>
    /// Evidence: Edge case - Sequence shorter than recognition sequence.
    /// </summary>
    [Test]
    public void Digest_SequenceShorterThanRecognition_ReturnsWholeSequence()
    {
        // Arrange: 4 bp sequence, EcoRI needs 6 bp
        var sequence = new DnaSequence("ACGT");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1), "No possible cut sites");
            Assert.That(fragments[0].Length, Is.EqualTo(sequence.Length));
        });
    }

    /// <summary>
    /// Evidence: Edge case — adjacent cut positions produce a 1-bp fragment.
    /// TaqI (T↓CGA) at position 3 cuts at 4; MboI (↓GATC) at position 5 cuts at 5.
    /// All invariants must hold: positive length, sum = original.
    /// </summary>
    [Test]
    public void Digest_AdjacentCutSites_ProducesSmallFragment()
    {
        // Arrange: TaqI cuts at 4, MboI cuts at 5 → adjacent cuts, 1-bp middle fragment
        var sequence = new DnaSequence("AAATCGATCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "TaqI", "MboI").ToList();

        // Assert: Adjacent cuts handled correctly with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two adjacent cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(12));
            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAT"), "Fragment 1: [0, 4)");
            Assert.That(fragments[1].Sequence, Is.EqualTo("C"), "Fragment 2: [4, 5) — 1 bp");
            Assert.That(fragments[2].Sequence, Is.EqualTo("GATCAAA"), "Fragment 3: [5, 12)");
            Assert.That(fragments.All(f => f.Length > 0), "All fragments have positive length");
        });
    }

    #endregion

    #region Circular Digest Tests (RESTR-DIGEST-001 / C5)

    // Evidence for circular-molecule digest (k sites → k fragments; origin-spanning join):
    // - Univ. of Illinois MolBio, "Restriction Digestion/Gel Electrophoresis Assignment 1"
    //   (https://www.life.illinois.edu/molbio/geldigest/assign1.html):
    //   "If you cut a circle once, you get one linear fragment ... cut it a second time to
    //   get 2 linear fragments."
    // - Addgene, "Plasmids 101: How to Verify Your Plasmid Using a Restriction Digest"
    //   (https://blog.addgene.org/plasmids-101-how-to-verify-your-plasmid):
    //   a single cutter linearizes the plasmid → 1 fragment; fragment count = cut-site count.
    // - Quora summary of the standard rule (linear k+1 vs circular k):
    //   https://www.quora.com/How-shall-I-calculate-DNA-fragments-produced-when-a-circular-linear-plasmid-is-digested-with-a-restriction-enzyme-having-N-sites
    // Rule: circular molecule with k cut sites → k fragments. The origin-spanning fragment is
    // the join of last-cut→end and start→first-cut, length = (SequenceLength − lastCut) + firstCut.

    /// <summary>
    /// Worked circular example (hand-traced). EcoRI (G↓AATTC) on the 18 bp circle
    /// "GAATTCAAAGAATTCAAA" has forward-strand sites at 0 and 9, cutting at positions 1 and 10.
    /// LINEAR (control): cuts {1,10} → 3 fragments [0,1)="G", [1,10)="AATTCAAAG", [10,18)="AATTCAAA".
    /// CIRCULAR: 2 cut sites → 2 fragments. Fragment 1 spans cut1→cut2 = [1,10)="AATTCAAAG" (len 9).
    /// Fragment 2 is the origin-spanning join of [10,18)+[0,1) = "AATTCAAA"+"G" = "AATTCAAAG"
    /// (len = (18−10)+1 = 9). This is a DISCRIMINATING test: the linear implementation returns 3
    /// fragments here, the circular path returns exactly 2.
    /// </summary>
    [Test]
    public void Digest_Circular_TwoSites_ProducesTwoFragments_LinearProducesThree()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "EcoRI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            // Discriminating: linear k+1 = 3 vs circular k = 2.
            Assert.That(linear, Has.Count.EqualTo(3), "Linear: 2 cuts → 3 fragments (k+1)");
            Assert.That(circular, Has.Count.EqualTo(2), "Circular: 2 cuts → 2 fragments (k)");

            // Default overload must match Linear (non-breaking).
            var defaultDigest = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();
            Assert.That(defaultDigest.Select(f => f.Sequence),
                Is.EqualTo(linear.Select(f => f.Sequence)),
                "Default Digest overload must equal Linear topology");

            // Exact circular fragments (sequence + length + start).
            Assert.That(circular[0].StartPosition, Is.EqualTo(1));
            Assert.That(circular[0].Length, Is.EqualTo(9));
            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCAAAG"), "Cut1→Cut2 fragment");
            Assert.That(circular[0].LeftEnzyme, Is.EqualTo("EcoRI"), "Circular: both flanks are cuts");
            Assert.That(circular[0].RightEnzyme, Is.EqualTo("EcoRI"));

            Assert.That(circular[1].StartPosition, Is.EqualTo(10));
            Assert.That(circular[1].Length, Is.EqualTo(9), "Origin-spanning len = (18−10)+1 = 9");
            Assert.That(circular[1].Sequence, Is.EqualTo("AATTCAAAG"), "Origin-spanning join [10,18)+[0,1)");
            Assert.That(circular[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(circular[1].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment sum invariant still holds on a circle.
            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(18),
                "Circular fragment lengths sum to molecule length");
        });
    }

    /// <summary>
    /// Circular, 0 cut sites → a single full-length uncut circular fragment.
    /// Distinct from linear-0 only in topology semantics; both return the whole molecule.
    /// </summary>
    [Test]
    public void Digest_Circular_NoSites_ReturnsSingleFullLengthFragment()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAA"); // 12 bp, no EcoRI site

        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular, Has.Count.EqualTo(1), "0 sites → 1 uncut circular fragment");
            Assert.That(circular[0].Length, Is.EqualTo(12));
            Assert.That(circular[0].Sequence, Is.EqualTo(sequence.Sequence));
            Assert.That(circular[0].LeftEnzyme, Is.Null);
            Assert.That(circular[0].RightEnzyme, Is.Null);
        });
    }

    /// <summary>
    /// Circular, 1 cut site → the plasmid is linearized into ONE full-length fragment.
    /// EcoRI on "AAAGAATTCAAA" (12 bp) cuts at position 4. Circular fragment is the
    /// origin-spanning join [4,12)+[0,4) = "AATTCAAA"+"AAAG" = "AATTCAAAAAAG" (len 12).
    /// DISCRIMINATING: linear here yields 2 fragments, circular yields exactly 1.
    /// </summary>
    [Test]
    public void Digest_Circular_SingleSite_LinearizesToOneFragment_LinearProducesTwo()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA"); // 12 bp, EcoRI cut at 4

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "EcoRI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(linear, Has.Count.EqualTo(2), "Linear: 1 cut → 2 fragments");
            Assert.That(circular, Has.Count.EqualTo(1), "Circular: 1 cut → 1 linearized fragment");

            Assert.That(circular[0].StartPosition, Is.EqualTo(4));
            Assert.That(circular[0].Length, Is.EqualTo(12), "Full-length linearized: (12−4)+4 = 12");
            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCAAAAAAG"),
                "Origin-spanning join [4,12)+[0,4)");
            Assert.That(circular[0].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(circular[0].RightEnzyme, Is.EqualTo("EcoRI"));
        });
    }

    /// <summary>
    /// Three cut sites on a circle → exactly 3 fragments (vs 4 linear), with the origin-spanning
    /// fragment joining last-cut→end and start→first-cut. EcoRI on the 27 bp circle
    /// "GAATTCAAAGAATTCAAAGAATTCAAA" has sites at 0, 9, 18 → cuts at 1, 10, 19.
    /// CIRCULAR fragments: [1,10)="AATTCAAAG" (9), [10,19)="AATTCAAAG" (9),
    /// origin-spanning [19,27)+[0,1)="AATTCAAA"+"G"="AATTCAAAG" (len (27−19)+1 = 9).
    /// </summary>
    [Test]
    public void Digest_Circular_ThreeSites_ProducesThreeFragments_WithOriginSpan()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAAGAATTCAAA"); // 27 bp

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "EcoRI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(linear, Has.Count.EqualTo(4), "Linear: 3 cuts → 4 fragments");
            Assert.That(circular, Has.Count.EqualTo(3), "Circular: 3 cuts → 3 fragments");

            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCAAAG"));
            Assert.That(circular[1].Sequence, Is.EqualTo("AATTCAAAG"));

            // Origin-spanning fragment.
            Assert.That(circular[2].StartPosition, Is.EqualTo(19));
            Assert.That(circular[2].Length, Is.EqualTo(9), "Origin-spanning len = (27−19)+1 = 9");
            Assert.That(circular[2].Sequence, Is.EqualTo("AATTCAAAG"), "Join [19,27)+[0,1)");

            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(27),
                "Circular fragment lengths sum to molecule length");
        });
    }

    /// <summary>
    /// API contract: the topology overload still rejects an empty enzyme list and null sequence.
    /// </summary>
    [Test]
    public void Digest_Circular_NoEnzymesOrNullSequence_Throws()
    {
        var sequence = new DnaSequence("GAATTCAAA");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular).ToList());
            Assert.Throws<ArgumentNullException>(
                () => RestrictionAnalyzer.Digest(null!, MoleculeTopology.Circular, "EcoRI").ToList());
        });
    }

    /// <summary>
    /// Independent validator test (RESTR-DIGEST-001 re-validation). The existing circular tests
    /// all use tandem-repeat sequences in which every fragment — including the origin-spanning
    /// one — happens to read "AATTCAAAG", so a mis-implemented wrap-around join (wrong slice
    /// offset, or head+tail swapped) could coincidentally still pass the sequence assertion.
    /// This case uses a NON-repetitive sequence whose origin-spanning fragment ("AATTCAACG")
    /// is unique: it does not appear as any contiguous substring of the molecule, so it can only
    /// arise from the correct tail([lastCut,len)) + head([0,firstCut)) join.
    ///
    /// Hand-derived independently (EcoRI cuts G↓AATTC, confirmed Wikipedia "EcoRI"):
    /// "CGAATTCTTTGAATTCAA" (18 bp), GAATTC starts at 1 and 10 → cuts at {2, 11}.
    ///   Fragment 1: [2,11)  = "AATTCTTTG" (len 9)
    ///   Fragment 2: origin-spanning [11,18)+[0,2) = "AATTCAA"+"CG" = "AATTCAACG"
    ///               (len = (18−11)+2 = 9)
    /// DISCRIMINATING: linear yields k+1 = 3 fragments; circular yields exactly k = 2.
    /// Fragment-count rule sourced: Univ. of Illinois MolBio "Restriction Digestion/Gel
    /// Electrophoresis" ("If you cut a circle once, you get one linear fragment ... cut it a
    /// second time to get 2 linear fragments"); Addgene Plasmids 101 (single cutter linearizes).
    /// </summary>
    [Test]
    public void Digest_Circular_TwoSites_OriginSpanningSequenceIsUnique_NotTandemRepeat()
    {
        var sequence = new DnaSequence("CGAATTCTTTGAATTCAA"); // 18 bp, NON-repetitive

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "EcoRI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(linear, Has.Count.EqualTo(3), "Linear: 2 cuts → 3 fragments (k+1)");
            Assert.That(circular, Has.Count.EqualTo(2), "Circular: 2 cuts → 2 fragments (k)");

            // Cut1→Cut2 fragment.
            Assert.That(circular[0].StartPosition, Is.EqualTo(2));
            Assert.That(circular[0].Length, Is.EqualTo(9));
            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCTTTG"), "Cut1→Cut2 fragment");

            // Origin-spanning fragment: its sequence is UNIQUE and only forms from tail+head.
            Assert.That(circular[1].StartPosition, Is.EqualTo(11));
            Assert.That(circular[1].Length, Is.EqualTo(9), "Origin-spanning len = (18−11)+2 = 9");
            Assert.That(circular[1].Sequence, Is.EqualTo("AATTCAACG"),
                "Origin-spanning join is tail [11,18) + head [0,2), not any internal slice");

            // The origin-spanning sequence must not be confusable with the other fragment.
            Assert.That(circular[0].Sequence, Is.Not.EqualTo(circular[1].Sequence),
                "Distinct-content fragments guard against tandem-repeat self-confirmation");

            // The wrap sequence is genuinely absent as a contiguous substring of the molecule.
            Assert.That(sequence.Sequence.Contains("AATTCAACG"), Is.False,
                "Origin-spanning content can only come from a correct wrap-around join");

            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(18),
                "Circular fragment lengths sum to molecule length");
        });
    }

    #endregion

    #region Circular / Overhang Re-validation (RESTR-DIGEST-001, 2026-10-14)

    // Reference oracle for this region: Biopython 1.88 Bio.Restriction (Restriction.py).
    // Every expected value below was produced in-session by running Biopython 1.88 and is quoted
    // with the exact call. Biopython reports 1-based top-strand cut positions, so a Biopython cut
    // at p corresponds to Seqeron CutPosition p-1.
    //
    // NOTE on the oracle: Biopython's own catalyse(..., linear=False) is unsafe as a fragment
    // oracle because AbstractCut._drop (Restriction.py L542-568) folds an out-of-range circular
    // cut back by +/-length in place WITHOUT re-sorting, and catalyse (L1146-1176) then indexes
    // r[0]/r[-1] assuming ascending order. See
    // Digest_Circular_BiopythonUnsortedCutCase_FragmentsStillPartitionMolecule below. Where that
    // happens the oracle used is sorted(search(..., linear=False)).

    /// <summary>
    /// DEFECT REGRESSION (claim 1 of the RESTR-FIND-001 hand-off, reproduced 2026-10-14):
    /// <c>DigestCircular</c> called the LINEAR site finder, so a circular digest could not see an
    /// origin-spanning cut. Before the fix,
    /// <c>Digest(new DnaSequence("TTCAAAAAAAGAA"), MoleculeTopology.Circular, "EcoRI")</c> returned
    /// 1 fragment with StartPosition 0, Length 13 and both flanking enzymes null — an UNCUT circle
    /// — although the plasmid carries an origin-spanning GAATTC at position 10.
    ///
    /// Reference (Biopython 1.88):
    ///   Bio.Restriction.EcoRI.search(Seq("TTCAAAAAAAGAA"))                == []   (linear)
    ///   Bio.Restriction.EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False)  == [12] (1-based)
    ///   Bio.Restriction.EcoRI.catalyse(Seq("TTCAAAAAAAGAA"), linear=False)
    ///       == (Seq("AATTCAAAAAAAG"),)
    /// So the circular digest must linearize the plasmid into ONE 13 bp fragment starting at the
    /// cut (0-based 11) and wrapping past the origin; the LINEAR digest of the same sequence is
    /// correctly uncut.
    /// </summary>
    [Test]
    public void Digest_Circular_OriginSpanningSite_IsCut_MatchesBiopython()
    {
        var sequence = new DnaSequence("TTCAAAAAAAGAA"); // 13 bp; GAATTC spans the origin at 10

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "EcoRI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            // Linear control: Biopython EcoRI.search(...) == [] -> no cut -> whole molecule.
            Assert.That(linear, Has.Count.EqualTo(1), "Linear: the site does not exist on a linear molecule");
            Assert.That(linear[0].Sequence, Is.EqualTo("TTCAAAAAAAGAA"));
            Assert.That(linear[0].LeftEnzyme, Is.Null);
            Assert.That(linear[0].RightEnzyme, Is.Null);

            // Circular: Biopython catalyse(..., linear=False) == ("AATTCAAAAAAAG",).
            Assert.That(circular, Has.Count.EqualTo(1), "Circular: origin-spanning site linearizes the plasmid");
            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCAAAAAAAG"),
                "Biopython EcoRI.catalyse(Seq(\"TTCAAAAAAAGAA\"), linear=False) == (Seq(\"AATTCAAAAAAAG\"),)");
            Assert.That(circular[0].StartPosition, Is.EqualTo(11),
                "Biopython search(..., linear=False) == [12] (1-based) -> 0-based cut 11");
            Assert.That(circular[0].Length, Is.EqualTo(13));
            Assert.That(circular[0].LeftEnzyme, Is.EqualTo("EcoRI"), "Regression: was null (uncut circle)");
            Assert.That(circular[0].RightEnzyme, Is.EqualTo("EcoRI"), "Regression: was null (uncut circle)");
        });
    }

    /// <summary>
    /// Second consequence of claim 1: on a closed circle the linear-only cut-feasibility filter
    /// must not apply, so an edge cutter cuts at every site. Sau3AI is ^GATC (CutPositionForward 0).
    ///
    /// Reference (Biopython 1.88):
    ///   Sau3AI.search(Seq("GATCGATCGATC"))                == [5]
    ///   Sau3AI.catalyse(Seq("GATCGATCGATC"))              == ("GATC", "GATCGATC")
    ///   Sau3AI.search(Seq("GATCGATCGATC"), linear=False)   == [1, 5, 9]
    ///   Sau3AI.catalyse(Seq("GATCGATCGATC"), linear=False) == ("GATC", "GATC", "GATC")
    /// Before the fix the circular digest returned 1 uncut fragment of 12 bp.
    /// </summary>
    [Test]
    public void Digest_Circular_EdgeCutter_CutsEverySite_MatchesBiopython()
    {
        var sequence = new DnaSequence("GATCGATCGATC"); // 12 bp

        var linear = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, "Sau3AI").ToList();
        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "Sau3AI").ToList();

        Assert.Multiple(() =>
        {
            // Linear: Biopython drops the two edge cuts -> one internal cut at 0-based 4.
            Assert.That(linear.Select(f => f.Sequence), Is.EqualTo(new[] { "GATC", "GATCGATC" }),
                "Biopython Sau3AI.catalyse(Seq(\"GATCGATCGATC\")) == (\"GATC\", \"GATCGATC\")");

            // Circular: 3 cuts {0, 4, 8} -> 3 fragments.
            Assert.That(circular.Select(f => f.Sequence), Is.EqualTo(new[] { "GATC", "GATC", "GATC" }),
                "Biopython Sau3AI.catalyse(..., linear=False) == (\"GATC\", \"GATC\", \"GATC\")");
            Assert.That(circular.Select(f => f.StartPosition), Is.EqualTo(new[] { 0, 4, 8 }),
                "Biopython search(..., linear=False) == [1, 5, 9] (1-based)");
            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(12));
        });
    }

    /// <summary>
    /// Circular DOUBLE digest with one origin-spanning site. Reference (Biopython 1.88):
    ///   RestrictionBatch([EcoRI, BamHI]).search(Seq("TTCGGATCCAAAAAAAGAA"), linear=False)
    ///       == {EcoRI: [18], BamHI: [5]}
    /// i.e. 0-based top-strand cuts {4, 17}. Circular -> 2 fragments:
    ///   [4,17)            = "GATCCAAAAAAAG"  (13, BamHI -> EcoRI)
    ///   [17,19) + [0,4)   = "AA" + "TTCG" = "AATTCG" (6, EcoRI -> BamHI, origin-spanning, last)
    /// Before the fix the circular digest saw only the BamHI cut and returned the single
    /// 19 bp fragment "GATCCAAAAAAAGAATTCG".
    /// </summary>
    [Test]
    public void Digest_Circular_DoubleDigest_OriginSpanningSite_MatchesBiopython()
    {
        var sequence = new DnaSequence("TTCGGATCCAAAAAAAGAA"); // 19 bp

        var circular = RestrictionAnalyzer.Digest(
            sequence, MoleculeTopology.Circular, "EcoRI", "BamHI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular, Has.Count.EqualTo(2), "2 cuts on a circle -> 2 fragments");
            Assert.That(circular[0].StartPosition, Is.EqualTo(4), "BamHI cut: Biopython 5 (1-based)");
            Assert.That(circular[0].Sequence, Is.EqualTo("GATCCAAAAAAAG"));
            Assert.That(circular[0].LeftEnzyme, Is.EqualTo("BamHI"));
            Assert.That(circular[0].RightEnzyme, Is.EqualTo("EcoRI"));

            Assert.That(circular[1].StartPosition, Is.EqualTo(17), "EcoRI cut: Biopython 18 (1-based)");
            Assert.That(circular[1].Sequence, Is.EqualTo("AATTCG"), "Origin-spanning join [17,19)+[0,4)");
            Assert.That(circular[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(circular[1].RightEnzyme, Is.EqualTo("BamHI"));

            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(19));
        });
    }

    /// <summary>
    /// Circular multi-enzyme edge-cutter double digest. Reference (Biopython 1.88):
    ///   RestrictionBatch([MboI, DpnI]).search(Seq("GATCGATCGATC"), linear=False)
    ///       == {MboI: [1, 5, 9], DpnI: [3, 7, 11]}
    /// i.e. 0-based cuts {0, 2, 4, 6, 8, 10} -> 6 fragments of 2 bp each:
    /// GA, TC, GA, TC, GA and the origin-spanning TC. Before the fix the circular digest returned
    /// 4 fragments ("TC", "GA", "TCGA", "TCGA"), because the linear filter dropped MboI's edge cuts.
    /// </summary>
    [Test]
    public void Digest_Circular_MultiEnzymeEdgeCutters_MatchesBiopython()
    {
        var sequence = new DnaSequence("GATCGATCGATC"); // 12 bp

        var circular = RestrictionAnalyzer.Digest(
            sequence, MoleculeTopology.Circular, "MboI", "DpnI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular.Select(f => f.Sequence),
                Is.EqualTo(new[] { "GA", "TC", "GA", "TC", "GA", "TC" }));
            Assert.That(circular.Select(f => f.StartPosition),
                Is.EqualTo(new[] { 0, 2, 4, 6, 8, 10 }));
            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(12));
            Assert.That(circular.Select(f => f.FragmentNumber), Is.EqualTo(new[] { 1, 2, 3, 4, 5, 6 }));
        });
    }

    /// <summary>
    /// The reference implementation's own circular bug, locked as a Seqeron guarantee.
    /// Biopython 1.88 <c>AbstractCut._drop</c> (Restriction.py L542-568) folds an out-of-range
    /// circular cut back by ±length in place WITHOUT re-sorting the list, and
    /// <c>catalyse</c> (L1146-1176) then indexes <c>r[0]</c>/<c>r[-1]</c> as if it were ascending:
    ///   AluI.search(Seq("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"), linear=False)   == [11, 2]   (unsorted)
    ///   AluI.catalyse(Seq("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"), linear=False)
    ///       == ("CTGTTATAGCTAAGAGCGAACCGGAGTATAGCTGTTATAG", "")   -> 40 bases for a 31 bp molecule
    /// Measured in-session: 55 of 4000 wrap-planted circular cases (1.4%) trigger it.
    /// Seqeron sorts its cut positions, so the correct answer is cuts {1, 10} (sorted Biopython
    /// [2, 11] minus 1) -> 2 fragments that partition the 31 bp molecule exactly.
    /// </summary>
    [Test]
    public void Digest_Circular_BiopythonUnsortedCutCase_FragmentsStillPartitionMolecule()
    {
        var sequence = new DnaSequence("GCTGTTATAGCTAAGAGCGAACCGGAGTATA"); // 31 bp; AluI = AG^CT

        var circular = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Circular, "AluI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular, Has.Count.EqualTo(2), "2 cuts on a circle -> 2 fragments");
            Assert.That(circular[0].StartPosition, Is.EqualTo(1));
            Assert.That(circular[0].Sequence, Is.EqualTo("CTGTTATAG"));
            Assert.That(circular[1].StartPosition, Is.EqualTo(10));
            Assert.That(circular[1].Sequence, Is.EqualTo("CTAAGAGCGAACCGGAGTATAG"),
                "Origin-spanning join [10,31)+[0,1)");
            Assert.That(circular.Sum(f => f.Length), Is.EqualTo(31),
                "Unlike Biopython's catalyse here, the fragments must partition the molecule exactly");
            Assert.That(circular, Has.None.Matches<DigestFragment>(f => f.Length == 0),
                "No empty fragment (Biopython emits one here)");
        });
    }

    /// <summary>
    /// Claim 2 of the hand-off: <c>FragmentNumber</c> used to be the boundary loop index, so a
    /// skipped zero-length piece would have broken the contiguous numbering. Verified this session
    /// that a zero-length piece is UNREACHABLE through the public API — the linear cut-feasibility
    /// filter (<c>CutsAreInternal</c>) confines every reported forward-strand cut to
    /// <c>[1, Length-1]</c>, and the positions are distinct — so the "k cuts -> k+1 fragments"
    /// invariant holds. This test enforces both halves on the densest cut patterns available:
    /// every fragment is non-empty, the numbering is 1..n contiguous, and the count is k+1
    /// (linear) / k (circular) for the k distinct forward-strand cut positions actually reported.
    /// </summary>
    [Test]
    public void Digest_FragmentNumbering_IsContiguous_AndNoZeroLengthFragments()
    {
        // Dense / adjacent-cut cases: 4-cutters on a 4-mer-repeat molecule, plus edge cutters.
        (string Seq, string[] Enzymes)[] cases =
        [
            ("GATCGATCGATC", ["MboI", "DpnI"]),
            ("GATCGATCGATC", ["Sau3AI"]),
            ("AGCTAGCTAGCT", ["AluI"]),
            ("CCGGCCGGCCGG", ["MspI", "HpaII"]),
            ("GATCGATCGATC", ["MboI", "DpnI", "TaqI", "AluI"]),
            ("GAATTCGAATTC", ["EcoRI"]),
            ("CTGCAGCTGCAG", ["PstI"]),
        ];

        Assert.Multiple(() =>
        {
            foreach (var (seqText, enzymes) in cases)
            {
                var sequence = new DnaSequence(seqText);

                foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
                {
                    var fragments = RestrictionAnalyzer.Digest(sequence, topology, enzymes).ToList();
                    string label = $"{seqText} / {string.Join("+", enzymes)} / {topology}";

                    int cuts = enzymes
                        .SelectMany(e => RestrictionAnalyzer.FindSites(sequence, e, topology))
                        .Where(s => s.IsForwardStrand)
                        .Select(s => s.CutPosition)
                        .Distinct()
                        .Count();

                    int expected = topology == MoleculeTopology.Linear
                        ? cuts + 1
                        : Math.Max(cuts, 1);

                    Assert.That(fragments, Has.Count.EqualTo(expected), $"fragment count: {label}");
                    Assert.That(fragments.Select(f => f.FragmentNumber),
                        Is.EqualTo(Enumerable.Range(1, fragments.Count).ToList()),
                        $"FragmentNumber must be 1..n contiguous: {label}");
                    Assert.That(fragments.All(f => f.Length > 0), Is.True,
                        $"no zero-length fragment: {label}");
                    Assert.That(fragments.All(f => f.Length == f.Sequence.Length), Is.True,
                        $"Length must equal Sequence.Length: {label}");
                    Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(sequence.Length),
                        $"fragment lengths must sum to the molecule length: {label}");
                    Assert.That(fragments.Select(f => f.StartPosition),
                        Is.Ordered.Ascending, $"StartPosition ascending: {label}");
                }
            }
        });
    }

    /// <summary>
    /// Claim 3 of the hand-off: the digest dropped the reverse-strand site records, leaving every
    /// fragment implicitly blunt with nothing recorded about its ends. Dropping the reverse-strand
    /// record is CORRECT — it carries the bottom-strand cut of the same palindromic site, and
    /// Biopython's <c>catalyse</c> likewise partitions the top strand at <c>search()</c>'s
    /// top-strand cuts — but the bottom-strand cut is what makes the end sticky, so it is now
    /// reported per fragment end.
    ///
    /// Reference (Biopython 1.88): EcoRI.ovhgseq == "AATT", is_5overhang() == True;
    /// PstI.ovhgseq == "TGCA", is_3overhang() == True; EcoRV.ovhgseq == "", is_blunt() == True;
    /// AluI.ovhgseq == "", is_blunt() == True; BamHI.ovhgseq == "GATC"; Sau3AI.ovhgseq == "GATC".
    /// Verified in-session for all 39 built-in enzymes (39/39 agreement on both the overhang
    /// sequence and the 5'/3'/blunt classification).
    /// </summary>
    [Test]
    public void Digest_FragmentEnds_ReportSourcedOverhangTypeAndSequence()
    {
        // EcoRI G^AATT_C on "AAAGAATTCAAA": Biopython catalyse == ("AAAG", "AATTCAAA").
        var ecoRI = RestrictionAnalyzer.Digest(new DnaSequence("AAAGAATTCAAA"), "EcoRI").ToList();

        // PstI C_TGCA^G on "CTGCAGAAACTGCAG": Biopython catalyse == ("CTGCA", "GAAACTGCA", "G").
        var pstI = RestrictionAnalyzer.Digest(new DnaSequence("CTGCAGAAACTGCAG"), "PstI").ToList();

        // EcoRV GAT^ATC (blunt) on "AAAGATATCAAA": cut at 0-based 6.
        var ecoRV = RestrictionAnalyzer.Digest(new DnaSequence("AAAGATATCAAA"), "EcoRV").ToList();

        Assert.Multiple(() =>
        {
            // 5' overhang (EcoRI, AATT).
            Assert.That(ecoRI.Select(f => f.Sequence), Is.EqualTo(new[] { "AAAG", "AATTCAAA" }));
            Assert.That(ecoRI[0].RightOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(ecoRI[0].RightOverhang, Is.EqualTo("AATT"), "Biopython EcoRI.ovhgseq == \"AATT\"");
            Assert.That(ecoRI[1].LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(ecoRI[1].LeftOverhang, Is.EqualTo("AATT"));

            // 3' overhang (PstI, TGCA).
            Assert.That(pstI.Select(f => f.Sequence), Is.EqualTo(new[] { "CTGCA", "GAAACTGCA", "G" }),
                "Biopython PstI.catalyse(Seq(\"CTGCAGAAACTGCAG\")) == (\"CTGCA\", \"GAAACTGCA\", \"G\")");
            Assert.That(pstI[0].RightOverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(pstI[0].RightOverhang, Is.EqualTo("TGCA"), "Biopython PstI.ovhgseq == \"TGCA\"");
            Assert.That(pstI[1].LeftOverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(pstI[1].LeftOverhang, Is.EqualTo("TGCA"));
            Assert.That(pstI[1].RightOverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(pstI[2].LeftOverhang, Is.EqualTo("TGCA"));

            // Blunt (EcoRV): type Blunt, empty overhang - NOT null, which means "no cut here".
            Assert.That(ecoRV, Has.Count.EqualTo(2));
            Assert.That(ecoRV[0].RightOverhangType, Is.EqualTo(OverhangType.Blunt));
            Assert.That(ecoRV[0].RightOverhang, Is.Empty, "Biopython EcoRV.ovhgseq == \"\"");
            Assert.That(ecoRV[1].LeftOverhangType, Is.EqualTo(OverhangType.Blunt));
            Assert.That(ecoRV[1].LeftOverhang, Is.Empty);

            // An original terminus of a LINEAR molecule was not created by a cut -> null.
            Assert.That(ecoRI[0].LeftOverhangType, Is.Null, "molecule 5' terminus is not a cut end");
            Assert.That(ecoRI[0].LeftOverhang, Is.Null);
            Assert.That(ecoRI[^1].RightOverhangType, Is.Null, "molecule 3' terminus is not a cut end");
            Assert.That(ecoRI[^1].RightOverhang, Is.Null);
        });
    }

    /// <summary>
    /// On a circle every fragment end is a cut, so both ends carry an enzyme AND an overhang; an
    /// uncut circle has no ends at all and therefore no overhangs. Also locks the ligation
    /// invariant: at one cut, the upstream fragment's right overhang is the reverse complement of
    /// the downstream fragment's left overhang (identical for a palindromic enzyme — which is why
    /// the fragments of a single-enzyme digest re-ligate).
    /// </summary>
    [Test]
    public void Digest_Circular_EveryEndIsACut_UncutCircleHasNoEnds()
    {
        // Two EcoRI cuts on an 18 bp circle (existing worked example).
        var cut = RestrictionAnalyzer.Digest(
            new DnaSequence("GAATTCAAAGAATTCAAA"), MoleculeTopology.Circular, "EcoRI").ToList();

        // No EcoRI site at all.
        var uncut = RestrictionAnalyzer.Digest(
            new DnaSequence("AAAAAAAAAAAA"), MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(cut, Has.Count.EqualTo(2));
            foreach (var f in cut)
            {
                Assert.That(f.LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
                Assert.That(f.LeftOverhang, Is.EqualTo("AATT"));
                Assert.That(f.RightOverhangType, Is.EqualTo(OverhangType.FivePrime));
                Assert.That(f.RightOverhang, Is.EqualTo("AATT"));
            }

            // Ligation invariant across each cut (EcoRI's AATT is self-complementary).
            Assert.That(cut[0].RightOverhang,
                Is.EqualTo(DnaSequence.GetReverseComplementString(cut[1].LeftOverhang!)));

            Assert.That(uncut, Has.Count.EqualTo(1), "0 cuts -> one uncut circular fragment");
            Assert.That(uncut[0].LeftOverhangType, Is.Null, "an uncut circle has no ends");
            Assert.That(uncut[0].LeftOverhang, Is.Null);
            Assert.That(uncut[0].RightOverhangType, Is.Null);
            Assert.That(uncut[0].RightOverhang, Is.Null);
        });
    }

    /// <summary>
    /// The overhang window is taken from the MOLECULE, so it wraps correctly when a cut sits close
    /// enough to the origin that the enzyme's two strand cuts straddle it. On the 13 bp circle
    /// "TTCAAAAAAAGAA" the origin-spanning EcoRI site at 10 has its top-strand cut at 11 and its
    /// bottom-strand cut at 11+4 = 15 == 2 (mod 13), so the overhang window is
    /// seq[11], seq[12], seq[0], seq[1] = "A","A","T","T" = "AATT" — Biopython EcoRI.ovhgseq.
    /// </summary>
    [Test]
    public void Digest_Circular_OverhangWindow_WrapsPastTheOrigin()
    {
        var fragments = RestrictionAnalyzer.Digest(
            new DnaSequence("TTCAAAAAAAGAA"), MoleculeTopology.Circular, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1));
            Assert.That(fragments[0].LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(fragments[0].LeftOverhang, Is.EqualTo("AATT"),
                "Window wraps: seq[11],seq[12],seq[0],seq[1]");
            Assert.That(fragments[0].RightOverhang, Is.EqualTo("AATT"));
        });
    }

    /// <summary>
    /// Claim 4 of the hand-off: <c>GetDigestSummary</c> had no topology overload, so a circular
    /// digest summary was unreachable. The additive overload must agree with the corresponding
    /// <c>Digest</c> call, and the topology-free overload must stay exactly Linear.
    ///
    /// EcoRI on the 18 bp circle "GAATTCAAAGAATTCAAA": Biopython search(..., linear=False)
    /// == [2, 11] -> cuts {1, 10}. Linear: 3 fragments, sizes [9, 8, 1]. Circular: 2 fragments,
    /// sizes [9, 9].
    /// </summary>
    [Test]
    public void GetDigestSummary_TopologyOverload_MatchesDigest_AndDefaultsToLinear()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        var linear = RestrictionAnalyzer.GetDigestSummary(sequence, MoleculeTopology.Linear, "EcoRI");
        var circular = RestrictionAnalyzer.GetDigestSummary(sequence, MoleculeTopology.Circular, "EcoRI");
        var legacy = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        Assert.Multiple(() =>
        {
            Assert.That(linear.TotalFragments, Is.EqualTo(3), "Linear: 2 cuts -> 3 fragments");
            Assert.That(linear.FragmentSizes, Is.EqualTo(new[] { 9, 8, 1 }));
            Assert.That(linear.Topology, Is.EqualTo(MoleculeTopology.Linear));

            Assert.That(circular.TotalFragments, Is.EqualTo(2), "Circular: 2 cuts -> 2 fragments");
            Assert.That(circular.FragmentSizes, Is.EqualTo(new[] { 9, 9 }));
            Assert.That(circular.LargestFragment, Is.EqualTo(9));
            Assert.That(circular.SmallestFragment, Is.EqualTo(9));
            Assert.That(circular.AverageFragmentSize, Is.EqualTo(9.0).Within(1e-9));
            Assert.That(circular.EnzymesUsed, Is.EqualTo(new[] { "EcoRI" }));
            Assert.That(circular.Topology, Is.EqualTo(MoleculeTopology.Circular));

            // Addgene: the fragment sizes sum to the original length in BOTH topologies.
            Assert.That(linear.FragmentSizes.Sum(), Is.EqualTo(18));
            Assert.That(circular.FragmentSizes.Sum(), Is.EqualTo(18));

            // Backward compatibility: the topology-free overload is exactly Linear. (Compared
            // field by field: DigestSummary's list members are List<T>, which the compiler-
            // generated record equality compares by reference, not by content.)
            Assert.That(legacy.TotalFragments, Is.EqualTo(linear.TotalFragments));
            Assert.That(legacy.FragmentSizes, Is.EqualTo(linear.FragmentSizes),
                "GetDigestSummary(seq, enzymes) == Linear overload");
            Assert.That(legacy.LargestFragment, Is.EqualTo(linear.LargestFragment));
            Assert.That(legacy.SmallestFragment, Is.EqualTo(linear.SmallestFragment));
            Assert.That(legacy.AverageFragmentSize, Is.EqualTo(linear.AverageFragmentSize).Within(1e-9));
            Assert.That(legacy.EnzymesUsed, Is.EqualTo(linear.EnzymesUsed));
            Assert.That(legacy.Topology, Is.EqualTo(MoleculeTopology.Linear),
                "the topology-free overload must default to Linear");
        });
    }

    /// <summary>
    /// A circular summary over an ORIGIN-SPANNING site is the case that was unreachable before:
    /// EcoRI on the 13 bp circle "TTCAAAAAAAGAA" -> 1 fragment of 13 bp, whereas the linear
    /// summary of the same molecule is also 1 fragment of 13 bp but for the opposite reason
    /// (uncut). The discriminating molecule is the 19 bp circle "TTCGGATCCAAAAAAAGAA" digested
    /// with EcoRI + BamHI: linear sees only BamHI (cut 4) -> sizes [15, 4]; circular also sees the
    /// origin-spanning EcoRI (cut 17) -> sizes [13, 6].
    /// </summary>
    [Test]
    public void GetDigestSummary_Circular_SeesOriginSpanningSite()
    {
        var sequence = new DnaSequence("TTCGGATCCAAAAAAAGAA"); // 19 bp

        var linear = RestrictionAnalyzer.GetDigestSummary(
            sequence, MoleculeTopology.Linear, "EcoRI", "BamHI");
        var circular = RestrictionAnalyzer.GetDigestSummary(
            sequence, MoleculeTopology.Circular, "EcoRI", "BamHI");

        Assert.Multiple(() =>
        {
            Assert.That(linear.TotalFragments, Is.EqualTo(2));
            Assert.That(linear.FragmentSizes, Is.EqualTo(new[] { 15, 4 }));

            Assert.That(circular.TotalFragments, Is.EqualTo(2));
            Assert.That(circular.FragmentSizes, Is.EqualTo(new[] { 13, 6 }),
                "Biopython batch search(..., linear=False) == {EcoRI: [18], BamHI: [5]}");
            Assert.That(circular.FragmentSizes.Sum(), Is.EqualTo(19));
        });
    }

    /// <summary>
    /// Several enzymes cutting at the SAME coordinate: the boundary is one cut, not two, and the
    /// fragment keeps one representative enzyme name (the last seen) with that enzyme's overhang.
    /// Reference (Biopython 1.88):
    ///   RestrictionBatch([MspI, HpaII]).search(Seq("CCGGAAACCGG"))      == {MspI: [2,9], HpaII: [2,9]}
    ///   RestrictionBatch([MboI, Sau3AI, DpnI]).search(Seq("GATCAAAGATC"), linear=False)
    ///       == {MboI: [1,8], Sau3AI: [1,8], DpnI: [3,10]}
    /// i.e. 0-based cuts {1, 8} and {0, 2, 7, 9} — the duplicate coordinate must be merged.
    /// </summary>
    [Test]
    public void Digest_EnzymesSharingACutCoordinate_ProduceOneBoundary()
    {
        // MspI and HpaII are the same site (CCGG) with the same offsets, so every cut coincides.
        var shared = RestrictionAnalyzer.Digest(
            new DnaSequence("CCGGAAACCGG"), MoleculeTopology.Linear, "MspI", "HpaII").ToList();

        // MboI/Sau3AI (^GATC) and DpnI (GA^TC) share the site but cut 2 bases apart.
        var mixed = RestrictionAnalyzer.Digest(
            new DnaSequence("GATCAAAGATC"), MoleculeTopology.Circular, "MboI", "Sau3AI", "DpnI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(shared.Select(f => f.Sequence), Is.EqualTo(new[] { "C", "CGGAAAC", "CGG" }),
                "2 merged cuts -> 3 linear fragments, not 5");
            Assert.That(shared.Select(f => f.StartPosition), Is.EqualTo(new[] { 0, 1, 8 }));
            Assert.That(shared.Sum(f => f.Length), Is.EqualTo(11));
            Assert.That(shared[1].LeftOverhang, Is.EqualTo("CG"), "Biopython MspI/HpaII ovhgseq == \"CG\"");

            Assert.That(mixed.Select(f => f.Sequence), Is.EqualTo(new[] { "GA", "TCAAA", "GA", "TC" }),
                "4 distinct cuts on a circle -> 4 fragments");
            Assert.That(mixed.Select(f => f.StartPosition), Is.EqualTo(new[] { 0, 2, 7, 9 }));
            Assert.That(mixed.Sum(f => f.Length), Is.EqualTo(11));
            // ^GATC leaves a 5' GATC overhang; GA^TC is blunt.
            Assert.That(mixed[0].LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(mixed[0].LeftOverhang, Is.EqualTo("GATC"));
            Assert.That(mixed[0].RightOverhangType, Is.EqualTo(OverhangType.Blunt));
            Assert.That(mixed[0].RightOverhang, Is.Empty);
        });
    }

    #endregion

    #region GetDigestSummary Tests

    /// <summary>
    /// Evidence: Summary invariants verified with exact values.
    /// Invariants #8 (sorted descending), #9 (bounds), #10 (enzyme list).
    /// EcoRI sites at positions 0 and 9 → fragments of length 1, 9, 8.
    /// </summary>
    [Test]
    public void GetDigestSummary_ReturnsCorrectSummaryWithInvariants()
    {
        // Arrange: 2 EcoRI sites → fragments [1, 9, 8], sorted desc: [9, 8, 1]
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        // Act
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        // Assert: All invariants with exact values
        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalFragments, Is.EqualTo(3));

            // Invariant #8: Fragment sizes sorted descending — exact values
            Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 9, 8, 1 }),
                "Fragment sizes sorted descending");

            // Invariant #9: Exact bounds
            Assert.That(summary.LargestFragment, Is.EqualTo(9));
            Assert.That(summary.SmallestFragment, Is.EqualTo(1));
            Assert.That(summary.AverageFragmentSize, Is.EqualTo(6.0));

            // Invariant #10: Exact enzyme list
            Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "EcoRI" }));
        });
    }

    #endregion

    #region Restriction Map Tests

    /// <summary>
    /// Evidence: CreateMap returns complete map data with exact values.
    /// Invariants #12, #13: Unique cutters, site count (forward-strand only).
    /// Two EcoRI sites → not a unique cutter.
    /// </summary>
    [Test]
    public void CreateMap_ReturnsCorrectMapWithAllFields()
    {
        // Arrange: Two EcoRI sites at positions 0 and 9 (palindromic → both strands per site)
        var sequence = new DnaSequence("GAATTCAAAGAATTC"); // 15 bp

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        // Assert: Exact values for all fields
        Assert.Multiple(() =>
        {
            Assert.That(map.SequenceLength, Is.EqualTo(15));
            Assert.That(map.TotalSites, Is.EqualTo(2), "Forward-strand sites only (Invariant #13)");
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.EqualTo(new[] { 0, 0, 9, 9 }),
                "Both strands: 2 forward + 2 reverse at positions 0 and 9");
            Assert.That(map.UniqueCutters, Is.Empty, "EcoRI cuts twice — not a unique cutter");
            Assert.That(map.NonCutters, Is.Empty);
        });
    }

    /// <summary>
    /// Evidence: Implementation - Unique cutters are enzymes with exactly one forward-strand site.
    /// Invariant #12: UniqueCutters list contains enzymes with one forward-strand site.
    /// For palindromic enzymes, both strands match at the same position — only forward counts.
    /// </summary>
    [Test]
    public void CreateMap_IdentifiesUniqueCutters()
    {
        // Arrange: Sequence where each enzyme cuts only once (forward strand)
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Assert: Invariant #12 — unique cutters correctly identified for palindromic sites
        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2), "Two forward-strand sites");
            Assert.That(map.UniqueCutters, Contains.Item("EcoRI"), "EcoRI is a unique cutter");
            Assert.That(map.UniqueCutters, Contains.Item("BamHI"), "BamHI is a unique cutter");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Non-cutters are enzymes with zero sites.
    /// Invariant #11: NonCutters list contains enzymes with no sites.
    /// </summary>
    [Test]
    public void CreateMap_IdentifiesNonCutters()
    {
        // Arrange: Sequence without any recognition sites
        var sequence = new DnaSequence("AAAAAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Assert: Invariant #11
        Assert.Multiple(() =>
        {
            Assert.That(map.NonCutters, Contains.Item("EcoRI"), "EcoRI is non-cutter");
            Assert.That(map.NonCutters, Contains.Item("BamHI"), "BamHI is non-cutter");
            Assert.That(map.TotalSites, Is.EqualTo(0), "Zero total sites");
        });
    }

    /// <summary>
    /// Evidence: Implementation - TotalSites counts forward-strand sites only.
    /// Invariant #13: Avoid double-counting palindromic sites.
    /// </summary>
    [Test]
    public void CreateMap_TotalSites_CountsForwardStrandOnly()
    {
        // Arrange: Single EcoRI site (palindromic, appears on both strands)
        var sequence = new DnaSequence("AAAAAAGAATTCAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        // Assert: Invariant #13 - should count only forward strand
        Assert.That(map.TotalSites, Is.EqualTo(1),
            "Should count forward-strand sites only, not both strands");
    }

    /// <summary>
    /// Evidence: API contract - Null sequence must throw.
    /// </summary>
    [Test]
    public void CreateMap_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.CreateMap(null!, "EcoRI"));
    }

    /// <summary>
    /// Evidence: API contract — When no enzymes specified, CreateMap searches all known enzymes.
    /// Sequence contains EcoRI (pos 9), BamHI (pos 30), and GATC-cutters (Sau3AI, MboI, DpnI at pos 31).
    /// </summary>
    [Test]
    public void CreateMap_NoEnzymesSpecified_SearchesAll()
    {
        // Arrange: Known EcoRI and BamHI sites (BamHI's GGATCC contains GATC for other cutters)
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA"); // 45 bp

        // Act: No enzymes specified — should search all
        var map = RestrictionAnalyzer.CreateMap(sequence);

        // Assert: Finds multiple enzyme families
        Assert.Multiple(() =>
        {
            Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), Is.True, "Should find EcoRI site");
            Assert.That(map.SitesByEnzyme.ContainsKey("BamHI"), Is.True, "Should find BamHI site");
            Assert.That(map.SitesByEnzyme.ContainsKey("Sau3AI"), Is.True,
                "Should find Sau3AI (GATC within GGATCC)");
            Assert.That(map.TotalSites, Is.EqualTo(5),
                "5 forward-strand sites: EcoRI, BamHI, Sau3AI, MboI, DpnI");
            // RESTR-FIND-001 review (2026-10, F19): with no enzyme names the considered set is
            // the whole built-in table, so every built-in enzyme without a site IS a non-cutter
            // (TestSpec invariant #11: "Enzymes with zero sites listed in NonCutters"). The
            // previous assertion (NonCutters empty) locked an implementation artefact — the
            // filter was `enzymeNames.Where(...)` over an empty array — that contradicted it.
            Assert.That(map.NonCutters, Is.Not.Empty,
                "every built-in enzyme with no site in this sequence is a non-cutter");
            Assert.That(map.NonCutters, Does.Not.Contain("EcoRI"), "EcoRI cut the sequence");
            Assert.That(map.SitesByEnzyme.Count + map.NonCutters.Count,
                Is.EqualTo(RestrictionAnalyzer.Enzymes.Count),
                "cutters and non-cutters partition the considered set (the whole table)");
        });
    }

    #endregion

    #region Compatibility Tests

    /// <summary>
    /// Evidence: Wikipedia (Sticky and blunt ends) — "blunt ends are always compatible with each other."
    /// Invariant #14: All blunt-end enzymes are compatible.
    /// </summary>
    [Test]
    public void AreCompatible_BluntEnzymes_AreCompatible()
    {
        // EcoRV and SmaI are both blunt cutters (Wikipedia examples table: * = blunt ends)
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI");

        Assert.That(compatible, Is.True, "Blunt enzymes should be compatible (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia (Restriction enzyme) — BamHI and BglII both produce 5' GATC overhangs.
    /// BamHI: 5'---G↓GATCC---3' / 3'---CCTAG↓G---5' → 5' overhang GATC
    /// BglII: 5'---A↓GATCT---3' / 3'---TCTAG↓A---5' → 5' overhang GATC
    /// Invariant #15: Same overhang type + same sequence = compatible.
    /// </summary>
    [Test]
    public void AreCompatible_SameOverhang_AreCompatible_BamHI_BglII()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("BamHI", "BglII");

        Assert.That(compatible, Is.True, "BamHI/BglII both produce 5' GATC overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — SalI and XhoI both produce 5' TCGA overhangs.
    /// SalI:  5'---G↓TCGAC---3' / 3'---CAGCT↓G---5' → 5' overhang TCGA
    /// XhoI:  5'---C↓TCGAG---3' / 3'---GAGCT↓C---5' → 5' overhang TCGA
    /// </summary>
    [Test]
    public void AreCompatible_SameOverhang_AreCompatible_SalI_XhoI()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("SalI", "XhoI");

        Assert.That(compatible, Is.True, "SalI/XhoI both produce 5' TCGA overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — MboI, Sau3AI, BamHI, BglII all produce 5' GATC overhangs.
    /// MboI/Sau3AI: 5'---↓GATC---3' / 3'---CTAG↓---5' → 5' overhang GATC
    /// </summary>
    [Test]
    [TestCase("MboI", "BamHI")]
    [TestCase("Sau3AI", "BglII")]
    [TestCase("MboI", "Sau3AI")]
    public void AreCompatible_GATC_OverhangFamily_AllCompatible(string enzyme1, string enzyme2)
    {
        bool compatible = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);

        Assert.That(compatible, Is.True,
            $"{enzyme1}/{enzyme2} both produce 5' GATC overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — NheI, XbaI, SpeI, AvrII all produce 5' CTAG overhangs.
    /// NheI:  5'---G↓CTAGC---3' → 5' overhang CTAG
    /// XbaI:  5'---T↓CTAGA---3' → 5' overhang CTAG
    /// SpeI:  5'---A↓CTAGT---3' → 5' overhang CTAG
    /// AvrII: 5'---C↓CTAGG---3' → 5' overhang CTAG
    /// </summary>
    [Test]
    [TestCase("NheI", "XbaI")]
    [TestCase("NheI", "SpeI")]
    [TestCase("XbaI", "AvrII")]
    public void AreCompatible_CTAG_OverhangFamily_AllCompatible(string enzyme1, string enzyme2)
    {
        bool compatible = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);

        Assert.That(compatible, Is.True,
            $"{enzyme1}/{enzyme2} both produce 5' CTAG overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia (Sticky and blunt ends) — "overhangs have to be complementary
    /// in order for the ligase to work." A 5' overhang CANNOT ligate with a 3' overhang,
    /// even when the overhang sequence string is the same palindrome, because both
    /// single-stranded extensions would be on the same strand and cannot base-pair.
    ///
    /// HindIII: 5'---A↓AGCTT---3' / 3'---TTCGA↓A---5' → 5' overhang AGCT
    /// SacI:    5'---GAGCT↓C---3' / 3'---C↓TCGAG---5' → 3' overhang AGCT
    /// </summary>
    [Test]
    public void AreCompatible_CrossTypeOverhang_NotCompatible_HindIII_SacI()
    {
        // HindIII produces 5' overhang AGCT; SacI produces 3' overhang AGCT
        // Same overhang string but different types → NOT compatible
        bool compatible = RestrictionAnalyzer.AreCompatible("HindIII", "SacI");

        Assert.That(compatible, Is.False,
            "HindIII (5' AGCT) and SacI (3' AGCT) are NOT compatible — different overhang types (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — SphI produces 3' overhang CATG, NcoI produces 5' overhang CATG.
    /// SphI:  5'---GCATG↓C---3' / 3'---C↓GTACG---5' → 3' overhang CATG
    /// NcoI:  5'---C↓CATGG---3' / 3'---GGTAC↓C---5' → 5' overhang CATG
    /// Same sequence CATG but different types → NOT compatible.
    /// </summary>
    [Test]
    public void AreCompatible_CrossTypeOverhang_NotCompatible_SphI_NcoI()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("SphI", "NcoI");

        Assert.That(compatible, Is.False,
            "SphI (3' CATG) and NcoI (5' CATG) are NOT compatible — different overhang types (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — Different overhang sequences cannot ligate.
    /// EcoRI (5' AATT overhang) and PstI (3' TGCA overhang): different type AND different sequence.
    /// </summary>
    [Test]
    public void AreCompatible_DifferentOverhangs_NotCompatible()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "PstI");

        Assert.That(compatible, Is.False, "Different overhang enzymes should not be compatible");
    }

    /// <summary>
    /// Evidence: API contract — Unknown enzyme returns false (no exception).
    /// </summary>
    [Test]
    public void AreCompatible_UnknownEnzyme_ReturnsFalse()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "UnknownEnzyme");

        Assert.That(compatible, Is.False, "Unknown enzyme should return false, not throw");
    }

    /// <summary>
    /// Evidence: Mathematical property — Compatibility is symmetric.
    /// Invariant #16: AreCompatible(A, B) == AreCompatible(B, A).
    /// </summary>
    [Test]
    [TestCase("BamHI", "BglII")]
    [TestCase("EcoRV", "SmaI")]
    [TestCase("EcoRI", "PstI")]
    [TestCase("HindIII", "SacI")]
    [TestCase("SphI", "NcoI")]
    public void AreCompatible_IsSymmetric(string enzyme1, string enzyme2)
    {
        bool forward = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);
        bool reverse = RestrictionAnalyzer.AreCompatible(enzyme2, enzyme1);

        Assert.That(forward, Is.EqualTo(reverse),
            $"AreCompatible({enzyme1}, {enzyme2}) should equal AreCompatible({enzyme2}, {enzyme1})");
    }

    /// <summary>
    /// Evidence: Wikipedia (Restriction enzyme) — Verified enzyme cut positions from Examples table.
    /// EcoRI:  G↓AATTC (forward=1, reverse=5), 5' overhang AATT
    /// BamHI:  G↓GATCC (forward=1, reverse=5), 5' overhang GATC
    /// PstI:   CTGCA↓G (forward=5, reverse=1), 3' overhang TGCA
    /// SmaI:   CCC↓GGG (forward=3, reverse=3), Blunt
    /// NotI:   GC↓GGCCGC (forward=2, reverse=6), 5' overhang GGCC
    /// </summary>
    [Test]
    [TestCase("EcoRI", "GAATTC", 1, 5, false, Description = "Wikipedia: G↓AATTC, 5' overhang")]
    [TestCase("BamHI", "GGATCC", 1, 5, false, Description = "Wikipedia: G↓GATCC, 5' overhang")]
    [TestCase("HindIII", "AAGCTT", 1, 5, false, Description = "Wikipedia: A↓AGCTT, 5' overhang")]
    [TestCase("PstI", "CTGCAG", 5, 1, false, Description = "Wikipedia: CTGCA↓G, 3' overhang")]
    [TestCase("SmaI", "CCCGGG", 3, 3, true, Description = "Wikipedia: CCC↓GGG, Blunt")]
    [TestCase("EcoRV", "GATATC", 3, 3, true, Description = "Wikipedia: GAT↓ATC, Blunt")]
    [TestCase("AluI", "AGCT", 2, 2, true, Description = "Wikipedia: AG↓CT, Blunt")]
    [TestCase("HaeIII", "GGCC", 2, 2, true, Description = "Wikipedia: GG↓CC, Blunt")]
    [TestCase("NotI", "GCGGCCGC", 2, 6, false, Description = "Wikipedia: GC↓GGCCGC, 5' overhang")]
    [TestCase("TaqI", "TCGA", 1, 3, false, Description = "Wikipedia: T↓CGA, 5' overhang")]
    [TestCase("Sau3AI", "GATC", 0, 4, false, Description = "Wikipedia: ↓GATC, 5' overhang")]
    [TestCase("KpnI", "GGTACC", 5, 1, false, Description = "Wikipedia: GGTAC↓C, 3' overhang")]
    [TestCase("SacI", "GAGCTC", 5, 1, false, Description = "Wikipedia: GAGCT↓C, 3' overhang")]
    [TestCase("SphI", "GCATGC", 5, 1, false, Description = "Wikipedia: GCATG↓C, 3' overhang")]
    [TestCase("XbaI", "TCTAGA", 1, 5, false, Description = "Wikipedia: T↓CTAGA, 5' overhang")]
    [TestCase("SalI", "GTCGAC", 1, 5, false, Description = "Wikipedia: G↓TCGAC, 5' overhang")]
    [TestCase("BglII", "AGATCT", 1, 5, false, Description = "Wikipedia: A↓GATCT, 5' overhang")]
    [TestCase("SpeI", "ACTAGT", 1, 5, false, Description = "Wikipedia: A↓CTAGT, 5' overhang")]
    [TestCase("ScaI", "AGTACT", 3, 3, true, Description = "Wikipedia: AGT↓ACT, Blunt")]
    [TestCase("StuI", "AGGCCT", 3, 3, true, Description = "Wikipedia: AGG↓CCT, Blunt")]
    public void EnzymeDatabase_MatchesWikipediaData(
        string name, string recognition, int cutForward, int cutReverse, bool isBlunt)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"Enzyme {name} should exist in database");
            Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(recognition),
                $"{name} recognition sequence (Wikipedia)");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(cutForward),
                $"{name} forward cut position (Wikipedia)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(cutReverse),
                $"{name} reverse cut position (Wikipedia)");
            Assert.That(enzyme.IsBluntEnd, Is.EqualTo(isBlunt),
                $"{name} blunt end status (Wikipedia)");
        });
    }

    /// <summary>
    /// Evidence: FindCompatibleEnzymes returns known pairs.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_FindsKnownPairs()
    {
        // Act
        var compatiblePairs = RestrictionAnalyzer.FindCompatibleEnzymes().ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(compatiblePairs, Has.Count.GreaterThan(0), "Should find compatible pairs");

            // Check that BamHI/BglII pair is found (known compatible pair from Wikipedia)
            bool hasBamHIBglII = compatiblePairs.Any(p =>
                (p.Enzyme1 == "BamHI" && p.Enzyme2 == "BglII") ||
                (p.Enzyme1 == "BglII" && p.Enzyme2 == "BamHI"));
            Assert.That(hasBamHIBglII, Is.True, "BamHI/BglII should be listed as compatible");

            // HindIII (5' AGCT) and SacI (3' AGCT) must NOT appear as compatible
            bool hasHindIIISacI = compatiblePairs.Any(p =>
                (p.Enzyme1 == "HindIII" && p.Enzyme2 == "SacI") ||
                (p.Enzyme1 == "SacI" && p.Enzyme2 == "HindIII"));
            Assert.That(hasHindIIISacI, Is.False,
                "HindIII/SacI must NOT be listed — different overhang types (Wikipedia)");
        });
    }

    /// <summary>
    /// Evidence: All returned pairs should be actually compatible.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_AllReturnedPairsAreActuallyCompatible()
    {
        // Act
        var compatiblePairs = RestrictionAnalyzer.FindCompatibleEnzymes().Take(10).ToList();

        // Assert: Verify each pair
        foreach (var pair in compatiblePairs)
        {
            bool compatible = RestrictionAnalyzer.AreCompatible(pair.Enzyme1, pair.Enzyme2);
            Assert.That(compatible, Is.True,
                $"{pair.Enzyme1} and {pair.Enzyme2} should be compatible per AreCompatible()");
        }
    }

    #endregion
}
