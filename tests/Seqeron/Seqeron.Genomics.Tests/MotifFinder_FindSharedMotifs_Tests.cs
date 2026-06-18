// MOTIF-SHARED-001 — Shared Motifs via fixed-length word enumeration with matching-sequence quorum
// Evidence: docs/Evidence/MOTIF-SHARED-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-SHARED-001.md
// Source: RSAT oligo-analysis manual (matching sequences = number of input sequences containing
//         >=1 occurrence of the oligonucleotide); Das & Dai (2007) BMC Bioinformatics 8(S7):S21
//         (word-enumeration family, exact matching, quorum across sequences);
//         van Helden, Andre, Collado-Vides (1998) J Mol Biol 281(5):827-842.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical test class for MOTIF-SHARED-001: shared fixed-length words across multiple DNA
/// sequences. Verifies <see cref="MotifFinder.FindSharedMotifs(System.Collections.Generic.IEnumerable{DnaSequence}, int, int)"/>
/// against the RSAT "matching sequences" definition (number of input sequences containing at least
/// one exact occurrence) and the word-enumeration quorum (Das &amp; Dai 2007).
/// </summary>
[TestFixture]
public class MotifFinder_FindSharedMotifs_Tests
{
    // Hand-traced dataset (Evidence §Test Datasets):
    //   S0 = ATGATG (ATG at 0 and 3), S1 = ATGCCC (ATG at 0, CCC at 3), S2 = CCCGGG (CCC at 0, GGG at 3)
    private static DnaSequence[] TraceSet() => new[]
    {
        new DnaSequence("ATGATG"),
        new DnaSequence("ATGCCC"),
        new DnaSequence("CCCGGG"),
    };

    #region FindSharedMotifs — MUST

    // M1 — "ATG" occurs in S0 and S1 -> matching sequences {0,1}. RSAT matching-sequence definition.
    [Test]
    public void FindSharedMotifs_QuorumWord_ReturnsExactSequenceIndices()
    {
        var seqs = TraceSet();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var atg = shared.Single(m => m.Sequence == "ATG");
        Assert.That(atg.SequenceIndices.OrderBy(i => i), Is.EqualTo(new[] { 0, 1 }),
            "\"ATG\" occurs in sequences 0 (ATGATG) and 1 (ATGCCC); its matching-sequence set is exactly {0,1} per the RSAT matching-sequence definition.");
    }

    // M2 — "CCC" occurs in S1 and S2 -> matching sequences {1,2}.
    [Test]
    public void FindSharedMotifs_SecondQuorumWord_ReturnsExactSequenceIndices()
    {
        var seqs = TraceSet();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var ccc = shared.Single(m => m.Sequence == "CCC");
        Assert.That(ccc.SequenceIndices.OrderBy(i => i), Is.EqualTo(new[] { 1, 2 }),
            "\"CCC\" occurs in sequences 1 (ATGCCC) and 2 (CCCGGG); matching-sequence set is exactly {1,2}.");
    }

    // M3 — "ATG" appears twice within S0 but contributes 1 to its matching-sequence count.
    // RSAT: "at least one occurrence" / "only the first occurrence of each sequence is taken into consideration".
    [Test]
    public void FindSharedMotifs_WordRepeatedInOneSequence_CountsThatSequenceOnce()
    {
        var seqs = TraceSet();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var atg = shared.Single(m => m.Sequence == "ATG");
        Assert.Multiple(() =>
        {
            Assert.That(atg.SequenceIndices.Count(i => i == 0), Is.EqualTo(1),
                "ATG occurs at positions 0 and 3 in sequence 0, but matching-sequence counting records sequence 0 exactly once.");
            Assert.That(atg.SequenceIndices.Count, Is.EqualTo(2),
                "Total matching sequences for ATG is 2 (sequences 0 and 1), not the 3 raw occurrences.");
        });
    }

    // M4 — "GGG" occurs only in S2 -> below quorum (1 < 2) -> excluded. Quorum criterion (Das & Dai).
    [Test]
    public void FindSharedMotifs_WordBelowQuorum_IsExcluded()
    {
        var seqs = TraceSet();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        Assert.That(shared.Any(m => m.Sequence == "GGG"), Is.False,
            "\"GGG\" occurs only in sequence 2 (matching sequences = 1), which is below the quorum of 2, so it must not be reported.");
    }

    // M5 — Prevalence = matchingSequences / totalSequences = 2/3 exactly.
    [Test]
    public void FindSharedMotifs_Prevalence_EqualsMatchingSequencesOverTotal()
    {
        var seqs = TraceSet();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var atg = shared.Single(m => m.Sequence == "ATG");
        Assert.That(atg.Prevalence, Is.EqualTo(2.0 / 3.0).Within(1e-10),
            "ATG matches 2 of 3 input sequences, so prevalence = 2/3 (INV-04).");
    }

    // M6 — Exact matching: ACGT vs ACTT differ at position 2, so no length-4 word is shared.
    // Das & Dai: "no variations allowed within an oligonucleotide" (INV-05).
    [Test]
    public void FindSharedMotifs_ExactMatchingOnly_NoSharedWordForSubstitutedSequences()
    {
        var seqs = new[] { new DnaSequence("ACGT"), new DnaSequence("ACTT") };

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 4, minSequences: 2).ToList();

        Assert.That(shared, Is.Empty,
            "ACGT and ACTT differ at one position; with exact (non-degenerate) matching they share no length-4 word, so no shared motif is reported.");
    }

    #endregion

    #region FindSharedMotifs — SHOULD

    // S1 — A sequence shorter than k yields no length-k window, so no 3-mer can be shared with it.
    [Test]
    public void FindSharedMotifs_SequenceShorterThanK_ContributesNoWords()
    {
        var seqs = new[] { new DnaSequence("AT"), new DnaSequence("ATGAAT") };

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        Assert.That(shared, Is.Empty,
            "Sequence 0 (\"AT\", length 2) has no length-3 window, so it shares no 3-mer with sequence 1; nothing reaches the quorum of 2.");
    }

    // S2 — Quorum equal to the full set: ATG present in all three -> reported with all three indices.
    [Test]
    public void FindSharedMotifs_FullQuorum_ReturnsWordPresentInAllSequences()
    {
        var seqs = new[]
        {
            new DnaSequence("ATGAAA"),
            new DnaSequence("CCATGC"),
            new DnaSequence("GGGATG"),
        };

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 3).ToList();

        var atg = shared.Single(m => m.Sequence == "ATG");
        Assert.Multiple(() =>
        {
            Assert.That(atg.SequenceIndices.OrderBy(i => i), Is.EqualTo(new[] { 0, 1, 2 }),
                "ATG occurs in all three sequences, so at a full quorum its matching-sequence set is {0,1,2}.");
            Assert.That(atg.Prevalence, Is.EqualTo(1.0).Within(1e-10),
                "ATG matches all 3 of 3 sequences -> prevalence = 1.0.");
        });
    }

    // S3 — INV-01: every reported word has length exactly k.
    [Test]
    public void FindSharedMotifs_AllReportedWords_HaveLengthK()
    {
        var seqs = TraceSet();
        const int k = 3;

        var shared = MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: 2).ToList();

        Assert.That(shared, Is.Not.Empty, "Trace set must yield at least one shared 3-mer.");
        Assert.That(shared.All(m => m.Sequence.Length == k), Is.True,
            "Only fixed-length-k windows are enumerated, so every reported word has length exactly k (INV-01).");
    }

    #endregion

    #region FindSharedMotifs — COULD (validation / edge)

    // C1 — Empty collection -> no shared motifs.
    [Test]
    public void FindSharedMotifs_EmptyCollection_ReturnsEmpty()
    {
        var shared = MotifFinder.FindSharedMotifs(Array.Empty<DnaSequence>(), k: 3, minSequences: 2).ToList();

        Assert.That(shared, Is.Empty,
            "With no input sequences there are no words to enumerate, so the result is empty.");
    }

    // C2 — k < 1 -> ArgumentOutOfRangeException (contract).
    [Test]
    public void FindSharedMotifs_KLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MotifFinder.FindSharedMotifs(TraceSet(), k: 0, minSequences: 2).ToList(),
            "Word length k must be >= 1; k = 0 is invalid.");
    }

    // C3 — null collection -> ArgumentNullException (contract).
    [Test]
    public void FindSharedMotifs_NullCollection_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => MotifFinder.FindSharedMotifs(null!, k: 3, minSequences: 2).ToList(),
            "A null sequence collection is invalid input.");
    }

    // C2b — minSequences < 1 -> ArgumentOutOfRangeException (contract; quorum must be >= 1).
    [Test]
    public void FindSharedMotifs_MinSequencesLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MotifFinder.FindSharedMotifs(TraceSet(), k: 3, minSequences: 0).ToList(),
            "The matching-sequence quorum must be >= 1; minSequences = 0 is invalid.");
    }

    #endregion
}
