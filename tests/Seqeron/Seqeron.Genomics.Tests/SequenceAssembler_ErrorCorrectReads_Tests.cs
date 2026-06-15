// ASSEMBLY-CORRECT-001 — K-mer spectrum (two-sided) read error correction
// Evidence: docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-CORRECT-001.md
// Source: Liu, Schmidt & Maskell (2013). Musket. Bioinformatics 29(3):308-315.
//         Kelley, Schatz & Salzberg (2010). Quake. Genome Biol 11:R116.

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

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
}
