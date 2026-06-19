using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the K-mer area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-COUNT-001 — k-mer counting (K-mer).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 32.
///
/// API under test (KmerAnalyzer.CountKmers):
///   CountKmers(sequence, k) → Dictionary&lt;kmer,count&gt; over the (|seq| − k + 1) sliding
///   windows of length k (case-insensitive; empty when k &gt; |seq|).
///
/// Relations (derived from the definition, NOT from output):
///   • INV (conservation): the counts sum to the number of length-k windows, |seq| − k + 1.
///   • MON  (occurrence subsumption): every occurrence of a (k+1)-mer w implies an
///          occurrence of its length-k prefix AND its length-k suffix at the same/adjacent
///          position, so count_{k+1}(w) ≤ count_k(prefix) and ≤ count_k(suffix); and the
///          TOTAL instance count strictly decreases by exactly 1 when k grows by 1.
///   • INV (reverse / canonical): reversing the sequence maps each k-mer occurrence to its
///          reverse, so count_s(w) = count_{reverse(s)}(reverse(w)); hence the CANONICAL
///          histogram keyed by min(w, reverse(w)) is identical for s and reverse(s).
///
/// ── Reconciliation with the checklist wording ──
///   The checklist row reads "MON: k+1 → ≤ distinct k-mers". The literal claim that the
///   number of DISTINCT k-mers is non-increasing in k is FALSE in general — e.g. "AACAAC"
///   has 2 distinct 1-mers {A,C} but 3 distinct 2-mers {AA,AC,CA}. The rigorously TRUE
///   monotone facts are the per-k-mer occurrence subsumption (a longer word cannot occur
///   more often than its prefix/suffix) and the strict −1 step in the TOTAL instance count;
///   those are what is asserted here, rather than rubber-stamping a false relation.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class KmerMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    private static string Reverse(string s) => new(s.Reverse().ToArray());

    /// <summary>Canonical histogram keyed by min(kmer, reverse(kmer)) — invariant under sequence reversal.</summary>
    private static Dictionary<string, int> CanonicalHistogram(Dictionary<string, int> counts)
    {
        var hist = new Dictionary<string, int>();
        foreach (var (kmer, count) in counts)
        {
            string canon = string.CompareOrdinal(kmer, Reverse(kmer)) <= 0 ? kmer : Reverse(kmer);
            hist[canon] = hist.GetValueOrDefault(canon) + count;
        }
        return hist;
    }

    /// <summary>Uppercase bodies (repetitive, homopolymer, and fixed-seed random).</summary>
    private static IEnumerable<string> KmerBodies()
    {
        yield return "ACGTACGTACGT";
        yield return "AAAAAAAA";
        yield return "AACAACAACAAC";
        yield return "GATTACAGATTACA";
        yield return RandomDna(60);
        yield return RandomDna(100);
    }

    private static readonly int[] KValues = { 1, 2, 3, 4, 6 };

    #endregion

    #region INV — counts sum to the number of length-k windows (|seq| − k + 1)

    [Test]
    [Description("INV: the total k-mer instance count equals |seq| − k + 1 (the number of sliding windows), and 0 when k exceeds the sequence length.")]
    public void CountKmers_TotalInstances_EqualsWindowCount()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                var counts = KmerAnalyzer.CountKmers(body, k);
                int expected = body.Length >= k ? body.Length - k + 1 : 0;

                counts.Values.Sum().Should().Be(expected,
                    because: $"every length-{k} window contributes exactly one k-mer instance, so the counts sum to |seq| − k + 1 = {expected}");
            }
        }
    }

    #endregion

    #region MON — extending k by 1 cannot increase a k-mer's count; total instances drop by 1

    [Test]
    [Description("MON: every (k+1)-mer occurs no more often than its length-k prefix or suffix (occurrence subsumption), and the total instance count strictly decreases by exactly 1 as k grows by 1.")]
    public void CountKmers_ExtendingK_SubsumesOccurrences_AndTotalDropsByOne()
    {
        foreach (var body in KmerBodies())
        {
            foreach (int k in KValues)
            {
                if (body.Length < k + 1)
                    continue;

                var ck = KmerAnalyzer.CountKmers(body, k);
                var ck1 = KmerAnalyzer.CountKmers(body, k + 1);

                foreach (var (w, count) in ck1)
                {
                    ck[w[..k]].Should().BeGreaterThanOrEqualTo(count,
                        because: $"every occurrence of the (k+1)-mer '{w}' contains an occurrence of its length-{k} prefix, so the prefix occurs at least as often");
                    ck[w[1..]].Should().BeGreaterThanOrEqualTo(count,
                        because: $"every occurrence of '{w}' contains an occurrence of its length-{k} suffix, so the suffix occurs at least as often");
                }

                ck1.Values.Sum().Should().Be(ck.Values.Sum() - 1,
                    because: "there is exactly one fewer length-(k+1) window than length-k window, so the total instance count drops by 1");
            }
        }
    }

    #endregion

    #region INV — reversing the sequence maps each k-mer to its reverse (canonical histogram preserved)

    [Test]
    [Description("INV: reversing the sequence maps each k-mer to its reverse with the same count, so the canonical histogram (keyed by min(kmer, reverse(kmer))) is identical.")]
    public void CountKmers_ReverseSequence_PreservesCanonicalHistogram()
    {
        foreach (var body in KmerBodies())
        {
            string reversed = Reverse(body);

            foreach (int k in KValues)
            {
                var forward = KmerAnalyzer.CountKmers(body, k);
                var rev = KmerAnalyzer.CountKmers(reversed, k);

                rev.Count.Should().Be(forward.Count,
                    because: "the bijection kmer ↔ reverse(kmer) pairs distinct k-mers, so the number of distinct k-mers is preserved under reversal");

                foreach (var (w, count) in forward)
                {
                    rev.GetValueOrDefault(Reverse(w)).Should().Be(count,
                        because: $"every occurrence of '{w}' in s corresponds to an occurrence of '{Reverse(w)}' in reverse(s), so their counts match");
                }

                CanonicalHistogram(rev).Should().BeEquivalentTo(CanonicalHistogram(forward),
                    because: "min(kmer, reverse(kmer)) is invariant under reversal, so the canonical k-mer histogram is the same for s and reverse(s)");
            }
        }
    }

    #endregion
}
