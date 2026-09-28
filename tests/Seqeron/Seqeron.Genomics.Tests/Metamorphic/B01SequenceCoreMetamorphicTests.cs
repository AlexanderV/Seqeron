namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the behaviours introduced by review batch B01
/// (docs/Validation/review-2026-09/B01.md): IUPAC validation / IndexOfInvalid*, the linguistic-complexity
/// suffix-tree path (F4), Shannon entropy over U (F5), GC-skew/replication-origin relations (F6) and
/// AnalyzeGcContent with RNA (F7). Fixed-seed random inputs; every relation is exact.
///
/// Test Units: SEQ-VALID-001, SEQ-COMPLEX-001, SEQ-ENTROPY-001, SEQ-GCSKEW-001, SEQ-REPLICATION-001,
/// SEQ-GC-ANALYSIS-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Composition")]
public class B01SequenceCoreMetamorphicTests
{
    private const int Seed = 20260928;

    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static string ReverseComplementAcgt(string s) =>
        new(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' }).ToArray());

    #region SEQ-VALID-001 — concatenation relations

    /// <summary>
    /// MR: validity of x+y is the conjunction of validities, and IndexOfInvalid(x+y) =
    /// IndexOfInvalid(x) if x is invalid, else |x| + IndexOfInvalid(y) (or −1).
    /// </summary>
    [Test]
    public void Validation_Concatenation_ComposesExactly()
    {
        var rng = new Random(Seed);
        const string alphabet = "ACGTUacgtuNRYnry-.Xſ";
        for (int trial = 0; trial < 500; trial++)
        {
            string x = Random(rng, alphabet, rng.Next(0, 12));
            string y = Random(rng, alphabet, rng.Next(0, 12));
            string xy = x + y;

            xy.AsSpan().IsValidIupacDna().Should().Be(x.AsSpan().IsValidIupacDna() && y.AsSpan().IsValidIupacDna(), xy);
            xy.AsSpan().IsValidIupacRna().Should().Be(x.AsSpan().IsValidIupacRna() && y.AsSpan().IsValidIupacRna(), xy);

            static int Compose(int ix, int iy, int lenX) => ix >= 0 ? ix : iy >= 0 ? lenX + iy : -1;
            xy.AsSpan().IndexOfInvalidDna().Should().Be(
                Compose(x.AsSpan().IndexOfInvalidDna(), y.AsSpan().IndexOfInvalidDna(), x.Length), xy);
            xy.AsSpan().IndexOfInvalidRna().Should().Be(
                Compose(x.AsSpan().IndexOfInvalidRna(), y.AsSpan().IndexOfInvalidRna(), x.Length), xy);
        }
    }

    /// <summary>
    /// MR: transcription T→U maps a valid DNA string onto a valid RNA string (and IUPAC-DNA onto IUPAC-RNA),
    /// and DnaSequence.Transcribe agrees with the RnaSequence validity gate.
    /// </summary>
    [Test]
    public void Validation_Transcription_MapsDnaValidityToRnaValidity()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 200; trial++)
        {
            string dna = Random(rng, "ACGTacgt", rng.Next(1, 40));
            string rna = dna.Replace('T', 'U').Replace('t', 'u');
            dna.AsSpan().IsValidDna().Should().BeTrue();
            rna.AsSpan().IsValidRna().Should().BeTrue(rna);
            rna.AsSpan().IsValidDna().Should().Be(!dna.Contains('T') && !dna.Contains('t'), rna);

            string iupac = Random(rng, "ACGTRYSWKMBDHVN", rng.Next(1, 40));
            iupac.AsSpan().IsValidIupacDna().Should().BeTrue();
            iupac.Replace('T', 'U').AsSpan().IsValidIupacRna().Should().BeTrue();

            new DnaSequence(dna).Transcribe().Should().Be(rna.ToUpperInvariant());
        }
    }

    #endregion

    #region SEQ-COMPLEX-001 — suffix-tree path relations on long sequences

    /// <summary>
    /// MR: all-length LC (m = N, suffix-tree path) is invariant under reversal, reverse complement and ACGT
    /// relabeling; string and DnaSequence overloads agree. Long inputs (500–3000 nt) only reachable in
    /// practice through the linear-time suffix-tree path.
    /// </summary>
    [Test]
    public void LinguisticComplexity_LongSequences_SymmetryRelations()
    {
        var rng = new Random(Seed + 2);
        foreach (int n in new[] { 500, 1200, 3000 })
        {
            string s = Random(rng, "ACGT", n);
            // Low-complexity variant: tandem repeat with a few point mutations.
            char[] rep = string.Concat(Enumerable.Repeat("ACGTTG", n / 6 + 1)).Substring(0, n).ToCharArray();
            for (int k = 0; k < 5; k++) rep[rng.Next(n)] = "ACGT"[rng.Next(4)];
            foreach (string seq in new[] { s, new string(rep) })
            {
                double lc = SequenceComplexity.CalculateLinguisticComplexity(seq, int.MaxValue);
                lc.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(1.0);
                SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(seq), int.MaxValue).Should().Be(lc);
                SequenceComplexity.CalculateLinguisticComplexity(new string(seq.Reverse().ToArray()), int.MaxValue).Should().Be(lc);
                SequenceComplexity.CalculateLinguisticComplexity(ReverseComplementAcgt(seq), int.MaxValue).Should().Be(lc);
                string relabeled = new(seq.Select(c => c switch { 'A' => 'G', 'G' => 'T', 'T' => 'C', _ => 'A' }).ToArray());
                SequenceComplexity.CalculateLinguisticComplexity(relabeled, int.MaxValue).Should().Be(lc);
                SequenceComplexity.CalculateLinguisticComplexity(seq.ToLowerInvariant(), int.MaxValue).Should().Be(lc);
            }
        }
    }

    /// <summary>
    /// MR: for m ≥ N the result no longer depends on m (V_i = 0 does not exist beyond N; m is capped at N),
    /// so m = N, N+1, and int.MaxValue give the same value; and m = 12 (hash) vs m = 13 (suffix tree)
    /// differ exactly by adding V_13 and V_max,13 to the sums.
    /// </summary>
    [Test]
    public void LinguisticComplexity_PathBoundary_ConsistentSums()
    {
        var rng = new Random(Seed + 3);
        for (int trial = 0; trial < 40; trial++)
        {
            int n = rng.Next(14, 200);
            string s = trial % 2 == 0 ? Random(rng, "ACGT", n) : Random(rng, "AC", n);
            double atN = SequenceComplexity.CalculateLinguisticComplexity(s, n);
            SequenceComplexity.CalculateLinguisticComplexity(s, n + 1).Should().Be(atN);
            SequenceComplexity.CalculateLinguisticComplexity(s, int.MaxValue).Should().Be(atN);

            // Independent V_i for i ≤ 13 and the closed-form V_max,i.
            long obs12 = 0, pos12 = 0;
            for (int i = 1; i <= 12; i++)
            {
                obs12 += Enumerable.Range(0, n - i + 1).Select(j => s.Substring(j, i)).Distinct().Count();
                pos12 += Math.Min(1L << (2 * i), n - i + 1);
            }
            long v13 = Enumerable.Range(0, n - 13 + 1).Select(j => s.Substring(j, 13)).Distinct().Count();
            long max13 = Math.Min(1L << 26, n - 13 + 1);

            SequenceComplexity.CalculateLinguisticComplexity(s, 12).Should().Be((double)obs12 / pos12);
            SequenceComplexity.CalculateLinguisticComplexity(s, 13).Should().Be((double)(obs12 + v13) / (pos12 + max13));
        }
    }

    #endregion

    #region SEQ-ENTROPY-001 — Shannon entropy relations

    /// <summary>
    /// MR: H is invariant under T→U transcription, case change, repetition (s·k), shuffling, and insertion of
    /// non-nucleotide symbols; always within [0, 2].
    /// </summary>
    [Test]
    public void ShannonEntropy_InvarianceRelations()
    {
        var rng = new Random(Seed + 4);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Random(rng, trial % 3 == 0 ? "AT" : "ACGT", rng.Next(1, 80));
            double h = SequenceComplexity.CalculateShannonEntropy(s);
            h.Should().BeInRange(0, 2.0 + 1e-12);

            SequenceComplexity.CalculateShannonEntropy(s.Replace('T', 'U')).Should().Be(h, s);
            SequenceComplexity.CalculateShannonEntropy(s.Replace('T', 'U').ToLowerInvariant()).Should().Be(h, s);
            SequenceComplexity.CalculateShannonEntropy(s + s + s).Should().BeApproximately(h, 1e-12, s);
            string shuffled = new(s.OrderBy(_ => rng.Next()).ToArray());
            SequenceComplexity.CalculateShannonEntropy(shuffled).Should().BeApproximately(h, 1e-12, s);
            string noisy = string.Concat(s.Select(c => rng.Next(4) == 0 ? c + "N-R" : c.ToString()));
            SequenceComplexity.CalculateShannonEntropy(noisy).Should().Be(h, noisy);
            SequenceComplexity.CalculateShannonEntropy(new DnaSequence(s)).Should().Be(h, s);
        }
    }

    #endregion

    #region SEQ-GCSKEW-001 / SEQ-REPLICATION-001 / SEQ-GC-ANALYSIS-001

    /// <summary>
    /// MR: prepending k·w bases of A/T (skew-neutral, whole windows) keeps every cumulative value, shifts the
    /// cumulative profile by k windows (preceded by k zero points), and shifts a non-zero replication origin by k·w.
    /// </summary>
    [Test]
    public void GcSkew_PrependNeutralWindows_ShiftsProfileAndOrigin()
    {
        var rng = new Random(Seed + 5);
        for (int trial = 0; trial < 100; trial++)
        {
            int w = rng.Next(1, 10);
            int k = rng.Next(1, 5);
            string s = Random(rng, "ACGT", rng.Next(0, 150));
            string flank = Random(rng, "AT", k * w);
            string t = flank + s;

            var a = GcSkewCalculator.CalculateCumulativeGcSkew(s, w).ToList();
            var b = GcSkewCalculator.CalculateCumulativeGcSkew(t, w).ToList();
            b.Count.Should().Be(a.Count + k);
            b.Take(k).Should().OnlyContain(p => p.GcSkew == 0 && p.CumulativeGcSkew == 0);
            for (int i = 0; i < a.Count; i++)
            {
                b[k + i].CumulativeGcSkew.Should().Be(a[i].CumulativeGcSkew);
                b[k + i].Position.Should().Be(a[i].Position + k * w);
            }

            var oa = GcSkewCalculator.PredictReplicationOrigin(s);
            var ob = GcSkewCalculator.PredictReplicationOrigin(t);
            ob.OriginSkew.Should().Be(oa.OriginSkew);
            ob.TerminusSkew.Should().Be(oa.TerminusSkew);
            ob.PredictedOrigin.Should().Be(oa.PredictedOrigin == 0 ? 0 : oa.PredictedOrigin + k * w);
            ob.PredictedTerminus.Should().Be(oa.PredictedTerminus == 0 ? 0 : oa.PredictedTerminus + k * w);
        }
    }

    /// <summary>
    /// MR: swapping G↔C negates every windowed and cumulative GC skew and exchanges origin and terminus
    /// (Skew_i → −Skew_i turns the first argmin into the first argmax).
    /// </summary>
    [Test]
    public void GcSkew_GcSwap_NegatesProfileAndExchangesOriginTerminus()
    {
        var rng = new Random(Seed + 6);
        for (int trial = 0; trial < 100; trial++)
        {
            string s = Random(rng, "ACGT", rng.Next(1, 150));
            string swapped = new(s.Select(c => c == 'G' ? 'C' : c == 'C' ? 'G' : c).ToArray());
            int w = rng.Next(1, 10);

            var a = GcSkewCalculator.CalculateCumulativeGcSkew(s, w).Select(p => p.CumulativeGcSkew).ToList();
            var b = GcSkewCalculator.CalculateCumulativeGcSkew(swapped, w).Select(p => p.CumulativeGcSkew).ToList();
            b.Should().Equal(a.Select(x => -x));

            var oa = GcSkewCalculator.PredictReplicationOrigin(s);
            var ob = GcSkewCalculator.PredictReplicationOrigin(swapped);
            ob.PredictedOrigin.Should().Be(oa.PredictedTerminus);
            ob.PredictedTerminus.Should().Be(oa.PredictedOrigin);
            ob.OriginSkew.Should().Be(-oa.TerminusSkew);
        }
    }

    /// <summary>
    /// MR: AnalyzeGcContent GC metrics are invariant under T→U transcription and case change; the
    /// DnaSequence and string overloads agree on DNA input.
    /// </summary>
    [Test]
    public void AnalyzeGcContent_TranscriptionAndCase_PreserveGcMetrics()
    {
        var rng = new Random(Seed + 7);
        for (int trial = 0; trial < 100; trial++)
        {
            string dna = Random(rng, "ACGT", rng.Next(1, 200));
            int w = rng.Next(1, 30), step = rng.Next(1, 30);
            var d = GcSkewCalculator.AnalyzeGcContent(dna, w, step);
            var viaSeq = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(dna), w, step);
            foreach (var r in new[]
                     {
                         GcSkewCalculator.AnalyzeGcContent(dna.Replace('T', 'U'), w, step),
                         GcSkewCalculator.AnalyzeGcContent(dna.Replace('T', 'U').ToLowerInvariant(), w, step),
                         viaSeq,
                     })
            {
                r.OverallGcContent.Should().Be(d.OverallGcContent);
                r.OverallGcSkew.Should().Be(d.OverallGcSkew);
                r.GcContentVariance.Should().Be(d.GcContentVariance);
                r.WindowedGcContent.Select(x => x.GcContent).Should().Equal(d.WindowedGcContent.Select(x => x.GcContent));
            }
            viaSeq.OverallAtSkew.Should().Be(d.OverallAtSkew);
        }
    }

    #endregion
}
