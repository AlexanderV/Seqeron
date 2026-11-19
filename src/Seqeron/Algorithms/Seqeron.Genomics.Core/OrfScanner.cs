namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// One open reading frame found on a single strand by <see cref="OrfScanner"/>.
    /// </summary>
    /// <param name="Start">0-based index of the first base of the first codon (the START codon), in the scanned strand.</param>
    /// <param name="End">0-based exclusive end: the index after the last base of the ORF. For a terminated ORF the
    /// last codon is the STOP codon (<c>End − 3</c> is its first base); for an open-ended ORF it is the last
    /// complete codon before the strand end. <c>(End − Start)</c> is always a multiple of three.</param>
    /// <param name="IsTerminated">True when the ORF ends with a STOP codon; false when it runs to the strand end.</param>
    public readonly record struct OrfRegion(int Start, int End, bool IsTerminated)
    {
        /// <summary>Length in nucleotides, including the STOP codon of a terminated ORF.</summary>
        public int Length => End - Start;

        /// <summary>Reading-frame offset of the ORF within its strand (0, 1 or 2).</summary>
        public int FrameOffset => Start % 3;
    }

    /// <summary>
    /// Single-pass, linear-time scanner for open reading frames (ORFs) on one strand: the shared core of the
    /// START→STOP ORF finders. Policy choices that distinguish the published ORF definitions are parameters:
    /// which codons are STARTs, whether nested in-frame STARTs are reported, and whether an ORF that runs to the
    /// strand end without a STOP is reported.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Model. For every frame offset f ∈ {0, 1, 2} the strand is read in non-overlapping codons
    /// S[f..f+3), S[f+3..f+6), … (a trailing partial codon is ignored). A codon is a STOP when the
    /// genetic code translates it as '*' (<see cref="GeneticCode.IsStopCodon(string)"/>: an IUPAC-ambiguous codon
    /// only if every concrete codon it stands for is a stop, e.g. TAR) and a START when the supplied predicate
    /// accepts it (only unambiguous A/C/G/T codons are ever STARTs, as in EMBOSS <c>ajTrnCodonstrTypeC</c> and
    /// NCBI <c>COrf</c>). An ORF is a START codon through the first in-frame STOP codon, inclusive.
    /// </para>
    /// <para>
    /// Nested STARTs. With <c>nestedStarts = true</c> every START that is followed by an in-frame STOP yields an ORF,
    /// so in-frame STARTs sharing one STOP are all reported (Rosalind ORF "all distinct candidate proteins";
    /// NCBI <c>COrf::FindOrfs(longest_orfs = false)</c>). With <c>false</c> only the most upstream START after the
    /// previous STOP (or the frame start) opens an ORF (NCBI ORFfinder default "longest ORF per STOP"; EMBOSS
    /// getorf <c>-find 1</c>: "nested in-frame STARTs do not open a new ORF").
    /// </para>
    /// <para>
    /// Open ends. With <c>includeOpenEnded = true</c> the STARTs still pending at the end of a frame are reported as
    /// open (3'-partial) ORFs ending at the last complete codon (<see cref="OrfRegion.IsTerminated"/> = false) —
    /// NCBI ORFfinder marks them partial, EMBOSS getorf <c>WriteORF(start, pos+2)</c>. Otherwise they are dropped
    /// (a reading without a STOP is not a complete ORF; Rosalind translates "until a stop codon").
    /// </para>
    /// <para>
    /// Complexity: O(n) time per frame, O(number of pending STARTs) extra space (cleared at every STOP), plus the
    /// output. Regions are yielded frame by frame (offset 0, 1, 2), in ascending start position within a frame.
    /// </para>
    /// </remarks>
    public static class OrfScanner
    {
        private const int CodonLength = 3;
        private const int FramesPerStrand = 3;
        private const string Bases = "ACGT";

        /// <summary>
        /// Finds the ORFs of one strand.
        /// </summary>
        /// <param name="strand">The strand to scan, a nucleotide string (A/C/G/T/U, either case; other IUPAC symbols make
        /// a codon ambiguous: it can end an ORF only if every concrete codon is a stop, and never starts one).</param>
        /// <param name="geneticCode">The genetic code whose STOP codons end ORFs (default: Standard, NCBI table 1).</param>
        /// <param name="isStartCodon">Predicate on an upper-case DNA codon spelled with A/C/G/T (e.g. <c>"ATG"</c>) that
        /// selects the START codons; <see langword="null"/> makes every non-STOP (sense) codon a START
        /// (NCBI ORFfinder "any sense codon": stop-to-stop regions).</param>
        /// <param name="nestedStarts">Report every in-frame START (true) or only the most upstream one per STOP (false).</param>
        /// <param name="includeOpenEnded">Also report ORFs that reach the strand end without a STOP codon.</param>
        /// <param name="allowDualCodingStops">false (default): a genetic code whose stop codons also code for an amino acid
        /// (NCBI tables 27, 28, 31) is rejected. true: the codons the table lists as stops (<see cref="GeneticCode.StopCodons"/>)
        /// end ORFs regardless (NCBI ORFfinder 0.4.3 behaviour, measured on the binary for <c>-g 27/28/31</c>).</param>
        /// <param name="cancellationToken">Polled while the (lazy) scan advances, every 65 536 codons; a cancelled token throws
        /// <see cref="OperationCanceledException"/> from the enumeration.</param>
        /// <returns>The ORFs, frame offset 0 first, each frame in ascending start order.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="strand"/> is null.</exception>
        /// <exception cref="ArgumentException">The genetic code has dual-coding stop codons (NCBI tables 27, 28, 31) and
        /// <paramref name="allowDualCodingStops"/> is false: no codon unambiguously ends an ORF there (see
        /// <see cref="Translator.FindOrfs(DnaSequence, GeneticCode?, int, bool)"/>).</exception>
        public static IEnumerable<OrfRegion> FindRegions(string strand, GeneticCode? geneticCode = null,
            Func<string, bool>? isStartCodon = null, bool nestedStarts = true, bool includeOpenEnded = false,
            bool allowDualCodingStops = false, CancellationToken cancellationToken = default)
        {
            ArgumentNullException.ThrowIfNull(strand);
            var code = geneticCode ?? GeneticCode.Standard;
            if (!allowDualCodingStops && Translator.HasDualCodingStopCodons(code))
            {
                throw new ArgumentException(
                    $"ORF finding cannot be used with genetic code table {code.TableNumber} " +
                    "because its stop codons also code for an amino acid (context-dependent termination).",
                    nameof(geneticCode));
            }

            var (stopTable, startTable) = BuildCodonTables(code, isStartCodon);
            return Scan(strand, code, stopTable, startTable, nestedStarts, includeOpenEnded, cancellationToken);
        }

        /// <summary>
        /// The location of an ORF for the cross-frame nesting filter <see cref="FindNestedSpans"/>: <c>Start</c> and <c>End</c> are
        /// the 0-based inclusive coordinates of the ORF on the <b>forward</b> strand (stop codon included for a terminated ORF).
        /// </summary>
        public readonly record struct OrfSpan(int Start, int End, bool IsReverseStrand);

        /// <summary>
        /// NCBI ORFfinder 0.4.3 "Ignore nested ORFs (completely placed within another)" (<c>-n true</c>), measured on the binary
        /// (A1-18): an ORF is <b>nested</b> when its forward-strand span <c>[Start, End]</c> lies completely inside the span of
        /// <i>another</i> ORF of the list — any reading frame, either strand, same or different stop codon (non-strict
        /// containment: a shared start or a shared end counts). Two ORFs with exactly the same span can only come from opposite
        /// strands; then the <b>plus-strand ORF is nested</b> and the minus-strand ORF is kept. The reference ORFs are the ORFs of the
        /// list itself (those that passed the length filter and strand choice), so containment is transitive and the result does not
        /// depend on the order of the list.
        /// </summary>
        /// <param name="spans">The ORFs' forward-strand spans.</param>
        /// <returns>For each input span, true when it is nested in another one and must be dropped. O(n log n).</returns>
        public static bool[] FindNestedSpans(IReadOnlyList<OrfSpan> spans)
        {
            ArgumentNullException.ThrowIfNull(spans);
            var nested = new bool[spans.Count];
            var order = new int[spans.Count];
            for (int i = 0; i < order.Length; i++)
            {
                order[i] = i;
            }

            // Start ascending, End descending, minus strand before plus strand on an identical span: every ORF that can contain
            // another one is visited first, and an identical plus-strand span is nested in the minus-strand one.
            Array.Sort(order, (a, b) =>
            {
                var x = spans[a];
                var y = spans[b];
                int c = x.Start.CompareTo(y.Start);
                if (c != 0) return c;
                c = y.End.CompareTo(x.End);
                if (c != 0) return c;
                return y.IsReverseStrand.CompareTo(x.IsReverseStrand);
            });

            int maxEnd = int.MinValue;
            foreach (int i in order)
            {
                if (spans[i].End <= maxEnd)
                {
                    nested[i] = true;
                }
                else
                {
                    maxEnd = spans[i].End;
                }
            }

            return nested;
        }

        // 64-entry lookup (index = 16·b1 + 4·b2 + b3, bases in ACGT order) of the codon classes, built once per call from
        // the canonical GeneticCode (no private stop/start table).
        private static (bool[] Stop, bool[] Start) BuildCodonTables(GeneticCode code, Func<string, bool>? isStartCodon)
        {
            var stop = new bool[64];
            var start = new bool[64];
            for (int i = 0; i < 64; i++)
            {
                string codon = string.Concat(Bases[i >> 4], Bases[(i >> 2) & 3], Bases[i & 3]);
                stop[i] = code.IsStopCodon(codon);
                start[i] = !stop[i] && (isStartCodon is null || isStartCodon(codon));
            }

            return (stop, start);
        }

        private static IEnumerable<OrfRegion> Scan(string seq, GeneticCode code, bool[] stopTable, bool[] startTable,
            bool nestedStarts, bool includeOpenEnded, CancellationToken cancellationToken)
        {
            var pending = new List<int>();
            int polled = 0;
            for (int frame = 0; frame < FramesPerStrand; frame++)
            {
                pending.Clear();
                int i = frame;
                for (; i + CodonLength <= seq.Length; i += CodonLength)
                {
                    if ((++polled & 0xFFFF) == 0)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                    }

                    int index = CodonIndex(seq, i);
                    bool isStop = index >= 0 ? stopTable[index] : IsAmbiguousStop(code, seq, i);
                    if (isStop)
                    {
                        // Starts in ascending order, so the first (outermost) one is the longest ORF of this STOP.
                        foreach (int start in pending)
                        {
                            yield return new OrfRegion(start, i + CodonLength, true);
                        }

                        pending.Clear();
                    }
                    else if (index >= 0 && startTable[index] && (nestedStarts || pending.Count == 0))
                    {
                        pending.Add(i);
                    }
                }

                // i is now the first codon position that no longer fits: open ORFs end at the last complete codon.
                if (includeOpenEnded)
                {
                    foreach (int start in pending)
                    {
                        yield return new OrfRegion(start, i, false);
                    }
                }
            }
        }

        // Index of an A/C/G/T/U codon in the 64-entry tables, −1 when any base is another symbol.
        private static int CodonIndex(string seq, int i)
        {
            int a = BaseCode(seq[i]);
            int b = BaseCode(seq[i + 1]);
            int c = BaseCode(seq[i + 2]);
            return (a | b | c) < 0 ? -1 : (a << 4) | (b << 2) | c;
        }

        private static int BaseCode(char c) => c switch
        {
            'A' or 'a' => 0,
            'C' or 'c' => 1,
            'G' or 'g' => 2,
            'T' or 't' or 'U' or 'u' => 3,
            _ => -1,
        };

        // An IUPAC-ambiguous codon ends an ORF only when every concrete codon is a stop (e.g. TAR, TRA in table 1);
        // symbols outside IUPAC are never stops.
        private static bool IsAmbiguousStop(GeneticCode code, string seq, int i)
        {
            string codon = seq.Substring(i, CodonLength);
            foreach (char c in codon)
            {
                if (!IupacHelper.IsNucleotideCode(char.ToUpperInvariant(c)))
                {
                    return false;
                }
            }

            return code.IsStopCodon(codon);
        }
    }
}
