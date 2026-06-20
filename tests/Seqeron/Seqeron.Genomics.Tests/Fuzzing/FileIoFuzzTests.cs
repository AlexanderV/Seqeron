using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the FileIO area — FASTA parsing (PARSE-FASTA-001) and FASTQ
/// parsing (PARSE-FASTQ-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, truncated, random and injected input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no out-of-bounds indexing leaking from an internal Substring, no
/// NullReferenceException, no OutOfMemoryException or quadratic blow-up on a huge
/// header, and no silent emission of a corrupt record. Every input must resolve to
/// EITHER a well-defined, theory-correct parse, OR a DOCUMENTED, intentional
/// validation exception. A raw, unhandled runtime exception, a hang, or a silently
/// corrupt record on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing": parsers of file formats are
///   THE prime fuzz target because they sit on the untrusted-input boundary.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PARSE-FASTA-001 — FASTA parsing (FileIO)
/// Checklist: docs/checklists/03_FUZZING.md, row 64. This is the first FileIO-area
/// fuzz unit and the single highest-value fuzz target in the checklist.
/// Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — binary garbage / random characters fed as FASTA text.
///   • TF  = Truncated Fields — a record with no trailing newline; a header with no
///           sequence; a defline-less sequence block (the leading record header is
///           "truncated" away).
///   • MC  = Malformed Content — missing '>' marker, empty sequence body.
///   • INJ = Injection — embedded NUL bytes and other control / unicode characters
///           inside the defline and the sequence body.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The FASTA-parsing contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// FASTA is a line-oriented format: each record is a defline beginning with '>'
/// (its first whitespace-delimited token is the identifier, the remainder is the
/// optional description) followed by one or more sequence lines.[Lipman &amp; Pearson
/// 1985; NCBI BLAST FASTA guidance]
///   — docs/algorithms/FileIO/FASTA_Parsing.md §1, §2.1, §2.2.
///
/// API entry points (FastaParser.cs):
///   • IEnumerable&lt;FastaEntry&gt; FastaParser.Parse(string fastaContent)
///   • IEnumerable&lt;FastaEntry&gt; FastaParser.ParseFile(string filePath)
/// Both delegate to the SAME line-oriented state machine `ParseReader`, so the
/// string surface fully exercises the parsing contract; one file-path test below
/// additionally pins the `ParseFile` StreamReader path (a common real-world entry).
///
/// CRITICAL ENUMERATION NOTE: `Parse` is a `yield`-based iterator — NOTHING runs,
/// and no exception can surface, until the caller enumerates. Every test therefore
/// materializes the result with `.ToList()` so that any work (and any documented
/// throw) actually happens inside the assertion.
///
/// Documented, repository-specific behavior the tests pin (FASTA_Parsing.md §3.3,
/// §5.2, §5.3, §6.1; FastaParser.cs; DnaSequence.cs):
///   • DEFLINE REQUIRED (INV-01): an entry is emitted only after a '>' defline has
///     been seen AND at least one sequence character collected. A sequence block
///     with NO leading '>' therefore produces ZERO entries — the parser drops it
///     silently rather than crashing (header stays null). The fuzz test pins the
///     empty result, NOT a throw. (FASTA_Parsing.md §6.1 row 1, §2.2.)
///   • HEADER-ONLY DROPPED: a defline with no following sequence line is NOT
///     yielded, because emission is gated on `sequenceBuilder.Length > 0`
///     (FastaParser.cs lines 117 / 135). So "empty sequence" → record skipped, not
///     an empty-sequence record and not a crash. (FASTA_Parsing.md §5.2, §6.1.)
///   • DNA-ONLY MATERIALIZATION: each entry's payload is built with
///     `new DnaSequence(sequence)`, which upper-cases then REJECTS any character
///     outside A/C/G/T with an `ArgumentException` (DnaSequence.cs lines 112–124).
///     Therefore binary garbage / NUL bytes / non-DNA symbols that land in a
///     SEQUENCE body surface as a DOCUMENTED `ArgumentException` during enumeration
///     — never an unhandled IndexOutOfRange / NullReference / KeyNotFound. This is
///     the documented "rejected with ArgumentException" contract of §3.3, and the
///     fuzz tests assert exactly that exception type, never a weaker "NotThrow".
///   • WHITESPACE STRIPPED: inside sequence lines every `char.IsWhiteSpace` char is
///     dropped before materialization, so blank lines contribute nothing and a
///     record split across many lines concatenates cleanly (INV-02).
///   • HEADERS ARE NOT ALPHABET-VALIDATED: NUL bytes / control / unicode chars in a
///     DEFLINE are preserved verbatim (only the leading '>' is stripped and the
///     remainder trimmed and split on the first space/tab). A 1 MB defline is just a
///     long string token — O(n) Substring + Trim, no quadratic blow-up, no OOM.
///   • NO-TRAILING-NEWLINE: `TextReader.ReadLine` returns the final unterminated
///     line, so the last record is still parsed — the key real-world robustness case.
///   • NULL / EMPTY / WHITESPACE-ONLY input → no entries (`yield break` on
///     `string.IsNullOrWhiteSpace`, FastaParser.cs line 19).
///
/// Determinism note: random-byte payloads use a LOCAL fixed-seed `new Random(seed)`
/// (no shared static Rng), so every run is byte-for-byte reproducible. The
/// huge-header and binary-garbage tests carry `[CancelAfter]` to convert any
/// pathological hang into a deterministic failure rather than a stuck run.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-FASTQ-001 — FASTQ parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 65. The second FileIO-area fuzz
/// unit. Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — fixed-seed random characters fed into the sequence
///           AND quality fields (and a whole random-byte blob with no '@').
///   • TF  = Truncated Fields — a record that ends after only 1 / 2 / 3 of the 4
///           canonical lines (file truncated mid-record). THE key boundary case:
///           the classic IndexOutOfRange-on-a-truncated-4-line-record trap.
///   • MC  = Malformed Content — header line missing its leading '@'; separator
///           line missing its leading '+'; quality length ≠ sequence length.
///   • INJ = Injection — NUL bytes / control chars inside the quality string.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The FASTQ-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// A canonical FASTQ record is FOUR logical lines: an '@'-header, a sequence, a
/// '+' separator, and a quality string whose length equals the sequence length;
/// quality is Phred-encoded printable ASCII (33..126).[Cock et al. 2009; NCBI SRA]
///   — docs/algorithms/FileIO/FASTQ_Parsing.md §2.1, §2.2, INV-01, INV-02.
///
/// API entry points (FastqParser.cs):
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.Parse(string content, QualityEncoding)
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.ParseFile(string path, QualityEncoding)
///   • IEnumerable&lt;FastqRecord&gt; FastqParser.Parse(TextReader, QualityEncoding)
/// All three delegate to the SAME `Parse(TextReader,...)` line state machine, so the
/// string surface fully exercises the contract; one file-path test additionally pins
/// the `ParseFile` StreamReader path. Like FASTA, `Parse` is a `yield`-based iterator,
/// so every test materializes with `.ToList()` to force the work to actually run.
///
/// CRITICAL CONTRACT DIVERGENCE FROM FASTA — TOLERANT, NOT STRICT.
/// FASTQ_Parsing.md §1, §5.2, §5.3, §6.2 state IN WRITING that this parser favors
/// "tolerant record assembly over strict format validation" and that "malformed
/// records can be skipped or partially assembled rather than rejected with an
/// error" — there is explicitly NO strict FASTQ validator. The DnaSequence alphabet
/// gate that made FASTA throw does NOT apply here: a `FastqRecord` stores the
/// sequence as a raw `string` (FastqParser.cs line 19–24), never materializing a
/// DnaSequence, so non-DNA / garbage sequence content does NOT throw. THEREFORE the
/// theory-correct contract these fuzz tests pin is the documented one:
///   • The BUSINESS GUARANTEE — the parser must NEVER crash, hang, throw an
///     undocumented runtime exception (IndexOutOfRange on a truncated 4-line
///     record!, NullReference, OOM), or loop forever on malformed/truncated/random/
///     injected input. EVERY fuzz test below asserts `NotThrow` + a pinned, exact
///     structural outcome, so the documented tolerance can never silently drift into
///     a crash, and a partially-assembled record's shape is nailed down precisely.
/// Documented behaviors pinned per target (FASTQ_Parsing.md §3.3, §4.1, §5.2):
///   • MISSING '@' HEADER (MC): a line not starting with '@' is SKIPPED
///     (FastqParser.cs line 96–97). A whole "record" whose header lacks '@' is
///     dropped — the parser does not crash and does not invent a record.
///   • MISSING '+' SEPARATOR (MC): the sequence loop accumulates every line until
///     one starts with '+' (line 104). With NO '+', the sequence absorbs the lines
///     that should have been the separator AND the quality, and the quality string
///     comes out EMPTY — a partially-assembled, NOT-rejected record (the documented
///     §5.2 tolerance). We pin that exact shape; the point is it does not crash.
///   • QUALITY-LENGTH ≠ SEQUENCE-LENGTH (MC): the quality loop reads lines only
///     until `qualityBuilder.Length >= sequence.Length` (line 115). A SHORT quality
///     run at EOF yields a record whose QualityString is shorter than the Sequence
///     (qual-len &lt; seq-len, NOT rejected); a quality line LONGER than the sequence
///     is appended whole (Trim never truncates), yielding qual-len &gt; seq-len. The
///     documented contract is tolerant assembly, so we pin BOTH mismatched shapes as
///     no-throw rather than asserting a rejection the parser explicitly does not do.
///   • INVALID QUALITY CHARS (INJ): `DecodeQualityScores` is `Math.Max(0, c-offset)`
///     for every char (line 178–181) — any byte, incl. NUL / control / chars outside
///     33..126, decodes (clamped at 0) without validation or crash. Pinned: a record
///     with control-char quality parses, QualityString preserved verbatim, scores
///     all ≥ 0 and never throwing.
///   • TRUNCATED RECORD (TF): every read loop guards `ReadLine() != null`, and
///     `DecodeQualityScores("")` returns an empty array (line 172–173), so a record
///     truncated to 1 / 2 / 3 lines yields a record with empty/short sequence and/or
///     empty quality and NEVER an IndexOutOfRange. This is the highest-value TF case.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`, line 73).
/// Determinism note: the random-byte FASTQ tests use a LOCAL fixed-seed
/// `new Random(seed)` and carry `[CancelAfter]`, exactly as the FASTA tests do.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-BED-001 — BED parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 66. The third FileIO-area fuzz
/// unit. Fuzz strategies exercised for THIS unit:
///   • TF  = Truncated Fields — a line with FEWER than the 3 required columns; a
///           record cut off after only chrom (or chrom+start). THE classic
///           IndexOutOfRange-on-`fields[1]`/`fields[2]` trap on a short split.
///   • MC  = Malformed Content — non-numeric chromStart/chromEnd ("abc"); an
///           interval with start > end (a negative-length feature).
///   • BE  = Boundary Exploitation — negative coordinates (-1), `int.MaxValue`
///           and an overflowing coordinate (a value larger than `int.MaxValue`).
///   • INJ = Injection — extra interior tabs producing empty fields; a tab
///           embedded inside what should be a single field; NUL / control bytes
///           in the chrom name; a fixed-seed random-byte blob.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes; BE = boundary
///   values 0, -1, MaxInt; INJ = injection).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The BED-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// BED is a TAB-delimited interval format whose required core is `chrom`,
/// `chromStart` (0-based, inclusive) and `chromEnd` (0-based, EXCLUSIVE), with
/// optional name/score/strand/… up to BED12.[UCSC BED FAQ; Wikipedia BED]
///   — docs/algorithms/FileIO/BED_Parsing.md §1, §2.2, §3.3, §6.1.
///
/// API entry points (BedParser.cs):
///   • IEnumerable&lt;BedRecord&gt; BedParser.Parse(string content, BedFormat)
///   • IEnumerable&lt;BedRecord&gt; BedParser.ParseFile(string path, BedFormat)
///   • IEnumerable&lt;BedRecord&gt; BedParser.Parse(TextReader, BedFormat)
/// All three delegate to the SAME `Parse(TextReader,…)` line state machine, so the
/// string surface fully exercises the contract; one file-path test additionally
/// pins the `ParseFile` StreamReader path. Like FASTA/FASTQ, `Parse` is a
/// `yield`-based iterator — every test materializes with `.ToList()` so any work
/// (and any throw) actually runs inside the assertion.
///
/// CRITICAL CONTRACT — SKIP-THE-BAD-LINE, NEVER THROW, NEVER CORRUPT.
/// The BED parser is a SKIPPING validator: a malformed data line is dropped and
/// parsing continues; it does NOT throw a parse exception and does NOT emit a
/// corrupt record. The business guarantee these fuzz tests pin is therefore that
/// EVERY malformed/truncated/injected line resolves to EITHER a well-defined,
/// theory-correct `BedRecord` OR a clean skip — NEVER an unhandled
/// FormatException / IndexOutOfRange / OverflowException / hang. Each fuzz test
/// asserts `NotThrow` PLUS a pinned, exact structural outcome so the documented
/// behavior can never silently drift into a crash or a corrupt interval.
/// Documented behaviors pinned per target (BED_Parsing.md §3.3, §6.1;
/// BedParser.cs ParseLine):
///   • FEWER THAN 3 COLUMNS (TF): `ParseLine` splits on TAB, falls back to
///     whitespace split, and returns `null` when still `&lt; 3` fields
///     (BedParser.cs lines 161–169). In `Auto` mode `Parse` also skips lines whose
///     field count is `&lt; 3` BEFORE `ParseLine` (lines 137–139). So `fields[1]` /
///     `fields[2]` are NEVER indexed on a short line → the classic IndexOutOfRange
///     trap CANNOT fire. Pinned: a 1- or 2-column line yields ZERO records.
///   • NON-NUMERIC chromStart/chromEnd (MC): coordinates are read with
///     `int.TryParse` (lines 173–176); a non-numeric token makes TryParse return
///     `false` → the line returns `null` (skipped). This is the documented
///     "Invalid numeric coordinates → Line skipped" of §6.1 — NOT an unhandled
///     `FormatException` from an `int.Parse`. Pinned: "abc" coords → ZERO records.
///   • start &gt; end (MC): an explicit `if (chromStart &gt; chromEnd) return null`
///     guard (lines 179–180) rejects a negative-length interval. This is the KEY
///     correctness check — a robust parser must catch start&gt;end, not emit a
///     feature with `Length &lt; 0`. Pinned: a start&gt;end line yields ZERO records,
///     and `chromStart == chromEnd` (a zero-length insertion point) is ACCEPTED
///     per §6.1 / INV-02.
///   • NEGATIVE COORDINATES (BE): `int.TryParse(NumberStyles.Integer)` ACCEPTS a
///     leading '-', and the parser has NO explicit non-negative guard. Per the
///     repository contract (BED_Parsing.md §3.3 + INV-02) the ONLY coordinate
///     rejection rules are "not an integer" and "start &gt; end" — negativity by
///     itself is NOT a documented rejection. So a line with a negative `chromStart`
///     whose `start &lt;= end` PARSES, producing a record with a negative
///     `ChromStart`. We pin the DOCUMENTED behavior exactly (parse, negative start
///     preserved), and separately pin that a negative-start line whose start &gt; end
///     is still rejected by the start&gt;end guard. (We assert the repo's real
///     contract, not an idealized "reject all negatives" the parser does not do.)
///   • int.MaxValue / OVERFLOW (BE): `chromStart = int.MaxValue`, `chromEnd =
///     int.MaxValue` parses (a zero-length feature at the boundary, start==end). A
///     coordinate LARGER than `int.MaxValue` overflows `int.TryParse` → TryParse
///     returns `false` → line skipped, NOT an unhandled `OverflowException`.
///   • TAB INJECTION (INJ): the parser splits the line on '\t' (line 161). Extra
///     interior tabs create EMPTY string fields but never break the column index;
///     a tab embedded where a single field was expected simply shifts the column
///     count. Pinned: injected tabs are handled deterministically (parsed with the
///     shifted columns, or skipped on the consistency/field-count rule) with NO
///     IndexOutOfRange on the split.
///   • RANDOM BYTES (INJ): a fixed-seed random-byte blob must be consumed
///     deterministically — no crash, no hang. Carries `[CancelAfter]`; the LOCAL
///     fixed seed makes the run byte-for-byte reproducible.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`,
///     BedParser.cs lines 104–105); `track `/`browser `/`#` lines are skipped.
/// Determinism note: the random-byte BED test uses a LOCAL fixed-seed
/// `new Random(seed)` and carries `[CancelAfter]`, exactly as the FASTA/FASTQ
/// tests do.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-VCF-001 — VCF parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 67. The fourth FileIO-area fuzz
/// unit. Fuzz strategies exercised for THIS unit:
///   • TF  = Truncated Fields — a data line with FEWER than the 8 mandatory columns
///           (truncated mid-record). THE classic IndexOutOfRange-on-`fields[1]`…
///           `fields[7]` trap on a short TAB split — the highest-value VCF fuzz case.
///   • MC  = Malformed Content — a missing `#CHROM` column-header line; a missing
///           `##fileformat` meta line; a non-numeric / non-integer `POS`; invalid
///           genotypes in the GT sample field ("9/9" referencing a non-existent
///           allele, malformed "1|", non-numeric ".|x").
///   • INJ = Injection — a 1 MB REF and a 1 MB ALT allele (huge alleles) that must be
///           consumed in linear time with no OOM / quadratic blow-up.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The VCF-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// VCF is a TAB-delimited text format: `##` meta lines (one of which is
/// `##fileformat=…`), a `#CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT samples…]`
/// column-header line, then one TAB-separated data line per variant locus. The first
/// 8 columns are mandatory; POS is a 1-based integer (INV-01, INV-02).[SAMtools
/// hts-specs VCF; Wikipedia VCF]
///   — docs/algorithms/FileIO/VCF_Parsing.md §1, §2.1, §2.2, §2.4.
///
/// API entry points (VcfParser.cs):
///   • IEnumerable&lt;VcfRecord&gt; VcfParser.Parse(string content)
///   • IEnumerable&lt;VcfRecord&gt; VcfParser.ParseFile(string filePath)
///   • (VcfHeader, IEnumerable&lt;VcfRecord&gt;) VcfParser.ParseWithHeader(string content)
/// `Parse`/`ParseFile` are `yield`-based iterators, so every test materializes with
/// `.ToList()` to force the parse (and any throw) to actually run inside the assertion;
/// one file-path test additionally pins the `ParseFile` File.ReadAllText path.
///
/// CRITICAL CONTRACT — TOLERANT PARSE, SKIP THE BAD LINE, NEVER CRASH OR CORRUPT.
/// The VCF parser defers most semantic validation to callers (VCF_Parsing.md §1,
/// §5.3, §6.2): it does the minimum structural validation (≥ 8 columns, integer POS)
/// and stores everything else as strings. The business guarantee these fuzz tests pin
/// is therefore that EVERY malformed / truncated / injected input resolves to EITHER a
/// well-defined, theory-correct `VcfRecord` OR a clean skip — NEVER an unhandled
/// IndexOutOfRange (on a short data line!), FormatException (on a non-integer POS),
/// OutOfMemory (on a huge allele), or hang. Each fuzz test asserts `NotThrow` PLUS a
/// pinned, exact structural outcome so the documented behavior can never silently
/// drift into a crash or a corrupt variant record.
/// Documented behaviors pinned per target (VCF_Parsing.md §3.3, §5.2, §6.1;
/// VcfParser.cs Parse / ParseLine):
///   • FEWER THAN 8 COLUMNS (TF): `ParseLine` splits on TAB and returns `null` when
///     `fields.Length &lt; 8` (VcfParser.cs lines 301–303) BEFORE indexing `fields[1]`…
///     `fields[7]`, so the classic IndexOutOfRange trap CANNOT fire. Pinned: a 4- or
///     7-column data line yields ZERO records (the documented "fewer than 8 columns →
///     skipped" of §6.1). The KEY TF boundary case.
///   • NON-INTEGER `POS` (MC): POS is read with `int.TryParse` (line 306–307); a
///     non-numeric / float token makes TryParse return `false` → the line returns
///     `null` (skipped). This is the documented "Non-integer POS → skipped" of §6.1,
///     NOT an unhandled `FormatException` from an `int.Parse`. Pinned: ZERO records.
///   • MISSING `#CHROM` HEADER (MC): `Parse` does NOT require a `#CHROM` line — it
///     skips `##`/`#` lines and parses every non-`#` line as data (lines 130–146). With
///     no `#CHROM` line `sampleNames` stays null, so FORMAT/sample columns are NOT
///     parsed (the `fields.Length &gt; 9 &amp;&amp; sampleNames != null` gate, line 324), but the
///     8 mandatory columns STILL parse. The DOCUMENTED behavior is therefore: data
///     lines parse, samples are silently dropped — NOT a rejection and NOT a crash. We
///     pin that exact shape (records present, `Samples == null`) rather than an
///     idealized "reject when #CHROM missing" the parser does not implement.
///   • MISSING `##fileformat` (MC): in `Parse(string)` ALL `##` lines are skipped
///     uniformly — `##fileformat` carries NO special status, so its absence has ZERO
///     effect on record parsing. (`ParseWithHeader` separately DEFAULTS `FileFormat`
///     to `"VCFv4.3"` when no `##fileformat=` line is seen — VcfParser.cs line 176.)
///     The repo contract is "parse anyway / default the header", NOT "reject". We pin
///     BOTH: `Parse` ignores the missing meta line and still parses the records, and
///     `ParseWithHeader` yields the default `VCFv4.3` FileFormat.
///   • INVALID GENOTYPES (MC): the GT sample field is stored as a RAW string in the
///     sample dictionary (ParseLine lines 331–337) with NO allele-range validation at
///     parse time, so a GT like "9/9" (allele index with no matching ALT), a malformed
///     "1|" (trailing separator), or a non-numeric ".|x" is parsed and STORED verbatim
///     — never a crash. The downstream zygosity helpers (`IsHet`/`IsHomAlt`/
///     `IsHomRef`, lines 529–563) then interpret them DETERMINISTICALLY: they split on
///     `/`|`|`, and a "." allele or an out-of-shape genotype simply yields false rather
///     than throwing. We pin that the invalid GTs parse, survive verbatim in the sample
///     dict, and that every zygosity helper returns a defined bool without crashing.
///   • HUGE ALLELES (INJ/OVF): a 1 MB REF and a 1 MB ALT string. REF is stored as a
///     plain string and ALT is `Split(',')` once (line 311) — both O(n) over the allele
///     length, no quadratic re-scan and no unbounded allocation. The record parses with
///     the full-length alleles intact. `[CancelAfter]` converts any pathological stall
///     into a deterministic failure. The point: huge alleles complete without OOM.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`, line 118–119).
/// Determinism note: every VCF fuzz input below is a FIXED literal (no randomness);
/// the huge-allele test carries `[CancelAfter]`, exactly as the FASTA/FASTQ/BED tests.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-GFF-001 — GFF parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 68. The fifth FileIO-area fuzz
/// unit, targeting `Seqeron.Genomics.IO.GffParser` — the FULL file parser, DISTINCT
/// from the Annotation-area `GenomeAnnotator.ParseGff3` (ANNOT-GFF-001, row 31).
/// Fuzz strategies exercised for THIS unit:
///   • TF  = Truncated Fields — a data line with FEWER than the 8 mandatory columns
///           (a row truncated to seqid/source/type/start only). THE classic
///           IndexOutOfRange-on-`fields[3]`…`fields[7]` trap on a short TAB split.
///   • MC  = Malformed Content — non-integer / float `start`/`end`; malformed
///           attributes (a `key` with no `=`, a trailing `;`, an empty `key=`,
///           duplicate keys).
///   • BE  = Boundary Exploitation — a NEGATIVE coordinate (`-5`) and `start > end`
///           (a reversed interval).
///   • INJ = Injection — an INVALID strand char (not one of `+`/`-`/`.`/`?`);
///           PERCENT-ENCODED reserved characters in attribute values (`%3D`=`=`,
///           `%2C`=`,`, `%3B`=`;`, `%09`=tab) — the KEY INJ target: a parser that
///           splits column 9 on `;`/`=` MUST handle encoded separators correctly.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GFF-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// GFF3/GTF is a TAB-delimited NINE-column annotation format —
/// `seqid source type start end score strand phase attributes` — with 1-based,
/// fully-closed coordinates; strand is one of `+`/`-`/`.`/`?`; GFF3 reserved
/// characters in attribute values are PERCENT-encoded (`%3D`, `%2C`, `%09`, …).
/// [Sequence Ontology GFF3 v1.26; UCSC; GTF2.2]
///   — docs/algorithms/FileIO/GFF_Parsing.md §2.1, §2.2, §2.4 (INV-01…INV-04).
///
/// API entry points (GffParser.cs):
///   • IEnumerable&lt;GffRecord&gt; GffParser.Parse(string content, GffFormat)
///   • IEnumerable&lt;GffRecord&gt; GffParser.ParseFile(string path, GffFormat)
///   • IEnumerable&lt;GffRecord&gt; GffParser.Parse(TextReader, GffFormat)
/// All three delegate to the SAME `Parse(TextReader,…)` line state machine, so the
/// string surface fully exercises the contract; one file-path test additionally pins
/// the `ParseFile` StreamReader path. `Parse` is a `yield`-based iterator — every test
/// materializes with `.ToList()` so any work (and any throw) actually runs inside the
/// assertion.
///
/// CRITICAL CONTRACT — TOLERANT/SKIP, NEVER CRASH, NEVER CORRUPT.
/// `GffParser` does the minimum STRUCTURAL validation (≥ 8 columns, integer
/// start/end) and stores everything else as-is (GFF_Parsing.md §3.3, §6.1). The
/// business guarantee these fuzz tests pin is that EVERY malformed / truncated /
/// injected line resolves to EITHER a well-defined, theory-correct `GffRecord` OR a
/// clean skip — NEVER an unhandled IndexOutOfRange (on a short line!),
/// FormatException (on a non-integer start/end), or hang. Each fuzz test asserts
/// `NotThrow` PLUS a pinned, exact structural outcome so the documented behavior can
/// never silently drift into a crash or a corrupt feature. (NOTE: ANNOT-GFF-001 found
/// a real FormatException/IndexOutOfRange bug in the OTHER GFF parser. This FileIO
/// parser is hardened against that class of bug — `fields.Length &lt; 8` guards the
/// short-line trap, and `int.TryParse` guards the non-integer-coordinate trap — and
/// the fuzz tests below PIN that hardening so it cannot regress.)
/// Documented behaviors pinned per target (GFF_Parsing.md §3.3, §6.1;
/// GffParser.cs ParseLine / ParseAttributes / UnescapeGff):
///   • FEWER THAN 8 COLUMNS (TF): `ParseLine` splits on TAB and returns `null` when
///     `fields.Length &lt; 8` (GffParser.cs lines 121–123) BEFORE indexing `fields[3]`…
///     `fields[7]`, so the classic IndexOutOfRange trap CANNOT fire. Pinned: a 4-column
///     row yields ZERO records (the documented "fewer than 8 columns → skipped" of
///     §6.1). An 8-column row (no attributes) IS accepted with an empty attribute dict.
///     The KEY TF boundary case.
///   • NON-INTEGER start/end (MC): start/end are read with `int.TryParse`
///     (lines 129–132); a non-numeric / float token makes TryParse return `false` →
///     the line returns `null` (skipped). This is the documented "start/end not an
///     integer → skipped" of §3.3 — NOT an unhandled `FormatException` from a
///     `int.Parse`. Pinned: a `"abc"` / `"10.5"` start yields ZERO records.
///   • NEGATIVE COORDINATE (BE): `int.TryParse(NumberStyles.Integer)` ACCEPTS a
///     leading `-`, and `ParseLine` has NO non-negative guard AND NO `start &gt; end`
///     guard (unlike the BED parser). Per the repo contract (§3.3 lists ONLY the
///     "&lt; 8 columns" and "non-integer start/end" rejections), a negative coordinate
///     PARSES, producing a record with a negative `Start`/`End`. We pin the DOCUMENTED
///     behavior exactly (parse, negative preserved), not an idealized rejection the
///     parser does not implement.
///   • start &gt; end (BE): there is NO `start &gt; end` rejection in `ParseLine`, so a
///     reversed interval PARSES (Start &gt; End) rather than being skipped — the
///     INV-02 `start &lt;= end` is a property of VALID input, not a parse-time guard.
///     We pin the documented behavior (the reversed record is emitted, no crash),
///     distinct from BED where start&gt;end is explicitly rejected.
///   • INVALID STRAND (INJ): `strand = fields[6].Length &gt; 0 ? fields[6][0] : '.'`
///     (line 141) — the parser takes the FIRST char of column 7 with NO membership
///     check against `+`/`-`/`.`/`?`. So an invalid strand like `"X"` is PRESERVED
///     verbatim as `'X'` (not rejected, not normalized), and an EMPTY column 7 becomes
///     `'.'`. We pin BOTH: an invalid strand parses with the char preserved, and an
///     empty strand defaults to `'.'` — the documented "store as-is" behavior, never a
///     crash.
///   • MALFORMED ATTRIBUTES (MC): column 9 is split on `;` with
///     `RemoveEmptyEntries`, then each part split on the FIRST `=` with the
///     `eqIdx &gt; 0` guard (lines 183–193). A part with NO `=` is SILENTLY DROPPED
///     (eqIdx == -1); a TRAILING `;` produces no empty entry (RemoveEmptyEntries); an
///     empty `=value` (eqIdx == 0) is dropped (the key would be empty); a DUPLICATE
///     key is LAST-WINS (`attributes[key] = value`). We pin each of these exact shapes
///     — the parser never crashes on malformed column 9 and resolves it
///     deterministically.
///   • PERCENT-ENCODING (INJ — THE KEY TARGET): seqid/source/type AND every GFF3
///     attribute key/value pass through `UnescapeGff` = `Uri.UnescapeDataString`
///     (lines 125–127, 189–190, 199–202). CRUCIALLY the split on `;`/`=` happens on
///     the RAW (still-encoded) text BEFORE unescaping, so an ENCODED separator
///     (`%3D` = `=`, `%3B` = `;`, `%2C` = `,`, `%09` = tab) inside an attribute VALUE
///     survives the split and is then correctly decoded into the value — a parser that
///     decoded first would mis-split. We pin: `Name=val%3Due` → value `"val=ue"` (the
///     encoded `=` does NOT create a spurious split), `%2C` → `,`, `%09` → tab. A
///     MALFORMED percent-sequence (`%ZZ`, a lone trailing `%`) is passed through
///     VERBATIM by `Uri.UnescapeDataString` (verified: it does NOT throw), so a broken
///     escape can never crash the parser. This is the contrast with the BED/VCF units:
///     GFF is the one format with mandatory percent-decoding, and the INJ target is
///     precisely that the `;`/`=` splitter must not be fooled by encoded separators.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`, line 71–72);
///     `##`/`#` directive and comment lines are skipped.
/// Determinism note: every GFF fuzz input below is a FIXED literal (no randomness),
/// so each run is byte-for-byte reproducible; no `[CancelAfter]` is needed because no
/// input is large or random — the parser is a single linear pass over short literals.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: PARSE-GENBANK-001 — GenBank parsing (FileIO)
/// ═══════════════════════════════════════════════════════════════════════════
/// Checklist: docs/checklists/03_FUZZING.md, row 69. The sixth FileIO-area fuzz
/// unit, targeting `Seqeron.Genomics.IO.GenBankParser` — the GenBank flat-file
/// parser. Fuzz strategies exercised for THIS unit:
///   • RB  = Random Bytes — a fixed-seed random-byte blob (no `LOCUS` keyword) fed
///           as GenBank text.
///   • TF  = Truncated Fields — a record TRUNCATED mid-FEATURES (a feature line cut
///           off with no location; a qualifier line with no `=value`); and a record
///           with NO `//` terminator (the file ends mid-record). THE key loop/hang
///           boundary case.
///   • MC  = Malformed Content — a record with NO `LOCUS` line (the keyword that
///           gates parsing); an INVALID feature location (`abc..xyz`, an unbalanced
///           `complement(` with no closing paren).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GenBank-parsing contract under test (the DOCUMENTED, repository-specific one)
/// ───────────────────────────────────────────────────────────────────────────
/// A GenBank flat file is a set of keyworded sections starting in column 1 —
/// `LOCUS`, `DEFINITION`, `ACCESSION`, `VERSION`, `KEYWORDS`, `SOURCE`/`ORGANISM`,
/// `REFERENCE`, a `FEATURES` table (feature key + INSDC location + `/qualifier`
/// lines), an `ORIGIN` sequence block — and each record is TERMINATED by `//`.
/// [NCBI GenBank Sample Record; NCBI GenBank Overview; INSDC Feature Table]
///   — docs/algorithms/FileIO/GenBank_Parsing.md §2.1, §2.2, INV-01.
///
/// API entry points (GenBankParser.cs):
///   • IEnumerable&lt;GenBankRecord&gt; GenBankParser.Parse(string content)
///   • IEnumerable&lt;GenBankRecord&gt; GenBankParser.ParseFile(string filePath)
/// Both `Parse`/`ParseFile` are `yield`-based iterators, so every test materializes
/// with `.ToList()` to force the parse (and any throw) to actually run inside the
/// assertion; one file-path test additionally pins the `ParseFile` File.ReadAllText
/// path.
///
/// CRITICAL CONTRACT — TERMINATOR HANDLING IS A SPLIT, NOT A `while not //` LOOP.
/// The single highest-value fuzz concern for a GenBank parser is an unbounded
/// `while (line != "//")` read loop that hangs when the terminator is missing. THIS
/// parser has NO such loop: `Parse` does `content.Split(new[]{"\n//"}, ...)` over the
/// WHOLE text in one O(n) pass (GenBankParser.cs line 97), so a record with NO `//`
/// simply becomes ONE block that is still parsed — it can never loop forever. The
/// no-terminator test below carries `[CancelAfter]` as a loop/hang tripwire and pins
/// that the partial record IS emitted (not dropped, not a hang). The business
/// guarantee these fuzz tests pin is that EVERY malformed / truncated / random input
/// resolves to EITHER a well-defined, theory-correct `GenBankRecord` OR a clean skip
/// — NEVER an unhandled IndexOutOfRange (on a truncated feature/qualifier!),
/// NullReference (on a record missing sections), or hang (on a missing `//`).
/// Documented behaviors pinned per target (GenBank_Parsing.md §3.3, §5.2, §6.1;
/// GenBankParser.cs Parse / ParseRecord / ParseFeatures / ParseLocation):
///   • MISSING `LOCUS` (MC): `Parse` only calls `ParseRecord` for a trimmed block
///     that `StartsWith("LOCUS")` (line 102). A block with NO `LOCUS` line is SKIPPED
///     — the documented "parsed only when the block begins with LOCUS" of §3.3 / §5.2.
///     Pinned: a record body without a `LOCUS` keyword yields ZERO records, and a
///     two-record input where only the second has `LOCUS` yields exactly ONE record.
///   • TRUNCATED FEATURES (TF): `ParseFeatures` (lines 438–513) is defensive — it
///     guards every line index (`line.Length > 5`, `line.Length > 21`), so a feature
///     line CUT OFF with no location yields `currentLocation == ""` (a zeroed
///     `Location` via `ParseLocation("")`), and a qualifier line with NO `=` is stored
///     as `qualifiers[name] = "true"` (line 492). No `fields[n]` / Substring index can
///     ever go out of bounds → the classic IndexOutOfRange-on-a-truncated-record trap
///     CANNOT fire. Pinned: a record truncated mid-FEATURES parses, the partial
///     feature is emitted with a zeroed location, and no crash.
///   • INVALID LOCATION (MC): `ParseLocation` delegates to
///     `SequenceFormatHelper.ParseLocationParts`, which is REGEX-driven — it finds
///     `\d+(?:\.\.\d+)?` ranges (SequenceFormatHelper.cs line 57) and ignores
///     everything else. So an invalid location like `abc..xyz` matches NO ranges →
///     a zeroed `Location` (Start==0, End==0, no Parts) with the RawLocation preserved
///     verbatim; an unbalanced `complement(` with no closing paren still only scans for
///     digit ranges → no IndexOutOfRange, no mismatched-paren crash. `IsComplement` is
///     a `StartsWith("complement(")` flag, so the broken `complement(` is flagged but
///     not dereferenced. Pinned: invalid locations parse to a zeroed location with the
///     raw string retained, never an IndexOutOfRange.
///   • NO `//` TERMINATOR (TF — THE KEY LOOP/HANG CASE): a record whose text ends
///     mid-ORIGIN with no `//`. Because parsing is a `Split` not a loop, the whole
///     text is one block beginning with `LOCUS`, so the partial record IS parsed and
///     emitted (locus + whatever sections were present). Pinned: the unterminated
///     record is emitted with the correct locus, NOT dropped and NOT a hang. Carries
///     `[CancelAfter]` as the loop tripwire.
///   • RANDOM BYTES (RB): a fixed-seed random-byte blob that does not begin with
///     `LOCUS` is a single block that fails the `StartsWith("LOCUS")` gate → ZERO
///     records, no crash, no hang. Carries `[CancelAfter]`; the LOCAL fixed seed makes
///     the run byte-for-byte reproducible.
///   • NULL / EMPTY input → no records (`yield break` on `IsNullOrEmpty`, lines 93–94).
/// Determinism note: the random-byte GenBank test uses a LOCAL fixed-seed
/// `new Random(seed)` and the random-byte + no-terminator tests carry `[CancelAfter]`,
/// exactly as the FASTA/FASTQ/BED tests do.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
// I/O-bound fixture: the file-path tests share a per-fixture `_tempDir` set up in
// [SetUp] and removed in [TearDown]. Under the assembly-wide
// `Parallelizable(ParallelScope.Children)` policy a single fixture instance is shared
// across parallel test methods, so a concurrent [TearDown] can delete the temp dir
// another test is mid-write into (DirectoryNotFoundException). Per the documented
// house rule in Parallelization.cs ("I/O-bound fixtures should be marked
// [NonParallelizable]") this fixture runs its tests serially.
[NonParallelizable]
public class FileIoFuzzTests
{
    #region Helpers — temp-dir management for the ParseFile (file-path) surface

    private string _tempDir = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(),
            "SeqeronFastaFuzz_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_tempDir);
    }

    [TearDown]
    public void TearDown()
    {
        try
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, recursive: true);
        }
        catch (IOException)
        {
            // Best-effort cleanup; a leaked temp file must not fail the suite.
        }
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp file and returns its path, so the file-path tests control byte
    /// layout exactly (in particular: whether a trailing newline exists).
    /// </summary>
    private string WriteTempFasta(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".fasta");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.fastq` file and returns its path, so the FASTQ file-path tests
    /// control byte layout exactly (in particular: truncation mid-record).
    /// </summary>
    private string WriteTempFastq(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".fastq");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.bed` file and returns its path, so the BED file-path tests control
    /// byte layout exactly (in particular: TAB layout and truncation/injection).
    /// </summary>
    private string WriteTempBed(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".bed");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.vcf` file and returns its path, so the VCF file-path tests control
    /// byte layout exactly (in particular: a truncated final data line / missing header).
    /// </summary>
    private string WriteTempVcf(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".vcf");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.gff3` file and returns its path, so the GFF file-path tests control
    /// byte layout exactly (in particular: a truncated final data line / TAB layout).
    /// </summary>
    private string WriteTempGff(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".gff3");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    /// <summary>
    /// Writes <paramref name="content"/> verbatim (no added trailing newline) to a
    /// fresh temp `.gb` file and returns its path, so the GenBank file-path tests
    /// control byte layout exactly (in particular: a record with NO `//` terminator,
    /// or a record truncated mid-FEATURES).
    /// </summary>
    private string WriteTempGenBank(string content)
    {
        string path = Path.Combine(_tempDir, "in_" + Guid.NewGuid().ToString("N") + ".gb");
        File.WriteAllText(path, content, new UTF8Encoding(encoderShouldEmitUTF8Identifier: false));
        return path;
    }

    #endregion

    #region PARSE-FASTA-001 — FASTA parsing

    #region Positive sanity — a valid multi-record FASTA parses to the correct ids/sequences

    /// <summary>
    /// Positive control: a well-formed two-record FASTA (the worked example from
    /// FASTA_Parsing.md §7.1) must parse to exactly two entries with the correct
    /// id / description / sequence. If this fails, the fuzz suite is meaningless —
    /// it proves the parser is actually wired up and the happy path is intact.
    /// </summary>
    [Test]
    public void Parse_ValidMultiRecord_ParsesToCorrectEntries()
    {
        const string fasta =
            ">seq1 First sequence\n" +
            "ATGCATGC\n" +
            ">seq2 Second sequence\n" +
            "GGCCTTAA\n";

        List<FastaEntry> entries = FastaParser.Parse(fasta).ToList();

        entries.Should().HaveCount(2, "both well-formed records must be emitted");

        entries[0].Id.Should().Be("seq1");
        entries[0].Description.Should().Be("First sequence");
        entries[0].Sequence.Sequence.Should().Be("ATGCATGC");

        entries[1].Id.Should().Be("seq2");
        entries[1].Description.Should().Be("Second sequence");
        entries[1].Sequence.Sequence.Should().Be("GGCCTTAA");
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a multi-line, multi-record FASTA
    /// written to disk (sequence wrapped across lines, exactly as real files do)
    /// parses via `ParseFile` to the same correct entries. This pins the StreamReader
    /// path and the multi-line concatenation invariant (INV-02).
    /// </summary>
    [Test]
    public void ParseFile_ValidMultiLineMultiRecord_ParsesToCorrectEntries()
    {
        const string fasta =
            ">chrA alpha contig\n" +
            "ATGC\n" +
            "ATGC\n" +
            ">chrB\n" +
            "GGGGCCCC\n";
        string path = WriteTempFasta(fasta);

        List<FastaEntry> entries = FastaParser.ParseFile(path).ToList();

        entries.Should().HaveCount(2);
        entries[0].Id.Should().Be("chrA");
        entries[0].Description.Should().Be("alpha contig");
        entries[0].Sequence.Sequence.Should().Be("ATGCATGC",
            "wrapped sequence lines must concatenate (INV-02)");
        entries[1].Id.Should().Be("chrB");
        entries[1].Description.Should().BeNull("a description-less defline yields a null Description");
        entries[1].Sequence.Sequence.Should().Be("GGGGCCCC");
    }

    #endregion

    #region MC — Missing '>': a defline-less sequence block yields NO entries (no crash)

    /// <summary>
    /// MC / TF: input that is a bare sequence block with NO leading '>' defline.
    /// Per INV-01 (FASTA_Parsing.md §2.2, §6.1 row 1) an entry is emitted only after
    /// a defline has been seen, so `header` never leaves null and the buffered
    /// characters are dropped. The documented, theory-correct outcome is ZERO
    /// entries — a clean parse, NOT an exception and certainly NOT a crash. We pin
    /// the empty result exactly so this silent-drop contract cannot drift.
    /// </summary>
    [Test]
    public void Parse_MissingHeaderMarker_YieldsNoEntries()
    {
        const string headerless = "ATGCATGCATGC\nGGCCTTAAGGCC\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(headerless).ToList();

        act.Should().NotThrow("a defline-less block is dropped, not a crash");
        entries.Should().BeEmpty(
            "no entry can be emitted without a preceding '>' defline (INV-01)");
    }

    #endregion

    #region TF — Empty sequence: a header with no sequence line is dropped (no crash)

    /// <summary>
    /// TF / MC: a defline immediately followed by another defline (and a final
    /// defline with nothing after it) — i.e. headers with EMPTY sequence bodies.
    /// Emission is gated on `sequenceBuilder.Length > 0` (FastaParser.cs line 117 /
    /// 135), so header-only records are DROPPED (FASTA_Parsing.md §5.2, §6.1). The
    /// only record with sequence content survives; no empty-sequence record is
    /// emitted and nothing crashes. We verify the surviving record AND the dropping.
    /// </summary>
    [Test]
    public void Parse_HeaderWithEmptySequence_DropsEmptyRecords()
    {
        const string fasta =
            ">empty1\n" +     // header, then immediately another header → empty seq
            ">withSeq data\n" +
            "ACGT\n" +
            ">empty2\n";      // trailing header, no sequence at all

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("empty-sequence records are dropped, not a crash");
        entries.Should().ContainSingle("only the one record with sequence content survives")
            .Which.Id.Should().Be("withSeq");
        entries[0].Sequence.Sequence.Should().Be("ACGT");
    }

    #endregion

    #region RB — Binary garbage / random bytes: documented ArgumentException, never a crash

    /// <summary>
    /// RB: a valid defline followed by a sequence body of fixed-seed RANDOM bytes
    /// (the full 0x00–0xFF range as chars). Because the body is materialized via
    /// `new DnaSequence(...)`, any non-A/C/G/T character is REJECTED with a
    /// DOCUMENTED `ArgumentException` (DnaSequence.cs lines 112–124; FASTA_Parsing.md
    /// §3.3). That is the contract: random garbage is rejected deterministically, it
    /// is NEVER an unhandled IndexOutOfRange / NullReference / OutOfMemory and never a
    /// silently-corrupt record. We assert the exact exception type. `[CancelAfter]`
    /// guards against any pathological hang while scanning the garbage.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Parse_BinaryGarbageSequenceBody_ThrowsDocumentedArgumentException(
        CancellationToken token)
    {
        var rng = new Random(0xFA57A);             // local fixed seed — fully reproducible
        var sb = new StringBuilder(">garbage record\n");
        for (int i = 0; i < 4096; i++)
            sb.Append((char)rng.Next(0, 256));     // arbitrary bytes incl. control chars
        sb.Append('\n');
        string fasta = sb.ToString();

        // Enumerate eagerly so the DnaSequence materialization (and its throw) runs.
        var act = () => FastaParser.Parse(fasta).ToList();

        act.Should().Throw<ArgumentException>(
            "non-A/C/G/T sequence content is rejected by DnaSequence, the documented contract")
           .And.Message.Should().Contain("nucleotide");
        token.IsCancellationRequested.Should().BeFalse("parsing garbage must not hang");
    }

    /// <summary>
    /// RB / robustness companion: pure random bytes that do NOT begin with '>' must
    /// be handled deterministically. The first random byte is almost never '>' (and
    /// we force it to be a non-'>' byte), so with no defline the whole block is a
    /// defline-less sequence buffer that is never emitted (INV-01) → ZERO entries,
    /// no crash, no hang. This pins the "random bytes that aren't even FASTA" case.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Parse_PureRandomBytesNoDefline_YieldsNoEntriesAndDoesNotCrash(
        CancellationToken token)
    {
        var rng = new Random(0xB10C);
        var sb = new StringBuilder();
        for (int i = 0; i < 8192; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '>') c = 'x';   // guarantee no defline is ever introduced
            if (c == '\n' || c == '\r') c = 'y'; // keep it a single buffered "line"
            sb.Append(c);
        }

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("defline-less random bytes are buffered and dropped, not a crash");
        entries.Should().BeEmpty("with no '>' line no entry can ever be emitted (INV-01)");
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    #endregion

    #region INJ — Embedded NUL bytes: handled deterministically

    /// <summary>
    /// INJ: NUL bytes (\0) injected into BOTH the defline and the sequence body.
    /// • In the DEFLINE, '\0' is NOT alphabet-validated — it is preserved verbatim
    ///   as part of the header token (only '>' is stripped, the rest trimmed/split).
    ///   So the parser must not crash on a NUL in the header.
    /// • In the SEQUENCE body, '\0' is not whitespace, so it reaches DnaSequence and
    ///   is REJECTED with the documented `ArgumentException` (it is not A/C/G/T).
    /// We pin BOTH halves: a NUL-bearing sequence throws the documented exception,
    /// while a NUL purely in the header (with a clean sequence) parses cleanly and
    /// the NUL survives in the id. Deterministic either way — never an unhandled crash.
    /// </summary>
    [Test]
    public void Parse_NullBytesInSequence_ThrowsDocumentedArgumentException()
    {
        const string fastaWithNulInSeq = ">rec\nAC\0GT\n";

        var act = () => FastaParser.Parse(fastaWithNulInSeq).ToList();

        act.Should().Throw<ArgumentException>(
            "a NUL byte in the sequence body is a non-DNA char rejected by DnaSequence");
    }

    /// <summary>
    /// INJ companion: a NUL byte injected only into the DEFLINE (header), with a
    /// clean A/C/G/T sequence, must parse deterministically — the header is not
    /// alphabet-validated, so the record is emitted and the NUL survives verbatim in
    /// the parsed id token. No crash, no silent loss of the record.
    /// </summary>
    [Test]
    public void Parse_NullByteInHeaderOnly_ParsesAndPreservesNul()
    {
        const string fasta = ">id\0X clean\nACGT\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("a NUL in the header is preserved, not validated, so it must not crash");
        entries.Should().ContainSingle();
        entries[0].Id.Should().Be("id\0X", "the first whitespace-delimited token, NUL included");
        entries[0].Sequence.Sequence.Should().Be("ACGT");
    }

    #endregion

    #region OVF — Huge header (1 MB defline): completes in linear time, no OOM/quadratic blow-up

    /// <summary>
    /// OVF / TF: a single defline whose header is 1 MB long (a header-only record,
    /// no sequence). The parser strips the leading '>', trims, and splits on the
    /// first space/tab — all O(n) over the header length, no quadratic re-scan, no
    /// unbounded allocation. Because the record has no sequence body it is dropped
    /// (INV / header-only rule) → ZERO entries. The point of THIS test is that the
    /// huge header is processed without OutOfMemory or a quadratic stall, guarded by
    /// `[CancelAfter]`. We add a sequence after the giant header to ALSO prove the
    /// parser recovers and still emits the following real record.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Parse_HugeHeader_CompletesWithoutBlowUpAndRecovers(CancellationToken token)
    {
        const int headerLen = 1_000_000;                 // 1 MB defline
        string hugeHeader = new string('H', headerLen);
        string fasta = ">" + hugeHeader + "\n" +         // header-only → dropped
                       ">realId real desc\n" +
                       "ACGTACGT\n";

        List<FastaEntry> entries = new();
        var act = () => entries = FastaParser.Parse(fasta).ToList();

        act.Should().NotThrow("a 1 MB header must be processed in linear time, no OOM/quadratic blow-up");
        entries.Should().ContainSingle("the header-only giant record is dropped; the real record survives")
            .Which.Id.Should().Be("realId");
        entries[0].Sequence.Sequence.Should().Be("ACGTACGT");
        token.IsCancellationRequested.Should().BeFalse("a huge header must not stall the parser");
    }

    #endregion

    #region TF — No trailing newline: the last record is still parsed (key real-world case)

    /// <summary>
    /// TF: a FASTA whose final record has NO trailing newline — extremely common in
    /// real files and a classic off-by-one trap. `TextReader.ReadLine` returns the
    /// last unterminated line, and the end-of-input flush emits the final buffered
    /// record (FastaParser.cs lines 135–138). So the last record MUST still parse.
    /// We assert it over BOTH surfaces: the in-memory `Parse(string)` and the
    /// `ParseFile` StreamReader path (with a temp file written WITHOUT a trailing
    /// newline), since the file path is where this case bites in practice.
    /// </summary>
    [Test]
    public void Parse_NoTrailingNewline_StillParsesLastRecord()
    {
        const string noNewline = ">a first\nACGT\n>b second\nTTTT"; // no final '\n'

        List<FastaEntry> entries = FastaParser.Parse(noNewline).ToList();

        entries.Should().HaveCount(2, "the unterminated final record must still be emitted");
        entries[1].Id.Should().Be("b");
        entries[1].Description.Should().Be("second");
        entries[1].Sequence.Sequence.Should().Be("TTTT");
    }

    /// <summary>
    /// TF over the FILE-PATH surface: the same no-trailing-newline case written to
    /// disk and read back via `ParseFile`, pinning that the StreamReader flush emits
    /// the final unterminated record exactly as the string surface does.
    /// </summary>
    [Test]
    public void ParseFile_NoTrailingNewline_StillParsesLastRecord()
    {
        const string noNewline = ">x\nAAAA\n>y last\nCCCC"; // file ends mid-record
        string path = WriteTempFasta(noNewline);

        List<FastaEntry> entries = FastaParser.ParseFile(path).ToList();

        entries.Should().HaveCount(2);
        entries[1].Id.Should().Be("y");
        entries[1].Description.Should().Be("last");
        entries[1].Sequence.Sequence.Should().Be("CCCC",
            "ParseFile must flush the final unterminated record");
    }

    #endregion

    #endregion

    #region PARSE-FASTQ-001 — FASTQ parsing

    #region Positive sanity — a valid multi-record FASTQ parses to correct id/seq/quality

    /// <summary>
    /// Positive control: a well-formed two-record FASTQ must parse to exactly two
    /// records with the correct id / description / sequence / quality string, and the
    /// per-record quality-string length must equal the sequence length (INV-01). If
    /// this fails the whole FASTQ fuzz suite is meaningless — it proves the happy path
    /// is wired up before we throw garbage at it.
    /// </summary>
    [Test]
    public void ParseFastq_ValidMultiRecord_ParsesToCorrectRecords()
    {
        const string fastq =
            "@SEQ1 first read\n" +
            "ACGTACGT\n" +
            "+\n" +
            "IIIIIIII\n" +
            "@SEQ2\n" +
            "GGCCTTAA\n" +
            "+\n" +
            "!!!!!!!!\n";

        List<FastqParser.FastqRecord> records = FastqParser.Parse(fastq).ToList();

        records.Should().HaveCount(2, "both well-formed records must be emitted");

        records[0].Id.Should().Be("SEQ1");
        records[0].Description.Should().Be("first read");
        records[0].Sequence.Should().Be("ACGTACGT");
        records[0].QualityString.Should().Be("IIIIIIII");
        records[0].QualityScores.Should().HaveCount(records[0].Sequence.Length,
            "valid FASTQ has one quality score per base (INV-01)");
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0));

        records[1].Id.Should().Be("SEQ2");
        records[1].Description.Should().Be("", "a description-less header yields an empty Description");
        records[1].Sequence.Should().Be("GGCCTTAA");
        records[1].QualityString.Should().Be("!!!!!!!!");
        records[1].QualityScores.Should().HaveCount(8);
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a valid two-record FASTQ written
    /// to disk parses via `ParseFile` to the same correct records, pinning the
    /// StreamReader path (the common real-world entry point).
    /// </summary>
    [Test]
    public void ParseFileFastq_ValidMultiRecord_ParsesToCorrectRecords()
    {
        const string fastq =
            "@readA alpha\n" +
            "ACGT\n" +
            "+\n" +
            "IIII\n" +
            "@readB\n" +
            "TTTT\n" +
            "+\n" +
            "####\n";
        string path = WriteTempFastq(fastq);

        List<FastqParser.FastqRecord> records = FastqParser.ParseFile(path).ToList();

        records.Should().HaveCount(2);
        records[0].Id.Should().Be("readA");
        records[0].Description.Should().Be("alpha");
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Should().Be("IIII");
        records[1].Id.Should().Be("readB");
        records[1].Sequence.Should().Be("TTTT");
        records[1].QualityString.Should().Be("####");
    }

    #endregion

    #region MC — Missing '@' header: the malformed block is skipped, never a crash

    /// <summary>
    /// MC: the first "record" has a header line that does NOT begin with '@' (the '@'
    /// is missing). Per FastqParser.cs line 96–97 every line that does not start with
    /// '@' is SKIPPED, so the entire malformed block (header-less + its seq/+/qual
    /// lines, none of which start with '@') is dropped and NO record is invented for
    /// it. The following well-formed record is still parsed. The documented tolerant
    /// contract is "skip, do not crash" — we pin the surviving record exactly and
    /// assert no throw, so the skip cannot silently become a crash.
    /// </summary>
    [Test]
    public void ParseFastq_MissingAtHeaderMarker_SkipsMalformedBlockNoCrash()
    {
        const string fastq =
            "SEQ_no_at\n" +   // header missing '@' → skipped
            "ACGT\n" +        // these lines also lack '@' → skipped
            "+\n" +
            "IIII\n" +
            "@good\n" +       // a clean record after the garbage
            "GGCC\n" +
            "+\n" +
            "JJJJ\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a header without '@' is skipped, not a crash (tolerant parser)");
        records.Should().ContainSingle("only the well-formed '@good' record can be emitted")
            .Which.Id.Should().Be("good");
        records[0].Sequence.Should().Be("GGCC");
        records[0].QualityString.Should().Be("JJJJ");
    }

    #endregion

    #region MC — Missing '+' separator: sequence absorbs the rest, quality empty, no crash

    /// <summary>
    /// MC: a record with NO '+' separator line. The sequence-accumulation loop
    /// (FastqParser.cs line 104) reads lines until one begins with '+'; with no such
    /// line it absorbs BOTH the intended sequence AND the intended quality line, then
    /// hits EOF, leaving the quality string EMPTY. This is the documented §5.2
    /// "partially assembled rather than rejected" behavior — a tolerant parse, NOT a
    /// crash and NOT an exception. We pin the exact partially-assembled shape so the
    /// documented tolerance is nailed down and can never drift into an
    /// IndexOutOfRange/NullReference.
    /// </summary>
    [Test]
    public void ParseFastq_MissingPlusSeparator_AbsorbsIntoSequenceNoCrash()
    {
        const string fastq =
            "@rec\n" +
            "ACGT\n" +   // intended sequence
            "IIII\n";    // intended quality — but with no '+' it is absorbed too

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a missing '+' separator is tolerated (partial assembly), not a crash");
        records.Should().ContainSingle("the record is partially assembled, not rejected");
        records[0].Id.Should().Be("rec");
        records[0].Sequence.Should().Be("ACGTIIII",
            "with no '+' the sequence loop absorbs the would-be quality line too (§5.2)");
        records[0].QualityString.Should().BeEmpty(
            "nothing is left for the quality loop after the sequence loop hits EOF");
        records[0].QualityScores.Should().BeEmpty("an empty quality string decodes to no scores");
    }

    #endregion

    #region MC — Quality length ≠ sequence length: tolerant partial assembly, never a crash

    /// <summary>
    /// MC (qual-len &lt; seq-len): the quality run is SHORTER than the sequence and the
    /// input ends before the quality loop fills up. The loop (line 115) stops at EOF,
    /// so the record carries a QualityString shorter than its Sequence — the parser
    /// does NOT reject the corruption (documented tolerant behavior, §5.2/§6.2), it
    /// emits the mismatched record. The business guarantee under fuzzing is that this
    /// must not crash and must not produce more scores than quality chars. We pin the
    /// exact short-quality shape so the documented tolerance is explicit and stable.
    /// </summary>
    [Test]
    public void ParseFastq_QualityShorterThanSequence_TolerantNoCrash()
    {
        const string fastq =
            "@rec\n" +
            "ACGTACGT\n" +  // 8-base sequence
            "+\n" +
            "II\n";          // only 2 quality chars, then EOF

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("a too-short quality run is tolerated at EOF, not a crash");
        records.Should().ContainSingle();
        records[0].Sequence.Should().Be("ACGTACGT");
        records[0].QualityString.Should().Be("II",
            "the quality loop stops at EOF with quality shorter than the sequence (no rejection)");
        records[0].QualityString.Length.Should().BeLessThan(records[0].Sequence.Length,
            "this is the documented qual-len < seq-len corruption the tolerant parser accepts");
        records[0].QualityScores.Should().HaveCount(records[0].QualityString.Length,
            "scores are decoded per quality char, never more than the quality length");
    }

    /// <summary>
    /// MC (qual-len &gt; seq-len): the quality LINE is LONGER than the sequence. The
    /// quality loop appends the whole trimmed line in one shot (Trim never truncates),
    /// so `qualityBuilder.Length` overshoots the sequence length and the QualityString
    /// comes out LONGER than the Sequence. Again tolerant assembly, not rejection — we
    /// pin the over-long shape and assert no crash. The following clean record proves
    /// the parser recovers after the over-long line is fully consumed.
    /// </summary>
    [Test]
    public void ParseFastq_QualityLongerThanSequence_TolerantNoCrash()
    {
        const string fastq =
            "@rec1\n" +
            "ACGT\n" +          // 4-base sequence
            "+\n" +
            "IIIIIIIIIIII\n" +  // 12 quality chars (over-long, appended whole)
            "@rec2\n" +
            "TTTT\n" +
            "+\n" +
            "JJJJ\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("an over-long quality line is tolerated, not a crash");
        records.Should().HaveCount(2, "the parser recovers and still emits the following record");
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Length.Should().BeGreaterThan(records[0].Sequence.Length,
            "the whole over-long quality line is appended (Trim does not truncate)");
        records[1].Id.Should().Be("rec2", "the parser recovers after the over-long quality line");
        records[1].QualityString.Should().Be("JJJJ");
    }

    #endregion

    #region INJ — Invalid / control quality chars: decoded without validation, never a crash

    /// <summary>
    /// INJ: quality characters OUTSIDE the canonical Phred printable range 33..126 —
    /// here NUL (\0) and other control bytes injected into the quality string (with a
    /// matching-length sequence). `DecodeQualityScores` is `Math.Max(0, c - offset)`
    /// for every char (line 178–181): it never validates the alphabet, so control /
    /// out-of-range chars decode (clamped at ≥ 0) without throwing. The documented
    /// contract is "decode, do not reject"; the fuzz guarantee is no crash. We pin
    /// that the record parses, the QualityString is preserved verbatim (NUL included),
    /// and every decoded score is ≥ 0.
    /// </summary>
    [Test]
    public void ParseFastq_InvalidControlCharQuality_DecodesWithoutCrash()
    {
        // 4-base sequence; 4 quality chars, all control bytes below ASCII 33.
        const string fastq = "@rec\nACGT\n+\n\0\n";

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow("out-of-range/control quality chars are decoded, not validated (no crash)");
        records.Should().ContainSingle();
        records[0].Sequence.Should().Be("ACGT");
        records[0].QualityString.Should().Be("\0",
            "control-char quality is preserved verbatim, not sanitized");
        records[0].QualityScores.Should().HaveCount(4);
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0),
            "DecodeQualityScores clamps every value to ≥ 0, even for sub-offset bytes");
    }

    #endregion

    #region TF — Truncated record (1/2/3 of 4 lines): never IndexOutOfRange, defined result

    /// <summary>
    /// TF — THE key boundary case: a file truncated after only 1, 2, or 3 of the four
    /// canonical FASTQ lines. The classic trap is an IndexOutOfRange when code assumes
    /// four lines are always present. Every read loop here guards `ReadLine() != null`
    /// and `DecodeQualityScores("")` returns an empty array (line 172–173), so each
    /// truncation yields a well-DEFINED record (empty/short sequence and/or empty
    /// quality) and NEVER an unhandled exception. We assert no-throw plus the exact
    /// shape for each truncation point.
    /// </summary>
    [Test]
    public void ParseFastq_TruncatedRecords_NeverIndexOutOfRange()
    {
        // 1 line only: header, then EOF.
        const string oneLine = "@only_header\n";
        // 2 lines: header + sequence, then EOF (no '+', no quality).
        const string twoLines = "@hdr\nACGT\n";
        // 3 lines: header + sequence + '+' separator, then EOF (no quality).
        const string threeLines = "@hdr\nACGT\n+\n";

        List<FastqParser.FastqRecord> r1 = new(), r2 = new(), r3 = new();
        var act1 = () => r1 = FastqParser.Parse(oneLine).ToList();
        var act2 = () => r2 = FastqParser.Parse(twoLines).ToList();
        var act3 = () => r3 = FastqParser.Parse(threeLines).ToList();

        act1.Should().NotThrow("a header-only truncation must not IndexOutOfRange");
        act2.Should().NotThrow("a header+sequence truncation must not IndexOutOfRange");
        act3.Should().NotThrow("a header+sequence+'+' truncation must not IndexOutOfRange");

        // 1 line: header consumed; sequence loop hits EOF → empty seq; quality empty.
        r1.Should().ContainSingle();
        r1[0].Id.Should().Be("only_header");
        r1[0].Sequence.Should().BeEmpty();
        r1[0].QualityString.Should().BeEmpty();
        r1[0].QualityScores.Should().BeEmpty();

        // 2 lines: no '+', so the sequence absorbs "ACGT"; quality empty at EOF.
        r2.Should().ContainSingle();
        r2[0].Sequence.Should().Be("ACGT");
        r2[0].QualityString.Should().BeEmpty();

        // 3 lines: '+' stops the sequence loop at "ACGT"; quality loop hits EOF → empty.
        r3.Should().ContainSingle();
        r3[0].Sequence.Should().Be("ACGT");
        r3[0].QualityString.Should().BeEmpty(
            "the quality loop hits EOF immediately, leaving an empty (not crashed) quality");
        r3[0].QualityScores.Should().BeEmpty();
    }

    /// <summary>
    /// TF over the FILE-PATH surface: a multi-record FASTQ whose FINAL record is
    /// truncated mid-record (only 3 of the 4 lines, no trailing newline) is written to
    /// disk and read via `ParseFile`. The StreamReader flush must emit the complete
    /// first record AND the truncated final record without an IndexOutOfRange — the
    /// real-world shape where this bug actually bites.
    /// </summary>
    [Test]
    public void ParseFileFastq_FinalRecordTruncated_StillParsesNoCrash()
    {
        const string fastq =
            "@full\n" +
            "ACGT\n" +
            "+\n" +
            "IIII\n" +
            "@truncated\n" +
            "GGCC\n" +
            "+";              // file ends here: no quality line, no trailing newline
        string path = WriteTempFastq(fastq);

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.ParseFile(path).ToList();

        act.Should().NotThrow("a truncated final record on disk must not crash ParseFile");
        records.Should().HaveCount(2, "both the full and the truncated records are emitted");
        records[0].Id.Should().Be("full");
        records[0].QualityString.Should().Be("IIII");
        records[1].Id.Should().Be("truncated");
        records[1].Sequence.Should().Be("GGCC");
        records[1].QualityString.Should().BeEmpty("the truncated record has no quality at EOF");
    }

    #endregion

    #region RB — Random bytes: handled deterministically, no crash, no hang

    /// <summary>
    /// RB: a valid '@' header followed by fixed-seed RANDOM bytes (full 0x00–0xFF
    /// range) in the SEQUENCE/quality body. Unlike FASTA, the FASTQ sequence is stored
    /// as a raw string and never alphabet-validated, so random garbage does NOT throw —
    /// the documented tolerant contract. The fuzz guarantee is that the parser consumes
    /// the garbage deterministically: no crash, no hang, and any decoded quality scores
    /// are ≥ 0. `[CancelAfter]` converts any pathological hang into a deterministic
    /// failure; the LOCAL fixed seed makes the run byte-for-byte reproducible.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseFastq_RandomByteBody_DoesNotCrashOrHang(CancellationToken token)
    {
        var rng = new Random(0xFA5701);            // local fixed seed — fully reproducible
        var sb = new StringBuilder("@garbage record\n");
        for (int i = 0; i < 4096; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '\n' || c == '\r') c = 'N';   // keep it inside one logical body line
            sb.Append(c);
        }
        sb.Append('\n');
        string fastq = sb.ToString();

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(fastq).ToList();

        act.Should().NotThrow(
            "random bytes in a FASTQ body are tolerated (raw string, no alphabet gate), not a crash");
        records.Should().ContainSingle("one '@' header → one (tolerantly assembled) record");
        records[0].QualityScores.Should().AllSatisfy(s => s.Should().BeGreaterThanOrEqualTo(0),
            "every decoded score is clamped to ≥ 0 regardless of the byte");
        token.IsCancellationRequested.Should().BeFalse("parsing garbage must not hang");
    }

    /// <summary>
    /// RB companion: a whole blob of fixed-seed random bytes with NO '@' header line.
    /// Every line fails the `StartsWith('@')` gate (line 96–97) and is skipped, so the
    /// parser yields ZERO records — no crash, no hang. We force out '@', '\n' and '\r'
    /// so the blob is genuinely one header-less garbage line. `[CancelAfter]` guards
    /// against any stall while scanning the garbage.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseFastq_PureRandomBytesNoHeader_YieldsNoRecordsNoCrash(CancellationToken token)
    {
        var rng = new Random(0xB10C);
        var sb = new StringBuilder();
        for (int i = 0; i < 8192; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '@') c = 'x';                  // guarantee no header is ever introduced
            if (c == '\n' || c == '\r') c = 'y';    // keep it a single header-less line
            sb.Append(c);
        }

        List<FastqParser.FastqRecord> records = new();
        var act = () => records = FastqParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("header-less random bytes are skipped line by line, not a crash");
        records.Should().BeEmpty("with no '@' line no record can ever be emitted");
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    #endregion

    #endregion

    #region PARSE-BED-001 — BED parsing

    #region Positive sanity — a valid BED parses to the correct chrom/start/end intervals

    /// <summary>
    /// Positive control: a well-formed multi-record BED6 (TAB-separated) must parse to
    /// exactly the right chrom / start / end (and name/score/strand) records, with the
    /// derived `Length = chromEnd - chromStart` (INV-01) and `chromStart &lt;= chromEnd`
    /// (INV-02). If this fails the BED fuzz suite is meaningless — it proves the parser
    /// is wired up and the happy path is intact before we throw garbage at it.
    /// — BED_Parsing.md §2.2, §6.1.
    /// </summary>
    [Test]
    public void ParseBed_ValidMultiRecord_ParsesToCorrectIntervals()
    {
        const string bed =
            "chr1\t100\t200\tfeatA\t500\t+\n" +
            "chr2\t0\t50\tfeatB\t1000\t-\n";

        List<BedParser.BedRecord> records = BedParser.Parse(bed).ToList();

        records.Should().HaveCount(2, "both well-formed BED lines must be emitted");

        records[0].Chrom.Should().Be("chr1");
        records[0].ChromStart.Should().Be(100);
        records[0].ChromEnd.Should().Be(200);
        records[0].Length.Should().Be(100, "Length = chromEnd - chromStart (INV-01)");
        records[0].Name.Should().Be("featA");
        records[0].Score.Should().Be(500);
        records[0].Strand.Should().Be('+');

        records[1].Chrom.Should().Be("chr2");
        records[1].ChromStart.Should().Be(0, "the first base in a chromosome is numbered 0 (0-based)");
        records[1].ChromEnd.Should().Be(50);
        records[1].Strand.Should().Be('-');
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a valid BED3 written to disk parses
    /// via `ParseFile` to the same correct intervals, pinning the StreamReader path (the
    /// common real-world entry point).
    /// </summary>
    [Test]
    public void ParseFileBed_ValidBed3_ParsesToCorrectIntervals()
    {
        const string bed =
            "chrX\t10\t20\n" +
            "chrX\t30\t45\n";
        string path = WriteTempBed(bed);

        List<BedParser.BedRecord> records = BedParser.ParseFile(path).ToList();

        records.Should().HaveCount(2);
        records[0].Chrom.Should().Be("chrX");
        records[0].ChromStart.Should().Be(10);
        records[0].ChromEnd.Should().Be(20);
        records[1].ChromStart.Should().Be(30);
        records[1].ChromEnd.Should().Be(45);
        records[1].Length.Should().Be(15);
    }

    #endregion

    #region MC — start > end: the negative-length interval is rejected (skipped), never emitted

    /// <summary>
    /// MC — THE key correctness check: a line whose `chromStart &gt; chromEnd` describes a
    /// negative-length feature, which is invalid (INV-02). `ParseLine` has an explicit
    /// `if (chromStart &gt; chromEnd) return null` guard (BedParser.cs lines 179–180), so
    /// the bad line is SKIPPED — the parser must NOT emit a record with `Length &lt; 0`.
    /// The valid neighbours on either side prove the parser recovers and keeps parsing.
    /// We also pin that `chromStart == chromEnd` (a zero-length insertion point) is
    /// ACCEPTED per §6.1 / INV-02 — the boundary is `start &lt;= end`, not `start &lt; end`.
    /// </summary>
    [Test]
    public void ParseBed_StartGreaterThanEnd_RejectsNegativeLengthInterval()
    {
        const string bed =
            "chr1\t100\t200\n" +   // valid
            "chr1\t500\t300\n" +   // start > end → invalid, must be skipped
            "chr1\t700\t700\n" +   // start == end → zero-length insertion point, valid
            "chr1\t800\t900\n";    // valid

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow("a start>end line is skipped, not a crash");
        records.Should().HaveCount(3,
            "the start>end line is dropped; the three start<=end lines survive");
        records.Should().NotContain(r => r.ChromStart > r.ChromEnd,
            "no negative-length interval may ever be emitted (INV-02)");
        records.Should().NotContain(r => r.Length < 0,
            "a robust parser must never produce a feature with Length < 0");
        records.Should().Contain(r => r.ChromStart == 700 && r.ChromEnd == 700,
            "a zero-length (start==end) insertion point is valid and must survive");
    }

    #endregion

    #region MC — Non-numeric coordinates: TryParse fails → line skipped, never a FormatException

    /// <summary>
    /// MC: `chromStart` and/or `chromEnd` are non-numeric ("abc", "12x", a float). The
    /// parser reads them with `int.TryParse` (BedParser.cs lines 173–176), so a
    /// non-numeric token makes TryParse return `false` and the line returns `null`
    /// (the documented "Invalid numeric coordinates → Line skipped", §6.1). The KEY
    /// guarantee is that this is a clean skip, NEVER an UNHANDLED `FormatException`
    /// that a naive `int.Parse` would raise. We pin that every non-numeric-coord line
    /// is dropped and only the clean line survives.
    /// </summary>
    [Test]
    public void ParseBed_NonNumericCoordinates_SkippedNotFormatException()
    {
        const string bed =
            "chr1\tabc\t200\n" +    // non-numeric start
            "chr1\t100\txyz\n" +    // non-numeric end
            "chr1\t12.5\t200\n" +   // float (not an int) start
            "chr1\t100\t200\n";     // the one clean line

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow(
            "non-numeric coordinates are skipped via int.TryParse, never an unhandled FormatException");
        records.Should().ContainSingle("only the line with integer coordinates parses")
            .Which.ChromStart.Should().Be(100);
        records[0].ChromEnd.Should().Be(200);
    }

    #endregion

    #region BE — Negative coordinates: documented accept-when-start<=end, reject-when-start>end

    /// <summary>
    /// BE (boundary value -1): a line with a NEGATIVE `chromStart` whose `start &lt;= end`.
    /// `int.TryParse(NumberStyles.Integer)` accepts a leading '-', and the parser has NO
    /// non-negative guard; per the repository contract (BED_Parsing.md §3.3 + INV-02) the
    /// ONLY coordinate rejections are "not an integer" and "start &gt; end" — negativity by
    /// itself is NOT a documented rejection. So this line PARSES with a negative
    /// `ChromStart`. We pin the DOCUMENTED behavior exactly (parse, negative start
    /// preserved) rather than asserting an idealized "reject all negatives" rule the
    /// parser does not implement. A negative-start line whose start &gt; end is STILL
    /// rejected by the start&gt;end guard — pinned in the companion below.
    /// </summary>
    [Test]
    public void ParseBed_NegativeStartWithStartLeEnd_ParsesPerDocumentedContract()
    {
        const string bed = "chr1\t-5\t10\n";

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow("a negative coordinate is parsed, not a crash");
        records.Should().ContainSingle(
            "negativity alone is not a documented rejection; only non-integer or start>end reject");
        records[0].ChromStart.Should().Be(-5, "int.TryParse accepts the leading '-' and it is preserved");
        records[0].ChromEnd.Should().Be(10);
    }

    /// <summary>
    /// BE companion: a NEGATIVE start that is ALSO greater than end (start = -1, end = -9)
    /// must still be rejected by the explicit `start &gt; end` guard — proving negative
    /// coordinates do not bypass the negative-length check. The valid line on either side
    /// survives.
    /// </summary>
    [Test]
    public void ParseBed_NegativeStartGreaterThanEnd_StillRejected()
    {
        const string bed =
            "chr1\t10\t20\n" +     // valid
            "chr1\t-1\t-9\n" +     // start(-1) > end(-9) → invalid, skipped
            "chr1\t30\t40\n";      // valid

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow();
        records.Should().HaveCount(2, "the negative start>end line is rejected by the start>end guard");
        records.Should().NotContain(r => r.ChromStart > r.ChromEnd);
        records.Should().NotContain(r => r.Length < 0);
    }

    #endregion

    #region BE — int.MaxValue and overflow: boundary parses; overflow is skipped, not OverflowException

    /// <summary>
    /// BE (MaxInt boundary): `chromStart = chromEnd = int.MaxValue` is a zero-length
    /// feature at the integer boundary (start == end) and PARSES. A coordinate value
    /// LARGER than `int.MaxValue` overflows `int.TryParse`, which returns `false`
    /// (NOT an `OverflowException` like `int.Parse` would raise), so that line is
    /// cleanly SKIPPED. We pin both: the boundary line survives with the exact
    /// `int.MaxValue` coordinates, and the overflowing line is dropped without a crash.
    /// </summary>
    [Test]
    public void ParseBed_MaxIntAndOverflow_BoundaryParsesOverflowSkipped()
    {
        string bed =
            "chr1\t" + int.MaxValue + "\t" + int.MaxValue + "\n" +  // start==end at boundary → valid
            "chr1\t0\t99999999999999999999\n" +                     // > int.MaxValue → overflow → skip
            "chr1\t1\t2\n";                                          // valid neighbour

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow(
            "an overflowing coordinate is skipped via int.TryParse, never an unhandled OverflowException");
        records.Should().HaveCount(2, "the overflowing line is dropped; the two valid lines survive");
        records.Should().Contain(r => r.ChromStart == int.MaxValue && r.ChromEnd == int.MaxValue,
            "the int.MaxValue boundary (zero-length) feature parses exactly");
        records.Should().Contain(r => r.ChromStart == 1 && r.ChromEnd == 2);
    }

    #endregion

    #region TF — Fewer than 3 columns: no IndexOutOfRange on the split, line skipped

    /// <summary>
    /// TF — THE classic IndexOutOfRange trap: lines with FEWER than the 3 required
    /// columns (chrom only; chrom+start only). A naive parser that does `fields[1]` /
    /// `fields[2]` would throw IndexOutOfRange. Here `Parse` skips lines with `&lt; 3`
    /// fields in Auto mode (BedParser.cs lines 137–139) and `ParseLine` also returns
    /// `null` when `&lt; 3` (lines 168–169), so `fields[1]`/`fields[2]` are NEVER indexed
    /// on a short line. We pin that 1- and 2-column lines yield ZERO records and never
    /// crash; the valid 3-column line proves recovery.
    /// </summary>
    [Test]
    public void ParseBed_FewerThanThreeColumns_NeverIndexOutOfRange()
    {
        const string bed =
            "chr1\n" +          // 1 column → too few
            "chr1\t100\n" +     // 2 columns → too few
            "chr1\t100\t200\n"; // 3 columns → valid

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow("a short line must never IndexOutOfRange on fields[1]/fields[2]");
        records.Should().ContainSingle("only the 3-column line is a valid record")
            .Which.ChromStart.Should().Be(100);
        records[0].ChromEnd.Should().Be(200);
    }

    /// <summary>
    /// TF over the FILE-PATH surface: a truncated 2-column final line on disk (no
    /// trailing newline) read via `ParseFile` must not IndexOutOfRange — the real-world
    /// shape where a file truncated mid-line bites. The full first line still parses.
    /// </summary>
    [Test]
    public void ParseFileBed_TruncatedFinalLine_NoCrash()
    {
        const string bed =
            "chr1\t100\t200\n" +
            "chr2\t300";          // file ends mid-line: only 2 columns, no newline
        string path = WriteTempBed(bed);

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.ParseFile(path).ToList();

        act.Should().NotThrow("a truncated final line on disk must not crash ParseFile");
        records.Should().ContainSingle("only the complete first line is a valid record")
            .Which.Chrom.Should().Be("chr1");
    }

    #endregion

    #region INJ — Tab injection: extra/interior tabs handled deterministically, no IndexOutOfRange

    /// <summary>
    /// INJ: a line with EXTRA interior tabs between the required fields
    /// (`chr1\t\t100\t200` → an empty field appears between chrom and start). The split
    /// on '\t' produces an empty string field, shifting the coordinate columns so
    /// `int.TryParse("")` fails on `fields[1]` and the line is cleanly SKIPPED — NOT an
    /// IndexOutOfRange and NOT a corrupt record. The clean 3-column line FIRST sets the
    /// Auto-mode expected field count to 3, so it parses; the 4-field injected-tab line
    /// is then dropped by BOTH the field-count-consistency rule AND the empty-start
    /// guard. The point is that injected tabs are absorbed deterministically by the
    /// column split with no crash and no corrupt interval.
    /// </summary>
    [Test]
    public void ParseBed_ExtraInteriorTabs_HandledDeterministicallyNoCrash()
    {
        const string bed =
            "chr2\t100\t200\n" +    // clean line first → Auto sets expectedFieldCount = 3
            "chr1\t\t100\t200\n";   // injected empty field → 4 fields + start="" → skipped

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow("injected interior tabs must not IndexOutOfRange on the split");
        records.Should().ContainSingle("only the clean 3-column line parses")
            .Which.Chrom.Should().Be("chr2");
        records[0].ChromStart.Should().Be(100);
        records.Should().NotContain(r => r.Length < 0,
            "an injected-tab line must never yield a corrupt interval");
    }

    /// <summary>
    /// INJ: NUL / control bytes injected into the CHROM name (a valid coordinate pair).
    /// The chrom name is NOT alphabet-validated — it is taken verbatim as `fields[0]` —
    /// so the parser must not crash and must preserve the bytes in `Chrom`. We pin that
    /// the record parses with the NUL preserved and correct coordinates, deterministically.
    /// </summary>
    [Test]
    public void ParseBed_NullByteInChromName_ParsesAndPreservesByte()
    {
        const string bed = "ch\0r1\t100\t200\n";

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(bed).ToList();

        act.Should().NotThrow("a NUL in the chrom name is preserved, not validated, so it must not crash");
        records.Should().ContainSingle();
        records[0].Chrom.Should().Be("ch\0r1", "the chrom token is taken verbatim, NUL included");
        records[0].ChromStart.Should().Be(100);
        records[0].ChromEnd.Should().Be(200);
    }

    #endregion

    #region INJ — Random bytes: consumed deterministically, no crash, no hang

    /// <summary>
    /// INJ / robustness: a fixed-seed RANDOM-byte blob (full 0x00–0xFF range) fed as BED
    /// text. Whatever line structure the bytes happen to form, every line resolves to a
    /// clean skip (too few fields / non-numeric coords / start>end / inconsistent field
    /// count) or a well-formed record — the parser must consume the garbage
    /// deterministically with NO crash and NO hang, and must never emit a negative-length
    /// interval. `[CancelAfter]` converts any pathological hang into a deterministic
    /// failure; the LOCAL fixed seed makes the run byte-for-byte reproducible.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseBed_RandomByteBlob_DoesNotCrashOrHangOrCorrupt(CancellationToken token)
    {
        var rng = new Random(0xB3D);               // local fixed seed — fully reproducible
        var sb = new StringBuilder();
        for (int i = 0; i < 8192; i++)
            sb.Append((char)rng.Next(0, 256));     // arbitrary bytes incl. control + newlines

        List<BedParser.BedRecord> records = new();
        var act = () => records = BedParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("random bytes must be consumed deterministically, not a crash");
        records.Should().NotContain(r => r.ChromStart > r.ChromEnd,
            "no garbage line may ever produce a negative-length interval (INV-02)");
        records.Should().NotContain(r => r.Length < 0);
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    #endregion

    #endregion

    #region PARSE-VCF-001 — VCF parsing

    #region Positive sanity — a valid VCF parses to the correct variant records

    /// <summary>
    /// Positive control: a well-formed VCF (meta lines + `#CHROM` header with one
    /// sample + two data lines) must parse to exactly two variants with the correct
    /// CHROM / POS / ID / REF / ALT / QUAL / FILTER / INFO and the per-sample GT. If
    /// this fails the VCF fuzz suite is meaningless — it proves the parser is wired up
    /// and the happy path is intact before we throw malformed input at it.
    /// — VCF_Parsing.md §2.2, §3.2.
    /// </summary>
    [Test]
    public void ParseVcf_ValidMultiRecord_ParsesToCorrectVariants()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n" +
            "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30\tGT\t0/1\n" +
            "chr2\t250\t.\tAT\tA,ATG\t.\t.\t.\tGT\t1/2\n";

        List<VcfParser.VcfRecord> records = VcfParser.Parse(vcf).ToList();

        records.Should().HaveCount(2, "both well-formed data lines must be emitted");

        records[0].Chrom.Should().Be("chr1");
        records[0].Pos.Should().Be(100, "POS is the 1-based integer coordinate (INV-01)");
        records[0].Id.Should().Be("rs1");
        records[0].Ref.Should().Be("A");
        records[0].Alt.Should().Equal("G");
        records[0].Qual.Should().Be(50.0);
        records[0].Filter.Should().Equal("PASS");
        records[0].Info.Should().ContainKey("DP").WhoseValue.Should().Be("30");
        records[0].Samples.Should().NotBeNull();
        VcfParser.GetGenotype(records[0], 0).Should().Be("0/1");

        records[1].Chrom.Should().Be("chr2");
        records[1].Pos.Should().Be(250);
        records[1].Id.Should().Be(".");
        records[1].Ref.Should().Be("AT");
        records[1].Alt.Should().Equal("A", "ATG");
        records[1].Qual.Should().BeNull("QUAL '.' is normalized to null (§6.1)");
        records[1].Filter.Should().BeEmpty("FILTER '.' is normalized to an empty array (§6.1)");
        VcfParser.GetGenotype(records[1], 0).Should().Be("1/2");
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: the same valid VCF written to disk
    /// parses via `ParseFile` to the same correct variants, pinning the
    /// `File.ReadAllText` path (the common real-world entry point).
    /// </summary>
    [Test]
    public void ParseFileVcf_ValidMultiRecord_ParsesToCorrectVariants()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chrX\t1000\t.\tC\tT\t99.5\tPASS\t.\n" +
            "chrX\t2000\t.\tG\tA\t.\tq10\tAF=0.25\n";
        string path = WriteTempVcf(vcf);

        List<VcfParser.VcfRecord> records = VcfParser.ParseFile(path).ToList();

        records.Should().HaveCount(2);
        records[0].Chrom.Should().Be("chrX");
        records[0].Pos.Should().Be(1000);
        records[0].Ref.Should().Be("C");
        records[0].Alt.Should().Equal("T");
        records[0].Qual.Should().Be(99.5);
        records[1].Pos.Should().Be(2000);
        records[1].Filter.Should().Equal("q10");
        records[1].Info.Should().ContainKey("AF").WhoseValue.Should().Be("0.25");
    }

    #endregion

    #region TF — Truncated data line (< 8 columns): line skipped, never IndexOutOfRange

    /// <summary>
    /// TF — THE key boundary case: data lines truncated to FEWER than the 8 mandatory
    /// columns (4 columns; 7 columns). The classic trap is an IndexOutOfRange when code
    /// indexes `fields[1]`…`fields[7]` assuming 8 columns are always present. `ParseLine`
    /// returns `null` when `fields.Length &lt; 8` (VcfParser.cs lines 301–303) BEFORE any
    /// indexing, so the short line is cleanly SKIPPED — the documented "fewer than 8
    /// columns → skipped" of §6.1. We pin that the short lines yield ZERO records and
    /// never crash; the full 8-column line on either side proves the parser recovers.
    /// </summary>
    [Test]
    public void ParseVcf_TruncatedDataLine_SkippedNeverIndexOutOfRange()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\trs1\tA\n" +            // 4 columns → too few, skipped
            "chr1\t200\trs2\tA\tG\t.\t.\n" +   // 7 columns → too few, skipped
            "chr1\t300\trs3\tA\tG\t.\t.\t.\n"; // 8 columns → valid

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow(
            "a short data line must never IndexOutOfRange on fields[1]..fields[7]");
        records.Should().ContainSingle("only the complete 8-column line is a valid record")
            .Which.Pos.Should().Be(300);
        records[0].Chrom.Should().Be("chr1");
    }

    /// <summary>
    /// TF over the FILE-PATH surface: a multi-record VCF whose FINAL data line is
    /// truncated mid-record (only 3 columns, no trailing newline) is written to disk and
    /// read via `ParseFile`. The full first data line still parses and the truncated
    /// final line is cleanly skipped — the real-world shape where a file truncated
    /// mid-line bites — with NO IndexOutOfRange.
    /// </summary>
    [Test]
    public void ParseFileVcf_TruncatedFinalDataLine_NoCrash()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30\n" +
            "chr2\t200\trs2";   // file ends mid-line: only 3 columns, no newline
        string path = WriteTempVcf(vcf);

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.ParseFile(path).ToList();

        act.Should().NotThrow("a truncated final data line on disk must not crash ParseFile");
        records.Should().ContainSingle("only the complete first data line is a valid record")
            .Which.Pos.Should().Be(100);
    }

    #endregion

    #region MC — Non-integer POS: TryParse fails → line skipped, never a FormatException

    /// <summary>
    /// MC: data lines whose `POS` is non-numeric ("abc") or a float ("12.5"). POS is
    /// read with `int.TryParse` (VcfParser.cs lines 306–307), so a non-integer token
    /// makes TryParse return `false` and the line returns `null` (the documented
    /// "Non-integer POS → skipped", §6.1). The KEY guarantee is that this is a clean
    /// skip, NEVER an UNHANDLED `FormatException` a naive `int.Parse` would raise. We
    /// pin that every non-integer-POS line is dropped and only the integer-POS line
    /// survives.
    /// </summary>
    [Test]
    public void ParseVcf_NonIntegerPos_SkippedNotFormatException()
    {
        const string vcf =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\tabc\trs1\tA\tG\t.\t.\t.\n" +   // non-numeric POS → skipped
            "chr1\t12.5\trs2\tA\tG\t.\t.\t.\n" +  // float POS (not an int) → skipped
            "chr1\t400\trs3\tA\tG\t.\t.\t.\n";    // the one valid integer POS

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow(
            "non-integer POS is skipped via int.TryParse, never an unhandled FormatException");
        records.Should().ContainSingle("only the line with an integer POS parses")
            .Which.Pos.Should().Be(400);
    }

    #endregion

    #region MC — Missing #CHROM header: data lines still parse, samples silently dropped

    /// <summary>
    /// MC: a VCF with `##` meta lines and data lines but NO `#CHROM` column-header line.
    /// `Parse` does NOT require `#CHROM` — it skips `##`/`#` lines and parses every
    /// non-`#` line as data (VcfParser.cs lines 130–146). With no `#CHROM` line the
    /// `sampleNames` array stays null, so the FORMAT/sample columns are NOT materialized
    /// (the `fields.Length &gt; 9 &amp;&amp; sampleNames != null` gate, line 324), but the 8
    /// mandatory columns STILL parse. The DOCUMENTED behavior is therefore: data lines
    /// parse and samples are silently dropped — NOT a rejection and NOT a crash. We pin
    /// that exact shape (records present, `Samples == null`) rather than an idealized
    /// "reject when #CHROM missing" the parser does not implement.
    /// </summary>
    [Test]
    public void ParseVcf_MissingChromHeader_StillParsesDataLinesNoSamples()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
            // NO #CHROM line at all — the FORMAT+sample columns have no sample names
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30\tGT\t0/1\n" +
            "chr2\t200\trs2\tC\tT\t60\tPASS\t.\tGT\t1/1\n";

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow("a missing #CHROM header is tolerated, not a crash");
        records.Should().HaveCount(2,
            "the 8 mandatory columns parse even without a #CHROM header line");
        records[0].Chrom.Should().Be("chr1");
        records[0].Pos.Should().Be(100);
        records[0].Alt.Should().Equal("G");
        records.Should().AllSatisfy(r => r.Samples.Should().BeNull(
            "with no #CHROM header there are no sample names, so sample columns are not parsed"));
        VcfParser.GetGenotype(records[0], 0).Should().BeNull(
            "no samples were parsed, so there is no genotype to read");
    }

    #endregion

    #region MC — Missing ##fileformat: Parse ignores it; ParseWithHeader defaults to VCFv4.3

    /// <summary>
    /// MC: a VCF with NO `##fileformat=` meta line. In `Parse(string)` every `##` line
    /// is skipped uniformly — `##fileformat` carries no special status — so its absence
    /// has ZERO effect and the records still parse (VcfParser.cs lines 130–146). The
    /// repository contract is "parse anyway", NOT "reject". We pin that `Parse` still
    /// emits the records, and separately that `ParseWithHeader` DEFAULTS the
    /// `FileFormat` to `"VCFv4.3"` when no `##fileformat=` line is present (line 176) —
    /// a defined default, not a crash and not a rejection.
    /// </summary>
    [Test]
    public void ParseVcf_MissingFileFormat_StillParsesAndDefaultsHeader()
    {
        const string vcf =
            // NO ##fileformat line — straight to the #CHROM header
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30\n";

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow("a missing ##fileformat meta line is ignored, not a crash");
        records.Should().ContainSingle("the data line parses regardless of the missing ##fileformat")
            .Which.Pos.Should().Be(100);

        // ParseWithHeader supplies the documented default fileformat when none is present.
        var (header, headerRecords) = VcfParser.ParseWithHeader(vcf);
        header.FileFormat.Should().Be("VCFv4.3",
            "a missing ##fileformat defaults to VCFv4.3, the documented behavior (§5.2)");
        headerRecords.Should().ContainSingle();
    }

    #endregion

    #region MC — Invalid genotypes: stored verbatim, interpreted deterministically, no crash

    /// <summary>
    /// MC: the GT sample field carries INVALID genotypes — "9/9" (an allele index with
    /// no matching ALT), a malformed "1|" (trailing phase separator, no second allele),
    /// and a non-numeric ".|x". The parser stores GT as a RAW string in the sample dict
    /// with NO allele-range validation at parse time (ParseLine lines 331–337), so every
    /// invalid GT is parsed and preserved verbatim — never a crash. The downstream
    /// zygosity helpers (`IsHet`/`IsHomAlt`/`IsHomRef`, lines 529–563) then interpret
    /// them DETERMINISTICALLY: they split on `/`|`|` and a missing/out-of-shape genotype
    /// simply yields a defined bool rather than throwing. We pin that the invalid GTs
    /// parse, survive verbatim, and that every zygosity helper returns without crashing.
    /// </summary>
    [Test]
    public void ParseVcf_InvalidGenotypes_StoredVerbatimAndInterpretedDeterministically()
    {
        const string vcf =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
            "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t9/9\n" +   // allele index with no such ALT
            "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t1|\n" +    // malformed: trailing separator
            "chr1\t300\t.\tA\tG\t.\t.\t.\tGT\t.|x\n";    // non-numeric / missing allele

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow("invalid genotypes are stored as raw strings, not validated → no crash");
        records.Should().HaveCount(3, "every data line parses; the bad GT is not a rejection");

        // The raw GT strings survive verbatim in the sample dictionary.
        VcfParser.GetGenotype(records[0], 0).Should().Be("9/9");
        VcfParser.GetGenotype(records[1], 0).Should().Be("1|");
        VcfParser.GetGenotype(records[2], 0).Should().Be(".|x");

        // Every zygosity helper must return a DEFINED bool, deterministically, no throw.
        var zygosityProbe = () =>
        {
            foreach (var r in records)
            {
                _ = VcfParser.IsHomRef(r, 0);
                _ = VcfParser.IsHomAlt(r, 0);
                _ = VcfParser.IsHet(r, 0);
            }
        };
        zygosityProbe.Should().NotThrow(
            "zygosity helpers must interpret invalid genotypes deterministically, never crash");

        // "9/9": two equal non-zero, non-'.' alleles → homozygous-alt by the helper's rule;
        // never het (the two alleles are equal).
        VcfParser.IsHomAlt(records[0], 0).Should().BeTrue();
        VcfParser.IsHet(records[0], 0).Should().BeFalse();
        // ".|x": a '.' allele is "missing" → het/hom-alt are false, never a crash.
        VcfParser.IsHet(records[2], 0).Should().BeFalse(
            "a missing '.' allele cannot determine heterozygosity");
        VcfParser.IsHomAlt(records[2], 0).Should().BeFalse();
    }

    #endregion

    #region INJ/OVF — Huge alleles (1 MB REF + 1 MB ALT): linear time, no OOM/blow-up

    /// <summary>
    /// INJ / OVF: a single data line whose REF is 1 MB and whose ALT is a different 1 MB
    /// string (huge alleles). REF is stored as a plain string and ALT is `Split(',')`
    /// exactly once (VcfParser.cs line 311) — both O(n) over the allele length, no
    /// quadratic re-scan and no unbounded allocation. The record must parse with the
    /// full-length alleles intact, in linear time, with NO OutOfMemory and NO stall. We
    /// add a normal data line after the giant one to ALSO prove the parser recovers and
    /// still emits the following real record. `[CancelAfter]` converts any pathological
    /// stall into a deterministic failure.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void ParseVcf_HugeAlleles_CompletesWithoutBlowUpAndRecovers(CancellationToken token)
    {
        const int alleleLen = 1_000_000;                 // 1 MB REF and 1 MB ALT
        string hugeRef = new string('A', alleleLen);
        string hugeAlt = new string('C', alleleLen);
        string vcf =
            "##fileformat=VCFv4.3\n" +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
            "chr1\t100\t.\t" + hugeRef + "\t" + hugeAlt + "\t.\t.\t.\n" +
            "chr2\t200\trs2\tA\tG\t50\tPASS\t.\n";   // normal record after the giant one

        List<VcfParser.VcfRecord> records = new();
        var act = () => records = VcfParser.Parse(vcf).ToList();

        act.Should().NotThrow("a 1 MB REF/ALT must parse in linear time, no OOM/quadratic blow-up");
        records.Should().HaveCount(2, "the giant record parses and the following record survives");
        records[0].Ref.Length.Should().Be(alleleLen, "the full-length REF is preserved intact");
        records[0].Alt.Should().ContainSingle().Which.Length.Should().Be(alleleLen,
            "the full-length ALT is preserved intact");
        records[1].Pos.Should().Be(200, "the parser recovers and still emits the following record");
        token.IsCancellationRequested.Should().BeFalse("a huge allele must not stall the parser");
    }

    #endregion

    #endregion

    #region PARSE-GFF-001 — GFF parsing (FileIO)

    #region Positive sanity — a valid GFF3 parses to correct features with decoded attributes

    /// <summary>
    /// Positive control: a well-formed GFF3 (a `##gff-version 3` directive plus two
    /// 9-column data rows, one with a percent-encoded attribute value) must parse to
    /// exactly two `GffRecord`s with the correct seqid/source/type/coords/strand/phase
    /// and DECODED attributes. If this fails the whole GFF fuzz suite is meaningless —
    /// it proves the happy path (incl. percent-decoding) is wired up before we throw
    /// garbage at it.
    /// </summary>
    [Test]
    public void ParseGff_ValidGff3_ParsesToCorrectFeatures()
    {
        const string gff =
            "##gff-version 3\n" +
            "chr1\tEnsembl\tgene\t1000\t9000\t.\t+\t.\tID=gene1;Name=EDEN\n" +
            "chr1\tEnsembl\tmRNA\t1050\t9000\t0.95\t-\t0\tID=mrna1;Parent=gene1;Note=test%2Cnote\n";

        List<GffParser.GffRecord> records = GffParser.Parse(gff).ToList();

        records.Should().HaveCount(2, "both well-formed 9-column rows must be emitted");

        records[0].Seqid.Should().Be("chr1");
        records[0].Source.Should().Be("Ensembl");
        records[0].Type.Should().Be("gene");
        records[0].Start.Should().Be(1000);
        records[0].End.Should().Be(9000);
        records[0].Score.Should().BeNull("a '.' score is normalized to null (§6.1)");
        records[0].Strand.Should().Be('+');
        records[0].Phase.Should().BeNull("a '.' phase is normalized to null (§6.1)");
        records[0].Attributes["ID"].Should().Be("gene1");
        records[0].Attributes["Name"].Should().Be("EDEN");

        records[1].Type.Should().Be("mRNA");
        records[1].Score.Should().Be(0.95);
        records[1].Strand.Should().Be('-');
        records[1].Phase.Should().Be(0);
        records[1].Attributes["Parent"].Should().Be("gene1");
        records[1].Attributes["Note"].Should().Be("test,note",
            "the percent-encoded '%2C' in the attribute value decodes to a comma (INV-04)");
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: a valid GFF3 written to disk parses
    /// via `ParseFile` to the same correct records, pinning the StreamReader path (the
    /// common real-world entry point) and that an 8-column row (no attributes) is
    /// accepted with an empty attribute dictionary.
    /// </summary>
    [Test]
    public void ParseFileGff_ValidGff3_ParsesToCorrectFeatures()
    {
        const string gff =
            "##gff-version 3\n" +
            "scaf1\t.\texon\t5\t25\t.\t+\t.\tID=exon1\n" +
            "scaf1\t.\tCDS\t5\t25\t.\t+\t0\n";   // 8 columns, no attribute field
        string path = WriteTempGff(gff);

        List<GffParser.GffRecord> records = GffParser.ParseFile(path).ToList();

        records.Should().HaveCount(2);
        records[0].Type.Should().Be("exon");
        records[0].Attributes["ID"].Should().Be("exon1");
        records[1].Type.Should().Be("CDS");
        records[1].Phase.Should().Be(0);
        records[1].Attributes.Should().BeEmpty(
            "an 8-column row (no column 9) gets an empty attribute dictionary (§5.2)");
    }

    #endregion

    #region TF — Fewer than 8 columns: skipped, never IndexOutOfRange

    /// <summary>
    /// TF: a data line truncated to FEWER than the 8 mandatory columns (here 4:
    /// seqid/source/type/start). `ParseLine` splits on TAB and returns `null` when
    /// `fields.Length &lt; 8` (GffParser.cs lines 121–123) BEFORE ever indexing
    /// `fields[3]`…`fields[7]`, so the classic IndexOutOfRange trap CANNOT fire — the
    /// exact class of bug ANNOT-GFF-001 found in the OTHER GFF parser. Pinned: the
    /// 4-column row yields ZERO records and the following full 9-column row still
    /// parses, proving skip-and-recover, never a crash.
    /// </summary>
    [Test]
    public void ParseGff_FewerThanEightColumns_SkippedNoCrash()
    {
        const string gff =
            "chr1\tsrc\tgene\t100\n" +   // only 4 columns — truncated mid-record
            "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow(
            "a <8-column row is rejected on the field-count guard, never IndexOutOfRange");
        records.Should().ContainSingle("only the full 9-column row can be emitted")
            .Which.Attributes["ID"].Should().Be("g1");
    }

    #endregion

    #region MC — Non-integer / float start or end: skipped, never FormatException

    /// <summary>
    /// MC: a row whose `start` is non-numeric ("abc") and a row whose `end` is a FLOAT
    /// ("200.5"). Both coordinates are read with `int.TryParse` (lines 129–132); a
    /// non-integer token makes TryParse return `false`, so the line returns `null`
    /// (skipped) — the documented "start/end not an integer → skipped" of §3.3, NOT an
    /// unhandled `FormatException` from a `int.Parse` (the bug class ANNOT-GFF-001 hit).
    /// Pinned: both malformed-coordinate rows yield ZERO records, the valid row between
    /// them survives, and nothing throws.
    /// </summary>
    [Test]
    public void ParseGff_NonIntegerCoordinate_SkippedNoCrash()
    {
        const string gff =
            "chr1\tsrc\tgene\tabc\t200\t.\t+\t.\tID=bad1\n" +    // non-numeric start
            "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=ok\n" +      // valid
            "chr1\tsrc\tgene\t100\t200.5\t.\t+\t.\tID=bad2\n";   // float end

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow(
            "non-integer start/end is rejected via int.TryParse, never a FormatException");
        records.Should().ContainSingle("only the integer-coordinate row survives")
            .Which.Attributes["ID"].Should().Be("ok");
    }

    #endregion

    #region BE — Negative coordinate and start > end: parsed per the documented contract

    /// <summary>
    /// BE: a row with a NEGATIVE `start` (`-5`). `int.TryParse(NumberStyles.Integer)`
    /// accepts a leading `-`, and `ParseLine` has NO non-negative guard — §3.3 lists
    /// ONLY "&lt; 8 columns" and "non-integer start/end" as rejection rules. So the
    /// negative coordinate PARSES, producing a record with a negative `Start`. We pin
    /// the DOCUMENTED behavior exactly (parse, negative preserved), NOT an idealized
    /// "reject negatives" the FileIO parser does not implement — never a crash.
    /// </summary>
    [Test]
    public void ParseGff_NegativeCoordinate_ParsesWithNegativeStart()
    {
        const string gff = "chr1\tsrc\tgene\t-5\t200\t.\t+\t.\tID=neg\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("a negative coordinate is accepted by int.TryParse, no crash");
        records.Should().ContainSingle("there is no non-negative guard in ParseLine (§3.3)");
        records[0].Start.Should().Be(-5, "the negative start is preserved verbatim");
        records[0].End.Should().Be(200);
    }

    /// <summary>
    /// BE: a row with `start &gt; end` (a REVERSED interval, 500 &gt; 100). Unlike the BED
    /// parser, `GffParser.ParseLine` has NO `start &gt; end` rejection — INV-02
    /// (`start &lt;= end`) is a property of VALID input, not a parse-time guard. So the
    /// reversed interval PARSES (Start &gt; End) rather than being skipped. We pin the
    /// DOCUMENTED behavior (the reversed record is emitted, coordinates preserved),
    /// distinct from BED where start&gt;end is explicitly rejected — never a crash.
    /// </summary>
    [Test]
    public void ParseGff_StartGreaterThanEnd_ParsesReversedInterval()
    {
        const string gff = "chr1\tsrc\tgene\t500\t100\t.\t+\t.\tID=rev\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("there is no start>end guard in the GFF parser, so it parses");
        records.Should().ContainSingle();
        records[0].Start.Should().Be(500);
        records[0].End.Should().Be(100,
            "the reversed interval is preserved as-is (no start>end rejection, unlike BED)");
    }

    #endregion

    #region INJ — Invalid strand: first char preserved verbatim, empty defaults to '.'

    /// <summary>
    /// INJ: column 7 (strand) holding an INVALID value — one NOT in the spec set
    /// `+`/`-`/`.`/`?` ("X"), and an EMPTY strand field. `strand =
    /// fields[6].Length &gt; 0 ? fields[6][0] : '.'` (line 141) takes the FIRST char with
    /// NO membership check, so "X" is PRESERVED verbatim as `'X'` (not rejected, not
    /// normalized) and an empty column 7 becomes `'.'`. We pin BOTH shapes — the
    /// documented "store as-is" behavior — and assert no throw, so the lack of strand
    /// validation can never silently drift into a crash.
    /// </summary>
    [Test]
    public void ParseGff_InvalidStrand_PreservedVerbatimNoCrash()
    {
        const string gff =
            "chr1\tsrc\tgene\t100\t200\t.\tX\t.\tID=invalidStrand\n" +  // invalid strand 'X'
            "chr1\tsrc\tgene\t100\t200\t.\t\t.\tID=emptyStrand\n";       // empty strand column

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("an invalid/empty strand is stored as-is, never a crash");
        records.Should().HaveCount(2, "neither row is rejected — strand is not validated");
        records[0].Strand.Should().Be('X',
            "an invalid strand is preserved verbatim (first char, no membership check)");
        records[1].Strand.Should().Be('.',
            "an empty strand column defaults to '.' (the §2.2 undefined-strand symbol)");
    }

    #endregion

    #region MC — Malformed attributes: handled deterministically, never a crash

    /// <summary>
    /// MC: column 9 holding malformed attributes — a bare key with NO `=`
    /// (`flagonly`), a TRAILING `;`, an empty `=value` (no key), and a DUPLICATE key
    /// (`Name=first;Name=second`). The GFF3 attribute parser splits on `;` with
    /// `RemoveEmptyEntries`, then on the FIRST `=` with the `eqIdx &gt; 0` guard
    /// (GffParser.cs lines 183–193). So: the no-`=` part is silently DROPPED; the
    /// trailing `;` yields no empty entry; the empty `=value` (eqIdx == 0) is dropped;
    /// the duplicate key is LAST-WINS. We pin each exact outcome — the parser resolves
    /// malformed column 9 deterministically and never crashes (the documented tolerant
    /// behavior, §3.3).
    /// </summary>
    [Test]
    public void ParseGff_MalformedAttributes_HandledDeterministicallyNoCrash()
    {
        // flagonly: no '='. =orphan: empty key. Name appears twice → last wins.
        // Trailing ';' after the last pair must not create an empty/garbage entry.
        const string gff =
            "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1;flagonly;=orphan;Name=first;Name=second;\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("malformed attributes are parsed deterministically, never a crash");
        records.Should().ContainSingle();
        var attrs = records[0].Attributes;

        attrs["ID"].Should().Be("g1", "a well-formed key=value still parses");
        attrs.Should().NotContainKey("flagonly",
            "a part with no '=' is dropped by the eqIdx>0 guard");
        attrs.Should().NotContainKey("",
            "an empty '=value' (eqIdx==0) is dropped — no empty-key entry");
        attrs["Name"].Should().Be("second",
            "a duplicate key is last-wins (attributes[key] = value)");
    }

    /// <summary>
    /// MC companion: an EMPTY / whitespace-only attribute column (`.` or blank where
    /// callers sometimes place a placeholder). `ParseAttributes` early-returns an empty
    /// dictionary on `IsNullOrWhiteSpace` (lines 161–162), and a lone `.` simply has no
    /// `=` so contributes nothing. We pin that the record still parses with an empty
    /// attribute dictionary — no crash, no spurious attribute.
    /// </summary>
    [Test]
    public void ParseGff_EmptyAttributeColumn_YieldsEmptyAttributesNoCrash()
    {
        const string gff = "chr1\tsrc\tgene\t100\t200\t.\t+\t.\t.\n";  // attr col is just "."

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("a placeholder '.' attribute column must not crash");
        records.Should().ContainSingle();
        records[0].Attributes.Should().BeEmpty(
            "a lone '.' has no '=' so contributes no attribute");
    }

    #endregion

    #region INJ — Percent-encoding (THE key target): encoded separators decoded correctly

    /// <summary>
    /// INJ (THE KEY TARGET): attribute VALUES containing PERCENT-ENCODED reserved
    /// characters — `%3D` (=`=`), `%2C` (=`,`), `%3B` (=`;`), `%09` (=tab). The GFF3
    /// attribute parser splits column 9 on `;`/`=` on the RAW, still-encoded text
    /// BEFORE unescaping (GffParser.cs lines 183–193) and only THEN runs
    /// `UnescapeGff` = `Uri.UnescapeDataString` on each key/value. So an ENCODED
    /// separator inside a value survives the split intact and decodes correctly — a
    /// parser that decoded FIRST would mis-split on the embedded `=`/`;`. We pin:
    ///   • `Note=a%3Db`  → value `"a=b"`   (encoded `=` does NOT split the pair)
    ///   • `List=x%2Cy`  → value `"x,y"`   (encoded `,` decodes to a comma)
    ///   • `Semi=p%3Bq`  → value `"p;q"`   (encoded `;` does NOT split into two attrs)
    ///   • `Tab=u%09v`   → value `"u\tv"`  (encoded tab decodes inside the value)
    /// This is the central GFF3 robustness invariant (INV-04) and the exact behavior a
    /// naive `;`/`=` splitter gets wrong; we prove the parser gets it right.
    /// </summary>
    [Test]
    public void ParseGff_PercentEncodedSeparatorsInAttributeValues_DecodedCorrectly()
    {
        const string gff =
            "chr1\tsrc\tgene\t100\t200\t.\t+\t.\t" +
            "Note=a%3Db;List=x%2Cy;Semi=p%3Bq;Tab=u%09v\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow("percent-encoded separators must be handled, never a crash");
        records.Should().ContainSingle("the encoded ';'/'=' must NOT split the row into extra attrs");
        var attrs = records[0].Attributes;

        attrs.Should().HaveCount(4,
            "exactly four attributes — encoded separators do not create spurious entries");
        attrs["Note"].Should().Be("a=b", "%3D decodes to '=' AFTER the split on '='");
        attrs["List"].Should().Be("x,y", "%2C decodes to a comma");
        attrs["Semi"].Should().Be("p;q", "%3B decodes to ';' AFTER the split on ';'");
        attrs["Tab"].Should().Be("u\tv", "%09 decodes to a tab inside the value");
    }

    /// <summary>
    /// INJ companion: seqid/source/type columns also pass through `UnescapeGff`
    /// (lines 125–127), so a percent-encoded seqid (`chr%3A1` → `chr:1`) decodes; AND a
    /// MALFORMED percent-sequence (a lone `%`, `%ZZ`, a truncated `%2`) is passed
    /// through VERBATIM by `Uri.UnescapeDataString` (verified: it does NOT throw on a
    /// broken escape). We pin BOTH: a valid encoded seqid decodes, and a malformed
    /// escape in an attribute value survives verbatim without crashing the parser — so
    /// a broken escape can never become an unhandled exception.
    /// </summary>
    [Test]
    public void ParseGff_MalformedPercentEscape_PassedThroughVerbatimNoCrash()
    {
        const string gff =
            "chr%3A1\tsrc\tgene\t100\t200\t.\t+\t.\tBad=50%off;Trunc=ab%2;Hex=%ZZ\n";

        List<GffParser.GffRecord> records = new();
        var act = () => records = GffParser.Parse(gff).ToList();

        act.Should().NotThrow(
            "a malformed percent-escape is passed through verbatim by UnescapeDataString, never a crash");
        records.Should().ContainSingle();
        records[0].Seqid.Should().Be("chr:1", "a valid encoded seqid '%3A' decodes to ':'");
        var attrs = records[0].Attributes;
        attrs["Bad"].Should().Be("50%off", "a lone '%' that is not a valid escape is preserved verbatim");
        attrs["Trunc"].Should().Be("ab%2", "a truncated '%2' escape is preserved verbatim");
        attrs["Hex"].Should().Be("%ZZ", "an invalid-hex '%ZZ' escape is preserved verbatim");
    }

    #endregion

    #region Robustness — null/empty input yields no records, never a crash

    /// <summary>
    /// Robustness: null, empty, and directive/comment-only input must each yield ZERO
    /// records via the `IsNullOrEmpty` early-return (line 71–72) and the `##`/`#`
    /// skip (lines 96–111), never a crash. Pins the documented empty-input contract
    /// (§6.1) so a future refactor cannot turn empty input into an exception.
    /// </summary>
    [Test]
    public void ParseGff_NullEmptyAndDirectiveOnly_YieldNoRecordsNoCrash()
    {
        List<GffParser.GffRecord> nullRecords = new();
        var actNull = () => nullRecords = GffParser.Parse((string)null!).ToList();
        actNull.Should().NotThrow("null input is guarded by IsNullOrEmpty");
        nullRecords.Should().BeEmpty();

        GffParser.Parse("").ToList().Should().BeEmpty("empty input yields no records");

        const string directiveOnly =
            "##gff-version 3\n" +
            "# a comment line\n" +
            "##sequence-region chr1 1 1000\n";
        List<GffParser.GffRecord> dirRecords = new();
        var actDir = () => dirRecords = GffParser.Parse(directiveOnly).ToList();
        actDir.Should().NotThrow("directive/comment-only input must not crash");
        dirRecords.Should().BeEmpty("##/# lines are skipped — no data rows means no records");
    }

    #endregion

    #endregion

    #region PARSE-GENBANK-001 — GenBank parsing

    // A valid, fully-formed GenBank record (LOCUS + header + FEATURES + ORIGIN + //),
    // used by the positive-sanity tests and as the well-formed scaffold the fuzz
    // tests deliberately damage.
    private const string ValidGenBankRecord =
        "LOCUS       TEST001                  40 bp    DNA     linear   UNK 01-JAN-2024\n" +
        "DEFINITION  Test sequence for fuzzing.\n" +
        "ACCESSION   TEST001\n" +
        "VERSION     TEST001.1\n" +
        "KEYWORDS    test; fuzz.\n" +
        "SOURCE      Homo sapiens\n" +
        "  ORGANISM  Homo sapiens\n" +
        "            Eukaryota; Metazoa; Chordata.\n" +
        "FEATURES             Location/Qualifiers\n" +
        "     gene            1..20\n" +
        "                     /gene=\"testGene\"\n" +
        "     CDS             10..40\n" +
        "                     /product=\"test protein\"\n" +
        "ORIGIN      \n" +
        "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n" +
        "//\n";

    #region Positive sanity — a valid GenBank record parses to correct locus/features/sequence

    /// <summary>
    /// Positive control: a well-formed GenBank record (LOCUS + header + FEATURES +
    /// ORIGIN, terminated by `//`) must parse to exactly one record with the correct
    /// locus name, declared length, molecule type, the two FEATURES entries with their
    /// locations/qualifiers, and the uppercased ORIGIN sequence. If this fails the
    /// whole GenBank fuzz suite is meaningless — it proves the happy path is wired up
    /// before we throw garbage at it.
    /// </summary>
    [Test]
    public void ParseGenBank_ValidRecord_ParsesToCorrectLocusFeaturesAndSequence()
    {
        List<GenBankParser.GenBankRecord> records = GenBankParser.Parse(ValidGenBankRecord).ToList();

        records.Should().ContainSingle("a single well-formed record must be emitted");
        var rec = records[0];

        rec.Locus.Should().Be("TEST001");
        rec.SequenceLength.Should().Be(40, "the declared LOCUS length");
        rec.MoleculeType.Should().Be("DNA");
        rec.Topology.Should().Be("linear");
        rec.Definition.Should().Be("Test sequence for fuzzing.");
        rec.Accession.Should().Be("TEST001");

        rec.Features.Should().HaveCount(2, "the gene and CDS features");
        rec.Features[0].Key.Should().Be("gene");
        rec.Features[0].Location.Start.Should().Be(1);
        rec.Features[0].Location.End.Should().Be(20);
        rec.Features[0].Qualifiers["gene"].Should().Be("testGene");
        rec.Features[1].Key.Should().Be("CDS");
        rec.Features[1].Qualifiers["product"].Should().Be("test protein");

        rec.Sequence.Should().Be("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            "ORIGIN letters are extracted, digits/spaces stripped, uppercased");
    }

    /// <summary>
    /// Positive control over the FILE-PATH surface: the same valid record written to
    /// disk parses via `ParseFile` to the same single record, pinning the
    /// File.ReadAllText path (the common real-world entry point).
    /// </summary>
    [Test]
    public void ParseFileGenBank_ValidRecord_ParsesToCorrectRecord()
    {
        string path = WriteTempGenBank(ValidGenBankRecord);

        List<GenBankParser.GenBankRecord> records = GenBankParser.ParseFile(path).ToList();

        records.Should().ContainSingle();
        records[0].Locus.Should().Be("TEST001");
        records[0].Sequence.Should().Be("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        records[0].Features.Should().HaveCount(2);
    }

    #endregion

    #region MC — Missing LOCUS: a record with no LOCUS line is skipped (no crash)

    /// <summary>
    /// MC: a record body that has DEFINITION / FEATURES / ORIGIN sections but NO
    /// `LOCUS` line. Per GenBank_Parsing.md §3.3 (`Parse` parses a block "only when the
    /// resulting block begins with LOCUS", GenBankParser.cs line 102), a block that does
    /// NOT start with `LOCUS` is SKIPPED — the documented behavior, NOT a crash and NOT
    /// a record with a null locus. We pin ZERO records exactly so this gate cannot drift.
    /// </summary>
    [Test]
    public void ParseGenBank_MissingLocusLine_YieldsNoRecordsNoCrash()
    {
        const string noLocus =
            "DEFINITION  A record with no LOCUS keyword.\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene            1..20\n" +
            "ORIGIN      \n" +
            "        1 acgtacgtac\n" +
            "//\n";

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(noLocus).ToList();

        act.Should().NotThrow("a block not starting with LOCUS is skipped, not a crash");
        records.Should().BeEmpty("only blocks beginning with LOCUS are parsed (§3.3)");
    }

    /// <summary>
    /// MC companion: a two-record input where the FIRST block lacks `LOCUS` and the
    /// SECOND is a valid LOCUS record. The first is skipped, the second parses — pinning
    /// that a missing-LOCUS block does not poison the parse of the records around it.
    /// </summary>
    [Test]
    public void ParseGenBank_MissingLocusFollowedByValid_ParsesOnlyTheValidRecord()
    {
        string twoRecords =
            "DEFINITION  Orphan block with no LOCUS.\n" +
            "ORIGIN      \n" +
            "        1 acgtacgtac\n" +
            "//\n" +
            ValidGenBankRecord;

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(twoRecords).ToList();

        act.Should().NotThrow();
        records.Should().ContainSingle("only the LOCUS-bearing block is parsed")
            .Which.Locus.Should().Be("TEST001");
    }

    #endregion

    #region TF — Truncated features: a record cut off mid-FEATURES parses, never an IndexOutOfRange

    /// <summary>
    /// TF: a record whose FEATURES table is TRUNCATED — a feature line cut off with no
    /// location, AND a qualifier line with no `=value`, AND the file ends right after
    /// the qualifier (no ORIGIN, no `//`). `ParseFeatures` guards every line index
    /// (`line.Length > 5`, `line.Length > 21`, GenBankParser.cs lines 457/480/496), so:
    ///   • a feature key with no location → `currentLocation == ""` → a zeroed
    ///     `Location` via `ParseLocation("")`;
    ///   • a qualifier with no `=` → stored as value `"true"` (line 492).
    /// No `fields[n]` / Substring index can go out of bounds → the classic
    /// IndexOutOfRange-on-a-truncated-record trap CANNOT fire. We pin that the partial
    /// feature is emitted with a zeroed location and the bare qualifier, and that
    /// nothing crashes.
    /// </summary>
    [Test]
    public void ParseGenBank_TruncatedFeatures_ParsesPartialFeatureNoCrash()
    {
        // LOCUS present (so the block is parsed), FEATURES table cut off mid-feature:
        // a feature key with NO location, then a bare qualifier with no '=value',
        // then EOF (no ORIGIN, no '//').
        const string truncated =
            "LOCUS       TRUNC001                  40 bp    DNA     linear   UNK\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene\n" +                                  // feature key, NO location
            "                     /pseudo\n";                // qualifier, NO '=value', then EOF

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(truncated).ToList();

        act.Should().NotThrow("a truncated FEATURES table must not IndexOutOfRange");
        records.Should().ContainSingle("the LOCUS block is still parsed");
        var rec = records[0];
        rec.Locus.Should().Be("TRUNC001");
        rec.Features.Should().ContainSingle("the partial feature is emitted, not dropped");
        rec.Features[0].Key.Should().Be("gene");
        rec.Features[0].Location.Start.Should().Be(0, "a location-less feature yields a zeroed location");
        rec.Features[0].Location.End.Should().Be(0);
        rec.Features[0].Location.Parts.Should().BeEmpty();
        rec.Features[0].Qualifiers.Should().ContainKey("pseudo");
        rec.Features[0].Qualifiers["pseudo"].Should().Be("true",
            "a valueless qualifier is stored as \"true\" (line 492)");
    }

    #endregion

    #region MC — Invalid feature locations: parsed to a zeroed location, never an IndexOutOfRange

    /// <summary>
    /// MC: feature locations that are syntactically INVALID —
    ///   • `abc..xyz` (no digits at all);
    ///   • `complement(1..` (an unbalanced `complement(` with no closing paren).
    /// `ParseLocation` delegates to the REGEX-driven `ParseLocationParts`, which scans
    /// for `\d+(?:\.\.\d+)?` ranges and ignores everything else (SequenceFormatHelper.cs
    /// line 57). So `abc..xyz` matches NO range → a zeroed `Location` (Start==0, End==0,
    /// no Parts) with the RawLocation preserved verbatim; the unbalanced `complement(1..`
    /// flags `IsComplement` (a `StartsWith` check) and extracts the lone `1` as a single
    /// position — there is NO paren-matching that could throw on the missing `)`. We pin
    /// that invalid locations parse without an IndexOutOfRange and retain their raw text.
    /// </summary>
    [Test]
    public void ParseGenBank_InvalidFeatureLocations_ParseToZeroedLocationNoCrash()
    {
        const string invalidLoc =
            "LOCUS       BADLOC001                 40 bp    DNA     linear   UNK\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene            abc..xyz\n" +                  // no digits → no ranges
            "                     /gene=\"junk\"\n" +
            "     CDS             complement(1..\n" +            // unbalanced complement(
            "                     /product=\"trunc\"\n" +
            "ORIGIN      \n" +
            "        1 acgtacgtac\n" +
            "//\n";

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(invalidLoc).ToList();

        act.Should().NotThrow("invalid locations must parse to a zeroed location, never IndexOutOfRange");
        records.Should().ContainSingle();
        var features = records[0].Features;
        features.Should().HaveCount(2);

        // abc..xyz → no numeric ranges → zeroed location, raw string retained.
        features[0].Key.Should().Be("gene");
        features[0].Location.Start.Should().Be(0);
        features[0].Location.End.Should().Be(0);
        features[0].Location.Parts.Should().BeEmpty("a digitless location matches no ranges");
        features[0].Location.RawLocation.Should().Be("abc..xyz", "the raw location string is preserved");

        // complement(1.. → IsComplement flagged, the lone digit '1' extracted, no crash.
        features[1].Key.Should().Be("CDS");
        features[1].Location.IsComplement.Should().BeTrue("StartsWith(\"complement(\") flags it");
        features[1].Location.RawLocation.Should().Be("complement(1..",
            "the unbalanced raw location is preserved, no paren-match crash");
    }

    #endregion

    #region TF — No `//` terminator: the partial record is still emitted, never a hang

    /// <summary>
    /// TF — THE KEY LOOP/HANG CASE. A record that ends mid-ORIGIN with NO `//`
    /// terminator. The classic hang bug is a `while (line != "//")` loop with no EOF
    /// guard; this parser has NONE — `Parse` does `content.Split(new[]{"\n//"}, ...)`
    /// over the whole text in one O(n) pass (GenBankParser.cs line 97), so a record with
    /// no `//` is simply ONE block beginning with `LOCUS` and is STILL parsed. We pin
    /// that the unterminated record IS emitted with the correct locus and sequence (not
    /// dropped), over BOTH the string and `ParseFile` surfaces. `[CancelAfter]` is the
    /// loop/hang tripwire: a regression to an unbounded read loop would time out here
    /// rather than hang the run.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseGenBank_NoTerminator_EmitsPartialRecordAndDoesNotHang(CancellationToken token)
    {
        // A valid record body, but the file is cut off mid-ORIGIN with NO '//'.
        const string noTerminator =
            "LOCUS       NOTERM001                 40 bp    DNA     linear   UNK\n" +
            "DEFINITION  Record with no terminator.\n" +
            "FEATURES             Location/Qualifiers\n" +
            "     gene            1..20\n" +
            "                     /gene=\"openGene\"\n" +
            "ORIGIN      \n" +
            "        1 acgtacgtac gtacgtacgt";   // ends mid-sequence, NO trailing '//'

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(noTerminator).ToList();

        act.Should().NotThrow("a missing // must not crash");
        records.Should().ContainSingle("the unterminated record is still parsed (Split, not a loop)");
        records[0].Locus.Should().Be("NOTERM001");
        records[0].Features.Should().ContainSingle().Which.Qualifiers["gene"].Should().Be("openGene");
        records[0].Sequence.Should().Be("ACGTACGTACGTACGTACGT",
            "the partial ORIGIN sequence is still extracted");
        token.IsCancellationRequested.Should().BeFalse("a missing // terminator must not hang");
    }

    /// <summary>
    /// TF over the FILE-PATH surface: the same no-`//` record written to disk and read
    /// back via `ParseFile`, pinning that the File.ReadAllText path emits the partial
    /// record exactly as the string surface does. Carries `[CancelAfter]` as the
    /// loop/hang tripwire.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseFileGenBank_NoTerminator_EmitsPartialRecordAndDoesNotHang(CancellationToken token)
    {
        const string noTerminator =
            "LOCUS       FNOTERM01                 20 bp    DNA     linear   UNK\n" +
            "ORIGIN      \n" +
            "        1 acgtacgtac";   // file ends mid-record, NO '//'
        string path = WriteTempGenBank(noTerminator);

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.ParseFile(path).ToList();

        act.Should().NotThrow("ParseFile on a //-less file must not crash");
        records.Should().ContainSingle().Which.Locus.Should().Be("FNOTERM01");
        records[0].Sequence.Should().Be("ACGTACGTAC");
        token.IsCancellationRequested.Should().BeFalse("a missing // terminator must not hang ParseFile");
    }

    #endregion

    #region RB — Random bytes / null / empty: handled deterministically, never a crash or hang

    /// <summary>
    /// RB: a fixed-seed random-byte blob that does NOT begin with `LOCUS`. The whole
    /// blob is a single `Split` block; it fails the `StartsWith("LOCUS")` gate (line
    /// 102) → ZERO records, no crash, no hang. We guarantee the first bytes are not
    /// `"LOCUS"` and strip embedded `\n//` so the random bytes stay a single block. The
    /// LOCAL fixed seed makes the run byte-for-byte reproducible; `[CancelAfter]` is the
    /// hang tripwire.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ParseGenBank_RandomBytesNoLocus_YieldsNoRecordsAndDoesNotHang(CancellationToken token)
    {
        var rng = new Random(0x6B6E);          // local fixed seed — fully reproducible
        var sb = new StringBuilder("X");       // guarantee the block never starts with "LOCUS"
        for (int i = 0; i < 8192; i++)
        {
            char c = (char)rng.Next(0, 256);
            if (c == '/') c = 'z';             // keep it a single block (no '\n//' split / no // marker)
            sb.Append(c);
        }

        List<GenBankParser.GenBankRecord> records = new();
        var act = () => records = GenBankParser.Parse(sb.ToString()).ToList();

        act.Should().NotThrow("a random-byte block that does not start with LOCUS must not crash");
        records.Should().BeEmpty("only blocks beginning with LOCUS are parsed — random bytes yield none");
        token.IsCancellationRequested.Should().BeFalse("parsing random bytes must not hang");
    }

    /// <summary>
    /// Robustness: null and empty input must each yield ZERO records via the
    /// `IsNullOrEmpty` early-return (GenBankParser.cs lines 93–94), never a crash. Pins
    /// the documented empty-input contract (§6.1) so a future refactor cannot turn empty
    /// input into an exception.
    /// </summary>
    [Test]
    public void ParseGenBank_NullAndEmpty_YieldNoRecordsNoCrash()
    {
        List<GenBankParser.GenBankRecord> nullRecords = new();
        var actNull = () => nullRecords = GenBankParser.Parse(null!).ToList();
        actNull.Should().NotThrow("null input is guarded by IsNullOrEmpty");
        nullRecords.Should().BeEmpty();

        GenBankParser.Parse("").ToList().Should().BeEmpty("empty input yields no records");
    }

    #endregion

    #endregion
}
