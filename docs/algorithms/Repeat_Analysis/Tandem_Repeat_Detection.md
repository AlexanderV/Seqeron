# Tandem Repeat Detection

## Algorithm Overview

**Test Unit ID:** REP-TANDEM-001
**Area:** Repeat Analysis
**Complexity:** O(n²) for FindTandemRepeats

---

## Definition

A **tandem repeat** is a pattern of DNA where one or more nucleotides are repeated consecutively in a sequence. The repetitions are directly adjacent to each other, e.g., `ATTCGATTCGATTCG` contains the unit `ATTCG` repeated three times.

### Terminology (from Wikipedia)

| Term | Definition | Unit Length |
|------|------------|-------------|
| **Microsatellite (STR)** | Short Tandem Repeat | 1–6 bp (sometimes 1–10 bp) |
| **Minisatellite** | Variable Number Tandem Repeat (VNTR) | 10–60 bp |
| **Macrosatellite** | Large tandem repeats | ~1,000+ bp |
| **Dinucleotide repeat** | Two nucleotides repeated | 2 bp (e.g., ACACAC) |
| **Trinucleotide repeat** | Three nucleotides repeated | 3 bp (e.g., CAGCAGCAG) |

### Classification by Unit Length

| Type | Unit Length | Example |
|------|-------------|---------|
| Mononucleotide | 1 bp | AAAAA |
| Dinucleotide | 2 bp | ATATAT |
| Trinucleotide | 3 bp | CAGCAGCAG |
| Tetranucleotide | 4 bp | GATAGATAGATA |
| Pentanucleotide | 5 bp | AATGGAATGG |
| Hexanucleotide | 6 bp | TTAGGGTTAGGG |

---

## Biological Significance

### Medical Relevance

Tandem repeats are implicated in:
- **Trinucleotide repeat disorders** (e.g., Huntington's disease - CAG expansions, Fragile X syndrome)
- **Microsatellite instability** in cancers (particularly colorectal cancer)
- ~8% of the human genome consists of tandem repeats
- Linked to >50 human diseases

### Applications

- **Forensic DNA profiling**: STRs are standard markers (CODIS uses 13+ loci)
- **Paternity testing**: Highly polymorphic nature enables kinship analysis
- **Population genetics**: Used for diversity studies
- **Cancer diagnosis**: Microsatellite instability detection

---

## Mutation Mechanism

Tandem repeats mutate primarily through **replication slippage**:
1. DNA polymerase slips during replication of repetitive sequences
2. Template strand loops cause repeat loss
3. New strand loops cause repeat gain
4. Rate: ~1 slippage per 1,000 generations (3 orders of magnitude higher than point mutations)

---

## Algorithm Details

### FindTandemRepeats (Canonical)

**Location:** `GenomicAnalyzer.FindTandemRepeats`

**Signature:**
```csharp
IEnumerable<TandemRepeat> FindTandemRepeats(
    DnaSequence sequence,
    int minUnitLength = 2,
    int minRepetitions = 2)
```

**Algorithm:**
1. For each unit length from `minUnitLength` to `sequence.Length / minRepetitions`:
2. For each starting position where a repeat could fit:
3. Extract candidate unit and count consecutive occurrences
4. If repetitions ≥ `minRepetitions`, yield result
5. Skip to end of detected tandem to avoid overlapping reports

**Complexity:** O(n² × m) where n = sequence length, m = max unit length

**Output:** `TandemRepeat` struct containing:
- `Unit`: The repeated pattern
- `Position`: Start position (0-based)
- `Repetitions`: Number of consecutive repeats
- `TotalLength`: Unit.Length × Repetitions
- `FullSequence`: The complete repeated sequence

### GetTandemRepeatSummary (Summary/Wrapper)

**Location:** `RepeatFinder.GetTandemRepeatSummary`

**Signature:**
```csharp
TandemRepeatSummary GetTandemRepeatSummary(
    DnaSequence sequence,
    int minRepeats = 3)
```

**Behavior:** Delegates to `FindMicrosatellites` (1-6 bp units) and aggregates results.

**Output:** `TandemRepeatSummary` record containing:
- Total repeat count
- Total bases covered
- Percentage of sequence
- Counts by type (mono-, di-, tri-, tetra-nucleotide)
- Longest repeat identified
- Most frequent unit

---

## Implementation Notes

### Current Implementation

- Uses brute-force O(n²) approach (simple, correct, suitable for moderate sequences)
- Reports non-overlapping results by skipping to end of detected repeats
- Does not use suffix tree optimization (ASSUMPTION: direct string comparison chosen for clarity)

### Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| Empty sequence | Return empty enumerable |
| No repeats found | Return empty enumerable |
| Entire sequence is one repeat | Return single result covering full sequence |
| Overlapping repeat patterns | Implementation skips after detection (non-overlapping output) |
| minRepetitions < 2 | Minimum valid value is 2 (at least one repetition) |
| Unit length > sequence/minReps | No results possible |

### Invariants

1. `repetitions >= minRepetitions` for all results
2. `unit.Length >= minUnitLength` for all results
3. `position + unit.Length × repetitions <= sequence.Length`
4. `TotalLength == unit.Length × repetitions`
5. Results represent actual consecutive occurrences in sequence

---

## Sources

1. **Wikipedia - Tandem repeat**: https://en.wikipedia.org/wiki/Tandem_repeat
   - Definition, terminology, detection methods, biological significance

2. **Wikipedia - Microsatellite (Short tandem repeat)**: https://en.wikipedia.org/wiki/Microsatellite
   - Classification (1-6 bp), mutation mechanisms (slippage), applications (forensics, disease)
   - Richard et al. (2008) - Comparative genomics of DNA repeats

3. **Richard GF, Kerrest A, Dujon B (2008)**: "Comparative genomics and molecular dynamics of DNA repeats in eukaryotes" - Microbiology and Molecular Biology Reviews 72(4):686-727
   - Authoritative review on STR structure, function, and variation

4. **Benson G (1999)**: Tandem Repeats Finder (TRF) algorithm - conceptual reference for detection approaches

---

## Version History

| Date | Version | Author | Changes |
|------|---------|--------|---------|
| 2026-01-22 | 1.0 | Algorithm QA | Initial documentation |
