# Test Specification: <TEST_UNIT_ID>

**Test Unit ID:** <ID>
**Area:** <Area>
**Algorithm:** <n>
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** <date>

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

<!-- Every source MUST have a verifiable link. DOI preferred, PubMed/PMC acceptable. -->

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| | | | | |

### 1.2 Key Evidence Points

1. <point> — <source citation>

### 1.3 Documented Corner Cases

<!-- From evidence. If none found: "No authoritative sources explicitly specify corner cases for <X>." -->

### 1.4 Known Failure Modes / Pitfalls

1. <failure mode> — <source citation>

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| | | | |

<!-- Type values: -->
<!-- **Canonical** — deep evidence-based testing -->
<!-- **Delegate** — smoke verification only (1–2 tests proving delegation) -->
<!-- **Internal** — tested indirectly via canonical methods -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | <statement> | Yes/No | <source or **ASSUMPTION**> |

---

## 4. Test Cases

<!-- For complex inputs (sequences, matrices): describe briefly here, -->
<!-- put full data in test code or a data file, reference it. -->

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | <n> | <what and why> | <outcome> | <source or **ASSUMPTION**> |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- <where tests were found, file paths>

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| | | |

<!-- Status values: ✅ Covered, ⚠ Weak, ❌ Missing, 🔁 Duplicate -->

### 5.3 Consolidation Plan

- **Canonical file:** `<path>` — <what goes here>
- **Remove:** <file/tests and reason>

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| | | |

### 5.5 Phase 7 Work Queue

<!-- Before writing code, list every ❌ and ⚠ from §5.2. -->
<!-- Implement each, mark ✅ Done or ⛔ Blocked. All must be resolved. -->

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| | | | | |

**Total items:** <N>
**✅ Done:** <N> | **⛔ Blocked:** <N> | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

<!-- After Phase 7, re-audit every row from §5.2. ALL must be ✅. -->
<!-- If any remain ❌ or ⚠, the Test Unit cannot be ☑ Complete. -->

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| | | |

---

## 6. Assumption Register

**Total assumptions:** <N>

<!-- List every ASSUMPTION referenced in this spec. -->
<!-- Zero assumptions is the goal; each one is a risk. -->

| # | Assumption | Used In |
|---|-----------|---------|
| | | |

---

## 7. Open Questions / Decisions

<!-- Unresolved items. "None" only if genuinely nothing is open. -->

1. <question or decision needed>