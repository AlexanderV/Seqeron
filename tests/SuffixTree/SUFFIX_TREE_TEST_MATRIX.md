# Suffix Tree Test Matrix

## Scope
This matrix tracks `ISuffixTree` contract coverage for:
- `SuffixTree.Tests` (in-memory implementation)
- `SuffixTree.Persistent.Tests` (persistent implementation: heap + mapped file)

Status labels:
- `Strong`: functional + oracle/parity/contract assertions
- `Medium`: functional checks exist, but oracle depth is limited

## ISuffixTree Contract Coverage

| API | In-Memory Coverage | Persistent Coverage | Strength | Notes |
|---|---|---|---|---|
| `Text` | `Core/BuildTests`, `Core/InvariantTests` | `Parity/ParityTests`, `Safety/DisposeBehaviorTests` | Strong | Includes disposed-state contract. |
| `NodeCount` / `LeafCount` / `MaxDepth` / `IsEmpty` | `Core/InvariantTests`, `Core/BuildTests`, `Properties/PropertyBasedTests` | `Parity/*`, `Validation/*`, `Safety/DisposeBehaviorTests` | Strong | Includes invariant and parity checks. |
| `Contains(string)` | `Search/ContainsTests`, `Algorithms/*`, `Properties/PropertyBasedTests` | `Parity/*`, `Validation/*`, `Format/*`, `Safety/*` | Strong | Cross-checked with `Count` and `FindAll`. |
| `Contains(ReadOnlySpan<char>)` | `Search/ContainsTests`, `Algorithms/SuffixTreePropertyTests` | `Validation/SpanOverloadContractTests` | Strong | Span/string parity on fixed + random inputs. |
| `FindAllOccurrences(string)` | `Search/FindAllOccurrencesTests`, `Algorithms/BruteForceVerificationTests` | `Parity/*`, `Validation/SpanOverloadContractTests`, `Format/*` | Strong | Includes brute-force oracle and parity. |
| `FindAllOccurrences(ReadOnlySpan<char>)` | `Search/FindAllOccurrencesTests`, `Algorithms/SuffixTreePropertyTests` | `Validation/SpanOverloadContractTests` | Strong | Span/string parity + brute-force contract. |
| `CountOccurrences(string)` | `Search/CountOccurrencesTests`, `Core/InvariantTests`, `Properties/PropertyBasedTests` | `Parity/*`, `Validation/*`, `Format/*` | Strong | Consistency with `FindAll` is enforced. |
| `CountOccurrences(ReadOnlySpan<char>)` | `Search/CountOccurrencesTests`, `Algorithms/SuffixTreePropertyTests` | `Validation/SpanOverloadContractTests` | Strong | Span/string parity + random checks. |
| `LongestRepeatedSubstring()` | `Algorithms/LongestRepeatedSubstringTests`, `Algorithms/BruteForceVerificationTests`, `Properties/PropertyBasedTests` | `Parity/*`, `Format/*`, `Serialization/*` | Strong | Includes maximality + repetition invariants. |
| `LongestRepeatedSubstringMemory()` | `Algorithms/LongestRepeatedSubstringTests`, `Core/BuildTests` | `Validation/LcsAndEnumerationContractTests` | Strong | Interface default method validated on persistent trees. |
| `GetAllSuffixes()` | `Algorithms/SuffixEnumerationTests`, `Properties/PropertyBasedTests` | `Validation/LcsAndEnumerationContractTests`, `Validation/TreeContractTests` | Strong | Compared to reference + enumeration contract. |
| `EnumerateSuffixes()` | `Algorithms/SuffixEnumerationTests` | `Validation/LcsAndEnumerationContractTests`, `Validation/TreeContractTests` | Strong | Equality with `GetAllSuffixes` is enforced. |
| `LongestCommonSubstring(string)` | `Algorithms/LongestCommonSubstringTests`, `Algorithms/LongestCommonSubstringOracleTests` | `Parity/*`, `Validation/SpanOverloadContractTests`, `Validation/LcsAndEnumerationContractTests` | Strong | Length oracle + parity with reference implementation. |
| `LongestCommonSubstring(ReadOnlySpan<char>)` | `Algorithms/LongestCommonSubstringOracleTests` | `Validation/SpanOverloadContractTests`, `Validation/LcsAndEnumerationContractTests` | Strong | Span/string consistency asserted. |
| `LongestCommonSubstringInfo(string)` | `Algorithms/LongestCommonSubstringTests`, `Algorithms/LongestCommonSubstringOracleTests` | `Validation/LcsAndEnumerationContractTests` | Strong | Position validity and consistency with `FindAllLCS`. |
| `FindAllLongestCommonSubstrings(string)` | `Algorithms/LongestCommonSubstringTests`, `Algorithms/LongestCommonSubstringOracleTests` | `Validation/LcsAndEnumerationContractTests`, `Safety/DisposeBehaviorTests` | Strong | Returned positions validated against reported substring. |
| `PrintTree()` | `Core/DiagnosticsTests` | `Validation/DiagnosticsContractTests` | Strong | Verifies depth/indent format, deterministic ordering, node/leaf line cardinality, empty-tree contract. |
| `Traverse(ISuffixTreeVisitor)` | `Core/DiagnosticsTests` | `Validation/TreeContractTests`, `Parity/TopologyParityTests`, `Format/HybridTransitionZoneTests` | Strong | Depth semantics + branch balancing + topology parity. |
| `FindExactMatchAnchors(...)` | `Algorithms/ExactMatchAnchorTests` | `Parity/*`, `Format/*`, `Serialization/*` | Strong | Validated for semantics and cross-implementation parity. |

## Persistent-Specific Additional API

| API | Coverage | Strength | Notes |
|---|---|---|---|
| `LongestCommonSubstringInfo(ReadOnlySpan<char>)` | `Validation/LcsAndEnumerationContractTests` | Strong | Span/string parity + reference parity. |
| Hybrid transition/jump-table behavior | `Format/HybridTransitionZoneTests` | Strong | Boundary sweeps, round-trip load, cross-zone suffix-link regression. |
| Serialization/hash contracts | `Serialization/*` | Strong | Logical hash parity, truncation guards, import/export contracts. |
| Storage lifecycle/safety | `Safety/*` | Strong | Dispose behavior, concurrency, read-only provider, overflow/guard checks. |

## Residual Risk
- `PrintTree()` still intentionally avoids full golden-snapshot locking to keep formatter evolution low-friction; structure and invariants are now asserted directly.
