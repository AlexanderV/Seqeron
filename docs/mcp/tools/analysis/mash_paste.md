# mash_paste

Mash 2.3 'mash paste', byte for byte (35 random runs equal): one sketch file of the sketches of several. The parameters (k, sketch size, seed, alphabet, strand mode, case) are those of the first file; a later file with another alphabet, seed, k or strand mode or a smaller sketch size is skipped (the program warns on the error stream), a larger sketch is reduced to the first file's size; hash counts of read-set sketches are not kept.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_paste` |
| **Method ID** | `MashCommands.Paste` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash paste', byte for byte (35 random runs equal): one sketch file of the sketches of several. The parameters (k, sketch size, seed, alphabet, strand mode, case) are those of the first file; a later file with another alphabet, seed, k or strand mode or a smaller sketch size is skipped (the program warns on the error stream), a larger sketch is reduced to the first file's size; hash counts of read-set sketches are not kept.

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sketches` | array<string> | Yes | The sketch files: 'msh:' + base64 each. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `msh` | string | msh: + base64 of the pasted file. |
| `bytes` | integer | Size of the file. |
| `sketches` | integer | Number of sketches. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (the program's message for incompatible options or a sketch file that is not a .msh). |

## Examples

### Example 1: Paste two sketch files

**User Prompt:**
> Combine these two sketch files.

**Expected Tool Call:**
```json
{"tool": "mash_paste", "arguments": {"sketches": ["msh:...", "msh:..."]}}
```

**Response:**
```json
{"msh": "msh:...", "bytes": 300, "sketches": 2}
```

## See Also

- [mash_dist](mash_dist.md) - distances
- [mash_screen](mash_screen.md) - containment
