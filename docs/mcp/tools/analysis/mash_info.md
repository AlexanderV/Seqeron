# mash_info

Mash 2.3 'mash info', byte for byte (120 sketch files x 5 output modes equal): the header (hash function and seed, k-mer size with the hash width, the alphabet with '(canonical)' / '(case-sensitive)', the target sketch size, the number of sketches) and the table of the sketches (hashes, length, id, comment; an empty field is '-'), or with headerOnly only the header, with tabular the tab separated table, with counts the hash count histogram of each read-set sketch, with dump the JSON dump of the sketches and hashes. The options headerOnly / tabular / counts / dump exclude each other as in the program (the errors are the program's messages).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `mash_info` |
| **Method ID** | `MashCommands.Info` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Mash 2.3 'mash info', byte for byte (120 sketch files x 5 output modes equal): the header (hash function and seed, k-mer size with the hash width, the alphabet with '(canonical)' / '(case-sensitive)', the target sketch size, the number of sketches) and the table of the sketches (hashes, length, id, comment; an empty field is '-'), or with headerOnly only the header, with tabular the tab separated table, with counts the hash count histogram of each read-set sketch, with dump the JSON dump of the sketches and hashes. The options headerOnly / tabular / counts / dump exclude each other as in the program (the errors are the program's messages).

## Core Documentation Reference

- Algorithm: [Mash_Commands.md](../../../algorithms/Alignment/Mash_Commands.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sketch` | string | Yes | The sketch file: 'msh:' + base64 (the output of mash_sketch / mash_paste). |
| `headerOnly` | boolean | No | -H: only the header. (default: false) |
| `tabular` | boolean | No | -t: tabular output (no header). (default: false) |
| `counts` | boolean | No | -c: hash count histograms (only for read-set sketches). (default: false) |
| `dump` | boolean | No | -d: dump the sketches as JSON. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The text mash info prints. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (the program's message for incompatible options or a sketch file that is not a .msh). |

## Examples

### Example 1: Show the header and the sketches

**User Prompt:**
> What is in this sketch file?

**Expected Tool Call:**
```json
{"tool": "mash_info", "arguments": {"sketch": "msh:..."}}
```

**Response:**
```json
{"text": "Header:\n..."}
```

## See Also

- [mash_dist](mash_dist.md) - distances
- [mash_screen](mash_screen.md) - containment
