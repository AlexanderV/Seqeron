# A2-12b-s6b run 2 — stopped by the user, notes for the fix

Run 2 of `fetch_dalilite_snapshots.py` (package as of `f115ac3e0`) started 2026-11-26 01:20 local and was **stopped by the user at 18:20** (`docker stop`, ~17 h) while in stage D. The archive `returned/A2-12b-s6b-out.run2.zip` (2.4 MB) is the `out/` at the stop, with one exception: `out/found/.tmp` was removed. That file was a 150 MB partial `dump.tar.gz`, and with it the zip was 153 MB, over GitHub's 100 MB limit. The original zip and the partial file are kept outside the repo. `report.txt` / `summary.json` do not exist (written only at the end of the script); `requests.log` / `fetch.log` are complete up to the stop.

## What run 2 established

- **Stages A + B** (82 requests, ~15 min): 24 distinct lead pages, **0 candidate links**; CDX of the 20 bases. `CDX http://ekhidna.biocenter.helsinki.fi/dali/` returned **exactly 5000 lines = the `limit=5000` cap: truncated**, entries beyond it were never seen.
- **Stage C done**: `STAGE C 1548 URLs: 9 interesting, 2 FAILED`. 1539 HEAD probes answered 404. **None of the 77 release names × 20 bases is archived** (nearest capture to 2008). The 9 "interesting" ones are exactly the 9 CDX candidates (below), not name hits.
- **Kept** (gzip, distinct sha256):
  - `01_DaliLite.tar.gz` ← `ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v3/DaliLite.tar.gz`, capture 20170709235707, 1 367 393 B, sha256 `85f6dd8d…ef81e6`; top dir `DaliLite_3.3/`, its `CHANGES` lists 3.3, 3.2, 3.1, 3.0, 2.4.5, 2.4.4 (version history only, not older code).
  - `02_DaliLite.v4.tar.gz` ← `ekhidna.biocenter.helsinki.fi/dali/DaliLite.v4.tar.gz`, capture 20241231122421, 800 664 B, sha256 `6267acef…0cb5cb`; top dir `DaliLite.v4.1/`.
  - The text inspection of the kept archives (`out/inspect/`, `report.txt`) never ran.
- **FAILED (2)**: `ekhidna2…/dali/AF-Digest.tar.gz` and `ekhidna2…/dali/pdb_and_scope.tar`: every HEAD timed out at 90 s on a 302 (the archive never finished redirecting). Both are database files, not DaliLite releases.
- **Stage D at the stop**: `v3/DaliLite.tar.gz` for all 11 other years → 11 × `KEPT duplicate` (always the same 2017 capture); `DAT.tar.gz` 2002 = full retry cycle, 2003–2013 = 404; `dump.tar.gz` 2002 in its retry cycle. Not reached: the rest of dump, GTG, AF-Digest ×2, v4, benchmark, pdb_and_scope, the second ask of the 2 FAILED.

## Where the 17 hours went

38 GET timeouts × 900 s (~9.5 h) + 14 HEAD timeouts × 90 s + 66 back-offs totalling 293 min. Practically all of it on 7 multi-GB DALI database files that the CDX filter admitted: `dali_lite/downloads/v3/{DAT,dump,GTG}.tar.gz`, `dali/AF-Digest.tar.gz` (both hosts), `ekhidna2…/dali/{benchmark.tar.gz,pdb_and_scope.tar}`. The archive serves them at ~0.35 MB/s; each GET got 250–520 MB before `--max-time 900`, then was retried 6 more times (~2.3 h per file, per year).

## Defects and proposed fixes (`fetch_dalilite_snapshots.py`)

1. **Candidate filter admits any archive** (lines 158 and 170): `ARCHIVE_EXT.search(x) or re.search('dali', x) and ARCHIVE_EXT.search(x)`. The first branch alone accepts every `.tar.gz/.tar/.gz` under the bases; the second is redundant. Require the basename to look like a DaliLite release, e.g. `re.search(r'dali.?lite', basename, re.I) and ARCHIVE_EXT.search(basename)`, and list the rejected ones in the log.
2. **No size guard.** The real releases are 0.8–1.4 MB. Add `--max-filesize 20M` to the GET in `probe()` (curl rc 63 → a final, non-transient status `TOO-LARGE`, logged with the size), and/or read `Content-Length` from the HEAD (`once()` already captures the headers with `-D -` but discards them) and skip the GET above the cap.
3. **Timeouts with data flowing are retried as network faults** (`curl()`, `transient = rc != 0 or …`). rc 28 with `http=200` and `bytes>0` means slow/huge, not unreachable: do not retry it (or once). HEAD rc 28 on a 302 after 90 s: at most one retry.
4. **Stage D re-probes every candidate unconditionally**: `interesting = [u for u in urls if u in candidates or …]`. Only URLs that answered 2xx at 2008 and were kept (or at least not FAILED/TOO-LARGE) should go to D.
5. **Years resolve to the same capture.** `web/<year>id_/` returns the nearest capture, so all 11 years of `v3/DaliLite.tar.gz` downloaded the same 20170709235707 file. Either skip the GET when the HEAD's final URL (capture timestamp) was already fetched, or drop the year sweep and fetch each distinct capture listed by CDX (the timestamps are already in `out/cdx/*.txt`) once.
6. **CDX truncation**: `limit=5000` was hit for `ekhidna.biocenter.helsinki.fi/dali/`. Paginate (`showResumeKey=true` / `resumeKey`, or `page=`) or query narrower prefixes, and log "TRUNCATED" when the cap is reached.
7. **Nothing on the console for hours**: `in_container.sh` sends all script output to `out/fetch.log`, so the user's terminal (and `run_console.log`, 0 bytes) looks hung. Tee the log to stdout, or print a progress line per stage / every N requests to stderr.
8. **No summary on interruption**: `report.txt` / `summary.json` and the inspection are written only at the very end. Write the SUMMARY line and each KEPT archive's inspection as they happen, and handle SIGTERM/SIGINT (`docker stop` sends SIGTERM to bash as PID 1 and SIGKILLs after 10 s; run with `--init`, `exec` the python, or install a handler that flushes the summary).
9. **`pack()` zips the partial download**: exclude `out/found/.tmp` (`zip -qr … out -x 'out/found/.tmp'`), so that an interrupted GET cannot push the archive past the 100 MB GitHub limit.

## What is still open for the question

Stage C already answers it for the 77 names × 20 bases (nothing archived under those names). Open: the truncated CDX of `ekhidna…/dali/` (defect 6), and running the text inspection of the two kept archives. A run 3 with fixes 1–6 would take about 2 hours instead of days: ~82 requests for A + B, CDX pagination, ~1550 HEAD probes at 2.5 s (~1.5 h), small GETs only.
