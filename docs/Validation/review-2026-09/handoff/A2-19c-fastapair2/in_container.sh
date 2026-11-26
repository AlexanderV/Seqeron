#!/usr/bin/env bash
# Runs inside the image: for every set of inputs/sets.tsv (1) `mafft <options> in.fa` WITHOUT --quiet with FASTA_4_MAFFT = a logging wrapper (rows, names, return code, standard error, tree / distance files, every raw fasta34 call: arguments, query, library, raw output),
# (2) dndfast7 [-A] run directly on the same sequences in an empty directory: its hat2, hat3 and standard error.
# RUN 4: only the sets SETS_FROM .. SETS_TO (default 254 .. 288) are run; every set runs in its own try/except (an error is recorded in results.tsv column 10 and in sets/<k>.error.txt, the run goes on); the raw standard OUTPUT of mafft
# is stored as it is (sets/<k>.stdout.gz: CLUSTAL, PHYLIP, --text, any non-FASTA output), as are all files the run leaves in the working directory (sets/<k>.files.json.gz).
set -uo pipefail
cd /work
rm -rf out/sets out/results.tsv out/progress.log; mkdir -p out/sets
{ mafft --version 2>&1 | head -1; command -v fasta34; sha256sum "$(command -v fasta34)"; uname -m; tail -3 /opt/fasta34.build.log; ls /opt/FASTA34_BUILD_FAILED 2>&1; ls -l /usr/lib/mafft/lib/mafft/dndfast7 2>&1; } > out/versions.txt 2>&1
export FASTA_4_MAFFT=/work/fasta34-log
python3 -I - <<'PY' &
import glob, gzip, hashlib, io, json, os, shutil, signal, subprocess, sys, tarfile, tempfile, traceback


TIMEOUT = int(os.environ.get('SET_TIMEOUT', '1200'))     # seconds per program run: a hang is recorded (rc 124), it never blocks the run
SETS_FROM = int(os.environ.get('SETS_FROM', '254'))       # run 4: sets 1 - 253 (complete runs in run 3, raw calls and direct runs in run 2) are NOT repeated
SETS_TO = int(os.environ.get('SETS_TO', '288'))
DIRECT_FROM = int(os.environ.get('DIRECT_FROM', str(SETS_FROM)))    # first set whose direct dndfast7 run is made
MAXFILE = 4 << 20


class Stop(Exception):
    pass


def on_term(signum, frame):
    raise Stop()


signal.signal(signal.SIGTERM, on_term)      # docker stop (run.sh after Ctrl-C): finish cleanly, the zip is written by run.sh afterwards


def run_logged(argv, cwd, stdin_path=None):
    """Runs argv in its own process group; a hang or a missing program is a recorded result (rc 124 / 127), never an exception of the script."""
    d = os.path.join(cwd, 'fastalog')
    shutil.rmtree(d, ignore_errors=True)
    env = dict(os.environ, FASTALOG=d)
    stdin = open(stdin_path, 'rb') if stdin_path else subprocess.DEVNULL
    try:
        pr = subprocess.Popen(argv, cwd=cwd, env=env, stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    except OSError as e:
        return subprocess.CompletedProcess(argv, 127, b'', ('@@cannot start %s: %s@@\n' % (argv[0], e)).encode()), d
    try:
        o, e = pr.communicate(timeout=TIMEOUT)
        p = subprocess.CompletedProcess(argv, pr.returncode, o, e)
    except subprocess.TimeoutExpired:
        os.killpg(pr.pid, signal.SIGKILL)
        o, e = pr.communicate()
        p = subprocess.CompletedProcess(argv, 124, o or b'', (e or b'') + b'\n@@TIMEOUT after %d s@@\n' % TIMEOUT)
    except BaseException:
        os.killpg(pr.pid, signal.SIGKILL)
        pr.communicate()
        raise
    finally:
        if stdin_path:
            stdin.close()
    return p, d


def pack_dir(d):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as t:
        for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            t.add(os.path.join(d, f), arcname=f)
    return buf.getvalue()


def read_text(path):
    return open(path, 'rb').read().decode('latin-1') if os.path.isfile(path) else ''


def rows_of(stdout):
    """FASTA rows / names of the standard output. Output that is not FASTA (CLUSTAL, PHYLIP, ...): lines before the first '>' header are not rows (the raw output is stored anyway), never an error."""
    out, names = [], []
    for l in stdout.decode('latin-1').splitlines():
        if l.startswith('>'):
            names.append(l[1:]); out.append('')
        elif l.strip() and out:
            out[-1] += l.strip()
    return out, names


esc = lambda s: s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')
results = open('out/results.tsv', 'a')                    # one row per set, written at once
progress = open('out/progress.log', 'a')
stopped = False
for line in open('inputs/sets.tsv'):
    cols = line.rstrip('\n').split('\t')
    k = cols[0]
    if not k.isdigit() or not (SETS_FROM <= int(k) <= SETS_TO):
        continue
    opt = cols[1] if len(cols) > 1 else ''
    print('set', k, opt, flush=True)
    progress.write('set %s %s\n' % (k, opt)); progress.flush()
    # defaults: whatever fails, the set is recorded (rc NA + the error) and the run goes on
    err, seqs_field = '', ''
    rc, out, names, stderr, tree, hat2, stdout = 'NA', [], [], b'', '', '', b''
    calls, direct, outfiles = pack_dir('/nonexistent'), {'rc': 'NA', 'stderr': '', 'hat2': '', 'hat3': ''}, {}
    try:
        seqs_field = cols[2]
        extra = json.loads(cols[3]) if len(cols) > 3 else {}      # files that the options name (new.fa, seed.fa, table.txt ...)
        seqs = seqs_field.split('|')
        with tempfile.TemporaryDirectory() as td:
            fa = os.path.join(td, 'in.fa')
            infa = extra.pop('in.fa', None) or ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))
            open(fa, 'w').write(infa)
            inputs = {'in.fa': infa}
            for fname, text in extra.items():
                open(os.path.join(td, fname), 'w').write(text)
                inputs[fname] = text
            p, d = run_logged(['mafft'] + opt.split() + [fa], td)
            rc, stderr, stdout = p.returncode, p.stderr, p.stdout
            tree = read_text(fa + '.tree')
            hat2 = read_text(fa + '.hat2')
            calls = pack_dir(d)
            # every file the run left in the working directory (outputs of the options: .tree, .map, hat2, pre, ...), not the inputs it was given
            for f in sorted(os.listdir(td)):
                fp = os.path.join(td, f)
                if os.path.isfile(fp) and os.path.getsize(fp) <= MAXFILE:
                    t = read_text(fp)
                    if inputs.get(f) != t:
                        outfiles[f] = t
            out, names = rows_of(stdout)
            # (2) dndfast7 directly
            wd = os.path.join(td, 'direct'); os.makedirs(wd)
            sw = ['-A'] if '--fastswpair' in opt else []
            exe = next((x for x in ('/usr/lib/mafft/lib/mafft/dndfast7', '/usr/lib/mafft/libexec/mafft/dndfast7', '/usr/libexec/mafft/dndfast7') if os.path.exists(x)), 'dndfast7')
            if int(k) >= DIRECT_FROM:
                q, d2 = run_logged([exe] + sw, wd, fa)
                direct = {'rc': q.returncode, 'stderr': q.stderr.decode('latin-1'), 'hat2': read_text(os.path.join(wd, 'hat2')), 'hat3': read_text(os.path.join(wd, 'hat3')), 'stdout': q.stdout.decode('latin-1')}
            else:
                direct = {'rc': 0, 'stderr': '', 'hat2': '', 'hat3': ''}      # not run again: the direct runs of these sets are in the archive of run 2
    except (Stop, KeyboardInterrupt):
        stopped = True
        break
    except Exception as e:
        err = '%s: %s' % (type(e).__name__, e)
        with open('out/sets/%s.error.txt' % k, 'w') as f:
            f.write(traceback.format_exc())
        print('  ERROR (recorded, the run goes on):', err, flush=True)
    for name, data in (('stderr.gz', stderr), ('stdout.gz', stdout)):
        with gzip.open('out/sets/%s.%s' % (k, name), 'wb') as f:
            f.write(data)
    with open('out/sets/%s.fastalog.tar.gz' % k, 'wb') as f:
        f.write(calls)
    with gzip.open('out/sets/%s.direct.json.gz' % k, 'wt', encoding='utf-8') as f:
        json.dump(direct, f)
    with gzip.open('out/sets/%s.files.json.gz' % k, 'wt', encoding='utf-8') as f:
        json.dump(outfiles, f)
    results.write('\t'.join([k, opt, seqs_field, str(rc), '|'.join(out), '|'.join(names), hashlib.sha256(stderr).hexdigest(), esc(tree), esc(hat2), esc(err)]) + '\n')
    results.flush()
print('interrupted' if stopped else 'done')
PY
PYPID=$!
trap 'kill -TERM $PYPID 2>/dev/null' TERM INT        # docker stop / Ctrl-C: the script finishes cleanly (the current set is dropped, results.tsv has every finished set)
wait $PYPID; wait $PYPID 2>/dev/null
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
