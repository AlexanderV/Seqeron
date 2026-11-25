"""Line-by-line Python transcription of Musket 1.1 (Liu, Schmidt & Maskell 2013) error correction.

Source transcribed: musket-1.1/paracorrect.cpp (successor, predecessor, selectAllMutations,
paraCorrectCoreTwoSides, paraCorrectCore, ParaCorrect), voting.cpp (voteErroneousRead,
fixErroneousBase), paracorrect.h (CorrectRegionPing/PongComp), kmer.cpp (2-bit encoding,
operator< = memcmp over the packed bytes, twin), main.cpp/thread.cpp (canonical k-mer counting,
one count per occurrence). libstdc++ std::push_heap/pop_heap transcribed for region order.

Flags:
  canonical   - Musket counts canonical k-mers (min(km, twin) under memcmp order); False = forward only.
  fix_quirks  - True: use the true leftmost covering k-mer in stage 1 for reads shorter than 2k-1 and
                the true preceding k-mers in stage-2 predecessor(); False = Musket 1.1 code verbatim.
Non-ACGT symbols are kept as their own symbols (Musket encodes them as A): their candidate order is
A,C,G,T on the representative strand.
"""
CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
DEC = 'ACGT'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def rc(s):
    return ''.join(COMP.get(c, c) for c in reversed(s))


def musket_bytes(s):
    k = len(s)
    out = []
    for g in range((k + 3) // 4):
        v = 0
        for j in range(4):
            p = 4 * g + j
            if p < k:
                v |= CODE.get(s[p], 0) << (2 * j)
        out.append(v)
    return out


def musket_less(a, b):
    return musket_bytes(a) < musket_bytes(b)


class Spectrum:
    def __init__(self, reads, k, canonical):
        self.k = k
        self.canonical = canonical
        self.cnt = {}
        for r in reads:
            for i in range(len(r) - k + 1):
                x = r[i:i + k]
                self.cnt[x] = self.cnt.get(x, 0) + 1

    def mult(self, x):
        c = self.cnt.get(x, 0)
        if self.canonical:
            y = rc(x)
            if y != x:
                c += self.cnt.get(y, 0)
        return c

    def rep(self, km):
        """(rep, revcomp) as in Musket: rep = km if km < twin else twin (palindrome -> twin, revcomp)."""
        if not self.canonical:
            return km, False
        tw = rc(km)
        if musket_less(km, tw):
            return km, False
        return tw, True


def select_all_mutations(sp, w, km, index):
    k = sp.k
    rep, revcomp = sp.rep(km)
    base_index = k - 1 - index if revcomp else index
    ch = rep[base_index]
    if ch in CODE:
        o = CODE[ch]
        order = [(o + 1) & 3, (o + 2) & 3, (o + 3) & 3]
    else:
        order = [0, 1, 2, 3]
    res = []
    for b in order:
        alt = rep[:base_index] + DEC[b] + rep[base_index + 1:]
        m = sp.mult(alt)
        if m >= w:
            res.append(((b ^ 3) if revcomp else b, m))
    return res


def solids_of(sp, w, seq):
    k, L = sp.k, len(seq)
    sol = [False] * L
    for i in range(L - k + 1):
        if sp.mult(''.join(seq[i:i + k])) >= w:
            for j in range(i, i + k):
                sol[j] = True
    return sol


def two_sides(sp, w, seq, fix_quirks):
    k, L = sp.k, len(seq)
    if L < k:
        return True, []
    sol = solids_of(sp, w, seq)
    if all(sol):
        return True, []
    s = ''.join(seq)
    lkm = s[:k]  # Musket 1.1 lkm, advanced by forwardBase(seq[ipos]) for ipos >= k
    for ipos in range(L):
        if ipos >= k or ipos > L - k:
            lkm = lkm[1:] + s[ipos]
        if sol[ipos]:
            continue
        rstart = min(ipos, L - k)
        bases = select_all_mutations(sp, w, s[rstart:rstart + k], ipos - rstart)
        if fix_quirks:
            lstart = max(0, ipos - k + 1)
            lbases = select_all_mutations(sp, w, s[lstart:lstart + k], ipos - lstart)
        else:
            # Musket 1.1 verbatim: ipos < k probes lkm (= k-mer at 0) at index ipos; from ipos >= k on, and
            # in the tail loop, lkm has been advanced by seq[ipos] and is probed at index k-1. In the tail
            # loop Musket advances lkm for every ipos > L-k even when ipos < k, which splices the k-mer
            # (only when L < 2k-1).
            lbases = select_all_mutations(sp, w, lkm, ipos if (ipos < k and ipos <= L - k) else k - 1)
        corr = [(ipos, b) for b, _ in bases for lb, _ in lbases if b == lb]
        if len(corr) == 1:
            return False, corr
    return False, []


# ---- libstdc++ heap -------------------------------------------------------------------------
def _push_heap(a, hole, top, value, comp):
    parent = (hole - 1) // 2
    while hole > top and comp(a[parent], value):
        a[hole] = a[parent]
        hole = parent
        parent = (hole - 1) // 2
    a[hole] = value


def push_heap(a, comp):
    _push_heap(a, len(a) - 1, 0, a[-1], comp)


def _adjust_heap(a, hole, length, value, comp):
    top = hole
    second = hole
    while second < (length - 1) // 2:
        second = 2 * (second + 1)
        if comp(a[second], a[second - 1]):
            second -= 1
        a[hole] = a[second]
        hole = second
    if (length & 1) == 0 and second == (length - 2) // 2:
        second = 2 * (second + 1)
        a[hole] = a[second - 1]
        hole = second - 1
    _push_heap(a, hole, top, value, comp)


def pop_heap(a, comp):
    n = len(a)
    if n > 1:
        value = a[n - 1]
        a[n - 1] = a[0]
        _adjust_heap(a, 0, n - 1, value, comp)


def ping(a, b):
    return (a[1] - a[0] + 1) - (b[1] - b[0] + 1) < 0


def pong(a, b):
    return (a[1] - a[0] + 1) - (b[1] - b[0] + 1) > 0


def successor(sp, w, seq, start, km, dist):
    # seq[start:] is the remaining sequence; km is the k-mer ending at start-1
    k = sp.k
    L = len(seq) - start
    if L < k or dist <= 0:
        return True
    end = min(L - k, dist - 1)
    for i in range(end + 1):
        km = km[1:] + seq[start + i + k - 1]
        if sp.mult(km) < w:
            return False
    return True


def predecessor(sp, w, seq, seqlen, km, dist):
    if seqlen <= 0 or dist <= 0:
        return True
    startpos = max(0, seqlen - dist)
    for i in range(seqlen - 1, startpos - 1, -1):
        km = seq[i] + km[:-1]
        if sp.mult(km) < w:
            return False
    return True


def one_side(sp, w, seq, max_err, stride, pingpong, fix_quirks):
    k, L = sp.k, len(seq)
    if L < k:
        return False, []
    comp = ping if pingpong else pong
    sol = [False] * L
    regions = []
    solid_region = False
    left = right = -1
    for ipos in range(L - k + 1):
        if sp.mult(''.join(seq[ipos:ipos + k])) >= w:
            if not solid_region:
                solid_region = True
                left = right = ipos
            else:
                right += 1
            for j in range(ipos, ipos + k):
                sol[j] = True
        else:
            if left >= 0:
                regions.append((left, right))
                push_heap(regions, comp)
                left = right = -1
            solid_region = False
    if solid_region and left >= 0:
        regions.append((left, right))
        push_heap(regions, comp)
    if not regions:
        return True, []
    out = []
    while regions:
        region = regions[0]
        pop_heap(regions, comp)
        regions.pop()
        left, right = region
        km = ''.join(seq[right:right + k])
        last = -1
        ncorr = 0
        corr = []
        for ipos in range(right + 1, L - k + 1):
            t = ipos + k - 1
            if sol[t]:
                break
            km = km[1:] + seq[t]
            if sp.mult(km) < w:
                bases = select_all_mutations(sp, w, km, k - 1)
                done = False
                if len(bases) == 1:
                    b = bases[0][0]
                    corr.append((t, b))
                    km = km[:-1] + DEC[b]
                    done = True
                else:
                    best = (0, 0)
                    for b, m in bases:
                        if successor(sp, w, seq, ipos + 1, km[:-1] + DEC[b], stride):
                            if best[1] < m:
                                best = (b, m)
                    if best[1] > 0:
                        corr.append((t, best[0]))
                        km = km[:-1] + DEC[best[0]]
                        done = True
                if done:
                    if last < 0:
                        last = t
                    if t - last < k:
                        ncorr += 1
                        if ncorr > max_err:
                            for _ in range(ncorr):
                                corr.pop()
                            break
                    else:
                        last = t
                        ncorr = 0
                    continue
                break
        out.extend(corr)
        last = corr[0][0] if corr else -1
        ncorr = 0
        corr = []
        if left > 0:
            km = ''.join(seq[left:left + k])
            for ipos in range(left - 1, -1, -1):
                if sol[ipos]:
                    break
                km = seq[ipos] + km[:-1]
                if sp.mult(km) < w:
                    bases = select_all_mutations(sp, w, km, 0)
                    done = False
                    if len(bases) == 1:
                        b = bases[0][0]
                        corr.append((ipos, b))
                        km = DEC[b] + km[1:]
                        done = True
                    else:
                        best = (0, 0)
                        plen = ipos if fix_quirks else ipos - 1
                        for b, m in bases:
                            if predecessor(sp, w, seq, plen, DEC[b] + km[1:], stride):
                                if best[1] < m:
                                    best = (b, m)
                        if best[1] > 0:
                            corr.append((ipos, best[0]))
                            km = DEC[best[0]] + km[1:]
                            done = True
                    if done:
                        if last < 0:
                            last = ipos
                        if last - ipos < k:
                            ncorr += 1
                            if ncorr > max_err:
                                for _ in range(ncorr):
                                    corr.pop()
                                break
                        else:
                            last = ipos
                            ncorr = 0
                        continue
                    break
        out.extend(corr)
    return False, out


def vote(sp, w, seq):
    k, L = sp.k, len(seq)
    if L < k:
        return 0, None
    votes = [[0] * 4 for _ in range(L)]
    error_free = True
    for ipos in range(L - k + 1):
        km = ''.join(seq[ipos:ipos + k])
        rep, revcomp = sp.rep(km)
        if sp.mult(rep) >= w:
            continue
        error_free = False
        for off in range(k):
            bi = k - off - 1 if revcomp else off
            ch = rep[bi]
            if ch in CODE:
                o = CODE[ch]
                order = [(o + 1) & 3, (o + 2) & 3, (o + 3) & 3]
            else:
                order = [0, 1, 2, 3]
            for b in order:
                alt = rep[:bi] + DEC[b] + rep[bi + 1:]
                if sp.mult(alt) >= w:
                    votes[ipos + off][(b ^ 3) if revcomp else b] += 1
    if error_free:
        return 0, votes
    return max(max(v) for v in votes), votes


def fix_by_votes(votes, max_vote, min_vote, k, L):
    out = []
    if max_vote < min_vote or L < k:
        return out
    for pos in range(L):
        alt = -1
        for b in range(4):
            if votes[pos][b] == max_vote:
                alt = b if alt == -1 else -2
        if alt >= 0:
            out.append((pos, alt))
    return out


def para_correct(sp, w, read, max_iters=2, num_errors=4, fix_quirks=True):
    seq = list(read)
    for _ in range(max_iters):
        ef, corr = two_sides(sp, w, seq, fix_quirks)
        if ef:
            return ''.join(seq)
        if not corr:
            break
        for i, b in corr:
            seq[i] = DEC[b]
    for nerr in range(1, num_errors + 1):
        pp = 1
        for _ in range(max_iters):
            ef, corr = one_side(sp, w, seq, nerr, num_errors - nerr + 1, pp & 1, fix_quirks)
            pp += 1
            if ef:
                return ''.join(seq)
            if not corr:
                break
            for i, b in corr:
                seq[i] = DEC[b]
        mv, votes = vote(sp, w, seq)
        if mv == 0:
            return ''.join(seq)
        for i, b in fix_by_votes(votes, mv, 3, sp.k, len(seq)):
            seq[i] = DEC[b]
    return ''.join(seq)


def correct_all(reads, k, w, canonical=True, max_iters=2, num_errors=4, fix_quirks=True):
    reads = [r.upper() for r in reads]
    sp = Spectrum(reads, k, canonical)
    return [para_correct(sp, w, r, max_iters, num_errors, fix_quirks) for r in reads]
