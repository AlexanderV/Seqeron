// Oracle harness for RobinHoodMap (B09 A5-10b3): the iteration order of robin_hood::unordered_map<uint32|uint64, uint32>
// of Mash 2.3's bundled robin_hood.h 3.9.1, after a deterministic operation stream (splitmix64).
// usage: robin_hood_order_harness < cases   (each line: bits mode seed n param)
#include <limits>
#include <cstring>
#include <stdexcept>
#include "robin_hood.h"
#include <cstdio>
#include <set>
#include <cstdint>
static uint64_t sm_state;
static uint64_t next64() {
    uint64_t z = (sm_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
template <typename K> static void run(int mode, uint64_t seed, long n, long param) {
    sm_state = seed;
    robin_hood::unordered_map<K, uint32_t> m;
    std::set<K> ordered;
    auto mask = [](uint64_t x) { return (K)x; };
    for (long i = 0; i < n; i++) {
        if (mode == 0) {                       // MinHashHeap::tryInsert with a bottom-param set; keys drawn from a domain of 4 * param
            K h = mask(next64() % (uint64_t)(param * 4));
            if (ordered.size() < (size_t)param || h < *ordered.rbegin()) {
                if (m.count(h) == 0) { m[h] = 1; ordered.insert(h); } else m[h] = m.at(h) + 1;
                if (m.size() > (size_t)param) { K top = *ordered.rbegin(); m.erase(top); ordered.erase(top); }
            }
        } else {                               // random inserts and erases over a key domain of param keys
            uint64_t r = next64();
            K h = mask((r >> 8) % (uint64_t)param * 0x9E3779B97F4A7C15ULL);
            if (r % 3 == 0) m.erase(h); else { if (m.count(h)) m[h] = m.at(h) + 1; else m[h] = (uint32_t)(r & 0xFF); }
        }
    }
    uint64_t fnv = 1469598103934665603ULL;
    for (auto it = m.begin(); it != m.end(); ++it) {
        fnv = (fnv ^ (uint64_t)it->first) * 1099511628211ULL;
        fnv = (fnv ^ (uint64_t)it->second) * 1099511628211ULL;
    }
    printf("%zu %llu\n", m.size(), (unsigned long long)fnv);
}
int main() {
    int bits, mode; unsigned long long seed; long n, param;
    while (scanf("%d %d %llu %ld %ld", &bits, &mode, &seed, &n, &param) == 5) {
        if (bits == 32) run<uint32_t>(mode, seed, n, param); else run<uint64_t>(mode, seed, n, param);
    }
}
