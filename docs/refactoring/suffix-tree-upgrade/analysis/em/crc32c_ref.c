/* Outside the shipped code (work package 1b.7e, gate G2): CRC-32C (Castagnoli) of a byte range of a file, a plain table implementation written from the definition
   (polynomial 0x82F63B78 reflected, initial value and final xor 0xFFFFFFFF), independent of the C# code. usage: crc32c_ref <file> <offset> <length>  -> prints the CRC as 8 hex digits */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: crc32c_ref <file> <offset> <length>\n"); return 2; }
    uint32_t table[256];
    for (uint32_t i = 0; i < 256; i++) { uint32_t c = i; for (int k = 0; k < 8; k++) c = (c & 1) ? (c >> 1) ^ 0x82F63B78u : c >> 1; table[i] = c; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("open"); return 1; }
    long long offset = atoll(argv[2]), length = atoll(argv[3]);
    if (fseeko(f, offset, SEEK_SET) != 0) { perror("seek"); return 1; }
    static unsigned char buf[1 << 20];
    uint32_t crc = 0xFFFFFFFFu;
    while (length > 0) {
        size_t want = length < (long long)sizeof buf ? (size_t)length : sizeof buf;
        size_t got = fread(buf, 1, want, f);
        if (got == 0) { fprintf(stderr, "short read\n"); return 1; }
        for (size_t i = 0; i < got; i++) crc = table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
        length -= (long long)got;
    }
    printf("%08x\n", crc ^ 0xFFFFFFFFu);
    return 0;
}
