/*
 * Speed tests for different sample generation algorithms
 */

#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>
#include <vector>
#include <getopt.h>
#include <string>
#include <algorithm>
#include "sobolmatrices.h"

// =============================================================================
// Utils
// =============================================================================

//#define RND_ZERO
#ifdef RND_ZERO
uint32_t rnd() { return 0; }
#else
std::random_device rd;                                                          // Will be used to obtain a seed for the random number engine
std::mt19937 rng(rd());
#define rnd() static_cast<uint32_t>(rng())
#endif

inline int randomBit() {
    static int counter = 0;
    static uint32_t bits;
    if (counter == 0) {
        bits = rnd();
        counter = 32;
    }
    counter--;
    int bit = bits & 1;
    bits >>= 1;
    return bit;
}

// =============================================================================
// Data structures
// =============================================================================

struct Vector {
    uint32_t x, y;
    inline uint32_t& operator[](uint32_t i) { return ((uint32_t *)this)[i]; };
    inline Vector operator^(const Vector   &v) { return {x ^ v.x, y ^ v.y}; };
    inline Vector operator>>(const int &i)     { return {x >> i , y >> i }; };
    inline void operator^= (const Vector &v) { x ^= v.x; y ^= v.y; };
};

typedef Vector Point;                                                           // An alias
typedef std::vector<Point> Points;
typedef void (*GenerationFunction) (Points &p);
typedef uint64_t (*InversionFunction) (
    int imgBitResolution, const Point &pixel, uint64_t pixelSampleIndex
);
typedef std::vector<uint32_t> TList;

// =============================================================================
// Generation Algorithms
// =============================================================================

// -----------------------------------------------------------------------------
// Random mt19937
// -----------------------------------------------------------------------------

inline void random(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        p[i] = {rnd(), rnd()};
    }
}

// -----------------------------------------------------------------------------
// Random, mainly for speed benchmarking
// -----------------------------------------------------------------------------

inline void randomFloat(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        float x = drand48();
        float y = drand48();
        p[i] = {*((uint32_t *)&x), *((uint32_t *)&y)};                          // A dummy assignment for speed benchmarking
    }
}

// -----------------------------------------------------------------------------
// van der Corput
// -----------------------------------------------------------------------------

inline uint32_t J(uint32_t i) {
    i = (i << 16) | (i >> 16);
    i = ((i & 0x00FF00FF) << 8) | ((i & 0xFF00FF00) >> 8);
    i = ((i & 0x0F0F0F0F) << 4) | ((i & 0xF0F0F0F0) >> 4);
    i = ((i & 0x33333333) << 2) | ((i & 0xCCCCCCCC) >> 2);
    i = ((i & 0x55555555) << 1) | ((i & 0xAAAAAAAA) >> 1);
    return i;
}

// -----------------------------------------------------------------------------
// SobolPBRT
// -----------------------------------------------------------------------------

inline uint32_t PascalPBRT(uint32_t i, uint32_t scramble = 0) {
    for (uint32_t v = 1 << 31; i != 0; i >>= 1, v ^= v >> 1) {
        if (i & 1) scramble ^= v;
    }
    return scramble;
}

inline void PBRT(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        p[i] = {J(i), PascalPBRT(i)};
    }
}

// -----------------------------------------------------------------------------
// I to X to Y
// -----------------------------------------------------------------------------

inline uint32_t P(uint32_t v) {                                                 // Input is a bit-reversed uint32_t
    v ^=  v               << 16;
    v ^= (v & 0x00FF00FF) <<  8;
    v ^= (v & 0x0F0F0F0F) <<  4;
    v ^= (v & 0x33333333) <<  2;
    v ^= (v & 0x55555555) <<  1;
    return v;
}

inline void I2X2Y(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = J(i);
        p[i] = {x, P(x)};
    }
}

inline Point sobol(uint32_t seqNo) {                                            // Will use this for individual points retrieval
    uint32_t x = J(seqNo);
    return {x, P(x)};
}

// -----------------------------------------------------------------------------
// I to (X, Y)
// -----------------------------------------------------------------------------

inline uint32_t PJ(uint32_t i) {                                                // Combine J and Sobol
    uint32_t j;
    j = (i & 0xFFFF0000); i = j ^ (j >> 16) ^ (i << 16);
    j = (i & 0xFF00FF00); i = j ^ (j >>  8) ^ ((i & 0x00FF00FF) << 8);
    j = (i & 0xF0F0F0F0); i = j ^ (j >>  4) ^ ((i & 0x0F0F0F0F) << 4);
    j = (i & 0xCCCCCCCC); i = j ^ (j >>  2) ^ ((i & 0x33333333) << 2);
    j = (i & 0xAAAAAAAA); i = j ^ (j >>  1) ^ ((i & 0x55555555) << 1);
    return i;
}

inline void I2XY(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        p[i] = {J(i), PJ(i)};
    }
}

// -----------------------------------------------------------------------------
// Morton: V to (X, Y)
// -----------------------------------------------------------------------------

inline uint32_t xi_u_from_VDC(uint32_t v) {                                     // Input is a bit-reversed uint32_t
    v ^= (v & 0x33333333) << 1;
    v ^= (v & 0x0F0F0F0F) << 2;
    v ^= (v & 0x00FF00FF) << 4;
    v ^= (v & 0x0000FFFF) << 8;
    return v;
}

inline uint32_t xiPlus_u_from_VDC(uint32_t v) {                                 // Input is a bit-reversed uint32_t
    v ^= (v & 0x55555555) << 1;
    v ^= (v & 0x33333333) << 1;
    v ^= (v & 0x0F0F0F0F) << 2;
    v ^= (v & 0x00FF00FF) << 4;
    v ^= (v & 0x0000FFFF) << 8;
    return v;
}

inline void N_V2XY(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = J(i);
        p[i] = {xi_u_from_VDC(x), xiPlus_u_from_VDC(x)};
    }
}

// -----------------------------------------------------------------------------
// Morton: I to X to Y
// -----------------------------------------------------------------------------

inline uint32_t xi_u_from_I(uint32_t i) {
    uint32_t j;
    i = ((i & 0x55555555) << 1) | ((i & 0xAAAAAAAA) >> 1);
    j = (i & 0xCCCCCCCC) >> 1; i = j ^ (j >> 1) ^ ((i & 0x33333333) <<  2);
    j = (i & 0xF0F0F0F0) >> 2; i = j ^ (j >> 2) ^ ((i & 0x0F0F0F0F) <<  4);
    j = (i & 0xFF00FF00) >> 4; i = j ^ (j >> 4) ^ ((i & 0x00FF00FF) <<  8);
    j = (i & 0xFFFF0000) >> 8; i = j ^ (j >> 8) ^ ( i               << 16);
    return i;
}

inline void N_I2X2Y(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = xi_u_from_I(i);
        p[i] = {x, P(x)};
    }
}

// -----------------------------------------------------------------------------
// Morton: I to (X,Y)
// -----------------------------------------------------------------------------

inline uint32_t xiPlus_u_from_I(uint32_t i) {                                       // Input is a bit-reversed uint32_t
    uint32_t j;
    j = (i & 0xAAAAAAAA)     ; i = j ^ (j >> 1) ^ ((i & 0x55555555) <<  1);
    j = (i & 0xCCCCCCCC) >> 1; i = j ^ (j >> 1) ^ ((i & 0x33333333) <<  2);
    j = (i & 0xF0F0F0F0) >> 2; i = j ^ (j >> 2) ^ ((i & 0x0F0F0F0F) <<  4);
    j = (i & 0xFF00FF00) >> 4; i = j ^ (j >> 4) ^ ((i & 0x00FF00FF) <<  8);
    j = (i & 0xFFFF0000) >> 8; i = j ^ (j >> 8) ^ ( i               << 16);
    return i;
}

inline void N_I2XY(Points &p) {
    for (uint32_t i = 0; i < p.size(); i++) {
        p[i] = {xi_u_from_I(i), xiPlus_u_from_I(i)};
    }
}

// -----------------------------------------------------------------------------
// Helmer implementation of Owen:
// -----------------------------------------------------------------------------

inline void helmer(Points &p) {
    static const uint32_t yXors[32] = {
        0x00000000, 0x00000001, 0x00000001, 0x00000007,
        0x00000001, 0x00000013, 0x00000015, 0x0000007f,
        0x00000001, 0x00000103, 0x00000105, 0x0000070f,
        0x00000111, 0x00001333, 0x00001555, 0x00007fff,
        0x00000001, 0x00010003, 0x00010005, 0x0007000f,
        0x00010011, 0x00130033, 0x00150055, 0x007f00ff,
        0x00010101, 0x01030303, 0x01050505, 0x070f0f0f,
        0x01111111, 0x13333333, 0x15555555, 0x7fffffff
    };
    int mMax = __builtin_ctz(p.size());
    p[0] = {0, 0};                                                              // Originally {rnd(), rnd()} for Owen scrambling
    for (int m = 0; m < mMax; m++) {                                            // Current resolution of net
        uint32_t N = 1 << m;                                                    // Current number of computed samples
        uint32_t msb = 0x80000000 >> m;                                         // Leading bit to subsequent octave
        for (uint32_t i = 0; i < N; i++) {                                      // Iterate through existing samples
            uint32_t X = i;                                                     // Sibling stratum of sample of p[N + i] along the x-axis
            uint32_t Y = i ^ yXors[m];                                          // Sibling stratum of sample of p[N + i] along the y-axis
            p[N + i].x = p[X].x ^ msb;                                          // Originally "p[X].x ^ ((rnd() | 0x80000000) >> m)" for Owen scramble
            p[N + i].y = p[Y].y ^ msb;                                          // Originally "p[Y].y ^ ((rnd() | 0x80000000) >> m)" for Owen scramble
        }
    }
}

// -----------------------------------------------------------------------------
// \xi
// -----------------------------------------------------------------------------

inline uint32_t ximul(uint32_t x) {
    return (x >> 1) ^ (x >> 2) ^ (x >> 4) ^ (x >> 8) ^ (x >> 16);
}

inline void xi(Points &p) {
    Vector u[4];
    Vector a = {rnd() | 0x80000000, rnd() | 0x80000000};
    Vector p0 = {rnd(), rnd()};
    Vector b = {ximul(a.x), ximul(a.y) ^ a.y};                                  // Oddly indexed columns; note adding a in y, equivalent to multiplying by 1 + 2xi
    u[0] = p0 ^ (p0 >> 1) ^ (p0 >> 16) ^ (p0 >> 17);                            // Obtained by solving {p0 = sum_{i=0}^{15} (u0 >> i)} for u0
    u[1] = u[0] ^ a;
    u[2] = u[0] ^ b;
    u[3] = u[1] ^ b;
    for (uint32_t i = 0; i < p.size(); i++) {                                   // Iterate through bit-reversed indices
        uint32_t seqNo = i;
        Point q = {0, 0};
        for (int depth = 0; depth < 16; depth++) {
            q ^= u[seqNo & 3] >> depth;                                             // Index the appropriate vector using the quadrant, seqNo & 3
            seqNo >>= 2;                                                            // Recursively subdivide the quadrant
        }
        p[i] = q;
    }
}

// =============================================================================
// GFaure: Scramble to optimize for m-bit point density
// =============================================================================

// -----------------------------------------------------------------------------
// JPJ, common
// -----------------------------------------------------------------------------

inline uint32_t JPJ(uint32_t v) {                                  //
    v ^=  v               >> 16;
    v ^= (v & 0xFF00FF00) >>  8;
    v ^= (v & 0xF0F0F0F0) >>  4;
    v ^= (v & 0xCCCCCCCC) >>  2;
    v ^= (v & 0xAAAAAAAA) >>  1;
    return v;
}

inline uint32_t JPJ(uint32_t v, int m) {
    return (JPJ(v >> (32 - m)) << (32 - m)) | (v & (0xFFFFFFFF >> m));          // Scramble leading m bits by
}

// -----------------------------------------------------------------------------
// Hammersley
// -----------------------------------------------------------------------------

inline void hammersley(Points &p) {
    int m = __builtin_ctz(p.size());
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = J(i);
        uint32_t y = P(x);
        p[i] = {JPJ(x, m), JPJ(y, m)};
    }
}

// -----------------------------------------------------------------------------
// LP
// -----------------------------------------------------------------------------

inline uint32_t G(uint32_t x, int m) {
    uint32_t v = JPJ(x >> (32 - m));
    v ^=  v >> 1;                                                               // Inverse of lower LP matrix
    return (v << (32 - m)) | (x & (0xFFFFFFFF >> m));
}

inline void LP(Points &p) {
    int m = __builtin_ctz(p.size());
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = J(i);
        uint32_t y = P(x);
        p[i] = {G(x, m), JPJ(y, m)};
    }
}

// -----------------------------------------------------------------------------
// mmd
//           Lx      Ly
//           1.....  1.....
//           .1....  11....
//           .11...  111...
//           ...1..  ...1..
//           ....1.  ..1.1.
//           .....1  .....1
// -----------------------------------------------------------------------------

inline uint32_t mmdX(uint32_t x, int m) {
    uint32_t v = JPJ(x >> (32 - m));
    int padding = (m - 6) >> 1;                                                 // Core matrix is 6x6, padded with identities on both sides
    v ^= (v & (0x10 << padding)) >> 1;                                          // <J U' J> has a single 1 below the fifth diagonal entry from bottom
    return (v << (32 - m)) | (x & (0xFFFFFFFF >> m));
}

inline uint32_t mmdY(uint32_t y, int m) {
    uint32_t v = JPJ(y >> (32 - m));
    int padding = (m - 6) >> 1;                                                 // Core matrix is 6x6, padded with identities on both sides
     v ^= ((v & (0x30 << padding)) >> 1) ^ ((v & (0x8 << padding)) >> 2);
    return (v << (32 - m)) | (y & (0xFFFFFFFF >> m));
}

inline void mmd(Points &p) {
    int m = __builtin_ctz(p.size());
    for (uint32_t i = 0; i < p.size(); i++) {
        uint32_t x = J(i);
        uint32_t y = P(x);
        p[i] = {mmdX(x, m), mmdY(y, m)};
    }
}

// =============================================================================
// Inversion Algorithms
// =============================================================================

// -----------------------------------------------------------------------------
// PBRT
// -----------------------------------------------------------------------------

inline uint64_t SobolIntervalToIndex(
    int m,                                                                      // Image plane pixel resolution per axis
    const Point &p,                                                             // Pixel coordinates
    uint64_t frame                                                              // Pixel sample index
) {
    if (m == 0) return 0;

    const uint32_t m2 = m << 1;
    uint64_t index = uint64_t(frame) << m2;

    uint64_t delta = 0;
    for (int c = 0; frame; frame >>= 1, ++c)
        if (frame & 1)  // Add flipped column m + c + 1.
            delta ^= VdCSobolMatrices[m - 1][c];

    // flipped b
    uint64_t b = (((uint64_t)((uint32_t)p.x) << m) | ((uint32_t)p.y)) ^ delta;

    for (int c = 0; b; b >>= 1, ++c)
        if (b & 1)  // Add column 2 * m - c.
            index ^= VdCSobolMatricesInv[m - 1][c];

    return index;
}

// -----------------------------------------------------------------------------
// I:P
// -----------------------------------------------------------------------------

inline uint32_t IP_inv (uint32_t v) {
    v ^= (v & 0x3333) <<  2;
    v ^= (v & 0x5555) <<  1;
    v ^= (v & 0x0F0F) <<  4;
    v ^= (v & 0x00FF) <<  8;
    v ^=  v           >> 16;
    return v;
}

inline uint32_t xi_u_inv(uint32_t v) {
    v ^= ((v & 0x0000FFFF) << 8) ^ ((v & 0x000000FF) << 16);
    v ^= ((v & 0x00FF00FF) << 4) ^ ((v & 0x000F000F) <<  8);
    v ^= ((v & 0x0F0F0F0F) << 2) ^ ((v & 0x03030303) <<  4);
    v ^= ((v & 0x33333333) << 1) ^ ((v & 0x11111111) <<  2);
    return v;
}

inline uint64_t Jxi_u_invJ(uint64_t v) {
    v ^= ((v & 0xFFFFFFFF00000000ll) >> 16) ^ ((v & 0xFFFF000000000000ll) >> 32);
    v ^= ((v & 0xFFFF0000FFFF0000ll) >>  8) ^ ((v & 0xFF000000FF000000ll) >> 16);
    v ^= ((v & 0xFF00FF00FF00FF00ll) >>  4) ^ ((v & 0xF000F000F000F000ll) >>  8);
    v ^= ((v & 0xF0F0F0F0F0F0F0F0ll) >>  2) ^ ((v & 0xC0C0C0C0C0C0C0C0ll) >>  4);
    v ^= ((v & 0xCCCCCCCCCCCCCCCCll) >>  1) ^ ((v & 0x8888888888888888ll) >>  2);
    return v;
}

inline uint64_t Jxi_uJ (uint64_t v) {
    v ^= (v & 0xCCCCCCCCCCCCCCCCllu) >>  1;
    v ^= (v & 0xF0F0F0F0F0F0F0F0llu) >>  2;
    v ^= (v & 0xFF00FF00FF00FF00llu) >>  4;
    v ^= (v & 0xFFFF0000FFFF0000llu) >>  8;
    v ^= (v & 0xFFFFFFFF00000000llu) >> 16;
    return v;
}


inline uint64_t invIP(int m, const Point &p, uint64_t sampleNo) {
    if (m == 0) return 0;
    int m2 = m << 1;
    uint32_t v = ((p.x << 16) | p.y) << (16 - m);                               // concat x:y as 16-bit fixed-point fractions
    sampleNo <<= m2;                                                            // Position sampleNo above pixel-indexing range of bits
    // -------------------------------------------------------------------------
    // Invert using [i:p]-inverse matrix.
    // This crude inversion gives the sample index of the pixel Origin
    // at 16-bit resolution, and will be refined subsequently
    // -------------------------------------------------------------------------
    v = IP_inv(v);
    // -------------------------------------------------------------------------
    // Translate to alpha-ordering so that all pixel samples have the same
    // suffix/prefix in normal/reversed order of bits.
    // -------------------------------------------------------------------------
    v = xi_u_inv(v);
    // -------------------------------------------------------------------------
    // Isolate the suffix/prefix bits for the pixel alpha index
    // -------------------------------------------------------------------------
    v = J(v);                                                          // Restore normal ordering of bits
    v &= ((1 << m2) - 1);                                                       // Isolate the 2*m-bits suffix
    // -------------------------------------------------------------------------
    // Adjust sampleNo to (I,P) ordering.
    // This is optional, ignoring it would reorder the samples within each
    // power-of-two prefix of the pixel samples.
    // We noticed that it does not change the rate, possibly thanks to
    // concurrent processing of the CPU.
    // -------------------------------------------------------------------------
    sampleNo = Jxi_u_invJ(sampleNo);
    sampleNo &= 0xFFFFFFFFFFFFFFFFllu << m2;                                    // Mask trailing bits.
    // -------------------------------------------------------------------------
    // Contatenate, bit-reverse, restore (I,P) ordering, and bit-reverse again
    // -------------------------------------------------------------------------
    return Jxi_uJ(v | sampleNo);                                                // Restore (I, P) ordering and
}

// -----------------------------------------------------------------------------
// Serial processing for evaluation
// -----------------------------------------------------------------------------

// double invert(
//     InversionFunction f, int imgRes, int log2spp,
//     std::vector<uint64_t> &list
// ) {
//     std::vector<int> list()
// }


// =============================================================================

const char *USAGE_MESSAGE = "Usage: %s [outputfile1 outputfile2 ..]\n"
"Options:\n"
"   -h                  Write this help message\n"
"   -a <algName>        sobol, sobolpbrt, xi, Helmer, random, drand48\n"
"   -r <m>              Benchmark generation rate using 2^m points\n"
"   -v                  Verify generation algorithms\n"
"   -R <imgRes>         Benchmark inversion\n"
"   -m <value>          log2 of number of produced points; default is 12.\n"
;

int main(int argc, char** argv) {
    int opt;                                                                    // For use by getopt, the command line options utility
    int m = 12;
    std::string fileName;
    int mRate = 0;
    int imgRes = 0;                                                             // Image resolution for inversion benchmarking
    std::string algName = "i2x2y";
    bool verifyGen = false;
    struct GenAlgorithm {
        std::string name;
        GenerationFunction function;
    };
    std::vector<GenAlgorithm> genAlgorithms{
        {"pbrt"      , PBRT         },
        {"i2x2y"     , I2X2Y        },
        {"i2xy"      , I2XY         },
        {"n-v2xy"    , N_V2XY       },
        {"n-i2x2y"   , N_I2X2Y      },
        {"n-i2xy"    , N_I2XY       },
        {"hammersley", hammersley   },
        {"lp"        , LP           },
        {"mmd"       , mmd         },
        {"helmer"    , helmer       },
        {"xi"        , xi           },
        {"random"    , random       },
        {"drand48"   , randomFloat  }
    };
    while ((opt = getopt(argc, argv, "m:va:r:hR:")) != -1) {
        switch (opt) {
            case 'r': mRate = atoi(optarg); break;
            case 'R': imgRes = atoi(optarg); break;
            case 'm': m = atoi(optarg); break;
            case 'v': verifyGen = true; break;
            case 'a': algName = optarg; break;
            case 'h':
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    transform(
        algName.begin(), algName.end(), algName.begin(), ::tolower
    );
    GenerationFunction netFunction = NULL;
    for (int i = 0; i < genAlgorithms.size(); i++) {
        if (algName == genAlgorithms[i].name) {
            netFunction = genAlgorithms[i].function;
            break;
        }
    }
    if (!netFunction) {
        fprintf(stderr, "Unknown algorithm name '%s'\n", algName);
        fprintf(stderr, USAGE_MESSAGE);
        exit(1);
    }
    uint32_t N = 1 << m;

    // =========================================================================
    // Validation and speed profiling of generation
    // =========================================================================
    if (mRate) {
        int m = mRate;
        uint32_t N = 1 << m;
        printf("Using (2^%d) = %d points\n", m, N);
        printf(
            "%10s %9s %4s %4s%s\n",
            "AlgName", "Time", "MSPP", "xRef",
            verifyGen ? " Same?" : ""
        );
        Points p(N);
        std::vector<uint32_t> X2Y;
        double refTime;
        for (int i = 0; i < genAlgorithms.size(); i++) {
            clock_t t0 = clock();
            genAlgorithms[i].function(p);
            clock_t t1 = clock();
            double totalTime = (double)(t1 - t0) / CLOCKS_PER_SEC;
            bool sameAsRef = true;                                              // Assume will end up true
            if (!i) {
                refTime = totalTime;                                            // Initialize reference time
                if (verifyGen) {
                    X2Y.resize(N);
                    for (int seqNo = 0; seqNo < N; seqNo++) {                       // Initialize reference mapping
                        X2Y[p[seqNo].x >> (32 - m)] = p[seqNo].y;
                    }
                }
            }
            else if (verifyGen) {
                for (int seqNo = 0; (seqNo < N) && sameAsRef; seqNo++) {        // Compare to reference mapping
                    uint32_t x = p[seqNo].x >> (32 - m);
                    if (x >= N) {
                        sameAsRef = false;
                        break;
                    }
                    sameAsRef = (X2Y[x] == p[seqNo].y);
                }
            }
            printf(
                "%10s %9.6f %4d %4.2f%5s\n",
                genAlgorithms[i].name.c_str(),
                totalTime, (int)(N/(totalTime * 1000000)),
                refTime / totalTime,
                verifyGen ? (sameAsRef ? " yes" : " no") : ""
            );
        }
    }

    // =========================================================================
    // Validation and speed profiling of inversion
    // =========================================================================

    if (imgRes) {
        int w = 1 << imgRes;
        printf("Using %d x %d image resolution\n", w, w);

        int x = rnd() & (w - 1);
        int y = rnd() & (w - 1);
        printf("x = %d, y = %d\n", x, y);
        for (int sampleNo = 0; sampleNo < 16; sampleNo++) {
            int iRef = SobolIntervalToIndex(imgRes, {x, y}, sampleNo);
            int iOurs = invIP(imgRes, {x, y}, sampleNo);
            Point pRef = sobol(iRef) >> (32 - imgRes);
            Point pOurs = sobol(iOurs) >> (32 - imgRes);
            printf(
                "%x: %8x --> (%4d, %4d) %s %8x --> (%4d, %4d)\n",
                sampleNo,
                iRef, pRef.x, pRef.y,
                iOurs == iRef ? "==" : "!=",
                iOurs, pOurs.x, pOurs.y
            );
        }

//         for (int y = w - 1; y >= 0; y--) {
//             for (int x = 0; x < w; x++) {
//                 Point p = {x, y};
//                 int i = SobolIntervalToIndex(imgRes, p, 0);
//                 int iOurs = invIP(imgRes, p, 0);
//                 printf("%3x/%3x%c", i, iOurs, x == w - 1? '\n' : ' ');
//             }
//         }

        clock_t t0 = clock();
        std::vector<uint32_t> seqNo(w * w);
        for (int yx2i = 0; yx2i < w * w; yx2i++) {
            Point p = {yx2i & (w-1), yx2i >> imgRes};
            int i = SobolIntervalToIndex(imgRes, p, 0);
            //int i = invIP(imgRes, p, 0);
            seqNo[yx2i] = i;
        }
        clock_t t1 = clock();
        double totalTime = (double)(t1 - t0) / CLOCKS_PER_SEC;
        printf(
            "Time = %9.6fs, Inversion rate = %4d MSPP\n",
            totalTime, (int)(w*w/(totalTime * 1000000))
        );
    }

    // =========================================================================
    // Save point sets
    // =========================================================================
    if (optind < argc) {                                                        // File names supplied?
        int n = argc - optind;
        double clocks = 0;
        Points p(N);
        std::vector<uint32_t> order(N);
        for (int i = 0; i < N; i++) { order[i] = i;  }                          // Natural ordering
        for (int fileNo = 0; fileNo < n; fileNo++) {
            clock_t t0 = clock();
            netFunction(p);
            clocks += clock() - t0;
            FILE *file = fopen(argv[optind + fileNo], "w");
            const double inv = 1.0 / (1ull << 32);
            fprintf(file, "%d\n", N);
            for (uint32_t seqNo = 0; seqNo < N; seqNo++) {
                uint32_t i = order[seqNo];
                fprintf(file, "%0.18f %0.18f\n", inv * p[i].x, inv * p[i].y);
            }
            fclose(file);
        }
        fprintf(
            stderr, "Average generation rate is %4dM points per second)\n",
            int(n * N * double(CLOCKS_PER_SEC) / (clocks * 1000000.0))
        );
    }
}
