#ifndef ARANS_ARANS_8_H
#define ARANS_ARANS_H

//process configuration
#ifdef ARANS_STATIC
#define STORAGE_SPEC static
#else
#define STORAGE_SPEC extern
#endif

//includes
#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>

//constants
#ifndef RATE_BITS
#define RATE_BITS 7                 //number of rate bits for adaption shift
#endif

#define CODE_BITS 21                //number of bits for coding
#define PROB_BITS 15                //number of bits for probability
#define CODE_NORM (1 << CODE_BITS)  //lower bound for normalization
#define PROB_SIZE (1 << PROB_BITS)  //total for probability factors
#define CHUNK_SIZE (1 << 13)        //number of bytes per chunk
#define ALIGN_SHIFT 7               //structure aligning

#define ALPH1_SIZE (1 << 2)         //number of characters in the alphabet 1
#define ALPH2_SIZE (1 << 6)         //number of characters in the alphabet 2
#define CDF1_SIZE (ALPH1_SIZE + 1)  //number of elements in cdf 1
#define CDF2_SIZE (ALPH2_SIZE + 1)  //number of elements in cdf 2

//implementation section
#ifdef __cplusplus
#define ALIGN16(A) A alignas(16)
#else
#define ALIGN16(A) alignas(16) A
#endif


//structs
struct Arans {
    uint16_t cdf1[CDF1_SIZE];
    uint16_t cdf2[ALPH1_SIZE][CDF2_SIZE];
};

struct Range {
    uint16_t start;
    uint16_t width;
};

static uint16_t ALIGN16(updateMtx[ALPH2_SIZE][ALPH2_SIZE]);


// Encoder

//public function declarations
STORAGE_SPEC void aransInit(struct Arans *);

STORAGE_SPEC size_t aransEncode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

//internal function declarations
static size_t encChunk(uint16_t *, uint16_t (*)[CDF2_SIZE], unsigned char *, size_t, const unsigned char *, size_t);

static int encPut(uint32_t *, unsigned char **, struct Range);

static size_t putOriginalSize(unsigned char *, size_t);

static int encFlush(const uint32_t *, unsigned char **, const unsigned char *);

static inline struct Range modRange(const uint16_t *, unsigned char);

static inline struct Range modSecondRange(const uint16_t (*)[CDF2_SIZE], unsigned char, unsigned char);

static inline void modUpdate(uint16_t *, unsigned char);

static inline void modSecondUpdate(uint16_t (*)[CDF2_SIZE], unsigned char, unsigned char);

//public functions
STORAGE_SPEC void aransInit(struct Arans *arans) {
    for (int i = 0; i < CDF1_SIZE; ++i)
        arans->cdf1[i] = i << (PROB_BITS - 8);

    for (int i = 0; i < ALPH1_SIZE; ++i)
        for (int j = 0; j < CDF2_SIZE; ++j)
            arans->cdf2[i][j] = j << (PROB_BITS - 8);

    for (int i = 0; i < ALPH2_SIZE; ++i)
        for (int j = 0; j < ALPH2_SIZE; ++j)
            updateMtx[i][j] = ((j + 1) + ((i < (j + 1)) ? PROB_SIZE - 1 : 0));
}

STORAGE_SPEC size_t
aransEncode(struct Arans *arans, unsigned char *out, size_t out_size, const unsigned char *in, size_t in_size) {
    size_t offset = putOriginalSize(out, in_size);
    unsigned char *out_cur = &out[offset];
    const unsigned char *in_cur = in;
    size_t out_rem = out_size - offset;
    size_t in_rem = in_size;
    size_t ret;

    uint16_t ALIGN16(cdf_align1[CDF1_SIZE + ALIGN_SHIFT]);
    uint16_t ALIGN16(cdf_align2[ALPH1_SIZE + ALIGN_SHIFT][CDF2_SIZE]);

    uint16_t *cdf1 = &cdf_align1[ALIGN_SHIFT];
    uint16_t (*cdf2)[CDF2_SIZE] = &cdf_align2[ALIGN_SHIFT];

    memcpy(cdf1, arans->cdf1, sizeof(arans->cdf1));
    memcpy(cdf2, arans->cdf2, sizeof(arans->cdf2));

    while (in_rem > CHUNK_SIZE) {
        if (!(ret = encChunk(cdf1, cdf2, out_cur, out_rem, in_cur, CHUNK_SIZE)))
            return 0;

        out_cur = &out_cur[ret];
        out_rem -= ret;
        in_cur = &in_cur[CHUNK_SIZE];
        in_rem -= CHUNK_SIZE;
    }

    if (!(ret = encChunk(cdf1, cdf2, out_cur, out_rem, in_cur, in_rem)))
        return 0;

    out_rem -= ret;
    memcpy(arans->cdf1, cdf1, sizeof(arans->cdf1));
    memcpy(arans->cdf2, cdf2, sizeof(arans->cdf2));
    return out_size - out_rem;
}

//internal functions
static size_t
encChunk(uint16_t *cdf1, uint16_t (*cdf2)[CDF2_SIZE], unsigned char *out, size_t out_size, const unsigned char *in,
         size_t in_size) {
    unsigned char *ptr = &out[out_size];
    struct Range range1[CHUNK_SIZE];
    struct Range range2[CHUNK_SIZE];
    uint32_t cod1 = CODE_NORM;
    uint32_t cod2 = CODE_NORM;

    for (size_t i = 0; i < in_size; ++i) {
        unsigned char n1 = in[i] >> 6;
        unsigned char n2 = in[i] & 0x3F;

        range1[i] = modRange(cdf1, n1);
        range2[i] = modSecondRange(cdf2, n1, n2);

        modUpdate(cdf1, n1);
        modSecondUpdate(cdf2, n1, n2);
    }

    for (size_t i = in_size; i > 0; --i) {
        if (encPut(&cod2, &ptr, range2[i - 1]))
            return 0;

        if (encPut(&cod1, &ptr, range1[i - 1]))
            return 0;
    }

    if (encFlush(&cod2, &ptr, out))
        return 0;

    if (encFlush(&cod1, &ptr, out))
        return 0;

    size_t size = &out[out_size] - ptr;
    memmove(out, ptr, size);
    return size;
}

static int encPut(uint32_t *c, unsigned char **pptr, struct Range range) {
    uint32_t x = *c;
    uint32_t x_max = range.width << (CODE_BITS - PROB_BITS + 8);

    if (x >= x_max) {
        unsigned char *ptr = *pptr;
        do {
            *--ptr = x;
            x >>= 8;
        } while (x >= x_max);
        *pptr = ptr;
    }

    *c = x + (PROB_SIZE - range.width) * (x / range.width) + range.start;
    return 0;
}

static size_t putOriginalSize(unsigned char *out, size_t in_size) {
    *out++ = in_size >> 24;
    *out++ = in_size >> 16;
    *out++ = in_size >> 8;
    *out++ = in_size;

    return 4;
}

static int encFlush(const uint32_t *c, unsigned char **pptr, const unsigned char *lim) {
    if (*pptr < &lim[4])
        return 1;

    *--*pptr = *c;
    *--*pptr = *c >> 8;
    *--*pptr = *c >> 16;
    *--*pptr = *c >> 24;

    return 0;
}

static inline struct Range modRange(const uint16_t *cdf, unsigned char c) {
    return (struct Range) {cdf[c], cdf[c + 1] - cdf[c]};
}

static inline struct Range modSecondRange(const uint16_t (*cdf)[CDF2_SIZE], unsigned char n1, unsigned char n2) {
    return (struct Range) {cdf[n1][n2], cdf[n1][n2 + 1] - cdf[n1][n2]};
}

static inline void modUpdate(uint16_t *cdf, unsigned char c) {
    for (int i = 1; i < CDF1_SIZE; ++i)
        cdf[i] = ((cdf[i] << RATE_BITS) + updateMtx[c][i - 1] - cdf[i]) >> RATE_BITS;
}

static inline void modSecondUpdate(uint16_t (*cdf)[CDF2_SIZE], unsigned char n1, unsigned char n2) {
    for (int i = 1; i < CDF2_SIZE; ++i)
        cdf[n1][i] = ((cdf[n1][i] << RATE_BITS) + updateMtx[n2][i - 1] - cdf[n1][i]) >> RATE_BITS;
}

// Encoder


// Decoder

// public function declarations
STORAGE_SPEC size_t aransDecode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

STORAGE_SPEC size_t aransGetOutFileSize(unsigned char *);

// internal function declarations
static size_t decChunk(uint16_t *, uint16_t (*)[CDF2_SIZE], unsigned char *, size_t, const unsigned char *, size_t);

static int decInit(uint32_t *, unsigned char **, const unsigned char *);

static int decPut(uint32_t *, unsigned char **, struct Range);

static inline uint16_t decGet(const uint32_t *);

static unsigned char modSymb(const uint16_t *, uint16_t);

static unsigned char modSecondSymb(const uint16_t (*)[CDF2_SIZE], unsigned char, uint16_t);

// public functions
STORAGE_SPEC size_t aransDecode(struct Arans *arans, unsigned char *out, const size_t out_size, const unsigned char *in,
                                const size_t in_size) {
    size_t offset = 4;
    unsigned char *out_cur = out;
    const unsigned char *in_cur = &in[offset];
    size_t out_rem = out_size;
    size_t in_rem = in_size - offset;
    size_t ret;

    uint16_t ALIGN16(cdf_align1[CDF1_SIZE + ALIGN_SHIFT]);
    uint16_t ALIGN16(cdf_align2[ALPH1_SIZE + ALIGN_SHIFT][CDF2_SIZE]);

    uint16_t *cdf1 = &cdf_align1[ALIGN_SHIFT];
    uint16_t (*cdf2)[CDF2_SIZE] = &cdf_align2[ALIGN_SHIFT];

    memcpy(cdf1, arans->cdf1, sizeof(arans->cdf1));
    memcpy(cdf2, arans->cdf2, sizeof(arans->cdf2));

    while (out_rem > CHUNK_SIZE) {
        if (!(ret = decChunk(cdf1, cdf2, out_cur, CHUNK_SIZE, in_cur, in_rem)))
            return 0;

        in_cur = &in_cur[ret];
        in_rem -= ret;
        out_cur = &out_cur[CHUNK_SIZE];
        out_rem -= CHUNK_SIZE;
    }

    if (!(ret = decChunk(cdf1, cdf2, out_cur, out_rem, in_cur, in_rem)))
        return 0;

    in_rem -= ret;
    memcpy(arans->cdf1, cdf1, sizeof(arans->cdf1));
    memcpy(arans->cdf2, cdf2, sizeof(arans->cdf2));
    return in_size - in_rem - offset;
}

STORAGE_SPEC size_t aransGetOutFileSize(unsigned char *in) {
    size_t size = *in++ << 24;
    size |= *in++ << 16;
    size |= *in++ << 8;
    size |= *in;

    return size;
}

// internal functions
static size_t decChunk(uint16_t *cdf1, uint16_t (*cdf2)[CDF2_SIZE], unsigned char *out, const size_t out_size,
                       const unsigned char *in,
                       const size_t in_size) {
    unsigned char *ptr = (unsigned char *) in;
    uint32_t cod1;
    uint32_t cod2;

    if (decInit(&cod1, &ptr, &in[in_size]))
        return 0;

    if (decInit(&cod2, &ptr, &in[in_size]))
        return 0;

    for (size_t i = 0; i < out_size; ++i) {
        unsigned char n1 = modSymb(cdf1, decGet(&cod1));
        unsigned char n2 = modSecondSymb(cdf2, n1, decGet(&cod2));

        struct Range range1 = modRange(cdf1, n1);
        struct Range range2 = modSecondRange(cdf2, n1, n2);

        modUpdate(cdf1, n1);
        modSecondUpdate(cdf2, n1, n2);

        if (decPut(&cod1, &ptr, range1))
            return 0;

        if (decPut(&cod2, &ptr, range2))
            return 0;

        out[i] = (n1 << 6) | n2;
    }

    if ((cod1 != CODE_NORM) || (cod2 != CODE_NORM))
        return 0;

    return ptr - in;
}

static int decInit(uint32_t *c, unsigned char **pptr, const unsigned char *lim) {
    if (*pptr > &lim[-4])
        return 1;

    *c = *(*pptr)++ << 24;
    *c |= *(*pptr)++ << 16;
    *c |= *(*pptr)++ << 8;
    *c |= *(*pptr)++;

    return 0;
}

static int decPut(uint32_t *c, unsigned char **pptr, const struct Range range) {
    uint32_t x = *c;
    x = range.width * (x >> PROB_BITS) + (x & (PROB_SIZE - 1)) - range.start;

    if (x < CODE_NORM) {
        unsigned char *ptr = *pptr;
        do {
            x = (x << 8) | *ptr++;
        } while (x < CODE_NORM);
        *pptr = ptr;
    }

    *c = x;
    return 0;
}

static inline uint16_t decGet(const uint32_t *c) {
    return *c & (PROB_SIZE - 1);
}

static unsigned char modSymb(const uint16_t *cdf, const uint16_t prb) {
    for (int i = 1; i < CDF1_SIZE; ++i)
        if (prb < cdf[i])
            return i - 1;
}

static unsigned char modSecondSymb(const uint16_t (*cdf)[CDF2_SIZE], unsigned char n1, const uint16_t prb) {
    for (int i = 1; i < CDF2_SIZE; ++i)
        if (prb < cdf[n1][i])
            return i - 1;
}

// Decoder

#endif //ARANS_ARANS_8_H