#ifndef ARANS_ARANS_8_H
#define ARANS_ARANS_8_H

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
#define RATE_BITS 8                 //number of rate bits for adaption shift
#endif

#define CODE_BITS 21                //number of bits for coding
#define PROB_BITS 15                //number of bits for probability
#define CODE_NORM (1 << CODE_BITS)  //lower bound for normalization
#define PROB_SIZE (1 << PROB_BITS)  //total for probability factors
#define CHUNK_SIZE (1 << 13)        //number of bytes per chunk
#define ALPH_SIZE (1 << 8)          //number of characters in the alphabet
#define CDF_SIZE (ALPH_SIZE + 1)    //number of elements in cdf
#define ALIGN_SHIFT 7               //structure aligning

//implementation section
#ifdef __cplusplus
#define ALIGN32(A) A alignas(32)
#else
#define ALIGN32(A) alignas(32) A
#endif


//structs
struct Arans {
    uint32_t cdf[CDF_SIZE];
};

struct Range {
    uint32_t start;
    uint32_t width;
};

static uint16_t ALIGN32(updateMtx[ALPH_SIZE][ALPH_SIZE]);


// Encoder

//public function declarations
STORAGE_SPEC void aransInit(struct Arans *);

STORAGE_SPEC size_t aransEncode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

//internal function declarations
static size_t encChunk(uint32_t *, unsigned char *, size_t, const unsigned char *, size_t);

static int encPut(uint32_t *, unsigned char **, struct Range);

static size_t putOriginalSize(unsigned char *, size_t);

static int encFlush(const uint32_t *, unsigned char **, const unsigned char *);

static inline struct Range modRange(const uint32_t *, unsigned char);

static inline void modUpdate(uint32_t *, unsigned char);

//public functions
STORAGE_SPEC void aransInit(struct Arans *arans) {
    for (int i = 0; i < CDF_SIZE; ++i)
        arans->cdf[i] = i << (PROB_BITS - 8);

    for (int i = 0; i < ALPH_SIZE; ++i)
        for (int j = 0; j < ALPH_SIZE; ++j)
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

    uint32_t ALIGN32(cdf_align[CDF_SIZE + ALIGN_SHIFT]);
    uint32_t *cdf = &cdf_align[ALIGN_SHIFT];
    memcpy(cdf, arans->cdf, sizeof(arans->cdf));

    while (in_rem > CHUNK_SIZE) {
        if (!(ret = encChunk(cdf, out_cur, out_rem, in_cur, CHUNK_SIZE)))
            return 0;

        out_cur = &out_cur[ret];
        out_rem -= ret;
        in_cur = &in_cur[CHUNK_SIZE];
        in_rem -= CHUNK_SIZE;
    }

    if (!(ret = encChunk(cdf, out_cur, out_rem, in_cur, in_rem)))
        return 0;

    out_rem -= ret;
    memcpy(arans->cdf, cdf, sizeof(arans->cdf));
    return out_size - out_rem;
}

//internal functions
static size_t encChunk(uint32_t *cdf, unsigned char *out, size_t out_size, const unsigned char *in, size_t in_size) {
    unsigned char *ptr = &out[out_size];
    struct Range range[CHUNK_SIZE];
    uint32_t cod = CODE_NORM;

    for (size_t i = 0; i < in_size; ++i) {
        unsigned char n = in[i];
        range[i] = modRange(cdf, n);
        modUpdate(cdf, n);
    }

    for (size_t i = in_size; i > 0; --i)
        if (encPut(&cod, &ptr, range[i - 1]))
            return 0;

    if (encFlush(&cod, &ptr, out))
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

static int encFlush(const uint32_t *c,
                    unsigned char **pptr,
                    const unsigned char *lim) {
    if (*pptr < &lim[4])
        return 1;

    *--*pptr = *c;
    *--*pptr = *c >> 8;
    *--*pptr = *c >> 16;
    *--*pptr = *c >> 24;

    return 0;
}

static inline struct Range modRange(const uint32_t *cdf, unsigned char c) {
    return (struct Range) {cdf[c], cdf[c + 1] - cdf[c]};
}

static inline void modUpdate(uint32_t *cdf, unsigned char c) {
    for (int i = 1; i < ALPH_SIZE; ++i)
        cdf[i] = ((cdf[i] << RATE_BITS) + updateMtx[c][i - 1] - cdf[i]) >> RATE_BITS;
}

// Encoder


// Decoder

// public function declarations
STORAGE_SPEC size_t aransDecode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

STORAGE_SPEC size_t aransGetOutFileSize(unsigned char *);

// internal function declarations
static size_t decChunk(uint32_t *, unsigned char *, size_t, const unsigned char *, size_t);

static int decInit(uint32_t *, unsigned char **, const unsigned char *);

static int decPut(uint32_t *, unsigned char **, struct Range);

static inline uint32_t decGet(const uint32_t *);

static unsigned char modSymb(const uint32_t *, uint32_t);

// public functions
STORAGE_SPEC size_t aransDecode(struct Arans *arans, unsigned char *out, const size_t out_size, const unsigned char *in,
                                const size_t in_size) {
    size_t offset = 4;
    unsigned char *out_cur = out;
    const unsigned char *in_cur = &in[offset];
    size_t out_rem = out_size;
    size_t in_rem = in_size - offset;
    size_t ret;

    uint32_t ALIGN32(cdf_align[CDF_SIZE + ALIGN_SHIFT]);
    uint32_t *cdf = &cdf_align[ALIGN_SHIFT];
    memcpy(cdf, arans->cdf, sizeof(arans->cdf));

    while (out_rem > CHUNK_SIZE) {
        if (!(ret = decChunk(cdf, out_cur, CHUNK_SIZE, in_cur, in_rem)))
            return 0;

        in_cur = &in_cur[ret];
        in_rem -= ret;
        out_cur = &out_cur[CHUNK_SIZE];
        out_rem -= CHUNK_SIZE;
    }

    if (!(ret = decChunk(cdf, out_cur, out_rem, in_cur, in_rem)))
        return 0;

    in_rem -= ret;
    memcpy(arans->cdf, cdf, sizeof(arans->cdf));
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
static size_t
decChunk(uint32_t *cdf, unsigned char *out, const size_t out_size, const unsigned char *in, const size_t in_size) {
    unsigned char *ptr = (unsigned char *) in;
    uint32_t cod;

    if (decInit(&cod, &ptr, &in[in_size]))
        return 0;

    for (size_t i = 0; i < out_size; ++i) {
        unsigned char n = modSymb(cdf, decGet(&cod));
        struct Range range = modRange(cdf, n);
        modUpdate(cdf, n);

        if (decPut(&cod, &ptr, range))
            return 0;

        out[i] = n;
    }

    if (cod != CODE_NORM)
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

static inline uint32_t decGet(const uint32_t *c) {
    return *c & (PROB_SIZE - 1);
}

static unsigned char modSymb(
        const uint32_t *cdf,
        const uint32_t prb) {
    for (int i = 1; i < CDF_SIZE; ++i)
        if (prb < cdf[i])
            return i - 1;
}

// Decoder

#endif //ARANS_ARANS_8_H