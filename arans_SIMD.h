#ifndef ARANS_ARANS_H
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
#include <immintrin.h>

//constants
#ifndef RATE_BITS
#define RATE_BITS 8 //number of rate bits for adaption shift
#endif

#define CODE_BITS 20                //number of bits for coding
#define PROB_BITS 15                //number of bits for probability
#define CODE_NORM (1 << CODE_BITS)  //lower bound for normalization
#define PROB_SIZE (1 << PROB_BITS)  //total for probability factors
#define CHUNK_SIZE (1 << 13)        //number of bytes per chunk
#define ALPH_SIZE (1 << 8)          //number of characters in the alphabet
#define CDF_SIZE (ALPH_SIZE + 1)    //number of characters in the alphabet
#define ALIGN_SHIFT 17              //structure aligning

//implementation section
#ifdef __cplusplus
#define ALIGN_ALPH_SIZE(A) A alignas(ALPH_SIZE)
#else
#define ALIGN_ALPH_SIZE(A) alignas(ALPH_SIZE) A
#endif


//structs
struct Arans {
    uint16_t cdf[CDF_SIZE];
};

struct range {
    uint16_t start;
    uint16_t width;
};

static uint16_t ALIGN_ALPH_SIZE(updateMtx[ALPH_SIZE][ALPH_SIZE]);


// Encoder

//public function declarations
STORAGE_SPEC void aransInit(struct Arans *);

STORAGE_SPEC size_t aransEncode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

//internal function declarations
static size_t encChunk(uint16_t *, unsigned char *, size_t, const unsigned char *, size_t);

static int encPut(uint32_t *, unsigned char **, struct range);

static size_t putOriginalSize(unsigned char *, size_t);

static int encFlush(const uint32_t *, unsigned char **, const unsigned char *);

static struct range modRange(const uint16_t *, unsigned char);

static void modUpdate(uint16_t *, unsigned char);

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
    uint16_t ALIGN_ALPH_SIZE(cdf_align[CDF_SIZE + ALIGN_SHIFT]);
    uint16_t *cdf = &cdf_align[ALIGN_SHIFT];
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
static size_t
encChunk(uint16_t *cdf, unsigned char *out, size_t out_size, const unsigned char *in, size_t in_size) {
    struct range rng_buff[CHUNK_SIZE];
    unsigned char *ptr = &out[out_size];
    uint32_t cod = CODE_NORM;

    for (size_t i = 0; i < in_size; ++i) {
        unsigned char n = in[i];
        rng_buff[i] = modRange(cdf, n);
        modUpdate(cdf, n);
    }

    for (size_t i = in_size; i > 0; --i)
        if (encPut(&cod, &ptr, rng_buff[i - 1]))
            return 0;

    if (encFlush(&cod, &ptr, out))
        return 0;

    size_t size = &out[out_size] - ptr;
    memmove(out, ptr, size);
    return size;
}

static int encPut(uint32_t *c, unsigned char **pptr, struct range rng) {
    uint32_t x = *c;
    uint32_t x_max = rng.width << (CODE_BITS - PROB_BITS + 8);

    if (x >= x_max) {
        unsigned char *ptr = *pptr;
        do {
            *--ptr = x;
            x >>= 8;
        } while (x >= x_max);
        *pptr = ptr;
    }

    *c = x + (PROB_SIZE - rng.width) * (x / rng.width) + rng.start;
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

static struct range modRange(const uint16_t *cdf, unsigned char c) {
    return (struct range) {cdf[c], cdf[c + 1] - cdf[c]};
}

static void modUpdate(uint16_t *cdf, unsigned char c) {
#ifndef __AVX2__
    for (int i = 1; i < ALPH_SIZE; ++i)
        cdf[i] = ((cdf[i] << RATE_BITS) + updateMtx[c][i - 1] - cdf[i]) >> RATE_BITS;
#else
    __m256i cdf1 = _mm256_load_si256((__m256i*)&cdf[1]);
    __m256i cdf17 = _mm256_load_si256((__m256i*)&cdf[17]);
    __m256i cdf33 = _mm256_load_si256((__m256i*)&cdf[33]);
    __m256i cdf49 = _mm256_load_si256((__m256i*)&cdf[49]);
    __m256i cdf65 = _mm256_load_si256((__m256i*)&cdf[65]);
    __m256i cdf81 = _mm256_load_si256((__m256i*)&cdf[81]);
    __m256i cdf97 = _mm256_load_si256((__m256i*)&cdf[97]);
    __m256i cdf113 = _mm256_load_si256((__m256i*)&cdf[113]);
    __m256i cdf129 = _mm256_load_si256((__m256i*)&cdf[129]);
    __m256i cdf145 = _mm256_load_si256((__m256i*)&cdf[145]);
    __m256i cdf161 = _mm256_load_si256((__m256i*)&cdf[161]);
    __m256i cdf177 = _mm256_load_si256((__m256i*)&cdf[177]);
    __m256i cdf193 = _mm256_load_si256((__m256i*)&cdf[193]);
    __m256i cdf209 = _mm256_load_si256((__m256i*)&cdf[209]);
    __m256i cdf225 = _mm256_load_si256((__m256i*)&cdf[225]);
    __m256i cdf241 = _mm256_load_si256((__m256i*)&cdf[241]);

    *(__m256i*)&cdf[1] = _mm256_add_epi16(cdf1, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][0], cdf1), RATE_BITS));
    *(__m256i*)&cdf[17] = _mm256_add_epi16(cdf17, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][16], cdf17), RATE_BITS));
    *(__m256i*)&cdf[33] = _mm256_add_epi16(cdf33, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][32], cdf33), RATE_BITS));
    *(__m256i*)&cdf[49] = _mm256_add_epi16(cdf49, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][48], cdf49), RATE_BITS));
    *(__m256i*)&cdf[65] = _mm256_add_epi16(cdf65, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][64], cdf65), RATE_BITS));
    *(__m256i*)&cdf[81] = _mm256_add_epi16(cdf81, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][80], cdf81), RATE_BITS));
    *(__m256i*)&cdf[97] = _mm256_add_epi16(cdf97, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][96], cdf97), RATE_BITS));
    *(__m256i*)&cdf[113] = _mm256_add_epi16(cdf113, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][112], cdf113), RATE_BITS));
    *(__m256i*)&cdf[129] = _mm256_add_epi16(cdf129, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][128], cdf129), RATE_BITS));
    *(__m256i*)&cdf[145] = _mm256_add_epi16(cdf145, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][144], cdf145), RATE_BITS));
    *(__m256i*)&cdf[161] = _mm256_add_epi16(cdf161, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][160], cdf161), RATE_BITS));
    *(__m256i*)&cdf[177] = _mm256_add_epi16(cdf177, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][176], cdf177), RATE_BITS));
    *(__m256i*)&cdf[193] = _mm256_add_epi16(cdf193, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][192], cdf193), RATE_BITS));
    *(__m256i*)&cdf[209] = _mm256_add_epi16(cdf209, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][208], cdf209), RATE_BITS));
    *(__m256i*)&cdf[225] = _mm256_add_epi16(cdf225, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][224], cdf225), RATE_BITS));
    *(__m256i*)&cdf[241] = _mm256_add_epi16(cdf241, _mm256_srai_epi16(_mm256_sub_epi16(*(__m256i*)&updateMtx[c][240], cdf241), RATE_BITS));
#endif
}

// Encoder


// Decoder

// public function declarations
STORAGE_SPEC size_t aransDecode(struct Arans *, unsigned char *, size_t, const unsigned char *, size_t);

STORAGE_SPEC size_t aransGetOutFileSize(unsigned char *);

// internal function declarations
static size_t decChunk(uint16_t *, unsigned char *, size_t, const unsigned char *, size_t);

static int decInit(uint32_t *, unsigned char **, const unsigned char *);

static int decPut(uint32_t *, unsigned char **, struct range);

static uint16_t decGet(const uint32_t *);

static unsigned char modSymb(const uint16_t *, uint16_t);

// public functions
STORAGE_SPEC size_t
aransDecode(struct Arans *arans, unsigned char *out, const size_t out_size, const unsigned char *in,
            const size_t in_size) {
    size_t offset = 4;
    unsigned char *out_cur = out;
    const unsigned char *in_cur = &in[offset];
    size_t out_rem = out_size;
    size_t in_rem = in_size - offset;
    size_t ret;
    uint16_t ALIGN_ALPH_SIZE(cdf_align[CDF_SIZE + ALIGN_SHIFT]);
    uint16_t *cdf = &cdf_align[ALIGN_SHIFT];
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
    size_t size = 0;
    size = *in++ << 24;
    size |= *in++ << 16;
    size |= *in++ << 8;
    size |= *in;

    return size;
}

// internal functions
static size_t
decChunk(uint16_t *cdf, unsigned char *out, const size_t out_size, const unsigned char *in, const size_t in_size) {
    unsigned char *ptr = (unsigned char *) in;
    uint32_t cod;

    if (decInit(&cod, &ptr, &in[in_size]))
        return 0;

    for (size_t i = 0; i < out_size; ++i) {
        unsigned char n = modSymb(cdf, decGet(&cod));
        struct range rng = modRange(cdf, n);
        modUpdate(cdf, n);

        if (decPut(&cod, &ptr, rng))
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

static int decPut(uint32_t *c, unsigned char **pptr, const struct range rng) {
    uint32_t x = *c;
    x = rng.width * (x >> PROB_BITS) + (x & (PROB_SIZE - 1)) - rng.start;

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

static uint16_t decGet(const uint32_t *c) {
    return *c & (PROB_SIZE - 1);
}

static unsigned char modSymb(const uint16_t *cdf, const uint16_t prb) {
    for (int i = 1; i < CDF_SIZE; ++i) {
        if (prb < cdf[i])
            return i - 1;
    }
}

// Decoder

#endif //ARANS_ARANS_H