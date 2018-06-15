#include <stdlib.h>
#include <assert.h>
#include <x86intrin.h>

#include <ccfec/ff_basic.h>
#include <ccfec/ff.h>

static const uint8_t* restrict gf2_8_inversion_table = NULL;
static const uint8_t* restrict gf2_8_inversion_table_raw = NULL;
static const uint8_t* restrict gf2_8_multiplication_table = NULL;
static const uint8_t* restrict gf2_8_multiplication_table_raw = NULL;

int SIMD_size()
{
    return 32;
}

int ceil_to_grid(const int arg)
{
    const int grid = SIMD_size();
    const int modulo = arg & (grid - 1);
    const int res = (modulo == 0 ? arg : arg + grid - modulo);
    return res;
}

uint8_t* ceil_to_grid_p(uint8_t* const arg_p)
{
    const intptr_t arg = (intptr_t)arg_p;
    const int grid = SIMD_size();
    const intptr_t modulo = arg & (grid - 1);
    const intptr_t res = (modulo == 0 ? arg : arg + grid - modulo);
    uint8_t* res_p = (uint8_t*)res;
    return res_p;
}

ff_unit* get_coef_p(uint8_t* arg)
{
    return arg;
}

int ff_row_size(const int symbol_size)
{
    return ceil_to_grid(symbol_size);
}

int ff_math_size(const int symbol_size)
{
    return ceil_to_grid(symbol_size) / SIMD_size();
}

int payload_size(const int k, const int symbol_size)
{
    return ceil_to_grid(k * sizeof(ff_unit) + ff_row_size(symbol_size + sizeof(int)));
}

const uint8_t* gf2_8_calculate_inversion_table()
{
    const int max = 256;
    uint8_t* const p = (uint8_t*)malloc(max + SIMD_size());
    uint8_t* const table = ceil_to_grid_p(p);
    gf2_8_inversion_table_raw = p;

    table[1] = 1;
    for (int i = 2; i < max; ++i)
    {
        table[i] = gf2_8_invert(i);
        assert(table[i] != 0);
        assert(table[i] != 1);
    }
    return table;
}

const uint8_t* gf2_8_calculate_multiplication_table()
{
    uint8_t* const p = (uint8_t*)malloc(256 * 16 * 2 * 2 + SIMD_size());
    uint8_t* const multiplication_table = ceil_to_grid_p(p);
    gf2_8_multiplication_table_raw = p;

    for (int i = 0; i < 256; ++i)
    {
        const int offset = i * 16 * 2;
        for (int j = 0; j < 16; ++j)
        {
            uint8_t const low = gf2_8_multiply(i, j);
            uint8_t const high = gf2_8_multiply(i, j << 4);

            multiplication_table[offset + j] = low;
            multiplication_table[offset + j + 16] = low;
            multiplication_table[offset + j + 8192] = high;
            multiplication_table[offset + j + 8192 + 16] = high;
        }
    }
    return multiplication_table;
}

void encode_unit_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length)
{
    // Create the speed up pointers
    const __m256i* src_p = (const __m256i*)src;
    __m256i* dst_p = (__m256i*)dst;
    for (int i = 0; i < length; ++i, ++dst_p, ++src_p)
    {
        // Load coefficients from dst (unaligned) and src (aligned)
        const __m256i c_src = _mm256_load_si256(src_p);
        __m256i c_dst = _mm256_loadu_si256(dst_p);
        // Do XOR
        c_dst = _mm256_xor_si256(c_dst, c_src);
        // Store result (unaligned)
        _mm256_storeu_si256(dst_p, c_dst);
    }
}

void encode_rich_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length, const ff_unit coef)
{
    const __m256i tableLow  = _mm256_load_si256((const __m256i*)(gf2_8_multiplication_table + coef * 16 * 2));
    const __m256i tableHigh = _mm256_load_si256((const __m256i*)(gf2_8_multiplication_table + coef * 16 * 2 + 8192));
        const __m256i mask = _mm256_set1_epi8(0x0f);
    const __m256i* src_p = (const __m256i*)src;
    __m256i* dst_p = (__m256i*)dst;
    for (int i = 0; i < length; ++i, ++dst_p, ++src_p)
    {
        __m256i c_src = _mm256_load_si256(src_p);
        __m256i low = _mm256_and_si256(c_src, mask);
        low = _mm256_shuffle_epi8(tableLow, low);

        __m256i high = _mm256_srli_epi64(c_src, 4);
        high = _mm256_and_si256(high, mask);
        high = _mm256_shuffle_epi8(tableHigh, high);
        c_src = _mm256_xor_si256(low, high);

        __m256i c_dst = _mm256_loadu_si256(dst_p);
        c_dst = _mm256_xor_si256(c_dst, c_src);

        _mm256_storeu_si256(dst_p, c_dst);
    }
}

void encode_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length, const ff_unit coef)
{
    if (coef == 1)
    {
        encode_unit_symbol(dst, src, ff_math_size(length));
    }
    else
    {
        encode_rich_symbol(dst, src, ff_math_size(length), coef);
    }
}

void normalise_symbol(uint8_t* restrict const dst, const int length_in, const ff_unit coef)
{
    const int length = ff_math_size(length_in);
    const uint8_t inv_coef = gf2_8_inversion_table[coef];
    const __m256i tableLow  = _mm256_load_si256((const __m256i*)(gf2_8_multiplication_table + inv_coef * 16 * 2));
    const __m256i tableHigh = _mm256_load_si256((const __m256i*)(gf2_8_multiplication_table + inv_coef * 16 * 2 + 8192));
        const __m256i mask = _mm256_set1_epi8(0x0f);
    __m256i* dst_p = (__m256i*)dst;
    for (int i = 0; i < length; ++i, ++dst_p)
    {
        __m256i full = _mm256_load_si256(dst_p);
        __m256i low = _mm256_and_si256(full, mask);
        low = _mm256_shuffle_epi8(tableLow, low);

        __m256i high = _mm256_srli_epi64(full, 4);
        high = _mm256_and_si256(high, mask);
        high = _mm256_shuffle_epi8(tableHigh, high);
        full = _mm256_xor_si256(low, high);
        _mm256_store_si256(dst_p, full);
    }
}

void init_ff()
{
    gf2_8_inversion_table = gf2_8_calculate_inversion_table();
    gf2_8_multiplication_table = gf2_8_calculate_multiplication_table();
}

void free_ff()
{
    free((void*)gf2_8_inversion_table_raw);
    free((void*)gf2_8_multiplication_table_raw);
}

const uint8_t * get_gf2_8_multiplication_table()
{
    return gf2_8_multiplication_table;
}

const uint8_t * get_gf2_8_inversion_table()
{
    return gf2_8_inversion_table;
}

int full_multiplication_table()
{
    return 0;
}
