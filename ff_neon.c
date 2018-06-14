#include <stdlib.h>
#include <assert.h>
#include <arm_neon.h>

#include <ccfec/ff_basic.h>
#include <ccfec/ff.h>

static const uint8_t* restrict gf2_8_inversion_table = NULL;
static const uint8_t* restrict gf2_8_inversion_table_raw = NULL;
static const uint8_t* restrict gf2_8_multiplication_table = NULL;
static const uint8_t* restrict gf2_8_multiplication_table_raw = NULL;

int SIMD_size()
{
    return 16;
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
    uint8_t* const p = (uint8_t*)malloc(256 * 16 * 2 + SIMD_size());
    uint8_t* const multiplication_table = ceil_to_grid_p(p);
    gf2_8_multiplication_table_raw = p;

    for (int i = 0; i < 256; ++i)
    {
        const int offset = i * 16;
        for (int j = 0; j < 16; ++j)
        {
            uint8_t const low = gf2_8_multiply(i, j);
            uint8_t const high = gf2_8_multiply(i, j << 4);

            multiplication_table[offset + j] = low;
            multiplication_table[offset + j + 4096] = high;
        }
    }
    return multiplication_table;
}

void encode_unit_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length)
{
    uint8_t* restrict dst_m = dst; // Mutable version of dst
    const uint8_t* restrict src_m = src; // Mutable version of src
    for (int i = 0; i < length; ++i, dst += SIMD_size(), src_m += SIMD_size())
    {
        // Load coefficients from dst_m (must be unaligned) and src_m (should be aligned)
        const uint8x16_t c_src = vld1q_u8(src_m);
        uint8x16_t c_dst = vld1q_u8(dst);
        c_dst = veorq_u8(c_dst, c_src);
        vst1q_u8(dst, c_dst);
    }
}

void encode_rich_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length, const ff_unit coef)
{
    // TODO: Implement
    abort();
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
    // TODO: Implement
    abort();
}

void init_ff()
{
    gf2_8_inversion_table = gf2_8_calculate_inversion_table();
    gf2_8_multiplication_table = gf2_8_calculate_multiplication_table();
}

void free_ff()
{
    free((void *)gf2_8_inversion_table_raw);
    free((void *)gf2_8_multiplication_table_raw);
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
