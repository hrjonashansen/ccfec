#include <stdlib.h>
#include <assert.h>

#include <ccfec/ff_basic.h>
#include <ccfec/ff.h>

static const uint8_t* restrict gf2_8_inversion_table = NULL;
static const uint8_t* restrict gf2_8_multiplication_table = NULL;

int SIMD_size()
{
    return 0;
}

uint8_t* ceil_to_grid_p(uint8_t* const arg_p)
{
    return arg_p;
}

ff_unit* get_coef_p(uint8_t* arg)
{
    return arg;
}

int ff_row_size(const int symbol_size)
{
    return symbol_size;
}

int ff_math_size(const int symbol_size)
{
    return symbol_size;
}

int payload_size(const int k, const int symbol_size)
{
    return k * (int)sizeof(ff_unit) + ff_row_size(symbol_size + sizeof(int));
}

const uint8_t* gf2_8_calculate_inversion_table()
{
    const int max = 256;
    uint8_t* const table = (uint8_t*)malloc(max);

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
    uint8_t* const multiplication_table = (uint8_t*)malloc(256 * 256);

    multiplication_table[0] = 0;
    for (uint8_t x = 1; x != 0; ++x)
    {
        const int offset = x * 256;
        multiplication_table[x] = 0;
        multiplication_table[offset] = 0;
        for (uint8_t y = 1; y != 0; ++y)
        {
            multiplication_table[offset + y] = gf2_8_multiply(x, y);
        }
    }
    return multiplication_table;
}

void encode_unit_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length)
{
    for (int i = 0; i < length; ++i)
    {
        dst[i] ^= src[i];
    }
}

void encode_rich_symbol(uint8_t* restrict const dst, const uint8_t* restrict const src, const int length, const ff_unit coef)
{
    assert(coef != 0);
    const ff_unit* restrict const offset = gf2_8_multiplication_table + (coef << 8);
    for (int i = 0; i < length; ++i)
    {
        dst[i] ^= offset[src[i]];
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

void normalise_symbol(uint8_t* restrict const dst, const int length, const ff_unit coef)
{
    assert(coef != 0);
    const ff_unit* const offset = gf2_8_multiplication_table + (gf2_8_inversion_table[coef] << 8);
    for (int i = 0; i < length; ++i)
    {
        dst[i] = offset[dst[i]];
    }
}

void init_ff()
{
    gf2_8_inversion_table = gf2_8_calculate_inversion_table();
    gf2_8_multiplication_table = gf2_8_calculate_multiplication_table();
}

void free_ff()
{
    free((void *)gf2_8_inversion_table);
    free((void *)gf2_8_multiplication_table);
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
    return 1;
}
