#include <stdlib.h>
#include <assert.h>

#include <ccfec/ff_basic.h>
#include <ccfec/ff.h>

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
    for (int i = 0; i < length; ++i)
    {
        dst[i] ^= gf2_8_multiply(coef, src[i]);
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
    const ff_unit inv_coef = gf2_8_invert(coef);
    for (int i = 0; i < length; ++i)
    {
        dst[i] = gf2_8_multiply(inv_coef, dst[i]);
    }
}

void init_ff()
{
}

void free_ff()
{
}

const uint8_t * get_gf2_8_multiplication_table()
{
    return NULL;
}

const uint8_t * get_gf2_8_inversion_table()
{
    return NULL;
}

int full_multiplication_table()
{
    return 0;
}
