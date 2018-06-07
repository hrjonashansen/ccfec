#pragma once

#include <stdint.h>

#ifdef __cplusplus
#define RESTRICT
#else
#define RESTRICT restrict
#endif

typedef uint8_t ff_unit;

void encode_symbol(const uint8_t* RESTRICT const s, uint8_t* RESTRICT const mem, const int length, const ff_unit coef);
void normalise_symbol(uint8_t* RESTRICT const mem, const int length, const ff_unit coef);
void init_ff();
void free_ff();
int ff_row_size(const int symbol_size);
int payload_size(const int k, const int symbol_size);
int SIMD_size();
intptr_t ceil_to_grid_p(intptr_t const arg);
ff_unit* get_coef_p(uint8_t* arg);
const uint8_t * get_gf2_8_multiplication_table();
const uint8_t * get_gf2_8_inversion_table();
int full_multiplication_table();

static inline intptr_t ceil_to_p(intptr_t const arg)
{
    const intptr_t modulo = arg & (8 - 1);
    const intptr_t res = (modulo == 0 ? arg : arg + 8 - modulo);
    return res;
}
