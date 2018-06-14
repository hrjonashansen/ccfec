#pragma once

#include <stdint.h>

#ifdef __cplusplus
#define RESTRICT
#else
#define RESTRICT restrict
#endif

typedef uint8_t ff_unit;

void encode_symbol(uint8_t* RESTRICT const dst, const uint8_t* RESTRICT const src, const int length, const ff_unit coef);
void normalise_symbol(uint8_t* RESTRICT const mem, const int length, const ff_unit coef);
void init_ff();
void free_ff();
int ff_row_size(const int symbol_size);
int payload_size(const int k, const int symbol_size);
int SIMD_size();
uint8_t* ceil_to_grid_p(uint8_t* const arg_p);
ff_unit* get_coef_p(uint8_t* arg);
const uint8_t * get_gf2_8_multiplication_table();
const uint8_t * get_gf2_8_inversion_table();
int full_multiplication_table();

static inline uint8_t* ceil_to_p(uint8_t* RESTRICT const arg_p)
{
    const int pointer_alignment = 8;
    const intptr_t arg = (intptr_t)arg_p;
    const intptr_t modulo = arg & (pointer_alignment - 1);
    const intptr_t res = (modulo == 0 ? arg : arg + pointer_alignment - modulo);
    uint8_t* const res_p = (uint8_t*)res;
    return res_p;
}
