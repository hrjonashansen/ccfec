#pragma once

#include <stdint.h>

#include "ff.h"

typedef struct encoder_t
{
    int k;
    int n;
    int symbol_size;
    int next;
    uint8_t* RESTRICT * data;
    ff_unit* RESTRICT * RESTRICT coef;
    uint8_t* RESTRICT status;
    uint8_t* RESTRICT p;
} encoder_t;

int encode(encoder_t* const encoder, uint8_t* RESTRICT const payload);
void init_encoder(encoder_t* const encoder, const int k, const int symbol_size, const int n);
void free_encoder(encoder_t* const encoder);
void reset_encoder(encoder_t* const encoder);
int set_symbol(encoder_t* const encoder, uint8_t* RESTRICT const payload, const int id);
int set_next_symbol_with_size(encoder_t* const encoder, uint8_t* RESTRICT const payload, const int size);

void init_rs_coef(encoder_t* const encoder);
void init_systematic_rs_coef(encoder_t* const encoder);
