#pragma once

#include <stdint.h>
#include <stdbool.h>

#include "ff.h"

typedef struct decoder_t
{
    int k;
    int symbol_size;
    int rank;
    uint8_t* RESTRICT * data;
    ff_unit* RESTRICT * RESTRICT coef;
    uint8_t* RESTRICT status;
    uint8_t* RESTRICT p;
} decoder_t;

typedef struct decode_info_t
{
    int index;
    bool innovative;
} decode_info_t;

decode_info_t decode(decoder_t* decoder, uint8_t* RESTRICT const payload);
decode_info_t decode_backward(decoder_t* decoder, uint8_t* RESTRICT const payload);
void init_decoder(decoder_t* const decoder, const int k, const int symbol_size);
void free_decoder(decoder_t* const decoder);
void reset_decoder(decoder_t* const decoder);
const uint8_t* get_symbol(const decoder_t* const decoder, const int i);
const uint8_t* get_symbol_with_size(const decoder_t* const decoder, const int i, int* const size_p);
bool symbol_decoded(const decoder_t* const decoder, const int i);
bool decoding_is_complete(const decoder_t* const decoder);

void print_status(const decoder_t* const decoder);
void print_matrix(const decoder_t* const decoder);
void print_vector(const decoder_t* const decoder, uint8_t* RESTRICT const payload);
