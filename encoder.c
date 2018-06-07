#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

#include <ccfec/ff_online.h>
#include <ccfec/encoder.h>

int encode(encoder_t* const encoder, uint8_t* RESTRICT const payload)
{
    const int k = encoder->k;
    const int symbol_size = encoder->symbol_size;
//    const int payload_size = k + symbol_size;
    int real_size = 0;

    ff_unit* RESTRICT const coef = encoder->coef[encoder->next];

    memset(payload + k * sizeof(ff_unit), 0, ff_row_size(symbol_size));
    memcpy(payload, coef, k * sizeof(ff_unit));

    for (int i = 0; i < k; ++i)
    {
        if (encoder->status[i])
        {
            const ff_unit c = coef[i];
            if (c == 0)
            {
                continue;
            }
            int size_temp;
            memcpy(&size_temp, encoder->data[i], sizeof(int));
            if (size_temp > real_size)
                real_size = size_temp;
            assert(size_temp <= symbol_size - sizeof(int));
            encode_symbol(encoder->data[i], payload + k * sizeof(ff_unit), size_temp + sizeof(int), c);
        }
    }

    ++encoder->next;

    assert(real_size > 0);

    return real_size + k * sizeof(ff_unit) + sizeof(int);
}

void init_encoder(encoder_t* const encoder, const int k, const int symbol_size, const int n)
{
    encoder->k = k;
    encoder->n = n;
    encoder->next = 0;
    encoder->symbol_size = symbol_size + sizeof(int);

    const int data_block_size = k * ff_row_size(encoder->symbol_size);
    const int coef_block_size = n * k * sizeof(ff_unit);
    const int status_block_size = k;
    const int data_pointer_size = k * sizeof(uint8_t*);
    const int coef_pointer_size = n * sizeof(uint8_t*);
    const int total_size = SIMD_size() + 8 + data_block_size + coef_block_size + status_block_size + data_pointer_size + coef_pointer_size;

    encoder->p = (uint8_t* RESTRICT)malloc(total_size);
    uint8_t* RESTRICT const p_alligned_data = (uint8_t* RESTRICT)ceil_to_grid_p((intptr_t)encoder->p);
    encoder->status = p_alligned_data + data_block_size + coef_block_size;
    uint8_t* RESTRICT const p_alligned_pointers = (uint8_t* RESTRICT)ceil_to_p((intptr_t)(encoder->status + status_block_size));
    encoder->data = (uint8_t* RESTRICT *)(p_alligned_pointers);
    encoder->coef = (uint8_t* RESTRICT * RESTRICT)(p_alligned_pointers + data_pointer_size);
    if (!encoder->p)
    {
        abort();
    }
    for (int i = 0; i < k; ++i)
    {
        encoder->data[i] = p_alligned_data + i * ff_row_size(encoder->symbol_size);
        memset(encoder->data[i] + encoder->symbol_size, 0, ff_row_size(encoder->symbol_size) - encoder->symbol_size);
    }
    for (int i = 0; i < n; ++i)
    {
        encoder->coef[i] = p_alligned_data + i * k * sizeof(ff_unit) + data_block_size;
    }
    memset(encoder->status, 0, encoder->k);
}

void init_rs_coef(encoder_t* const encoder)
{
    for (int i = 0; i < encoder->n; ++i)
    {
        encoder->coef[i][0] = 1;
        for (int u = 1; u < encoder->k; ++u)
        {
            encoder->coef[i][u] = gf2_8_multiply(i + 1, encoder->coef[i][u - 1]);
        }
    }
}

void init_systematic_rs_coef(encoder_t* const encoder)
{
    // Initialize the encoding matrix
    init_rs_coef(encoder);

    // Allocate space for the inversion matrix
    ff_unit** const matrix = (ff_unit**)malloc(encoder->n * sizeof(ff_unit*));
    ff_unit** const matrix_to_free = (ff_unit**)malloc(encoder->n * sizeof(ff_unit*));
    for (int i = 0; i < encoder->n; ++i)
    {
        matrix_to_free[i] = (ff_unit*)malloc(ff_row_size(2 * encoder->k * sizeof(ff_unit)) + SIMD_size());
        matrix[i] = (ff_unit*)ceil_to_grid_p((intptr_t)matrix_to_free[i]);
    }

    // Copy the top k x k matrix from the encoder (and add the identity matrix)
    for (int r = 0; r < encoder->k; ++r)
    {
        for (int c = 0; c < encoder->k; ++c)
        {
            matrix[r][c] = encoder->coef[r][c];
            matrix[r][c + encoder->k] = c == r ? 1 : 0;
        }
    }

    // Zero out the lower part of the matrix
    for (int r = encoder->k; r < encoder->n; ++r)
    {
        for (int c = 0; c < 2* encoder->k; ++c)
        {
            matrix[r][c] = 0;
        }
    }

    // Forward substitute (first part of inverting the left half of the matrix)
    for (int r = 0; r < encoder->k; ++r)
    {
        normalise_symbol(matrix[r], 2 * encoder->k, matrix[r][r]);
        for (int u = r + 1; u < encoder->k; ++u)
        {
            if (matrix[u][r])
            {
                encode_symbol(matrix[r], matrix[u], 2 * encoder->k, matrix[u][r]);
            }
        }
    }

    // Backward substitute (second part of inverting the left half of the matrix)
    for (int r = encoder->k - 1; r > 0 ; --r)
    {
        normalise_symbol(matrix[r], 2 * encoder->k, matrix[r][r]);
        for (int u = 0; u < r; ++u)
        {
            if (matrix[u][r])
            {
                encode_symbol(matrix[r], matrix[u], 2 * encoder->k, matrix[u][r]);
            }
        }
    }

    // Prepare the left half for multiplication, by setting the diagonal to zero (all other elements are already zero)
    for (int r = 0; r < encoder->k; ++r)
    {
        matrix[r][r] = 0;
    }

    // Matrix multiplication, multiply the original matrix (full length) with the inverse of the upper k x k matrix.
    for (int r = 0; r < encoder->n; ++r)
    {
        for (int c = 0; c < encoder->k; ++c)
        {
            for (int k = 0; k < encoder->k; ++k)
            {
                matrix[r][c] ^= gf2_8_multiply(encoder->coef[r][k], matrix[k][c + encoder->k]);
            }
        }
    }

    // Copy the new coefficients back
    for (int r = 0; r < encoder->n; ++r)
    {
        for (int c = 0; c < encoder->k; ++c)
        {
            encoder->coef[r][c] = matrix[r][c];
        }
    }

    // Free the space for the inversion matrix
    for (int i = 0; i < encoder->n; ++i)
    {
        free(matrix_to_free[i]);
    }
    free(matrix_to_free);
    free(matrix);
}

void free_encoder(encoder_t* const encoder)
{
    free(encoder->p);
}

void reset_encoder(encoder_t* const encoder)
{
    memset(encoder->status, 0, encoder->k);
    encoder->next = 0;
}

int set_symbol(encoder_t* const encoder, uint8_t* RESTRICT const payload, const int id)
{
    assert(encoder->status[id] == 0);
    const int size = encoder->symbol_size - sizeof(int);
    memcpy(encoder->data[id] + sizeof(int), payload, size);
    memcpy(encoder->data[id], &size, sizeof(int));
    encoder->status[id] = 1;
    return encoder->symbol_size;
}

int set_next_symbol_with_size(encoder_t* const encoder, uint8_t* RESTRICT const payload, const int size)
{
    for (int i = 0; i < encoder->k; ++i)
    {
        if (encoder->status[i] == 0)
        {
            assert(size > 0 && size <= encoder->symbol_size - sizeof(int));
            memcpy(encoder->data[i] + sizeof(int), payload, size);
            if (encoder->symbol_size != size)
            {
                memset(encoder->data[i] + size + sizeof(int), 0, encoder->symbol_size - size - sizeof(int));
            }
            memcpy(encoder->data[i], &size, sizeof(int));
            encoder->status[i] = 1;
            return i;
        }
    }
    assert(false);
    return encoder->k;
}
