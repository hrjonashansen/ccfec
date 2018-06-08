#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <ccfec/decoder.h>

#define STATUS_UNKNOWN 0
#define STATUS_KNOWN   1
#define STATUS_DECODED 2

void refresh_symbol_status_forward(decoder_t* const decoder, const int index)
{
    for (int i = index + 1; i < decoder->k; ++i)
    {
        const ff_unit c = get_coef_p(decoder->coef[index])[i];
        if (c != 0)
        {
            return;
        }
    }
    decoder->status[index] = STATUS_DECODED;
}

void refresh_symbol_status_backward(decoder_t* const decoder, const int index)
{
    for (int i = index - 1; i >= 0; --i)
    {
        const ff_unit c = get_coef_p(decoder->coef[index])[i];
        if (c != 0)
        {
            return;
        }
    }
    decoder->status[index] = STATUS_DECODED;
}

// TODO: Break up into smaller functions...
decode_info_t decode(decoder_t* decoder, uint8_t* RESTRICT const payload)
{
    decode_info_t info = {0, false};

    const int k = decoder->k;
    const int symbol_size = decoder->symbol_size;
    const int ps = payload_size(k, symbol_size - sizeof(int));

    for (int i = 0; i < k; ++i)
    {
        const ff_unit c = get_coef_p(payload)[i];
        if (c == 0)
        {
            continue;
        }
        else if (decoder->status[i] != STATUS_UNKNOWN)
//        else if (decoder->status[i] == STATUS_KNOWN)
        {
            encode_symbol(decoder->coef[i], payload, ps, c);
        }
//        else if (decoder->status[i] == STATUS_DECODED)
//        {
//            encode_symbol(decoder->coef[i], payload, *(int*)(decoder->data[i]) + k * sizeof(ff_unit) + sizeof(int), c);
//        }
        else if (info.innovative == false)
        {
            info.innovative = true;
            info.index = i;
        }
    }

    if (info.innovative == false)
    {
        return info;
    }
    ++decoder->rank;

    memcpy(decoder->coef[info.index], payload, k * sizeof(ff_unit));
    memcpy(decoder->data[info.index], payload + k * sizeof(ff_unit), symbol_size);
    decoder->status[info.index] = STATUS_KNOWN;

    const ff_unit c = get_coef_p(payload)[info.index];
    if (c != 1)
    {
        normalise_symbol(decoder->coef[info.index], ps, c);
    }

    refresh_symbol_status_forward(decoder, info.index);

    for (int i = 0; i < info.index; ++i)
    {
        if (decoder->status[i] == STATUS_KNOWN)
        {
            const ff_unit c = get_coef_p(decoder->coef[i])[info.index];
            if (c == 0)
            {
                continue;
            }
            encode_symbol(decoder->coef[info.index], decoder->coef[i], ps, c);
            refresh_symbol_status_forward(decoder, i);
        }
    }

    return info;
}

// TODO: Break up into smaller functions...
decode_info_t decode_backward(decoder_t* decoder, uint8_t* RESTRICT const payload)
{
    decode_info_t info = {0, false};

    const int k = decoder->k;
    const int symbol_size = decoder->symbol_size;
    const int ps = payload_size(k, symbol_size - sizeof(int));

    for (int i = k - 1; i >= 0; --i)
    {
        const ff_unit c = get_coef_p(payload)[i];
        if (c == 0)
        {
            continue;
        }
        else if (decoder->status[i] != STATUS_UNKNOWN)
        {
            encode_symbol(decoder->coef[i], payload, ps, c);
        }
        else if (info.innovative == false)
        {
            info.innovative = true;
            info.index = i;
        }
    }

    if (info.innovative == false)
    {
        return info;
    }
    ++decoder->rank;

    memcpy(decoder->coef[info.index], payload, k * sizeof(ff_unit));
    memcpy(decoder->data[info.index], payload + k * sizeof(ff_unit), symbol_size);
    decoder->status[info.index] = STATUS_KNOWN;

    const ff_unit c = get_coef_p(payload)[info.index];
    if (c != 1)
    {
        normalise_symbol(decoder->coef[info.index], ps, c);
    }

    refresh_symbol_status_backward(decoder, info.index);

    for (int i = info.index + 1; i < k; ++i)
    {
        if (decoder->status[i] == STATUS_KNOWN)
        {
            const ff_unit c = get_coef_p(decoder->coef[i])[info.index];
            if (c == 0)
            {
                continue;
            }
            encode_symbol(decoder->coef[info.index], decoder->coef[i], ps, c);
            refresh_symbol_status_backward(decoder, i);
        }
    }

    return info;
}

void init_decoder(decoder_t* const decoder, const int k, const int symbol_size)
{
    decoder->k = k;
    decoder->rank = 0;
    decoder->symbol_size = symbol_size + sizeof(int);

    const int coef_and_data_row_size = payload_size(k, decoder->symbol_size - sizeof(int));
    const int coef_and_data_block_size = k * coef_and_data_row_size;
    const int status_block_size = k;
    const int data_pointer_size = k * sizeof(uint8_t*);
    const int coef_pointer_size = k * sizeof(uint8_t*);
    const int total_size = SIMD_size() + 8 + coef_and_data_block_size + status_block_size + data_pointer_size + coef_pointer_size;

    decoder->p = (uint8_t* RESTRICT)malloc(total_size);
    if (!decoder->p)
    {
        abort();
    }
    uint8_t* RESTRICT const p_alligned_data = (uint8_t* RESTRICT)ceil_to_grid_p((intptr_t)decoder->p);
    decoder->status = p_alligned_data + coef_and_data_block_size;
    uint8_t* RESTRICT const p_alligned_pointers = (uint8_t* RESTRICT)ceil_to_p((intptr_t)(decoder->status + status_block_size));
    decoder->data = (uint8_t* RESTRICT *)(p_alligned_pointers);
    decoder->coef = (uint8_t* RESTRICT * RESTRICT)(p_alligned_pointers + data_pointer_size);
    for (int i = 0; i < k; ++i)
    {
        decoder->coef[i] = p_alligned_data + i * coef_and_data_row_size;
        decoder->data[i] = p_alligned_data + i * coef_and_data_row_size + k * sizeof(ff_unit);
        memset(decoder->data[i] + decoder->symbol_size, 0, ff_row_size(decoder->symbol_size) - decoder->symbol_size);
    }
    memset(decoder->status, 0, decoder->k);
}

const uint8_t* get_symbol(const decoder_t* const decoder, const int i)
{
    return decoder->data[i] + sizeof(int);
}

const uint8_t* get_symbol_with_size(const decoder_t* const decoder, const int i, int* const size_p)
{
    if (size_p)
    {
        memcpy(size_p, decoder->data[i], sizeof(int));
        assert(*size_p > 0 && *size_p <= decoder->symbol_size);
    }
    return decoder->data[i] + sizeof(int);
}

bool symbol_decoded(const decoder_t* const decoder, const int i)
{
    return decoder->status[i] == STATUS_DECODED;
}

bool decoding_is_complete(const decoder_t* const decoder)
{
    return decoder->rank == decoder->k;
}

void free_decoder(decoder_t* const decoder)
{
    free(decoder->p);
}

void reset_decoder(decoder_t* const decoder)
{
    memset(decoder->status, 0, decoder->k);
    decoder->rank = 0;
}

void print_status(const decoder_t* const decoder)
{
    for (int i = 0; i < decoder->k; ++i)
    {
        switch (decoder->status[i])
        {
        case STATUS_UNKNOWN:
            printf("%d: Unknown\n", i);
            break;
        case STATUS_KNOWN:
            printf("%d: Known\n", i);
            break;
        case STATUS_DECODED:
            printf("%d: Decoded\n", i);
            break;
        default:
            printf("%d: ERROR\n", i);
            abort();
        }
    }
}

void print_matrix(const decoder_t* const decoder)
{
    for (int i = 0; i < decoder->k; ++i)
    {
        printf("%d: ", i);
        for (int u = 0; u < decoder->k; ++u)
        {
            printf("%3d ", decoder->coef[i][u]);
        }
        printf("\n");
    }
}

void print_vector(const decoder_t* const decoder, uint8_t* RESTRICT const payload)
{
    const int k = decoder->k;

    for (int i = 0; i < k; ++i)
    {
        const ff_unit c = get_coef_p(payload)[i];
        printf(" %3d", c);
    }
    printf("\n");
}
