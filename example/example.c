#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include <ccfec/ff.h>
#include <ccfec/encoder.h>
#include <ccfec/decoder.h>

int main()
{
    init_ff(); // Initialize the finite field arithmetic

    const int k = 8, n = 16, symbol_size = 1600, N = 100000;
    uint8_t payload[payload_size(k, symbol_size)];

    encoder_t e;
    decoder_t d;
    init_encoder(&e, k, symbol_size, n);
    init_rs_coef(&e);
//    init_systematic_rs_coef(&e);
    init_decoder(&d, k, symbol_size);

    uint8_t data[symbol_size * k];

    for (int i = 0; i < k; ++i)
    {
        memset(data + i * symbol_size, i + 0xf, symbol_size);
    }

    for (int u = 0; u < N; ++u)
    {
        for (int i = 0; i < k; ++i)
        {
            set_symbol(&e, data + i * symbol_size, i);
        }

        for (int i = 0; i < n; ++i)
        {
            encode(&e, payload);
//            print_vector(&d, payload);
            decode(&d, payload);
//            decode_backward(&d, payload);
        }

        reset_encoder(&e);
        reset_decoder(&d);
    }
    free_encoder(&e);
    free_decoder(&d);
    free_ff();
}
