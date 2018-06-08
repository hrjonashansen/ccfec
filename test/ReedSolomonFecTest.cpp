#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include <boost/test/unit_test.hpp>

extern "C"
{
#include <ccfec/ff.h>
#include <ccfec/encoder.h>
#include <ccfec/decoder.h>
}

BOOST_AUTO_TEST_SUITE(reedsolomonfec)

BOOST_AUTO_TEST_CASE(ff_init_free)
{
    init_ff(); // Initialize the finite field arithmetic
    free_ff();
}

BOOST_AUTO_TEST_CASE(ff_tables)
{
    if (not full_multiplication_table())
    {
        return;
    }

    init_ff();

    const uint8_t * const invert = get_gf2_8_inversion_table();

    for (int a = 0; a < 256; ++a)
    {
        const uint8_t * const multiply_a = get_gf2_8_multiplication_table() + (a << 8);
        BOOST_CHECK_EQUAL(a, multiply_a[1]);
        BOOST_CHECK_EQUAL(0, multiply_a[0]);

        for (int b = 1; b < 256; ++b)
        {
            uint8_t product = multiply_a[b];
            const uint8_t * const multiply_product = get_gf2_8_multiplication_table() + (product << 8);

            BOOST_CHECK_EQUAL(a, multiply_product[invert[b]]);
        }
    }

    free_ff();
}

BOOST_AUTO_TEST_CASE(encoder_init_free)
{
    const int k = 8;
    const int n = 16;
    const int symbol_size = 1280;
    encoder_t e;

    init_encoder(&e, k, symbol_size, n);
    free_encoder(&e);
}

BOOST_AUTO_TEST_CASE(encoder_init_coef_free)
{
    const int k = 8;
    const int n = 16;
    const int symbol_size = 1280;
    encoder_t e;

    init_ff();

    init_encoder(&e, k, symbol_size, n);
    init_rs_coef(&e);
    free_encoder(&e);

    free_ff();
}

BOOST_AUTO_TEST_CASE(encoder_init_systematic_coef_free)
{
    const int k = 8;
    const int n = 16;
    const int symbol_size = 1280;
    encoder_t e;

    init_ff();

    init_encoder(&e, k, symbol_size, n);
    init_systematic_rs_coef(&e);
    free_encoder(&e);

    free_ff();
}

BOOST_AUTO_TEST_CASE(decoder_init_free)
{
    const int k = 8;
    const int symbol_size = 1280;
    decoder_t d;

    init_decoder(&d, k, symbol_size);
    free_decoder(&d);
}

void test_reed_solomon(const int k, const int symbol_size, const bool systematic, const bool lossy)
{
    const int n = 2 * k;
    const int N = 2;

    std::vector<uint8_t> payload(payload_size(k, symbol_size));

    encoder_t e;
    decoder_t d;
    init_encoder(&e, k, symbol_size, n);
    if (systematic)
    {
        init_systematic_rs_coef(&e);
    }
    else
    {
        init_rs_coef(&e);
    }
    init_decoder(&d, k, symbol_size);

    std::vector<uint8_t> data(symbol_size * k);

    for (int i = 0; i < k; ++i)
    {
        std::fill_n(data.begin() + i * symbol_size, symbol_size, i + 0xf);
    }

    for (int u = 0; u < N; ++u)
    {
        for (int i = 0; i < k; ++i)
        {
            set_symbol(&e, data.data() + i * symbol_size, i);
        }

        for (int i = 0; i < n; ++i)
        {
            encode(&e, payload.data());

            if (lossy and (i & 1))
            {
                continue; 
            }

            decode(&d, payload.data());
        }

        for (int i = 0; i < k; ++i)
        {
            const uint8_t * const original_symbol = data.data() + i * symbol_size;
            const uint8_t * const decoded_symbol = get_symbol(&d, i);
            BOOST_CHECK_EQUAL_COLLECTIONS(original_symbol, original_symbol + symbol_size, decoded_symbol, decoded_symbol + symbol_size);
        }

        reset_encoder(&e);
        reset_decoder(&d);
    }
    free_encoder(&e);
    free_decoder(&d);
}

BOOST_AUTO_TEST_CASE(encoder_decoder_full)
{
    init_ff();

    for (int k = 1; k < 16; ++k)
    {
        for (int symbol_size = 1; symbol_size < 1600; symbol_size += k)
        {
            test_reed_solomon(k, symbol_size, false, false);
            test_reed_solomon(k, symbol_size, true, false);
        }
    }

    for (int k = 1; k < 16; k += 2)
    {
        test_reed_solomon(k, 101, true, true);
    }

    free_ff();
}

BOOST_AUTO_TEST_SUITE_END() // reedsolomonfec
