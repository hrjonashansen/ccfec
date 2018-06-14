#include <assert.h>

#include <ccfec/ff_basic.h>

uint8_t degree(uint8_t A)
{
    uint8_t res = 0;
    A >>= 1;
    while (A != 0)
    {
        A >>= 1;
        ++res;
    }
    return res;
}

uint8_t gf2_8_invert(uint8_t u)
{
    assert(u);

    uint8_t v = 0x1D;
    uint8_t g1 = 1, g2 = 0;
    int j = degree(u) - 8;

    while (u != 1)
    {
        if (j < 0)
        {
            // Swap u and v
            const uint8_t tmp_1 = u;
            u = v;
            v = tmp_1;

            // Swap g1 and g2
            const uint8_t tmp_2 = g1;
            g1 = g2;
            g2 = tmp_2;

            j = -j;
        }
        u ^= v << j;
        g1 ^= g2 << j;
        j = degree(u) - degree(v);
    }
    return g1;
}

uint8_t gf2_8_multiply(uint8_t lhs, uint8_t rhs)
{
    uint8_t res = 0;
    for (uint8_t i = 0; i < 8; ++i)
    {
        if (rhs & 1)
        {
            res ^= lhs;
        }
        // Detect if x^8 term is about to be generated
        const uint8_t carry = lhs & 0x80;
        lhs <<= 1;
        if (carry)
        {
            // Replace x^8 with x^4 + x^3 + x^2 + 1
            lhs ^= 0x1D;
        }
        rhs >>= 1;
    }
    return res;
}
