#if defined(__AVX2__)
#include "ff_avx2.c"
#elif defined(__SSSE3__)
#include "ff_ssse3.c"
#elif defined(__ARM_NEON_FP)
#include "ff_neon.c"
#elif CCFEC_FF_ONLINE
#include "ff_online.c"
#else
#include "ff_simple.c"
#endif
