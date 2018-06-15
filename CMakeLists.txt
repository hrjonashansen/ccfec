cmake_minimum_required(VERSION 3.0)

project(ccfec)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")

if (NOT AROS_BUILD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -flto -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -march=native -flto -Wall -Wextra -pedantic")
endif()

add_library(ccfec STATIC
    encoder.c
    decoder.c
    ff.c
    ff_basic.c
)

target_include_directories(ccfec PUBLIC "include")

if (NOT AROS_BUILD)
    if (CCFEC_TEST_ALL_FF_IMPLEMENTATIONS)
        include(CheckCXXSourceRuns)

        # Check whether or not SSSE3 instructions are available
        check_cxx_source_runs("
            #include <immintrin.h>
            int main()
            {
                const int src[4] = {0, 1, 2, 3};
                int dst[4];
                const __m128i a = _mm_load_si128(reinterpret_cast<const __m128i*>(src));
                const __m128i b = _mm_add_epi32(a, a);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(dst), b);
                for (int i = 0; i < 4; i++)
                {
                    if (src[i] * 2 != dst[i])
                    {
                        return 1;
                    }
                }
                return 0;
            }"
            SSSE3_INSTRUCTIONS_ENABLED)

        # Check whether or not AVX2 instructions are available
        check_cxx_source_runs("
            #include <immintrin.h>
            int main()
            {
                const int src[8] = {0, 1, 2, 3, 4, 5, 6, 7};
                int dst[8];
                const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src));
                const __m256i b = _mm256_add_epi32(a, a);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst), b);
                for (int i = 0; i < 8; i++)
                {
                    if (src[i] * 2 != dst[i])
                    {
                        return -1;
                    }
                }
                return 0;
            }"
            AVX2_INSTRUCTIONS_ENABLED)

        # Check whether or not NEON instructions are available (When compiled on x86/x64 this code requires "-march=armv7-a -mfpu=neon -mfloat-abi=hard")
        check_cxx_source_runs("
            #include <arm_neon.h>
            int main()
            {
                const int32_t src[8] = {0, 1, 2, 3};
                int32_t dst[8];
                const int32x4_t a = vld1q_s32(src);
                const int32x4_t b = vaddq_s32(a, a);
                vst1q_s32(dst, b);
                for (int i = 0; i < 4; i++)
                {
                    if (src[i] * 2 != dst[i])
                    {
                        return -1;
                    }
                }
                return 0;
            }"
            NEON_INSTRUCTIONS_ENABLED)

        add_library(ccfec-online STATIC
            encoder.c
            decoder.c
            ff_online.c
            ff_basic.c
        )
        target_include_directories(ccfec-online PUBLIC "include")

        add_library(ccfec-simple STATIC
            encoder.c
            decoder.c
            ff_simple.c
            ff_basic.c
        )
        target_include_directories(ccfec-simple PUBLIC "include")

        if (SSSE3_INSTRUCTIONS_ENABLED)
            add_library(ccfec-ssse3 STATIC
                encoder.c
                decoder.c
                ff_ssse3.c
                ff_basic.c
            )
            target_include_directories(ccfec-ssse3 PUBLIC "include")
        endif()

        if (AVX2_INSTRUCTIONS_ENABLED)
            add_library(ccfec-avx2 STATIC
                encoder.c
                decoder.c
                ff_avx2.c
                ff_basic.c
            )
            target_include_directories(ccfec-avx2 PUBLIC "include")
        endif()

        if (NEON_INSTRUCTIONS_ENABLED)
            add_library(ccfec-neon STATIC
                encoder.c
                decoder.c
                ff_neon.c
                ff_basic.c
            )
            target_include_directories(ccfec-neon PUBLIC "include")
        endif()
    endif() # CCFEC_TEST_ALL_FF_IMPLEMENTATIONS

    add_subdirectory(example)
    add_subdirectory(test)
else()
    add_subdirectory(test EXCLUDE_FROM_ALL)
endif()
