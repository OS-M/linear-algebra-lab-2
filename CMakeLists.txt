cmake_minimum_required(VERSION 3.20)
project(linear_algebra_lab_2)

add_compile_options(
        "$<$<CONFIG:RELEASE>:-O3>"
        "$<$<CONFIG:DEBUG>:-O0>"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_INCLUDE_CURRENT_DIR true)

# ThreadSanitizer
set(TSAN_FLAGS
        "-fsanitize=thread"
        CACHE STRING "Flags used by the C++ compiler during ThreadSanitizer builds."
        FORCE)

# AddressSanitize
set(ASAN_FLAGS
        "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during AddressSanitizer builds."
        FORCE)

# LeakSanitizer
set(LSAN_FLAGS
        "-fsanitize=leak -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during LeakSanitizer builds."
        FORCE)

# MemorySanitizer
set(MSAN_FLAGS
        "-fsanitize=memory -fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer"
        CACHE STRING "Flags used by the C++ compiler during MemorySanitizer builds."
        FORCE)

# UndefinedBehaviour
set(UBSAN_FLAGS
        "-fsanitize=undefined"
        CACHE STRING "Flags used by the C++ compiler during UndefinedBehaviourSanitizer builds."
        FORCE)

# Build Types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${ASAN_FLAGS}")

add_executable(linear_algebra_lab_2 main.cpp Matrix/matrix.h Algebra/gauss.h Algebra/gauss_back_substitution.h Algebra/euclidean_norm.h Algebra/qr_hessenberg.h Algebra/qr_algorithm.h Algebra/rotations.h Algebra/power_iteration_method.h Algebra/qr_decompose.h Algebra/minimal_square_problem.h Algebra/eigenvalues.h Algebra/frobenius_form.h Algebra/polynomial.h Algebra/danilevski_eigenvalues.h Algebra/polynomial_roots.h)
