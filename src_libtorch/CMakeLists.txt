cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(torch_loader LANGUAGES CXX)

find_package(Torch REQUIRED)
find_package(ZLIB REQUIRED)

add_library(cascadia_runtime STATIC
  CascadiaPipeline.cpp
  CascadiaModelConfig.cpp
  CascadiaMzmlReader.cpp
  CascadiaPeptideTokenizer.cpp
  CascadiaSequenceForward.cpp
  CheckpointModelLoader.cpp
  InferenceParameterInspector.cpp
)
set_target_properties(cascadia_runtime PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_compile_options(cascadia_runtime PUBLIC ${TORCH_CXX_FLAGS})
target_link_libraries(cascadia_runtime PUBLIC ${TORCH_LIBRARIES} ZLIB::ZLIB)

# When a Conda environment is active, GCC can accidentally pick Conda's ld
# while using the system compiler and glibc. Prefer system binutils for this
# system-compiler build to avoid GLIBC_PRIVATE link errors.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER MATCHES "^/usr/bin/"
   AND EXISTS "/usr/bin/ld")
  target_link_options(cascadia_runtime PUBLIC "-B/usr/bin/")
endif()

option(CASCADIA_ENABLE_SANITIZERS "Enable native sanitizers (ASan+UBSan on Linux, UBSan with prebuilt macOS Torch)" OFF)
if(CASCADIA_ENABLE_SANITIZERS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    if(APPLE)
      # The local ASan build stalled at startup with prebuilt macOS Torch.
      # Linux CI exercises ASan+UBSan; local macOS checks use UBSan.
      set(CASCADIA_SANITIZER_FLAG -fsanitize=undefined)
    else()
      set(CASCADIA_SANITIZER_FLAG -fsanitize=address,undefined)
    endif()
    target_compile_options(cascadia_runtime PUBLIC ${CASCADIA_SANITIZER_FLAG} -fno-omit-frame-pointer)
    target_link_options(cascadia_runtime PUBLIC ${CASCADIA_SANITIZER_FLAG})
  else()
    message(FATAL_ERROR "CASCADIA_ENABLE_SANITIZERS requires GCC or Clang")
  endif()
endif()

target_include_directories(cascadia_runtime PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(torch_loader main.cpp)
target_link_libraries(torch_loader PRIVATE cascadia_runtime)
target_compile_features(torch_loader PRIVATE cxx_std_17)
add_executable(megapx_sequence sequence_main.cpp)
target_link_libraries(megapx_sequence PRIVATE cascadia_runtime)
target_compile_features(megapx_sequence PRIVATE cxx_std_17)

include(CTest)
if(BUILD_TESTING)
  add_executable(cascadia_native_tests ../tests/cascadia_native_tests.cpp)
  target_link_libraries(cascadia_native_tests PRIVATE cascadia_runtime)
  target_compile_features(cascadia_native_tests PRIVATE cxx_std_17)
  target_compile_definitions(cascadia_native_tests PRIVATE CASCADIA_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../tests/fixtures")
  add_test(NAME cascadia_native_tests COMMAND cascadia_native_tests)
  set_tests_properties(cascadia_native_tests PROPERTIES TIMEOUT 60)
  find_package(Python3 3.9 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME cascadia_benchmark_tests COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_cascadia_benchmark.py)
    set_tests_properties(cascadia_benchmark_tests PROPERTIES TIMEOUT 30)
  endif()
endif()
