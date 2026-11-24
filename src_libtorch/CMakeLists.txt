cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(torch_loader LANGUAGES CXX)

option(CASCADIA_DOWNLOAD_LIBTORCH "Download the Linux CPU SDK if Torch is missing" ON)
find_package(Torch CONFIG QUIET)
if(NOT Torch_FOUND)
  string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" torch_processor)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND torch_processor MATCHES "^(x86_64|amd64)$")
    set(torch_url "https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip")
    set(torch_hash 137a842d1cf1e9196b419390133a1623ef92f8f84dc7a072f95ada684f394afd)
    set(torch_package_root ".")
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND torch_processor MATCHES "^(aarch64|arm64)$")
    # The wheel contains the C++ SDK under torch/ and its bundled libraries alongside it.
    set(torch_url "https://download-r2.pytorch.org/whl/torch-2.0.1-cp310-cp310-manylinux2014_aarch64.whl")
    set(torch_hash 359bfaad94d1cda02ab775dc1cc386d585712329bb47b8741607ef6ef4950747)
    set(torch_package_root "torch")
  else()
    message(FATAL_ERROR "No Torch SDK found for ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}. "
      "Automatic downloads support Linux x86_64 and ARM64. Set Torch_DIR to the directory containing TorchConfig.cmake.")
  endif()
  if(NOT CASCADIA_DOWNLOAD_LIBTORCH)
    message(FATAL_ERROR "No Torch SDK found and CASCADIA_DOWNLOAD_LIBTORCH is OFF in this build directory. "
      "Reconfigure with -U Torch_DIR -DCASCADIA_DOWNLOAD_LIBTORCH=ON to enable the download, "
      "or set Torch_DIR to an installed SDK directory containing TorchConfig.cmake.")
  endif()

  include(FetchContent)
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
  endif()
  message(STATUS "Fetching LibTorch 2.0.1 CPU (Linux ${torch_processor})")
  FetchContent_Declare(cascadia_libtorch
    URL "${torch_url}"
    DOWNLOAD_NAME libtorch-2.0.1-${torch_processor}.zip
    URL_HASH SHA256=${torch_hash}
    TLS_VERIFY ON)
  FetchContent_MakeAvailable(cascadia_libtorch)
  find_package(Torch CONFIG REQUIRED PATHS "${cascadia_libtorch_SOURCE_DIR}/${torch_package_root}" NO_DEFAULT_PATH)
endif()
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
