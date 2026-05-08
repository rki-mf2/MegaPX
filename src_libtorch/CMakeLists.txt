cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(torch_loader LANGUAGES CXX)

find_package(Torch REQUIRED)

add_executable(torch_loader
  main.cpp
  CheckpointModelLoader.cpp
)
target_compile_features(torch_loader PRIVATE cxx_std_20)
target_compile_options(torch_loader PRIVATE ${TORCH_CXX_FLAGS})
target_link_libraries(torch_loader PRIVATE ${TORCH_LIBRARIES})

# When a Conda environment is active, GCC can accidentally pick Conda's ld
# while using the system compiler and glibc. Prefer system binutils for this
# system-compiler build to avoid GLIBC_PRIVATE link errors.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER MATCHES "^/usr/bin/"
   AND EXISTS "/usr/bin/ld")
  target_link_options(torch_loader PRIVATE "-B/usr/bin/")
endif()
