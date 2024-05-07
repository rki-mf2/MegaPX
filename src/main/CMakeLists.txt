#include lyra for argument parsing
set(Main_SRC main.cpp)
file(GLOB MAIN_HDRS "*.hpp")
add_executable (${PROJECT_NAME} ${Main_SRC} ${MAIN_HDRS})
find_package (Threads)
target_link_libraries (${PROJECT_NAME} pthread mutate seqan3::seqan3 hibf_seqan fm_index_mx)
target_compile_options(${PROJECT_NAME} PRIVATE -O2)


