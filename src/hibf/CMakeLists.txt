set (seqan3_git_tag "main")
 
message (STATUS "Fetching SeqAn3 ${seqan3_git_tag}:")
 
include (FetchContent)
FetchContent_Declare (
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG "${seqan3_git_tag}"
)
 
FetchContent_MakeAvailable (seqan3_fetch_content)

FetchContent_Declare (
    hibf_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/hibf"
    GIT_TAG "main")
option (INSTALL_HIBF "" OFF)
FetchContent_MakeAvailable (hibf_fetch_content)

file(GLOB IBFSRCS "*.cpp")
file(GLOB IBFHDRS "*.hpp" "*.h")

find_package(OpenMP REQUIRED)

add_library (hibf_seqan STATIC ${IBFSRCS} ${IBFHDRS})

target_include_directories (hibf_seqan PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hibf_seqan PUBLIC OpenMP::OpenMP_CXX seqan3::seqan3 seqan::hibf mutate)

install(TARGETS hibf_seqan DESTINATION lib COMPONENT libraries)

install(FILES ${IBFHDRS} DESTINATION include COMPONENT headers)
