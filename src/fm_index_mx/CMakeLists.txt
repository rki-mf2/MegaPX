file(GLOB FMINDEXMXSRCS "*.cpp")
file(GLOB FMINDEXMXHDRS "*.hpp" "*.h")

find_package(OpenMP REQUIRED)

add_library (fm_index_mx STATIC ${FMINDEXMXSRCS} ${FMINDEXMXHDRS})


target_include_directories (fm_index_mx PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fm_index_mx PUBLIC OpenMP::OpenMP_CXX seqan3::seqan3 mutate)

install(TARGETS fm_index_mx DESTINATION lib COMPONENT libraries)

install(FILES ${FMINDEXMXHDRS} DESTINATION include COMPONENT headers)