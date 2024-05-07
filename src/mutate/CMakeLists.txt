file(GLOB MUTATESRCS "*.cpp")
file(GLOB MUTATEHDRS "*.hpp" "*.h")

find_package(OpenMP REQUIRED)

add_library (mutate STATIC ${MUTATESRCS} ${MUTATEHDRS})

#target_link_libraries (Mutate seqan3::seqan3)

target_include_directories (mutate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mutate PUBLIC OpenMP::OpenMP_CXX seqan3::seqan3 hibf_seqan )

install(TARGETS mutate DESTINATION lib COMPONENT libraries)

install(FILES ${MUTATEHDRS} DESTINATION include COMPONENT headers)
