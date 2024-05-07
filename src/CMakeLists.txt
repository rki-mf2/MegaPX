cmake_minimum_required (VERSION 3.16)

project (megaX CXX C)

add_compile_options(-DDEBUG)
add_compile_options(-DBUILD_SHARED_LIBS=ON)
add_compile_options(-Wno-format)


if (WIN32)
  set(CMAKE_GENERATORS "MinGW Makefiles")
  set(CMAKE_ASM_NASM_OBJECT_FORMAT win32)
  add_compile_options(-D_WIN32_WINNT=0x600)
  add_definitions(-D_WIN32_WINNT=0x600)
else()
  set(CMAKE_GENERATORS "Unix Makefiles")
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fopenmp")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")


if (WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/Apps/${PROJECT_NAME}")
    find_file(libbz2_DLL "libbz2-1__.dll")
    install(FILES ${libbz2_DLL} DESTINATION bin COMPONENT applications)
else()
    set(CMAKE_INSTALL_PREFIX "/usr/local/${PROJECT_NAME}")
endif()

# add subdirectories
#add_subdirectory(seqan3)
add_subdirectory(main)
add_subdirectory(mutate)
add_subdirectory(hibf)
add_subdirectory(fm_index_mx)
