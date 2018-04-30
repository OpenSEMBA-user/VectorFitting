cmake_minimum_required(VERSION 2.8)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

add_sources(. SRCS)

add_library(vectorfitting STATIC ${SRCS})