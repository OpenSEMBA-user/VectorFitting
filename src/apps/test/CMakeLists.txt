cmake_minimum_required(VERSION 2.8)

if (GTEST_FOUND)
    project(vectorfitting_test CXX)
    include_directories(${GTEST_INCLUDE_DIRS}
                        ${CMAKE_CURRENT_LIST_DIR})
                        
    add_subdirectories(./ ./obj)
    add_sources(. SRCS)

    add_executable(vectorfitting_test ${SRCS})
    target_link_libraries(vectorfitting_test vectorfitting
                                             ${GTEST_BOTH_LIBRARIES}) 
endif()
