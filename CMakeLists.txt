cmake_minimum_required(VERSION 3.20)
project(EstimatorsForMeasuresOfPresortedness CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(${CMAKE_SOURCE_DIR}/external/matplotplusplus/source)


add_executable(EstimatorsForMeasuresOfPresortedness
        main.cpp
)


target_include_directories(EstimatorsForMeasuresOfPresortedness PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/Analysis
        ${CMAKE_SOURCE_DIR}/src/Data
)


include(FetchContent)
set(gtest_force_shared_crt ON CACHE BOOL "Use shared CRT" FORCE)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

enable_testing()
include(GoogleTest)

add_executable(EstimatorsForMeasuresOfPresortednessTest
        tests/test_runner.cpp
        tests/PresortednessMeasures.h
        tests/SamplerTest.h
        tests/DataGeneratorTest.h
)



target_include_directories(EstimatorsForMeasuresOfPresortednessTest PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/tests
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${GTEST_INC_CANDIDATES}
)
target_include_directories(EstimatorsForMeasuresOfPresortednessTest PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/Analysis
        ${CMAKE_SOURCE_DIR}/src/Data
)

target_link_libraries(EstimatorsForMeasuresOfPresortednessTest PRIVATE
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(EstimatorsForMeasuresOfPresortednessTest)
