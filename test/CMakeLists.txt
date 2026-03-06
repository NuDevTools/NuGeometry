# Load Catch2 information into variables
include(FetchContent)
FetchContent_GetProperties(catch2)
include(${catch2_SOURCE_DIR}/contrib/Catch.cmake)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2 trompeloeil)

add_executable(nugeometry-testsuite
    # Files with tests
    test_translations.cc
    test_rotations.cc
    test_element.cc
    test_material.cc
    test_shape.cc
    test_volume.cc
    test_parser.cc
    test_bvh.cc
    bench_segments.cc
)
target_link_libraries(nugeometry-testsuite PRIVATE project_options project_warnings catch_main 
                                           PUBLIC geom spdlog::spdlog)

target_compile_definitions(nugeometry-testsuite PRIVATE
    NUGEOM_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

include(CTest)
catch_discover_tests(nugeometry-testsuite)
