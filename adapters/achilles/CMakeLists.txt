# Achilles adapter. Disabled by default; enable with:
#   cmake -B build -DNUGEOM_ENABLE_ACHILLES=ON -DAchilles_DIR=<path-to-achilles-install/cmake>
#
# Requires an Achilles install that exposes the targets
#   achilles::Achilles   (core library)
# or, equivalently, the subdirectory form: pointing ACHILLES_SOURCE_DIR at
# a source tree whose top CMakeLists declares those targets.

option(NUGEOM_ENABLE_ACHILLES "Build the Achilles GeneratorInterface adapter" OFF)

if(NOT NUGEOM_ENABLE_ACHILLES)
    return()
endif()

find_package(Achilles CONFIG QUIET)
if(NOT Achilles_FOUND)
    message(FATAL_ERROR
        "NUGEOM_ENABLE_ACHILLES=ON but Achilles was not found. "
        "Pass -DAchilles_DIR=<path> or add Achilles to CMAKE_PREFIX_PATH.")
endif()

add_library(nugeom_achilles SHARED
    AchillesAdapter.cc
)
target_include_directories(nugeom_achilles PUBLIC
    "${CMAKE_SOURCE_DIR}"
)
target_link_libraries(nugeom_achilles
    PRIVATE project_options
    PUBLIC  geom achilles::Achilles)

# End-to-end driver: Achilles events through the NuGeometry detector.
add_executable(achilles_geom_driver achilles_driver.cc)
target_link_libraries(achilles_geom_driver PRIVATE project_options nugeom_achilles)

# Diagnostic: Achilles' sigma_tot spline vs the integral of its own generated
# weights, at fixed energy.  See mode_discrepancy_findings.md.
add_executable(xsec_scan xsec_scan.cc)
target_link_libraries(xsec_scan PRIVATE project_options nugeom_achilles)
