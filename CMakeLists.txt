cmake_minimum_required(VERSION 3.12)

# Policy to address @foo@ variable expansion
if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
endif()

# Set the project name and basic settings
project(NUGEOMETRY LANGUAGES CXX VERSION 1.0.0)
include(CMake/StandardProjectSettings.cmake)

# Link this 'library' to set the c++ standard / compile-time options requested
# Additionally, link to get include and external dependencies
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)
target_include_directories(project_options INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Windows options
if(MSVC) 
    target_compile_definitions(project_options INTERFACE _USE_MATH_DEFINES)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# Standard compiler warnings
include(CMake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitizer options if supported by compiler
include(CMake/Sanitizers.cmake)
enable_sanitizers(project_options)

# Allow for static analysis options
include(CMake/StaticAnalyzers.cmake)

# Build options
option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)
option(ENABLE_TESTING "Enable Test Builds" OFF)
option(ENABLE_INTERACTIVE "Enable interactively moving around" OFF)

# Include external libraries
include(CMake/CPM.cmake)
add_subdirectory(external)

# Very basic PCH example
option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
    # This sets a global PCH parameter, each project will build its own PCH,
    # which is a good idea if any #define's change
    # These should be headers included in many places
    target_precompile_headers(project_options INTERFACE <vector> <string> <map>)
endif()

# Configure system specific variables for loading shared libraries
set(LIBPREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(LIBSUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
set(BUILDPATH ${CMAKE_CURRENT_BINARY_DIR})
set(INSTALLPATH ${CMAKE_INSTALL_PREFIX})
configure_file(include/geom/System.hh.in ${CMAKE_CURRENT_SOURCE_DIR}/include/geom/System.hh)

# Interactive
if(ENABLE_INTERACTIVE)
    find_package(SFML 2.5 COMPONENTS system window graphics REQUIRED)
endif()

# Testing
if(ENABLE_TESTING)
    if(COVERAGE) 
        target_compile_options(project_options INTERFACE --coverage) 
        target_link_libraries(project_options INTERFACE gcov)
    endif()
    target_compile_definitions(project_options INTERFACE TESTING)
    enable_testing()
    message(STATUS "Building Tests.")
    add_subdirectory(test)
endif()

if(ENABLE_FUZZING)
    message(STATUS 
        "Building Fuzz Tests, using fuzzing sanitizer https://www.llvm.org/docs/LibFuzzer.html"
    )
    add_subdirectory(fuzz_test)
endif()

# Add python support
if(ENABLE_PYTHON)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

# Main code
add_subdirectory(src)

# Install CMake find_package files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/CMake/nugeometry-config.cmake.in
    ${CMAKE_BINARY_DIR}/CMake/nugeometry-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nugeometry
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/CMake/nugeometry-config-version.cmake
    VERSION ${NUGEOMETRY_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/CMake/nugeometry-config.cmake
        ${CMAKE_BINARY_DIR}/CMake/nugeometry-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nugeometry
)
