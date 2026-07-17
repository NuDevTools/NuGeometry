# Make imported targets from find_package (used by CPMFindPackage below) visible
# in sibling directories like src/, not just here in external/.
set(CMAKE_FIND_PACKAGE_TARGETS_GLOBAL TRUE)

# Prefer an already-available yaml-cpp (e.g. the one Achilles installs) so the
# Achilles adapter does not load two different yaml-cpp ABIs. Pinned to match
# Achilles' external/CMakeLists.txt.
CPMFindPackage(
    NAME yaml-cpp
    GITHUB_REPOSITORY jbeder/yaml-cpp
    GIT_TAG 4861d04
    OPTIONS
        "YAML_CPP_BUILD_TESTS OFF"
        "YAML_CPP_BUILD_CONTRIB OFF"
        "YAML_CPP_BUILD_TOOLS OFF"
        "YAML_CPP_INSTALL ON"
)
# yaml-cpp exports yaml-cpp::yaml-cpp when found and yaml-cpp when freshly built;
# expose NuGeom's yaml::cpp alias over whichever exists.
if(NOT TARGET yaml::cpp)
    add_library(nugeom_yaml_cpp INTERFACE)
    if(TARGET yaml-cpp::yaml-cpp)
        target_link_libraries(nugeom_yaml_cpp INTERFACE yaml-cpp::yaml-cpp)
    else()
        target_link_libraries(nugeom_yaml_cpp INTERFACE yaml-cpp)
    endif()
    add_library(yaml::cpp ALIAS nugeom_yaml_cpp)
endif()

CPMAddPackage(
    NAME pugixml
    GITHUB_REPOSITORY zeux/pugixml
    VERSION 1.16
)

# Version floors so find_package rejects an older parallel install (e.g. an
# ~/.local fmt 11 / spdlog 1.14) and resolves to the same system fmt 12 /
# spdlog 1.17 that the Achilles libraries were built against. Mismatched ABIs
# here load two copies of the library and crash at runtime.
CPMFindPackage(
    NAME fmt
    VERSION 12
    GITHUB_REPOSITORY fmtlib/fmt
    GIT_TAG 12.1.0
    OPTIONS
        "FMT_INSTALL ON"
)
if(NOT TARGET fmt::fmt)
    add_library(fmt::fmt ALIAS fmt)
endif()

CPMFindPackage(
    NAME spdlog
    VERSION 1.16
    GITHUB_REPOSITORY gabime/spdlog
    GIT_TAG v1.16.0
    OPTIONS
        "SPDLOG_COMPILED_LIB ON"
        "SPDLOG_FMT_EXTERNAL ON"
        "SPDLOG_BUILD_SHARED ON"
        "SPDLOG_INSTALL ON"
)
if(NOT TARGET spdlog::spdlog)
    add_library(spdlog::spdlog ALIAS spdlog)
endif()

CPMAddPackage(
    NAME cli11
    VERSION 2.6.1
    GITHUB_REPOSITORY CLIUtils/CLI11
    OPTIONS
)
add_library(CLI11::CLI11 ALIAS CLI11)

#pybind11
CPMAddPackage(
    NAME pybind11
    GITHUB_REPOSITORY pybind/pybind11
    GIT_TAG v2.9.2
)

# NuHepMC convention helpers.  With NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON it
# fetches and builds HepMC3 itself, exposing HepMC3::HepMC3 transitively.
# Override install() locally to swallow NuHepMC_CPPUtils' install(EXPORT)
# rules, which would otherwise error on the CPM-fetched spdlog target not
# being part of any export set.  We never install these deps anyway.
function(install)
endfunction()
CPMAddPackage(
    NAME NuHepMC_CPPUtils
    GITHUB_REPOSITORY NuHepMC/cpputils
    GIT_TAG v1-RC19
    OPTIONS
        "NuHepMC_CPPUtils_ENABLE_TESTS OFF"
        "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
)

if(ENABLE_INTERACTIVE)
    CPMAddPackage(
        NAME imgui
        GITHUB_REPOSITORY ocornut/imgui
        GIT_TAG v1.91.6
    )
    CPMAddPackage(
        NAME imgui-sfml
        GITHUB_REPOSITORY SFML/imgui-sfml
        GIT_TAG v3.0
        OPTIONS
            "IMGUI_DIR ${imgui_SOURCE_DIR}"
            "IMGUI_SFML_FIND_SFML OFF"
            "IMGUI_SFML_IMGUI_DEMO ON"
    )
endif()

# Install testing framework
if(ENABLE_TESTING OR ENABLE_FUZZING)
    # Catch2
    CPMAddPackage(
        NAME catch2
        GITHUB_REPOSITORY catchorg/Catch2
        VERSION 2.13.9
        OPTIONS
            "CATCH_ENABLE_COVERAGE"
    )
    add_library(Catch2::Catch2 ALIAS Catch2)

    # Trompeloeil
    CPMAddPackage(
        NAME trompeloeil
        GITHUB_REPOSITORY rollbear/trompeloeil
        VERSION 37 
    )
endif()
