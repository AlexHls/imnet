cmake_minimum_required(VERSION 3.16)
project(imnet VERSION 0.1.0 LANGUAGES C CXX)
include(CTest)
include(GNUInstallDirs)

set(IMNET_NETWORK_BACKEND "YANN" CACHE STRING "Network backend: YANN or NUPPN")
set_property(CACHE IMNET_NETWORK_BACKEND PROPERTY STRINGS YANN NUPPN)
string(TOUPPER "${IMNET_NETWORK_BACKEND}" IMNET_NETWORK_BACKEND)
if(NOT IMNET_NETWORK_BACKEND MATCHES "^(YANN|NUPPN)$")
    message(FATAL_ERROR "IMNET_NETWORK_BACKEND must be YANN or NUPPN")
endif()

if(IMNET_NETWORK_BACKEND STREQUAL "YANN")
    set(IMNET_YANN_SOURCE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/yann"
        CACHE PATH "Path to a Yann source checkout")
    get_filename_component(IMNET_YANN_SOURCE_DIR "${IMNET_YANN_SOURCE_DIR}" ABSOLUTE)
    foreach(YANN_REQUIRED_FILE
            network.c network_integrate.c network_solver.c network_nse.c
            helm_eos.c utilities.c network.h network_integrate.h
            network_solver.h network_nse.h helm_eos.h utilities.h)
        if(NOT EXISTS "${IMNET_YANN_SOURCE_DIR}/${YANN_REQUIRED_FILE}")
            message(FATAL_ERROR
                "Set IMNET_YANN_SOURCE_DIR to a Yann source directory; missing ${YANN_REQUIRED_FILE}")
        endif()
    endforeach()
    set(YANN_DIR "${IMNET_YANN_SOURCE_DIR}")
elseif(IMNET_NETWORK_BACKEND STREQUAL "NUPPN")
    enable_language(Fortran)
    set(IMNET_NUPPN_SOURCE_DIR "" CACHE PATH "Path to a NuPPN source checkout")
    if(NOT IMNET_NUPPN_SOURCE_DIR)
        message(FATAL_ERROR "Set IMNET_NUPPN_SOURCE_DIR when IMNET_NETWORK_BACKEND=NUPPN")
    endif()
    get_filename_component(IMNET_NUPPN_SOURCE_DIR "${IMNET_NUPPN_SOURCE_DIR}" ABSOLUTE)
    set(IMNET_NUPPN_RUN_DIR
        "${IMNET_NUPPN_SOURCE_DIR}/frames/ppn/run_template"
        CACHE PATH "NuPPN run directory containing ppn_*.input and isotopedatabase.txt")
    find_program(IMNET_MAKE_PROGRAM NAMES gmake make REQUIRED)
endif()

# C and C++ standards
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# On macOS, add Homebrew paths to CMAKE_PREFIX_PATH for GLFW
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH /opt/homebrew/opt/glfw)
endif()

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ============================================================================
# Dependency: GSL (GNU Scientific Library)
# ============================================================================
find_package(GSL REQUIRED)
find_package(OpenGL REQUIRED)
if(APPLE AND NOT TARGET OpenGL::OpenGL)
    add_library(OpenGL::OpenGL INTERFACE IMPORTED)
    target_link_libraries(OpenGL::OpenGL INTERFACE "-framework OpenGL")
endif()

# ============================================================================
# Sub-projects
# ============================================================================
add_subdirectory(external/imgui)

add_library(implot STATIC
    external/implot/implot.cpp
    external/implot/implot_items.cpp
)

target_include_directories(implot PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/implot
)

target_link_libraries(implot PUBLIC imgui)

if(MSVC)
    target_compile_options(implot PRIVATE /W4)
else()
    target_compile_options(implot PRIVATE
        -Wall -Wextra -Wpedantic -Wno-unused-parameter
    )
endif()

# ============================================================================
# Define paths
# ============================================================================
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

# ============================================================================
# Subdirectories
# ============================================================================
add_subdirectory(src)
