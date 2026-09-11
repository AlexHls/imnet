# ============================================================================
# Source files
# ============================================================================

# C++ application sources
set(IMNET_SOURCES
    main.cpp
    network_wrapper.cpp
    file_io.cpp
    app_state.cpp
    ui_theme.cpp
    ui_main_window.cpp
    ui_nuclide_chart.cpp
)

# ============================================================================
# Define library: libnetwork
# ============================================================================
if(IMNET_NETWORK_BACKEND STREQUAL "YANN")
    # C network core library
    set(YANN_SOURCES
        ${YANN_DIR}/network.c
        ${YANN_DIR}/network_integrate.c
        ${YANN_DIR}/network_solver.c
        ${YANN_DIR}/network_nse.c
        ${YANN_DIR}/helm_eos.c
        ${YANN_DIR}/utilities.c
    )

    add_library(libnetwork STATIC ${YANN_SOURCES})

    # Set include directories for C core
    target_include_directories(libnetwork PUBLIC ${YANN_DIR})

    # Link GSL for C core
    target_link_libraries(libnetwork PUBLIC GSL::gsl GSL::gslcblas)

    if(MSVC)
        target_compile_options(libnetwork PRIVATE /W4)
    else()
        target_compile_options(libnetwork PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wno-unused-parameter
            -Wno-sign-compare
        )
    endif()

    # Necessary C defines for YANN
    target_compile_definitions(libnetwork PRIVATE
        EOS_HELM
        USE_GSL
        NETWORK_T_MAX=1e10
        _POSIX_C_SOURCE=200112L
    )
else()
    set(NUPPN_FRAME_SOURCE_DIR "${IMNET_NUPPN_SOURCE_DIR}/frames/ppn/source")
    set(NUPPN_PHYSICS_SOURCE_DIR "${IMNET_NUPPN_SOURCE_DIR}/physics/source")
    set(NUPPN_SOLVER_SOURCE_DIR "${IMNET_NUPPN_SOURCE_DIR}/solver/source")
    set(NUPPN_UTILS_SOURCE_DIR "${IMNET_NUPPN_SOURCE_DIR}/utils/source")
    set(NUPPN_SUPERLU_DIR "${NUPPN_SOLVER_SOURCE_DIR}/SuperLU_5.0")

    add_custom_target(nuppn_external
        COMMAND "${IMNET_MAKE_PROGRAM}" -C "${NUPPN_FRAME_SOURCE_DIR}" "PPN=${IMNET_NUPPN_SOURCE_DIR}"
        WORKING_DIRECTORY "${NUPPN_FRAME_SOURCE_DIR}"
        VERBATIM
    )

    add_library(libnetwork STATIC nuppn_c_api.F90)
    add_dependencies(libnetwork nuppn_external)
    target_include_directories(libnetwork PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NUPPN_FRAME_SOURCE_DIR}
        ${NUPPN_PHYSICS_SOURCE_DIR}
        ${NUPPN_SOLVER_SOURCE_DIR}
        ${NUPPN_SOLVER_SOURCE_DIR}/slu_f90
        ${NUPPN_UTILS_SOURCE_DIR}
    )
    target_compile_definitions(libnetwork PRIVATE PPN)
    target_link_libraries(libnetwork PUBLIC
        ${NUPPN_FRAME_SOURCE_DIR}/frame_knobs.o
        ${NUPPN_SOLVER_SOURCE_DIR}/slu_f90/slu.o
        ${NUPPN_SOLVER_SOURCE_DIR}/slu_f90/slu_f90.o
        ${NUPPN_SUPERLU_DIR}/FORTRAN/c_fortran_dgssv.o
        ${NUPPN_SOLVER_SOURCE_DIR}/libsolver.a
        ${NUPPN_PHYSICS_SOURCE_DIR}/libphysics.a
        ${NUPPN_UTILS_SOURCE_DIR}/libutils.a
        ${NUPPN_SUPERLU_DIR}/lib/libsuperlu_5.0.a
        blas
        m
    )
endif()

# ============================================================================
# Define executable: imnet
# ============================================================================
add_executable(imnet ${IMNET_SOURCES})

# Set include directories
target_include_directories(imnet PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link libraries
target_link_libraries(imnet PRIVATE
    libnetwork
    implot
)

target_compile_definitions(imnet PRIVATE
    IMNET_VERSION="${PROJECT_VERSION}"
)
if(IMNET_NETWORK_BACKEND STREQUAL "NUPPN")
    target_compile_definitions(imnet PRIVATE
        IMNET_USE_NUPPN
        IMNET_NUPPN_RUN_DIR="${IMNET_NUPPN_RUN_DIR}"
    )
endif()
if(APPLE)
    target_compile_definitions(imnet PRIVATE GL_SILENCE_DEPRECATION)
endif()

if(MSVC)
    target_compile_options(imnet PRIVATE /W4)
else()
    target_compile_options(imnet PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(BUILD_TESTING AND IMNET_NETWORK_BACKEND STREQUAL "YANN")
    add_executable(imnet_app_tests
        ${CMAKE_SOURCE_DIR}/tests/app_layer_tests.cpp
        file_io.cpp
        app_state.cpp
        network_wrapper.cpp
    )
    target_include_directories(imnet_app_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_link_libraries(imnet_app_tests PRIVATE libnetwork)
    target_compile_definitions(imnet_app_tests PRIVATE
        IMNET_TEST_DATA_DIR="${DATA_DIR}"
    )
    add_test(NAME imnet_app_tests COMMAND imnet_app_tests)
endif()

install(TARGETS imnet
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install data files
install(DIRECTORY ${DATA_DIR}/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/imnet/data
)
