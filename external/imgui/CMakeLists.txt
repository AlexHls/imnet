# ============================================================================
# imgui - Immediate Mode GUI Library
# Vendored sub-project for imyann
# ============================================================================

# Find dependencies
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLFW REQUIRED glfw3)
find_package(OpenGL REQUIRED)

add_library(imgui STATIC
    imgui.cpp
    imgui_demo.cpp
    imgui_draw.cpp
    imgui_tables.cpp
    imgui_widgets.cpp
    backends/imgui_impl_glfw.cpp
    backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/backends
    ${GLFW_INCLUDE_DIRS}
)

target_link_directories(imgui PUBLIC ${GLFW_LIBRARY_DIRS})

target_link_libraries(imgui PUBLIC
    ${GLFW_LIBRARIES}
    OpenGL::OpenGL
)

# Compile options
if(NOT MSVC)
    target_compile_options(imgui PRIVATE -Wno-unused-parameter -Wno-deprecated-declarations)
endif()

# Export library for main project
set(IMGUI_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
