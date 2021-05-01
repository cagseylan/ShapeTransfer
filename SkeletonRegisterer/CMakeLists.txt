cmake_minimum_required(VERSION 3.1)
project(register)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"	  ON)

find_package(LIBIGL REQUIRED QUIET)

# Libigl requires a modern C++ compiler that supports c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/include/eigen3/ -std=c++11 -O3") 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
