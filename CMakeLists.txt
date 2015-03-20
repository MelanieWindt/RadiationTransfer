cmake_minimum_required(VERSION 2.8)

project(shortChar)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
add_subdirectory(meshProcessor)

add_executable(shortChar shortChar.cpp coloring.cpp LebedevQuad.cpp matrixSolution.cpp)
target_link_libraries(shortChar mesh3d)
