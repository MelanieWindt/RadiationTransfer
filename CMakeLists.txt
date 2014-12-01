cmake_minimum_required(VERSION 2.8)

project(shortChar)

add_subdirectory(meshProcessor)

add_executable(shortChar shortChar.cpp coloring.cpp LebedevQuad.cpp matrixSolution.cpp)
target_link_libraries(shortChar mesh3d)
