cmake_minimum_required(VERSION 2.8)

project(delone)

add_subdirectory(meshProcessor)

add_executable(delone delone.cpp coloring.cpp LebedevQuad.cpp matrixSolution.cpp)
target_link_libraries(delone mesh3d)
