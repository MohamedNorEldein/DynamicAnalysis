cmake_minimum_required(VERSION 3.10)  # Specify the minimum CMake version required

# Set the project name and version
project(MyProject)


set(ARCH "x64")  # 64-bit Windows
set(CMAKE_BUILD_TYPE Release)
# Add an executable target
#add_library(MyProject SHARED vec2D.cpp DynamicAnalysis.cpp)

#add_executable(MyProject main.cpp vec2D.cpp DynamicAnalysis_new.cpp  )
add_executable(MyProject main.cpp Frame.cpp structure.cpp)

