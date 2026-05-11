cmake_minimum_required(VERSION 3.10)
project(minimiseurAgregats)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(includes)

set(SOURCES 
    src/main.cpp 
    src/fonctions.cpp 
    src/maths.cpp
    src/minimisation.cpp 
    src/placement.cpp 
    src/struct.cpp 
)

add_executable(minimiseur ${SOURCES})