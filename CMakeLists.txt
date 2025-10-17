cmake_minimum_required(VERSION 3.28)
project(dipole_FEM_Cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(dipole_FEM_Cpp main.cpp
        dipoleAntenna.cpp
        dipoleAntenna.h
        ImpedanceMatcher.cpp
        ImpedanceMatcher.h
        DisplayResults.cpp
        DisplayResults.h
)
