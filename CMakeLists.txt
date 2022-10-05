cmake_minimum_required(VERSION 3.21)
project(Task)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(Matrix)

add_executable(
        Task main.cpp
        Matrix/Matrix.cpp Matrix/Matrix.h
        Matrix/Column.cpp Matrix/Column.h
        Matrix/Row.cpp Matrix/Row.h
        Matrix/Pair.cpp Matrix/Pair.h)

target_link_libraries(Task Matrix)
