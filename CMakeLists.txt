cmake_minimum_required(VERSION 3.10)
project(iucount)

set(CMAKE_CXX_STANDARD 14)

add_executable(iucount main.cpp)

target_link_libraries(iucount hts)
