# project name
project(phoxtrot)

# minimum cmake version
cmake_minimum_required(VERSION 2.8.11)

include_directories("${PROJECT_SOURCE_DIR}/src")

add_executable(phoxtrot src/main.cpp)
