cmake_minimum_required(VERSION 3.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

project(DashMiniMD VERSION 1.0
        DESCRIPTION "Dash implementation of miniMD"
        LANGUAGES CXX)

add_executable(minimd main.cpp)

target_include_directories(minimd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_compile_options(minimd PRIVATE -Wall)

