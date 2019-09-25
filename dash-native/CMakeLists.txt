cmake_minimum_required(VERSION 3.1)

project(DashMiniMD VERSION 1.0
        DESCRIPTION "Dash implementation of miniMD"
        LANGUAGES CXX)

add_executable(minimd main.cpp)

target_include_directories(minimd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)
