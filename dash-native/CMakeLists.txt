# Note: simply set CXX to the DASH compiler wrapper

cmake_minimum_required(VERSION 3.9)

project(DashMiniMD VERSION 1.0
        DESCRIPTION "Dash implementation of miniMD"
        LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(DASH_ROOT $ENV{HOME}/opt/dash-0.3.0-dev)
set(MEPHISTO_ROOT $ENV{HOME}/src/mephisto)

find_package(Boost REQUIRED)
find_package(MPI REQUIRED)
find_package(Threads REQUIRED)


add_executable(minimd main.cpp)

# target_compile_definitions(minimd PUBLIC DASH_RELEASE DASH_HAVE_STD_TRIVIALLY_COPYABLE DASH_ENV_HOST_SYSTEM_ID=default DASH_MPI_IMPL_ID=mpich)
target_include_directories(minimd PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
        )
target_include_directories(minimd SYSTEM PRIVATE
        ${DASH_ROOT}/include
        ${MEPHISTO_ROOT}/alpaka/include
        ${DASH_ROOT}/patterns/include
        ${DASH_ROOT}/dash/include
        ${Boost_INCLUDE_DIRS}
        )

target_link_libraries(minimd dash-mpi dart-mpi dart-base MPI::MPI_CXX Threads::Threads)
target_link_directories(minimd PUBLIC ${DASH_ROOT}/lib)
target_compile_options(minimd PRIVATE -Wall)
