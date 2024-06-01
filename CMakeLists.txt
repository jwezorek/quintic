cmake_minimum_required(VERSION 3.16)

project(quintic VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.80 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

add_executable(quintic
    src/laguerre.cpp
    src/main.cpp
    src/quintic.hpp
)