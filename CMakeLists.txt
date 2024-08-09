cmake_minimum_required(VERSION 3.11)

project(PyrTex)
set(PROJECT_VERSION 1.0.0)

include_directories("lib")

add_executable(
    PyrTex
    main.cpp)