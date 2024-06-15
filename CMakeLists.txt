cmake_minimum_required(VERSION 3.14)

project(prefix-codes CXX)

if (MSVC OR CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL MSVC)
    add_compile_options(/W4)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall)
endif()

add_library(prefix-codes STATIC
    "src/main.cpp"
    "src/utils.cpp"
    "src/utils.h"
)

macro(add_program name)
    add_executable(${name} "src/${name}.cpp")
    target_link_libraries(${name} PRIVATE prefix-codes)
endmacro()

add_program(01-naive-huffman)
add_program(02-canonical-huffman)
add_program(03-fast-huffman)
add_program(04-fast-lengths)
add_program(05-dumb-length-limiting)
add_program(06-smart-length-limiting)
