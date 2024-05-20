cmake_minimum_required(VERSION 3.14)

project(prefix-codes CXX)

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
