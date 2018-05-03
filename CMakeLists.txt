cmake_minimum_required(VERSION 3.9)
project(LKH_2 C)

set(CMAKE_C_STANDARD 99)

#将第一个目录参数中包含的所有源文件名赋值给第二个变量
AUX_SOURCE_DIRECTORY(SRC SRC_DIR)

add_executable(LKH_2 ${SRC_DIR})
include_directories(SRC/INCLUDE)
#add_executable(TEST test.c)