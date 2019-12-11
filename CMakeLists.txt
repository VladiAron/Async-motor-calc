cmake_minimum_required(VERSION 3.14)
project(motorCalc C)

set(CMAKE_C_STANDARD 99)

add_executable(motorCalc main.c)

target_link_libraries(motorCalc m)
#target_link_options(motorCalc PUBLIC -lm)