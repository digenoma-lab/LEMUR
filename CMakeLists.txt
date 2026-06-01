cmake_minimum_required(VERSION 3.14)
project(MergeBedMethyl VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(merge_bedmethyl_lib
    src/locus.cpp
    src/reader.cpp
    src/format.cpp
    src/merger.cpp
    src/options.cpp
)
target_include_directories(merge_bedmethyl_lib PUBLIC include)

add_executable(merge_bedmethyl src/main.cpp)
target_link_libraries(merge_bedmethyl PRIVATE merge_bedmethyl_lib)

install(TARGETS merge_bedmethyl RUNTIME DESTINATION bin)

enable_testing()
add_subdirectory(tests)
