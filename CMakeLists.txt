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

add_library(impute_methylation_lib
    src/impute/tsv_io.cpp
    src/impute/haplotype_target.cpp
    src/impute/beta_binomial.cpp
    src/impute/evaluate.cpp
)
target_include_directories(impute_methylation_lib PUBLIC include)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(impute_methylation_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(merge_bedmethyl src/main.cpp)
target_link_libraries(merge_bedmethyl PRIVATE merge_bedmethyl_lib impute_methylation_lib)

add_executable(impute_methylation src/impute_main.cpp)
target_link_libraries(impute_methylation PRIVATE impute_methylation_lib)

add_executable(evaluate src/evaluate_main.cpp)
target_link_libraries(evaluate PRIVATE impute_methylation_lib)

install(TARGETS merge_bedmethyl impute_methylation evaluate RUNTIME DESTINATION bin)

enable_testing()
add_subdirectory(tests)
