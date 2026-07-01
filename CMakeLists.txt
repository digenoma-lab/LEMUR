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

add_library(dml_lib
    src/dml/metadata.cpp
    src/dml/model.cpp
    src/dml/fdr.cpp
    src/dml/tsv_stream.cpp
    src/dml/util.cpp
    src/dml/call_dmr.cpp
)
target_include_directories(dml_lib PUBLIC include)

add_executable(dml src/dml_main.cpp)
target_link_libraries(dml PRIVATE dml_lib)

add_executable(call_dmr src/call_dmr_main.cpp)
target_link_libraries(call_dmr PRIVATE dml_lib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(impute_methylation_lib PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(dml_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(merge_bedmethyl src/main.cpp)
target_link_libraries(merge_bedmethyl PRIVATE merge_bedmethyl_lib impute_methylation_lib)

add_executable(impute_methylation src/impute_main.cpp)
target_link_libraries(impute_methylation PRIVATE impute_methylation_lib)

add_executable(smooth_methylation src/smooth_main.cpp)
target_link_libraries(smooth_methylation PRIVATE impute_methylation_lib)

add_executable(evaluate src/evaluate_main.cpp)
target_link_libraries(evaluate PRIVATE impute_methylation_lib)

add_executable(tracts src/tracts_main.cpp)

add_executable(la-dml src/la_dml_main.cpp)
if(OpenMP_CXX_FOUND)
    target_link_libraries(la-dml PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS merge_bedmethyl impute_methylation smooth_methylation evaluate dml call_dmr tracts la-dml RUNTIME DESTINATION bin)

enable_testing()
add_subdirectory(tests)
