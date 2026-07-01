set(TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)

add_test(
    NAME merge_two_samples
    COMMAND merge_bedmethyl -c 3 -s 2 ${CMAKE_CURRENT_BINARY_DIR}/test_out.tsv
            S1 ${TEST_DATA_DIR}/S1_hp1.bedmethyl ${TEST_DATA_DIR}/S1_hp2.bedmethyl
            S2 ${TEST_DATA_DIR}/S2_hp1.bedmethyl ${TEST_DATA_DIR}/S2_hp2.bedmethyl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

set_tests_properties(merge_two_samples PROPERTIES
    PASS_REGULAR_EXPRESSION "Wrote 1 rows"
)

add_test(
    NAME merge_output_format
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${CMAKE_CURRENT_BINARY_DIR}/test_out.tsv
        ${CMAKE_CURRENT_SOURCE_DIR}/expected/merge_two_samples.tsv
)
set_tests_properties(merge_output_format PROPERTIES DEPENDS merge_two_samples)

add_test(
    NAME merge_two_samples_sample_mode
    COMMAND merge_bedmethyl --sample -c 3 -s 2 ${CMAKE_CURRENT_BINARY_DIR}/test_out_sample.tsv
            S1 ${TEST_DATA_DIR}/S1_combined.bedmethyl
            S2 ${TEST_DATA_DIR}/S2_combined.bedmethyl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(merge_two_samples_sample_mode PROPERTIES
    PASS_REGULAR_EXPRESSION "Wrote 1 rows")

add_test(
    NAME merge_output_format_sample_mode
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${CMAKE_CURRENT_BINARY_DIR}/test_out_sample.tsv
        ${CMAKE_CURRENT_SOURCE_DIR}/expected/merge_two_samples_sample_mode.tsv
)
set_tests_properties(merge_output_format_sample_mode PROPERTIES
    DEPENDS merge_two_samples_sample_mode)

add_test(
    NAME impute_stream_tiny
    COMMAND impute_methylation -a 1 -b 1
            ${TEST_DATA_DIR}/tiny.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_imputed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME impute_output_columns
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_imputed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_hap1_imputed.tsv
)
set_tests_properties(impute_output_columns PROPERTIES DEPENDS impute_stream_tiny)

add_test(
    NAME impute_stream_tiny_sample
    COMMAND impute_methylation --sample -a 1 -b 1
            ${TEST_DATA_DIR}/tiny_sample.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_imputed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME impute_output_columns_sample
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_imputed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_sample_imputed.tsv
)
set_tests_properties(impute_output_columns_sample PROPERTIES DEPENDS impute_stream_tiny_sample)

add_test(
    NAME impute_stream_tiny_sample_counts_cov
    COMMAND impute_methylation --sample --counts-cov -a 1 -b 1
            ${TEST_DATA_DIR}/tiny_sample.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_counts_cov_imputed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME impute_output_columns_sample_counts_cov
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_counts_cov_imputed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_sample_counts_cov_imputed.tsv
)
set_tests_properties(impute_output_columns_sample_counts_cov PROPERTIES
    DEPENDS impute_stream_tiny_sample_counts_cov)

add_test(
    NAME smooth_stream_tiny
    COMMAND smooth_methylation -a 1 -b 1
            ${TEST_DATA_DIR}/tiny.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_smoothed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME smooth_output_columns
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_smoothed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_hap1_smoothed.tsv
)
set_tests_properties(smooth_output_columns PROPERTIES DEPENDS smooth_stream_tiny)

add_test(
    NAME smooth_stream_tiny_sample
    COMMAND smooth_methylation --sample -a 1 -b 1
            ${TEST_DATA_DIR}/tiny_sample.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_smoothed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME smooth_output_columns_sample
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_smoothed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_sample_smoothed.tsv
)
set_tests_properties(smooth_output_columns_sample PROPERTIES DEPENDS smooth_stream_tiny_sample)

add_test(
    NAME smooth_stream_tiny_sample_counts_cov
    COMMAND smooth_methylation --sample --counts-cov -a 1 -b 1
            ${TEST_DATA_DIR}/tiny_sample.tsv
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_counts_cov_smoothed.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(
    NAME smooth_output_columns_sample_counts_cov
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample_counts_cov_smoothed.tsv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/tiny_sample_counts_cov_smoothed.tsv
)
set_tests_properties(smooth_output_columns_sample_counts_cov PROPERTIES
    DEPENDS smooth_stream_tiny_sample_counts_cov)

add_test(
    NAME evaluate_tiny
    COMMAND evaluate -c S1.hap1_counts -m 0.2 -s 42 -w 200 -n 5
            ${TEST_DATA_DIR}/tiny.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(evaluate_tiny PROPERTIES
    PASS_REGULAR_EXPRESSION "S1\\.hap1_frac_imputed.*MSE="
)

add_test(
    NAME evaluate_cohort_tiny
    COMMAND evaluate -chr chr1 -m 0.2 -s 42 -w 200 -n 1 -j 2
            -o ${CMAKE_CURRENT_BINARY_DIR}/tiny.cohort.eval.tsv
            ${TEST_DATA_DIR}/tiny.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(evaluate_cohort_tiny PROPERTIES
    PASS_REGULAR_EXPRESSION "Evaluated 2 haplotype columns"
)

add_test(
    NAME evaluate_tiny_sample
    COMMAND evaluate --sample -c S1.counts -m 0.2 -s 42 -w 200 -n 5
            ${TEST_DATA_DIR}/tiny_sample.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(evaluate_tiny_sample PROPERTIES
    PASS_REGULAR_EXPRESSION "S1\\.frac_imputed.*MSE="
)

add_test(
    NAME evaluate_cohort_tiny_sample
    COMMAND evaluate --sample -chr chr1 -m 0.2 -s 42 -w 200 -n 1 -j 2
            -o ${CMAKE_CURRENT_BINARY_DIR}/tiny_sample.cohort.eval.tsv
            ${TEST_DATA_DIR}/tiny_sample.tsv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(evaluate_cohort_tiny_sample PROPERTIES
    PASS_REGULAR_EXPRESSION "Evaluated 1 sample columns"
)

add_test(
    NAME dml_tiny_cohort
    COMMAND dml --sample -j 1
            ${TEST_DATA_DIR}/dml_cohort.tsv
            ${TEST_DATA_DIR}/dml_metadata.csv
            ${CMAKE_CURRENT_BINARY_DIR}/dml_cohort.out.csv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(dml_tiny_cohort PROPERTIES
    PASS_REGULAR_EXPRESSION "tested 5 CpGs"
)

add_test(
    NAME dml_output_format
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/dml_cohort.out.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/dml_cohort.csv
)
set_tests_properties(dml_output_format PROPERTIES DEPENDS dml_tiny_cohort)

add_test(
    NAME call_dmr_tiny
    COMMAND call_dmr --sample -j 1
            ${TEST_DATA_DIR}/dml_call_dmr.csv
            ${CMAKE_CURRENT_BINARY_DIR}/call_dmr.out.csv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_tests_properties(call_dmr_tiny PROPERTIES
    PASS_REGULAR_EXPRESSION "DMRs found: 1"
)

add_test(
    NAME call_dmr_output_format
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/call_dmr.out.csv
            ${CMAKE_CURRENT_SOURCE_DIR}/expected/call_dmr.csv
)
set_tests_properties(call_dmr_output_format PROPERTIES DEPENDS call_dmr_tiny)
