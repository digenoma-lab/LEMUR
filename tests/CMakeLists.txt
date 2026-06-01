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
