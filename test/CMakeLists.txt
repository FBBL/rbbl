cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# custom function for testing
function(my_add_test name source dependency1 dependency2 expected_response)

  add_executable(${name} 
  	"${source}")

  target_link_libraries(${name} 
  	PRIVATE ${dependency1} ${dependency2}
  	)

  add_test(test_${name} ${name} 
  	WORKING_DIRECTORY ${OUTPUT_TEST_DIR}
  	)

  set_tests_properties(test_${name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${expected_response}"
  	)

endfunction()


link_directories(${BINARY_DIR})

# test20_005
my_add_test(20_005 "${TEST_DIR}/test20_005.c" m rbbl "Test Passed")
# test22_005_bruteforce
my_add_test(22_005_bruteforce "${TEST_DIR}/test22_005_bruteforce.c" m rbbl "Test Passed")
# test40_005_bruteforce
my_add_test(40_005_bruteforce "${TEST_DIR}/test40_005_bruteforce.c" m rbbl "Test Passed")
