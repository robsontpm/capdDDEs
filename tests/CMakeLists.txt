find_package(Boost REQUIRED COMPONENTS unit_test_framework)

file(GLOB TEST_SOURCES "*.cpp")

foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME "${TEST_SRC}" NAME_WE)

    # Check if file has "main" to avoid conflict with BOOST_TEST_MODULE which provides main
    # Actually, Boost Test sources usually define BOOST_TEST_MODULE and include unit_test.hpp,
    # which generates main.
    # If we are refactoring existing files that HAD main, we will remove it.

    add_executable(${TEST_NAME} "${TEST_SRC}")

    target_link_libraries(${TEST_NAME} PRIVATE capdDDEs Boost::unit_test_framework)
    target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/programs/utils) # In case they need headers from there

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
