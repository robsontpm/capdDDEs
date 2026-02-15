find_package(Boost REQUIRED)

# automatically scan this folder for all .cpp files
# TODO: rethink, if there should be separated subfolders? When there is a lot of tests...
file(GLOB TEST_SOURCES "*.cpp")
foreach(TEST_SRC ${TEST_SOURCES})
    get_filename_component(TEST_NAME "${TEST_SRC}" NAME_WE)

    add_executable(${TEST_NAME} "${TEST_SRC}")

    target_link_libraries(${TEST_NAME} PRIVATE capdDDEs)
    target_include_directories(${TEST_NAME} PRIVATE ${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/programs/utils) # In case they need headers from there

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
