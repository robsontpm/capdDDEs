# programs/utils/CMakeLists.txt

add_programs_from_subdir(benchmark)
add_programs_from_subdir(converter)

# Tests require Boost
find_package(Boost QUIET COMPONENTS unit_test_framework)

if(Boost_FOUND)
    message(STATUS "Boost unit_test_framework found. Enabling tests.")
#    add_programs_from_subdir(tests LIBS Boost::unit_test_framework)
else()
    message(STATUS "Boost unit_test_framework not found. Skipping tests.")
endif()
