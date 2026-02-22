find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include(BoostTestDiscoverTests)

function(add_simple_boost_test target_name)
  add_executable(${target_name})
  target_sources(${target_name}
    PRIVATE
      "${target_name}.cpp")
  target_link_libraries(${target_name}
    PRIVATE
      # Boost::unit_test_framework
      capdDDEs::capdDDEs
  )
  target_include_directories(${target_name}
    PRIVATE
      ${Boost_INCLUDE_DIRS}
  )
  boosttest_discover_tests(${target_name})
endfunction()

add_simple_boost_test(capd-ddes-DDECommon)
add_simple_boost_test(capd-ddes-storage-BasicDoubleton)
add_simple_boost_test(cubeset_tests)
add_simple_boost_test(legacy_tests)
