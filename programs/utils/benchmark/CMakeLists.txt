add_executable(benchmark)

target_sources(benchmark
  PRIVATE
    benchmark.cpp
)

target_link_libraries(benchmark
  PRIVATE
    capdDDEs::capdDDEs
)

add_executable(benchmark2)

target_sources(benchmark2
  PRIVATE
    benchmark2.cpp
)

target_link_libraries(benchmark2
  PRIVATE
    capdDDEs::capdDDEs
)
