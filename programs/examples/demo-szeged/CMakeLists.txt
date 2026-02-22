add_executable(demo-szeged)

target_sources(demo-szeged
  PRIVATE
    "1_nonrig.cpp"
    "2_manipulations.cpp"
    "3_rigorous_timemap.cpp"
    "4_periodic_candidate.cpp"
    "5_rigorous_poincare.cpp"
    "setup.cpp"
)

target_link_libraries(demo-szeged
  PRIVATE
    capdDDEs::capdDDEs
)
