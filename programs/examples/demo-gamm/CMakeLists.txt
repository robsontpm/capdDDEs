add_executable(demo-gamm)

target_sources(demo-gamm
  PRIVATE
    "numeric.cpp"
    "setup.cpp"
    "validated.cpp"
)

target_link_libraries(demo-gamm
  PRIVATE
    capdDDEs::capdDDEs
)
