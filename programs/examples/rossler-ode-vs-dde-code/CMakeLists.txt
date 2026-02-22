add_executable(rossler-ode-vs-dde-code)

target_sources(command-line-helper
  PRIVATE
    "capd-original-code.cpp"
    "dde-vs-ode-code.cpp"
    "gnuplot.cpp"
    "setup.cpp"
)

target_link_libraries(command-line-helper
  PRIVATE
    capdDDEs::capdDDEs
)
