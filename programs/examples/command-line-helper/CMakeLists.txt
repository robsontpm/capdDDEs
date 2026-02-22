add_executable(command-line-helper)

target_sources(command-line-helper
  PRIVATE
    "cmd-parser.cpp"
)

target_link_libraries(command-line-helper
  PRIVATE
    capdDDEs::capdDDEs
)
