add_executable(user-program-template)

target_sources(user-program-template
  PRIVATE
    "program.cpp"
    "setup.cpp"
)

target_link_libraries(user-program-template
  PRIVATE
    capdDDEs::capdDDEs
)
