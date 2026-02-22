add_executable(mackey-glass-stable-periodic)

target_sources(mackey-glass-stable-periodic
  PRIVATE
    "nonrig-find.cpp"
    "nonrig-setup.cpp"
    "rig-prove.cpp"
    "rig-setup.cpp"
)

target_link_libraries(mackey-glass-stable-periodic
  PRIVATE
    capdDDEs::capdDDEs
)
