add_executable(demo-elninio)

target_sources(demo-elninio
  PRIVATE
    "draw-nonrig.cpp"
    "setup-nonrig.cpp"
)

target_link_libraries(demo-elninio
  PRIVATE
    capdDDEs::capdDDEs
)
