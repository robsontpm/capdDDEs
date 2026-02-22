function(add_simple_converter target_name)
  add_executable(${target_name})
  target_sources(${target_name}
    PRIVATE
      "${target_name}.cpp"
  )
  target_link_libraries(${target_name}
    PRIVATE
      capdDDEs::capdDDEs)
endfunction()

add_simple_converter(convmatrix)
add_simple_converter(convvector)
add_simple_converter(crmatrix)
add_simple_converter(gencoords)
add_simple_converter(growvector)
add_simple_converter(invmatrixtest)
add_simple_converter(matrixcmp)
add_simple_converter(midmatrix)
add_simple_converter(mulvector)
add_simple_converter(rmatrix)
add_simple_converter(splitvector)
add_simple_converter(statmatrix)
add_simple_converter(vectorcmp)
add_simple_converter(vectorhull)
