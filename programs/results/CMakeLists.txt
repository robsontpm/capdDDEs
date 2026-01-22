# programs/results/CMakeLists.txt

# Iterate over subdirectories and add programs
file(GLOB children RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
foreach(child ${children})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${child})
        add_programs_from_subdir(${child})
    endif()
endforeach()
