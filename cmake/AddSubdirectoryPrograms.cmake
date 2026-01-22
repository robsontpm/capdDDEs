# cmake/AddSubdirectoryPrograms.cmake

function(add_programs_from_subdir SUBDIR)
    # Parse extra arguments (libraries to link)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs LIBS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Get absolute path to the subdirectory
    get_filename_component(SUBDIR_ABS "${SUBDIR}" ABSOLUTE)

    # Check if directory exists
    if(NOT IS_DIRECTORY "${SUBDIR_ABS}")
        message(STATUS "Directory ${SUBDIR} does not exist or is not a directory. Skipping.")
        return()
    endif()

    # Find all .cpp files
    file(GLOB CPP_FILES "${SUBDIR_ABS}/*.cpp")

    set(SOURCES "")
    set(EXECUTABLES "")

    foreach(FILE ${CPP_FILES})
        file(READ "${FILE}" FILE_CONTENT)
        # Simple heuristic to check for main()
        # Checks for " int main" or " void main" or "\nmain(" or " main("
        # We look for standard signatures.
        if("${FILE_CONTENT}" MATCHES "[ \t\n](int|void)[ \t\n]+main[ \t\n]*\\(")
            list(APPEND EXECUTABLES "${FILE}")
        else()
            list(APPEND SOURCES "${FILE}")
        endif()
    endforeach()

    if(NOT EXECUTABLES)
        return()
    endif()

    foreach(SOURCE ${SOURCES})
        message(STATUS "Using extra source file (no main): ${SOURCE}")
    endforeach()    

    # Process executables
    foreach(EXE_SOURCE ${EXECUTABLES})
        message(STATUS "Adding entry for executable: ${EXE_SOURCE}")
        get_filename_component(EXE_NAME "${EXE_SOURCE}" NAME_WE)

        # Create a unique target name by sanitizing the subdirectory path
        file(RELATIVE_PATH REL_SUBDIR "${CMAKE_SOURCE_DIR}" "${SUBDIR_ABS}")
        string(REPLACE "/" "_" SUBDIR_SANITIZED "${REL_SUBDIR}")
        set(TARGET_NAME "${SUBDIR_SANITIZED}_${EXE_NAME}")

        # Add executable
        add_executable(${TARGET_NAME} "${EXE_SOURCE}" ${SOURCES})

        # Set the output name to match the file name
        set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME "${EXE_NAME}")

        # Link libraries
        target_link_libraries(${TARGET_NAME} PRIVATE capdDDEs ${CAPD_LIBRARIES} ${ARG_LIBS})

        # Add compile options if any from CAPD
        if(CAPD_COMPILE_OPTIONS)
            target_compile_options(${TARGET_NAME} PRIVATE ${CAPD_COMPILE_OPTIONS})
        endif()
        if(CAPD_COMPILE_DEFINITIONS)
            target_compile_definitions(${TARGET_NAME} PRIVATE ${CAPD_COMPILE_DEFINITIONS})
        endif()

        # Add include directories (current subdirectory)
        target_include_directories(${TARGET_NAME} PRIVATE "${SUBDIR_ABS}")

        # Set output directory to match the original layout: ./bin/ inside the source subdirectory
        set_target_properties(${TARGET_NAME} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        )
    endforeach()
endfunction()
