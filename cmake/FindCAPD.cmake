# FindCAPD.cmake

# Try to find capd-config
# Prioritize bundled CAPD build directory
find_program(CAPD_CONFIG_EXECUTABLE
    NAMES capd-config
    PATHS
        ${CMAKE_SOURCE_DIR}/bin/capd_build/bin
        ${CAPD_DIR}/bin
    NO_DEFAULT_PATH
    DOC "Path to bundled or specified capd-config executable"
)

# Fallback to system paths only if not found above and not strictly forbidden
if(NOT CAPD_CONFIG_EXECUTABLE)
    find_program(CAPD_CONFIG_EXECUTABLE
        NAMES capd-config
        PATHS
            /usr/local/bin
            /usr/bin
        DOC "Path to system capd-config executable"
    )
endif()

if(CAPD_CONFIG_EXECUTABLE)
    message(STATUS "Using CAPD from: ${CAPD_CONFIG_EXECUTABLE}")

    execute_process(
        COMMAND ${CAPD_CONFIG_EXECUTABLE} --cflags
        OUTPUT_VARIABLE CAPD_CFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${CAPD_CONFIG_EXECUTABLE} --libs
        OUTPUT_VARIABLE CAPD_LDFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Process CFLAGS
    # We use separate_arguments to handle spaces correctly
    separate_arguments(CAPD_CFLAGS_LIST UNIX_COMMAND "${CAPD_CFLAGS}")

    set(CAPD_INCLUDE_DIRS "")
    set(CAPD_COMPILE_OPTIONS "")
    set(CAPD_COMPILE_DEFINITIONS "")

    foreach(ARG ${CAPD_CFLAGS_LIST})
        if(ARG MATCHES "^-I(.+)")
            list(APPEND CAPD_INCLUDE_DIRS "${CMAKE_MATCH_1}")
        elseif(ARG MATCHES "^-D(.+)")
            list(APPEND CAPD_COMPILE_DEFINITIONS "${CMAKE_MATCH_1}")
        else()
            list(APPEND CAPD_COMPILE_OPTIONS "${ARG}")
        endif()
    endforeach()

    # Process LDFLAGS
    separate_arguments(CAPD_LIBRARIES UNIX_COMMAND "${CAPD_LDFLAGS}")

endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CAPD
    REQUIRED_VARS CAPD_CONFIG_EXECUTABLE CAPD_LIBRARIES
    FOUND_VAR CAPD_FOUND
)

mark_as_advanced(CAPD_CONFIG_EXECUTABLE)
