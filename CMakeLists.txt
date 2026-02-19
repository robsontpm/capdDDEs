cmake_minimum_required(VERSION 3.13)
project(capdDDEs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_COVERAGE "Enable code coverage" OFF)
if(ENABLE_COVERAGE)
    message(STATUS "Enabling code coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find CAPD
find_package(CAPD REQUIRED)
message(STATUS "CAPD Found: ${CAPD_LIBRARIES}")
message(STATUS "CAPD Includes: ${CAPD_INCLUDE_DIRS}")

# --- Build capdDDEs library ---

# Read source list
file(READ "src/libcapdddes.txt" LIB_SOURCES_RAW)
string(REPLACE " " ";" LIB_SOURCES_LIST "${LIB_SOURCES_RAW}")
string(REPLACE "\n" ";" LIB_SOURCES_LIST "${LIB_SOURCES_LIST}")

# Read VPATH
file(READ "src/libvpath.txt" VPATH_RAW)
string(REPLACE " " ";" VPATH_LIST "${VPATH_RAW}")
string(REPLACE "\n" ";" VPATH_LIST "${VPATH_LIST}")

# Resolve source paths
set(LIB_SOURCES "")
foreach(SRC_NAME ${LIB_SOURCES_LIST})
    string(STRIP "${SRC_NAME}" SRC_NAME)
    if(NOT SRC_NAME)
        continue()
    endif()

    set(FOUND_SRC FALSE)
    foreach(PATH ${VPATH_LIST})
        string(STRIP "${PATH}" PATH)
        if(EXISTS "${CMAKE_SOURCE_DIR}/src/${PATH}/${SRC_NAME}.cpp")
            list(APPEND LIB_SOURCES "src/${PATH}/${SRC_NAME}.cpp")
            set(FOUND_SRC TRUE)
            break()
        endif()
    endforeach()

    if(NOT FOUND_SRC)
        message(FATAL_ERROR "Could not find source file: ${SRC_NAME}.cpp in VPATH: ${VPATH_LIST}")
    endif()
endforeach()

add_library(capdDDEs STATIC ${LIB_SOURCES})

target_include_directories(capdDDEs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${CAPD_INCLUDE_DIRS}
)

if(CAPD_COMPILE_OPTIONS)
    target_compile_options(capdDDEs PUBLIC ${CAPD_COMPILE_OPTIONS})
endif()
if(CAPD_COMPILE_DEFINITIONS)
    target_compile_definitions(capdDDEs PUBLIC ${CAPD_COMPILE_DEFINITIONS})
endif()

target_link_libraries(capdDDEs PUBLIC ${CAPD_LIBRARIES})

# --- Subdirectories ---

# Helper function
include(AddSubdirectoryPrograms)

option(BUILD_PROGRAMS "Build programs" ON)

if(BUILD_PROGRAMS)
    add_subdirectory(programs/utils)
    add_subdirectory(programs/results)
    add_subdirectory(programs/examples)
endif()

find_package(Boost QUIET COMPONENTS unit_test_framework)
if(Boost_FOUND)
    enable_testing()
    add_subdirectory(tests)
else()
    message(WARNING "Boost unit_test_framework not found. Tests will be skipped.")
endif()
