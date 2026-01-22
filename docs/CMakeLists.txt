cmake_minimum_required(VERSION 3.13)
project(MyCapdDDEsProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set path to capdDDEs root directory
# Assuming this CMakeLists.txt is used in a structure where capdDDEs is available.
# Change this path to point to your capdDDEs checkout.
set(CAPDDDES_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
get_filename_component(CAPDDDES_ROOT "${CAPDDDES_ROOT}" ABSOLUTE)

# Add cmake module path to find CAPD using the script provided in capdDDEs
list(APPEND CMAKE_MODULE_PATH "${CAPDDDES_ROOT}/cmake")

# Find CAPD
find_package(CAPD REQUIRED)

# Include directories
include_directories("${CAPDDDES_ROOT}/include")
include_directories(${CAPD_INCLUDE_DIRS})

# Library path
# Assuming build directory is "build" inside root.
# If you used a different build directory, update this path.
link_directories("${CAPDDDES_ROOT}/build")

# Create executable
# Replace my_program.cpp with your source file
# add_executable(my_program my_program.cpp)

# Link libraries
# target_link_libraries(my_program capdDDEs ${CAPD_LIBRARIES})

message(STATUS "Configured to use capdDDEs from ${CAPDDDES_ROOT}")
message(STATUS "Please uncomment add_executable and target_link_libraries lines to build your program.")
