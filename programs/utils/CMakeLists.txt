# programs/utils/CMakeLists.txt

add_programs_from_subdir(benchmark)
add_programs_from_subdir(issue_bad_alloc)

include(CheckCXXSourceCompiles)

# Check for CAPD multiprecision support
set(CMAKE_REQUIRED_INCLUDES ${CAPD_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${CAPD_LIBRARIES})
# We construct flags string from list
string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CAPD_COMPILE_OPTIONS}")
foreach(DEF ${CAPD_COMPILE_DEFINITIONS})
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -D${DEF}")
endforeach()

# Ensure C++17
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++17")

check_cxx_source_compiles("
#include <capd/mpcapdlib.h>
int main() {
    capd::MpInterval a;
    return 0;
}
" CAPD_HAS_MULTIPRECISION)

if(CAPD_HAS_MULTIPRECISION)
    message(STATUS "CAPD multiprecision support detected. Enabling 'converter' programs.")
    add_programs_from_subdir(converter)
else()
    message(STATUS "CAPD multiprecision support not found (MpInterval check failed). Skipping 'converter' programs.")
endif()

# Tests require Boost
find_package(Boost QUIET COMPONENTS unit_test_framework)

if(Boost_FOUND)
    message(STATUS "Boost unit_test_framework found.")
else()
    message(STATUS "Boost unit_test_framework not found.")
endif()
