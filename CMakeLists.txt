cmake_minimum_required(VERSION 3.15)
project(reshala VERSION 0.1.0 LANGUAGES CXX)

# ----------------------------------------
# Build options
# ----------------------------------------
option(BUILD_TESTS "Build unit tests" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")

# Optional: unit tests with GoogleTest
if(BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    FetchContent_MakeAvailable(googletest)
endif()

# ----------------------------------------
# Library
# ----------------------------------------
set(LIB_SOURCES
    src/io/lp/lp_reader.cpp
    src/linalg/dense_vector.cpp
    src/linalg/sparse_vector.cpp
    src/linalg/sparse_matrix.cpp
    src/linalg/operators.cpp
    src/model/milp_model.cpp
    src/lp/dual_simplex.cpp
)
set(EXE_SOURCES
    src/main.cpp
)

add_library(${PROJECT_NAME}-lib SHARED ${LIB_SOURCES})
add_executable(${PROJECT_NAME}-cli ${EXE_SOURCES})
target_include_directories(${PROJECT_NAME}-lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(${PROJECT_NAME}-cli PUBLIC ${PROJECT_NAME}-lib )
#target_compile_definitions(milp_solver PRIVATE -DMILP_BUILD)

# ----------------------------------------
# Examples
# ----------------------------------------
if(BUILD_EXAMPLES)
    add_executable(simple_mip examples/simple_mip.cpp)
    target_link_libraries(simple_mip milp_solver)
endif()

# ----------------------------------------
# Tests
# ----------------------------------------
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ----------------------------------------
# Installation (optional)
# ----------------------------------------
install(TARGETS ${PROJECT_NAME}-cli
    EXPORT reshalaTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(DIRECTORY include/ DESTINATION include)
