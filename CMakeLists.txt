cmake_minimum_required(VERSION 3.15)
project(reshala VERSION 0.1.0 LANGUAGES CXX)

# ----------------------------------------
# Build options
# ----------------------------------------
option(BUILD_TESTS "Build test executable" ON)
option(BUILD_EXAMPLES "Build examples" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")

# ----------------------------------------
# Library
# ----------------------------------------

add_executable(${PROJECT_NAME}-cli ${EXE_SOURCES} src/main.cpp)
add_library(${PROJECT_NAME}-lib SHARED ${LIB_SOURCES})
add_subdirectory(src)

target_link_libraries(${PROJECT_NAME}-cli PUBLIC ${PROJECT_NAME}-lib )
target_include_directories(${PROJECT_NAME}-lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

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
