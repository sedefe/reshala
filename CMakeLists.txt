cmake_minimum_required(VERSION 3.15)
project(reshala VERSION 0.1.0 LANGUAGES CXX)

# ----------------------------------------
# Build options
# ----------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer -DDEBUG")

option(ENABLE_SANITIZERS "Enable Address and UndefinedBehavior sanitizers" ON)
if(ENABLE_SANITIZERS AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Enabling sanitizers (ASan + UBSan)")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
endif()

# ----------------------------------------
# Library
# ----------------------------------------

add_executable(${PROJECT_NAME}-cli src/main.cpp)
add_library(${PROJECT_NAME}-lib SHARED)
add_subdirectory(src)

target_link_libraries(${PROJECT_NAME}-cli PRIVATE ${PROJECT_NAME}-lib)
target_include_directories(${PROJECT_NAME}-lib
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# ----------------------------------------
# Tests
# ----------------------------------------
option(BUILD_TESTS "Build test executable" ON)
if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

# ----------------------------------------
# Installation (optional)
# ----------------------------------------
install(TARGETS ${PROJECT_NAME}-cli ${PROJECT_NAME}-lib
    EXPORT reshalaTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(DIRECTORY include/ DESTINATION include)
