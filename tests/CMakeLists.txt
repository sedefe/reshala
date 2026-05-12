add_executable(${PROJECT_NAME}-test run_tests.cpp)

target_link_libraries(${PROJECT_NAME}-test PRIVATE ${PROJECT_NAME}-lib)

target_include_directories(${PROJECT_NAME}-test
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

enable_testing()
add_test(NAME ${PROJECT_NAME}-test 
    COMMAND ${PROJECT_NAME}-test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
