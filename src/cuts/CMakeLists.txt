add_subdirectory(generators)

target_sources(${PROJECT_NAME}-lib PRIVATE
    cut.cpp
    cutter.cpp
)
