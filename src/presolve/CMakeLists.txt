target_sources(${PROJECT_NAME}-lib PRIVATE
    bounder.cpp
    model_tracker.cpp
    presolve.cpp
    transforms.cpp
)

add_subdirectory(rules)
