target_sources(${PROJECT_NAME}-lib PRIVATE
    bnb.cpp
    milp.cpp
)

add_subdirectory(branching)
