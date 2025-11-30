set(RELATION_SRC
        candidate_space.cpp
        labeled_graph.cpp
        relation.cpp
        graph.cpp)

add_library(relation SHARED ${RELATION_SRC})
target_link_libraries(relation utility)