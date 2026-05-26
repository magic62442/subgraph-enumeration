set(PLAN_SRC
       decomposition.cpp estimator.cpp optimizer.cpp
        all_td.cpp
        all_td.h)

add_library(plan SHARED ${PLAN_SRC})
target_link_libraries(plan relation)
