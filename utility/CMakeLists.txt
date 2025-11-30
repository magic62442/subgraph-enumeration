set(UTILITY_SRC
    command_parser.cpp utils.cpp command.cpp compute_set_intersection.cpp td/preprocessing.h clique.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/automorphism)

find_library(NAUTY_LIB nauty PATHS ${CMAKE_CURRENT_SOURCE_DIR}/automorphism)
if (NOT NAUTY_LIB)
    message(FATAL_ERROR "Could not find nauty library")
endif ()

#link_directories(${CMAKE_CURRENT_SOURCE_DIR}/automorphism)

add_library(utility SHARED  ${UTILITY_SRC})

target_link_libraries(utility PUBLIC td ${NAUTY_LIB})