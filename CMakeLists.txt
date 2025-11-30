cmake_minimum_required(VERSION 3.14)
project(subgraph-enumeration)

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    set(HOMEBREW_PREFIX "/opt/homebrew")
    set(OpenMP_ROOT "${HOMEBREW_PREFIX}/opt/libomp")
endif()

if(APPLE)
    message(STATUS "Configuring for macOS")
    message(STATUS "Configuring OpenMP for macOS (Homebrew)")
    include_directories("${OpenMP_ROOT}/include")
    link_libraries("${OpenMP_ROOT}/lib/libomp.dylib")
    add_compile_options(-Xpreprocessor -fopenmp)
    set(CMAKE_OSX_ARCHITECTURES "arm64")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    # OpenMP 配置
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found (非macOS): ${OpenMP_CXX_FLAGS}")
        include_directories(${OpenMP_INCLUDE_DIRS})
        add_compile_options(${OpenMP_CXX_FLAGS})
        link_libraries(${OpenMP_CXX_LIBRARIES})
    endif()
endif()

# 发布版本的额外优化
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(
        ${PROJECT_SOURCE_DIR}/config
        ${PROJECT_SOURCE_DIR}/utility
        ${PROJECT_SOURCE_DIR}/relation
        ${PROJECT_SOURCE_DIR}/plan
        ${PROJECT_SOURCE_DIR}/matching
)

add_subdirectory(config)
add_subdirectory(relation)
add_subdirectory(plan)
add_subdirectory(matching)
add_subdirectory(utility)
add_subdirectory(utility/td)
add_subdirectory(executable)

include_directories(${PROJECT_SOURCE_DIR}/utility/td/glpk_lib/include)
link_directories(${PROJECT_SOURCE_DIR}/utility/td/glpk_lib/lib)

add_executable(graph_mining_system main.cpp)
target_link_libraries(graph_mining_system matching)