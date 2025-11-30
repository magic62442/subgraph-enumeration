set(MATCH_SRC
        join.cpp parallel.cpp)

add_library(matching SHARED ${MATCH_SRC})
target_link_libraries(matching plan)

# Link OpenMP
if(APPLE)
    # macOS 使用 homebrew 的 libomp
    target_link_libraries(matching "/opt/homebrew/opt/libomp/lib/libomp.dylib")
else()
    if(OpenMP_CXX_FOUND)
        target_link_libraries(matching OpenMP::OpenMP_CXX)
    endif()
endif()