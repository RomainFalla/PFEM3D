add_subdirectory(mesh)

add_subdirectory(solver)

set(SRCS
main.cpp)

find_path(JSON_INCLUDE_DIRS NAMES "nlohmann/json.hpp")
message(STATUS "Found nlohmann json: " ${JSON_INCLUDE_DIRS})
if(NOT JSON_INCLUDE_DIRS)
    message(FATAL_ERROR "nlohmann json not found!")
endif()

add_executable(pfem ${SRCS})
target_include_directories(pfem SYSTEM
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${JSON_INCLUDE_DIRS})
target_link_libraries(pfem PUBLIC pfemSolver)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfem PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfem PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfem PRIVATE /W4 /WX /wd4251)
endif()
