add_subdirectory(mesh)

add_subdirectory(simulation)

set(SRCS
main.cpp
sharedLib_defines.h)

add_executable(pfem ${SRCS})
target_include_directories(pfem SYSTEM
                           PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${JSON_INCLUDE_DIRS})
target_link_libraries(pfem PUBLIC pfemSimulation)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    target_compile_options(pfem PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES CLANG)
    target_compile_options(pfem PRIVATE -Wall -Wextra -pedantic-errors -Wold-style-cast -Wnull-dereference -Wshadow)
elseif(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    target_compile_options(pfem PRIVATE /W4 /WX /wd4251)
endif()
