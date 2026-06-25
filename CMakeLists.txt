cmake_minimum_required(VERSION 3.11)
project(nazgul)



set(CMAKE_BUILD_TYPE Release)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)

add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE libnazgul)