cmake_minimum_required(VERSION 3.6)

project(hamming CXX)

find_package(OpenMP)

add_executable (hamming src/new_hamming.cpp)
add_executable (fragment_hamming src/fragment_hamming.cpp)
add_executable (homology src/homology.cpp)

target_link_libraries(fragment_hamming OpenMP::OpenMP_CXX)
target_link_libraries(hamming OpenMP::OpenMP_CXX)
target_link_libraries(homology OpenMP::OpenMP_CXX)
