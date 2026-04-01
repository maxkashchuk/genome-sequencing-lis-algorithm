cmake_minimum_required(VERSION 3.11)

project(fasta_lis_sequencing VERSION 3.1.0
                  LANGUAGES CXX
                  DESCRIPTION "Sequencing of FASTA data using LIS algorithm")

add_subdirectory(bioparser)

add_executable(fasta_lis_sequencing src/main.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include")

target_link_libraries(${PROJECT_NAME} PRIVATE bioparser)