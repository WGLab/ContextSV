cmake_minimum_required(VERSION 3.16)
project(ContextSV)

set(CMAKE_CXX_STANDARD 14)

# Set directories
set(INCL_DIR ./include)
set(SRC_DIR ./src)

# Specify project and external library directories
include_directories(${INCL_DIR} /usr/local/include)
link_directories(usr/local/lib)

# Set source files
set(
        SOURCE_FILES
        ${INCL_DIR}/cli.h
        ${INCL_DIR}/bam_reader.h
        ${SRC_DIR}/main.cpp
        ${SRC_DIR}/cli.cpp
        ${SRC_DIR}/bam_reader.cpp)

add_executable(contextsv ${SOURCE_FILES})

target_link_libraries(contextsv PUBLIC hts)  # Link HTS
# For Linux development: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:HTSInstallPath
