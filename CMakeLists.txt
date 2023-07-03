cmake_minimum_required(VERSION 3.16)
project(ContextSV)

set(CMAKE_CXX_STANDARD 14)

# Set directories
set(INCL_DIR ./include)
set(SRC_DIR ./src)

#set(HTSLIB_INCLUDE_DIR "/home/perdomoj/software/lib/include")
#set(HTSLIB_LIBRARIES "/home/perdomoj/software/lib/lib")

#set(HTSLIB_INCLUDE_DIR "/home/perdomoj/miniconda3/envs/lrst_py39/include")
#set(HTSLIB_LIBRARIES "/home/perdomoj/miniconda3/envs/lrst_py39/lib")

set(HTSLIB_INCLUDE_DIR "/home/perdomoj/software/htslib/include")
set(HTSLIB_LIBRARIES "/home/perdomoj/software/htslib/lib")

# Specify project and external library directories
#include_directories(${INCL_DIR} /usr/local/include ${HTSLIB_INCLUDE_DIR})
#link_directories(usr/local/lib/ ${HTSLIB_LIBRARIES})

include_directories(${INCL_DIR} ${HTSLIB_INCLUDE_DIR})
link_directories(${HTSLIB_LIBRARIES} "/home/perdomoj//miniconda3/envs/contextsvenv/lib")

# Set source files
set(
        SOURCE_FILES
        ${INCL_DIR}/cli.h
        ${INCL_DIR}/integrative_caller.h
	${INCL_DIR}/common.h
        ${INCL_DIR}/cnv_caller.h
        ${INCL_DIR}/khmm.h
        ${INCL_DIR}/kc.h
        ${SRC_DIR}/cli.cpp
        ${SRC_DIR}/integrative_caller.cpp
	${SRC_DIR}/common.cpp
        ${SRC_DIR}/cnv_caller.cpp
        ${SRC_DIR}/khmm.cpp
        ${SRC_DIR}/kc.cpp
        ${SRC_DIR}/main.cpp)

# Add executable
add_executable(contextsv ${SOURCE_FILES})

target_include_directories(contextsv PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(contextsv PUBLIC hts)  # Link HTS

