cmake_minimum_required(VERSION 3.16)
project(ContextSV)

set(CMAKE_CXX_STANDARD 14)

# Set directories
set(INCL_DIR ./include)
set(SRC_DIR ./src)

# Find HTS
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTS REQUIRED htslib)

# Set HTS directories
set(HTSLIB_INCLUDE_DIR ${HTS_INCLUDE_DIRS})
set(HTSLIB_LIBRARIES ${HTS_LIBRARIES})

# set(HTSLIB_INCLUDE_DIR "/home/perdomoj/software/htslib/include")
# set(HTSLIB_LIBRARIES "/home/perdomoj/software/htslib/lib")

include_directories(${INCL_DIR} ${HTSLIB_INCLUDE_DIR})
link_directories(${HTSLIB_LIBRARIES} "/home/perdomoj//miniconda3/envs/contextsvenv/lib")

# Set source files
set(
        SOURCE_FILES
        ${INCL_DIR}/cli.h
        ${INCL_DIR}/integrative_caller.h
	${INCL_DIR}/common.h
	${INCL_DIR}/common.h
        ${INCL_DIR}/cnv_caller.h
        ${INCL_DIR}/khmm.h
        ${INCL_DIR}/kc.h
        ${SRC_DIR}/cli.cpp
        ${SRC_DIR}/integrative_caller.cpp
	${SRC_DIR}/common.cpp
	${SRC_DIR}/common.cpp
        ${SRC_DIR}/cnv_caller.cpp
        ${SRC_DIR}/khmm.cpp
        ${SRC_DIR}/kc.cpp
        ${SRC_DIR}/main.cpp)

# Add executable
add_executable(contextsv ${SOURCE_FILES})
target_include_directories(contextsv PRIVATE ${HTSLIB_INCLUDE_DIR})
target_link_libraries(contextsv PUBLIC hts)  # Link HTS

