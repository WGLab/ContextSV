cmake_minimum_required(VERSION 3.16)
project(ContextSV)

set(CMAKE_CXX_STANDARD 14)

# Set project directories
set(INCL_DIR ./include)
set(SRC_DIR ./src)

include_directories(${INCL_DIR})

# Set HTSLib directories
include_directories(${PROJECT_SOURCE_DIR}/include/htslib)
# include_directories(${PROJECT_SOURCE_DIR}/lib)

link_directories(${PROJECT_SOURCE_DIR}/lib)

# Set the OpenSSL directories
include_directories(${PROJECT_SOURCE_DIR}/include/openssl)
# set(HTSLIB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/htslib)
# set(HTSLIB_LIBRARIES ${CMAKE_CURRENT_SOURCE_DIR}/lib/htslib)

# Set include and link directories
# link_directories(${HTSLIB_LIBRARIES})
# include_directories(${INCL_DIR} ${HTSLIB_INCLUDE_DIR})
# link_directories(${HTSLIB_LIBRARIES})


# Set source files
set(
        SOURCE_FILES
        ${INCL_DIR}/contextsv.h
        ${SRC_DIR}/contextsv.cpp
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

# Link HTSLib
# target_include_directories(contextsv PUBLIC ${HTSLIB_INCLUDE_DIR})
# target_link_libraries(contextsv ${HTSLIB_LIBRARIES})

target_link_libraries(contextsv hts)
