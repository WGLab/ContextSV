cmake_minimum_required(VERSION 3.16)
project(ContextSV)

set(CMAKE_CXX_STANDARD 14)

# Set directories
set(INCL_DIR ./include)
set(SRC_DIR ./src)

# Include the directories and now your cpp files will recognize your headers
include_directories(${INCL_DIR})

# Set source files
set(
        SOURCE_FILES
        ${INCL_DIR}/CLI.h
        ${SRC_DIR}/main.cpp
        ${SRC_DIR}/CLI.cpp)
#[[
set(
        SOURCE_FILES
        main.cpp
        ${APP_DOMAIN}/person.h
        ${APP_DOMAIN}/person.cpp
        ${APP_CALC}/add.h
        ${APP_CALC}/add.cpp
)
]]
#add_executable(contextsv ContextSV.cpp)

add_executable(contextsv ${SOURCE_FILES})
