cmake_minimum_required(VERSION 3.25)

set(SXWNL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)

file(GLOB SXWNL_CORE_SRCS "${SXWNL_SRC_DIR}/*.cpp")

add_library(sxwnl_capi STATIC
    sxwnl_capi.cpp
    sxwnl_eclipse_map.cpp
    ${SXWNL_CORE_SRCS}
)

target_include_directories(sxwnl_capi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SXWNL_SRC_DIR}
)

set_target_properties(sxwnl_capi PROPERTIES
    CXX_STANDARD 17
    C_STANDARD 11
    POSITION_INDEPENDENT_CODE ON
)
