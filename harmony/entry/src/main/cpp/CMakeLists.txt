cmake_minimum_required(VERSION 3.16)
project(sxwnl_napi)

set(CMAKE_CXX_STANDARD 17)

set(CAPI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../capi)
set(SXWNL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../../../src)

file(GLOB SXWNL_CORE_SRCS CONFIGURE_DEPENDS "${SXWNL_SRC_DIR}/*.cpp")

add_library(sxwnl_napi SHARED
    napi_init.cpp
    napi_calendar.cpp
    napi_bazi.cpp
    napi_eclipse_map.cpp
    ${CAPI_DIR}/sxwnl_capi.cpp
    ${CAPI_DIR}/sxwnl_eclipse_map.cpp
    ${SXWNL_CORE_SRCS}
)

target_include_directories(sxwnl_napi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CAPI_DIR}
    ${SXWNL_SRC_DIR}
)

# 全局强制包含 long double 兼容垫片: 绕开 arm64 (128 位 binary128) libm 缺陷
# (如 sinl 对某些参数返回 0, 导致真太阳时/天文计算错误)。详见 ld_compat.h。
target_compile_options(sxwnl_napi PRIVATE
    -include ${CMAKE_CURRENT_SOURCE_DIR}/ld_compat.h
)

target_link_libraries(sxwnl_napi PUBLIC libace_napi.z.so)
