cmake_minimum_required(VERSION 3.22)
project(sxwnl_jni)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 项目根目录：android/app/src/main/cpp -> libsxwnl/  (cpp -> main -> src -> app -> android -> libsxwnl)
set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../../../..)

# 源码目录（已确认存在）
set(SRC_DIR ${PROJECT_ROOT}/src)
set(CAPI_DIR ${PROJECT_ROOT}/capi)

# 自动收集核心 C++ 源文件
file(GLOB SXWNL_CORE_SRCS "${SRC_DIR}/*.cpp")

add_library(sxwnl_jni SHARED
        jni_bridge.cpp
        jni_eclipse_map.cpp
        ${CAPI_DIR}/sxwnl_capi.cpp
        ${CAPI_DIR}/sxwnl_eclipse_map.cpp
        ${SXWNL_CORE_SRCS}
)

target_include_directories(sxwnl_jni PRIVATE
        ${CAPI_DIR}
        ${SRC_DIR}
)

# 全局强制包含 long double 兼容垫片: 绕开安卓 arm64 bionic 的 binary128 libm 缺陷
# (如 sinl 对某些参数返回 0, 导致真太阳时/天文计算错误)。详见 ld_compat.h。
target_compile_options(sxwnl_jni PRIVATE
        -include ${CMAKE_CURRENT_SOURCE_DIR}/ld_compat.h
)

target_link_libraries(sxwnl_jni log)