# 阴历月名(古历"十三月/后九月"与普通闰月)验证测试
# 通过导出的 C API (sxwnl_capi) 端到端校验真实显示逻辑, 故需一并编译 capi 源。
cmake_minimum_required(VERSION 3.25)

add_executable(lunar_month_name_test
    lunar_month_name_test.cpp
    ${PROJECT_SOURCE_DIR}/capi/sxwnl_capi.cpp
    ${PROJECT_SOURCE_DIR}/capi/sxwnl_eclipse_map.cpp
)

target_link_libraries(lunar_month_name_test PRIVATE sxwnl)

target_include_directories(lunar_month_name_test PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/capi
)

set_target_properties(lunar_month_name_test PROPERTIES
    CXX_STANDARD 17
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
