# 起运/交运 算法验证用例
cmake_minimum_required(VERSION 3.25)
add_executable(qiyun_jiaoyun_test qiyun_jiaoyun_test.cpp)
target_link_libraries(qiyun_jiaoyun_test PRIVATE sxwnl)
target_include_directories(qiyun_jiaoyun_test PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_target_properties(qiyun_jiaoyun_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
