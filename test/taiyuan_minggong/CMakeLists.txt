# 胎元 / 命宫 —— 规格符合性测试
cmake_minimum_required(VERSION 3.25)
add_executable(taiyuan_minggong_test taiyuan_minggong_test.cpp)
target_link_libraries(taiyuan_minggong_test PRIVATE sxwnl)
target_include_directories(taiyuan_minggong_test PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_target_properties(taiyuan_minggong_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
