# 八字排盘 —— 规格符合性 / 黑盒不变量测试
cmake_minimum_required(VERSION 3.25)
add_executable(bazi_invariants_test bazi_invariants_test.cpp)
target_link_libraries(bazi_invariants_test PRIVATE sxwnl)
target_include_directories(bazi_invariants_test PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_target_properties(bazi_invariants_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
