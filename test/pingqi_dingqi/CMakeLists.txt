# 项目根目录下的CMakeLists.txt
cmake_minimum_required(VERSION 3.25)
# 平气/定气 八字模式在节气交接时刻的自洽性测试
add_executable(pingqi_dingqi_test pingqi_dingqi_test.cpp)
target_link_libraries(pingqi_dingqi_test PRIVATE sxwnl)
target_include_directories(pingqi_dingqi_test PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_target_properties(pingqi_dingqi_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
