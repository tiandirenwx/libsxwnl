# 历谱口径节气 (Day::getLiPuJieQiName, 对应权威 sxwnl ob.Ljq) 回归用例
cmake_minimum_required(VERSION 3.25)
add_executable(lipu_jieqi_test lipu_jieqi_test.cpp)
target_link_libraries(lipu_jieqi_test PRIVATE sxwnl)
target_include_directories(lipu_jieqi_test PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_target_properties(lipu_jieqi_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
