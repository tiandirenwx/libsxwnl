# test/CMakeLists.txt

# 收集test目录下的所有.cpp文件作为测试程序的源文件
file(GLOB TEST_SRC_FILES "*.cpp")

# 添加测试程序的可执行文件
add_executable(test_sxwnl ${TEST_SRC_FILES})

# 链接sxwnl库到测试程序
target_link_libraries(test_sxwnl PRIVATE sxwnl)

# 添加库的头文件路径到测试程序的头文件搜索路径
# 注意这里的${PROJECT_SOURCE_DIR}/src是假设您的头文件在项目的src目录下
target_include_directories(test_sxwnl PRIVATE ${PROJECT_SOURCE_DIR}/src)

# 确保使用相对路径引用库
set_target_properties(test_sxwnl PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)