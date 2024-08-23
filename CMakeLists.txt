# 项目根目录下的CMakeLists.txt
cmake_minimum_required(VERSION 3.25)
project(sxwnl)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# 设置库和可执行文件的输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# 添加src目录到构建
add_subdirectory(src)

# 添加test目录到构建
add_subdirectory(test)
