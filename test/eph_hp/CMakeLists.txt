# test/eph_hp/CMakeLists.txt
# ===========================
# 高精度天体位置计算的"备件"测试目录.
#
# 入库原则: 只入库脚本 / C++ / CMake / Python / README.
#   - data/*                       通过 fetch_data.sh 按需下载, 不入库
#   - vsop87d_earth_hp.h           由 tools/vsop87d_to_header.py 生成, 不入库
#   - iau2000a_nut_lon_hp.h        由 tools/iau2000a_to_header.py 生成, 不入库
#   - vsop2013_emb_reference.csv   由 tools/vsop2013_eval.py 生成, 不入库
#
# 构建流程:
#   1. cd test/eph_hp && ./fetch_data.sh    # 一次性下数据 + 生成 header
#   2. cmake --build build --target eph_hp_test
#
# 若数据/header 缺失, CMake 会在 configure 时跳过对应 target 并打印提示,
# 主库 libsxwnl 和其他测试完全不受影响.

cmake_minimum_required(VERSION 3.25)

# ---------------------------------------------------------------------------
# 用到的路径
# ---------------------------------------------------------------------------
set(HP_DIR         ${CMAKE_CURRENT_SOURCE_DIR})
set(HP_DATA_DIR    ${HP_DIR}/data)
set(HP_TOOLS_DIR   ${HP_DIR}/tools)
set(VSOP87D_EAR    ${HP_DATA_DIR}/VSOP87D.ear)
set(IAU2000A_TXT   ${HP_DATA_DIR}/IAU2000A_tab5.3a.txt)
set(VSOP87D_HDR    ${HP_DIR}/vsop87d_earth_hp.h)
set(IAU2000A_HDR   ${HP_DIR}/iau2000a_nut_lon_hp.h)

# ---------------------------------------------------------------------------
# 检查数据/生成产物就绪情况
# ---------------------------------------------------------------------------
set(HP_READY TRUE)

# 优先检查是否已有生成的 header. 若缺失但原始数据在, 自动补生成一次.
if(NOT EXISTS ${VSOP87D_HDR} OR NOT EXISTS ${IAU2000A_HDR})
    if(EXISTS ${VSOP87D_EAR} AND EXISTS ${IAU2000A_TXT})
        find_program(PYTHON3_EXECUTABLE NAMES python3 python)
        if(PYTHON3_EXECUTABLE)
            message(STATUS "HP: 原始数据已下载但 header 未生成, 现在自动生成 ...")
            if(NOT EXISTS ${VSOP87D_HDR})
                execute_process(
                    COMMAND ${PYTHON3_EXECUTABLE}
                            ${HP_TOOLS_DIR}/vsop87d_to_header.py
                            ${VSOP87D_EAR} ${VSOP87D_HDR}
                    RESULT_VARIABLE _rv)
                if(NOT _rv EQUAL 0)
                    message(WARNING "HP: 生成 vsop87d_earth_hp.h 失败")
                endif()
            endif()
            if(NOT EXISTS ${IAU2000A_HDR})
                execute_process(
                    COMMAND ${PYTHON3_EXECUTABLE}
                            ${HP_TOOLS_DIR}/iau2000a_to_header.py
                            ${IAU2000A_TXT} ${IAU2000A_HDR}
                    RESULT_VARIABLE _rv)
                if(NOT _rv EQUAL 0)
                    message(WARNING "HP: 生成 iau2000a_nut_lon_hp.h 失败")
                endif()
            endif()
        else()
            message(STATUS "HP: 未找到 python3, 无法自动生成 header")
            set(HP_READY FALSE)
        endif()
    else()
        set(HP_READY FALSE)
    endif()
endif()

# 二次检查: 生成后 header 是否真的到位了
if(NOT EXISTS ${VSOP87D_HDR} OR NOT EXISTS ${IAU2000A_HDR})
    set(HP_READY FALSE)
endif()

# ---------------------------------------------------------------------------
# 主 HP 测试 (Test A/B/C/D). 数据/header 就绪才定义 target.
# ---------------------------------------------------------------------------
if(HP_READY)
    add_executable(eph_hp_test
        eph_hp_test.cpp
        eph_hp.cpp
    )
    target_link_libraries(eph_hp_test PRIVATE sxwnl)
    target_include_directories(eph_hp_test PRIVATE
        ${HP_DIR}
        ${PROJECT_SOURCE_DIR}/src
    )
    set_target_properties(eph_hp_test PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
    message(STATUS "HP: eph_hp_test 已启用 (VSOP87D + IAU 2000A 全项对比)")
else()
    message(STATUS
        "HP: eph_hp_test 已跳过 (缺少数据/header)\n"
        "     运行以下命令下载并生成:\n"
        "         cd test/eph_hp && ./fetch_data.sh")
endif()

# ---------------------------------------------------------------------------
# 可选: DE441/DE440 二进制星历测试, 依赖 IMCCE calceph 库
# ---------------------------------------------------------------------------
find_package(calceph CONFIG QUIET)
if(calceph_FOUND AND HP_READY)
    message(STATUS "HP: calceph 已找到, 启用 de_hp_test 目标")
    add_executable(de_hp_test
        de_hp_test.cpp
        de_hp.cpp
        eph_hp.cpp    # 用 E_Lon_HP 做 DE 结果的量级对照
    )
    target_link_libraries(de_hp_test PRIVATE sxwnl calceph)
    target_include_directories(de_hp_test PRIVATE
        ${HP_DIR}
        ${PROJECT_SOURCE_DIR}/src
    )
    set_target_properties(de_hp_test PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    )
    # 把星历文件路径通过宏传给 C++, 免去写死或环境变量
    target_compile_definitions(de_hp_test PRIVATE
        DE_EPHEMERIS_PATH="${HP_DATA_DIR}/de440s.bsp"
    )
elseif(calceph_FOUND)
    message(STATUS "HP: calceph 已找到但 HP 数据缺失, de_hp_test 跳过")
else()
    message(STATUS "HP: calceph 未安装, de_hp_test 跳过 "
                   "(brew install calceph / apt install libcalceph-dev)")
endif()
