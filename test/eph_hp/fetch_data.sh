#!/usr/bin/env bash
# fetch_data.sh
# ==============
# 一键把 HP 备件的**全部**数据从 IMCCE / IERS / NASA 官方源下载, 并把
# C++ 需要的 header (vsop87d_earth_hp.h / iau2000a_nut_lon_hp.h) 生成好.
# 也可以顺带生成 vsop2013_emb_reference.csv (VSOP2013 参考评估输出).
#
# 仓库入库原则:
#   * 只入库脚本 / C++/CMake/Python 源码 / README / .cursor 规则
#   * 全部 data/* 和生成的 .h / .csv 都通过本脚本按需产出
#
# 用法:
#   ./fetch_data.sh                    # 默认: 下小文件 + vsop2013 + de440, 生成 header
#   ./fetch_data.sh small              # 只下 4 个小参考文件, 生成 header (不下大文件)
#   ./fetch_data.sh vsop2013           # 只下 VSOP2013 大文件
#   ./fetch_data.sh de440              # 只下 DE440 星历
#   ./fetch_data.sh de441              # 下完整 3.2 GB 的 DE441
#   ./fetch_data.sh gen                # 只跑 Python 生成 header/CSV (数据要已存在)
#
# 环境:
#   - 需要 curl (macOS/Linux 自带)
#   - 生成 header 需要 python3

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
TOOLS_DIR="$SCRIPT_DIR/tools"
mkdir -p "$DATA_DIR"

CURL_OPTS="--fail --show-error --location -C - --retry 3"

# --- 通用下载函数: 幂等, 已存在且 >= min_bytes 就跳过 ---
fetch() {
    local url="$1"
    local out="$2"
    local min_bytes="$3"
    if [ -f "$out" ]; then
        local sz
        sz=$(wc -c <"$out" | tr -d ' ')
        if [ "$sz" -ge "$min_bytes" ]; then
            echo "  [SKIP] $(basename "$out") 已存在 ($sz bytes)"
            return 0
        fi
        echo "  [重下] $(basename "$out") 存在但太小 ($sz < $min_bytes)"
    fi
    echo "  [下载] $url"
    curl $CURL_OPTS -o "$out.tmp" "$url" && mv "$out.tmp" "$out"
    echo "  [完成] $(basename "$out") -> $out"
}

# --- 4 个小参考文件 ---
fetch_small() {
    echo "== 小参考文件 (400 KB 级, IMCCE / IERS) =="
    fetch "https://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87D.ear" \
          "$DATA_DIR/VSOP87D.ear" \
          300000
    fetch "https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3a.txt" \
          "$DATA_DIR/IAU2000A_tab5.3a.txt" \
          60000
    fetch "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013.f" \
          "$DATA_DIR/VSOP2013.f" \
          10000
    fetch "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013.ctl" \
          "$DATA_DIR/VSOP2013.ctl" \
          30000
}

fetch_vsop2013() {
    echo "== VSOP2013 (地-月轨道要素, 33 MB) =="
    fetch "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013p3.dat" \
          "$DATA_DIR/VSOP2013p3.dat" \
          30000000
}

fetch_de440() {
    echo "== JPL DE440s (二进制星历, 1849~2150 CE, 33 MB) =="
    fetch "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp" \
          "$DATA_DIR/de440s.bsp" \
          30000000
}

fetch_de441() {
    echo "== JPL DE441 (二进制星历, -13000~17191 CE, ~3.2 GB) =="
    echo "  注意: 文件极大, 会占用较长下载时间和磁盘空间"
    fetch "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp" \
          "$DATA_DIR/de441.bsp" \
          3000000000
}

# --- Python 生成 C++ header 和 CSV 参考 ---
run_generators() {
    echo "== 生成 C++ header 与参考 CSV =="
    if ! command -v python3 >/dev/null 2>&1; then
        echo "  [错误] 未找到 python3, 无法生成 header. 跳过." >&2
        return 1
    fi

    if [ -f "$DATA_DIR/VSOP87D.ear" ]; then
        echo "  [生成] vsop87d_earth_hp.h"
        python3 "$TOOLS_DIR/vsop87d_to_header.py" \
                "$DATA_DIR/VSOP87D.ear" \
                "$SCRIPT_DIR/vsop87d_earth_hp.h"
    else
        echo "  [跳过] VSOP87D.ear 未下载, 无法生成 vsop87d_earth_hp.h" >&2
    fi

    if [ -f "$DATA_DIR/IAU2000A_tab5.3a.txt" ]; then
        echo "  [生成] iau2000a_nut_lon_hp.h"
        python3 "$TOOLS_DIR/iau2000a_to_header.py" \
                "$DATA_DIR/IAU2000A_tab5.3a.txt" \
                "$SCRIPT_DIR/iau2000a_nut_lon_hp.h"
    else
        echo "  [跳过] IAU2000A_tab5.3a.txt 未下载, 无法生成 iau2000a_nut_lon_hp.h" >&2
    fi

    if [ -f "$DATA_DIR/VSOP2013p3.dat" ]; then
        echo "  [生成] vsop2013_emb_reference.csv (-2000 ~ +8000, 500 年步长)"
        python3 "$TOOLS_DIR/vsop2013_eval.py" \
                "$DATA_DIR/VSOP2013p3.dat" \
                --csv "$SCRIPT_DIR/vsop2013_emb_reference.csv" \
                --range -2000 8000 500
    else
        echo "  [跳过] VSOP2013p3.dat 未下载, 无法生成参考 CSV"
    fi
}

# --- 主流程 ---
if [ $# -eq 0 ]; then
    fetch_small
    fetch_vsop2013
    fetch_de440
    run_generators
    echo ""
    echo "默认下载 + 生成完成. 3.2 GB 的 DE441 完整版请单独运行:"
    echo "  $0 de441"
    exit 0
fi

for arg in "$@"; do
    case "$arg" in
        small)    fetch_small ;;
        vsop2013) fetch_vsop2013 ;;
        de440)    fetch_de440 ;;
        de441)    fetch_de441 ;;
        gen)      run_generators ;;
        all)
            fetch_small
            fetch_vsop2013
            fetch_de440
            run_generators
            ;;
        *)
            echo "未知选项: $arg" >&2
            echo "可用: small | vsop2013 | de440 | de441 | gen | all | (无参数=默认)" >&2
            exit 2
            ;;
    esac
done

# 如果下载了原始数据但没显式跑生成, 自动补一遍
case " $* " in
    *" gen "*|*" all "*) ;;
    *)
        if [ -f "$DATA_DIR/VSOP87D.ear" ] || [ -f "$DATA_DIR/IAU2000A_tab5.3a.txt" ] \
           || [ -f "$DATA_DIR/VSOP2013p3.dat" ]; then
            echo ""
            run_generators
        fi
        ;;
esac
