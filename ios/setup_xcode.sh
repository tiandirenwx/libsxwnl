#!/bin/bash
# ═══════════════════════════════════════════════════════════
# 寿星万年历 - macOS/iOS Xcode 项目快速设置脚本
# 用法: cd libsxwnl/ios && bash setup_xcode.sh
# ═══════════════════════════════════════════════════════════

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_DIR="$SCRIPT_DIR/SxwnlCalendar"

echo "╔══════════════════════════════════════════╗"
echo "║  寿星万年历 - macOS 项目设置             ║"
echo "╚══════════════════════════════════════════╝"
echo ""

# Step 1: 验证 C++ 核心在 macOS 上可编译
echo "➤ Step 1: 编译验证 C API..."
cd "$ROOT_DIR/capi"

if command -v g++ &>/dev/null; then
    CXX=g++
elif command -v clang++ &>/dev/null; then
    CXX=clang++
else
    echo "  ✗ 未找到 C++ 编译器，请安装 Xcode Command Line Tools:"
    echo "    xcode-select --install"
    exit 1
fi

$CXX -std=c++17 -I../src -I. test_capi.cpp sxwnl_capi.cpp ../src/*.cpp -o test_capi_macos 2>&1
if [ $? -eq 0 ]; then
    echo "  ✓ C++ 编译成功"
    echo ""
    echo "  运行功能测试..."
    ./test_capi_macos
    echo ""
    rm -f test_capi_macos
else
    echo "  ✗ 编译失败，请检查错误信息"
    exit 1
fi

# Step 2: 创建 Xcode 项目目录结构
echo "➤ Step 2: 准备 Xcode 项目文件..."

# 创建源文件符号链接（避免复制）
XCODE_SRC="$PROJECT_DIR/CppSources"
mkdir -p "$XCODE_SRC"

# 链接 src/ 下所有 .cpp / .h (含 world_map* 等)
for f in "$ROOT_DIR/src"/*.cpp "$ROOT_DIR/src"/*.h; do
    ln -sf "$f" "$XCODE_SRC/" 2>/dev/null || true
done
# 链接 capi/ 下所有 .cpp / .h (sxwnl_capi.*, sxwnl_eclipse_map.* 等)
for f in "$ROOT_DIR/capi"/*.cpp "$ROOT_DIR/capi"/*.h; do
    ln -sf "$f" "$XCODE_SRC/" 2>/dev/null || true
done

echo "  ✓ C++ 源码已链接到 $XCODE_SRC"

# Step 3: 输出 Xcode 手动配置指引
echo ""
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "  ✓ C++ 核心编译验证通过！"
echo ""
echo "  接下来请在 Xcode 中完成以下步骤："
echo ""
echo "  1. 打开 Xcode → File → New → Project"
echo "     → Multiplatform → App → Product Name: SxwnlCalendar"
echo "     → 保存到: $SCRIPT_DIR/"
echo ""
echo "  2. 删除 Xcode 自动生成的 ContentView.swift"
echo ""
echo "  3. 将以下文件拖入项目 (Add to target ✓):"
echo "     ├── SxwnlCalendar/SxwnlCalendarApp.swift"
echo "     ├── SxwnlCalendar/Bridge/SxwnlBridge.swift"
echo "     ├── SxwnlCalendar/Views/ (所有 .swift)"
echo "     ├── SxwnlCalendar/Models/ (所有 .swift)"
echo "     └── SxwnlCalendar/CppSources/ (所有 .cpp 和 .h)"
echo ""
echo "  4. 第一次添加 .cpp 文件时，Xcode 会提示创建 Bridging Header"
echo "     → 选择 Create，然后编辑内容为:"
echo '     #import "sxwnl_capi.h"'
echo ""
echo "  5. Build Settings 搜索并设置:"
echo "     → C++ Language Dialect: C++17"
echo "     → Header Search Paths 添加:"
echo "       $XCODE_SRC"
echo ""
echo "  6. 选择 My Mac 作为目标 → ⌘B 编译 → ⌘R 运行"
echo ""
echo "═══════════════════════════════════════════════════════════"
