#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# 寿星万年历 — iOS / macOS 一键打包脚本
#
# 用法:
#   scripts/build_ios.sh                  # 默认: 模拟器 .app
#   scripts/build_ios.sh sim              # iOS 模拟器 .app
#   scripts/build_ios.sh device           # iOS 真机 .ipa (archive + export)
#   scripts/build_ios.sh mac              # macOS .app (Mac Catalyst/原生)
#   scripts/build_ios.sh all              # sim + device + mac
#
# 环境要求:
#   - 已安装 Xcode 和 Command Line Tools (xcode-select --install)
#   - device 打包: Xcode 已登录开发者账号, Apple Developer 团队 ID=9TKPYKJA67
#   - 真机包导出配置: ios/ExportOptions.plist (可按需修改 method/teamID)
#
# 产物输出:
#   <repo>/dist/sxwnl-ios-sim.app/        模拟器 (拖入模拟器即可)
#   <repo>/dist/sxwnl-ios-device.ipa      真机包 (用 Xcode/Apple Configurator 安装)
#   <repo>/dist/sxwnl-macos.app/          macOS 桌面版
#
# 清理策略:
#   - 默认每次打包前清理: dist/ 旧产物 + ios/build/ci 中本次目标的 DerivedData/archive/export
#   - xcodebuild 命令显式使用 'clean build' / 'clean archive', 确保中间产物全新
#   - 设置环境变量 SXWNL_NO_CLEAN=1 可跳过, 适合本地增量调试
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

# ── 颜色 ──
RED='\033[1;31m'; GREEN='\033[1;32m'; YELLOW='\033[1;33m'; BLUE='\033[1;34m'; NC='\033[0m'
info()  { printf "${BLUE}[INFO]${NC} %s\n" "$*"; }
ok()    { printf "${GREEN}[ OK ]${NC} %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
err()   { printf "${RED}[ERR ]${NC} %s\n" "$*" >&2; }
die()   { err "$*"; exit 1; }

# ── 路径 ──
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
IOS_DIR="$REPO_ROOT/ios"
PROJECT="$IOS_DIR/SxwnlCalendar/SxwnlCalendar.xcodeproj"
SCHEME="SxwnlCalendar"
EXPORT_PLIST="$IOS_DIR/ExportOptions.plist"
DIST_DIR="$REPO_ROOT/dist"
BUILD_DIR="$IOS_DIR/build/ci"

[[ -d "$PROJECT" ]] || die "找不到 Xcode 项目: $PROJECT"
command -v xcodebuild >/dev/null 2>&1 || die "未找到 xcodebuild, 请安装 Xcode 与 Command Line Tools"
mkdir -p "$DIST_DIR" "$BUILD_DIR"

# ── 参数 ──
TARGET="${1:-sim}"
case "$TARGET" in
    sim|device|mac|all) ;;
    -h|--help)
        awk '/^set -euo/{exit} NR>1' "$0"; exit 0 ;;
    *) die "未知参数: $TARGET (支持: sim | device | mac | all)" ;;
esac

# ── SDK / Runtime 可用性检测 ──
have_ios_sim_runtime() {
    xcrun simctl list runtimes 2>/dev/null | grep -qE "^iOS[[:space:]]+[0-9]"
}

have_ios_device_dest() {
    xcodebuild -project "$PROJECT" -scheme "$SCHEME" -showdestinations 2>&1 \
        | awk '/Available destinations/{flag=1;next} /Ineligible destinations/{flag=0} flag' \
        | grep -q "platform:iOS,"
}

have_macos_dest() {
    xcodebuild -project "$PROJECT" -scheme "$SCHEME" -showdestinations 2>&1 \
        | awk '/Available destinations/{flag=1;next} /Ineligible destinations/{flag=0} flag' \
        | grep -q "platform:macOS"
}

hint_install_platform() {
    err "  -> 打开 Xcode > Settings > Components, 下载 $1 后重试"
    err "  -> 或命令行: xcodebuild -downloadPlatform $2"
}

# ── 清理 ──
# 在每个 build_xxx 入口调用, 删除该目标的 DerivedData / archive / dist 旧产物
# SXWNL_NO_CLEAN=1 可跳过 (本地反复调试时加快速度)
should_clean() { [[ "${SXWNL_NO_CLEAN:-0}" != "1" ]]; }

clean_paths() {
    # 用法: clean_paths "目标名" path1 path2 ...
    local label="$1"; shift
    if ! should_clean; then
        warn "已设 SXWNL_NO_CLEAN=1, 跳过 $label 清理 (增量构建)"
        return 0
    fi
    info "清理 $label:"
    local p
    for p in "$@"; do
        if [[ -e "$p" ]]; then
            printf "    rm -rf %s\n" "$p"
            rm -rf "$p"
        fi
    done
}

# ── 工具函数: 包一层 xcodebuild, 失败时立即 die ──
#   注: 必须显式调 die, 因为 bash 中函数 return 非零在管道中
#       不会触发外层 set -e (尤其与 || / && 组合时)
xc() {
    local log="$BUILD_DIR/xcodebuild.log"
    info "xcodebuild $*"
    set +e
    xcodebuild "$@" 2>&1 | tee "$log" >/dev/null
    local status=${PIPESTATUS[0]}
    set -e
    if [[ $status -ne 0 ]]; then
        err "xcodebuild 失败 (exit=$status); 最近 40 行日志:"
        tail -40 "$log" | sed 's/^/    /' >&2
        die "请修正上述错误后重试"
    fi
}

# ── 1) 模拟器 .app ──
build_sim() {
    info "═══ 构建 iOS 模拟器版本 ═══"
    if ! have_ios_sim_runtime; then
        err "未检测到任何 iOS Simulator 运行时 (xcrun simctl list runtimes 为空)"
        hint_install_platform "iOS Simulator Platform" "iOS"
        return 1
    fi

    local derived="$BUILD_DIR/sim-DerivedData"
    clean_paths "iOS Simulator" "$derived" "$DIST_DIR/sxwnl-ios-sim.app"

    # 即使保留 derived (NO_CLEAN), 也强制 xcodebuild 走 'clean build' 流程,
    # 确保中间链接/资源/编译信息每次都重做.
    xc -project "$PROJECT" \
       -scheme "$SCHEME" \
       -configuration Debug \
       -sdk iphonesimulator \
       -destination 'generic/platform=iOS Simulator' \
       -derivedDataPath "$derived" \
       CODE_SIGNING_ALLOWED=NO CODE_SIGNING_REQUIRED=NO \
       clean build

    local app_src
    app_src="$(find "$derived/Build/Products" -maxdepth 3 -name "${SCHEME}.app" -type d | head -1)"
    [[ -n "$app_src" ]] || die "找不到模拟器 .app 产物"
    local dst="$DIST_DIR/sxwnl-ios-sim.app"
    rm -rf "$dst"; cp -R "$app_src" "$dst"
    ok ".app -> $dst"
    info "安装到模拟器: xcrun simctl install booted \"$dst\""
}

# ── 2) 真机 .ipa (archive -> export) ──
build_device() {
    info "═══ 构建 iOS 真机版本 (archive + export) ═══"
    if ! have_ios_device_dest; then
        err "未检测到 iOS Device 可用 destination (平台未完整安装)"
        hint_install_platform "iOS Platform" "iOS"
        return 1
    fi
    [[ -f "$EXPORT_PLIST" ]] || die "找不到导出配置: $EXPORT_PLIST"

    local derived="$BUILD_DIR/device-DerivedData"
    local archive="$BUILD_DIR/SxwnlCalendar.xcarchive"
    local export_dir="$BUILD_DIR/export"

    # 这三个都必须清理: 旧 archive 会被 exportArchive 复用; 旧 derived 会被增量编译复用;
    # 旧 export 会让 cp 拿到上次的 .ipa
    clean_paths "iOS Device" "$derived" "$archive" "$export_dir" "$DIST_DIR/sxwnl-ios-device.ipa"

    # archive 命令显式 'clean archive', 并把 DerivedData 锁定到 BUILD_DIR (默认会写 ~/Library/.../DerivedData)
    xc -project "$PROJECT" \
       -scheme "$SCHEME" \
       -configuration Release \
       -sdk iphoneos \
       -destination 'generic/platform=iOS' \
       -archivePath "$archive" \
       -derivedDataPath "$derived" \
       -allowProvisioningUpdates \
       clean archive

    xc -exportArchive \
       -archivePath "$archive" \
       -exportPath "$export_dir" \
       -exportOptionsPlist "$EXPORT_PLIST" \
       -allowProvisioningUpdates

    local ipa
    ipa="$(find "$export_dir" -maxdepth 2 -name "*.ipa" | head -1)"
    [[ -n "$ipa" ]] || die "找不到 .ipa 产物 (检查 ExportOptions.plist 配置)"
    local dst="$DIST_DIR/sxwnl-ios-device.ipa"
    cp -f "$ipa" "$dst"
    ok ".ipa -> $dst ($(du -h "$dst" | cut -f1))"
    info "安装方式: Xcode > Window > Devices, 或 Apple Configurator 2 拖入"
}

# ── 3) macOS .app ──
# 签名策略:
#   - 默认: 跳过签名 (CODE_SIGNING_ALLOWED=NO), 适合本机运行 / 内部分发
#   - 设置 SXWNL_MAC_SIGN=1: 用自动签名 (需要 Xcode 已登录开发者账号)
build_mac() {
    info "═══ 构建 macOS 桌面版 ═══"
    if ! have_macos_dest; then
        err "未检测到 macOS 可用 destination"
        hint_install_platform "macOS Platform" "macOS"
        return 1
    fi

    local derived="$BUILD_DIR/mac-DerivedData"
    clean_paths "macOS" "$derived" "$DIST_DIR/sxwnl-macos.app"

    local sign_args=()
    if [[ "${SXWNL_MAC_SIGN:-0}" == "1" ]]; then
        info "使用 Xcode 自动签名 (SXWNL_MAC_SIGN=1)"
        sign_args=(-allowProvisioningUpdates CODE_SIGN_STYLE=Automatic)
    else
        info "跳过签名 (仅本机运行用; 如需分发设置 SXWNL_MAC_SIGN=1)"
        sign_args=(CODE_SIGNING_ALLOWED=NO CODE_SIGNING_REQUIRED=NO CODE_SIGN_IDENTITY='')
    fi

    xc -project "$PROJECT" \
       -scheme "$SCHEME" \
       -configuration Release \
       -destination 'platform=macOS' \
       -derivedDataPath "$derived" \
       "${sign_args[@]}" \
       clean build

    local app_src
    app_src="$(find "$derived/Build/Products" -maxdepth 3 -name "${SCHEME}.app" -type d | head -1)"
    [[ -n "$app_src" ]] || die "找不到 macOS .app 产物"
    local dst="$DIST_DIR/sxwnl-macos.app"
    rm -rf "$dst"; cp -R "$app_src" "$dst"
    ok ".app -> $dst"
    info "直接运行: open \"$dst\""
    if [[ "${SXWNL_MAC_SIGN:-0}" != "1" ]]; then
        warn "未签名包首次运行可能被 Gatekeeper 拦截:"
        warn "  -> 右键 .app > 打开, 或: xattr -dr com.apple.quarantine \"$dst\""
    fi
}

# ── 执行 ──
echo "╔═══════════════════════════════════════╗"
echo "║  寿星万年历 - iOS/macOS 打包 ($TARGET)"
echo "╚═══════════════════════════════════════╝"

FAIL=0
case "$TARGET" in
    sim)    build_sim    || FAIL=1 ;;
    device) build_device || FAIL=1 ;;
    mac)    build_mac    || FAIL=1 ;;
    all)
        build_sim    || FAIL=1
        build_device || FAIL=1
        build_mac    || FAIL=1
        ;;
esac

if [[ $FAIL -ne 0 ]]; then
    err "部分目标打包失败 (见上方日志)"
    exit 1
fi
ok "iOS/macOS 打包完成 (dist/ 目录)"
