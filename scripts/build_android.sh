#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# 寿星万年历 — Android 一键打包脚本
#
# 用法:
#   scripts/build_android.sh                # 默认: debug + release APK (本地分发)
#   scripts/build_android.sh both           # debug + release APK (本地分发, 同上)
#   scripts/build_android.sh debug          # 仅 Debug APK
#   scripts/build_android.sh release        # 仅 Release APK (自动签名)
#   scripts/build_android.sh aab            # Release AAB (Google Play)
#   scripts/build_android.sh all            # debug + release APK + AAB
#
# 签名优先级 (release/aab):
#   1) <repo>/android/keystore.properties (storeFile/storePassword/keyAlias/keyPassword)
#   2) 环境变量 SXWNL_KEYSTORE_PATH / SXWNL_KEYSTORE_PASSWORD /
#                SXWNL_KEY_ALIAS / SXWNL_KEY_PASSWORD
#   3) 回退到 ~/.android/debug.keystore (适合开发/内测; 不可上架商店)
#
# 产物输出:
#   <repo>/dist/sxwnl-android-<variant>.{apk,aab}
#
# 清理策略:
#   - 默认每次打包前都会执行: ./gradlew clean + 删除 dist/sxwnl-android-*.{apk,aab}
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
ANDROID_DIR="$REPO_ROOT/android"
DIST_DIR="$REPO_ROOT/dist"

[[ -f "$ANDROID_DIR/gradlew" ]] || die "找不到 gradlew: $ANDROID_DIR/gradlew"
mkdir -p "$DIST_DIR"

# ── 参数 ──
TARGET="${1:-both}"
case "$TARGET" in
    debug|release|both|aab|all) ;;
    -h|--help)
        awk '/^set -euo/{exit} NR>1' "$0"; exit 0 ;;
    *) die "未知参数: $TARGET (支持: debug | release | both | aab | all)" ;;
esac

# ── 校验签名 (release/aab) ──
check_signing() {
    local ksprops="$ANDROID_DIR/keystore.properties"
    if [[ -f "$ksprops" ]]; then
        ok "使用项目级签名: $ksprops"
    elif [[ -n "${SXWNL_KEYSTORE_PATH:-}" ]]; then
        ok "使用环境变量签名: $SXWNL_KEYSTORE_PATH"
    else
        local debug_ks="$HOME/.android/debug.keystore"
        if [[ -f "$debug_ks" ]]; then
            warn "未配置正式签名, 回退到 debug.keystore (开发/内测可用, 不可上架商店)"
            warn "  -> 自定义签名: 创建 $ksprops 或导出 SXWNL_KEYSTORE_PATH 等环境变量"
        else
            warn "未找到 debug.keystore, gradle 将在首次构建时自动创建"
        fi
    fi
}

# ── 清理 ──
# 默认每次打包前清理: gradle clean + 删除 dist 中本平台旧产物
# 设置 SXWNL_NO_CLEAN=1 可跳过 (本地反复调试时加快速度)
clean_all() {
    if [[ "${SXWNL_NO_CLEAN:-0}" == "1" ]]; then
        warn "已设 SXWNL_NO_CLEAN=1, 跳过 clean (增量构建)"
        return 0
    fi
    info "清理: 删除 dist/ 中 Android 旧产物"
    rm -f "$DIST_DIR"/sxwnl-android-*.apk "$DIST_DIR"/sxwnl-android-*.aab 2>/dev/null || true

    info "清理: ./gradlew clean"
    (cd "$ANDROID_DIR" && ./gradlew --console=plain clean) \
        || die "gradle clean 失败"

    # 额外清理 jni 中间产物 (有时 CMake 缓存会卡住 NDK 升级或脚本变更)
    if [[ -d "$ANDROID_DIR/app/.cxx" ]]; then
        info "清理: app/.cxx (NDK CMake 缓存)"
        rm -rf "$ANDROID_DIR/app/.cxx"
    fi
    ok "清理完成"
}

# ── 构建 ──
# 支持单 task 或多 task: run_gradle assembleDebug
#                       run_gradle assembleDebug assembleRelease
run_gradle() {
    info "执行: ./gradlew $*"
    (cd "$ANDROID_DIR" && ./gradlew --console=plain "$@")
}

VERSION_TAG="$(date +%Y%m%d-%H%M)"

publish_apk() {
    local variant="$1"      # debug | release
    local src_dir="$ANDROID_DIR/app/build/outputs/apk/$variant"
    local src_file=""
    if [[ -f "$src_dir/app-$variant.apk" ]]; then
        src_file="$src_dir/app-$variant.apk"
    elif [[ -f "$src_dir/app-${variant}-unsigned.apk" ]]; then
        src_file="$src_dir/app-${variant}-unsigned.apk"
        warn "产物未签名: $src_file"
    else
        die "找不到 APK: $src_dir/app-*.apk"
    fi
    local dst="$DIST_DIR/sxwnl-android-$variant.apk"
    cp -f "$src_file" "$dst"
    ok "APK -> $dst ($(du -h "$dst" | cut -f1))"
}

publish_aab() {
    local src="$ANDROID_DIR/app/build/outputs/bundle/release/app-release.aab"
    [[ -f "$src" ]] || die "找不到 AAB: $src"
    local dst="$DIST_DIR/sxwnl-android-release.aab"
    cp -f "$src" "$dst"
    ok "AAB -> $dst ($(du -h "$dst" | cut -f1))"
}

# ── 执行 ──
echo "╔═══════════════════════════════════════╗"
echo "║  寿星万年历 - Android 打包 ($TARGET)"
echo "╚═══════════════════════════════════════╝"

clean_all

case "$TARGET" in
    debug)
        run_gradle "assembleDebug"
        publish_apk debug
        ;;
    release)
        check_signing
        run_gradle "assembleRelease"
        publish_apk release
        ;;
    both)
        # 本地分发常用: 同一次 gradle 调用里串联 debug + release, 共享
        # 一次 configure/JNI/资源处理, 比分两次跑明显省时.
        check_signing
        run_gradle assembleDebug assembleRelease
        publish_apk debug
        publish_apk release
        ;;
    aab)
        check_signing
        run_gradle bundleRelease
        publish_aab
        ;;
    all)
        check_signing
        run_gradle assembleDebug assembleRelease bundleRelease
        publish_apk debug
        publish_apk release
        publish_aab
        ;;
esac

ok "Android 打包完成 (dist/ 目录)"
