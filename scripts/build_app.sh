#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# 寿星万年历 — 一键打包总入口 (Android + iOS)
#
# 用法:
#   scripts/build_app.sh                    # 默认: android debug
#   scripts/build_app.sh android [...]      # 仅 Android, 后续参数转给 build_android.sh
#   scripts/build_app.sh ios [...]          # 仅 iOS,    后续参数转给 build_ios.sh
#   scripts/build_app.sh all                # android release + iOS sim + macOS
#
# 示例:
#   scripts/build_app.sh android release    # Android release APK
#   scripts/build_app.sh android aab        # Android Play Store AAB
#   scripts/build_app.sh ios sim            # iOS 模拟器 .app
#   scripts/build_app.sh ios device         # iOS 真机 IPA
#   scripts/build_app.sh ios mac            # macOS .app
#
# 产物统一输出到 <repo>/dist/
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── 颜色 ──
GREEN='\033[1;32m'; YELLOW='\033[1;33m'; BLUE='\033[1;34m'; NC='\033[0m'
info() { printf "${BLUE}[INFO]${NC} %s\n" "$*"; }
ok()   { printf "${GREEN}[ OK ]${NC} %s\n" "$*"; }
warn() { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }

usage() { sed -n '2,18p' "$0"; }

PLATFORM="${1:-default}"
shift || true

case "$PLATFORM" in
    -h|--help|help)
        usage; exit 0 ;;
    android)
        bash "$SCRIPT_DIR/build_android.sh" "$@"
        ;;
    ios)
        bash "$SCRIPT_DIR/build_ios.sh" "$@"
        ;;
    all)
        info "═══ 全平台打包 ═══"
        bash "$SCRIPT_DIR/build_android.sh" release || warn "Android release 失败, 继续"
        bash "$SCRIPT_DIR/build_ios.sh"      sim    || warn "iOS sim 失败, 继续"
        if [[ "${SXWNL_BUILD_DEVICE:-0}" == "1" ]]; then
            bash "$SCRIPT_DIR/build_ios.sh"  device || warn "iOS device 失败, 继续"
        else
            warn "已跳过 iOS device 包 (需要开发者签名, 设置 SXWNL_BUILD_DEVICE=1 开启)"
        fi
        bash "$SCRIPT_DIR/build_ios.sh"      mac    || warn "macOS 失败, 继续"
        ;;
    default)
        info "═══ 默认: 仅 Android debug (跨平台请显式指定 ios/all) ═══"
        bash "$SCRIPT_DIR/build_android.sh" debug
        ;;
    *)
        warn "未知平台: $PLATFORM"; usage; exit 1 ;;
esac

echo
ok "全部完成, 产物位置: $REPO_ROOT/dist/"
if [[ -d "$REPO_ROOT/dist" ]]; then
    ls -lh "$REPO_ROOT/dist/" | tail -n +2 | awk '{printf "      %s  %s\n", $5, $NF}'
fi
