#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# 八字页共享资源同步 — 将 assets/bazi/ 链接到三端资源目录
#
# 用法:
#   scripts/sync_bazi_assets.sh           # 同步 android + harmony + ios
#   scripts/sync_bazi_assets.sh android
#   scripts/sync_bazi_assets.sh harmony
#   scripts/sync_bazi_assets.sh ios
#
# 各端在编译前自动调用 (无需手动执行):
#   Android  — app/build.gradle.kts preBuild
#   iOS      — setup_xcode.sh --link-only / Xcode Run Script
#   Harmony  — entry/hvigorfile.ts 加载时 (复制, 因 hvigor 不识别符号链接)
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
SHARED="$ROOT_DIR/assets/bazi"

[[ -d "$SHARED" ]] || {
    echo "[sync_bazi_assets] 缺少目录: $SHARED" >&2
    exit 1
}

link_file() {
    local target_rel="$1"
    local link_path="$2"
    mkdir -p "$(dirname "$link_path")"
    ln -sf "$target_rel" "$link_path"
}

copy_file() {
    local src="$1"
    local dest="$2"
    mkdir -p "$(dirname "$dest")"
    rm -f "$dest"
    cp -f "$src" "$dest"
}

sync_android() {
    local drawable="$ROOT_DIR/android/app/src/main/res/drawable"
    local font="$ROOT_DIR/android/app/src/main/res/font"
    local rel="../../../../../../assets/bazi"

    for f in "$SHARED"/bz_*.png "$SHARED"/bz_paper.jpg; do
        [[ -e "$f" ]] || continue
        link_file "$rel/$(basename "$f")" "$drawable/$(basename "$f")"
    done
    link_file "$rel/WenYue.otf" "$font/wen_yue.otf"
    echo "  ✓ Android 共享资源已链接"
}

sync_harmony() {
    local media="$ROOT_DIR/harmony/entry/src/main/resources/base/media"
    local rawfile="$ROOT_DIR/harmony/entry/src/main/resources/rawfile"

    mkdir -p "$rawfile" "$media"
    # hvigor 资源打包不跟随符号链接, 必须复制实体文件
    for f in "$SHARED"/bz_*.png "$SHARED"/bz_paper.jpg; do
        [[ -e "$f" ]] || continue
        copy_file "$f" "$media/$(basename "$f")"
    done
    copy_file "$SHARED/WenYue.otf" "$rawfile/WenYue.otf"
    echo "  ✓ Harmony 共享资源已复制"
}

sync_ios() {
    local fonts="$ROOT_DIR/ios/SxwnlCalendar/SxwnlCalendar/Fonts"
    local xcassets="$ROOT_DIR/ios/SxwnlCalendar/SxwnlCalendar/Assets.xcassets"
    local rel_font="../../../../assets/bazi"
    local rel_img="../../../../../../assets/bazi"

    mkdir -p "$fonts"
    link_file "$rel_font/WenYue.otf" "$fonts/WenYue.otf"

    for f in "$SHARED"/bz_*.png; do
        [[ -e "$f" ]] || continue
        local name
        name="$(basename "$f" .png)"
        link_file "$rel_img/$(basename "$f")" "$xcassets/${name}.imageset/$(basename "$f")"
    done
    link_file "$rel_img/bz_paper.jpg" "$xcassets/bz_paper.imageset/bz_paper.jpg"
    echo "  ✓ iOS 共享资源已链接"
}

TARGET="${1:-all}"
case "$TARGET" in
    android) sync_android ;;
    ios) sync_ios ;;
    harmony) sync_harmony ;;
    all)
        sync_android
        sync_harmony
        sync_ios
        ;;
    *)
        echo "用法: $0 [android|ios|harmony|all]" >&2
        exit 1
        ;;
esac
