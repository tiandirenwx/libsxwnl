#!/usr/bin/env bash
# install_to_device.sh
# 一键把寿星万年历鸿蒙版构建并安装到已连接的真机/模拟器
# 前置条件：
#   1. 已在 DevEco Studio 中配置好 signingConfigs（Project Structure → Signing Configs → Automatically generate signature）
#   2. 已用 USB 连接 HarmonyOS NEXT 真机，并启用 USB 调试
#   3. macOS 上安装了 DevEco Studio

set -euo pipefail

# ---------- 配置 ----------
BUNDLE_NAME="com.sxwnl.calendar"
ABILITY_NAME="EntryAbility"
MODULE_NAME="entry"
PRODUCT_NAME="default"

DEVECO_APP="/Applications/DevEco-Studio.app"
HDC="${DEVECO_APP}/Contents/sdk/default/openharmony/toolchains/hdc"
HVIGORW="${DEVECO_APP}/Contents/tools/hvigor/bin/hvigorw"
NODE_HOME="${DEVECO_APP}/Contents/tools/node"
# hvigor 必须的 SDK 根目录（下面应包含 default/openharmony 和 default/hms）
DEVECO_SDK_HOME="${DEVECO_SDK_HOME:-${DEVECO_APP}/Contents/sdk}"

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HAP_OUT="${PROJECT_ROOT}/${MODULE_NAME}/build/${PRODUCT_NAME}/outputs/${PRODUCT_NAME}/${MODULE_NAME}-${PRODUCT_NAME}-signed.hap"
HAP_UNSIGNED="${PROJECT_ROOT}/${MODULE_NAME}/build/${PRODUCT_NAME}/outputs/${PRODUCT_NAME}/${MODULE_NAME}-${PRODUCT_NAME}-unsigned.hap"

# ---------- 颜色输出 ----------
RED='\033[1;31m'; GREEN='\033[1;32m'; YELLOW='\033[1;33m'; BLUE='\033[1;34m'; NC='\033[0m'
info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
ok()    { echo -e "${GREEN}[ OK ]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
err()   { echo -e "${RED}[ERR ]${NC} $*" >&2; }
die()   { err "$*"; exit 1; }

# ---------- 0. 环境检查 ----------
info "项目根目录: ${PROJECT_ROOT}"
[[ -x "$HDC" ]]     || die "找不到 hdc，请确认已安装 DevEco Studio: ${HDC}"
[[ -x "$HVIGORW" ]] || die "找不到 hvigorw: ${HVIGORW}"
[[ -d "$NODE_HOME" ]] || die "找不到 DevEco 自带 node: ${NODE_HOME}"
[[ -d "$DEVECO_SDK_HOME/default/openharmony" ]] \
  || die "DEVECO_SDK_HOME 无效: ${DEVECO_SDK_HOME}（期望其下存在 default/openharmony）"
info "DEVECO_SDK_HOME = ${DEVECO_SDK_HOME}"

# ---------- 1. 检查设备 ----------
info "检测已连接设备..."
DEVICES=$("$HDC" list targets 2>/dev/null | grep -v "^\[Empty\]" | grep -v "^$" || true)
if [[ -z "$DEVICES" ]]; then
  err "没有检测到任何设备"
  echo "请确认："
  echo "  1) 数据线已连接手机"
  echo "  2) 手机已开启 USB 调试（设置 → 系统和更新 → 开发人员选项 → USB 调试）"
  echo "  3) 手机上已点击 \"允许此电脑调试\""
  exit 1
fi
DEVICE_COUNT=$(echo "$DEVICES" | wc -l | tr -d ' ')
ok "检测到 ${DEVICE_COUNT} 台设备："
echo "$DEVICES" | sed 's/^/      /'

# 设备选择策略：
#   1) 若环境变量 HDC_TARGET 指定了，则用它
#   2) 否则优先选真机（非 127.0.0.1 / 非 localhost 开头的设备）
#   3) 没真机才回落到模拟器
if [[ -n "${HDC_TARGET:-}" ]]; then
  TARGET_DEVICE="$HDC_TARGET"
  info "通过 HDC_TARGET 指定设备: ${TARGET_DEVICE}"
else
  REAL_DEVICE=$(echo "$DEVICES" | grep -vE "^(127\.0\.0\.1:|localhost:)" | head -1 || true)
  if [[ -n "$REAL_DEVICE" ]]; then
    TARGET_DEVICE="$REAL_DEVICE"
    ok "自动选中真机: ${TARGET_DEVICE}"
  else
    TARGET_DEVICE=$(echo "$DEVICES" | head -1)
    warn "未发现真机，回落到: ${TARGET_DEVICE}（模拟器/远程通道）"
    read -r -p "继续部署到模拟器吗？(y/N) " yn
    [[ "$yn" =~ ^[Yy]$ ]] || exit 1
  fi
fi
HDC_RUN=("$HDC" -t "$TARGET_DEVICE")
info "使用设备: ${TARGET_DEVICE}"

# ---------- 2. 检查签名配置 ----------
info "检查签名配置..."
if ! grep -q "\"signingConfigs\":\s*\[" "${PROJECT_ROOT}/build-profile.json5" 2>/dev/null \
   || grep -Pzoq "\"signingConfigs\"\s*:\s*\[\s*\]" "${PROJECT_ROOT}/build-profile.json5" 2>/dev/null; then
  warn "harmony/build-profile.json5 中的 signingConfigs 看起来还是空的。"
  warn "请先在 DevEco Studio 完成自动签名："
  warn "  File → Project Structure → Project → Signing Configs"
  warn "  勾选 \"Automatically generate signature\" → Sign In"
  warn "完成后再次运行此脚本。"
  echo
  read -r -p "是否继续尝试构建？(y/N) " yn
  [[ "$yn" =~ ^[Yy]$ ]] || exit 1
fi

# ---------- 3. 构建 ----------
info "调用 hvigorw 构建 hap (assembleHap, mode=module, product=${PRODUCT_NAME})..."
export NODE_HOME
export DEVECO_SDK_HOME
# 把 DevEco 自带 node 放到 PATH 最前，避免系统里其它 node 版本干扰
export PATH="${NODE_HOME}/bin:${PATH}"
cd "$PROJECT_ROOT"

# ─── 3.1 硬清理：彻底干掉旧构建产物 + CMake 缓存 + hvigor 任务缓存 ───
#
# 为什么不只用 `hvigorw clean`：
#   - hvigorw clean 对 hap/ts 增量缓存有效, 但对 native cpp 的 CMake
#     CACHE/Ninja .ninja_deps 经常清不干净, 改了 .cpp/.h 后 .so 不重编
#   - hvigor 自身的任务级缓存 (.hvigor) 也会让某些 ts 增量编译错过 .ets 改动
#
# 由 SKIP_CLEAN=1 可跳过 (例如调试脚本时反复跑省时间).
if [[ "${SKIP_CLEAN:-0}" == "1" ]]; then
  warn "SKIP_CLEAN=1, 跳过硬清理 (仅用于调试)"
else
  info "硬清理旧构建产物 (设 SKIP_CLEAN=1 可跳过)..."
  # hvigorw clean 先跑一次, 让它正常释放自己持有的 lock
  "$HVIGORW" clean --no-daemon 2>/dev/null || warn "hvigorw clean 失败, 继续物理删除"
  CLEAN_PATHS=(
    "${PROJECT_ROOT}/build"
    "${PROJECT_ROOT}/.hvigor"
    "${PROJECT_ROOT}/.cxx"
    "${PROJECT_ROOT}/${MODULE_NAME}/build"
    "${PROJECT_ROOT}/${MODULE_NAME}/.hvigor"
    "${PROJECT_ROOT}/${MODULE_NAME}/.cxx"
    "${PROJECT_ROOT}/${MODULE_NAME}/src/main/cpp/build"
  )
  for p in "${CLEAN_PATHS[@]}"; do
    if [[ -e "$p" ]]; then
      rm -rf "$p" && info "  rm -rf ${p#$PROJECT_ROOT/}"
    fi
  done
  ok "硬清理完成"
fi

"$HVIGORW" assembleHap --mode module -p product=${PRODUCT_NAME} --no-daemon

# ---------- 4. 检查产物 ----------
if [[ -f "$HAP_OUT" ]]; then
  ok "构建成功: ${HAP_OUT}"
elif [[ -f "$HAP_UNSIGNED" ]]; then
  die "只构建出未签名 hap (${HAP_UNSIGNED})，HarmonyOS NEXT 真机无法安装未签名包。请先在 DevEco Studio 配置 Signing Configs。"
else
  die "构建产物未生成，请检查上方 hvigorw 输出"
fi

# ---------- 5. 卸载旧版（忽略未安装错误） ----------
info "卸载旧版（如果有）..."
"${HDC_RUN[@]}" uninstall "$BUNDLE_NAME" 2>&1 | sed 's/^/      /' || true

# ---------- 6. 安装 ----------
info "安装 hap 到设备..."
"${HDC_RUN[@]}" install -r "$HAP_OUT"

# ---------- 7. 启动 Ability ----------
info "启动 ${BUNDLE_NAME}/${ABILITY_NAME} ..."
"${HDC_RUN[@]}" shell aa start -a "$ABILITY_NAME" -b "$BUNDLE_NAME" || warn "自动启动失败，请手动在手机桌面点击应用图标"

ok "完成！请在手机上查看 \"寿星万年历\"。"
