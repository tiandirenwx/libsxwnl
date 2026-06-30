#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# 寿星万年历 — Android 正式签名 keystore 生成脚本
#
# 用法:
#   scripts/gen_android_keystore.sh                  # 交互式生成
#   scripts/gen_android_keystore.sh --non-interactive  # 用默认值生成 (需先设环境变量)
#
# 产物:
#   <repo>/android/keystore/sxwnl-release.jks   # keystore 文件 (已 .gitignore)
#   <repo>/android/keystore.properties          # gradle 读取的配置 (已 .gitignore)
#
# 重要:
#   - 生成后请把 .jks 文件和密码备份到安全的地方 (1Password / 加密硬盘 / 离线存储)
#   - 一旦丢失, 你的 App 将无法再发布升级 (用户必须卸载重装)
#   - .jks 和 keystore.properties 都不会被提交到 git

# CI 环境怎么用
# 不要把 keystore.properties 提交进 git，CI 上改用环境变量（已支持）：

# export SXWNL_KEYSTORE_PATH=/path/to/sxwnl-release.jks
# export SXWNL_KEYSTORE_PASSWORD=xxx
# export SXWNL_KEY_ALIAS=sxwnl
# export SXWNL_KEY_PASSWORD=xxx
# scripts/build_android.sh release
# 或者非交互式生成：
# SXWNL_KEYSTORE_PASSWORD=xxx scripts/gen_android_keystore.sh --non-interactive
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

RED='\033[1;31m'; GREEN='\033[1;32m'; YELLOW='\033[1;33m'; BLUE='\033[1;34m'; NC='\033[0m'
info()  { printf "${BLUE}[INFO]${NC} %s\n" "$*"; }
ok()    { printf "${GREEN}[ OK ]${NC} %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC} %s\n" "$*"; }
err()   { printf "${RED}[ERR ]${NC} %s\n" "$*" >&2; }
die()   { err "$*"; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
ANDROID_DIR="$REPO_ROOT/android"
KS_DIR="$ANDROID_DIR/keystore"
KS_FILE="$KS_DIR/sxwnl-release.jks"
KS_PROPS="$ANDROID_DIR/keystore.properties"

command -v keytool >/dev/null 2>&1 \
    || die "keytool 未找到, 请先安装 JDK 17+ 并把 keytool 加入 PATH"

# ── 已存在保护 ──
if [[ -f "$KS_FILE" ]]; then
    warn "keystore 已存在: $KS_FILE"
    read -r -p "是否覆盖? 覆盖后旧签名将永久丢失, 已发布的 App 无法再升级 [y/N]: " ans
    [[ "${ans:-N}" =~ ^[Yy]$ ]] || die "取消"
    rm -f "$KS_FILE"
fi
mkdir -p "$KS_DIR"

# ── 交互式收集信息 ──
NON_INTERACTIVE="${1:-}"

if [[ "$NON_INTERACTIVE" == "--non-interactive" ]]; then
    KEY_ALIAS="${SXWNL_KEY_ALIAS:-sxwnl}"
    STORE_PASSWORD="${SXWNL_KEYSTORE_PASSWORD:?需要环境变量 SXWNL_KEYSTORE_PASSWORD}"
    KEY_PASSWORD="${SXWNL_KEY_PASSWORD:-$STORE_PASSWORD}"
    VALIDITY_DAYS="${SXWNL_VALIDITY:-10950}"  # 30 年
    CN="${SXWNL_DNAME_CN:-Shouxing Wannianli}"
    OU="${SXWNL_DNAME_OU:-Dev}"
    O="${SXWNL_DNAME_O:-Sxwnl}"
    L="${SXWNL_DNAME_L:-Unknown}"
    ST="${SXWNL_DNAME_ST:-Unknown}"
    C="${SXWNL_DNAME_C:-CN}"
else
    echo "═══ 生成 Android 正式签名 keystore ═══"
    echo "(直接回车使用方括号内默认值)"
    echo ""

    read -r -p "Key 别名 [sxwnl]: " KEY_ALIAS
    KEY_ALIAS="${KEY_ALIAS:-sxwnl}"

    while true; do
        read -r -s -p "keystore 密码 (至少 6 位): " STORE_PASSWORD; echo
        [[ ${#STORE_PASSWORD} -ge 6 ]] || { warn "至少 6 位, 重输"; continue; }
        read -r -s -p "再输入一次确认: " STORE_PASSWORD2; echo
        [[ "$STORE_PASSWORD" == "$STORE_PASSWORD2" ]] || { warn "两次不一致, 重输"; continue; }
        break
    done

    read -r -p "key 密码 (回车则与 keystore 密码相同): " -s KEY_PASSWORD; echo
    KEY_PASSWORD="${KEY_PASSWORD:-$STORE_PASSWORD}"
    if [[ "$KEY_PASSWORD" != "$STORE_PASSWORD" ]]; then
        read -r -s -p "再输入一次 key 密码确认: " KEY_PASSWORD2; echo
        [[ "$KEY_PASSWORD" == "$KEY_PASSWORD2" ]] || die "两次 key 密码不一致"
    fi

    read -r -p "有效期 (天) [10950 ≈ 30 年]: " VALIDITY_DAYS
    VALIDITY_DAYS="${VALIDITY_DAYS:-10950}"

    echo ""
    echo "─── 证书 DName (可随便填, 但建议填真实信息便于识别) ───"
    read -r -p "CN  名字/常用名     [Shouxing Wannianli]: " CN; CN="${CN:-Shouxing Wannianli}"
    read -r -p "OU  部门            [Dev]: "                OU; OU="${OU:-Dev}"
    read -r -p "O   组织            [Sxwnl]: "              O;  O="${O:-Sxwnl}"
    read -r -p "L   城市            [Unknown]: "            L;  L="${L:-Unknown}"
    read -r -p "ST  省份            [Unknown]: "            ST; ST="${ST:-Unknown}"
    read -r -p "C   国家代码 (2 位) [CN]: "                 C;  C="${C:-CN}"
fi

# ── 调用 keytool 生成 ──
info "生成 keystore: $KS_FILE"
keytool -genkeypair -v \
    -keystore "$KS_FILE" \
    -storetype PKCS12 \
    -alias "$KEY_ALIAS" \
    -keyalg RSA -keysize 2048 \
    -validity "$VALIDITY_DAYS" \
    -storepass "$STORE_PASSWORD" \
    -keypass "$KEY_PASSWORD" \
    -dname "CN=$CN, OU=$OU, O=$O, L=$L, ST=$ST, C=$C" \
    >/dev/null

ok "keystore 已生成"

# ── 写 keystore.properties ──
# 路径以 rootProject (android/) 为基准
REL_PATH="keystore/sxwnl-release.jks"
cat > "$KS_PROPS" <<EOF
# 自动生成于 $(date "+%Y-%m-%d %H:%M:%S")
# !! 此文件包含密码, 已加入 .gitignore, 不要提交到 git !!
storeFile=$REL_PATH
storePassword=$STORE_PASSWORD
keyAlias=$KEY_ALIAS
keyPassword=$KEY_PASSWORD
EOF
chmod 600 "$KS_PROPS"
ok "配置已写入: $KS_PROPS (权限 600)"

# ── 摘要 ──
echo ""
info "证书指纹:"
keytool -list -v -keystore "$KS_FILE" -storepass "$STORE_PASSWORD" -alias "$KEY_ALIAS" \
    2>/dev/null | grep -E "SHA1|SHA256|MD5|有效期|Valid" | sed 's/^/        /'

echo ""
ok "完成! 现在可以执行: scripts/build_android.sh release"
echo ""
warn "请立即把以下两个文件备份到安全的地方:"
warn "  1) $KS_FILE"
warn "  2) $KS_PROPS  (或仅记下密码)"
warn "丢失后 App 将无法升级!"
