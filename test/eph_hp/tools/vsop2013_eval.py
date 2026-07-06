#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vsop2013_eval.py
================
VSOP2013 参考评估器 (Python 实现, 对标 IMCCE 官方 VSOP2013.f).

**用途**:
    在把 VSOP2013 集成到 C++ 之前, 用 Python 先跑参考值, 校验数据完整性,
    并产出对比表. 供将来的 C++ 实现做单元测试 ground truth.

**数据源**:
    - data/VSOP2013p3.dat  : 地球-月亮质心 (EMB) 轨道要素完整级数
    - data/VSOP2013.f      : IMCCE 官方 Fortran 参考实现 (校对用)
    - data/VSOP2013.ctl    : 官方精度/范围参数

**输入格式** (每 section 起始行, `format (9x,3i3,i7)`):
    "VSOP2013  ip iv it  nt"    表示第 ip 号天体的第 iv 变量, 时间幂 it, 共 nt 项
    ip: 1..9 = Me, Ve, EMB, Ma, J, Sa, U, Ne, Pl
    iv: 1..6 = A, L, K, H, Q, P (Poincaré 椭圆根数)
    it: 时间幂 t^it, 单位为儒略千年 (a1000 = 365250 天)

**每项行** (`format (6x,4i3,1x,5i3,1x,4i4,1x,i6,1x,3i3,2a24)`):
    17 个整数乘子 + 2 个双精度浮点 (S 正弦项, K 余弦项)
    17 个基本参数为: 9 大行星 + 5 小行星/冥王星 + 3 月球 Delaunay (D, F, l)

**级数计算公式** (VSOP2013.f L307-333):
    对每个变量 iv = 1..6:
        r[iv] = 0
        对每个时间幂 it = 0..20:
            对每一项 n:
                ARG = Σ_j iphi[j] · ci0[j]  +  t · Σ_j iphi[j] · ci1[j]
                r[iv] += t^it · (S[n]·sin(ARG) + K[n]·cos(ARG))
    // L 分量额外补充线性平均运动:
    r[2] = mod( r[2] + freqpla[ip] * t,  2π )

**从椭圆根数到位置** (未在本文件实现, 作为 TODO 供将来 C++ 扩展):
    (a, λ, k, h, q, p) → 解 Kepler → (x, y, z)_orbital → 旋转 → (X, Y, Z)_ecliptic

用法:
    python3 vsop2013_eval.py data/VSOP2013p3.dat  --verify
        只解析并验证格式.
    python3 vsop2013_eval.py data/VSOP2013p3.dat  --at 2020.5 2021.0
        评估某几年的 6 个 Poincaré 根数, 输出到 stdout.
    python3 vsop2013_eval.py data/VSOP2013p3.dat  --csv out.csv --range -1000 3000 500
        批量输出参考值 CSV, 供将来 C++ 测试对照.
"""
import sys
import math
import argparse
from pathlib import Path


# --- 17 个基本参数 (与 VSOP2013.f 的 ci0/ci1 完全一致, radians / rad·cy^-1) ---
CI0 = [
    0.4402608631669e1,   # 1  Mercury
    0.3176134461576e1,   # 2  Venus
    0.1753470369433e1,   # 3  Earth-Moon Barycenter
    0.6203500014141e1,   # 4  Mars
    0.4091360003050e1,   # 5  Vesta
    0.1713740719173e1,   # 6  Iris
    0.5598641292287e1,   # 7  Bamberga
    0.2805136360408e1,   # 8  Ceres
    0.2326989734620e1,   # 9  Pallas
    0.5995461070350e0,   # 10 Jupiter
    0.8740185101070e0,   # 11 Saturn
    0.5481225395663e1,   # 12 Uranus
    0.5311897933164e1,   # 13 Neptune
    0.0,                 # 14 (Pluto reserved, ci0 slot=0)
    5.19846640063,       # 15 Moon D
    1.62790513602,       # 16 Moon F
    2.35555563875,       # 17 Moon l
]

CI1 = [
    0.2608790314068555e5,   # Mercury mean motion (rad / kyear)
    0.1021328554743445e5,   # Venus
    0.6283075850353215e4,   # Earth-Moon Barycenter
    0.3340612434145457e4,   # Mars
    0.1731170452721855e4,   # Vesta
    0.1704450855027201e4,   # Iris
    0.1428948917844273e4,   # Bamberga
    0.1364756513629990e4,   # Ceres
    0.1361923207632842e4,   # Pallas
    0.5296909615623250e3,   # Jupiter
    0.2132990861084880e3,   # Saturn
    0.7478165903077800e2,   # Uranus
    0.3813297222612500e2,   # Neptune
    0.3595362285049309e0,   # Pluto
    77713.7714481804,       # Moon D
    84334.6615717837,       # Moon F
    83286.9142477147,       # Moon l
]

# 用于 L 分量的额外线性项 (freqpla)
FREQPLA = {
    1: 0.2608790314068555e5,   # Mercury
    2: 0.1021328554743445e5,   # Venus
    3: 0.6283075850353215e4,   # EMB
    4: 0.3340612434145457e4,   # Mars
    5: 0.5296909615623250e3,   # Jupiter
    6: 0.2132990861084880e3,   # Saturn
    7: 0.7478165903077800e2,   # Uranus
    8: 0.3813297222612500e2,   # Neptune
    9: 0.2533566020437000e2,   # Pluto
}

VAR_NAMES = {1: "A", 2: "L", 3: "K", 4: "H", 5: "Q", 6: "P"}
DPI = 2.0 * math.pi
A1000 = 365250.0  # 天/千年
J2000_JD = 2451545.0

TWO_PI = 2.0 * math.pi


def parse_file(path: Path):
    """
    解析 VSOP2013p<ip>.dat 文件. 返回:
        header_info: dict{ (iv, it): term_count }
        terms:       list of dicts { 'iv', 'it', 'iphi': [17 int], 's', 'k' }
    """
    header_info = {}
    terms = []
    current_iv = None
    current_it = None
    remaining = 0
    with open(path, "r", encoding="ascii") as fp:
        for line in fp:
            if line.strip().startswith("VSOP2013"):
                # e.g. " VSOP2013  3  1  0  32658    EARTH-MOON VARIABLE A   *T*00"
                # Fortran format 1001: (9x, 3i3, i7)
                # 前 9 列跳过, 接着 3 个 i3 是 ip iv it, 再 i7 是 nt
                ip = int(line[9:12])
                iv = int(line[12:15])
                it = int(line[15:18])
                nt = int(line[18:25])
                current_iv, current_it, remaining = iv, it, nt
                header_info[(iv, it)] = nt
                continue
            if current_iv is None or remaining == 0:
                continue
            # 数据行 format 1002:
            #   (6x, 4i3, 1x, 5i3, 1x, 4i4, 1x, i6, 1x, 3i3, 2a24)
            # 用固定列宽解析(注意 Fortran D 记法 => 手工用 e 记法解析)
            #   位置: 6-17(4i3), 18(空), 19-33(5i3), 34(空), 35-50(4i4), 51(空),
            #         52-57(i6), 58(空), 59-67(3i3), 68-91(a24), 92-115(a24)
            try:
                # 17 个整数乘子, 按 Fortran 定长切
                iphi = []
                # 4 i3
                for k in range(4):
                    iphi.append(int(line[6 + k * 3 : 9 + k * 3]))
                # 5 i3
                for k in range(5):
                    iphi.append(int(line[19 + k * 3 : 22 + k * 3]))
                # 4 i4
                for k in range(4):
                    iphi.append(int(line[35 + k * 4 : 39 + k * 4]))
                # 1 i6
                iphi.append(int(line[52:58]))
                # 3 i3
                for k in range(3):
                    iphi.append(int(line[59 + k * 3 : 62 + k * 3]))
                assert len(iphi) == 17

                # S 和 K: 24 字符各一个, Fortran D 记法 (e.g. "  0.1000001017641000 +01")
                # 转换 " D" 或 " +" 为标准 e 记法
                s_raw = line[68:92]
                k_raw = line[92:116] if len(line) >= 116 else line[92:].rstrip("\n")
                s_val = _fortran_d(s_raw)
                k_val = _fortran_d(k_raw)
            except (ValueError, IndexError):
                continue
            terms.append(
                {"iv": current_iv, "it": current_it, "iphi": iphi, "s": s_val, "k": k_val}
            )
            remaining -= 1
    return header_info, terms


def _fortran_d(s: str) -> float:
    """
    Fortran double: '  0.1234567890123456 +01'  → 0.1234567890123456e+01
    分隔符可能为 ' +' / ' -' / 'D+' / 'D-'.
    """
    s = s.strip()
    for marker in (" D+", " D-", " +", " -", "D+", "D-"):
        i = s.find(marker)
        if i > 0:
            mant = s[:i]
            expo = s[i + 1 :].replace("D", "").replace(" ", "")
            return float(mant) * (10.0 ** int(expo))
    return float(s)


def evaluate(terms, ip: int, jd_tdb: float):
    """
    在给定 TDB 儒略日下求 (A, L, K, H, Q, P) 6 个 Poincaré 根数.
    对应 VSOP2013.f 的 subroutine 主体.
    """
    t1 = (jd_tdb - J2000_JD) / A1000  # 千年 (从 J2000 起)
    tn = [1.0] + [0.0] * 20
    for i in range(1, 21):
        tn[i] = tn[i - 1] * t1

    r = [0.0] * 7  # 1..6 有效
    for term in terms:
        iv = term["iv"]
        it = term["it"]
        aa = 0.0
        bb = 0.0
        for j in range(17):
            n = term["iphi"][j]
            if n != 0:
                aa += n * CI0[j]
                bb += n * CI1[j]
        arg = aa + bb * t1
        r[iv] += tn[it] * (term["s"] * math.sin(arg) + term["k"] * math.cos(arg))

    # L 分量补充线性平均运动
    r[2] = (r[2] + FREQPLA[ip] * t1) % TWO_PI
    if r[2] < 0:
        r[2] += TWO_PI

    return r[1:7]  # 返回 [A, L, K, H, Q, P]


def year_to_jd_tdb(year: float) -> float:
    """
    近似: 把年份中点 (以 Julian 年计) 转 JD (TDB).
    对精度不严的比较测试足够(误差 ~1 天,不影响千年尺度对比).
    """
    return J2000_JD + (year - 2000.0) * 365.25


def cmd_verify(path: Path):
    print(f"验证 {path.name} ...")
    hdr, terms = parse_file(path)
    total = 0
    for (iv, it), nt in sorted(hdr.items()):
        got = sum(1 for t in terms if t["iv"] == iv and t["it"] == it)
        status = "OK" if got == nt else f"MISMATCH ({got}/{nt})"
        print(f"  iv={iv}({VAR_NAMES[iv]}) it={it:>2d}   expected={nt:>6d}  got={got:>6d}  {status}")
        total += nt
    print(f"总项数: {total}")
    return 0 if len(terms) == total else 1


def cmd_at(path: Path, years, ip: int = 3):
    hdr, terms = parse_file(path)
    print(f"# 天体 ip={ip}, VSOP2013 Poincaré 根数 (a, λ, k, h, q, p)")
    print(f"# {'Year':>8}  {'a (AU)':>16}  {'λ (rad)':>18}  {'k':>14}  {'h':>14}  {'q':>14}  {'p':>14}")
    for y in years:
        jd = year_to_jd_tdb(y)
        A, L, K, H, Q, P = evaluate(terms, ip, jd)
        print(f"  {y:>8.2f}  {A:>16.10f}  {L:>18.10f}  {K:>14.10e}  {H:>14.10e}  {Q:>14.10e}  {P:>14.10e}")
    return 0


def cmd_csv(path: Path, csv_path: Path, start: int, end: int, step: int, ip: int = 3):
    hdr, terms = parse_file(path)
    with open(csv_path, "w") as f:
        f.write("year,jd_tdb,A,L,K,H,Q,P\n")
        y = start
        while y <= end:
            jd = year_to_jd_tdb(y)
            A, L, K, H, Q, P = evaluate(terms, ip, jd)
            f.write(f"{y},{jd},{A},{L},{K},{H},{Q},{P}\n")
            y += step
    print(f"已写出 {csv_path}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_file", help="VSOP2013p3.dat 等文件路径")
    ap.add_argument("--verify", action="store_true", help="只验证格式解析")
    ap.add_argument("--at", nargs="+", type=float, help="求给定年份的 6 根数")
    ap.add_argument("--csv", help="输出 CSV")
    ap.add_argument("--range", nargs=3, type=int, help="年份范围 START END STEP")
    ap.add_argument("--ip", type=int, default=3, help="天体索引 (3=EMB)")
    args = ap.parse_args()

    p = Path(args.data_file)
    if not p.exists():
        print(f"文件不存在: {p}", file=sys.stderr)
        print(
            "提示: VSOP2013p3.dat 是 33 MB 大文件, 不随仓库分发.\n"
            "      运行 test/eph_hp/fetch_data.sh vsop2013 一键下载.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.verify:
        sys.exit(cmd_verify(p))
    if args.at:
        sys.exit(cmd_at(p, args.at, args.ip))
    if args.csv:
        if not args.range:
            print("--csv 需要配合 --range START END STEP", file=sys.stderr)
            sys.exit(2)
        start, end, step = args.range
        sys.exit(cmd_csv(p, Path(args.csv), start, end, step, args.ip))

    print(__doc__)


if __name__ == "__main__":
    main()
