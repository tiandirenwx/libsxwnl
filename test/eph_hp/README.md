# test/eph_hp/ — 高精度天体位置计算"备件"

本目录存放**将来升级天文精度**用的原始数据获取脚本、转换工具、参考实现和对比测试。

**核心约束**:
1. 此目录所有代码/数据**不进入 `libsxwnl` 主库**,只在测试时编译。主库 `src/` 零改动。
2. **只入库脚本 + 源码**:所有原始数据、生成的 C++ header、生成的 CSV 都**不进 git**,通过 `fetch_data.sh` 从上游按需下载并生成。

---

## 快速开始(60 秒)

```bash
# 1. 拉数据 + 生成 header (只需一次, ~500 KB, 秒级完成)
cd test/eph_hp && ./fetch_data.sh small

# 2. 从项目根目录构建并运行
cd ../..
cmake -S . -B build            # cmake 自动检测到数据, 启用 eph_hp_test
cmake --build build --target eph_hp_test
./build/bin/eph_hp_test        # 跑 Test A/B/C/D
```

**前置依赖**:`curl` + `python3`(生成 header 用)+ `cmake >= 3.25` + C++20 编译器。

如果 `cmake configure` 看到:
```
-- HP: eph_hp_test 已启用 (VSOP87D + IAU 2000A 全项对比)
```
就说明数据到位了。看到"已跳过"就是 `fetch_data.sh` 还没跑。

---

## 一、入库文件清单(git 追踪)

```
test/eph_hp/
├── README.md                       [入库] 本文件
├── CMakeLists.txt                  [入库] 独立测试目标 (含 calceph 可选分支)
├── fetch_data.sh                   [入库] 一键下载 + 生成脚本
│
├── eph_hp.h / eph_hp.cpp           [入库] HP 参考实现 (VSOP87D 全项 + IAU 2000A 章动)
├── eph_hp_test.cpp                 [入库] Test A/B/C/D 测试
├── de_hp.h / de_hp.cpp             [入库] calceph 包装 (DE441/DE440)
├── de_hp_test.cpp                  [入库] Test F 测试 (可选)
│
├── data/
│   └── .gitkeep                    [入库] 占位, 保证 fetch_data.sh 可写入
│
└── tools/                          [入库] 数据转换 / 评估工具 (Python)
    ├── vsop87d_to_header.py        VSOP87D → C++ header
    ├── iau2000a_to_header.py       IAU 2000A → C++ header
    └── vsop2013_eval.py            VSOP2013 Python 解析器 + 评估器
```

## 二、被 `.gitignore` 排除的产出(全部由 `fetch_data.sh` 生成)

| 文件 | 大小 | 来源 | 何时需要 |
|---|---|---|---|
| `data/VSOP87D.ear` | 318 KB | IMCCE | 生成 `vsop87d_earth_hp.h` |
| `data/IAU2000A_tab5.3a.txt` | 142 KB | IERS | 生成 `iau2000a_nut_lon_hp.h` |
| `data/VSOP2013.f` | 15 KB | IMCCE | VSOP2013 官方 Fortran 参考(仅文档用) |
| `data/VSOP2013.ctl` | 35 KB | IMCCE | VSOP2013 参数控制(仅文档用) |
| `data/VSOP2013p3.dat` | 33 MB | IMCCE | Python VSOP2013 评估器 |
| `data/de440s.bsp` | 33 MB | NASA | DE440 星历(Test F 用) |
| `data/de441.bsp` | 3.2 GB | NASA | DE441 完整星历(-13000~17191 CE) |
| `vsop87d_earth_hp.h` | 97 KB | Python 生成 | `eph_hp_test` C++ 编译 |
| `iau2000a_nut_lon_hp.h` | 92 KB | Python 生成 | `eph_hp_test` C++ 编译 |
| `vsop2013_emb_reference.csv` | 3 KB | Python 生成 | 交叉验证参考 |

---

## 三、`fetch_data.sh` 分场景用法

| 场景 | 命令 | 下载 | 生成 |
|---|---|---|---|
| **只跑基础 HP 测试**(Test A/B/C/D) | `./fetch_data.sh small` | 4 个小文件 (~500 KB) | 2 个 header |
| **跑 VSOP2013 Python 评估** | `./fetch_data.sh vsop2013` | +VSOP2013p3.dat (33 MB) | +CSV 参考 |
| **跑 DE440 星历对比**(Test F) | `./fetch_data.sh de440` | +de440s.bsp (33 MB) | — |
| **跑 DE441 长时段**(> 2150 CE) | `./fetch_data.sh de441` | +de441.bsp (3.2 GB) | — |
| **默认全套**(除 DE441) | `./fetch_data.sh` | 全部小 + vsop2013 + de440 | 全部 |
| **只跑生成**(数据已下过) | `./fetch_data.sh gen` | — | 2 个 header + CSV |

脚本特性:
- **幂等**:已存在且大小合理会跳过
- **断点续传**:`curl -C -` 支持
- **自动补生成**:下载完小文件后自动跑 Python 生成 header,无需手动 `gen`

---

## 四、构建流程详解

### 4.1 主 HP 测试(Test A/B/C/D)

```bash
cd test/eph_hp && ./fetch_data.sh small
cd ../..
cmake -S . -B build
cmake --build build --target eph_hp_test
./build/bin/eph_hp_test
```

CMake 智能检测:
- 数据/header 齐全 → 定义 `eph_hp_test` target,打印 `HP: eph_hp_test 已启用`
- 数据/header 缺失 → 跳过 target,打印 `HP: eph_hp_test 已跳过`,主库正常构建
- 原始数据在但 header 未生成 → **configure 阶段自动调用 Python 补生成**,零手工干预

### 4.2 DE440/DE441 星历测试(可选,需 calceph)

```bash
# 安装 calceph 库
brew install calceph                              # macOS
# apt install libcalceph-dev                     # Debian/Ubuntu

# 下星历数据
cd test/eph_hp && ./fetch_data.sh de440
cd ../..

# reconfigure 让 CMake 检测到 calceph
cmake -S . -B build
cmake --build build --target de_hp_test
./build/bin/de_hp_test
```

`de_hp_test` 三重优雅退化:
- calceph 未装 → CMake 不定义 target,打印安装提示
- calceph 已装,但 HP 数据缺失 → 跳过 target
- 全部就绪但 `de440s.bsp` 缺失 → 运行时打印 SKIP 提示,`exit 0`(CI 友好)

### 4.3 VSOP2013 Python 评估

```bash
cd test/eph_hp && ./fetch_data.sh vsop2013

# 验证数据完整性
python3 tools/vsop2013_eval.py data/VSOP2013p3.dat --verify

# 在指定年份评估
python3 tools/vsop2013_eval.py data/VSOP2013p3.dat --at 2000 2020 3000 5000

# 生成 CSV 用作 C++ 未来实现的对照参考
python3 tools/vsop2013_eval.py data/VSOP2013p3.dat --csv out.csv --range -1000 5000 500
```

---

## 五、常见问题排查

### 问:cmake 打印 "HP: eph_hp_test 已跳过 (缺少数据/header)"?
**答**:还没跑过 `fetch_data.sh`。执行 `cd test/eph_hp && ./fetch_data.sh small` 后再 `cmake -S . -B build` 重新配置。

### 问:`./fetch_data.sh` 报网络错误?
**答**:IMCCE / IERS / NASA 三个上游都在国外。可用代理:
```bash
export https_proxy=http://127.0.0.1:7890
./fetch_data.sh small
```
或按第九节手动 `curl` 下载后放到 `data/` 目录,再跑 `./fetch_data.sh gen`。

### 问:生成 header 报 "python3: command not found"?
**答**:装 Python 3(macOS 通常自带,Linux `apt install python3`)。无需第三方 pip 包。

### 问:每次编译是不是都要重下数据?
**答**:不需要。数据一旦下载并生成 header 就永久存在于 `data/` 和 `test/eph_hp/*.h`(都被 `.gitignore` 排除,不受 `git clean` 影响,除非用 `git clean -fdX`)。

### 问:CI 环境怎么用?
**答**:在 CI 脚本里加一行 `cd test/eph_hp && ./fetch_data.sh small` 即可(500 KB 下载,秒级)。DE440/VSOP2013 大文件按需。

---

## 六、当前实测结果摘要 (2026-07)

### VSOP87D 全项 vs 现有截断版
- Test A (1900-2024): 13/13 通过,最大差 0.025″ (0.6 秒)
- Test B (-4000~9999): 最大差 0.33″ (8 秒),出现在 Y=-4000
- Test C (S_aLon): Y=9999 差 32 秒

### IAU 2000A 章动 vs 现有中精度章动
- 近现代 (1900-2020): 3-25 mas 差异
- 远古/远期 (Y=9999): 最大 76 mas
- 符合"中精度 ~10 mas / 全项 ~0.2 mas"的量级预期

### DE440 vs VSOP87D-HP (Test F)
- 1900-2150 各年份 `|dL|` 与理论岁差(~50″/年)完全吻合(< 5″ 残差)
- **两大独立物理引擎的深度交叉验证:DE440 与 VSOP87D 精度都到 sub-arcsec**

---

## 七、将来"启用" HP 数据(切换到生产库)

**不要修改本目录里的文件**,而是:

### 启用 VSOP87D 全项(最小改动)
1. `./fetch_data.sh small` 生成 `vsop87d_earth_hp.h`
2. 把该文件里的 `XL0_0_HP[]` 数组内容复制到 `src/eph_data.h`(替换 `XL0_0`)
3. 保留 `XL0[]` 指针数组不变
4. 重新构建、跑 `pingqi_dingqi_test` 全时段回归

### 启用 IAU 2000A 章动
1. `./fetch_data.sh small` 生成 `iau2000a_nut_lon_hp.h`
2. 拷贝该文件到 `src/`
3. 拷贝 `eph_hp.cpp` 里 `nutationLon_HP()` 和 14 个 fundamental arguments 到 `src/eph.cpp`
4. 让 `nutation2()` / `nutationLon2()` 内部改调用 HP 版本
5. 跑回归

### 启用 VSOP2013(重大工程)
- 需要在 C++ 里实现:VSOP2013 series 求和 + Kepler equation 求解 + BCRS/ICRS→date 坐标变换
- 建议先用 Python 参考评估器输出大量 (JD, L, B, R) 参考值 CSV,然后写 C++ 实现对照

### 启用 DE441(需第三方运行时)
- libsxwnl 需要在打包时链接 calceph 或类似库
- 涉及所有平台(Android/iOS/HarmonyOS)的依赖管理,不推荐
- 更适合在**服务端**跑,作为 gold standard 提供 API

---

## 八、与主库解耦保证

- `src/CMakeLists.txt` 使用 `file(GLOB SRC_FILES "*.cpp")`,只扫描 `src/*.cpp`,不会扫到 `test/eph_hp/*.cpp`
- `test/CMakeLists.txt` 只在 test 构建时才 `add_subdirectory(eph_hp)`
- calceph 依赖只出现在 `test/eph_hp/CMakeLists.txt` 的可选分支,主库看不到
- 数据缺失时 CMake 直接跳过 HP target,主库 build 完全无影响
- Android/iOS/HarmonyOS 端 `libsxwnl` 二进制体积、依赖**完全不变**

---

## 九、上游数据源(手动 curl 备份)

若 `fetch_data.sh` 因网络问题失败,可手动下载:

```bash
# 小文件 (~500 KB)
curl -o data/VSOP87D.ear             https://ftp.imcce.fr/pub/ephem/planets/vsop87/VSOP87D.ear
curl -o data/IAU2000A_tab5.3a.txt    https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3a.txt
curl -o data/VSOP2013.f              https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013.f
curl -o data/VSOP2013.ctl            https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013.ctl

# 大文件
curl -o data/VSOP2013p3.dat          https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013p3.dat
curl -o data/de440s.bsp              https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
# DE441 完整版 (3.2 GB), 只在需要长时段时:
# curl -o data/de441.bsp             https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441.bsp

# 数据到位后生成 header:
./fetch_data.sh gen
```

---

## 十、精度上限与实用建议

Y=9999 时各源的精度对比:

| 源 | 视黄经精度 (arcsec) | 时间等效精度 (秒) |
|---|---|---|
| 现有 XL0_0 截断 | ~1.3 (基准) | ~32 |
| VSOP87D 全项 | ~0.05 | ~1.2 |
| VSOP2013 | ~0.001 | ~0.024 |
| JPL DE441 | ~0.0001 | ~0.0024 |
| **ΔT 物理下限** | (不适用) | **~3600 (约 1 小时)** |

**结论**:对**八字/命理**用途,VSOP87D 全项已远远超过 ΔT 的物理下限。VSOP2013 / DE441 只在**天文观测精算**(如日食贝塞尔要素、月食接触时刻)才值得启用。这批"备件"已足够满足未来 10-50 年内所有可预见的精度需求。
