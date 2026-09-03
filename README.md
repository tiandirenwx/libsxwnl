# libsxwnl
## 寿星万年历 C++ 核心 + 三端原生 App（HarmonyOS / Android / iOS）

### 说明

1. 本程序参考以下开源项目，仅限于兴趣研究，无任何商业化目的。
   - [sxwnl](https://github.com/sxwnl/sxwnl.git)
   - [sxwnl-cpp](https://github.com/sxwnl/sxwnl-cpp.git)
   - [sxtwl_cpp](https://github.com/yuangu/sxtwl_cpp.git)

   在此向许剑伟老师及上述作者致敬！

2. 程序修正了上述代码中的一些小问题，新增以下能力：
   - 支持平气法排盘（计算冬至仍不完善，精度控制在百年内）
   - 支持公元前 -721 到 9999 的农历转换（春节定界）
   - 修正历史上各种特殊历法，做了适配
   - **新增统一 C API 层 + 三端原生 UI**（HarmonyOS / Android / iOS）

3. 程序还有未测试完善的地方，欢迎大家批评指正！

---

## 三端 App 功能矩阵

四个 tab，三端 100% 功能对齐：

| Tab | 鸿蒙 (ArkUI) | Android (Compose) | iOS (SwiftUI) |
|---|:-:|:-:|:-:|
| 月历 | ✅ | ✅ | ✅ |
| 年历 | ✅ | ✅ | ✅ |
| 八字（含反推 + 简洁/专业版） | ✅ | ✅ | ✅ |
| 日月食 / 日月出没 | ✅ | ✅ | ✅ |

## 仓库结构

```
libsxwnl/
├── src/                    # C++20 核心库
├── capi/                   # 统一 C API（三端共用）
│   ├── sxwnl_capi.h        # 头文件
│   └── sxwnl_capi.cpp      # 实现
├── assets/bazi/            # 八字页共享资源 (字体/生肖图, 三端编译时自动链接)
├── harmony/                # 鸿蒙 App (ArkTS / ArkUI)
├── android/                # Android App (Kotlin / Compose) → 见 android/README.md
├── ios/                    # iOS App (SwiftUI)            → 见 ios/README.md
├── test/                   # C++ 单元测试
├── ARCHITECTURE.md         # 整体架构详解
└── README.md
```

## 快速构建

| 平台 | 命令 |
|---|---|
| C++ 核心 | `cmake -B build && cmake --build build` |
| Android | `cd android && ./gradlew :app:assembleDebug` |
| iOS | 打开 `ios/SxwnlCalendar/SxwnlCalendar.xcodeproj`，⌘R |
| 鸿蒙 | DevEco Studio 打开 `harmony/` |

各平台详细要求请参阅子目录 README。

## 架构总览

详见 [`ARCHITECTURE.md`](ARCHITECTURE.md)。简要：

```
平台 UI (ArkTS / Compose / SwiftUI)
    ↓
平台桥接 (N-API / JNI / Swift-C)
    ↓
统一 C API  (capi/sxwnl_capi.{h,cpp})
    ↓
libsxwnl C++20 核心（零外部依赖，星历数据内嵌）
```

核心层完全平台无关，新增平台只需补一层薄桥接 + UI。
