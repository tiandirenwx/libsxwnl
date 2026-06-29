# 寿星万年历 - 跨平台架构设计

## 整体分层架构

```
┌──────────────────────────────────────────────────────────────────┐
│                      各平台原生 UI 层                             │
│                                                                  │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐       │
│  │  ArkTS/ArkUI │  │ Kotlin/      │  │ SwiftUI          │       │
│  │  (HarmonyOS) │  │ Compose(Droid│  │ (iOS/macOS)      │       │
│  └──────┬───────┘  └──────┬───────┘  └────────┬─────────┘       │
│         │                 │                    │                  │
├─────────┼─────────────────┼────────────────────┼─────────────────┤
│         │    各平台薄桥接层（≈100 行/平台）       │                  │
│         │                 │                    │                  │
│  ┌──────┴───────┐  ┌──────┴───────┐  ┌────────┴─────────┐       │
│  │  N-API       │  │  JNI         │  │  Swift-C Bridge  │       │
│  │  (napi_*.cpp)│  │  (jni_*.cpp) │  │  (SxwnlBridge)   │       │
│  └──────┬───────┘  └──────┬───────┘  └────────┬─────────┘       │
│         │                 │                    │                  │
├─────────┴─────────────────┴────────────────────┴─────────────────┤
│                    公共 C API 层（唯一封装点）                       │
│                                                                  │
│              sxwnl_capi.h / sxwnl_capi.cpp                      │
│              ┌──────────────────────────────┐                    │
│              │ • 日历查询  sxwnl_get_day_*  │                    │
│              │ • 节气列表  sxwnl_get_jieqi_*│                    │
│              │ • 八字排盘  sxwnl_bazi_*     │                    │
│              │ • 日月食    sxwnl_*eclipse*  │                    │
│              │ • 升降计算  sxwnl_*rise*     │                    │
│              └──────────────────────────────┘                    │
│                                                                  │
├──────────────────────────────────────────────────────────────────┤
│                    libsxwnl C++17 核心                            │
│                                                                  │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌──────────────┐       │
│  │  Day    │  │BaziBase │  │  SSQ    │  │ Eph/XL 星历  │        │
│  │  sxtwl  │  │         │  │ (气朔)  │  │ EphRsgs(食)  │        │
│  │  JD     │  │         │  │         │  │ SunMoonRise  │        │
│  └─────────┘  └─────────┘  └─────────┘  └──────────────┘       │
│                                                                  │
│            零外部依赖 · 纯标准库 · 星历数据内嵌                      │
└──────────────────────────────────────────────────────────────────┘
```

## 内存管理策略

### 值类型（无需手动释放）
- `SxwnlDayInfo` — 栈上结构体，fill 后即可使用
- `SxwnlJieQiItem` — 同上
- `SxwnlPillar` / `SxwnlDaYun` — 同上

### Handle 类型（需 create/free 配对）
- `SxwnlBazi` — `sxwnl_bazi_create()` 分配，`sxwnl_bazi_free()` 释放
  - getter 返回的 `const char*` 生命周期绑定到 handle，handle 释放前有效

### 动态字符串（需 sxwnl_string_free）
- `sxwnl_search_eclipses()` 等返回 `char*`
- 调用方使用后必须调用 `sxwnl_string_free()` 释放

### 各平台桥接层的内存安全
- **N-API (鸿蒙)**：在函数内创建 C 对象，转为 napi_value 后立即释放
- **JNI (Android)**：在函数内创建 C 对象，转为 jobject 后立即释放
- **Swift (Apple)**：使用 defer { free() } 确保释放，SxwnlBridge 单例封装

## 目录结构

```
libsxwnl/
├── src/                    # C++17 核心库（原有）
├── capi/                   # 公共 C API 封装层
│   ├── sxwnl_capi.h       # C 头文件（三端共用）
│   ├── sxwnl_capi.cpp     # 实现
│   ├── CMakeLists.txt      # 独立编译配置
│   └── test_capi.cpp       # 编译验证
├── harmony/                # 鸿蒙应用
│   └── entry/src/main/
│       ├── cpp/            # N-API 桥接
│       └── ets/            # ArkTS UI
├── android/                # Android 应用
│   └── app/src/main/
│       ├── cpp/            # JNI 桥接
│       └── java/           # Kotlin UI
├── ios/                    # iOS/macOS 应用
│   └── SxwnlCalendar/
│       ├── Bridge/         # Swift-C 桥接
│       └── Views/          # SwiftUI
└── ARCHITECTURE.md         # 本文件
```

## 扩展指南

### 新增功能（如择日、风水）

1. **C++ 核心**：在 `src/` 中实现算法
2. **C API**：在 `sxwnl_capi.h` 中声明接口，`sxwnl_capi.cpp` 中实现
3. **各平台桥接**：各约 10-20 行转换代码
4. **UI**：各平台独立实现界面

### 新增平台

只需实现该平台的桥接层（调用 capi）和 UI，核心代码零改动。

## 三端 App 实施状态

| 模块 | C API | 鸿蒙 | Android | iOS |
|---|:-:|:-:|:-:|:-:|
| 月历（DayInfo + JieQi + MoonPhase + RTS） | ✅ | ✅ | ✅ | ✅ |
| 年历（Year overview / 月节气 / 朔日干支） | ✅ | ✅ | ✅ | ✅ |
| 八字排盘（公历 / 农历 / 反推） | ✅ | ✅ | ✅ | ✅ |
| 八字结果（简洁版 + 专业版） | ✅ | ✅ | ✅ | ✅ |
| 八字保存截图到相册 | — | ✅ | — (可选扩展) | ✅ |
| 日月食搜索 + 日月出没 | ✅ | ✅ | ✅ | ✅ |
| 食地图 / 行星事件 / 星表 | ✅ | — | — | — |

> 表中后三项 C API 已就绪，但鸿蒙端 UI 未实现，故 Android / iOS 暂未对齐——纯下游接入工作，未来按需求增量补 UI 即可。

## 三端构建一览

| 平台 | 工具链 | 主要构建命令 |
|---|---|---|
| C++ 核心 | CMake 3.22+, Clang/GCC C++17 | `cmake -B build && cmake --build build` |
| Android | JDK 17, Gradle 8.5, AGP 8.2.2, NDK 25.1.8937393, CMake 3.22.1 | `cd android && ./gradlew :app:assembleDebug` |
| iOS | Xcode 15+, iOS 16+ | 打开 Xcode 工程，⌘R |
| 鸿蒙 | DevEco Studio | 通过 IDE 构建 |
