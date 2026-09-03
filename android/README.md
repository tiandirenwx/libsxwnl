# 寿星万年历 — Android 版

基于 `libsxwnl` C++20 核心库的 Android 原生万年历应用，使用 **Kotlin + Jetpack Compose Material3**。

## 功能（与鸿蒙端一致）

- **月历 (CalendarScreen)** — 公农历对照、节气、月相 (朔/望/上弦/下弦)、节假日、底部 RTS / 月事件、点击日详情 Sheet
- **年历 (YearCalendarScreen)** — 年度概览、闰月、月度节气列表、朔日干支、可前后翻年
- **八字 (BaziScreen)**
  - 输入：姓名、性别、出生地经纬度、定气/平气夏至/平气冬至历法、真太阳时
  - 时间选择 (BaziDateTimePicker)：公历 / 农历 / 四柱反推 三 tab；滚轮选择月日时分；反推支持时辰未知
  - 结果 (BaziResultView)：简洁版「宣纸命书」 + 专业版命盘表（主星/天干/地支/藏干/纳音/星运/自坐/空亡/神煞 + 司令 + 五行旺衰 + 大运/流年）
- **日月食 (EclipseScreen)** — 日食/月食搜索、日月出没 (含经纬度输入)

## 架构

```
Kotlin/Compose (Material3) UI
    ↓
SxwnlNative.kt (Kotlin) ↔ JNI (jni_bridge.cpp + jni_eclipse_map.cpp)
    ↓
公共 C API (capi/sxwnl_capi.{h,cpp})  ← 三端共享
    ↓
libsxwnl C++20 核心 (src/*.cpp)
```

## 构建

### 工具链要求

| 工具 | 版本 |
|---|---|
| JDK | 17 (Temurin 推荐) |
| Gradle | 8.5 (项目 wrapper 已锁定) |
| AGP | 8.2.2 |
| Kotlin | 1.9.22 |
| Compose BOM | 2024.02 |
| Android SDK | 34 (compileSdk / targetSdk) |
| NDK | 25.1.8937393 |
| CMake | 3.22.1 |
| minSdk | 26 (Android 8.0) |

### 命令行构建

```bash
cd android

# Debug APK
./gradlew :app:assembleDebug

# Release APK (需要签名配置)
./gradlew :app:assembleRelease

# 安装到已连接设备
./gradlew :app:installDebug
```

构建产物：`android/app/build/outputs/apk/debug/app-debug.apk`

### Android Studio

打开 `android/` 目录即可。首次同步会自动下载所需 SDK Platforms 和 NDK。

## 工程目录

```
android/
├── app/
│   ├── build.gradle.kts
│   └── src/main/
│       ├── AndroidManifest.xml
│       ├── cpp/                    # JNI 桥接 + CMakeLists
│       │   ├── CMakeLists.txt
│       │   ├── jni_bridge.cpp
│       │   └── jni_eclipse_map.cpp
│       ├── java/com/sxwnl/calendar/
│       │   ├── MainActivity.kt     # 4 tab 主导航
│       │   ├── bridge/             # SxwnlNative.kt (JNI 声明)
│       │   ├── data/               # Models + CalendarRepository
│       │   ├── util/               # BaziCalc, YearUtil
│       │   └── ui/
│       │       ├── components/     # WheelPicker / BaziDateTimePicker / BaziResultView
│       │       ├── screens/        # CalendarScreen / YearCalendarScreen / BaziScreen / EclipseScreen
│       │       └── theme/          # Color, Dimens, Theme
│       └── res/
│           ├── drawable/           # adaptive icon foreground (☯)
│           ├── mipmap-anydpi-v26/  # adaptive launcher
│           └── values/
├── build.gradle.kts                # plugins 仅声明
├── settings.gradle.kts             # 仓库 + 模块包含
├── gradle.properties
└── gradle/wrapper/                 # Gradle 8.5 wrapper
```

## C++ 核心如何被链接

`app/src/main/cpp/CMakeLists.txt` 通过相对路径直接引用仓库根的 `src/*.cpp` 与 `capi/sxwnl_capi.cpp`，生成 `libsxwnl_jni.so`。无需先编译 libsxwnl 静态库。
