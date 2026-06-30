# 寿星万年历 — iOS / iPadOS 版

基于 `libsxwnl` C++17 核心库的 Apple 平台万年历应用，使用 **SwiftUI**。

## 功能（与鸿蒙端一致）

- **月历 (CalendarView)** — 公农历对照、节气、月相 (朔/望/上弦/下弦)、节假日、底部 RTS / 月事件、点击日详情 Sheet
- **年历 (YearCalendarView)** — 年度概览、闰月、月度节气列表、朔日干支、可前后翻年
- **八字 (BaziView)**
  - 输入：姓名、性别、出生地经纬度、定气/平气夏至/平气冬至历法、真太阳时
  - 时间选择 (BaziDateTimePicker)：公历 / 农历 / 四柱反推 三 tab；wheel picker 选择月日时分；反推支持时辰未知
  - 结果 (BaziResultView)：简洁版「宣纸命书」 + 专业版命盘表 + 大运/流年；支持 `ImageRenderer` 保存为图片到相册
- **日月食 (EclipseView)** — 日食/月食搜索、日月出没 (含经纬度输入)

## 架构

```
SwiftUI
    ↓
SxwnlBridge.swift  ←  Bridging Header (sxwnl_capi.h)
    ↓
公共 C API (capi/sxwnl_capi.{h,cpp})  ← 三端共享
    ↓
libsxwnl C++17 核心 (src/*.cpp)
```

## 构建

### 工具链要求

| 工具 | 版本 |
|---|---|
| Xcode | 15+ |
| iOS Deployment Target | 16.0+ |
| Swift | 5.9+ |
| C++ Dialect | C++17 |

### Xcode 工程

`ios/SxwnlCalendar/SxwnlCalendar.xcodeproj` 已配置好。克隆后直接打开，⌘R 运行即可。

编译前会自动执行 **Link C++ Sources** 构建阶段，将 `src/`、`capi/` 符号链接到 `SxwnlCalendar/CppSources/`（无需手动运行 `setup_xcode.sh`）。命令行打包时 `scripts/build_ios.sh` 也会同步链接。

### 工程关键配置

| 设置 | 值 |
|---|---|
| Objective-C Bridging Header | `SxwnlCalendar/Bridge/SxwnlCalendar-Bridging-Header.h` |
| Header Search Paths | `$(PROJECT_DIR)/../../capi` `$(PROJECT_DIR)/../../src` |
| C++ Language Dialect | GNU++17 |
| C++ and Objective-C Interop | C++ / Objective-C++ |
| Targeted Device Family | 1,2 (iPhone, iPad) |

### 命令行 typecheck（不实际编译产物，仅做静态校验）

```bash
SDK=$(xcrun --sdk iphoneos --show-sdk-path)
swiftc -typecheck -sdk "$SDK" -target arm64-apple-ios16.0 \
  -import-objc-header ios/SxwnlCalendar/SxwnlCalendar/Bridge/SxwnlCalendar-Bridging-Header.h \
  -I capi \
  $(find ios/SxwnlCalendar/SxwnlCalendar -name '*.swift')
```

## 工程目录

```
ios/SxwnlCalendar/
├── SxwnlCalendar.xcodeproj/
└── SxwnlCalendar/
    ├── SxwnlCalendarApp.swift      # 入口 + 4 tab 主导航
    ├── Bridge/
    │   ├── SxwnlCalendar-Bridging-Header.h
    │   └── SxwnlBridge.swift       # 全部 C API 的 Swift 封装
    ├── Models/                     # AppTheme, CalendarModels, BaziCalc, YearUtil
    └── Views/
        ├── CalendarView.swift
        ├── YearCalendarView.swift
        ├── BaziView.swift
        ├── BaziDateTimePicker.swift
        ├── BaziResultView.swift
        └── EclipseView.swift
```

## C++ 核心如何被链接

`CppSources/` 不入库（见根目录 `.gitignore`），由构建脚本在本地生成指向 `src/`、`capi/` 的符号链接。Xcode 每次编译前自动运行 `ios/setup_xcode.sh --link-only`；新增 `.cpp` 后重新编译即可被纳入。Bridging Header 暴露 C API 给 Swift。

## 共享 UI 资源 (八字页)

字体 `WenYue.otf`、12 张生肖图、`bz_paper.jpg` 仅在仓库根 `assets/bazi/` 保留一份。三端编译前由 `scripts/sync_bazi_assets.sh` 自动同步到各平台资源目录（Android/iOS 符号链接，鸿蒙复制实体文件；均由构建流程自动触发），换机器 clone 后直接编译即可，无需手动操作。

## 兼容性提示

- 由于 deployment target 是 iOS 16，`SwiftUI.ImageRenderer` 等 API 可正常使用。
- 真太阳时计算需要 SwiftUI 输入经纬度；默认值为北京 (116.4°E, 39.9°N)。
- 农历支持的天文年范围为 `-721 ~ 9999`（与核心库一致），公元前年份用 `B721` / `-720` 等输入格式。
