# 寿星万年历 - 鸿蒙版 (HarmonyOS NEXT)

基于 libsxwnl C++ 核心库的鸿蒙原生万年历应用。

## 功能

- **万年历**：公历/农历对照，节气标注，干支纪年，星座建除
- **八字排盘**：四柱排盘，大运流年，支持定气/平气多种历法
- **日月食工具**：日食月食搜索与详细计算

## 架构

```
ArkTS UI (ArkUI)
    ↓
N-API Bridge (napi_*.cpp)
    ↓
C API (sxwnl_capi.h) ← 三端共享
    ↓
libsxwnl (C++20 核心)
```

## 构建

需要 DevEco Studio 5.0+ 和 HarmonyOS SDK。

```bash
# 在 DevEco Studio 中打开 harmony/ 目录
# 或使用命令行
hvigorw assembleHap
```

## 目录结构

```
harmony/
├── entry/src/main/
│   ├── ets/          # ArkTS UI 代码
│   ├── cpp/          # N-API 桥接层
│   └── resources/    # 资源文件
├── build-profile.json5
└── oh-package.json5
```
