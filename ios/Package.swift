// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "SxwnlCalendar",
    platforms: [
        .iOS(.v17),
        .macOS(.v14)
    ],
    products: [
        .library(name: "SxwnlCore", targets: ["SxwnlCore"]),
    ],
    targets: [
        .target(
            name: "SxwnlCore",
            path: "capi",  // ← 现在是相对于 package 根目录的路径
            sources: ["sxwnl_capi.cpp"],
            publicHeadersPath: ".",
            cxxSettings: [
                .headerSearchPath("../src"),
                .define("SXWNL_BUILDING_CAPI"),
                .unsafeFlags(["-std=c++17"])
            ]
        ),
    ]
)
