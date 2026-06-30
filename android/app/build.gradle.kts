import java.util.Properties

plugins {
    id("com.android.application")
    id("org.jetbrains.kotlin.android")
}

// 编译前链接 assets/bazi/ 共享资源 (字体/生肖图/宣纸背景)
val repoRoot = rootProject.projectDir.parentFile!!
tasks.register<Exec>("syncBaziAssets") {
    workingDir(repoRoot)
    commandLine("bash", "scripts/sync_bazi_assets.sh", "android")
}
tasks.named("preBuild") { dependsOn("syncBaziAssets") }

// ── 读取签名配置 (优先级: keystore.properties > 环境变量 > debug.keystore) ──
val keystoreProps = Properties().apply {
    val f = rootProject.file("keystore.properties")
    if (f.exists()) f.inputStream().use { load(it) }
}

// 是否找到了真正的正式签名 (任意一种来源即可)
val hasReleaseKeystore: Boolean =
    keystoreProps.getProperty("storeFile") != null ||
        System.getenv("SXWNL_KEYSTORE_PATH") != null

fun resolveStoreFile(): java.io.File {
    keystoreProps.getProperty("storeFile")?.let { return rootProject.file(it) }
    System.getenv("SXWNL_KEYSTORE_PATH")?.let { return file(it) }
    return file("${System.getProperty("user.home")}/.android/debug.keystore")
}
fun resolveStorePassword(): String =
    keystoreProps.getProperty("storePassword")
        ?: System.getenv("SXWNL_KEYSTORE_PASSWORD")
        ?: "android"
fun resolveKeyAlias(): String =
    keystoreProps.getProperty("keyAlias")
        ?: System.getenv("SXWNL_KEY_ALIAS")
        ?: "androiddebugkey"
fun resolveKeyPassword(): String =
    keystoreProps.getProperty("keyPassword")
        ?: System.getenv("SXWNL_KEY_PASSWORD")
        ?: "android"

android {
    namespace = "com.sxwnl.calendar"
    compileSdk = 34

    defaultConfig {
        applicationId = "com.sxwnl.calendar"
        minSdk = 26
        targetSdk = 34
        versionCode = 1
        versionName = "1.0"

        ndk {
            abiFilters += listOf("arm64-v8a", "armeabi-v7a", "x86_64")
        }

        externalNativeBuild {
            cmake {
                cppFlags += "-std=c++17"
            }
        }
    }

    externalNativeBuild {
        cmake {
            path = file("src/main/cpp/CMakeLists.txt")
            version = "3.22.1"
        }
    }

    signingConfigs {
        create("release") {
            storeFile = resolveStoreFile()
            storePassword = resolveStorePassword()
            keyAlias = resolveKeyAlias()
            keyPassword = resolveKeyPassword()
        }
    }

    buildTypes {
        release {
            isMinifyEnabled = true
            proguardFiles(
                getDefaultProguardFile("proguard-android-optimize.txt"),
                "proguard-rules.pro"
            )
            signingConfig = signingConfigs.getByName("release")
        }
    }

    // ── 签名状态自检: 在配置阶段就把当前要用的 keystore 打印出来 ──
    // 同时, 如果要打 release 包却没有正式 keystore, 给出非常醒目的警告.
    gradle.taskGraph.whenReady {
        val ksFile = resolveStoreFile()
        val isReleaseBuild = allTasks.any {
            val n = it.name.lowercase()
            n.contains("release") && (n.startsWith("assemble") || n.startsWith("bundle"))
        }
        if (isReleaseBuild) {
            if (hasReleaseKeystore) {
                logger.lifecycle("[sxwnl] release 签名: ${ksFile.absolutePath} (alias=${resolveKeyAlias()})")
            } else {
                logger.warn("")
                logger.warn("╔═══════════════════════════════════════════════════════════════╗")
                logger.warn("║  [sxwnl] !! 警告: 正在用 debug.keystore 给 release 包签名 !!   ║")
                logger.warn("║                                                                ║")
                logger.warn("║  - 该包不可上架任何应用商店                                    ║")
                logger.warn("║  - 国产 ROM (小米/华为/OPPO/vivo) 可能直接拦截或卸载           ║")
                logger.warn("║  - 后续若换成正式签名, 用户必须卸载重装 (无法覆盖升级)         ║")
                logger.warn("║                                                                ║")
                logger.warn("║  生成正式 keystore: scripts/gen_android_keystore.sh            ║")
                logger.warn("╚═══════════════════════════════════════════════════════════════╝")
                logger.warn("")
            }
        }
    }

    compileOptions {
        sourceCompatibility = JavaVersion.VERSION_17
        targetCompatibility = JavaVersion.VERSION_17
    }

    kotlinOptions {
        jvmTarget = "17"
    }

    buildFeatures {
        compose = true
    }

    composeOptions {
        kotlinCompilerExtensionVersion = "1.5.8"
    }

    packaging {
        resources {
            excludes += "/META-INF/{AL2.0,LGPL2.1}"
        }
    }
}

dependencies {
    val composeBom = platform("androidx.compose:compose-bom:2024.02.00")
    implementation(composeBom)

    implementation("androidx.core:core-ktx:1.12.0")
    implementation("androidx.lifecycle:lifecycle-runtime-ktx:2.7.0")
    implementation("androidx.lifecycle:lifecycle-viewmodel-compose:2.7.0")
    implementation("androidx.activity:activity-compose:1.8.2")

    implementation("androidx.compose.ui:ui")
    implementation("androidx.compose.ui:ui-graphics")
    implementation("androidx.compose.ui:ui-tooling-preview")
    implementation("androidx.compose.material3:material3")
    implementation("androidx.compose.material:material-icons-extended")
    implementation("androidx.compose.foundation:foundation")

    implementation("androidx.navigation:navigation-compose:2.7.7")

    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-android:1.7.3")

    debugImplementation("androidx.compose.ui:ui-tooling")
    debugImplementation("androidx.compose.ui:ui-test-manifest")
}
