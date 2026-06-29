import java.util.Properties

plugins {
    id("com.android.application")
    id("org.jetbrains.kotlin.android")
}

// ── 读取签名配置 (优先级: keystore.properties > 环境变量 > debug.keystore) ──
val keystoreProps = Properties().apply {
    val f = rootProject.file("keystore.properties")
    if (f.exists()) f.inputStream().use { load(it) }
}

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
