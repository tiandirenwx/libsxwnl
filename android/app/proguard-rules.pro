# ═══════════════════════════════════════════════════════════
# 寿星万年历 - Android Release ProGuard 规则
# ═══════════════════════════════════════════════════════════

# Compose 默认规则由 AGP 自动注入, 这里只补充项目特有保留项

# 保留 JNI 回调的 data class (jni_bridge.cpp 使用反射创建并赋值字段)
-keep class com.sxwnl.calendar.data.** { *; }

# 保留 native 方法 (与 libsxwnl_capi 链接)
-keepclasseswithmembernames class * {
    native <methods>;
}

# 保留 Kotlin Coroutines 内部类
-keepnames class kotlinx.coroutines.internal.MainDispatcherFactory {}
-keepnames class kotlinx.coroutines.CoroutineExceptionHandler {}
