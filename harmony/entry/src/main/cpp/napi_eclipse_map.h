#pragma once
#include <napi/native_api.h>

napi_value NapiSearchSolarEclipses(napi_env env, napi_callback_info info);
napi_value NapiGetSolarEclipsePath(napi_env env, napi_callback_info info);
napi_value NapiGetLocalSolarEclipse(napi_env env, napi_callback_info info);
napi_value NapiSearchLunarEclipses(napi_env env, napi_callback_info info);
napi_value NapiGetLunarEclipseDetail(napi_env env, napi_callback_info info);
napi_value NapiGetWorldMapDitu0(napi_env env, napi_callback_info info);
napi_value NapiGetWorldMapData(napi_env env, napi_callback_info info);
