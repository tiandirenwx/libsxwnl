#pragma once
#include <napi/native_api.h>

napi_value NapiGetDayInfo(napi_env env, napi_callback_info info);
napi_value NapiGetMonthData(napi_env env, napi_callback_info info);
napi_value NapiLunarToSolar(napi_env env, napi_callback_info info);
napi_value NapiSolarToLunar(napi_env env, napi_callback_info info);
napi_value NapiGetJieQiList(napi_env env, napi_callback_info info);
napi_value NapiGetYearLeapMonth(napi_env env, napi_callback_info info);
napi_value NapiGetLunarMonths(napi_env env, napi_callback_info info);
napi_value NapiGetLunarMonthDays(napi_env env, napi_callback_info info);
napi_value NapiGetSolarMonthValidDays(napi_env env, napi_callback_info info);
napi_value NapiGetYearCalendar(napi_env env, napi_callback_info info);
napi_value NapiCalcDayRTS(napi_env env, napi_callback_info info);
