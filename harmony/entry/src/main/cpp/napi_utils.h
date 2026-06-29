#pragma once
#include <napi/native_api.h>
#include <string>
#include <vector>
#include <cstring>

// ─── napi value creation ───

inline napi_value js_str(napi_env env, const char *val) {
    napi_value r;
    if (!val) val = "";
    napi_create_string_utf8(env, val, strlen(val), &r);
    return r;
}

inline napi_value js_str(napi_env env, const std::string &val) {
    napi_value r;
    napi_create_string_utf8(env, val.c_str(), val.length(), &r);
    return r;
}

inline napi_value js_int(napi_env env, int32_t val) {
    napi_value r;
    napi_create_int32(env, val, &r);
    return r;
}

inline napi_value js_dbl(napi_env env, double val) {
    napi_value r;
    napi_create_double(env, val, &r);
    return r;
}

inline napi_value js_bool(napi_env env, bool val) {
    napi_value r;
    napi_get_boolean(env, val, &r);
    return r;
}

inline napi_value js_null(napi_env env) {
    napi_value r;
    napi_get_null(env, &r);
    return r;
}

// ─── fluent object builder ───

class NObj {
public:
    NObj(napi_env e) : e_(e) { napi_create_object(e_, &o_); }

    NObj &s(const char *k, const char *v)    { napi_set_named_property(e_, o_, k, js_str(e_, v)); return *this; }
    NObj &s(const char *k, const std::string &v) { return s(k, v.c_str()); }
    NObj &i(const char *k, int32_t v)        { napi_set_named_property(e_, o_, k, js_int(e_, v)); return *this; }
    NObj &d(const char *k, double v)         { napi_set_named_property(e_, o_, k, js_dbl(e_, v)); return *this; }
    NObj &b(const char *k, bool v)           { napi_set_named_property(e_, o_, k, js_bool(e_, v)); return *this; }
    NObj &v(const char *k, napi_value val)   { napi_set_named_property(e_, o_, k, val); return *this; }

    operator napi_value() const { return o_; }

private:
    napi_env e_;
    napi_value o_;
};

// ─── array builder ───

class NArr {
public:
    NArr(napi_env e, uint32_t sz = 0) : e_(e), n_(0) { napi_create_array_with_length(e_, sz, &a_); }

    NArr &push(napi_value v) { napi_set_element(e_, a_, n_++, v); return *this; }
    NArr &push(const char *v) { return push(js_str(e_, v)); }
    NArr &push(int32_t v)    { return push(js_int(e_, v)); }

    operator napi_value() const { return a_; }

private:
    napi_env e_;
    napi_value a_;
    uint32_t n_;
};

// ─── argument reader ───

class NArgs {
public:
    NArgs(napi_env e, napi_callback_info info, size_t n) : e_(e), c_(n), a_(n) {
        napi_get_cb_info(e, info, &c_, a_.data(), nullptr, nullptr);
    }

    int32_t     intAt(int i)  const { int32_t v=0;  if(i<(int)c_) napi_get_value_int32(e_, a_[i], &v); return v; }
    double      dblAt(int i)  const { double v=0;   if(i<(int)c_) napi_get_value_double(e_, a_[i], &v); return v; }
    bool        boolAt(int i) const { bool v=false;  if(i<(int)c_) napi_get_value_bool(e_, a_[i], &v); return v; }
    std::string strAt(int i)  const {
        if(i>=(int)c_) return "";
        size_t len=0; napi_get_value_string_utf8(e_, a_[i], nullptr, 0, &len);
        std::string s(len, '\0'); napi_get_value_string_utf8(e_, a_[i], s.data(), len+1, &len);
        return s;
    }

    // read from object argument
    napi_value prop(int i, const char *k) const {
        if(i>=(int)c_) return nullptr;
        napi_value v; napi_get_named_property(e_, a_[i], k, &v); return v;
    }
    int32_t     objInt(int i, const char *k)  const { int32_t v=0;  auto p=prop(i,k); if(p) napi_get_value_int32(e_, p, &v); return v; }
    double      objDbl(int i, const char *k)  const { double v=0;   auto p=prop(i,k); if(p) napi_get_value_double(e_, p, &v); return v; }
    bool        objBool(int i, const char *k) const { bool v=false;  auto p=prop(i,k); if(p) napi_get_value_bool(e_, p, &v); return v; }
    std::string objStr(int i, const char *k)  const {
        auto p = prop(i,k); if(!p) return "";
        size_t len=0; napi_get_value_string_utf8(e_, p, nullptr, 0, &len);
        std::string s(len, '\0'); napi_get_value_string_utf8(e_, p, s.data(), len+1, &len);
        return s;
    }

private:
    napi_env e_;
    size_t c_;
    std::vector<napi_value> a_;
};

// ─── export helper ───

struct NExport { const char *name; napi_callback fn; };

inline void napi_export_all(napi_env env, napi_value exports,
                            const NExport *e, size_t n) {
    for (size_t i = 0; i < n; i++) {
        napi_value fn;
        napi_create_function(env, e[i].name, NAPI_AUTO_LENGTH, e[i].fn, nullptr, &fn);
        napi_set_named_property(env, exports, e[i].name, fn);
    }
}
