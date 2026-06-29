#include "month_calendar.h"

#include <memory>
#include <stdexcept>
#include <string>

#include "JD.h"
#include "const.h"
#include "lunar_ob.h"
#include "sx_lang_zh.h"
#include "sxtwl.h"

namespace sxwnl
{

// ═══════════════════════════════════════════════════════════════
//  MonthCalendar
// ═══════════════════════════════════════════════════════════════
MonthCalendar::MonthCalendar(int solarYear, int solarMonth)
    : year_(solarYear), month_(solarMonth),
      w0_(0), totalWeeks_(0), d0_(0)
{
    if (solarMonth < 1 || solarMonth > 12)
        throw std::invalid_argument("MonthCalendar: month must be in [1,12]");

    // ── 月首/月末 ───────────────────────────────────────
    Time t{};
    t.h = 12; t.m = 0; t.s = 0.1;
    t.Y = solarYear; t.M = solarMonth; t.D = 1;
    d0_ = int2(JD::toJD(t)) - J2000;

    int ny = solarYear, nm = solarMonth + 1;
    if (nm > 12) { ++ny; nm = 1; }
    Time t2{}; t2.h = 12; t2.m = 0; t2.s = 0.1;
    t2.Y = ny; t2.M = nm; t2.D = 1;
    int dn = int2(JD::toJD(t2)) - J2000 - d0_;

    w0_ = (uint8_t)((d0_ + J2000 + 1 + 7000000) % 7);
    totalWeeks_ = (uint8_t)((w0_ + dn - 1) / 7 + 1);

    days_.reserve(dn);

    // ── 逐日填充 (复用 Day 已验证逻辑) ─────────────────
    for (int i = 0; i < dn; ++i)
    {
        std::shared_ptr<Day> d(sxtwl::fromSolar(solarYear, (uint8_t)solarMonth, i + 1));
        if (!d) throw std::runtime_error("MonthCalendar: failed to compute day");
        days_.push_back(CalendarDay{std::move(d), i});
    }

    // ── 年级别信息(春节为界), 全部由首日 Day 派生 ──────
    if (!days_.empty())
    {
        Day &first = *days_.front().day;
        yearGZ_    = first.getYearGZ(true);
        shengXiao_ = first.getShengXiao();
    }
    {
        LunarYear ly(solarYear);
        nianHao_ = ly.getNianHao();
    }
}

} // namespace sxwnl
