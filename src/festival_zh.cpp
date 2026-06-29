#include "festival_zh.h"
#include "day.h"
#include "sx_lang_zh.h"

namespace festival
{

// ═══════════════════════════════════════════════════════════════
//  公历定日节日表 (sFtv)
//  yearMin/yearMax 为 0 表示无限制(默认要求年份 >= 1850, 与上游一致)
// ═══════════════════════════════════════════════════════════════

struct SolarFixed
{
    uint8_t  month;
    uint8_t  day;
    char     type;
    int16_t  yearMin;  // 0 => 默认 >=1850
    int16_t  yearMax;  // 0 => 不设上限
    const char *name;
};

static constexpr SolarFixed SFTV[] = {
    // ─── 1 月 ───
    {1,  1, '#', 0,    0,    "元旦"},

    // ─── 2 月 ───
    {2,  2, 'I', 0,    0,    "世界湿地日"},
    {2, 10, '.', 0,    0,    "国际气象节"},
    {2, 14, 'I', 0,    0,    "情人节"},

    // ─── 3 月 ───
    {3,  1, '.', 0,    0,    "国际海豹日"},
    {3,  3, '.', 0,    0,    "全国爱耳日"},
    {3,  5, '.', 1963, 9999, "学雷锋纪念日"},
    {3,  8, 'I', 0,    0,    "妇女节"},
    {3, 12, 'I', 0,    0,    "植树节"},
    {3, 12, '.', 1925, 9999, "孙中山逝世纪念日"},
    {3, 14, '.', 0,    0,    "国际警察日"},
    {3, 15, 'I', 1983, 9999, "消费者权益日"},
    {3, 17, '.', 0,    0,    "中国国医节"},
    {3, 17, '.', 0,    0,    "国际航海日"},
    {3, 21, '.', 0,    0,    "世界森林日"},
    {3, 21, '.', 0,    0,    "消除种族歧视国际日"},
    {3, 21, '.', 0,    0,    "世界儿歌日"},
    {3, 22, 'I', 0,    0,    "世界水日"},
    {3, 23, 'I', 0,    0,    "世界气象日"},
    {3, 24, '.', 1982, 9999, "世界防治结核病日"},
    {3, 25, '.', 0,    0,    "全国中小学生安全教育日"},
    {3, 30, '.', 0,    0,    "巴勒斯坦国土日"},

    // ─── 4 月 ───
    {4,  1, 'I', 1564, 9999, "愚人节"},
    {4,  1, '.', 0,    0,    "全国爱国卫生运动月(四月)"},
    {4,  1, '.', 0,    0,    "税收宣传月(四月)"},
    {4,  7, 'I', 0,    0,    "世界卫生日"},
    {4, 22, 'I', 0,    0,    "世界地球日"},
    {4, 23, '.', 0,    0,    "世界图书和版权日"},
    {4, 24, '.', 0,    0,    "亚非新闻工作者日"},

    // ─── 5 月 ───
    {5,  1, '#', 1889, 9999, "劳动节"},
    {5,  4, 'I', 0,    0,    "青年节"},
    {5,  5, '.', 0,    0,    "碘缺乏病防治日"},
    {5,  8, '.', 0,    0,    "世界红十字日"},
    {5, 12, 'I', 0,    0,    "国际护士节"},
    {5, 15, 'I', 0,    0,    "国际家庭日"},
    {5, 17, '.', 0,    0,    "国际电信日"},
    {5, 18, '.', 0,    0,    "国际博物馆日"},
    {5, 20, '.', 0,    0,    "全国学生营养日"},
    {5, 23, '.', 0,    0,    "国际牛奶日"},
    {5, 31, 'I', 0,    0,    "世界无烟日"},

    // ─── 6 月 ───
    {6,  1, 'I', 1925, 9999, "国际儿童节"},
    {6,  5, '.', 0,    0,    "世界环境保护日"},
    {6,  6, '.', 0,    0,    "全国爱眼日"},
    {6, 17, '.', 0,    0,    "防治荒漠化和干旱日"},
    {6, 23, '.', 0,    0,    "国际奥林匹克日"},
    {6, 25, '.', 0,    0,    "全国土地日"},
    {6, 26, 'I', 0,    0,    "国际禁毒日"},

    // ─── 7 月 ───
    {7,  1, 'I', 1997, 9999, "香港回归纪念日"},
    {7,  1, 'I', 1921, 9999, "中共诞辰"},
    {7,  1, '.', 0,    0,    "世界建筑日"},
    {7,  2, '.', 0,    0,    "国际体育记者日"},
    {7,  7, 'I', 1937, 9999, "抗日战争纪念日"},
    {7, 11, 'I', 0,    0,    "世界人口日"},
    {7, 30, '.', 0,    0,    "非洲妇女日"},

    // ─── 8 月 ───
    {8,  1, 'I', 1927, 9999, "建军节"},
    {8,  8, '.', 0,    0,    "中国男子节(爸爸节)"},

    // ─── 9 月 ───
    {9,  3, 'I', 1945, 9999, "抗日战争胜利纪念"},
    {9,  8, '.', 1966, 9999, "国际扫盲日"},
    {9,  8, '.', 0,    0,    "国际新闻工作者日"},
    {9,  9, '.', 0,    0,    "毛泽东逝世纪念"},
    {9, 10, 'I', 0,    0,    "中国教师节"},
    {9, 14, '.', 0,    0,    "世界清洁地球日"},
    {9, 16, '.', 0,    0,    "国际臭氧层保护日"},
    {9, 18, 'I', 0,    0,    "九·一八事变纪念日"},
    {9, 20, '.', 0,    0,    "国际爱牙日"},
    {9, 27, '.', 0,    0,    "世界旅游日"},
    {9, 28, 'I', 0,    0,    "孔子诞辰"},

    // ─── 10 月 ───
    {10,  1, '#', 1949, 9999, "国庆节"},
    {10,  1, '.', 0,    0,    "世界音乐日"},
    {10,  1, '.', 0,    0,    "国际老人节"},
    {10,  2, '#', 1949, 9999, "国庆节假日"},
    {10,  2, '.', 0,    0,    "国际和平与民主自由斗争日"},
    {10,  3, '#', 1949, 9999, "国庆节假日"},
    {10,  4, '.', 0,    0,    "世界动物日"},
    {10,  6, '.', 0,    0,    "老人节"},
    {10,  8, '.', 0,    0,    "全国高血压日"},
    {10,  8, '.', 0,    0,    "世界视觉日"},
    {10,  9, '.', 0,    0,    "世界邮政日"},
    {10,  9, '.', 0,    0,    "万国邮联日"},
    {10, 10, 'I', 0,    0,    "辛亥革命纪念日"},
    {10, 10, '.', 0,    0,    "世界精神卫生日"},
    {10, 13, '.', 0,    0,    "世界保健日"},
    {10, 13, '.', 0,    0,    "国际教师节"},
    {10, 14, '.', 0,    0,    "世界标准日"},
    {10, 15, '.', 0,    0,    "国际盲人节(白手杖节)"},
    {10, 16, '.', 0,    0,    "世界粮食日"},
    {10, 17, '.', 0,    0,    "世界消除贫困日"},
    {10, 22, '.', 0,    0,    "世界传统医药日"},
    {10, 24, '.', 0,    0,    "联合国日"},
    {10, 31, '.', 0,    0,    "世界勤俭日"},

    // ─── 11 月 ───
    {11,  7, '.', 1917, 9999, "十月社会主义革命纪念日"},
    {11,  8, '.', 0,    0,    "中国记者日"},
    {11,  9, '.', 0,    0,    "全国消防安全宣传教育日"},
    {11, 10, '.', 0,    0,    "世界青年节"},
    {11, 11, '.', 0,    0,    "国际科学与和平周(本日所属的一周)"},
    {11, 12, '.', 0,    0,    "孙中山诞辰纪念日"},
    {11, 14, '.', 0,    0,    "世界糖尿病日"},
    {11, 17, '.', 0,    0,    "国际大学生节"},
    {11, 17, '.', 0,    0,    "世界学生节"},
    {11, 20, '.', 0,    0,    "彝族年"},
    {11, 21, '.', 0,    0,    "彝族年"},
    {11, 21, '.', 0,    0,    "世界问候日"},
    {11, 21, '.', 0,    0,    "世界电视日"},
    {11, 22, '.', 0,    0,    "彝族年"},
    {11, 29, '.', 0,    0,    "国际声援巴勒斯坦人民国际日"},

    // ─── 12 月 ───
    {12,  1, 'I', 1988, 9999, "世界艾滋病日"},
    {12,  3, '.', 0,    0,    "世界残疾人日"},
    {12,  5, '.', 0,    0,    "国际经济和社会发展志愿人员日"},
    {12,  8, '.', 0,    0,    "国际儿童电视日"},
    {12,  9, '.', 0,    0,    "世界足球日"},
    {12, 10, '.', 0,    0,    "世界人权日"},
    {12, 12, 'I', 0,    0,    "西安事变纪念日"},
    {12, 13, 'I', 0,    0,    "南京大屠杀(1937年)纪念日"},
    {12, 20, '.', 0,    0,    "澳门回归纪念"},
    {12, 21, '.', 0,    0,    "国际篮球日"},
    {12, 24, 'I', 0,    0,    "平安夜"},
    {12, 25, 'I', 0,    0,    "圣诞节"},
    {12, 26, '.', 0,    0,    "毛泽东诞辰纪念"},
};

// ═══════════════════════════════════════════════════════════════
//  公历周节日表 (wFtv)
//  weekN = 1-5; 5 表示"本月最后一个 X 曜日"
// ═══════════════════════════════════════════════════════════════
struct SolarWeekly
{
    uint8_t month;
    uint8_t weekN;     // 1-5 (5 = last)
    uint8_t weekday;   // 0-6 (0 = Sunday)
    char    type;
    const char *name;
};

static constexpr SolarWeekly WFTV[] = {
    { 1, 5, 0, 'I', "世界麻风日"},        // 1 月最后一个星期日
    { 5, 2, 0, '.', "国际母亲节"},        // 5 月第 2 个星期日
    { 5, 3, 0, 'I', "全国助残日"},        // 5 月第 3 个星期日
    { 6, 3, 0, '.', "父亲节"},            // 6 月第 3 个星期日
    { 7, 3, 0, '.', "被奴役国家周"},      // 7 月第 3 个星期日
    { 9, 3, 2, 'I', "国际和平日"},        // 9 月第 3 个星期二
    { 9, 4, 0, '.', "国际聋人节 世界儿童日"},
    { 9, 5, 0, 'I', "世界海事日"},        // 9 月最后一个星期日
    {10, 1, 1, '.', "国际住房日"},        // 10 月第 1 个星期一
    {10, 1, 3, 'I', "国际减轻自然灾害日(减灾日)"}, // 10 月第 1 个星期三
    {11, 4, 4, 'I', "感恩节"},            // 11 月第 4 个星期四
};

// ═══════════════════════════════════════════════════════════════
//  农历节日表
//   月份: 1=正月, 12=腊月
//   onlyIfNextIsZheng: 仅在"下一月为正月"时生效(除夕/小年)
//   名称按 '|' 分割为多个三元组: <type><name>
//    例: "#春节" 或 "#春节|.壮族歌墟节 苗族踩山节 达斡尔族卡钦"
// ═══════════════════════════════════════════════════════════════
struct LunarFixed
{
    uint8_t month;
    uint8_t day;
    bool    onlyIfNextIsZheng;
    bool    needDay30;  // 仅在 lunarMonthDays==30 时生效 (除夕之 30)
    bool    needDay29;  // 仅在 lunarMonthDays==29 时生效 (除夕之 29)
    const char *spec;   // "<type><name>[ | ...]"
};

static constexpr LunarFixed LFTV[] = {
    { 1,  1, false, false, false, "#春节"},
    { 1,  2, false, false, false, "#大年初二"},
    { 1, 15, false, false, false, "#元宵节|I上元节|.壮族歌墟节 苗族踩山节 达斡尔族卡钦"},
    { 1, 16, false, false, false, ".侗族芦笙节(至正月二十)"},
    { 1, 25, false, false, false, ".填仓节"},
    { 1, 29, false, false, false, ".送穷日"},

    { 2,  1, false, false, false, ".瑶族忌鸟节"},
    { 2,  2, false, false, false, "I春龙节(龙抬头)|.畲族会亲节"},
    { 2,  8, false, false, false, ".傈傈族刀杆节"},

    { 3,  3, false, false, false, "I北帝诞|.苗族黎族歌墟节"},
    { 3, 15, false, false, false, ".白族三月街(至三月二十)"},
    { 3, 23, false, false, false, "I天后诞 妈祖诞"},

    { 4,  8, false, false, false, "I牛王诞"},
    { 4, 18, false, false, false, ".锡伯族西迁节"},

    { 5,  5, false, false, false, "#端午节"},
    { 5, 13, false, false, false, "I关帝诞|.阿昌族泼水节"},
    { 5, 22, false, false, false, ".鄂温克族米阔鲁节"},
    { 5, 29, false, false, false, ".瑶族达努节"},

    { 6,  6, false, false, false, "I姑姑节 天贶节|.壮族祭田节 瑶族尝新节"},
    { 6, 24, false, false, false, ".火把节、星回节(彝、白、佤、阿昌、纳西、基诺族)"},

    { 7,  7, false, false, false, "I七夕(中国情人节,乞巧节,女儿节)"},
    { 7, 13, false, false, false, ".侗族吃新节"},
    { 7, 15, false, false, false, "I中元节 鬼节"},

    { 8, 15, false, false, false, "#中秋节"},

    { 9,  9, false, false, false, "I重阳节"},

    {10,  1, false, false, false, "I祭祖节(十月朝)"},
    {10, 15, false, false, false, "I下元节"},
    {10, 16, false, false, false, ".瑶族盘王节"},

    {12,  8, false, false, false, "I腊八节"},

    // 除夕 + 小年: 仅在"下一月为正月"时生效
    {12, 23, true,  false, false, "I小年"},
    {12, 30, true,  true,  false, "#除夕"},
    {12, 29, true,  false, true,  "#除夕"},
};

// ═══════════════════════════════════════════════════════════════
//  辅助: 拼接节日字符串(末尾追加空格分隔, 与上游显示一致)
// ═══════════════════════════════════════════════════════════════
static void appendByType(DayInfo &out, char type, const char *name)
{
    switch (type)
    {
    case TIER_HOLIDAY:
        out.holiday += name;
        out.holiday += ' ';
        out.isOffDay = true;
        break;
    case TIER_MAJOR:
        out.major += name;
        out.major += ' ';
        break;
    case TIER_MINOR:
        out.minor += name;
        out.minor += ' ';
        break;
    default:
        break;
    }
}

// 解析多类别串 "#春节|I上元节|.壮族..."
static void emitLunarSpec(DayInfo &out, const char *spec)
{
    const char *p = spec;
    while (*p)
    {
        char type = *p++;
        // 收集到下一个 '|' 或结尾
        const char *start = p;
        while (*p && *p != '|') ++p;
        // 复制 [start, p) 到临时
        std::string name(start, p - start);
        appendByType(out, type, name.c_str());
        if (*p == '|') ++p;
    }
}

// ═══════════════════════════════════════════════════════════════
//  公历定日节日
// ═══════════════════════════════════════════════════════════════
void appendSolarFixed(const DayContext &ctx, DayInfo &out)
{
    // 周末自动放假
    if (ctx.weekday == 0 || ctx.weekday == 6)
    {
        out.isOffDay = true;
    }

    for (const auto &f : SFTV)
    {
        if (f.month != ctx.month || f.day != ctx.day) continue;

        // 年份限制
        if (f.yearMin > 0 || f.yearMax > 0)
        {
            if (ctx.year < f.yearMin || ctx.year > f.yearMax) continue;
        }
        else
        {
            if (ctx.year < 1850) continue;
        }
        appendByType(out, f.type, f.name);
    }
}

// ═══════════════════════════════════════════════════════════════
//  公历周节日
//   规则: 本日是本月第 N 个 X 曜日;
//   "5" 表示"最后一个", 故位于本月最后一周的同曜日也同时匹配 5
// ═══════════════════════════════════════════════════════════════
void appendSolarWeekly(const DayContext &ctx, DayInfo &out)
{
    // 计算本日是本月第几个 (ctx.weekday) 曜日
    int w = ctx.weekiInMonth;
    if (ctx.weekday >= ctx.firstWeekday) ++w;

    int w2 = w;
    if (ctx.weekiInMonth == ctx.weekTotal - 1)
    {
        w2 = 5; // 最后一周, 允许匹配编号 5
    }

    for (const auto &f : WFTV)
    {
        if (f.month != ctx.month) continue;
        if (f.weekday != ctx.weekday) continue;
        if (f.weekN != w && f.weekN != w2) continue;
        appendByType(out, f.type, f.name);
    }
}

// ═══════════════════════════════════════════════════════════════
//  农历节日 + 节气
// ═══════════════════════════════════════════════════════════════
void appendLunarFixed(const DayContext &ctx, DayInfo &out)
{
    // 非闰月才匹配大部分农历节日; 闰月只匹配"除夕/小年"这类与下月名相关的
    for (const auto &f : LFTV)
    {
        if (f.month != ctx.lunarMonth || f.day != ctx.lunarDay) continue;

        // 除夕/小年: 只有下一月名为"正"才生效
        if (f.onlyIfNextIsZheng && !ctx.nextMonthIsZheng) continue;

        // 闰月排除(除"需要下月为正"那种, 它本身天然只在十二月生效)
        if (ctx.isLunarLeap && !f.onlyIfNextIsZheng) continue;

        // 除夕的 30/29 大小月限制
        if (f.needDay30 && ctx.lunarMonthDays != 30) continue;
        if (f.needDay29 && ctx.lunarMonthDays != 29) continue;

        emitLunarSpec(out, f.spec);
    }

    // 节气作为节日: 清明放假, 其余作为 B 类
    if (ctx.jieQiIdx >= 0 && ctx.jieQiIdx < 24)
    {
        const char *jqName = Jqmc[ctx.jieQiIdx];
        // 清明 = Jqmc[7]
        if (ctx.jieQiIdx == 7)
        {
            appendByType(out, TIER_HOLIDAY, jqName);
        }
        else
        {
            appendByType(out, TIER_MAJOR, jqName);
        }
    }
}

// ═══════════════════════════════════════════════════════════════
//  杂节: 数九 / 三伏 / 入梅 / 出梅
// ═══════════════════════════════════════════════════════════════
void appendMisc(const DayContext &ctx, DayInfo &out)
{
    // 数九: 自冬至开始, 共 81 天
    if (ctx.curDz >= 0 && ctx.curDz < 81)
    {
        int idx = ctx.curDz / 9;        // 0..8
        int rem = ctx.curDz % 9;        // 0..8
        const char *numCn = NumCn[idx + 1]; // 一..九
        if (rem == 0)
        {
            out.major += "『";
            out.major += numCn;
            out.major += "九』 ";
        }
        else
        {
            out.misc += numCn;
            out.misc += "九第";
            out.misc += std::to_string(rem + 1);
            out.misc += "天 ";
        }
    }

    // 三伏: 庚日规则
    // 初伏: 距夏至 [20, 30) 内的庚日
    // 中伏: 距夏至 [30, 40) 内的庚日
    // 末伏: 距立秋 [0, 10) 内的庚日
    bool isGeng = (ctx.dayGanIdx == 6); // 庚 = 第 7 个天干, 索引 6
    bool isBing = (ctx.dayGanIdx == 2); // 丙
    bool isWei  = (ctx.dayZhiIdx == 7); // 未

    if (isGeng)
    {
        if (ctx.curXz >= 20 && ctx.curXz < 30) out.major += "初伏 ";
        if (ctx.curXz >= 30 && ctx.curXz < 40) out.major += "中伏 ";
        if (ctx.curLq >=  0 && ctx.curLq < 10) out.major += "末伏 ";
    }

    // 入梅: 距芒种 [0, 10) 的丙日
    if (isBing && ctx.curMz >= 0 && ctx.curMz < 10) out.major += "入梅 ";
    // 出梅: 距小暑 [0, 12) 的未日
    if (isWei && ctx.curXs >= 0 && ctx.curXs < 12) out.major += "出梅 ";
}

// ═══════════════════════════════════════════════════════════════
//  主入口
// ═══════════════════════════════════════════════════════════════
DayInfo computeAll(const DayContext &ctx)
{
    DayInfo info;
    appendSolarFixed(ctx, info);
    appendSolarWeekly(ctx, info);
    appendLunarFixed(ctx, info);
    appendMisc(ctx, info);
    return info;
}

// 把 Day 上的相关字段抽出来组装为 DayContext.
// 这是 Day → festival 的唯一适配点; Day 类只暴露 getter, festival 模块独占适配逻辑。
DayContext fromDay(Day &d)
{
    DayContext ctx{};
    ctx.year         = d.getSolarYear();
    ctx.month        = d.getSolarMonth();
    ctx.day          = d.getSolarDay();
    ctx.weekday      = d.getWeek();
    ctx.weekiInMonth = d.getWeekIndex() - 1;   // getWeekIndex 从 1 起算
    ctx.firstWeekday = d.getFirstWeekDayOfMonth();
    ctx.weekTotal    = d.getTotalWeekNumsOfMonth();

    ctx.lunarMonth       = d.getLunarMonth();
    ctx.lunarDay         = d.getLunarDay();
    ctx.lunarMonthDays   = d.getLunarMonthDays();
    ctx.isLunarLeap      = d.isLunarLeap();
    ctx.nextMonthIsZheng = d.isNextLunarMonthZheng();

    ctx.jieQiIdx = d.hasJieQi() ? (int)d.getJieQi() : -1;

    ctx.curDz = d.getCurDz();
    ctx.curXz = d.getCurXz();
    ctx.curLq = d.getCurLq();
    ctx.curMz = d.getCurMz();
    ctx.curXs = d.getCurXs();

    GZ rgz = d.getDayGZ();
    ctx.dayGanIdx = rgz.tg;
    ctx.dayZhiIdx = rgz.dz;
    return ctx;
}

DayInfo computeAll(Day &d)
{
    return computeAll(fromDay(d));
}

} // namespace festival
