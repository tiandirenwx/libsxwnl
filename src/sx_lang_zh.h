#pragma once
#include <string>

static const char *Gan[] = {"甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"};
static const char *Zhi[] = {"子", "丑", "寅", "卯", "辰", "巳", "午", "未", "申", "酉", "戌", "亥"};
static const char *ShengXiao[] = {"鼠", "牛", "虎", "兔", "龙", "蛇", "马", "羊", "猴", "鸡", "狗", "猪"};
static const char *NumCn[] = {"零", "一", "二", "三", "四", "五", "六", "七", "八", "九", "十"}; //中文数字
static const char *Jqmc[] = {"冬至", "小寒", "大寒", "立春", "雨水", "惊蛰", "春分", "清明", "谷雨", "立夏", "小满", "芒种", "夏至", "小暑", "大暑", "立秋", "处暑", "白露", "秋分", "寒露", "霜降", "立冬", "小雪", "大雪"};
static const char *JieLing[ ]= {"大雪", "小寒", "立春", "惊蛰", "清明", "立夏", "芒种", "小暑", "立秋", "白露", "寒露", "立冬", "大雪"};
static const char *Ymc[] = { "正", "二", "三", "四", "五", "六", "七", "八", "九", "十","冬", "腊"}; //月名称,建寅
static const char *SYmc[] = { "壹", "贰", "叁", "肆", "伍", "陆", "柒", "捌", "玖", "拾","拾壹", "拾贰"}; //特殊月名称，区别正常月
static const char *Rmc[] = {"初一", "初二", "初三", "初四", "初五", "初六", "初七", "初八", "初九", "初十", "十一", "十二", "十三", "十四", "十五", "十六", "十七", "十八", "十九", "二十", "廿一", "廿二", "廿三", "廿四", "廿五", "廿六", "廿七", "廿八", "廿九", "三十"};
static const char *WeekCn[] = {"星期日", "星期一", "星期二", "星期三", "星期四", "星期五", "星期六"};
static const char *XiZ[] = {"摩羯", "水瓶", "双鱼", "白羊", "金牛", "双子", "巨蟹", "狮子", "处女", "天秤", "天蝎", "射手"};
static const char *RiJian12[] = {"建", "除", "满", "平", "定", "执", "破", "危", "成", "收", "开", "闭"};
static const char *BDLeapYueName[] = {"十三","后九"};
static const char *Yuedx[]  ={"小","大"};
static const char *DayAgo[] ={"昨天","今天","明天","天后","天前"};
static const char *Fangw[]  ={"南","东","北","西"};
static const char *DayTime[] ={"凌晨","早晨","上午","中午","下午","晚上","深夜"};
static const char *LeapYear[]  ={"  ","闰"};
static const char *LeapYear2[] ={"","闰"};

static const char *NaYinWuXing[]={"海中金","炉中火","大林木","路边土","剑峰金","山头火","涧下水","城墙土","白腊金","杨柳木","泉中水","屋上土","霹雷火","松柏木","长流水","沙中金","山下火","平地木","壁上土","金泊金","佛灯火","天河水","大泽土","钗钏金","桑松木", "大溪水","沙中土", "天上火","石榴木","大海水"};
static const char *XingQiName[]={"日", "一", "二", "三", "四", "五", "六"};

static const char *YueXiangName[]={"朔","上弦","望","下弦"};//月相名称表

static const char *ShiShen[] = {"比肩", "劫财", "食神", "伤官", "偏财", "正财", "七杀", "正官", "偏印", "正印", ""};

static const char *Gender[] = {"乾造", "坤造"};
static const char *LifaList[] = {"太初历", "四分历", "大明历", "戊寅元历", "麟德历", "正元历", "应天历", "崇天历", "淳祐历", "授时历","现代农历-定气法", "尤武伟子平历-定冬至","尤武伟子平历-定夏至"};
//地支藏干索引
static const  int gCangGan[][3] = {
    {9, 10, 10}, {5, 7, 9}, {0, 2, 4}, {1, 10, 10}, {1, 4, 9}, {2, 4, 6},
    {3, 5, 10}, {1, 3, 5}, {4, 6, 8}, {7, 10, 10}, {3, 4, 7}, {0, 8, 10}
    };

static const char *gSiZhuList[][5] = {
	//#甲己   乙庚   丙辛   丁壬   戊癸
	{"甲子", "丙子", "戊子", "庚子", "壬子"},
	{"乙丑", "丁丑", "己丑", "辛丑", "癸丑"},
	{"丙寅", "戊寅", "庚寅", "壬寅", "甲寅"},
	{"丁卯", "己卯", "辛卯", "癸卯", "乙卯"},
	{"戊辰", "庚辰", "壬辰", "甲辰", "丙辰"},
	{"己巳", "辛巳", "癸巳", "乙巳", "丁巳"},
	{"庚午", "壬午", "甲午", "丙午", "戊午"},
	{"辛未", "癸未", "乙未", "丁未", "己未"},
	{"壬申", "甲申", "丙申", "戊申", "庚申"},
	{"癸酉", "乙酉", "丁酉", "己酉", "辛酉"},
	{"甲戌", "丙戌", "戊戌", "庚戌", "壬戌"},
	{"乙亥", "丁亥", "己亥", "辛亥", "癸亥"}};

static const char *EphemerisList[] ={"地球","水星","金星","火星","木星","土星","天王星","海王星","冥王星","太阳","月亮"};
//std::array<std::string> xxName = {"地球","水星","金星","火星","木星","土星","天王星","海王星","冥王星"};


