#include "SSQ.h"
#include <cmath>
#include "const.h"
#include "eph.h"
#include "JD.h"


SSQ::SSQ()
{

	std::string  suoS = "", qiS = "";
	// 619-01-21开始16598个朔日修正表 d0=1947168
	suoS = "EqoFscDcrFpmEsF2DfFideFelFpFfFfFiaipqti1ksttikptikqckstekqttgkqttgkqteksttikptikq2fjstgjqttjkqttgkqt";
	suoS += "ekstfkptikq2tijstgjiFkirFsAeACoFsiDaDiADc1AFbBfgdfikijFifegF1FhaikgFag1E2btaieeibggiffdeigFfqDfaiBkF";
	suoS += "1kEaikhkigeidhhdiegcFfakF1ggkidbiaedksaFffckekidhhdhdikcikiakicjF1deedFhFccgicdekgiFbiaikcfi1kbFibef";
	suoS += "gEgFdcFkFeFkdcfkF1kfkcickEiFkDacFiEfbiaejcFfffkhkdgkaiei1ehigikhdFikfckF1dhhdikcfgjikhfjicjicgiehdik";
	suoS += "cikggcifgiejF1jkieFhegikggcikFegiegkfjebhigikggcikdgkaFkijcfkcikfkcifikiggkaeeigefkcdfcfkhkdgkegieid";
	suoS += "hijcFfakhfgeidieidiegikhfkfckfcjbdehdikggikgkfkicjicjF1dbidikFiggcifgiejkiegkigcdiegfggcikdbgfgefjF1";
	suoS += "kfegikggcikdgFkeeijcfkcikfkekcikdgkabhkFikaffcfkhkdgkegbiaekfkiakicjhfgqdq2fkiakgkfkhfkfcjiekgFebicg";
	suoS += "gbedF1jikejbbbiakgbgkacgiejkijjgigfiakggfggcibFifjefjF1kfekdgjcibFeFkijcfkfhkfkeaieigekgbhkfikidfcje";
	suoS += "aibgekgdkiffiffkiakF1jhbakgdki1dj1ikfkicjicjieeFkgdkicggkighdF1jfgkgfgbdkicggfggkidFkiekgijkeigfiski";
	suoS += "ggfaidheigF1jekijcikickiggkidhhdbgcfkFikikhkigeidieFikggikhkffaffijhidhhakgdkhkijF1kiakF1kfheakgdkif";
	suoS += "iggkigicjiejkieedikgdfcggkigieeiejfgkgkigbgikicggkiaideeijkefjeijikhkiggkiaidheigcikaikffikijgkiahi1";
	suoS += "hhdikgjfifaakekighie1hiaikggikhkffakicjhiahaikggikhkijF1kfejfeFhidikggiffiggkigicjiekgieeigikggiffig";
	suoS += "gkidheigkgfjkeigiegikifiggkidhedeijcfkFikikhkiggkidhh1ehigcikaffkhkiggkidhh1hhigikekfiFkFikcidhh1hit";
	suoS += "cikggikhkfkicjicghiediaikggikhkijbjfejfeFhaikggifikiggkigiejkikgkgieeigikggiffiggkigieeigekijcijikgg";
	suoS += "ifikiggkideedeijkefkfckikhkiggkidhh1ehijcikaffkhkiggkidhh1hhigikhkikFikfckcidhh1hiaikgjikhfjicjicgie";
	suoS += "hdikcikggifikigiejfejkieFhegikggifikiggfghigkfjeijkhigikggifikiggkigieeijcijcikfksikifikiggkidehdeij";
	suoS += "cfdckikhkiggkhghh1ehijikifffffkhsFngErD1pAfBoDd1BlEtFqA2AqoEpDqElAEsEeB2BmADlDkqBtC1FnEpDqnEmFsFsAFn";
	suoS += "llBbFmDsDiCtDmAB2BmtCgpEplCpAEiBiEoFqFtEqsDcCnFtADnFlEgdkEgmEtEsCtDmADqFtAFrAtEcCqAE1BoFqC1F1DrFtBmF";
	suoS += "tAC2ACnFaoCgADcADcCcFfoFtDlAFgmFqBq2bpEoAEmkqnEeCtAE1bAEqgDfFfCrgEcBrACfAAABqAAB1AAClEnFeCtCgAADqDoB";
	suoS += "mtAAACbFiAAADsEtBqAB2FsDqpFqEmFsCeDtFlCeDtoEpClEqAAFrAFoCgFmFsFqEnAEcCqFeCtFtEnAEeFtAAEkFnErAABbFkAD";
	suoS += "nAAeCtFeAfBoAEpFtAABtFqAApDcCGJ";

	//1645-09-23开始7567个节气修正表
	qiS = "FrcFs22AFsckF2tsDtFqEtF1posFdFgiFseFtmelpsEfhkF2anmelpFlF1ikrotcnEqEq2FfqmcDsrFor22FgFrcgDscFs22FgEe";
	qiS += "FtE2sfFs22sCoEsaF2tsD1FpeE2eFsssEciFsFnmelpFcFhkF2tcnEqEpFgkrotcnEqrEtFermcDsrE222FgBmcmr22DaEfnaF22";
	qiS += "2sD1FpeForeF2tssEfiFpEoeFssD1iFstEqFppDgFstcnEqEpFg11FscnEqrAoAF2ClAEsDmDtCtBaDlAFbAEpAAAAAD2FgBiBqo";
	qiS += "BbnBaBoAAAAAAAEgDqAdBqAFrBaBoACdAAf1AACgAAAeBbCamDgEifAE2AABa1C1BgFdiAAACoCeE1ADiEifDaAEqAAFe1AcFbcA";
	qiS += "AAAAF1iFaAAACpACmFmAAAAAAAACrDaAAADG0";



	SB = jieya(suoS);  //定朔修正表解压
	QB = jieya(qiS);   //定气修正表解压


	//朔直线拟合参数
	long double suoKBTmp[] = {
		1457698.231017,29.53067166, // -721-12-17 h=0.00032 古历·春秋
		1546082.512234,29.53085106, // -479-12-11 h=0.00053 古历·战国
		1640640.735300,29.53060000, // -221-10-31 h=0.01010 古历·秦汉
		1642472.151543,29.53085439, // -216-11-04 h=0.00040 古历·秦汉

		1683430.509300,29.53086148, // -104-12-25 h=0.00313 汉书·律历志(太初历)平气平朔
		1752148.041079,29.53085097, //   85-02-13 h=0.00049 后汉书·律历志(四分历)
		1807665.420323,29.53059851, //  237-02-12 h=0.00033 晋书·律历志(景初历)
		1883618.114100,29.53060000, //  445-01-24 h=0.00030 宋书·律历志(何承天元嘉历)
		1907360.704700,29.53060000, //  510-01-26 h=0.00030 宋书·律历志(祖冲之大明历)
		1936596.224900,29.53060000, //  590-02-10 h=0.01010 随书·律历志(开皇历)
		1939135.675300,29.53060000, //  597-01-24 h=0.00890 随书·律历志(大业历)
		1947168.00//  619-01-21
	};

	suoKB = new std::vector<long double>(suoKBTmp, suoKBTmp + sizeof(suoKBTmp) / sizeof(long double));

	long double qiKBTmp[] = {
		1640650.479938,15.21842500, // -221-11-09 h=0.01709 古历·秦汉
		1642476.703182,15.21874996, // -216-11-09 h=0.01557 古历·秦汉

		1683430.515601,15.218750011, // -104-12-25 h=0.01560 汉书·律历志(太初历)平气平朔 回归年=365.25000
		1752157.640664,15.218749978, //   85-02-23 h=0.01559 后汉书·律历志(四分历) 回归年=365.25000
		1807675.003759,15.218620279, //  237-02-22 h=0.00010 晋书·律历志(景初历) 回归年=365.24689
		1883627.765182,15.218612292, //  445-02-03 h=0.00026 宋书·律历志(何承天元嘉历) 回归年=365.24670
		1907369.128100,15.218449176, //  510-02-03 h=0.00027 宋书·律历志(祖冲之大明历) 回归年=365.24278
		1936603.140413,15.218425000, //  590-02-17 h=0.00149 随书·律历志(开皇历) 回归年=365.24220
		1939145.524180,15.218466998, //  597-02-03 h=0.00121 随书·律历志(大业历) 回归年=365.24321
		1947180.798300,15.218524844, //  619-02-03 h=0.00052 新唐书·历志(戊寅元历)平气定朔 回归年=365.24460
		1964362.041824,15.218533526, //  666-02-17 h=0.00059 新唐书·历志(麟德历) 回归年=365.24480
		1987372.340971,15.218513908, //  729-02-16 h=0.00096 新唐书·历志(大衍历,至德历) 回归年=365.24433
		1999653.819126,15.218530782, //  762-10-03 h=0.00093 新唐书·历志(五纪历) 回归年=365.24474
		2007445.469786,15.218535181, //  784-02-01 h=0.00059 新唐书·历志(正元历,观象历) 回归年=365.24484
		2021324.917146,15.218526248, //  822-02-01 h=0.00022 新唐书·历志(宣明历) 回归年=365.24463
		2047257.232342,15.218519654, //  893-01-31 h=0.00015 新唐书·历志(崇玄历) 回归年=365.24447
		2070282.898213,15.218425000, //  956-02-16 h=0.00149 旧五代·历志(钦天历) 回归年=365.24220
		2073204.872850,15.218515221, //  964-02-16 h=0.00166 宋史·律历志(应天历) 回归年=365.24437
		2080144.500926,15.218530782, //  983-02-16 h=0.00093 宋史·律历志(乾元历) 回归年=365.24474
		2086703.688963,15.218523776, // 1001-01-31 h=0.00067 宋史·律历志(仪天历,崇天历) 回归年=365.24457
		2110033.182763,15.218425000, // 1064-12-15 h=0.00669 宋史·律历志(明天历) 回归年=365.24220
		2111190.300888,15.218425000, // 1068-02-15 h=0.00149 宋史·律历志(崇天历) 回归年=365.24220
		2113731.271005,15.218515671, // 1075-01-30 h=0.00038 李锐补修(奉元历) 回归年=365.24438
		2120670.840263,15.218425000, // 1094-01-30 h=0.00149 宋史·律历志 回归年=365.24220
		2123973.309063,15.218425000, // 1103-02-14 h=0.00669 李锐补修(占天历) 回归年=365.24220
		2125068.997336,15.218477932, // 1106-02-14 h=0.00056 宋史·律历志(纪元历) 回归年=365.24347
		2136026.312633,15.218472436, // 1136-02-14 h=0.00088 宋史·律历志(统元历,乾道历,淳熙历) 回归年=365.24334
		2156099.495538,15.218425000, // 1191-01-29 h=0.00149 宋史·律历志(会元历) 回归年=365.24220
		2159021.324663,15.218425000, // 1199-01-29 h=0.00149 宋史·律历志(统天历) 回归年=365.24220
		2162308.575254,15.218461742, // 1208-01-30 h=0.00146 宋史·律历志(开禧历) 回归年=365.24308
		2178485.706538,15.218425000, // 1252-05-15 h=0.04606 淳祐历 回归年=365.24220
		2178759.662849,15.218445786, // 1253-02-13 h=0.00231 会天历 回归年=365.24270
		2185334.020800,15.218425000, // 1271-02-13 h=0.00520 宋史·律历志(成天历) 回归年=365.24220
		2187525.481425,15.218425000, // 1277-02-12 h=0.00520 本天历 回归年=365.24220
		2188621.191481,15.218437494, // 1280-02-13 h=0.00015 元史·历志(郭守敬授时历) 回归年=365.24250
		2322147.76// 1645-09-21
	};

	qiKB = new std::vector<long double>(qiKBTmp, qiKBTmp + sizeof(qiKBTmp) / sizeof(long double));
}

SSQ::~SSQ()
{
	delete suoKB;
	delete qiKB;
}

void str_replace(std::string & str, const std::string strsrc, const std::string strdst)
{
	std::string::size_type pos = 0;//位置
	std::string::size_type srclen = strsrc.size();//要替换的字符串大小
	std::string::size_type dstlen = strdst.size();//目标字符串大小
	while ((pos = str.find(strsrc, pos)) != std::string::npos)
	{
		str.replace(pos, srclen, strdst);
		pos += dstlen;
	}
}


std::string SSQ::jieya(std::string s) { //气朔解压缩
	std::string o = "0000000000", o2 = o + o;
	str_replace(s, "J", "00");
	str_replace(s, "I", "000");
	str_replace(s, "H", "0000");
	str_replace(s, "G", "00000");
	str_replace(s, "t", "02");
	str_replace(s, "s", "002");
	str_replace(s, "r", "0002");
	str_replace(s, "q", "00002");
	str_replace(s, "p", "000002");
	str_replace(s, "o", "0000002");
	str_replace(s, "n", "00000002");
	str_replace(s, "m", "000000002");
	str_replace(s, "l", "0000000002");
	str_replace(s, "k", "01");
	str_replace(s, "j", "0101");
	str_replace(s, "i", "001");
	str_replace(s, "h", "001001");
	str_replace(s, "g", "0001");
	str_replace(s, "f", "00001");
	str_replace(s, "e", "000001");
	str_replace(s, "d", "0000001");
	str_replace(s, "c", "00000001");
	str_replace(s, "b", "000000001");
	str_replace(s, "a", "0000000001");
	str_replace(s, "A", o2 + o2 + o2);
	str_replace(s, "B", o2 + o2 + o);
	str_replace(s, "C", o2 + o2);
	str_replace(s, "D", o2 + o);
	str_replace(s, "E", o2);
	str_replace(s, "F", o);
	return s;
}

int  SSQ::calc(long double jd, QSType qs)
{
	jd += J2000;
	int i;
	long double D;
	std::string n;
	std::vector<long double> B = *suoKB;
	long double pc = 14;
	//如果查的是气朔
	if (qs == QType)
	{
		B = *qiKB, pc = 7;
	}

	long double f1 = B[0] - pc, f2 = B[B.size() - 1] - pc, f3 = 2436935;   // 2436935为1960年1月1日JD

	if (jd < f1 || jd >= f3)
	{
		//平气朔表中首个之前，使用现代天文算法。1960.1.1以后，使用现代天文算法 (这一部分调用了qi_high和so_high,所以需星历表支持)
		if (qs == QType)
		{
			return  std::floor(qi_high(std::floor((jd + pc - 2451259) / 365.2422 * 24) * PI / 12) + 0.5); //2451259是1999.3.21,太阳视黄经为0,春分.定气计算
		}
		else
		{
			return  std::floor(so_high(std::floor((jd + pc - 2451551) / 29.5306) * PI * 2) + 0.5); //2451551是2000.1.7的那个朔日,黄经差为0.定朔计算
		}
	}

	if (jd >= f1 && jd < f2)
    { //平气或平朔
		for (i = 0; i < B.size(); i += 2)
        {
            if (jd + pc < B[i + 2])
            {
                break;
            }
        }

		D = B[i] + B[i + 1] * std::floor((jd + pc - B[i]) / B[i + 1]);
		D = std::floor(D + 0.5);
		if (D == 1683460)
        {
            D++; //如果使用太初历计算-103年1月24日的朔日,结果得到的是23日,这里修正为24日(实历)。修正后仍不影响-103的无中置闰。如果使用秦汉历，得到的是24日，本行D不会被执行。
        }
		return D - J2000;
	}

	if (jd >= f2 && jd < f3)
    { //定气或定朔
		if (qs == QType)
        {
			D = std::floor(qi_low(std::floor((jd + pc - 2451259) / 365.2422 * 24) * PI / 12) + 0.5); //2451259是1999.3.21,太阳视黄经为0,春分.定气计算
			n = QB.substr(std::floor((jd - f2) / 365.2422 * 24), 1); //找定气修正值
		}
		else {
			D = std::floor(so_low(std::floor((jd + pc - 2451551) / 29.5306) * PI * 2) + 0.5); //2451551是2000.1.7的那个朔日,黄经差为0.定朔计算
			n = SB.substr(std::floor((jd - f2) / 29.5306), 1); //找定朔修正值
		}
		if (n == "1")
        {
            return D + 1;
        }
		if (n == "2")
        {
            return D - 1;
        }

		return D;
	}
    return 0;
}


long double SSQ::qi_high(long double W)
{
	long double t = XL::S_aLon_t2(W) * 36525;
	t = t - dt_T(t) + 8.0 / 24.0;
	long double v = fmod(t + 0.5, 1) * 86400;
	if (v < 1200 || v >86400 - 1200) t = XL::S_aLon_t(W) * 36525 - dt_T(t) + 8.0 / 24.0;
	return  t;
}


long double SSQ::so_high(long double W)
{ //较高精度朔
	long double t = XL::MS_aLon_t2(W) * 36525;
	t = t - dt_T(t) + 8.0 / 24.0;
	long double v = fmod(t + 0.5, 1) * 86400;
	if (v < 1800 || v >86400 - 1800) t = XL::MS_aLon_t(W) * 36525 - dt_T(t) + 8.0 / 24.0;
	return  t;
}

long double SSQ::so_low(long double W) { //低精度定朔计算,在2000年至600，误差在2小时以内(仍比古代日历精准很多)
	long double v = 7771.37714500204;
	long double t = (W + 1.08472) / v, L;
	t -= (-0.0000331*t*t
		+ 0.10976 *cosl(0.785 + 8328.6914*t)
		+ 0.02224 *cosl(0.187 + 7214.0629*t)
		- 0.03342 *cosl(4.669 + 628.3076*t)) / v
		+ (32 * (t + 1.8)*(t + 1.8) - 20) / 86400 / 36525;
	return t * 36525 + 8.0 / 24.0;
}

long double SSQ::qi_low(long double W) { //最大误差小于30分钟，平均5分
	long double t, L, v = 628.3319653318;
	t = (W - 4.895062166) / v; //第一次估算,误差2天以内
	t -= (53 * t*t + 334116 * cosl(4.67 + 628.307585*t) + 2061 * cosl(2.678 + 628.3076*t)*t) / v / 10000000; //第二次估算,误差2小时以内

	L = 48950621.66 + 6283319653.318*t + 53 * t*t //平黄经
		+ 334166 * cosl(4.669257 + 628.307585*t) //地球椭圆轨道级数展开
		+ 3489 * cosl(4.6261 + 1256.61517*t) //地球椭圆轨道级数展开
		+ 2060.6 * cosl(2.67823 + 628.307585*t) * t  //一次泊松项
		- 994 - 834 * sinl(2.1824 - 33.75705*t); //光行差与章动修正

	t -= (L / 10000000 - W) / 628.332 + (32 * (t + 1.8)*(t + 1.8) - 20) / 86400 / 36525;
	return t * 36525 + 8.0 / 24.0;
}

bool isLeapYear(int year)
{
    int x = year % 4;
    int y = year % 100;
    int z = year % 400;

    return !x && (y || !z);
}
/*
 *
 *

中气不足话闰月
assistant
中气不足话闰月是指在中国农历中，当某一年的太阳黄经度数不足365.2422度时，需要增加一个月份来弥补这个差距，这个月份被称为闰月。闰月的出现是为了使农历与太阳历保持一定的同步，以便更好地反映季节的变化。

在农历中，每个月的开始是以月亮的新月为准的，所以每个农历月的天数约为29.53天，12个农历月的天数约为354.36天，而太阳历一年的天数约为365.2422天。为了弥补这两者之间的差距，大约每2.7年就需要增加一个闰月。

闰月的安排是根据天文学家观测月亮和太阳的运行情况来确定的，具体方法是：在19年中，安排7个闰月，这样每个19年的周期中，农历的月数总共为235个月，而太阳历的月数为228个月。通过这样的安排，农历与太阳历之间的差距就能够得到弥补，使农历更好地反映季节的变化。

 平气，定气:
平气和定气是中国古代天文学家用来描述太阳在黄道上运行的两个概念，它们都与二十四节气有关。

平气：平气是指太阳在黄道上平均运行的速度。古代天文学家将一年分为24个节气，每个节气的时间间隔相等，即平均每15天一个节气。这样，太阳在黄道上的平均运行速度为360°/365.2422 ≈ 0.9856°/天。这种假设下的太阳运行速度被称为平气。

定气：定气是指太阳在黄道上实际运行的速度。由于地球公转轨道呈椭圆形状，太阳在黄道上的实际运行速度是不均匀的。因此，实际上每个节气之间的时间间隔是不相等的。古代天文学家根据实际观测数据，计算出太阳在黄道上各个位置的实际运行速度，这种实际的太阳运行速度被称为定气。

总之，平气和定气都是描述太阳在黄道上运行速度的概念，其中平气是假设太阳在黄道上运行速度均匀，而定气是根据实际观测数据计算出的太阳在黄道上的实际运行速度。在二十四节气的安排上，采用定气能更准确地反映季节变化。
 */

//农历排月序计算,可定出农历,有效范围：两个冬至之间(冬至一 <= mDay < 冬至二)

/*
 * 我们不妨把“正月、二月、三月等”看作月建的别名。并且建子为十一，建丑为十二，建寅为正······由于某些原因，有时古代帝王随意更改了别名。本万年历也相应作了更改:
-721年 12 月 17 日至-479 年11 月春秋战国历采用 19年7，闰年的末月置闰并取名闰十三。年首为正月。
-221年 10月31日至-104年11 月秦汉历同样采用 19年7闰，闰年的末月置闰并取名后九月。年首为十月。
以上两行所述的年代，其月建别名计算方法:
月建序数 =(月积累数 - 置闰积累数 + c) mod 12，式中C为某一常数
008年01月15日至 023年12月02 日:建子为十二，建丑为正月，建寅为二月，其它顺推。
237年04月12 日至239年12月13日:建子为十二，建丑为正月，建寅为二月，其它顺推。
689 年12月 18日至701年01月 14日:建子为正，建寅为一月，其它不变761年12 月 02 日至 762年03 月30 日:建子为正，其它顺推。
由于上述月建别名的更改，造成排算的结果是 239 年 12.13 为十二月，240年01月 12 日变回原来的月建别名，这样240年01月 12 日建丑十二月，显然连续的出现两个十二月，本万年历中把上个月表示为“拾贰月”。同样的问题也出现在 23 年底。
237 年正月开始使用新的历法，造成前一月只有28 天。
 */
void SSQ::calcY(int jd) {
	std::vector<long double>& A = zhongqi_array_;
	std::vector<int>& B = sun_moon_hesuo_array_;  //中气表,日月合朔表(整日)
    std::vector<int>& Pe = zhongqi_pe_array_;
	int i, k;
	long double W, w;

	//该年的气
	W = int2((jd - 355 + 183) / 365.2422)*365.2422 + 355;  //355是2000.12冬至,得到较靠近jd的冬至估计值
	if (calc(W, QType) > jd)
	{
		W -= 365.2422;
	}

	//25个节气时刻(北京时间),从冬至开始到下一个冬至以后;
	A.clear();
	for (i = 0; i < 25; i++)
	{	
		int t = calc(W + 15.2184*i, QType);
		A.push_back(t);
	}

    //补算二气,确保一年中所有月份的“气”全部被计算在内
	Pe.clear();
	auto pe1 = calc(W - 15.2, QType);
    auto pe2 = calc(W - 30.4, QType);
    Pe.push_back(pe1);
    Pe.push_back(pe2);

    //今年"首朔"的日月黄经差w
	w = calc(A[0], SType); //求较靠近冬至的朔日
	if (w > A[0])
	{
		w -= 29.53;  //朔望月：每节气
	}

	//该年所有朔,包含14个月的始末
	B.clear();
	for (i = 0; i < 15; i++)
	{
		B.push_back( calc(w + 29.5306*i, SType) );
	}
		

	//月大小
	leap_month_ = 0;
	month_daxiao_array_.clear();
	month_order_array_.clear();
	month_name_array_.clear();
    specific_next_month_array_.clear();
	for (i = 0; i < 14; i++) {
		month_daxiao_array_.push_back(sun_moon_hesuo_array_[i + 1] - sun_moon_hesuo_array_[i] ); //月大小
		month_order_array_.push_back(i);  //月序初始化
		month_name_array_.push_back(0);
        specific_next_month_array_.push_back(false);
	}


	//-721年至-104年的后九月及月建问题,与朔有关，与气无关，19年7闰法”(-721年至-104年)，7个闰年均匀安插在19个年整数中，闰年的末月置为闰月
    /*
     * -721年12月17日至-479年11月春秋战国历采用19年7闰,闰年的末月置闰并取名闰十三。年首为正月。
       -221年10月31日至-104年11月秦汉历同样采用19年7闰，闰年的末月置闰并取名后九月。年首为十月。
       以上两行所述的年代，其月建别名计算方法:
       月建序数 =(月积累数 - 置闰积累数 + c) mod 12，式中C为某一常数
     */
    /*
     * 在中国，除了格列高里历（俗称阳历）之外，还有盛行千百年之久的农历法。这是一种特殊的阴阳历，并不是单纯的阴历。中国的百姓到现在仍然以它为依据，安排农事、渔业生产以及确定传统节日。

      农历是按朔望周期确定月份。月相朔（日月合朔）所在的日期为本月初一，下次朔的日期为下月初一。因为一个朔望月的周期是29.53天，所以分大月和小月，大月30天，小月29天。某月是“大”还是“小”，以及哪天是“朔日”，则根据太阳和月亮的真实位置来推算，古时称为“定朔”。

      农历的年以回归年为依据。为了和回归年的长度相似，农历使用增加闰月的方法（根据二十四节气制定），并将岁首调整到“雨水”所在的月初。农历一年12个月，一共是354日或者355日，平均19年有7个闰月，这样就保证了19年的农历与19年的回归年的长度基本相等。所以通常情况下，中国人的19岁、38岁、57岁及76岁时的阳历生日和农历生日会重合在一起。

      汉武帝太初元年（公元前104年）五月颁布的《太初历》，将含有雨水的月份定为正月，将这个月的初一定为岁首。因其更加科学地反映农业季节，除个别朝代有过短期改动外，一直沿用至今。
     * */

    /*
     * 代码首先确定了当前年份YY，然后根据不同历法的规则，计算了三个历年的首朔日（ns数组中的前三个元素），以及闰月名称（ns数组中的后三个元素）和月建（ns数组中的后六个元素）。

      接下来，代码遍历了14个朔日，计算每个月的积数（即该月第几个朔日），然后根据当前历法的规则，判断该月是否为闰月。如果不是闰月，则根据积数计算出该月的名称。如果是闰月，则根据闰月名称和积数计算出该月的名称。

      具体来说，代码计算闰月的方法是：首先找到当前月份所在的历年的首朔日（即ns数组中最后一个小于当前月份的元素），然后计算当前月份的积数。如果当前月份的积数小于12，则该月不是闰月，直接根据积数和月建计算出该月的名称。如果当前月份的积数等于12，则该月是闰月，名称为闰月名称。如果当前月份的积数大于12，则需要将积数减去1，然后再根据积数和月建计算出该月的名称。

      需要注意的是，这段代码只适用于-721年至-104年这段时间内的历法，对于其他时间段的历法，可能需要使用不同的计算方法。
     */
    //该行代码计算出当前历法的年份YY，其中this.ZQ[0]表示当前历法的立春时刻，int2函数表示取整操作，365.2422表示一年的平均天数，+10 +180表示将计算基准点从J2000调整到2000年春分前10天。
	int YY = int2((zhongqi_array_[0] + 10 + 180) / 365.2422) + 2000; //确定年份
	static const int yueIdx[13] =  {11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,1 };
	bd_ns_array_.clear();
	if (YY >= -721 && YY <= -104)
    {
        //创建一个新数组ns，用于存储历年首朔日、闰月名称和月建。yy变量用于存储当前遍历到的年份。使用for循环遍历3个历年
		int ns[12]={0,0,0,0,0,0,0,0,0,0,0,0};
		int yy;

        //根据不同的历法，计算当前历年的首朔日（存储在ns[i]中）、闰月名称（存储在ns[i+3]中）和月建（存储在ns[i+6]中）。
        // calc函数用于计算朔日，第一个参数是儒略日，第二个参数是计算方式（这里是'朔'，表示计算朔日）。
		for (i = 0; i < 3; i++)
        {
            //计算当前遍历到的年份yy。
			yy = YY + i - 1;
			//颁行历年首, 闰月名称, 月建
			if (yy >= -721)
            {
                ns[i] = calc(1457698 - J2000 + int2(0.342 + (yy + 721)*12.368422)*29.5306, SType);
                //ns[i + 3] = '十三';
                ns[i + 3] = 1;
                ns[i + 6] = 2;  //春秋历,ly为-722.12.17
                //leap_month_ = 2;
                //month_order_array_={0,1,1,2,3,4,5,6,7,8,9,10,11,12};
            }
            
			if (yy >= -479)
            {
                ns[i] = calc(1546083 - J2000 + int2(0.500 + (yy + 479)*12.368422)*29.5306, SType);
                //ns[i + 3] = '十三';
                ns[i + 3] = 1;
                ns[i + 6] = 2;  //战国历,ly为-480.12.11
                //leap_month_ = 2;
                //month_order_array_={0,1,1,2,3,4,5,6,7,8,9,10,11,12};
            }
            if (yy >= -220) {
                ns[i] = calc(1640641 - J2000 + int2(0.866 + (yy + 220)*12.369000)*29.5306, SType);
                //ns[i + 3] = '后九';  -220，10，19
                //-220年月序[("十二", "正", "二", "三", "四", "五", "六", "七", "八", "九", "后九"，"十"，"十一",)]
                ns[i + 3] = 9;
                ns[i + 6] = 11; //秦汉历,ly为-221.10.31
                //leap_month_ = 11;
                //month_order_array_={0,1,1,2,3,4,5,6,7,8,9,10,11,12};
            }
		}

		bd_ns_array_.push_back(ns[1]);
		bd_ns_array_.push_back(ns[2]);
		bd_ns_array_.push_back(ns[4]);

        //循环遍历14个朔日
		int nn, f1;
		for (i = 0; i < 14; i++)
        {
            int Dm = sun_moon_hesuo_array_[i] + J2000;
            //找到当前朔日所在的历年首朔日（即ns数组中最后一个小于当前朔日的元素）
			for (nn = 2; nn >= 0; nn--)
            {
                if (sun_moon_hesuo_array_[i] >= ns[nn])
                {
                    break;
                }
            }

            //月建序数 =(月积累数 - 置闰积累数 + c) mod 12，式中C为某一常数
            //计算当前朔日的积数（即该月第几个朔日）
			f1 = int2((sun_moon_hesuo_array_[i] - ns[nn] + 15) / 29.5306); //该月积数

            //根据积数和历年首朔日，判断当前朔日是否为闰月。如果是闰月，则使用闰月名称（存储在ns[nn+3]中）；
            // 如果不是闰月，则根据积数和月建（存储在ns[nn+6]中）计算出该月的名称
			if (f1 < 12)
            {
                month_order_array_[i] = (f1 + ns[nn + 6]) % 12;
            }
            else
            {
                month_order_array_[i] = ns[nn + 3]; //闰月
                leap_month_ = i;
            }
			month_name_array_[i] = yueIdx[month_order_array_[i]];

            //-221年10月31日，是后一个十月;-221-11-29: 是后一个冬月
            if (1640641 == Dm || 1640670 == Dm )
            {
                specific_next_month_array_[i] = true;
            }
		}

		return;
	}


	//无中气置闰法确定闰月,(气朔结合法,数据源需有冬至开始的的气和朔)
	if (B[13] <= A[24]) { //第13月的月末没有超过冬至(不含冬至),说明今年含有13个月
		for (i = 1; B[i + 1] > A[2 * i] && i < 13; i++); //在13个月中找第1个没有中气的月份
		leap_month_ = i;
		for (; i < 14; i++)
        {
            month_order_array_[i]--;
        }
	}

	//名称转换(月建别名)
    /*
     * 008年01月15日至 023年12月02 日:建子为十二，建丑为正月，建寅为二月，其它顺推。
       237年04月12 日至239年12月13日:建子为十二，建丑为正月，建寅为二月，其它顺推。
       689 年12月 18日至701年01月 14日:建子为正，建寅为一月，其它不变761年12 月 02 日至 762年03 月30 日:建子为正，其它顺推。
       由于上述月建别名的更改，造成排算的结果是 239 年 12.13 为十二月，240年01月 12 日变回原来的月建别名，这样240年01月 12 日建丑十二月，
       显然连续的出现两个十二月，本万年历中把上个月表示为“拾贰月”。同样的问题也出现在 23 年底。
       237 年正月开始使用新的历法，造成前一月只有28 天。
     */
	for (i = 0; i < 14; i++)
    {
		int Dm = sun_moon_hesuo_array_[i] + J2000, v2 = month_order_array_[i]; //Dm初一的儒略日,v2为月建序号
		int mc = v2 % 12; //月建对应的默认月名称：建子十一,建丑十二,建寅为正……
		if (Dm >= 1724360 && Dm <= 1729794)
        {
            mc = (v2 + 1) % 12; //  8.01.15至 23.12.02 建子为十二,其它顺推
        }
		else if (Dm >= 1807724 && Dm <= 1808699)
        {
            mc = (v2 + 1) % 12; //237.04.12至239.12.13 建子为十二,其它顺推
        }
		else if (Dm >= 1999349 && Dm <= 1999467)
        {
            mc = (v2 + 2) % 12; //761.12.02至762.03.30 建子为正月,其它顺推
        }
		else if (Dm >= 1973067 && Dm <= 1977052) //689-12-18~~700-11-15
        {
            if (v2 % 12 == 0)
            {
                mc = 2;
            }

            if (v2 == 2)
            {
                mc = 12; //689.12.18至700.11.15 建子为正月,建寅为一月,其它不变
                specific_next_month_array_[i] = true;
            }
        }

        //23-12-02 12:00:00,239-12-13 12:00:00
		if (Dm == 1729794 || Dm == 1808699)
        {
			mc = 1;
			//specific_next_month_array_[i] = true;
           // mc = '拾贰'; //239.12.13及23.12.02均为十二月,为避免两个连续十二月，此处改名
        }

        
        //1729823: 23-12-31 12:00:00; 1808729: 240-01-12 12:00:00
        if (Dm == 1729823 || Dm == 1808729)
        {
            //后一个十二月
            mc = 1;
            specific_next_month_array_[i] = true;
        }
	

        //1999497: 762-04-29 12:00:00, 1999526: 762-05-28 12:00:00; 1977112: 701-01-14 12:00:00
        if (Dm == 1999497 || Dm == 1999526 || Dm == 1977112)
        {
            //恢复建寅了，之前可能是建子，762年恢复了旧历
            specific_next_month_array_[i] = true;
        }

        month_order_array_[i] = mc;
		month_name_array_[i] = yueIdx[mc];
	}

	return;
}

std::vector<long double> SSQ::getZhongQi() const
 {
    return zhongqi_array_;
}

int SSQ::getLeap() const  
{
    return leap_month_;
}

std::vector<int> SSQ::getHS() const
{
    return sun_moon_hesuo_array_;
}

std::vector<int> SSQ::getYm() const
{
    return month_order_array_;
}

std::vector<int> SSQ::getYueName() const
{
    return month_name_array_;
}

//月大小
std::vector<int> SSQ::getDx() const
{
    return month_daxiao_array_;
}

//有些年份有两个特殊的农历月份
std::vector<bool> SSQ::getSpecificLunarMonth() const
{
    return specific_next_month_array_;
}

std::vector<int> SSQ::getBDNS() const
{
    return bd_ns_array_;
}

std::vector<int> SSQ::getZQPe() const {
    return zhongqi_pe_array_;
}
