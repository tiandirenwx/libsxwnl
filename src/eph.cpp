#include "eph.h"
#include "const.h"
#include <cmath>
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>

//=================================数学工具=========================================
// 对超过0-2PI的角度转为0-2PI;
long double rad2mrad(long double v)
{
	v = fmod(v, (2 * PI));
	if (v < 0)
	{
		return v + 2 * PI;
	}
	return v;
}

// 对超过-PI到PI的角度转为-PI到PI;
long double rad2rrad(long double v)
{
	v = fmod(v, (2 * PI));
	if (v <= -PI)
	{
		return v + 2 * PI;
	}
	if (v > PI)
	{
		return v - 2 * PI;
	}
	return v;
}

// 将弧度转为字串,保留2位
std::string rad2str(long double d, bool tim)
{
	return rad2strE(d, tim, 2);
}

//===============角度格式化==================
std::string rad2strE(long double d, bool flag, int ext)
{
	// 将弧度转为字串,ext为小数保留位数
	// flag=0输出格式示例: -23°59" 48.23"
	// flag=1输出格式示例:  18h 29m 44.52s
	std::string s = " ", w1 = "°", w2 = "\'", w3 = "\"";
	if (d < 0)
	{
		d = -d, s = "-";
	}
	if (flag)
	{
		d *= 12.0 / M_PI;
		w1 = "h", w2 = "m", w3 = "s";
	}
	else
	{
		d *= 180.0 / M_PI;
	}
	auto a = int64(d);
	d = (d - a) * 60;
	auto b = int64(d);
	d = (d - b) * 60;
	auto c = int64(d);

	long double Q = pow(10, ext);

	d = int64((d - c) * Q + 0.5);
	if (d >= Q)
	{
		d -= Q, c++;
	}
	if (c >= 60)
	{
		c -= 60, b++;
	}
	if (b >= 60)
	{
		b -= 60, a++;
	}

	std::string A, B, C, D;
	A = "   " + std::to_string(a);
	B = "0" + std::to_string(b);
	C = "0" + std::to_string(c);
	D = "00000" + std::to_string((int)d);
	s += A.substr(A.length() - 3, 3) + w1;
	s += B.substr(B.length() - 2, 2) + w2;
	s += C.substr(C.length() - 2, 2);
	if (ext)
	{
		s += "." + D.substr(D.length() - ext, ext) + w3;
	}
	return s;
}

// 将弧度转为字串,精确到分
std::string rad2str2(long double d)
{
	// 输出格式示例: -23°59"
	std::string s = "+", w1 = "°", w2 = "\'", w3 = "\"";
	if (d < 0)
	{
		d = -d, s = "-";
	}

	d *= 180 / M_PI;
	auto a = int64(d);
	auto b = int64((d - a) * 60 + 0.5);
	if (b >= 60)
	{
		b -= 60, a++;
	}
	std::string A = "   " + std::to_string(a), B = "0" + std::to_string(b);
	s += A.substr(A.length() - 3, 3) + w1;
	s += B.substr(B.length() - 2, 2) + w2;
	return s;
}

// 临界余数(a与最近的整倍数b相差的距离);
long double mod2(long double a, long double b)
{
	long double c = a / b;
	c -= floor(c);
	if (c > 0.5)
	{
		c -= 1;
	}
	return c * b;
}

std::string toFixed(long double z, int n)
{
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(n) << z;
	return oss.str();
}

std::string fill_str(std::string s, int n, std::string c)
{
	int len = s.length();
	for (int i = 0; i < n - len; i++)
	{
		s = c + s;
	}
	return s;
}

std::string dtoa_milo2(long double value, int precision, bool rightpad)
{
	std::ostringstream oss;
	if (std::signbit(value))
	{
		oss << '-';
		value = -value;
	}

	if (value == 0)
	{
		oss << '0';
		if (rightpad && precision > 0)
		{
			oss << '.';
			for (int i = 0; i < precision; ++i)
			{
				oss << '0';
			}
		}
	}
	else
	{
		oss << std::fixed << std::setprecision(precision) << value;
	}

	return oss.str();
}

std::string to_str(long double value, int precision, int length, bool right_align)
{
	std::string str = dtoa_milo2(value, precision, true);

	char fill[64] = "";
	size_t s_len = str.length();
	if (length > s_len)
	{
		size_t len = length - s_len;
		for (size_t i = 0; i < len; i++)
		{
			fill[i] = ' ';
		}
	}
	return right_align ? (std::string(fill) + str) : (str + std::string(fill));
}

std::string m2fm(long double v, int fx, int fs)
{
	std::string gn = "";
	if (v < 0)
	{
		v = -v, gn = "-";
	}
	auto f = int64(v / 60.0);
	long double m = v - f * 60;
	if (!fs)
	{
		return gn + std::to_string(f) + "\'" + toFixed(m, fx) + "\"";
	}
	if (fs == 1)
	{
		return gn + std::to_string(f) + "分" + toFixed(m, fx) + "秒";
	}
	if (fs == 2)
	{
		return gn + std::to_string(f) + "m" + toFixed(m, fx) + "s";
	}
	else
	{
		return "error";
	}
}

// 球面转直角坐标;
Vector3 llr2xyz(Vector3 JW)
{
	Vector3 r;
	long double J = JW[0], W = JW[1], R = JW[2];
	r[0] = R * cos(W) * cos(J);
	r[1] = R * cos(W) * sin(J);
	r[2] = R * sin(W);
	return r;
}

// 直角坐标转球;
Vector3 xyz2llr(Vector3 xyz)
{
	Vector3 r;
	long double x = xyz[0], y = xyz[1], z = xyz[2];
	r[2] = sqrt(x * x + y * y + z * z);
	r[1] = asin(z / r[2]);
	r[0] = rad2mrad(atan2(y, x));
	return r;
}

// 球面坐标旋转;
Vector3 llrConv(Vector3 JW, long double E)
{
	Vector3 r;
	long double J = JW[0], W = JW[1];
	r[0] = atan2(sin(J) * cos(E) - tan(W) * sin(E), cos(J));
	r[1] = asin(cos(E) * sin(W) + sin(E) * cos(W) * sin(J));
	r[2] = JW[2];
	r[0] = rad2mrad(r[0]);
	return r;
}

// 赤道坐标转为地平坐标;
Vector3 CD2DP(Vector3 z, long double L, long double fa, long double gst)
{
	// 转到相对于地平赤道分点的赤道坐标;
	Vector3 a = {z[0] + PI / 2 - gst - L, z[1], z[2]};
	a = llrConv(a, PI / 2 - fa);
	a[0] = rad2mrad(PI / 2 - a[0]);
	return a;
}

// 求角度差;
long double j1_j2(long double J1, long double W1, long double J2, long double W2)
{
	long double dJ = rad2rrad(J1 - J2), dW = W1 - W2;
	if (std::abs(dJ) < 1.0 / 1000 && std::abs(dW) < 1.0 / 1000)
	{
		dJ *= cos((W1 + W2) / 2.0);
		return sqrt(dJ * dJ + dW * dW);
	}
	return acos(sin(W1) * sin(W2) + cos(W1) * cos(W2) * cos(dJ));
}

// 日心球面转地心球面,Z星体球面坐标,A地球球面坐标;
// 本含数是通用的球面坐标中心平移函数,行星计算中将反复使用;
Vector3 h2g(Vector3 z, Vector3 a)
{
	a = llr2xyz(a); // 地球
	z = llr2xyz(z); // 星体
	z[0] -= a[0];
	z[1] -= a[1];
	z[2] -= a[2];
	return xyz2llr(z);
}

// 视差角(不是视差);
long double shiChaJ(long double gst, long double L, long double fa, long double J, long double W)
{
	long double H = gst + L - J; // 天体的时角;
	return rad2mrad(atan2(sin(H), tan(fa) * cos(W) - sin(W) * cos(H)));
}

//=================================deltat T计算=====================================
long double dt_ext(long double y, int jsd)
{
	long double dy = (y - 1820) / 100.0;
	return -20 + jsd * dy * dy;
} // 二次曲线外推

long double dt_calc(long double y)
{ // 计算世界时与原子时之差,传入年

	int dt_at_length = sizeof(dt_at) / sizeof(long double);

	long double y0 = dt_at[dt_at_length - 2]; // 表中最后一年
	long double t0 = dt_at[dt_at_length - 1]; // 表中最后一年的deltatT
	if (y >= y0)
	{
		int jsd = 31; // sjd是y1年之后的加速度估计。瑞士星历表jsd=31,NASA网站jsd=32,skmap的jsd=29
		if (y > y0 + 100)
		{
			return dt_ext(y, jsd);
		}
		long double v = dt_ext(y, jsd);		   // 二次曲线外推
		long double dv = dt_ext(y0, jsd) - t0; // ye年的二次外推与te的差
		return v - dv * (y0 + 100 - y) / 100.0;
	}
	int i;
	const long double *d = dt_at;
	for (i = 0; i < dt_at_length; i += 5)
	{
		if (y < d[i + 5])
		{
			break;
		}
	}

	long double t1 = (y - d[i]) / (d[i + 5] - d[i]) * 10, t2 = t1 * t1, t3 = t2 * t1;
	return d[i + 1] + d[i + 2] * t1 + d[i + 3] * t2 + d[i + 4] * t3;
}

// 传入儒略日(J2000起算),计算TD-UT(单位:日)
long double dt_T(long double t)
{
	return dt_calc(t / 365.2425 + 2000) / 86400.0;
}

//=================================岁差计算=========================================
// IAU1976岁差表

// t是儒略世纪数,sc是岁差量名称,mx是模型
long double prece(long double t, std::string sc, std::string mx)
{
	int i, tn = 1, c = 0, n;
	const long double *p = nullptr;
	if (mx == "IAU1976")
	{
		n = 4, p = preceTab_IAU1976;
	}
	if (mx == "IAU2000")
	{
		n = 6, p = preceTab_IAU2000;
	}
	if (mx == "P03")
	{
		n = 6, p = preceTab_P03;
	}
	auto isc = std::string("fi w  P  Q  E  x  pi II p  th Z  z ").find(sc + ' ') / 3;
	for (i = 0; i < n; i++, tn *= t)
	{
		c += p[isc * n + i] * tn;
	}
	return c / rad;
}

// 返回P03黄赤交角,t是世纪数
long double hcjj(long double t)
{
	long double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
	return (84381.4060 - 46.836769 * t - 0.0001831 * t2 + 0.00200340 * t3 - 5.76e-7 * t4 - 4.34e-8 * t5) / rad;
}

//=================================岁差旋转=========================================
// J2000赤道转Date赤道
Vector3 CDllr_J2D(long double t, Vector3 llr, std::string mx)
{
	long double Z = prece(t, "Z", mx) + llr[0];
	long double z = prece(t, "z", mx);
	long double th = prece(t, "th", mx);
	long double cosW = cos(llr[1]), cosH = cos(th);
	long double sinW = sin(llr[1]), sinH = sin(th);
	long double A = cosW * sin(Z);
	long double B = cosH * cosW * cos(Z) - sinH * sinW;
	long double C = sinH * cosW * cos(Z) + cosH * sinW;
	return Vector3{rad2mrad(atan2(A, B) + z), asin(C), llr[2]};
}
// Date赤道转J2000赤道
Vector3 CDllr_D2J(long double t, Vector3 llr, std::string mx)
{
	long double Z = -prece(t, "z", mx) + llr[0];
	long double z = -prece(t, "Z", mx);
	long double th = -prece(t, "th", mx);
	long double cosW = cos(llr[1]), cosH = cos(th);
	long double sinW = sin(llr[1]), sinH = sin(th);
	long double A = cosW * sin(Z);
	long double B = cosH * cosW * cos(Z) - sinH * sinW;
	long double C = sinH * cosW * cos(Z) + cosH * sinW;
	return Vector3{rad2mrad(atan2(A, B) + z), asin(C), llr[2]};
}

// 黄道球面坐标_J2000转Date分点,t为儒略世纪数
Vector3 HDllr_J2D(long double t, Vector3 llr, std::string mx)
{
	// J2000黄道旋转到Date黄道(球面对球面),也可直接由利用球面旋转函数计算,但交角接近为0时精度很低
	Vector3 r{llr[0], llr[1], llr[2]};
	r[0] += prece(t, "fi", mx);
	r = llrConv(r, prece(t, "w", mx));
	r[0] -= prece(t, "x", mx);
	r = llrConv(r, -prece(t, "E", mx));
	return r;
}

// 黄道球面坐标_Date分点转J2000,t为儒略世纪数
Vector3 HDllr_D2J(long double t, Vector3 llr, std::string mx)
{
	Vector3 r{llr[0], llr[1], llr[2]};
	r = llrConv(r, prece(t, "E", mx));
	r[0] += prece(t, "x", mx);
	r = llrConv(r, -prece(t, "w", mx));
	r[0] -= prece(t, "fi", mx);
	r[0] = rad2mrad(r[0]);
	return r;
}

//=================================章动计算=========================================
//==================================================================================

// 章动计算,t是J2000.0起算的儒略世纪数,zq表示只计算周期大于zq(天)的项
Vector2 nutation(long double t, long double zq)
{
	long double t2 = t * t, t3 = t2 * t, t4 = t3 * t; // t的二、三、四次方
	long double l = 485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3 - 0.00024470 * t4;
	long double l1 = 1287104.79305 + 129596581.0481 * t - 0.5532 * t2 - 0.000136 * t3 - 0.00001149 * t4;
	long double F = 335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3 + 0.00000417 * t4;
	long double D = 1072260.70369 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3 - 0.00003169 * t4;
	long double Om = 450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3 - 0.00005939 * t4;
	long double dL = 0, dE = 0, c;
	const long double *B = nuTab;
	long double q;
	for (int i = 0; i < 77 * 11; i += 11)
	{ // 周期项取和计算
		c = (B[i] * l + B[i + 1] * l1 + B[i + 2] * F + B[i + 3] * D + B[i + 4] * Om) / rad;
		if (zq)
		{ // 只算周期大于zq天的项
			q = 36526 * 2 * PI * rad / (1717915923.2178 * B[i] + 129596581.0481 * B[i + 1] + 1739527262.8478 * B[i + 2] + 1602961601.2090 * B[i + 3] + 6962890.5431 * B[i + 4]);
			if (q < zq)
				continue;
		}
		dL += (B[i + 5] + B[i + 6] * t) * sin(c) + B[i + 7] * cos(c);
		dE += (B[i + 8] + B[i + 9] * t) * cos(c) + B[i + 10] * sin(c);
	}
	return Vector2{dL /= 10000000 * rad, dE /= 10000000 * rad}; // 返回IAU2000B章动值, dL是黄经章动,dE是交角章动
}

// 本函数计算赤经章动及赤纬章动
Vector3 CDnutation(Vector3 z, long double E, long double dL, long double dE)
{
	Vector3 r(z);
	r[0] += (cos(E) + sin(E) * sin(z[0]) * tan(z[1])) * dL - cos(z[0]) * tan(z[1]) * dE; // 赤经章动
	r[1] += sin(E) * cos(z[0]) * dL + sin(z[0]) * dE;									 // 赤纬章动
	r[0] = rad2mrad(r[0]);
	return r;
}

Vector2 nutation2(long double t)
{ // 中精度章动计算,t是世纪数
	long double c, a, t2 = t * t;
	const long double *B = nutB;
	long double dL = 0, dE = 0;
	for (int i = 0; i < sizeof(nutB) / sizeof(long double); i += 5)
	{
		c = B[i] + B[i + 1] * t + B[i + 2] * t2;
		if (i == 0)
		{
			a = -1.742 * t;
		}
		else
		{
			a = 0;
		}
		dL += (B[i + 3] + a) * sin(c);
		dE += B[i + 4] * cos(c);
	}
	return Vector2{dL / 100.0 / rad, dE / 100.0 / rad}; // 黄经章动,交角章动
}

long double nutationLon2(long double t)
{ // 只计算黄经章动
	long double a, t2 = t * t, dL = 0;
	const long double *B = nutB;
	for (int i = 0; i < sizeof(nutB) / sizeof(long double); i += 5)
	{
		if (i == 0)
		{
			a = -1.742 * t;
		}
		else
		{
			a = 0;
		}
		dL += (B[i + 3] + a) * sin(B[i] + B[i + 1] * t + B[i + 2] * t2);
	}
	return dL / 100.0 / rad;
}

//=================================蒙气改正=========================================
//==================================================================================
long double MQC(long double h)
{
	// 大气折射,h是真高度
	return 0.0002967 / tan(h + 0.003138 / (h + 0.08919));
}
long double MQC2(long double ho)
{
	// 大气折射,ho是视高度
	return -0.0002909 / tan(ho + 0.002227 / (ho + 0.07679));
}

//=================================视差改正=========================================
//==================================================================================
Vector3 parallax(Vector3 z, long double H, long double fa, long double high)
{ // 视差修正
	// z赤道坐标,fa地理纬度,H时角,high海拔(千米)
	long double dw = 1;
	if (z[2] < 500)
	{
		dw = cs_AU;
	}
	z[2] *= dw;
	long double r0, x0, y0, z0, f = cs_ba, u = atan(f * tan(fa)), g = z[0] + H;
	r0 = cs_rEar * cos(u) + high * cos(fa);		// 站点与地地心向径的赤道投影长度
	z0 = cs_rEar * sin(u) * f + high * sin(fa); // 站点与地地心向径的轴向投影长度
	x0 = r0 * cos(g);
	y0 = r0 * sin(g);

	Vector3 s = llr2xyz(z);
	s[0] -= x0, s[1] -= y0, s[2] -= z0;
	s = xyz2llr(s);
	s[2] /= dw;
	return s;
}
//=================================星历数据=========================================
//==================================================================================

/********************************
8行星星历数据表,及数据表的计算
********************************/

// xt星体,zn坐标号,t儒略世纪数,n计算项数;
long double XL0_calc(int xt, int zn, long double t, int n)
{
	t /= 10; // 转为儒略千年数
	long double v = 0, tn = 1, c = 0;
	const long double *F = XL0[xt];
	long double n1, n2, N;
	long double n0;
	int pn = zn * 6 + 1;
	long double N0 = F[pn + 1] - F[pn]; // N0序列总数
	for (int i = 0; i < 6; i++, tn *= t)
	{
		n1 = F[pn + i], n2 = F[pn + 1 + i], n0 = n2 - n1;
		if (!n0)
		{
			continue;
		}
		if (n < 0)
		{
			N = n2; // 确定项数
		}
		else
		{
			N = int64(3 * n * n0 / N0 + 0.5) + n1;
			if (i)
			{
				N += 3;
			}
			if (N > n2)
			{
				N = n2;
			}
		}
		c = 0;
		for (int64_t j = n1; j < N; j += 3)
		{
			c += F[j] * cos(F[j + 1] + t * F[j + 2]);
		}
		v += c * tn;
	}
	v /= F[0];
	if (xt == 0)
	{										 // 地球
		long double t2 = t * t, t3 = t2 * t; // 千年数的各次方
		if (zn == 0)
		{
			v += (-0.0728 - 2.7702 * t - 1.1019 * t2 - 0.0996 * t3) / rad;
		}
		if (zn == 1)
		{
			v += (+0.0000 + 0.0004 * t + 0.0004 * t2 - 0.0026 * t3) / rad;
		}
		if (zn == 2)
		{
			v += (-0.0020 + 0.0044 * t + 0.0213 * t2 - 0.0250 * t3) / 1000000;
		}
	}
	else
	{ // 其它行星
		long double dv = XL0_xzb[(xt - 1) * 3 + zn];
		if (zn == 0)
		{
			v += -3 * t / rad;
		}
		if (zn == 2)
		{
			v += dv / 1000000;
		}
		else
		{
			v += dv / rad;
		}
	}
	return v;
}

// 返回冥王星J2000直角坐标(未实现)
Vector3 pluto_coord(long double t)
{
	long double c0 = PI / 180.0 / 100000.0;
	long double x = -1 + 2 * (t * 36525 + 1825394.5) / 2185000;
	long double T = t / 100000000;
	Vector3 r;

	for (int i = 0; i < 9; i++)
	{
		const long double *ob = XL0Pluto[i];
		int N = sizeof(XL0Pluto[i]) / sizeof(long double);
		long double v = 0;
		for (int j = 0; j < N; j += 3)
		{
			v += ob[j] * sin(ob[j + 1] * T + ob[j + 2] * c0);
		}
		if (i % 3 == 1)
		{
			v *= x;
		}
		if (i % 3 == 2)
		{
			v *= x * x;
		}
		r[int2(i / 3.0)] += v / 100000000.0;
	}
	r[0] += 9.922274 + 0.154154 * x;
	r[1] += 10.016090 + 0.064073 * x;
	r[2] += -3.947474 - 0.042746 * x;
	return r;
}

// xt星体,T儒略世纪数,TD
Vector3 p_coord(int xt, long double t, int n1, int n2, int n3)
{
	Vector3 z;
	if (xt < 8)
	{
		z[0] = XL0_calc(xt, 0, t, n1);
		z[1] = XL0_calc(xt, 1, t, n2);
		z[2] = XL0_calc(xt, 2, t, n3);
	}
	if (xt == 8)
	{ // 冥王星
		z = pluto_coord(t);
		z = xyz2llr(z);
		z = HDllr_J2D(t, z, "P03");
	}
	if (xt == 9)
	{ // 太阳
		z[0] = 0, z[1] = 0, z[2] = 0;
	}
	return z;
}

// 返回地球坐标,t为世纪数
Vector3 e_coord(long double t, int n1, int n2, int n3)
{
	Vector3 re;
	re[0] = XL0_calc(0, 0, t, n1);
	re[1] = XL0_calc(0, 1, t, n2);
	re[2] = XL0_calc(0, 2, t, n3);
	return re;
}

//=================================月亮星历--=======================================
//==================================================================================
// 计算月亮
long double XL1_calc(int zn, long double t, int n)
{
	const long double *const *ob = XL1[zn];
	long double F, N, v = 0, tn = 1, c = 0;
	long double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t, tx = t - 10;
	if (zn == 0)
	{
		v += (3.81034409 + 8399.684730072 * t - 3.319e-05 * t2 + 3.11e-08 * t3 - 2.033e-10 * t4) * rad; // 月球平黄经(弧度)
		v += 5028.792262 * t + 1.1124406 * t2 + 0.00007699 * t3 - 0.000023479 * t4 - 0.0000000178 * t5; // 岁差(角秒)
		if (tx > 0)
		{
			v += -0.866 + 1.43 * tx + 0.054 * tx * tx; // 对公元3000年至公元5000年的拟合,最大误差小于10角秒
		}
	}
	t2 /= 1e4, t3 /= 1e8, t4 /= 1e8;
	n *= 6;
	if (n < 0)
	{
		n = 3 /*ob[0].length*/;
	}
	for (int i = 0; i < 3 /*ob.length*/; i++, tn *= t)
	{
		const long double *F = ob[i];
		int F_len = 0;
		if (zn == 0 && i == 0)
		{
			F_len = sizeof(XL1_0_0) / sizeof(long double);
		}
		else if (zn == 0 && i == 1)
		{
			F_len = sizeof(XL1_0_1) / sizeof(long double);
		}
		else if (zn == 0 && i == 2)
		{
			F_len = sizeof(XL1_0_2) / sizeof(long double);
		}
		else if (zn == 1 && i == 0)
		{
			F_len = sizeof(XL1_1_0) / sizeof(long double);
		}
		else if (zn == 1 && i == 1)
		{
			F_len = sizeof(XL1_1_1) / sizeof(long double);
		}
		else if (zn == 1 && i == 2)
		{
			F_len = sizeof(XL1_1_2) / sizeof(long double);
		}
		else if (zn == 2 && i == 0)
		{
			F_len = sizeof(XL1_2_0) / sizeof(long double);
		}
		else if (zn == 2 && i == 1)
		{
			F_len = sizeof(XL1_2_1) / sizeof(long double);
		}
		else if (zn == 2 && i == 2)
		{
			F_len = sizeof(XL1_2_2) / sizeof(long double);
		}

		N = int64(n * F_len / 3.0 /*ob[0].length*/ + 0.5);
		if (i)
		{
			N += 6;
		}
		if (N >= F_len)
		{
			N = F_len;
		}
		c = 0;
		for (int j = 0; j < N; j += 6)
		{
			c += F[j] * cos(F[j + 1] + t * F[j + 2] + t2 * F[j + 3] + t3 * F[j + 4] + t4 * F[j + 5]);
		}
		v += c * tn;
	}
	if (zn != 2)
	{
		v /= rad;
	}
	return v;
};

// 返回月球坐标,n1,n2,n3为各坐标所取的项数
Vector3 m_coord(long double t, int n1, int n2, int n3)
{
	Vector3 re;
	re[0] = XL1_calc(0, t, n1);
	re[1] = XL1_calc(1, t, n2);
	re[2] = XL1_calc(2, t, n3);
	return re;
}

//=============================一些天文基本问题=====================================
//==================================================================================
// 返回朔日的编号,jd应在朔日附近，允许误差数天
int64_t suoN(long double jd)
{
	return int64((jd + 8) / 29.5306);
}

long double gxc_sunLon(long double t)
{																	  // 太阳光行差,t是世纪数
	long double v = -0.043126 + 628.301955 * t - 0.000002732 * t * t; // 平近点角
	long double e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;
	return (-20.49552 * (1 + e * cos(v))) / rad; // 黄经光行差
}

// 黄纬光行差
long double gxc_sunLat(long double t)
{
	return 0;
}

// 月球经度光行差,误差0.07"
long double gxc_moonLon(long double t)
{
	return -3.4E-6;
}

// 月球纬度光行差,误差0.006"
long double gxc_moonLat(long double t)
{
	return 0.063 * sin(0.057 + 8433.4662 * t + 0.000064 * t * t) / rad;
}

// 传入T是2000年首起算的日数(UT),dt是deltatT(日),精度要求不高时dt可取值为0																					 //传入T是2000年首起算的日数(UT),dt是deltatT(日),精度要求不高时dt可取值为0
long double pGST(long double T, long double dt)
{
	// 返回格林尼治平恒星时(不含赤经章动及非多项式部分),即格林尼治子午圈的平春风点起算的赤经
	long double t = (T + dt) / 36525, t2 = t * t, t3 = t2 * t, t4 = t3 * t;
	return pi2 * (0.7790572732640 + 1.00273781191135448 * T) // T是UT,下一行的t是力学时(世纪数)
		   + (0.014506 + 4612.15739966 * t + 1.39667721 * t2 - 0.00009344 * t3 + 0.00001882 * t4) / rad;
}

long double pGST2(long double jd)
{ // 传入力学时J2000起算日数，返回平恒星时
	long double dt = dt_T(jd);
	return pGST(jd - dt, dt);
}

// 太阳升降计算。jd儒略日(须接近L当地平午UT)，L地理经度，fa地理纬度，sj=-1升,sj=1降
long double sunShengJ(long double jd, long double L, long double fa, long double sj)
{
	jd = floor(jd + 0.5) - L / pi2;
	for (int i = 0; i < 2; i++)
	{
		long double T = jd / 36525, E = (84381.4060 - 46.836769 * T) / rad;	   // 黄赤交角
		long double t = T + (32 * (T + 1.8) * (T + 1.8) - 20) / 86400 / 36525; // 儒略世纪年数,力学时
		long double J = (48950621.66 + 6283319653.318 * t + 53 * t * t - 994 + 334166 * cos(4.669257 + 628.307585 * t) + 3489 * cos(4.6261 + 1256.61517 * t) + 2060.6 * cos(2.67823 + 628.307585 * t) * t) / 10000000;
		long double sinJ = sin(J), cosJ = cos(J);																						  // 太阳黄经以及它的正余弦值
		long double gst = (0.7790572732640 + 1.00273781191135448 * jd) * pi2 + (0.014506 + 4612.15739966 * T + 1.39667721 * T * T) / rad; // 恒星时(子午圈位置)
		long double A = atan2(sinJ * cos(E), cosJ);																						  // 太阳赤经
		long double D = asin(sin(E) * sinJ);																							  // 太阳赤纬
		long double cosH0 = (sin(-50 * 60 / rad) - sin(fa) * sin(D)) / (cos(fa) * cos(D));
		if (fabs(cosH0) >= 1)
		{
			return 0; // 太阳在地平线上的cos(时角)计算
		}
		jd += rad2rrad(sj * acos(cosH0) - (gst + L - A)) / 6.28; //(升降时角-太阳时角)/太阳速度
	}
	return jd; // 反回格林尼治UT
}

// 时差计算(高精度),t力学时儒略世纪数
long double pty_zty(long double t)
{
	long double t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
	long double L = (1753470142 + 628331965331.8 * t + 5296.74 * t2 + 0.432 * t3 - 0.1124 * t4 - 0.00009 * t5) / 1000000000 + PI - 20.5 / rad;

	long double E, dE, dL, f;
	Vector3 z;
	dL = -17.2 * sin(2.1824 - 33.75705 * t) / rad; // 黄经章
	dE = 9.2 * cos(2.1824 - 33.75705 * t) / rad;   // 交角章
	E = hcjj(t) + dE;							   // 真黄赤交角

	// 地球坐标
	z[0] = XL0_calc(0, 0, t, 50) + PI + gxc_sunLon(t) + dL;
	z[1] = -(2796 * cos(3.1987 + 8433.46616 * t) + 1016 * cos(5.4225 + 550.75532 * t) + 804 * cos(3.88 + 522.3694 * t)) / 1000000000;

	z = llrConv(z, E); // z太阳地心赤道坐标
	z[0] -= dL * cos(E);

	L = rad2rrad(L - z[0]);
	return L / pi2; // 单位是周(天)
}

long double pty_zty2(long double t)
{ // 时差计算(低精度),误差约在1秒以内,t力学时儒略世纪数
	long double L = (1753470142 + 628331965331.8 * t + 5296.74 * t * t) / 1000000000 + PI;
	Vector3 z;
	long double E = (84381.4088 - 46.836051 * t) / rad;
	z[0] = XL0_calc(0, 0, t, 5) + PI, z[1] = 0; // 地球坐标
	z = llrConv(z, E);							// z太阳地心赤道坐标
	L = rad2rrad(L - z[0]);
	return L / pi2; // 单位是周(天)
}

// 物件XL : 日月黄道平分点坐标、视坐标、速度、已知经度反求时间等方面的计算
namespace XL
{
	//=====================
	// 星历函数(日月球面坐标计算)

	// 地球经度计算,返回Date分点黄经,传入世纪数、取项数
	long double E_Lon(long double t, int n)
	{
		return XL0_calc(0, 0, t, n);
	}

	// 月球经度计算,返回Date分点黄经,传入世纪数,n是项数比例
	long double M_Lon(long double t, int n)
	{
		return XL1_calc(0, t, n);
	}

	//=========================
	long double E_v(long double t)
	{ // 地球速度,t是世纪数,误差小于万分3
		long double f = 628.307585 * t;
		return 628.332 + 21 * sin(1.527 + f) + 0.44 * sin(1.48 + f * 2) + 0.129 * sin(5.82 + f) * t + 0.00055 * sin(4.21 + f) * t * t;
	}

	// 月球速度计算,传入世经数
	long double M_v(long double t)
	{
		long double v = 8399.71 - 914 * sin(0.7848 + 8328.691425 * t + 0.0001523 * t * t); // 误差小于5%
		v -= 179 * sin(2.543 + 15542.7543 * t)											   // 误差小于0.3%
			 + 160 * sin(0.1874 + 7214.0629 * t) + 62 * sin(3.14 + 16657.3828 * t) + 34 * sin(4.827 + 16866.9323 * t) + 22 * sin(4.9 + 23871.4457 * t) + 12 * sin(2.59 + 14914.4523 * t) + 7 * sin(0.23 + 6585.7609 * t) + 5 * sin(0.9 + 25195.624 * t) + 5 * sin(2.32 - 7700.3895 * t) + 5 * sin(3.88 + 8956.9934 * t) + 5 * sin(0.49 + 7771.3771 * t);
		return v;
	}

	//=========================
	// 月日视黄经的差值
	long double MS_aLon(long double t, long double Mn, long double Sn)
	{
		return M_Lon(t, Mn) + gxc_moonLon(t) - (E_Lon(t, Sn) + gxc_sunLon(t) + PI);
	}

	// 太阳视黄经
	long double S_aLon(long double t, long double n)
	{
		return E_Lon(t, n) + nutationLon2(t) + gxc_sunLon(t) + PI; // 注意，这里的章动计算很耗时
	}

	//=========================
	// 已知地球真黄经求时间
	long double E_Lon_t(long double W)
	{
		long double t, v = 628.3319653318;
		t = (W - 1.75347) / v;
		v = E_v(t); // v的精度0.03%，详见原文
		t += (W - E_Lon(t, 10)) / v;
		v = E_v(t); // 再算一次v有助于提高精度,不算也可以
		t += (W - E_Lon(t, -1)) / v;
		return t;
	}

	// 已知真月球黄经求时间
	long double M_Lon_t(long double W)
	{
		long double t, v = 8399.70911033384;
		t = (W - 3.81034) / v;
		t += (W - M_Lon(t, 3)) / v;
		v = M_v(t); // v的精度0.5%，详见原文
		t += (W - M_Lon(t, 20)) / v;
		t += (W - M_Lon(t, -1)) / v;
		return t;
	};

	// 已知月日视黄经差求时间
	long double MS_aLon_t(long double W)
	{
		long double t, v = 7771.37714500204;
		t = (W + 1.08472) / v;
		t += (W - MS_aLon(t, 3, 3)) / v;
		v = M_v(t) - E_v(t); // v的精度0.5%，详见原文
		t += (W - MS_aLon(t, 20, 10)) / v;
		t += (W - MS_aLon(t, -1, 60)) / v;
		return t;
	}

	// 已知太阳视黄经反求时间
	long double S_aLon_t(long double W)
	{
		long double t, v = 628.3319653318;
		t = (W - 1.75347 - PI) / v;
		v = E_v(t); // v的精度0.03%，详见原文
		t += (W - S_aLon(t, 10)) / v;
		v = E_v(t); // 再算一次v有助于提高精度,不算也可以
		t += (W - S_aLon(t, -1)) / v;
		return t;
	};
	/****
	MS_aLon_t1:function(W){ //已知月日视黄经差求时间,高速低精度,误差不超过40秒
	  var t,v = 7771.37714500204;
	  t  = ( W + 1.08472               )/v;
	  t += ( W - this.MS_aLon(t, 3, 3) )/v;  v=this.M_v(t)-this.E_v(t);  //v的精度0.5%，详见原文
	  t += ( W - this.MS_aLon(t,50,20) )/v;
	  return t;
	},
	S_aLon_t1:function(W){ //已知太阳视黄经反求时间,高速低精度,最大误差不超过50秒,平均误差15秒
	  var t,v= 628.3319653318;
	  t  = ( W - 1.75347-Math.PI   )/v; v = 628.332 + 21*Math.sin( 1.527+628.307585*t );
	  t += ( W - this.S_aLon(t,3) )/v;
	  t += ( W - this.S_aLon(t,40))/v;
	  return t;
	},
	****/
	// 已知月日视黄经差求时间,高速低精度,误差不超过600秒(只验算了几千年)
	long double MS_aLon_t2(long double W)
	{
		long double t, v = 7771.37714500204;
		t = (W + 1.08472) / v;
		long double L, t2 = t * t;
		t -= (-0.00003309 * t2 + 0.10976 * cos(0.784758 + 8328.6914246 * t + 0.000152292 * t2) + 0.02224 * cos(0.18740 + 7214.0628654 * t - 0.00021848 * t2) - 0.03342 * cos(4.669257 + 628.307585 * t)) / v;
		L = M_Lon(t, 20) - (4.8950632 + 628.3319653318 * t + 0.000005297 * t * t + 0.0334166 * cos(4.669257 + 628.307585 * t) + 0.0002061 * cos(2.67823 + 628.307585 * t) * t + 0.000349 * cos(4.6261 + 1256.61517 * t) - 20.5 / rad);
		v = 7771.38 - 914 * sin(0.7848 + 8328.691425 * t + 0.0001523 * t * t) - 179 * sin(2.543 + 15542.7543 * t) - 160 * sin(0.1874 + 7214.0629 * t);
		t += (W - L) / v;
		return t;
	}

	long double S_aLon_t2(long double W)
	{ // 已知太阳视黄经反求时间,高速低精度,最大误差不超过600秒
		long double t, L, v = 628.3319653318;
		t = (W - 1.75347 - PI) / v;
		t -= (0.000005297 * t * t + 0.0334166 * cos(4.669257 + 628.307585 * t) + 0.0002061 * cos(2.67823 + 628.307585 * t) * t) / v;
		t += (W - E_Lon(t, 8) - PI + (20.5 + 17.2 * sin(2.1824 - 33.75705 * t)) / rad) / v;
		return t;
	}

	long double moonIll(long double t)
	{ // 月亮被照亮部分的比例
		long double t2 = t * t, t3 = t2 * t, t4 = t3 * t;
		long double D, M, m, a, dm = PI / 180;
		D = (297.8502042 + 445267.1115168 * t - 0.0016300 * t2 + t3 / 545868 - t4 / 113065000) * dm; // 日月平距角
		M = (357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000) * dm;				 // 太阳平近点
		m = (134.9634114 + 477198.8676313 * t + 0.0089970 * t2 + t3 / 69699 - t4 / 14712000) * dm;	 // 月亮平近点
		a = PI - D + (-6.289 * sin(m) + 2.100 * sin(M) - 1.274 * sin(D * 2 - m) - 0.658 * sin(D * 2) - 0.214 * sin(m * 2) - 0.110 * sin(D)) * dm;
		return (1 + cos(a)) / 2;
	}

	// 转入地平纬度及地月质心距离,返回站心视半径(角秒)
	long double moonRad(long double r, long double h)
	{
		return cs_sMoon / r * (1 + sin(h) * cs_rEar / r);
	}

	// 求月亮近点时间和距离,t为儒略世纪数力学时
	Vector2 moonMinR(long double t, bool min)
	{
		long double a = 27.55454988 / 36525, b;
		if (min)
		{
			b = -10.3302 / 36525;
		}
		else
		{
			b = 3.4471 / 36525;
		}
		t = b + a * int64((t - b) / a + 0.5); // 平近(远)点时间
		long double r1, r2, r3, dt;
		// 初算二次
		dt = 2.0 / 36525;
		r1 = XL1_calc(2, t - dt, 10);
		r2 = XL1_calc(2, t, 10);
		r3 = XL1_calc(2, t + dt, 10);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0;
		dt = 0.5 / 36525;
		r1 = XL1_calc(2, t - dt, 20);
		r2 = XL1_calc(2, t, 20);
		r3 = XL1_calc(2, t + dt, 20);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0;
		// 精算
		dt = 1200.0 / 86400.0 / 36525.0;
		r1 = XL1_calc(2, t - dt, -1);
		r2 = XL1_calc(2, t, -1);
		r3 = XL1_calc(2, t + dt, -1);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0;
		r2 += (r1 - r3) / (r1 + r3 - 2 * r2) * (r3 - r1) / 8.0;
		Vector2 re{t, r2};
		return re;
	}

	Vector2 moonNode(long double t, bool asc)
	{ // 月亮升交点
		long double a = 27.21222082 / 36525, b;
		if (asc)
		{
			b = 21.0 / 36525;
		}
		else
		{
			b = 35.0 / 36525;
		}
		t = b + a * int64((t - b) / a + 0.5); // 平升(降)交点时间
		long double w, v, w2, dt;
		dt = 0.5 / 36525;
		w = XL1_calc(1, t, 10);
		w2 = XL1_calc(1, t + dt, 10);
		v = (w2 - w) / dt;
		t -= w / v;
		dt = 0.05 / 36525;
		w = XL1_calc(1, t, 40);
		w2 = XL1_calc(1, t + dt, 40);
		v = (w2 - w) / dt;
		t -= w / v;
		w = XL1_calc(1, t, -1);
		t -= w / v;
		Vector2 re;
		re[0] = t;
		re[1] = XL1_calc(0, t, -1);
		return re;
	}

	// 地球近远点
	Vector2 earthMinR(long double t, bool min)
	{
		long double a = 365.25963586 / 36525, b;
		if (min)
		{
			b = 1.7 / 36525;
		}
		else
		{
			b = 184.5 / 36525;
		}
		t = b + a * int64((t - b) / a + 0.5); // 平近(远)点时间
		long double r1, r2, r3, dt;
		// 初算二次
		dt = 3.0 / 36525;
		r1 = XL0_calc(0, 2, t - dt, 10);
		r2 = XL0_calc(0, 2, t, 10);
		r3 = XL0_calc(0, 2, t + dt, 10);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0; // 误差几个小时
		dt = 0.2 / 36525;
		r1 = XL0_calc(0, 2, t - dt, 80);
		r2 = XL0_calc(0, 2, t, 80);
		r3 = XL0_calc(0, 2, t + dt, 80);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2; // 误差几分钟
		// 精算
		dt = 0.01 / 36525;
		r1 = XL0_calc(0, 2, t - dt, -1);
		r2 = XL0_calc(0, 2, t, -1);
		r3 = XL0_calc(0, 2, t + dt, -1);
		t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0; // 误差小于秒
		r2 += (r1 - r3) / (r1 + r3 - 2 * r2) * (r3 - r1) / 8.0;
		Vector2 re{t, r2};
		return re;
	}

	Vector2 daJu(int xt, long double t, bool dx)
	{ // 大距计算超底速算法, dx=1东大距,t儒略世纪TD
		long double a = 0, b, c[5] = {};
		if (xt == 1)
		{
			a = 115.8774777586 / 36525.0;
			long double arr[] = {2, 0.2, 0.01, 46, 87};
			memcpy(c, arr, sizeof(arr));
		} // 水星
		if (xt == 2)
		{
			a = 583.9213708245 / 36525.0;
			long double arr[] = {4, 0.2, 0.01, 382, 521};
			memcpy(c, arr, sizeof(arr));
		} // 金星
		if (dx)
		{
			b = c[3] / 36525.0;
		}
		else
		{
			b = c[4] / 36525.0;
		}
		t = b + a * int2((t - b) / a + 0.5); // 大距平时间
		long double dt, r1, r2, r3;
		int i;
		for (i = 0; i < 3; i++)
		{
			dt = c[i] / 36525.0;
			r1 = xingJJ(xt, t - dt, i);
			r2 = xingJJ(xt, t, i);
			r3 = xingJJ(xt, t + dt, i);
			t += (r1 - r3) / (r1 + r3 - 2 * r2) * dt / 2.0;
		}
		r2 += (r1 - r3) / (r1 + r3 - 2 * r2) * (r3 - r1) / 8.0;
		Vector2 re{t, r2};
		return re;
	}

	// 合冲日计算(视黄经合冲)
	Vector4 xingSP(int xt, long double t, int n, long double w0, long double ts, long double tp)
	{ // 行星太阳视黄经差与w0的差
		Vector3 a, p, s;
		a = p_coord(0, t - tp, n, n, n);  // 地球
		p = p_coord(xt, t - tp, n, n, n); // 行星
		s = p_coord(0, t - ts, n, n, n);
		s[0] += PI;
		s[1] = -s[1]; // 太阳
		p = h2g(p, a);
		Vector4 re{rad2rrad(p[0] - s[0] - w0), p[1] - s[1], s[2] * cs_Agx, p[2] * cs_Agx}; // 赤经差及光行时
		return re;
	}

	Vector2 xingHR(int xt, long double t, bool f)
	{ // xt星体号,t儒略世纪TD,f=1求冲(或下合)否则求合(或下合)
		Vector4 a, b;
		long double i, v, dt = 2e-5;
		long double w0 = PI, w1 = 0; // 合(或上合)时,日心黄经差为180，地心黄经差为0
		if (f)
		{			// 求冲(或下合)
			w0 = 0; // 日心黄经差
			if (xt > 2)
			{
				w1 = PI; // 地心黄经差(冲)
			}
		}
		v = pi2 / cs_xxHH[xt - 1] * 36525;
		if (xt > 2)
		{
			v = -v; // 行星相对地球的黄经平速度
		}
		for (i = 0; i < 6; i++)
		{
			t -= rad2rrad(XL0_calc(xt, 0, t, 8) - XL0_calc(0, 0, t, 8) - w0) / v; // 水星的平速度与真速度相差较多,所以多算几次
		}
		// 严格计算
		a = xingSP(xt, t, 8, w1, 0, 0);
		b = xingSP(xt, t + dt, 8, w1, 0, 0);
		v = (b[0] - a[0]) / dt;
		a = xingSP(xt, t, 40, w1, a[2], a[3]);
		t -= a[0] / v;
		a = xingSP(xt, t, -1, w1, a[2], a[3]);
		t -= a[0] / v;
		Vector2 re{t, a[1]};
		return re;
	}

	// 合月计算
	Vector4 xingMP(int xt, long double t, int n, long double E, Vector4 g)
	{ // 月亮行星视赤经差
		Vector3 a, p, m;
		a = p_coord(0, t - g[1], n, n, n);	// 地球
		p = p_coord(xt, t - g[1], n, n, n); // 行星
		m = m_coord(t - g[0], n, n, n);		// 月亮
		p = h2g(p, a);
		m[0] += g[2];
		p[0] += g[2];
		m = llrConv(m, E + g[3]);
		p = llrConv(p, E + g[3]);
		Vector4 re{rad2rrad(m[0] - p[0]), m[1] - p[1], m[2] / cs_GS / 86400 / 36525.0, p[2] / cs_GS / 86400 / 36525.0 * cs_AU}; // 赤经差及光行时
		return re;
	}

	Vector2 xingHY(int xt, long double t)
	{ // 行星合月(视赤经),t儒略世纪TD
		long double i, v, E;
		Vector4 d, d2, g = {0, 0, 0, 0};
		for (i = 0; i < 3; i++)
		{
			d = xingMP(xt, t, 8, 0.4091, g);
			t -= d[0] / 8192.0;
		}
		E = hcjj(t);
		Vector2 zd = nutation2(t);
		g = {d[2], d[3], zd[0], zd[1]}; // 光行时,章动

		d = xingMP(xt, t, 8, E, g);
		d2 = xingMP(xt, t + 1e-6, 8, E, g);
		v = (d2[0] - d[0]) / 1e-6; // 速度

		d = xingMP(xt, t, 30, E, g);
		t -= d[0] / v;
		d = xingMP(xt, t, -1, E, g);
		t -= d[0] / v;
		Vector2 re{t, d[1]};
		return re;
	}

	Vector3 xingLiu0(int xt, long double t, int n, long double gxs)
	{ // 行星的视坐标
		Vector3 a, z;
		Vector2 zd;
		long double E = hcjj(t);
		a = p_coord(0, t - gxs, n, n, n);  // 地球
		z = p_coord(xt, t - gxs, n, n, n); // 行星
		z = h2g(z, a);					   // 转到地心
		if (gxs)
		{					   // 如果计算了光行时，那么也计算章动
			zd = nutation2(t); // 章动计算
			z[0] += zd[0];
			E += zd[1];
		}
		z = llrConv(z, E);
		return z;
	}

	long double xingLiu(int xt, long double t, bool sn)
	{ // 留,sn=1顺留
		Vector3 y1, y2, y3;
		long double g;
		int i, n;
		// 先求冲(下合)
		long double hh = cs_xxHH[xt - 1] / 36525.0; // 会合周期
		long double v = pi2 / hh;
		if (xt > 2)
		{
			v = -v; // 行星相对地球的黄经平速度
		}
		for (i = 0; i < 6; i++)
		{
			t -= rad2rrad(XL0_calc(xt, 0, t, 8) - XL0_calc(0, 0, t, 8)) / v; // 水星的平速度与真速度相差较多,所以多算几次
		}

		long double tt[] = {5 / 36525.0, 1 / 36525.0, 0.5 / 36525.0, 2e-6, 2e-6};
		long double dt;
		long double tc0[] = {17.4, 28, 52, 82, 86, 88, 89, 90};
		long double tc = tc0[xt - 1] / 36525.0;

		if (sn)
		{
			if (xt > 2)
			{
				t -= tc;
			}
			else
			{
				t += tc;
			}
		} // 顺留
		else
		{
			if (xt > 2)
			{
				t += tc;
			}
			else
			{
				t -= tc;
			}
		} // 逆留
		for (i = 0; i < 4; i++)
		{
			dt = tt[i], n = 8, g = 0;
			if (i >= 3)
			{
				g = y2[2] * cs_Agx;
				n = -1;
			}
			y1 = xingLiu0(xt, t - dt, n, g);
			y2 = xingLiu0(xt, t, n, g);
			y3 = xingLiu0(xt, t + dt, n, g);
			t += (y1[0] - y3[0]) / (y1[0] + y3[0] - 2 * y2[0]) * dt / 2.0;
		}
		return t;
	}

	long double xingJJ(int xt, long double t, int jing)
	{
		Vector3 a, z;
		a = p_coord(0, t, 10, 10, 10);	// 地球
		z = p_coord(xt, t, 10, 10, 10); // 行星
		z = h2g(z, a);					// 转到地心
		// if (jing == 0) ;			//低精度
		if (jing == 1)
		{									// 中精度
			a = p_coord(0, t, 60, 60, 60);	// 地球
			z = p_coord(xt, t, 60, 60, 60); // 行星
			z = h2g(z, a);					// 转到地心
		}
		if (jing >= 2)
		{													// 高精度(补光行时)
			a = p_coord(0, t - a[2] * cs_Agx, -1, -1, -1);	// 地球
			z = p_coord(xt, t - z[2] * cs_Agx, -1, -1, -1); // 行星
			z = h2g(z, a);									// 转到地心
		}
		a[0] += PI, a[1] = -a[1]; // 太阳
		return j1_j2(z[0], z[1], a[0], a[1]);
	}

	std::string xingX(int xt, long double jd, double L, double fa)
	{
		// 行星计算,jd力学时
		// 基本参数计算

		long double T = jd / 36525.0;
		Vector2 zd = nutation2(T);
		long double dL = zd[0], dE = zd[1]; // 章动
		long double E = hcjj(T) + dE;		// 真黄赤交角

		long double gstPing = pGST2(jd);		 // 平恒星时
		long double gst = gstPing + dL * cos(E); // 真恒星时(不考虑非多项式部分)

		Vector3 z, a, z2, a2;
		std::string s;
		long double ra, rb, rc;
		int rfn = 8;
		if (xt == 10)
		{ // 月亮
			rfn = 2;
			// 求光行时并精确求出地月距
			a = e_coord(T, 15, 15, 15); // 地球
			z = m_coord(T, 1, 1, -1);
			ra = z[2];				  // 月亮
			T -= ra * cs_Agx / cs_AU; // 光行时计算

			// 求视坐标
			a2 = e_coord(T, 15, 15, 15); // 地球
			z = m_coord(T, -1, -1, -1);
			rc = z[2]; // 月亮
			// 求光行距
			a2 = h2g(a, a2);
			a2[2] *= cs_AU;
			z2 = h2g(z, a2);
			rb = z2[2];

			// 地心黄道及地心赤道
			z[0] = rad2mrad(z[0] + dL); // 补章动
			s += "视黄经 " + rad2str(z[0], false) + " 视黄纬 " + rad2str(z[1], false) + " 地心距 " + toFixed(ra, rfn) + "\n";

			z = llrConv(z, E); // 转到赤道坐标
			s += "视赤经 " + rad2str(z[0], true) + " 视赤纬 " + rad2str(z[1], false) + " 光行距 " + toFixed(rb, rfn) + "\n";
		}
		else if (xt < 10 && xt >= 0)
		{
			a = p_coord(0, T, -1, -1, -1);	// 地球
			z = p_coord(xt, T, -1, -1, -1); // 行星
			z[0] = rad2mrad(z[0]);
			s += "黄经一 " + rad2str(z[0], false) + " 黄纬一 " + rad2str(z[1], false) + " 向径一 " + toFixed(z[2], rfn) + "\n";

			// 地心黄道
			z = h2g(z, a);
			ra = z[2];		  // ra地心距
			T -= ra * cs_Agx; // 光行时

			// 重算坐标
			a2 = p_coord(0, T, -1, -1, -1);	 // 地球
			z2 = p_coord(xt, T, -1, -1, -1); // 行星

			z = h2g(z2, a);
			rb = z[2]; // rb光行距(在惯性系中看)
			z = h2g(z2, a2);
			rc = z[2];					// rc视距
			z[0] = rad2mrad(z[0] + dL); // 补章动
			s += "视黄经 " + rad2str(z[0], false) + " 视黄纬 " + rad2str(z[1], false) + " 地心距 " + toFixed(ra, rfn) + "\n";

			z = llrConv(z, E); // 转到赤道坐标
			s += "视赤经 " + rad2str(z[0], true) + " 视赤纬 " + rad2str(z[1], false) + " 光行距 " + toFixed(rb, rfn) + "\n";
		}
		long double sj = rad2rrad(gst + L - z[0]); // 得到天体时角
		z = parallax(z, sj, fa, 0);				   // 视差修正
		s += "站赤经 " + rad2str(z[0], true) + " 站赤纬 " + rad2str(z[1], false) + " 视距离 " + toFixed(rc, rfn) + "\n";

		z[0] += M_PI / 2 - gst - L;	   // 修正了视差的赤道坐标
		z = llrConv(z, M_PI / 2 - fa); // 转到时角坐标转到地平坐标
		z[0] = rad2mrad(M_PI / 2 - z[0]);

		if (z[1] > 0)
		{
			z[1] += MQC(z[1]); // 大气折射修正
		}
		s += "方位角 " + rad2str(z[0], false) + " 高度角 " + rad2str(z[1], false) + "\n";
		s += "恒星时 " + rad2str(rad2mrad(gstPing), true) + "(平) " + rad2str(rad2mrad(gst), true) + "(真)\n";

		return s;
	}
};

//========日月食计算使用的一些函数=============
COORDP lineEar(Vector3 P, Vector3 Q, long double gst)
{
	// 在分点坐标中求空间两点连线与地球的交点(靠近点P的交点),返回地标
	Vector3 p = llr2xyz(P), q = llr2xyz(Q);
	COORDP r = lineEll(p[0], p[1], p[2], q[0], q[1], q[2], cs_ba, cs_rEar);
	if (r.D < 0)
	{
		r.J = r.W = 100.0;
		return r;
	} // 反回100表示无解
	r.W = atan(r.z / cs_ba2 / sqrt(r.x * r.x + r.y * r.y));
	r.J = rad2rrad(atan2(r.y, r.x) - gst);
	return r;
}

COORDP lineEll(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2, long double e, long double r)
{
	// 求空间两点连线与地球的交点(靠近点x1的交点)
	long double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1, e2 = e * e, A, B, C, D, R, t;
	COORDP p = {};
	A = dx * dx + dy * dy + dz * dz / e2;
	B = x1 * dx + y1 * dy + z1 * dz / e2;
	C = x1 * x1 + y1 * y1 + z1 * z1 / e2 - r * r;
	p.D = B * B - A * C;
	if (p.D < 0)
	{
		return p; // 判别式小于0无解
	}
	D = sqrt(p.D);
	if (B < 0)
	{
		D = -D; // 只求靠近x1的交点
	}
	t = (-B + D) / A;
	p.x = x1 + dx * t, p.y = y1 + dy * t, p.z = z1 + dz * t;
	R = sqrt(dx * dx + dy * dy + dz * dz);
	p.R1 = R * fabs(t), p.R2 = R * fabs(t - 1); // R1,R2分别为x1,x2到交点的距离
	return p;
}

COORDP lineEar2(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2, long double e, long double r, Vector3 I)
{
	// I是贝塞尔坐标参数
	long double P = cos(I[1]), Q = sin(I[1]);
	long double X1 = x1, Y1 = P * y1 - Q * z1, Z1 = Q * y1 + P * z1;
	long double X2 = x2, Y2 = P * y2 - Q * z2, Z2 = Q * y2 + P * z2;
	COORDP p = lineEll(X1, Y1, Z1, X2, Y2, Z2, e, r);
	p.J = p.W = 100.0;
	if (p.D < 0)
	{
		return p;
	}
	p.J = rad2rrad(atan2(p.y, p.x) + I[0] - I[2]);
	p.W = atan(p.z / e / e / sqrt(p.x * p.x + p.y * p.y));
	return p;
}

NODE lineOvl(long double x1, long double y1, long double dx, long double dy, long double r, long double ba)
{
	long double A, B, C, D, L, t1, t2;
	NODE p = {};
	long double f = ba * ba;
	A = dx * dx + dy * dy / f;
	B = x1 * dx + y1 * dy / f;
	C = x1 * x1 + y1 * y1 / f - r * r;
	D = B * B - A * C;
	if (D < 0)
	{
		p.n = 0;
		return p;
	} // 判别式小于0无解
	if (!D)
	{
		p.n = 1;
	}
	else
	{
		p.n = 2;
	}
	D = sqrt(D);
	t1 = (-B + D) / A, t2 = (-B - D) / A;
	p.A = {x1 + dx * t1, y1 + dy * t1, 0};
	p.B = {x1 + dx * t2, y1 + dy * t2, 0};
	L = sqrt(dx * dx + dy * dy);
	p.R1 = L * fabs(t1); // x1到交点1的距离
	p.R2 = L * fabs(t2); // x1到交点2的距离
	return p;
}

NODE cirOvl(long double R, long double ba, long double R2, long double x0, long double y0)
{
	// 椭圆与圆的交点,R椭圆长半径,R2圆半径,x0,y0圆的圆心
	NODE re = {};
	long double d = sqrt(x0 * x0 + y0 * y0);
	long double sinB = y0 / d, cosB = x0 / d;
	long double cosA = (R * R + d * d - R2 * R2) / (2 * d * R);
	if (fabs(cosA) > 1)
	{
		re.n = 0;
		return re;
	} // 无解
	long double sinA = sqrt(1 - cosA * cosA);

	long double g, ba2 = ba * ba, C, S;
	int k;
	for (k = -1; k < 2; k += 2)
	{
		S = cosA * sinB + sinA * cosB * k;
		g = R - S * S * (1 / ba2 - 1) / 2;
		cosA = (g * g + d * d - R2 * R2) / (2 * d * g);
		if (fabs(cosA) > 1)
		{
			re.n = 0;
			return re;
		} // 无解
		sinA = sqrt(1 - cosA * cosA);
		C = cosA * cosB - sinA * sinB * k;
		S = cosA * sinB + sinA * cosB * k;
		if (k == 1)
		{
			re.A = {g * C, g * S, 0};
		}
		else
		{
			re.B = {g * C, g * S, 0};
		}
	}
	re.n = 2;
	return re;
}

_ECFAST ecFast(long double jd)
{
	// 快速日食搜索,jd为朔时间(J2000起算的儒略日数,不必很精确)
	_ECFAST re = {};
	long double t, t2, t3, t4;
	long double L, mB, mR, sR, vL, vB, vR;
	long double W = int64((jd + 8) / 29.5306) * PI * 2; // 合朔时的日月黄经差

	// 合朔时间计算,2000前+-4000年误差1小时以内，+-2000年小于10分钟
	t = (W + 1.08472) / 7771.37714500204; // 平朔时间
	re.jd = re.jdSuo = t * 36525;

	t2 = t * t, t3 = t2 * t, t4 = t3 * t;
	L = (93.2720993 + 483202.0175273 * t - 0.0034029 * t2 - t3 / 3526000 + t4 / 863310000) / 180 * PI;
	re.ac = 1, re.lx = "N";
	if (fabs(sin(L)) > 0.4)
	{
		return re; // 一般大于21度已不可能
	}

	t -= (-0.0000331 * t * t + 0.10976 * cos(0.785 + 8328.6914 * t)) / 7771;
	t2 = t * t;
	L = -1.084719 + 7771.377145013 * t - 0.0000331 * t2 +
		(22640 * cos(0.785 + 8328.6914 * t + 0.000152 * t2) + 4586 * cos(0.19 + 7214.063 * t - 0.000218 * t2) + 2370 * cos(2.54 + 15542.754 * t - 0.000070 * t2) + 769 * cos(3.1 + 16657.383 * t) + 666 * cos(1.5 + 628.302 * t) + 412 * cos(4.8 + 16866.93 * t) + 212 * cos(4.1 - 1114.63 * t) + 205 * cos(0.2 + 6585.76 * t) + 192 * cos(4.9 + 23871.45 * t) + 165 * cos(2.6 + 14914.45 * t) + 147 * cos(5.5 - 7700.39 * t) + 125 * cos(0.5 + 7771.38 * t) + 109 * cos(3.9 + 8956.99 * t) + 55 * cos(5.6 - 1324.18 * t) + 45 * cos(0.9 + 25195.62 * t) + 40 * cos(3.8 - 8538.24 * t) + 38 * cos(4.3 + 22756.82 * t) + 36 * cos(5.5 + 24986.07 * t) - 6893 * cos(4.669257 + 628.3076 * t) - 72 * cos(4.6261 + 1256.62 * t) - 43 * cos(2.67823 + 628.31 * t) * t + 21) / rad;

	t += (W - L) / (7771.38 - 914 * sin(0.7848 + 8328.691425 * t + 0.0001523 * t2) - 179 * sin(2.543 + 15542.7543 * t) - 160 * sin(0.1874 + 7214.0629 * t));
	re.jd = re.jdSuo = jd = t * 36525; // 朔时刻

	// 纬 52,15 (角秒)
	t2 = t * t / 10000, t3 = t2 * t / 10000;
	mB = 18461 * cos(0.0571 + 8433.46616 * t - 0.640 * t2 - 1 * t3) + 1010 * cos(2.413 + 16762.1576 * t + 0.88 * t2 + 25 * t3) + 1000 * cos(5.440 - 104.7747 * t + 2.16 * t2 + 26 * t3) + 624 * cos(0.915 + 7109.2881 * t + 0 * t2 + 7 * t3) + 199 * cos(1.82 + 15647.529 * t - 2.8 * t2 - 19 * t3) + 167 * cos(4.84 - 1219.403 * t - 1.5 * t2 - 18 * t3) + 117 * cos(4.17 + 23976.220 * t - 1.3 * t2 + 6 * t3) + 62 * cos(4.8 + 25090.849 * t + 2 * t2 + 50 * t3) + 33 * cos(3.3 + 15437.980 * t + 2 * t2 + 32 * t3) + 32 * cos(1.5 + 8223.917 * t + 4 * t2 + 51 * t3) + 30 * cos(1.0 + 6480.986 * t + 0 * t2 + 7 * t3) + 16 * cos(2.5 - 9548.095 * t - 3 * t2 - 43 * t3) + 15 * cos(0.2 + 32304.912 * t + 0 * t2 + 31 * t3) + 12 * cos(4.0 + 7737.590 * t) + 9 * cos(1.9 + 15019.227 * t) + 8 * cos(5.4 + 8399.709 * t) + 8 * cos(4.2 + 23347.918 * t) + 7 * cos(4.9 - 1847.705 * t) + 7 * cos(3.8 - 16133.856 * t) + 7 * cos(2.7 + 14323.351 * t);
	mB /= rad;

	// 距 106, 23 (千米)
	mR = 385001 + 20905 * cos(5.4971 + 8328.691425 * t + 1.52 * t2 + 25 * t3) + 3699 * cos(4.900 + 7214.06287 * t - 2.18 * t2 - 19 * t3) + 2956 * cos(0.972 + 15542.75429 * t - 0.66 * t2 + 6 * t3) + 570 * cos(1.57 + 16657.3828 * t + 3.0 * t2 + 50 * t3) + 246 * cos(5.69 - 1114.6286 * t - 3.7 * t2 - 44 * t3) + 205 * cos(1.02 + 14914.4523 * t - 1 * t2 + 6 * t3) + 171 * cos(3.33 + 23871.4457 * t + 1 * t2 + 31 * t3) + 152 * cos(4.94 + 6585.761 * t - 2 * t2 - 19 * t3) + 130 * cos(0.74 - 7700.389 * t - 2 * t2 - 25 * t3) + 109 * cos(5.20 + 7771.377 * t) + 105 * cos(2.31 + 8956.993 * t + 1 * t2 + 25 * t3) + 80 * cos(5.38 - 8538.241 * t + 2.8 * t2 + 26 * t3) + 49 * cos(6.24 + 628.302 * t) + 35 * cos(2.7 + 22756.817 * t - 3 * t2 - 13 * t3) + 31 * cos(4.1 + 16171.056 * t - 1 * t2 + 6 * t3) + 24 * cos(1.7 + 7842.365 * t - 2 * t2 - 19 * t3) + 23 * cos(3.9 + 24986.074 * t + 5 * t2 + 75 * t3) + 22 * cos(0.4 + 14428.126 * t - 4 * t2 - 38 * t3) + 17 * cos(2.0 + 8399.679 * t);
	mR /= 6378.1366;

	t = jd / 365250, t2 = t * t, t3 = t2 * t;
	// 误0.0002AU
	sR = 10001399 // 日地距离
		 + 167070 * cos(3.098464 + 6283.07585 * t) + 1396 * cos(3.0552 + 12566.1517 * t) + 10302 * cos(1.10749 + 6283.07585 * t) * t + 172 * cos(1.064 + 12566.152 * t) * t + 436 * cos(5.785 + 6283.076 * t) * t2 + 14 * cos(4.27 + 6283.08 * t) * t3;
	sR *= 1.49597870691 / 6378.1366 * 10;

	// 经纬速度
	t = jd / 36525.0;
	vL = 7771 // 月日黄经差速度
		 - 914 * sin(0.785 + 8328.6914 * t) - 179 * sin(2.543 + 15542.7543 * t) - 160 * sin(0.187 + 7214.0629 * t);
	vB = -755 * sin(0.057 + 8433.4662 * t) // 月亮黄纬速度
		 - 82 * sin(2.413 + 16762.1576 * t);
	vR = -27299 * sin(5.497 + 8328.691425 * t) - 4184 * sin(4.900 + 7214.06287 * t) - 7204 * sin(0.972 + 15542.75429 * t);
	vL /= 36525, vB /= 36525, vR /= 36525; // 每日速度

	long double gm = mR * sin(mB) * vL / sqrt(vB * vB + vL * vL), smR = sR - mR; // gm伽马值,smR日月距
	long double mk = 0.2725076, sk = 109.1222;
	long double f1 = (sk + mk) / smR, r1 = mk + f1 * mR; // tanf1半影锥角, r1半影半径
	long double f2 = (sk - mk) / smR, r2 = mk - f2 * mR; // tanf2本影锥角, r2本影半径
	long double b = 0.9972, Agm = fabs(gm), Ar2 = fabs(r2);
	long double fh2 = mR - mk / f2, h = Agm < 1 ? sqrt(1 - gm * gm) : 0; // fh2本影顶点的z坐标
	long double ls1, ls2, ls3, ls4;

	if (fh2 < h)
	{
		re.lx = "T";
	}
	else
	{
		re.lx = "A";
	}

	ls1 = Agm - (b + r1);
	if (fabs(ls1) < 0.016)
	{
		re.ac = 0; // 无食分界
	}
	ls2 = Agm - (b + Ar2);
	if (fabs(ls2) < 0.016)
	{
		re.ac = 0; // 偏食分界
	}
	ls3 = Agm - (b);
	if (fabs(ls3) < 0.016)
	{
		re.ac = 0; // 无中心食分界
	}
	ls4 = Agm - (b - Ar2);
	if (fabs(ls4) < 0.016)
	{
		re.ac = 0; // 有中心食分界(但本影未全部进入)
	}

	if (ls1 > 0)
	{
		re.lx = "N"; // 无日食
	}
	else if (ls2 > 0)
	{
		re.lx = "P"; // 偏食
	}
	else if (ls3 > 0)
	{
		re.lx += "0"; // 无中心
	}
	else if (ls4 > 0)
	{
		re.lx += "1"; // 有中心(本影未全部进入)
	}
	else
	{ // 本影全进入
		if (fabs(fh2 - h) < 0.019)
		{
			re.ac = 0;
		}
		if (fabs(fh2) < h)
		{
			long double dr = vR * h / vL / mR;
			long double H1 = mR - dr - mk / f2; // 入点影锥z坐标
			long double H2 = mR + dr - mk / f2; // 出点影锥z坐标
			if (H1 > 0)
			{
				re.lx = "H3"; // 环全全
			}
			if (H2 > 0)
			{
				re.lx = "H2"; // 全全环
			}
			if (H1 > 0 && H2 > 0)
			{
				re.lx = "H"; // 环全环
			}
			if (fabs(H1) < 0.019)
			{
				re.ac = 0;
			}
			if (fabs(H2) < 0.019)
			{
				re.ac = 0;
			}
		}
	}
	return re;
}

long double moonIll(long double t)
{
	// 月亮被照亮部分的比例
	long double t2 = t * t, t3 = t2 * t, t4 = t3 * t;
	long double D, M, m, a, dm = PI / 180;
	D = (297.8502042 + 445267.1115168 * t - 0.0016300 * t2 + t3 / 545868 - t4 / 113065000) * dm; // 日月平距角
	M = (357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000) * dm;				 // 太阳平近点
	m = (134.9634114 + 477198.8676313 * t + 0.0089970 * t2 + t3 / 69699 - t4 / 14712000) * dm;	 // 月亮平近点
	a = PI - D + (-6.289 * sin(m) + 2.100 * sin(M) - 1.274 * sin(D * 2 - m) - 0.658 * sin(D * 2) - 0.214 * sin(m * 2) - 0.110 * sin(D)) * dm;
	return (1 + cos(a)) / 2;
}

long double moonRad(long double r, long double h)
{
	// 转入地平纬度及地月质心距离,返回站心视半径(角秒)
	return cs_sMoon / r * (1 + sin(h) * cs_rEar / r);
}
