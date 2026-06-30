// g++ -std=c++17 -I../src -I../capi scripts/verify_bazi_shishen.cpp capi/sxwnl_capi.cpp src/*.cpp -o /tmp/verify_bazi
#include "sxwnl_capi.h"
#include <cstdio>
#include <cstring>

static void test(int hour) {
    SxwnlBaziParams bp{};
    strcpy(bp.name, "Test");
    bp.gender = true;
    bp.year = 2017;
    bp.month = 10;
    bp.day = 21;
    bp.hour = hour;
    bp.minute = 0;
    bp.second = 0;
    bp.longitude = 116.4;
    bp.latitude = 39.9;
    bp.lifa = 11;
    SxwnlBazi bz = sxwnl_bazi_create(&bp);
    if (!bz) {
        printf("hour=%d create failed\n", hour);
        return;
    }
    SxwnlPillar pls[4];
    sxwnl_bazi_get_si_zhu(bz, pls);
    printf("hour=%d SiZhu:", hour);
    for (int i = 0; i < 4; i++) printf(" %s%s", pls[i].tian_gan, pls[i].di_zhi);
    printf("\n");
    SxwnlBaziColumn cols[4];
    sxwnl_bazi_get_columns(bz, cols);
    const char *names[] = {"年", "月", "日", "时"};
    for (int i = 0; i < 4; i++) {
        printf("  %s柱: %s%s 主星=%s 藏干:", names[i], cols[i].gan, cols[i].zhi, cols[i].gan_shi_shen);
        for (int j = 0; j < cols[i].cang_gan_count; j++)
            printf(" %s(%s)", cols[i].cang_gan[j], cols[i].cang_gan_shi_shen[j]);
        printf("\n");
    }
    sxwnl_bazi_free(bz);
}

int main() {
    test(11);
    test(12);
    return 0;
}
