/* -*- coding: utf-8 -*-
 * test_eop.c —— IERS EOP C04 接入测试
 *
 * 用法: test_eop <eopc04 文件> [catalog.bin]
 *
 * 测试项:
 *   T1 加载: 行数/时间范围/最新一行数值
 *   T2 表列时刻插值 == 表值 (<=1e-12)
 *   T3 数据范围外 -> AC_W_EOP_RANGE 警告
 *   T4 ITRS->GCRS->ITRS 往返一致 (位置 <=1e-9 km, 速度 <=1e-9 km/s)
 *   T5 极移效应量化: 有/无极移的指向差 (应 ~ |PM| 量级, 0.1-0.5")
 *   T6 导出 eop_rotation_sample.csv 供 Python(erfa) 交叉验证
 *   T7 (可选, 给定星表) ac_apparent_catalog_itrs 冒烟测试
 *
 * 版本: v1.1.261115-5
 */
#include "apparent_catalog.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;

static void check(const char *name, int ok, const char *detail)
{
    printf("  [%s] %s  %s\n", ok ? "PASS" : "FAIL", name, detail ? detail : "");
    if (!ok) g_fail++;
}

static double norm3(const double v[3])
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

int main(int argc, char *argv[])
{
    ac_eop_t eop;
    int st, i, n_ep;
    double jd0, mjd;
    double r_i[3] = { 7000.0, 1000.0, -1500.0 };      /* km, ITRS 示例 */
    double v_i[3] = { -1.1, 6.9, 0.3 };               /* km/s */
    double r_g[3], v_g[3], r_b[3], v_b[3];
    double xp, yp, dut1, dx, dy;

    if (argc < 2) {
        fprintf(stderr, "用法: %s <eopc04 文件> [catalog.bin]\n", argv[0]);
        return 2;
    }

    printf("apparent_catalog %s (EOP test)\n", AC_VERSION_STRING);

    /* ---- T1 加载 ---- */
    st = ac_eop_load(argv[1], &eop);
    {
        char buf[256];
        const ac_eop_record_t *l = &eop.rec[eop.n - 1];
        snprintf(buf, sizeof(buf),
                 "N=%d, MJD %.1f..%.1f, 最新: x=%.6f\" y=%.6f\" DUT1=%.7fs "
                 "dX=%.6f\" dY=%.6f\"",
                 eop.n, eop.mjd_first, eop.mjd_last,
                 l->xp_as, l->yp_as, l->dut1_s, l->dx_as, l->dy_as);
        check("T1 ac_eop_load", st == AC_OK && eop.n > 20000, buf);
    }

    /* ---- T2 表列时刻插值 == 表值 ---- */
    mjd = eop.rec[eop.n - 1].mjd_utc;
    st = ac_eop_get(&eop, mjd, &xp, &yp, &dut1, &dx, &dy);
    {
        const ac_eop_record_t *l = &eop.rec[eop.n - 1];
        double d = fabs(xp - l->xp_as) + fabs(yp - l->yp_as)
                 + fabs(dut1 - l->dut1_s);
        check("T2 端点时刻插值", st == AC_OK && d < 1e-12, "diff<1e-12");
    }
    /* 中间表列时刻 */
    mjd = eop.rec[1000].mjd_utc;
    st = ac_eop_get(&eop, mjd, &xp, &yp, &dut1, &dx, &dy);
    {
        const ac_eop_record_t *r = &eop.rec[1000];
        double d = fabs(xp - r->xp_as) + fabs(yp - r->yp_as)
                 + fabs(dut1 - r->dut1_s);
        check("T2 内部表列时刻插值", st == AC_OK && d < 1e-12, "diff<1e-12");
    }

    /* ---- T3 越界外推（长期产品策略: DUT1 用 -LOD, 极移用线性趋势） ---- */
    {
        const ac_eop_record_t *l = &eop.rec[eop.n - 1];
        double d = 100.0, exp_dut1;
        st = ac_eop_get(&eop, eop.mjd_last + d, &xp, &yp, &dut1, &dx, &dy);
        exp_dut1 = l->dut1_s - l->lod_s * d;
        check("T3 越界返回警告", st == AC_W_EOP_RANGE, "外推值+警告");
        check("T3 DUT1 按 -LOD 外推",
              fabs(dut1 - exp_dut1) < 1e-12, "dut1=dut1_end-LOD*d");
        check("T3 极移线性趋势外推量级合理",
              fabs(xp - l->xp_as) < 0.5 && fabs(yp - l->yp_as) < 0.5,
              "100 天外推偏离端点 <0.5\"");
        st = ac_eop_get(&eop, eop.mjd_first - 5.0, &xp, &yp, &dut1, &dx, &dy);
        check("T3 首端越界返回警告", st == AC_W_EOP_RANGE, "端点保持");
    }

    /* ---- T4 往返一致 ---- */
    {
        double eps_r, eps_v;
        jd0 = ac_utc_to_jd(2026, 6, 1, 7.5);
        st = ac_itrs_to_gcrs(&eop, jd0, r_i, v_i, r_g, v_g);
        ac_gcrs_to_itrs(&eop, jd0, r_g, v_g, r_b, v_b);
        eps_r = fabs(r_b[0] - r_i[0]) + fabs(r_b[1] - r_i[1]) + fabs(r_b[2] - r_i[2]);
        eps_v = fabs(v_b[0] - v_i[0]) + fabs(v_b[1] - v_i[1]) + fabs(v_b[2] - v_i[2]);
        check("T4 ITRS<->GCRS 往返",
              st == AC_OK && eps_r < 1e-9 && eps_v < 1e-9,
              "位置差<1e-9 km, 速度差<1e-9 km/s");
        printf("      (2026-06-01 07:30 UTC: r_g=[%.6f %.6f %.6f] km, "
               "v_g=[%.9f %.9f %.9f] km/s)\n",
               r_g[0], r_g[1], r_g[2], v_g[0], v_g[1], v_g[2]);
    }

    /* ---- T5 极移效应 ---- */
    {
        double r_nopm[3], v_nopm[3], ang;
        int k;
        /* 暂存并清零极移（测完恢复, 避免污染共享记录数组） */
        double *xp_save = (double *)malloc((size_t)eop.n * sizeof(double));
        double *yp_save = (double *)malloc((size_t)eop.n * sizeof(double));
        for (k = 0; k < eop.n; k++) {
            xp_save[k] = eop.rec[k].xp_as;
            yp_save[k] = eop.rec[k].yp_as;
            eop.rec[k].xp_as = 0.0;
            eop.rec[k].yp_as = 0.0;
        }
        ac_itrs_to_gcrs(&eop, jd0, r_i, v_i, r_nopm, v_nopm);
        for (k = 0; k < eop.n; k++) {
            eop.rec[k].xp_as = xp_save[k];
            eop.rec[k].yp_as = yp_save[k];
        }
        free(xp_save);
        free(yp_save);
        /* 两方向夹角 */
        ang = atan2(sqrt(pow(r_g[1] * r_nopm[2] - r_g[2] * r_nopm[1], 2)
                       + pow(r_g[2] * r_nopm[0] - r_g[0] * r_nopm[2], 2)
                       + pow(r_g[0] * r_nopm[1] - r_g[1] * r_nopm[0], 2)),
                    r_g[0] * r_nopm[0] + r_g[1] * r_nopm[1] + r_g[2] * r_nopm[2]);
        ang *= (180.0 / 3.14159265358979323846) * 3600.0;
        check("T5 极移效应量级", ang > 0.05 && ang < 1.0,
              "有/无极移指向差(角秒)");
    }

    /* ---- T6 导出旋转样本（Python 交叉验证用） ---- */
    {
        double eps[6][3] = {
            {2026, 6, 1}, {2026, 6, 1}, {2026, 2, 15},
            {2024, 11, 27}, {2026, 10, 14}, {2026, 10, 14}
        };
        double hrs[6] = { 0.0, 7.5, 12.0, 23.999, 6.0, 18.0 };
        FILE *f = fopen("eop_rotation_sample.csv", "w");
        if (!f) { perror("eop_rotation_sample.csv"); return 1; }
        fprintf(f, "jd_utc,rgx,rgy,rgz,vgx,vgy,vgz\n");
        n_ep = 6;
        for (i = 0; i < n_ep; i++) {
            double jd = ac_utc_to_jd((int)eps[i][0], (int)eps[i][1],
                                     (int)eps[i][2], hrs[i]);
            ac_itrs_to_gcrs(&eop, jd, r_i, v_i, r_g, v_g);
            fprintf(f, "%.9f,%.9f,%.9f,%.9f,%.12f,%.12f,%.12f\n",
                    jd, r_g[0], r_g[1], r_g[2], v_g[0], v_g[1], v_g[2]);
        }
        fclose(f);
        printf("  [INFO] 旋转样本已导出: eop_rotation_sample.csv (%d 组)\n", n_ep);
    }

    /* ---- T7 星表 + ITRS 状态一步计算 ---- */
    if (argc >= 3) {
        ac_catalog_t cat;
        st = ac_load_catalog(argv[2], &cat);
        if (st == AC_OK) {
            double ra, dec, u[3];
            ac_breakdown_t brk1;
            ac_star_record_t *rec = &cat.stars[0];
            st = ac_apparent_catalog_itrs(&cat, &eop,
                    ac_utc_to_jd(2026, 6, 1, 0.0), r_i, v_i,
                    AC_DEFL_SUN_EARTH, NULL, NULL, NULL, NULL);
            check("T7 ac_apparent_catalog_itrs", st == AC_OK, "整表 OK");
            /* 单星视方向: 先转 GCRS 再算 */
            ac_itrs_to_gcrs(&eop, ac_utc_to_jd(2026, 6, 1, 0.0),
                            r_i, v_i, r_g, v_g);
            ac_apparent_star(rec, ac_utc_to_jd(2026, 6, 1, 0.0),
                             r_g, v_g, AC_DEFL_SUN_EARTH, u, &ra, &dec, &brk1);
            printf("      首星(source_id=%lld): 视方向 ra=%.6f dec=%.6f, "
                   "总改正=%.3f\"\n",
                   (long long)rec->source_id, ra, dec, brk1.total_as);
            ac_free_catalog(&cat);
        } else {
            check("T7 星表读取", 0, "读取失败");
        }
    }

    ac_eop_free(&eop);
    printf(g_fail ? "\nEOP 测试: FAIL (%d 项未过)\n" : "\nEOP 测试: 全部通过\n", g_fail);
    return g_fail ? 1 : 0;
}
