/* -*- coding: utf-8 -*-
 * test_apparent_catalog.c —— 标准星表 + 在轨视方向 示例/测试程序
 *
 * 用法:
 *   test_apparent_catalog <catalog.bin> [year month day hour_utc]
 *                        [rx ry rz km] [vx vy vz km/s] [ncsv]
 * 缺省: 2026-11-15 12:00 UTC, 550 km 圆轨道示例状态 (GCRS), 导出 2000 星。
 *
 * 输出:
 *   - 星表摘要 + 最亮 10 星明细（J2016 方向 vs 视方向, 各项改正量）
 *   - 全星表改正量统计（均值/中位数/P99/最大, 与误差分析表同口径）
 *   - c_sample_output.csv（前 ncsv 星, 供与 Python 参考实现逐星比对）
 *
 * 版本: v1.1.261115-5
 */
#include "apparent_catalog.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int cmp_d(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static void stats(const char *name, double *v, int n, const char *unit)
{
    double mean = 0.0, median, p99;
    int i;
    for (i = 0; i < n; i++) mean += v[i];
    mean /= n;
    qsort(v, (size_t)n, sizeof(double), cmp_d);
    median = (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    p99 = v[(int)(0.99 * (n - 1))];
    printf("  %-14s mean=%10.4f median=%10.4f P99=%10.4f max=%10.4f  %s\n",
           name, mean, median, p99, v[n - 1], unit);
}

int main(int argc, char *argv[])
{
    ac_catalog_t cat;
    int st, i, n, ncsv = 2000;
    int yr = 2026, mo = 11, dy = 15;
    double hr = 12.0;
    double r[3] = { 6928.137, 0.0, 0.0 };        /* 550 km 圆轨道示例 */
    double v[3] = { 0.0, -1.0065, 7.5219 };
    double jd, *ra = NULL, *dec = NULL;
    ac_breakdown_t *brk = NULL;
    double *tmp;
    clock_t t0;
    double dt_ms;

    if (argc < 2) {
        fprintf(stderr, "用法: %s <catalog.bin> [y m d h] [rx ry rz] [vx vy vz] [ncsv]\n",
                argv[0]);
        return 2;
    }
    if (argc >= 6) { yr = atoi(argv[2]); mo = atoi(argv[3]); dy = atoi(argv[4]); hr = atof(argv[5]); }
    if (argc >= 9) { r[0] = atof(argv[6]); r[1] = atof(argv[7]); r[2] = atof(argv[8]); }
    if (argc >= 12) { v[0] = atof(argv[9]); v[1] = atof(argv[10]); v[2] = atof(argv[11]); }
    if (argc >= 13) ncsv = atoi(argv[12]);

    printf("apparent_catalog %s (test)\n", AC_VERSION_STRING);

    st = ac_load_catalog(argv[1], &cat);
    if (st != AC_OK) {
        fprintf(stderr, "星表读取失败: %s (code %d)\n", argv[1], st);
        return 1;
    }
    n = cat.hdr.n_stars;
    printf("星表: %s  N=%d  G<=%.1f  源历元 JD%.1f → 传播历元 JD%.1f\n",
           argv[1], n, cat.hdr.g_max,
           cat.hdr.jd_epoch_source, cat.hdr.jd_epoch_propagated);

    jd = ac_utc_to_jd(yr, mo, dy, hr);
    if (!isfinite(jd)) {
        fprintf(stderr, "观测时刻非法\n");
        return 1;
    }
    printf("观测: %04d-%02d-%02d %05.2fh UTC (JD %.6f)\n", yr, mo, dy, hr, jd);
    printf("卫星: r=[%.3f %.3f %.3f] km, v=[%.4f %.4f %.4f] km/s, |v|=%.4f km/s (GCRS)\n",
           r[0], r[1], r[2], v[0], v[1], v[2],
           sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));

    ra  = (double *)malloc((size_t)n * sizeof(double));
    dec = (double *)malloc((size_t)n * sizeof(double));
    brk = (ac_breakdown_t *)malloc((size_t)n * sizeof(ac_breakdown_t));
    if (!ra || !dec || !brk) { fprintf(stderr, "内存不足\n"); return 1; }

    t0 = clock();
    st = ac_apparent_catalog(&cat, jd, r, v, AC_DEFL_SUN_EARTH, ra, dec, NULL, brk);
    dt_ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (st != AC_OK) {
        fprintf(stderr, "计算失败 code %d\n", st);
        return 1;
    }
    printf("整表计算: %.1f ms (%d 星, %.2f us/星)\n", dt_ms, n, 1e3 * dt_ms / n);

    /* 最亮 10 星（按 G 部分选择） */
    printf("\n最亮 10 星（J2016 星表方向 → 在轨视方向）:\n");
    printf("%19s %5s %13s %12s %13s %12s %9s %10s %9s %9s\n",
           "source_id", "G", "ra_2016", "dec_2016", "ra_app", "dec_app",
           "总改正\"", "自行+视差\"", "偏折mas", "光行差\"");
    {
        int idx[10];
        for (i = 0; i < 10 && i < n; i++) idx[i] = i;
        /* 简单选择: 前 10 个最小 G */
        {
            int k, m, mi;
            for (k = 0; k < n; k++) {
                m = k;
                for (i = 0; i < 10; i++) if (idx[i] == k) { m = -1; break; }
                if (m < 0) continue;
                for (i = 0; i < 10; i++)
                    if (cat.stars[k].g_mag < cat.stars[idx[i]].g_mag) { idx[i] = k; break; }
            }
            /* 冒泡排序这 10 个 */
            for (i = 0; i < 9; i++)
                for (m = 0; m < 9 - i; m++)
                    if (cat.stars[idx[m + 1]].g_mag < cat.stars[idx[m]].g_mag) {
                        mi = idx[m]; idx[m] = idx[m + 1]; idx[m + 1] = mi;
                    }
        }
        for (i = 0; i < 10 && i < n; i++) {
            int k = idx[i];
            printf("%19lld %5.2f %13.6f %12.6f %13.6f %12.6f %9.3f %10.3f %9.2f %9.3f\n",
                   (long long)cat.stars[k].source_id, cat.stars[k].g_mag,
                   cat.stars[k].ra_j2016_deg, cat.stars[k].dec_j2016_deg,
                   ra[k], dec[k], brk[k].total_as, brk[k].pm_parallax_as,
                   brk[k].deflection_mas, brk[k].aberration_as);
        }
    }

    /* 全星表统计（注意 stats 会排序, 复制后使用） */
    printf("\n全星表改正量统计 (%d 星):\n", n);
    tmp = (double *)malloc((size_t)n * sizeof(double));
    for (i = 0; i < n; i++) tmp[i] = brk[i].pm_parallax_as;
    stats("自行+视差", tmp, n, "arcsec");
    for (i = 0; i < n; i++) tmp[i] = brk[i].deflection_mas;
    stats("引力偏折", tmp, n, "mas");
    for (i = 0; i < n; i++) tmp[i] = brk[i].aberration_as;
    stats("光行差", tmp, n, "arcsec");
    for (i = 0; i < n; i++) tmp[i] = brk[i].total_as;
    stats("总改正", tmp, n, "arcsec");
    free(tmp);

    /* 导出交叉验证样本 */
    if (ncsv > n) ncsv = n;
    {
        FILE *f = fopen("c_sample_output.csv", "w");
        if (!f) { perror("c_sample_output.csv"); return 1; }
        fprintf(f, "index,source_id,ra_app_deg,dec_app_deg,total_as,"
                   "pm_px_as,defl_mas,abbr_as\n");
        for (i = 0; i < ncsv; i++)
            fprintf(f, "%d,%lld,%.12f,%.12f,%.9f,%.9f,%.9f,%.9f\n",
                    i, (long long)cat.stars[i].source_id, ra[i], dec[i],
                    brk[i].total_as, brk[i].pm_parallax_as,
                    brk[i].deflection_mas, brk[i].aberration_as);
        fclose(f);
        printf("\n交叉验证样本已导出: c_sample_output.csv (%d 星)\n", ncsv);
    }

    free(ra); free(dec); free(brk);
    ac_free_catalog(&cat);
    return 0;
}
