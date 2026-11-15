/* -*- coding: utf-8 -*-
 * eop_c04.c —— IERS EOP C04 (IAU2000) 加载/插值 + ITRS(ECEF) 状态支持
 *
 * 数据文件: IERS EOP 产品中心 (巴黎天文台) 发布的
 *   eopc04_IAU2000.62-now   （update_eop.py 下载/更新）
 * 格式: 日期(年 月 日) MJD x(") y(") UT1-UTC(s) LOD(s) dX(") dY(") ...误差列
 *
 * ITRS→GCRS 旋转 = erfaC2t06a（EOP 极移 + DUT1 + IAU2006/2000A 岁差章动,
 * 无旋转原点 CIO 体系）; 速度含地球自转项 v_g = Rᵀ(v_i + ω×r_i)。
 *
 * 版本: v1.1.261115-6
 */
#include "apparent_catalog.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "erfa.h"
#include "erfam.h"

#define AC_MJD0 2400000.5
#define AC_EARTH_W 7.292115146706979e-5   /* 地球自转角速度 rad/s (IERS 2010) */

/* ------------------------------------------------------------------ */
/* C04 文件解析                                                        */
/* ------------------------------------------------------------------ */
int ac_eop_load(const char *path, ac_eop_t *eop)
{
    FILE *f;
    char line[512];
    int cap = 4096, n = 0;

    if (!path || !eop) return AC_E_ARG;
    memset(eop, 0, sizeof(*eop));

    f = fopen(path, "r");
    if (!f) return AC_E_IO;

    eop->rec = (ac_eop_record_t *)malloc((size_t)cap * sizeof(ac_eop_record_t));
    if (!eop->rec) { fclose(f); return AC_E_IO; }

    while (fgets(line, sizeof(line), f)) {
        int iy, im, id, k;
        double mjd, x, y, dut1, lod, dx, dy;
        k = sscanf(line, "%d %d %d %lf %lf %lf %lf %lf %lf %lf",
                   &iy, &im, &id, &mjd, &x, &y, &dut1, &lod, &dx, &dy);
        if (k != 10) continue;                       /* 表头/注释行 */
        if (iy < 1962 || iy > 2100) continue;
        if (mjd < 30000.0 || mjd > 90000.0) continue;
        if (n >= cap) {
            cap *= 2;
            eop->rec = (ac_eop_record_t *)realloc(eop->rec,
                        (size_t)cap * sizeof(ac_eop_record_t));
            if (!eop->rec) { fclose(f); return AC_E_IO; }
        }
        eop->rec[n].mjd_utc = mjd;
        eop->rec[n].xp_as   = x;
        eop->rec[n].yp_as   = y;
        eop->rec[n].dut1_s  = dut1;
        eop->rec[n].lod_s   = lod;
        eop->rec[n].dx_as   = dx;
        eop->rec[n].dy_as   = dy;
        n++;
    }
    fclose(f);

    if (n < 2) { free(eop->rec); eop->rec = NULL; return AC_E_EOP; }

    /* 保险: 按 mjd 升序（C04 本身有序, 仅防御） */
    {
        int i, j;
        for (i = 1; i < n; i++) {
            ac_eop_record_t key = eop->rec[i];
            if (eop->rec[i - 1].mjd_utc > key.mjd_utc) {
                for (j = i - 1; j >= 0 && eop->rec[j].mjd_utc > key.mjd_utc; j--)
                    eop->rec[j + 1] = eop->rec[j];
                eop->rec[j + 1] = key;
            }
        }
    }

    eop->n = n;
    eop->mjd_first = eop->rec[0].mjd_utc;
    eop->mjd_last  = eop->rec[n - 1].mjd_utc;
    return AC_OK;
}

void ac_eop_free(ac_eop_t *eop)
{
    if (eop && eop->rec) {
        free(eop->rec);
        eop->rec = NULL;
    }
}

/* ------------------------------------------------------------------ */
/* 插值与外推（本项目固定使用长期产品 C04, 不使用 IERS 快速产品）        */
/* ------------------------------------------------------------------ */

/* 末端 [i_end-span+1, i_end] 窗口最小二乘拟合 y = a + b·(mjd − t0), 返回外推值 */
static double eop_fit_extrapolate(const double *mjd, const double *y,
                                  int i_end, int span, double d)
{
    int i0 = i_end - span + 1;
    int i, cnt = 0;
    double t0 = mjd[i_end];
    double sx = 0, sy = 0, sxx = 0, sxy = 0, den, b, a;

    if (i0 < 0) i0 = 0;
    for (i = i0; i <= i_end; i++) {
        double x = mjd[i] - t0;
        sx += x; sy += y[i]; sxx += x * x; sxy += x * y[i];
        cnt++;
    }
    den = (double)cnt * sxx - sx * sx;
    if (den == 0.0) return y[i_end];
    b = ((double)cnt * sxy - sx * sy) / den;
    a = sy / (double)cnt - b * (sx / (double)cnt);
    return a + b * d;
}

int ac_eop_get(const ac_eop_t *eop, double mjd_utc,
               double *xp_as, double *yp_as, double *dut1_s,
               double *dx_as, double *dy_as)
{
    int lo, hi, mid;
    double t, a0, a1;
    const ac_eop_record_t *r = eop->rec;

    if (!eop || !eop->rec || eop->n < 2) return AC_E_EOP;

    if (mjd_utc <= eop->mjd_first) {
        *xp_as = r[0].xp_as; *yp_as = r[0].yp_as; *dut1_s = r[0].dut1_s;
        *dx_as = r[0].dx_as; *dy_as = r[0].dy_as;
        return (mjd_utc < eop->mjd_first) ? AC_W_EOP_RANGE : AC_OK;
    }
    if (mjd_utc == eop->mjd_last) {
        const ac_eop_record_t *l = &r[eop->n - 1];
        *xp_as = l->xp_as; *yp_as = l->yp_as; *dut1_s = l->dut1_s;
        *dx_as = l->dx_as; *dy_as = l->dy_as;
        return AC_OK;
    }
    if (mjd_utc > eop->mjd_last) {
        /* 末端之外: 长期产品外推策略（返回 AC_W_EOP_RANGE 告知调用方）
           - DUT1: d(UT1-UTC)/dt = -LOD （解析关系）
           - 极移: 末端 30 天线性趋势
           - dX/dY: 毫角秒级, 端点保持 */
        const ac_eop_record_t *l = &r[eop->n - 1];
        double d = mjd_utc - eop->mjd_last;
        int i, i0 = eop->n - 30;
        if (i0 < 0) i0 = 0;
        double *mm = (double *)malloc((size_t)(eop->n - i0) * sizeof(double));
        double *xx = (double *)malloc((size_t)(eop->n - i0) * sizeof(double));
        double *yy = (double *)malloc((size_t)(eop->n - i0) * sizeof(double));
        if (!mm || !xx || !yy) {
            free(mm); free(xx); free(yy);
            *xp_as = l->xp_as; *yp_as = l->yp_as; *dut1_s = l->dut1_s;
            *dx_as = l->dx_as; *dy_as = l->dy_as;
            return AC_W_EOP_RANGE;
        }
        for (i = i0; i < eop->n; i++) {
            mm[i - i0] = r[i].mjd_utc;
            xx[i - i0] = r[i].xp_as;
            yy[i - i0] = r[i].yp_as;
        }
        *xp_as  = eop_fit_extrapolate(mm, xx, eop->n - 1 - i0, 30, d);
        *yp_as  = eop_fit_extrapolate(mm, yy, eop->n - 1 - i0, 30, d);
        *dut1_s = l->dut1_s - l->lod_s * d;
        *dx_as = l->dx_as; *dy_as = l->dy_as;
        free(mm); free(xx); free(yy);
        return AC_W_EOP_RANGE;
    }

    /* 二分: r[lo].mjd <= mjd_utc < r[hi].mjd */
    lo = 0;
    hi = eop->n - 1;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (r[mid].mjd_utc <= mjd_utc) lo = mid; else hi = mid;
    }
    a0 = r[lo].mjd_utc;
    a1 = r[hi].mjd_utc;
    t = (mjd_utc - a0) / (a1 - a0);

    *xp_as  = r[lo].xp_as  + t * (r[hi].xp_as  - r[lo].xp_as);
    *yp_as  = r[lo].yp_as  + t * (r[hi].yp_as  - r[lo].yp_as);
    *dut1_s = r[lo].dut1_s + t * (r[hi].dut1_s - r[lo].dut1_s);
    *dx_as  = r[lo].dx_as  + t * (r[hi].dx_as  - r[lo].dx_as);
    *dy_as  = r[lo].dy_as  + t * (r[hi].dy_as  - r[lo].dy_as);
    return AC_OK;
}

/* ------------------------------------------------------------------ */
/* ITRS <-> GCRS                                                       */
/* ------------------------------------------------------------------ */
int ac_itrs_to_gcrs(const ac_eop_t *eop, double jd_utc,
                    double r_itrs[3], double v_itrs[3],
                    double r_gcrs[3], double v_gcrs[3])
{
    double jd_tt, xp, yp, dut1, dx, dy, ut1, rc2t[3][3];
    double wvec[3] = { 0.0, 0.0, AC_EARTH_W };
    double wxr[3], vi[3];
    int i, st;

    if (!eop || !eop->rec) return AC_E_EOP;

    jd_tt = ac_jd_utc_to_jd_tt(jd_utc);
    if (!isfinite(jd_tt)) return AC_E_TIME;

    st = ac_eop_get(eop, jd_utc - AC_MJD0, &xp, &yp, &dut1, &dx, &dy);
    if (st < 0) return st;

    ut1 = jd_utc + dut1 / 86400.0;
    eraC2t06a(jd_tt, 0.0, ut1, 0.0,
              xp * ERFA_DAS2R, yp * ERFA_DAS2R, rc2t);   /* 天球→地固 */

    eraPxp(wvec, r_itrs, wxr);               /* ω × r, ITRS 中 */
    for (i = 0; i < 3; i++) vi[i] = v_itrs[i] + wxr[i];

    /* r_g = Rᵀ r_i,  v_g = Rᵀ vi */
    for (i = 0; i < 3; i++) {
        r_gcrs[i] = rc2t[0][i] * r_itrs[0]
                  + rc2t[1][i] * r_itrs[1]
                  + rc2t[2][i] * r_itrs[2];
        v_gcrs[i] = rc2t[0][i] * vi[0]
                  + rc2t[1][i] * vi[1]
                  + rc2t[2][i] * vi[2];
    }
    return st;
}

int ac_gcrs_to_itrs(const ac_eop_t *eop, double jd_utc,
                    double r_gcrs[3], double v_gcrs[3],
                    double r_itrs[3], double v_itrs[3])
{
    double jd_tt, xp, yp, dut1, dx, dy, ut1, rc2t[3][3];
    double wvec[3] = { 0.0, 0.0, AC_EARTH_W };
    double wxr[3];
    int i, st;

    if (!eop || !eop->rec) return AC_E_EOP;

    jd_tt = ac_jd_utc_to_jd_tt(jd_utc);
    if (!isfinite(jd_tt)) return AC_E_TIME;

    st = ac_eop_get(eop, jd_utc - AC_MJD0, &xp, &yp, &dut1, &dx, &dy);
    if (st < 0) return st;

    ut1 = jd_utc + dut1 / 86400.0;
    eraC2t06a(jd_tt, 0.0, ut1, 0.0,
              xp * ERFA_DAS2R, yp * ERFA_DAS2R, rc2t);

    /* r_i = R r_g,  v_i = R v_g − ω×r_i */
    for (i = 0; i < 3; i++) {
        r_itrs[i] = rc2t[i][0] * r_gcrs[0]
                  + rc2t[i][1] * r_gcrs[1]
                  + rc2t[i][2] * r_gcrs[2];
        v_itrs[i] = rc2t[i][0] * v_gcrs[0]
                  + rc2t[i][1] * v_gcrs[1]
                  + rc2t[i][2] * v_gcrs[2];
    }
    eraPxp(wvec, r_itrs, wxr);
    for (i = 0; i < 3; i++) v_itrs[i] -= wxr[i];
    return st;
}

/* ------------------------------------------------------------------ */
/* ITRS 卫星状态 → 在轨观测星表                                        */
/* ------------------------------------------------------------------ */
int ac_apparent_catalog_itrs(const ac_catalog_t *cat,
                             const ac_eop_t *eop,
                             double jd_utc,
                             double r_sat_itrs_km[3],
                             double v_sat_itrs_kms[3],
                             unsigned defl_flags,
                             double *ra_app_deg,
                             double *dec_app_deg,
                             double *unit_gcrs,
                             ac_breakdown_t *brk)
{
    double r_g[3], v_g[3];
    int st, st2;

    st = ac_itrs_to_gcrs(eop, jd_utc, r_sat_itrs_km, v_sat_itrs_kms,
                         r_g, v_g);
    if (st < 0) return st;

    st2 = ac_apparent_catalog(cat, jd_utc, r_g, v_g, defl_flags,
                              ra_app_deg, dec_app_deg, unit_gcrs, brk);
    if (st2 < 0) return st2;
    return st;    /* 透传 EOP 越界警告 */
}
