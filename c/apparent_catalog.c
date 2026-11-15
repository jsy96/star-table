/* -*- coding: utf-8 -*-
 * apparent_catalog.c —— 观测时刻 + 卫星状态矢量 → 在轨观测星表（GCRS 视方向）
 *
 * 算法链路逐条使用 ERFA(SOFA) 标准例程:
 *   eraEpv00 地球日心/质心 pv 星历    eraPmpx 自行+视差
 *   eraLd    引力光线偏折             eraAb   光行差
 *   eraCal2jd/eraJd2cal/eraDat 时间尺度 UTC→TT
 * 与 SOFA iauATCIQ 的区别仅在于不施加 bias-岁差-章动旋转（全天空共模旋转
 * 由姿态解算吸收），因此输出为 GCRS 视方向，可直接与星敏感器测量方向比对。
 *
 * 版权说明: 本文件调用 third_party/erfa 下原样收录的 ERFA 库子集
 * (LGPL-3.0, Copyright NumFOCUS Foundation, 源自 IAU SOFA 库)。
 *
 * 版本: v1.1.261115-6
 */
#include "apparent_catalog.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "erfa.h"

/* 编译期确认二进制布局 */
typedef char ac_check_rec[(sizeof(ac_star_record_t) == 80) ? 1 : -1];
typedef char ac_check_hdr[(sizeof(ac_catalog_header_t) == 80) ? 1 : -1];

#define AC_AU_KM   149597870.7        /* 1 AU, km */
#define AC_C_AUD   173.1446326846693  /* 光速, AU/day (与 ERFA_DC 一致) */
#define AC_M_EARTH 3.00348959632e-6   /* 地球质量, 太阳质量单位 */
#define AC_DR2D    57.2957795130823208767981548141051703  /* 弧度→度 */
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* MSVC C 模式无 M_PI */
#endif

/* ------------------------------------------------------------------ */
/* 内部: 观测上下文（对给定时刻+卫星状态只算一次）                      */
/* ------------------------------------------------------------------ */
typedef struct {
    double jd_tt;      /* TT 儒略日 */
    double pmt;        /* J2000 起儒略年数 (TDB≈TT), 供 eraPmpx */
    double pob[3];     /* SSB→观测者, AU */
    double vnat[3];    /* 观测者质心系速度, 自然单位 c=1 */
    double bm1;        /* sqrt(1-|v|²) */
    double ph[3];      /* 太阳→观测者, AU */
    double em;         /* |太阳→观测者|, AU */
    double r_au[3];    /* 地心→卫星, AU */
} ac_context_t;

/* 两单位矢量夹角 (rad) */
static double ac_vsep(double a[3], double b[3])
{
    double axb[3], s, c;
    eraPxp(a, b, axb);
    s = sqrt(eraPdp(axb, axb));
    c = eraPdp(a, b);
    return atan2(s, c);
}

/* ------------------------------------------------------------------ */
/* 时间                                                                */
/* ------------------------------------------------------------------ */
double ac_utc_to_jd(int year, int month, int day, double hour_utc)
{
    double djm0, djm1;
    if (year < 1972 || year > 2100) return NAN;  /* 闰秒表/星历适用域 */
    if (eraCal2jd(year, month, day, &djm0, &djm1) != 0) return NAN;
    return djm0 + djm1 + hour_utc / 24.0;
}

double ac_jd_utc_to_jd_tt(double jd_utc)
{
    int iy, im, id, j;
    double fd, dat;
    j = eraJd2cal(jd_utc, 0.0, &iy, &im, &id, &fd);
    if (j < 0) return NAN;
    j = eraDat(iy, im, id, fd, &dat);
    if (j < 0) return NAN;
    return jd_utc + (dat + 32.184) / 86400.0;
}

/* ------------------------------------------------------------------ */
/* 观测上下文初始化                                                    */
/* ------------------------------------------------------------------ */
static int ac_context_init(ac_context_t *ctx, double jd_utc,
                           const double r_sat_km[3], const double v_sat_kms[3])
{
    double pvh[2][3], pvb[2][3], vv;
    int i, j;

    if (!isfinite(jd_utc)) return AC_E_ARG;
    ctx->jd_tt = ac_jd_utc_to_jd_tt(jd_utc);
    if (!isfinite(ctx->jd_tt)) return AC_E_TIME;

    /* 地球日心/质心星历 (1900–2100) */
    j = eraEpv00(ctx->jd_tt, 0.0, pvh, pvb);
    if (j < 0) return AC_E_TIME;

    ctx->pmt = (ctx->jd_tt - AC_EPOCH_J2000) / 365.25;

    for (i = 0; i < 3; i++) {
        ctx->r_au[i] = r_sat_km[i] / AC_AU_KM;
        ctx->pob[i] = pvb[0][i] + ctx->r_au[i];           /* SSB→观测者, AU */
        ctx->vnat[i] = (pvb[1][i]                          /* AU/day */
                        + v_sat_kms[i] * 86400.0 / AC_AU_KM) / AC_C_AUD;
        ctx->ph[i] = pvh[0][i] + ctx->r_au[i];             /* 太阳→观测者, AU */
    }
    ctx->em = sqrt(eraPdp(ctx->ph, ctx->ph));

    vv = eraPdp(ctx->vnat, ctx->vnat);
    if (vv >= 1.0) return AC_E_ARG;      /* 速度超光速, 输入有误 */
    ctx->bm1 = sqrt(1.0 - vv);

    return AC_OK;
}

/* ------------------------------------------------------------------ */
/* 单星视方向                                                          */
/* ------------------------------------------------------------------ */
static void ac_star_direction(ac_context_t *ctx,
                              const ac_star_record_t *rec,
                              unsigned defl_flags,
                              double u_gcrs[3],
                              double *ra_deg, double *dec_deg,
                              ac_breakdown_t *brk)
{
    double p[3], p0[3], p_after_pm[3], e[3], em, dlim, em2;
    int i;

    /* 1) 自行传播 + 周年视差 → BCRS 坐标方向 */
    eraPmpx(rec->ra_j2000, rec->dec_j2000,
            rec->pmra_j2000, rec->pmdec_j2000,
            rec->parallax_as, rec->rv_kms,
            ctx->pmt, ctx->pob, p);
    for (i = 0; i < 3; i++) p_after_pm[i] = p[i];

    /* 2) 引力光线偏折 */
    if (defl_flags & AC_DEFL_SUN) {
        em = ctx->em;
        em2 = em * em;
        if (em2 < 1.0) em2 = 1.0;
        dlim = 1e-6 / em2;                        /* 与 SOFA iauLdsun 相同 */
        for (i = 0; i < 3; i++) e[i] = ctx->ph[i] / em;
        eraLd(1.0, p, p, e, em, dlim, p);         /* bm=1 太阳质量 */
    }
    if (defl_flags & AC_DEFL_EARTH) {
        em = sqrt(eraPdp(ctx->r_au, ctx->r_au));
        if (em > 0.0) {
            em2 = em * em;
            if (em2 < 1.0) em2 = 1.0;
            dlim = 1e-6 / em2;
            for (i = 0; i < 3; i++) e[i] = ctx->r_au[i] / em;
            eraLd(AC_M_EARTH, p, p, e, em, dlim, p);
        }
    }

    /* 3) 光行差（观测者质心系总速度） → GCRS 视方向 */
    eraAb(p, ctx->vnat, ctx->em, ctx->bm1, u_gcrs);

    /* 赤经赤纬 */
    if (ra_deg && dec_deg) {
        *ra_deg  = eraAnp(atan2(u_gcrs[1], u_gcrs[0])) * AC_DR2D;
        *dec_deg = atan2(u_gcrs[2],
                         sqrt(u_gcrs[0] * u_gcrs[0] + u_gcrs[1] * u_gcrs[1]))
                   * AC_DR2D;
    }

    /* 各项改正量（相对 J2016.0 星表原始方向） */
    if (brk) {
        double cd = cos(rec->dec_j2016_deg * (M_PI / 180.0));
        p0[0] = cos(rec->ra_j2016_deg * (M_PI / 180.0)) * cd;
        p0[1] = sin(rec->ra_j2016_deg * (M_PI / 180.0)) * cd;
        p0[2] = sin(rec->dec_j2016_deg * (M_PI / 180.0));
        brk->pm_parallax_as = ac_vsep(p0, p_after_pm) * (180.0 / M_PI) * 3600.0;
        brk->deflection_mas = ac_vsep(p_after_pm, p)     * (180.0 / M_PI) * 3.6e6;
        brk->aberration_as  = ac_vsep(p, u_gcrs)         * (180.0 / M_PI) * 3600.0;
        brk->total_as       = ac_vsep(p0, u_gcrs)        * (180.0 / M_PI) * 3600.0;
    }
}

/* ------------------------------------------------------------------ */
/* 公开 API                                                            */
/* ------------------------------------------------------------------ */
int ac_apparent_catalog(const ac_catalog_t *cat,
                        double jd_utc,
                        const double r_sat_gcrs_km[3],
                        const double v_sat_gcrs_kms[3],
                        unsigned defl_flags,
                        double *ra_app_deg,
                        double *dec_app_deg,
                        double *unit_gcrs,
                        ac_breakdown_t *brk)
{
    ac_context_t ctx;
    int i, j, n, st;

    if (!cat || !cat->stars || !r_sat_gcrs_km || !v_sat_gcrs_kms)
        return AC_E_ARG;
    if (defl_flags == 0) defl_flags = AC_DEFL_SUN_EARTH;

    st = ac_context_init(&ctx, jd_utc, r_sat_gcrs_km, v_sat_gcrs_kms);
    if (st != AC_OK) return st;

    n = cat->hdr.n_stars;
    for (i = 0; i < n; i++) {
        double u[3], ra = 0.0, dec = 0.0;
        ac_star_direction(&ctx, &cat->stars[i], defl_flags, u,
                          (ra_app_deg || dec_app_deg) ? &ra : NULL,
                          (ra_app_deg || dec_app_deg) ? &dec : NULL,
                          brk ? &brk[i] : NULL);
        if (ra_app_deg)  ra_app_deg[i] = ra;
        if (dec_app_deg) dec_app_deg[i] = dec;
        if (unit_gcrs)
            for (j = 0; j < 3; j++) unit_gcrs[3 * (size_t)i + j] = u[j];
    }
    return AC_OK;
}

int ac_apparent_star(const ac_star_record_t *rec,
                     double jd_utc,
                     const double r_sat_gcrs_km[3],
                     const double v_sat_gcrs_kms[3],
                     unsigned defl_flags,
                     double u_gcrs[3],
                     double *ra_app_deg,
                     double *dec_app_deg,
                     ac_breakdown_t *brk)
{
    ac_context_t ctx;
    int st;

    if (!rec || !r_sat_gcrs_km || !v_sat_gcrs_kms) return AC_E_ARG;
    if (defl_flags == 0) defl_flags = AC_DEFL_SUN_EARTH;

    st = ac_context_init(&ctx, jd_utc, r_sat_gcrs_km, v_sat_gcrs_kms);
    if (st != AC_OK) return st;

    ac_star_direction(&ctx, rec, defl_flags, u_gcrs,
                      ra_app_deg, dec_app_deg, brk);
    return AC_OK;
}

/* ------------------------------------------------------------------ */
/* 星表读写                                                            */
/* ------------------------------------------------------------------ */
int ac_load_catalog(const char *path, ac_catalog_t *cat)
{
    FILE *f;
    size_t n;

    if (!path || !cat) return AC_E_ARG;
    memset(cat, 0, sizeof(*cat));

    f = fopen(path, "rb");
    if (!f) return AC_E_IO;
    if (fread(&cat->hdr, sizeof(cat->hdr), 1, f) != 1) {
        fclose(f);
        return AC_E_IO;
    }
    if (memcmp(cat->hdr.magic, "STDCAT01", 8) != 0) {
        fclose(f);
        return AC_E_FORMAT;
    }
    if (cat->hdr.format_version != 1 ||
        cat->hdr.n_stars <= 0 || cat->hdr.n_stars > 100000000) {
        fclose(f);
        return AC_E_FORMAT;
    }
    n = (size_t)cat->hdr.n_stars;
    cat->stars = (ac_star_record_t *)malloc(n * sizeof(ac_star_record_t));
    if (!cat->stars) {
        fclose(f);
        return AC_E_IO;
    }
    if (fread(cat->stars, sizeof(ac_star_record_t), n, f) != n) {
        free(cat->stars);
        cat->stars = NULL;
        fclose(f);
        return AC_E_IO;
    }
    fclose(f);
    return AC_OK;
}

void ac_free_catalog(ac_catalog_t *cat)
{
    if (cat && cat->stars) {
        free(cat->stars);
        cat->stars = NULL;
    }
}
