/* -*- coding: utf-8 -*-
 * apparent_catalog.h —— 观测时刻 + 卫星状态矢量 → 在轨观测星表（GCRS 视方向）
 *
 * 物理链路（与 SOFA iauATCIQ 一致，不做 bias-岁差-章动旋转，姿态解算自行吸收）:
 *   eraPmpx   自行传播 + 周年视差（含 Roemer 光行时近似）
 *   eraLd     太阳/地球引力光线偏折（可扩展其它行星）
 *   eraAb     光行差（观测者质心系总速度 = 地球公转 + 卫星轨道速度）
 *   地球星历  eraEpv00（1900–2100, 速度精度 42 mm/s → 光行差误差 ~29 μas）
 *
 * 输入时间尺度 UTC（内部经 eraDat 闰秒表转 TT；TT≈TDB，差异 <0.1 μas）。
 * 卫星状态矢量须为 GCRS（地心赤道惯性系, ≈J2000 惯性系）; 单位 km、km/s。
 *
 * 星表二进制格式由 build_std_catalog.py 生成（Gaia DR3, J2016.0 历元源数据）。
 *
 * 版本: v1.0.261115-1
 */
#ifndef APPARENT_CATALOG_H
#define APPARENT_CATALOG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- 常量 ---- */
#define AC_VERSION_STRING   "v1.0.261115-1"
#define AC_EPOCH_GAIA_J2016 2457389.0   /* Gaia DR3 星表历元 J2016.0 (JD) */
#define AC_EPOCH_J2000      2451545.0   /* 星表已传播至 J2000.0 */

/* 错误码 */
#define AC_OK          0
#define AC_E_ARG      (-1)   /* 参数非法 */
#define AC_E_TIME     (-2)   /* 时刻超出 epv00/闰秒表有效范围 (1900–2100) */
#define AC_E_FORMAT   (-3)   /* 星表文件格式不符 */
#define AC_E_IO       (-4)   /* 星表文件读写错误 */

/* ---- 二进制星表: 记录 (80 字节, 与文件布局一一对应) ---- */
typedef struct {
    int64_t source_id;      /* Gaia DR3 source_id */
    double  ra_j2000;       /* ICRS 赤经 @J2000.0, rad（自 J2016.0 严格传播） */
    double  dec_j2000;      /* ICRS 赤纬 @J2000.0, rad */
    double  pmra_j2000;     /* 赤经自行 dRA/dt（不含 cosδ, ERFA 约定）, rad/yr */
    double  pmdec_j2000;    /* 赤纬自行 dDec/dt, rad/yr */
    double  parallax_as;    /* 视差, arcsec（ERFA 约定单位） */
    double  rv_kms;         /* 视向速度, km/s（无测值处为 0） */
    double  g_mag;          /* Gaia G 星等 */
    double  ra_j2016_deg;   /* 原始 Gaia DR3 赤经 @J2016.0, deg（误差分析基准） */
    double  dec_j2016_deg;  /* 原始 Gaia DR3 赤纬 @J2016.0, deg */
} ac_star_record_t;         /* 9*8 + 8 = 80 字节, 无填充 */

/* ---- 二进制星表: 文件头 (80 字节) ---- */
typedef struct {
    char     magic[8];          /* "STDCAT01" */
    int32_t  format_version;    /* 1 */
    int32_t  n_stars;
    double   g_max;             /* 生成时的星等截断 */
    double   jd_epoch_source;   /* 2457389.0 (J2016.0) */
    double   jd_epoch_propagated; /* 2451545.0 (J2000.0) */
    int32_t  reserved[10];
} ac_catalog_header_t;          /* 8+4+4+8+8+8+40 = 80 字节 */

typedef struct {
    ac_catalog_header_t hdr;
    ac_star_record_t   *stars;  /* n_stars 条记录（ac_load_catalog 内 malloc） */
} ac_catalog_t;

/* ---- 各项改正量（弧秒/毫角秒），供误差分析 ---- */
typedef struct {
    double pm_parallax_as;  /* 自行传播+视差 合计改正量, arcsec */
    double deflection_mas;  /* 引力光线偏折, mas */
    double aberration_as;   /* 光行差, arcsec */
    double total_as;        /* 相对 J2016.0 星表方向的总变化, arcsec（矢量合成） */
} ac_breakdown_t;

/* ---- 引力偏折体选择 ---- */
#define AC_DEFL_NONE   0
#define AC_DEFL_SUN    1
#define AC_DEFL_EARTH  2
#define AC_DEFL_SUN_EARTH (AC_DEFL_SUN | AC_DEFL_EARTH)  /* 缺省 */

/* ============================================================
 * 时间工具 (UTC)
 * ============================================================ */
/* 公历日期 + UTC 小时 → 儒略日 (UTC 尺度)。hour 可含小数秒。 */
double ac_utc_to_jd(int year, int month, int day, double hour_utc);

/* UTC 儒略日 → TT 儒略日（经 eraDat 闰秒表）。 */
double ac_jd_utc_to_jd_tt(double jd_utc);

/* ============================================================
 * 星表读写
 * ============================================================ */
/* 读取标准星表二进制文件（build_std_catalog.py 生成）。成功返回 AC_OK。 */
int ac_load_catalog(const char *path, ac_catalog_t *cat);

/* 释放 ac_load_catalog 分配的资源 */
void ac_free_catalog(ac_catalog_t *cat);

/* ============================================================
 * 核心: 在轨观测星表
 * ============================================================ */
/* 整星表计算。观测时刻 jd_utc（儒略日, UTC 尺度），卫星 GCRS 位置 km / 速度 km/s。
 * 四个输出指针均可为 NULL（不需要的不计算）:
 *   ra_app_deg/dec_app_deg : [n] 视赤经/视赤纬, deg, GCRS
 *   unit_gcrs              : [n][3] 视方向单位矢量, 行主序
 *   brk                    : [n] 各项改正量
 * 返回 AC_OK 或负错误码。defl_flags 取 AC_DEFL_* 组合，缺省用 AC_DEFL_SUN_EARTH。
 */
int ac_apparent_catalog(const ac_catalog_t *cat,
                        double jd_utc,
                        const double r_sat_gcrs_km[3],
                        const double v_sat_gcrs_kms[3],
                        unsigned defl_flags,
                        double *ra_app_deg,
                        double *dec_app_deg,
                        double *unit_gcrs,
                        ac_breakdown_t *brk);

/* 单星版本（重复初始化星历，仅适合零星调用） */
int ac_apparent_star(const ac_star_record_t *rec,
                     double jd_utc,
                     const double r_sat_gcrs_km[3],
                     const double v_sat_gcrs_kms[3],
                     unsigned defl_flags,
                     double u_gcrs[3],
                     double *ra_app_deg,
                     double *dec_app_deg,
                     ac_breakdown_t *brk);

#ifdef __cplusplus
}
#endif

#endif /* APPARENT_CATALOG_H */
