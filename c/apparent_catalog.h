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
#define AC_VERSION_STRING   "v1.1.261115-4"
#define AC_EPOCH_GAIA_J2016 2457389.0   /* Gaia DR3 星表历元 J2016.0 (JD) */
#define AC_EPOCH_J2000      2451545.0   /* 星表已传播至 J2000.0 */

/* 错误码 */
#define AC_OK          0
#define AC_E_ARG      (-1)   /* 参数非法 */
#define AC_E_TIME     (-2)   /* 时刻超出 epv00/闰秒表有效范围 (1900–2100) */
#define AC_E_FORMAT   (-3)   /* 星表文件格式不符 */
#define AC_E_IO       (-4)   /* 星表文件读写错误 */
#define AC_E_EOP      (-5)   /* EOP 文件读取/解析错误 */
#define AC_W_EOP_RANGE (+1)  /* 警告: 请求时刻超出 EOP 数据范围, 已用端点值 */

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

/* ============================================================
 * IERS EOP C04（IAU2000）接入与 ITRS/ECEF 支持
 * ============================================================ */

/* 单日 EOP 记录（每日 0h UTC） */
typedef struct {
    double mjd_utc;   /* 简化儒略日 (UTC) */
    double xp_as;     /* 极移 X, arcsec */
    double yp_as;     /* 极移 Y, arcsec */
    double dut1_s;    /* UT1-UTC, s */
    double lod_s;     /* 日长, s */
    double dx_as;     /* 天极偏差 dX, arcsec */
    double dy_as;     /* 天极偏差 dY, arcsec */
} ac_eop_record_t;

typedef struct {
    int n;
    ac_eop_record_t *rec;   /* 按 mjd 升序 */
    double mjd_first, mjd_last;
} ac_eop_t;

/* 解析 IERS EOP C04 (IAU2000) 文本文件（update_eop.py 下载的
 * eopc04_IAU2000.62-now，格式: 日期 MJD x y UT1-UTC LOD dX dY ...） */
int ac_eop_load(const char *path, ac_eop_t *eop);

void ac_eop_free(ac_eop_t *eop);

/* 按时刻线性插值 EOP（mjd_utc = JD_UTC − 2400000.5）。
 * 本项目固定使用长期产品 C04（不使用 IERS 快速产品）。C04 约月更、
 * 末端滞后约 30 天；超出数据范围时采用外推并返回 AC_W_EOP_RANGE 警告：
 *   DUT1: dut1(t0+Δ) = dut1(t0) − LOD(t0)·Δ （解析关系）
 *   极移: 末端 30 天线性趋势外推
 *   dX/dY: 毫角秒级, 端点保持
 * 外推残差量级见 update_eop.py 的留出回测（EOP_report.md）。 */
int ac_eop_get(const ac_eop_t *eop, double mjd_utc,
               double *xp_as, double *yp_as, double *dut1_s,
               double *dx_as, double *dy_as);

/* ITRS(地固/ECEF) → GCRS(惯性) 状态转换（位置 km、速度 km/s）。
 * 使用 EOP 极移 + DUT1 + IAU 2006/2000A 无旋转原点岁差章动
 * （erfaC2t06a 链路）；速度含地球自转修正 v_g = Rᵀ(v_i + ω×r_i)，
 * 忽略岁差章动变率（<4 μm/s）。 */
int ac_itrs_to_gcrs(const ac_eop_t *eop, double jd_utc,
                    double r_itrs[3], double v_itrs[3],
                    double r_gcrs[3], double v_gcrs[3]);

/* GCRS → ITRS 逆变换 */
int ac_gcrs_to_itrs(const ac_eop_t *eop, double jd_utc,
                    double r_gcrs[3], double v_gcrs[3],
                    double r_itrs[3], double v_itrs[3]);

/* 一步到位: ITRS 卫星状态 → 在轨观测星表（内部先转 GCRS） */
int ac_apparent_catalog_itrs(const ac_catalog_t *cat,
                             const ac_eop_t *eop,
                             double jd_utc,
                             double r_sat_itrs_km[3],
                             double v_sat_itrs_kms[3],
                             unsigned defl_flags,
                             double *ra_app_deg,
                             double *dec_app_deg,
                             double *unit_gcrs,
                             ac_breakdown_t *brk);

#ifdef __cplusplus
}
#endif

#endif /* APPARENT_CATALOG_H */
