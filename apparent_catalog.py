# -*- coding: utf-8 -*-
"""
apparent_catalog.py —— 观测时刻 + 卫星状态矢量 → 在轨观测星表（GCRS 视方向）

功能
----
输入观测时刻 (UTC) 与卫星状态矢量 (位置 km、速度 km/s，GCRS 惯性系或 ITRS 地固系)，
对标准星表 (Gaia DR3, 由 build_std_catalog.py 生成) 逐星计算观测瞬间的视方向：

    pmpx (自行传播 + 周年视差)  →  ld (太阳/地球/行星引力光线偏折)  →  ab (光行差)

链路与 SOFA iauATCIQ 完全一致（self_test() 逐微角秒校验），但**不做 bias-岁差-章动
旋转**：星敏感器姿态解算会吸收这种全天空共有的旋转，逐星修正项（自行/视差/偏折/
光行差）才必须逐星施加。

输出 GCRS 单位矢量与视赤经/视赤纬 (度)。TDB≈TT（差 <1.7 ms，对光行差影响 <0.1 μas）。

用法示例
--------
>>> import numpy as np, apparent_catalog as ac
>>> cat = ac.load_catalog("std_catalog_G9p0.npz")
>>> out = ac.apparent_catalog(cat, "2026-11-15 12:00:00",
...                           r_sat=[6928.137, 0, 0],      # km, GCRS
...                           v_sat=[0, -1.0065, 7.5219],  # km/s, GCRS
...                           include_breakdown=True)
>>> out["ra_app_deg"], out["dec_app_deg"], out["total_as"]

依赖: numpy, pyerfa   （uv run --no-project --with numpy --with pyerfa ...）
"""
from __future__ import annotations

from datetime import datetime, timezone
import argparse
import time

import numpy as np
import erfa

__version__ = "v1.1.261115-4"

DJ00 = 2451545.0            # J2000.0 (TDB)
GAIA_EPOCH_JD = 2457389.0   # Gaia DR3 星表历元 J2016.0
DJY = 365.25
AU_KM = erfa.DAU / 1000.0   # 1 AU = 149597870.7 km
C_AUD = erfa.DC             # 光速, AU/day

# erfa.plan94 的行星编号 与 引力体质量 (太阳质量单位, IAU 2015 标称值)
_BODY_NP = {"MERCURY": 1, "VENUS": 2, "EARTH": 3, "MARS": 4,
            "JUPITER": 5, "SATURN": 6, "URANUS": 7, "NEPTUNE": 8}
_BODY_MASS = {"SUN": 1.0, "MERCURY": 1.6601e-7, "VENUS": 2.4478e-6,
              "EARTH": 3.00348959632e-6, "MARS": 3.2271e-7,
              "JUPITER": 9.5479e-4, "SATURN": 2.858e-7,
              "URANUS": 4.3662e-8, "NEPTUNE": 5.1506e-8}


# ----------------------------------------------------------------------
# 时间与坐标
# ----------------------------------------------------------------------
def _utc_two_part(utc):
    """UTC → 两部儒略日 (UTC 尺度)。接受 datetime 或 (d1, d2)。"""
    if isinstance(utc, datetime):
        dt = utc if utc.tzinfo is None else utc.astimezone(timezone.utc)
        d1, d2 = erfa.dtf2d("UTC", dt.year, dt.month, dt.day,
                            dt.hour, dt.minute, dt.second + dt.microsecond * 1e-6)
        return d1, d2
    if isinstance(utc, (tuple, list, np.ndarray)) and len(utc) == 2:
        return float(utc[0]), float(utc[1])
    raise TypeError("utc 需为 datetime 或 (d1,d2) 两部儒略日(UTC)")


def utc_to_tt(utc):
    """UTC → TT 两部儒略日（TAI 中转，自动处理闰秒）。"""
    d1, d2 = _utc_two_part(utc)
    return erfa.taitt(*erfa.utctai(d1, d2))


def parse_utc(s):
    """解析 'YYYY-MM-DD HH:MM:SS[.fff]' (UTC) → datetime"""
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S.%f",
                "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"无法解析时刻字符串: {s!r}")


def itrs_to_gcrs(utc, r_itrs, v_itrs, dut1=0.0, xp=0.0, yp=0.0):
    """ITRS(地固) 位置 km / 速度 km/s → GCRS。

    dut1 = UT1-UTC (秒)。缺省 0 意味着旋转误差最大 ~0.9s×ω ≈ 13″，
    精确使用时必须给真实 DUT1（IERS 公报）；极移 xp/yp 缺省 0（<0.3″ 影响）。
    """
    d1, d2 = _utc_two_part(utc)
    ut1 = erfa.utcut1(d1, d2, dut1)
    tt = erfa.taitt(*erfa.utctai(d1, d2))
    rc2t = erfa.c2t06a(tt[0], tt[1], ut1[0], ut1[1], xp, yp)  # 天球→地固
    r = rc2t.T @ np.asarray(r_itrs, float)
    v = rc2t.T @ np.asarray(v_itrs, float)
    w = np.array([0.0, 0.0, erfa.W])                          # 地球自转角速度 rad/s
    return r, v + np.cross(w, r)


# ----------------------------------------------------------------------
# IERS EOP C04 (IAU2000)
# ----------------------------------------------------------------------
def load_eop_c04(path):
    """解析 IERS EOP C04 (IAU2000) 文本文件（update_eop.py 下载维护）。

    返回 dict: mjd_utc, xp_as, yp_as, dut1_s, lod_s, dx_as, dy_as (升序数组)。
    """
    rows = []
    with open(path, "r", encoding="ascii", errors="replace") as f:
        for line in f:
            t = line.split()
            if len(t) < 10:
                continue
            try:
                iy, im, id_ = int(t[0]), int(t[1]), int(t[2])
                vals = [float(x) for x in t[3:10]]
            except ValueError:
                continue
            mjd = vals[0]
            if not (1962 <= iy <= 2100 and 30000.0 < mjd < 90000.0):
                continue
            rows.append((mjd, *vals[1:]))
    if len(rows) < 2:
        raise ValueError(f"EOP 文件解析失败: {path}")
    arr = np.array(rows)  # mjd, xp, yp, dut1, lod, dx, dy
    return {"mjd_utc": arr[:, 0], "xp_as": arr[:, 1], "yp_as": arr[:, 2],
            "dut1_s": arr[:, 3], "lod_s": arr[:, 4],
            "dx_as": arr[:, 5], "dy_as": arr[:, 6]}


def eop_get(eop, mjd_utc):
    """线性插值取 EOP（与 C 版 ac_eop_get 一致, 长期产品外推策略）。

    末端越界时: DUT1 按 dut1−LOD·Δ 外推（解析关系）, 极移按末端 30 天
    线性趋势外推, dX/dY 端点保持; 首端越界端点保持。本项目固定使用
    长期产品 C04, 不使用 IERS 快速产品。
    """
    m = eop["mjd_utc"]
    if mjd_utc <= m[0]:
        return (eop["xp_as"][0], eop["yp_as"][0], eop["dut1_s"][0],
                eop["dx_as"][0], eop["dy_as"][0])
    if mjd_utc >= m[-1]:
        d = mjd_utc - m[-1]
        dut1 = eop["dut1_s"][-1] - eop["lod_s"][-1] * d
        span = 30
        m0, m1 = max(0, len(m) - span), len(m) - 1
        x = m[m0:m1 + 1] - m[m1]
        out = []
        for k in ("xp_as", "yp_as"):
            b, a = np.polyfit(x, eop[k][m0:m1 + 1], 1)
            out.append(a + b * d)              # a = 拟合在末端 x=0 处的值
        return (out[0], out[1], dut1, eop["dx_as"][-1], eop["dy_as"][-1])
    i = int(np.clip(np.searchsorted(m, mjd_utc), 1, len(m) - 1))
    t = (mjd_utc - m[i - 1]) / (m[i] - m[i - 1])
    out = []
    for k in ("xp_as", "yp_as", "dut1_s", "dx_as", "dy_as"):
        v = eop[k]
        out.append(v[i - 1] + t * (v[i] - v[i - 1]))
    return tuple(out)


def itrs_to_gcrs_eop(eop, utc, r_itrs, v_itrs):
    """用 IERS EOP (极移+DUT1) 做 ITRS→GCRS 状态转换（与 C 版一致）。

    返回 (r_gcrs_km, v_gcrs_kms)。v_g = Rᵀ(v_i + ω×r_i)。
    """
    d1, d2 = _utc_two_part(utc)
    tt = erfa.taitt(*erfa.utctai(d1, d2))
    jd_utc = d1 + d2
    mjd = jd_utc - 2400000.5
    xp, yp, dut1, dx, dy = eop_get(eop, mjd)
    ut1 = jd_utc + dut1 / 86400.0
    rc2t = erfa.c2t06a(tt[0], tt[1], ut1, 0.0,
                       xp * erfa.DAS2R, yp * erfa.DAS2R)
    r = np.asarray(r_itrs, float)
    v = np.asarray(v_itrs, float)
    w = np.array([0.0, 0.0, erfa.W])
    return rc2t.T @ r, rc2t.T @ (v + np.cross(w, r))


# ----------------------------------------------------------------------
# 观测者 (BCRS)
# ----------------------------------------------------------------------
def observer_bcrs(tt1, tt2, r_sat_eci, v_sat_eci):
    """观测者日心系/太阳系质心系状态。

    地球质心状态由 erfa.epv00 给出 (BCRS, AU, AU/day)，叠加卫星地心 GCRS 状态。
    GCRS→BCRS 的相对论尺度差 ~1e-10，对角秒/毫角秒级应用可忽略。
    返回 (pob, vob, ph_obs): 质心系位置 AU、质心系速度 AU/day、日心位置 AU。
    """
    pvh, pvb = erfa.epv00(tt1, tt2)
    r_au = np.asarray(r_sat_eci, float) / AU_KM
    v_au = np.asarray(v_sat_eci, float) * 86400.0 / AU_KM
    return pvb[0] + r_au, pvb[1] + v_au, pvh[0] + r_au


def apply_light_deflection(p, tt1, tt2, ph_obs, r_sat_eci, bodies):
    """逐体施加引力光线偏折 (erfa.ld, 首参为引力体太阳质量)。

    dlim 采用 SOFA iauLdsun 同款限制器 1e-6/max(em²,1)，与 atciq 链路一致。
    太阳在 1 AU 处最大 ~8.1 mas、临日 1.75″；地球从 LEO 看 ≤0.5 mas；
    木星临边 16 mas。多体按离观测者距离降序依次施加（ERFA ld 注 7）。
    """
    r_au = np.asarray(r_sat_eci, float) / AU_KM
    entries = []
    for body in bodies:
        b = str(body).upper()
        if b not in _BODY_MASS:
            raise ValueError(f"未知引力体 {body!r}，可选: {sorted(_BODY_MASS)}")
        if b == "SUN":
            e = ph_obs
        elif b == "EARTH":
            e = r_au
        else:
            pv_planet = erfa.plan94(tt1, tt2, _BODY_NP[b])[0]  # 日心行星位置 AU
            e = ph_obs - pv_planet                              # 行星→观测者
        em = float(np.linalg.norm(e))
        entries.append((em, _BODY_MASS[b], e / em))
    for em, mass, ehat in sorted(entries, key=lambda t: -t[0]):
        dlim = 1e-6 / max(em * em, 1.0)
        p = erfa.ld(mass, p, p, ehat, em, dlim)
    return p


# ----------------------------------------------------------------------
# 主函数
# ----------------------------------------------------------------------
def apparent_catalog(catalog, utc, r_sat, v_sat, frame="GCRS",
                     bodies=("SUN", "EARTH"), include_breakdown=False,
                     dut1=0.0, xp=0.0, yp=0.0):
    """观测时刻 + 卫星状态矢量 → 星表视方向 (GCRS)。

    参数
    ----
    catalog : dict
        load_catalog() 载入的标准星表（Gaia DR3, 已传播到 J2000.0）。
    utc : datetime | (d1, d2) | 'YYYY-MM-DD HH:MM:SS'
        观测时刻 (UTC)。
    r_sat, v_sat : (3,) array_like
        卫星位置 km / 速度 km/s。frame='GCRS'（缺省, 地心赤道惯性系）或
        'ITRS'（地固系，需 dut1 才精确）。
    bodies : tuple
        光线偏折包含的引力体，缺省 ('SUN','EARTH')。
    include_breakdown : bool
        True 时返回各项改正量大小（角秒），用于误差分析。

    返回
    ----
    dict: ra_app_deg, dec_app_deg, unit_vec_gcrs (N,3), g_mag, source_id,
          utc_tt_jd;   include_breakdown=True 时另含
          pm_parallax_as / deflection_mas / aberration_as / total_as
    """
    if isinstance(utc, str):
        utc = parse_utc(utc)
    tt1, tt2 = utc_to_tt(utc)

    r_sat = np.asarray(r_sat, float).reshape(3)
    v_sat = np.asarray(v_sat, float).reshape(3)
    if frame.upper() == "ITRS":
        r_sat, v_sat = itrs_to_gcrs(utc, r_sat, v_sat, dut1, xp, yp)
    elif frame.upper() != "GCRS":
        raise ValueError("frame 只支持 'GCRS' 或 'ITRS'")

    pob, vob, ph_obs = observer_bcrs(tt1, tt2, r_sat, v_sat)
    vnat = vob / C_AUD                       # 观测者质心系速度 (自然单位 c=1)
    bm1 = float(np.sqrt(1.0 - float(vnat @ vnat)))
    pmt = (tt1 + tt2 - DJ00) / DJY           # 自 J2000 起的儒略年数 (TDB≈TT)

    # 1) 自行传播 + 周年视差（含光行时近似，SOFA pmpx; 视差单位 arcsec）
    p_pm_px = erfa.pmpx(catalog["ra_j2000"], catalog["dec_j2000"],
                        catalog["pmra_j2000"], catalog["pmdec_j2000"],
                        catalog["parallax_as"], catalog["rv_kms"], pmt, pob)
    # 2) 引力光线偏折
    p_defl = apply_light_deflection(p_pm_px, tt1, tt2, ph_obs, r_sat, bodies)
    # 3) 光行差（质心系总速度 = 地球公转 + 卫星轨道速度）
    em_sun = float(np.linalg.norm(ph_obs))
    p_app = erfa.ab(p_defl, vnat, em_sun, bm1)

    c1, c2 = erfa.c2s(p_app)
    out = {"utc_tt_jd": tt1 + tt2,
           "ra_app_deg": np.degrees(erfa.anp(c1)),
           "dec_app_deg": np.degrees(c2),
           "unit_vec_gcrs": p_app,
           "g_mag": catalog["g_mag"],
           "source_id": catalog["source_id"]}
    if include_breakdown:
        p0 = erfa.s2c(np.deg2rad(catalog["ra_j2016_deg"]),
                      np.deg2rad(catalog["dec_j2016_deg"]))
        out["pm_parallax_as"] = np.degrees(erfa.sepp(p0, p_pm_px)) * 3600.0
        out["deflection_mas"] = np.degrees(erfa.sepp(p_pm_px, p_defl)) * 3.6e6
        out["aberration_as"] = np.degrees(erfa.sepp(p_defl, p_app)) * 3600.0
        out["total_as"] = np.degrees(erfa.sepp(p0, p_app)) * 3600.0
    return out


def load_catalog(npz_path):
    """载入 build_std_catalog.py 生成的标准星表 npz。"""
    with np.load(npz_path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files if k != "meta"}


# ----------------------------------------------------------------------
# 自检：与 SOFA apci13+atciq (地心观测者, 仅太阳偏折) 逐星比对
# ----------------------------------------------------------------------
def self_test(utc=None, n=200, tol_uas=1.0):
    """对随机星做地心观测者的视位置比对。

    SOFA atciq 输出 CIRS 坐标 = 本模块 GCRS 输出左乘 bias-岁差-章动矩阵。
    两者应逐星一致到 <1 μas，证明分步链路 == SOFA 标准链路。
    """
    utc = utc or datetime(2026, 11, 15, 12, 0, 0)
    rng = np.random.default_rng(42)
    ra = rng.uniform(0.0, 2.0 * np.pi, n)
    dec = np.arcsin(rng.uniform(-1.0, 1.0, n))
    pr = rng.normal(0.0, 200.0, n) / 3.6e6   # rad/yr
    pd = rng.normal(0.0, 200.0, n) / 3.6e6
    px = rng.uniform(0.0001, 0.768, n)       # arcsec (ERFA 约定)
    rv = rng.normal(0.0, 30.0, n)            # km/s
    cat = {"ra_j2000": ra, "dec_j2000": dec, "pmra_j2000": pr,
           "pmdec_j2000": pd, "parallax_as": px, "rv_kms": rv,
           "g_mag": np.zeros(n), "source_id": np.zeros(n, dtype=np.int64)}

    mine = apparent_catalog(cat, utc, np.zeros(3), np.zeros(3), bodies=("SUN",))

    tt = utc_to_tt(utc)
    astrom, _ = erfa.apci13(tt[0], tt[1])
    ri, di = erfa.atciq(ra, dec, pr, pd, px, rv, astrom)
    p_cirs = erfa.s2c(ri, di)                       # CIRS
    p_expect = p_cirs @ astrom["bpn"]               # 行矢量右乘 = bpn^T·p, 逆旋转回 GCRS
    ang = np.degrees(erfa.sepp(p_expect, mine["unit_vec_gcrs"])) * 3.6e6  # μas
    ok = float(ang.max()) < tol_uas
    print(f"[self_test] vs SOFA atciq: max|Δ| = {ang.max():.3f} μas "
          f"(median {np.median(ang):.3f}) -> {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------
# CLI 演示
# ----------------------------------------------------------------------
def _stats(x, unit):
    return (f"mean={np.mean(x):10.4f} median={np.median(x):10.4f} "
            f"P99={np.percentile(x, 99):10.4f} max={np.max(x):10.4f} {unit}")


def main():
    ap = argparse.ArgumentParser(
        description="观测时刻 + 卫星状态矢量 → 在轨观测星表 (GCRS 视方向)")
    ap.add_argument("--catalog", default="std_catalog_G9p0.npz", help="标准星表 npz")
    ap.add_argument("--utc", default="2026-11-15 12:00:00", help="UTC 观测时刻")
    ap.add_argument("--pos", nargs=3, type=float,
                    default=[6928.137, 0.0, 0.0],
                    help="卫星位置 km (缺省: 550 km 圆轨道示例)")
    ap.add_argument("--vel", nargs=3, type=float,
                    default=[0.0, -1.0065, 7.5219],
                    help="卫星速度 km/s")
    ap.add_argument("--frame", default="GCRS", choices=["GCRS", "ITRS"])
    ap.add_argument("--dut1", type=float, default=0.0, help="UT1-UTC 秒 (ITRS 时需要)")
    ap.add_argument("--bodies", nargs="+", default=["SUN", "EARTH"],
                    help="光线偏折引力体, 如 SUN EARTH JUPITER")
    ap.add_argument("--top", type=int, default=12, help="打印最亮 N 颗星的明细")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    print(f"apparent_catalog {__version__}")
    self_test()

    cat = load_catalog(args.catalog)
    print(f"\n星表: {args.catalog}  N={len(cat['g_mag'])} 颗 (G≤9, Gaia DR3)")
    r = np.array(args.pos)
    v = np.array(args.vel)
    print(f"观测时刻 UTC : {args.utc}")
    print(f"卫星位置     : [{r[0]:.3f}, {r[1]:.3f}, {r[2]:.3f}] km ({args.frame})")
    print(f"卫星速度     : [{v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f}] km/s, "
          f"|v|={np.linalg.norm(v):.4f} km/s")

    t0 = time.perf_counter()
    out = apparent_catalog(cat, args.utc, r, v, frame=args.frame,
                           bodies=tuple(args.bodies), include_breakdown=True)
    dt = time.perf_counter() - t0
    print(f"计算耗时     : {dt * 1000:.1f} ms  ({len(cat['g_mag'])} 星)")

    idx = np.argsort(out["g_mag"])[: args.top]
    cat16_ra = np.deg2rad(cat["ra_j2016_deg"][idx])
    cat16_dec = np.deg2rad(cat["dec_j2016_deg"][idx])
    print(f"\n最亮 {args.top} 颗（J2016 星表方向 vs 在轨视方向）:")
    print(f"{'source_id':>19} {'G':>5} {'ra_2016':>12} {'dec_2016':>11} "
          f"{'ra_app':>12} {'dec_app':>11} {'总改正\"':>9} {'自行+视差\"':>10} "
          f"{'偏折mas':>9} {'光行差\"':>9}")
    for k, i in enumerate(idx):
        ra16, d16 = erfa.anp(cat16_ra[k]), cat16_dec[k]
        print(f"{out['source_id'][i]:>19d} {out['g_mag'][i]:>5.2f} "
              f"{np.degrees(ra16):>12.6f} {np.degrees(d16):>11.6f} "
              f"{out['ra_app_deg'][i]:>12.6f} {out['dec_app_deg'][i]:>11.6f} "
              f"{out['total_as'][i]:>9.3f} {out['pm_parallax_as'][i]:>10.3f} "
              f"{out['deflection_mas'][i]:>9.2f} {out['aberration_as'][i]:>9.3f}")

    print(f"\n全星表改正量统计 ({len(cat['g_mag'])} 颗, 与误差分析表同口径):")
    print("  自行+视差  : " + _stats(out["pm_parallax_as"], "arcsec"))
    print("  引力偏折   : " + _stats(out["deflection_mas"], "mas"))
    print("  光行差     : " + _stats(out["aberration_as"], "arcsec"))
    print("  总改正     : " + _stats(out["total_as"], "arcsec"))


if __name__ == "__main__":
    main()
