# -*- coding: utf-8 -*-
"""
build_std_catalog.py —— Gaia DR3 摘要 → 标准星表（供 C/C++ 与 Python 使用）

输入 CSV 字段: ra, dec, parallax, pmra, pmdec, radial_velocity,
               phot_g_mean_mag, source_id   （Gaia DR3, 历元 J2016.0）

输出:
  1) std_catalog_G{X}.bin —— C/C++ 直接 fread 的定长二进制:
       文件头 80B: "STDCAT01" + format_version + n_stars + g_max
                   + jd_epoch_source + jd_epoch_propagated + reserved[10]
       记录   80B: source_id(i8) + 9×f8
                (ra/dec/pmra/pmdec @J2000.0 严格传播自 J2016.0,
                 parallax_mas, rv_kms, g_mag, ra/dec_j2016_deg)
  2) std_catalog_G{X}.npz —— 同内容 numpy 版（apparent_catalog.py 使用）

处理:
  - G 星等截断（缺省 9.0）
  - 自行/视差缺失(NaN)的行剔除并报告
  - 视向速度缺失置 0（透视加速度影响 μas 级）
  - eraPmsafe 把 (α,δ,μα*,μδ) 从 J2016.0 严格传播到 J2000.0（ERFA 标准星表
    假定 J2000 历元, 运行时 eraPmpx/erfa.pmpx 直接可用）

用法:
  uv run --no-project --with numpy --with pandas --with pyerfa \
      python -I build_std_catalog.py --input gaia_dr3_13.csv --gmax 9.0

版本: v1.0.261115-1
"""
import argparse
import json
import struct
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import erfa

__version__ = "v1.1.261115-2"

JD_J2000 = 2451545.0
JD_GAIA_J2016 = 2457389.0     # Gaia DR3 历元 J2016.0

# 与 c/apparent_catalog.h 的 ac_star_record_t 逐字段一致
REC_DTYPE = np.dtype([
    ("source_id", "<i8"),
    ("ra_j2000", "<f8"), ("dec_j2000", "<f8"),
    ("pmra_j2000", "<f8"), ("pmdec_j2000", "<f8"),
    ("parallax_as", "<f8"), ("rv_kms", "<f8"), ("g_mag", "<f8"),
    ("ra_j2016_deg", "<f8"), ("dec_j2016_deg", "<f8"),
])
assert REC_DTYPE.itemsize == 80

HEADER_FMT = "<8s2i3d10i"       # 与 ac_catalog_header_t 逐字段一致
assert struct.calcsize(HEADER_FMT) == 80


def main():
    ap = argparse.ArgumentParser(description="Gaia DR3 → 标准星表 (.bin/.npz)")
    ap.add_argument("--input", default="gaia_dr3_13.csv")
    ap.add_argument("--gmax", type=float, default=9.0, help="G 星等截断 (缺省 9.0)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    gtag = ("G" + str(args.gmax)).replace(".", "p")
    bin_path = f"{args.outdir}/std_catalog_{gtag}.bin"
    npz_path = f"{args.outdir}/std_catalog_{gtag}.npz"

    print(f"build_std_catalog {__version__}")
    print(f"读取 {args.input} ...")
    df = pd.read_csv(args.input)

    n0 = len(df)
    g = df["phot_g_mean_mag"].to_numpy(float)
    keep = np.isfinite(g) & (g <= args.gmax)
    for col in ("ra", "dec", "parallax", "pmra", "pmdec"):
        keep &= np.isfinite(df[col].to_numpy(float))
    n1 = int(keep.sum())
    print(f"总行数 {n0}, G<={args.gmax} 且自行/视差完备: {n1} 颗"
          f"（剔除 {n0 - n1} 行）")

    ra = np.deg2rad(df["ra"].to_numpy(float)[keep])
    dec = np.deg2rad(df["dec"].to_numpy(float)[keep])
    px = df["parallax"].to_numpy(float)[keep]                 # mas
    pmra = df["pmra"].to_numpy(float)[keep]                   # mas/yr (含 cosδ)
    pmdec = df["pmdec"].to_numpy(float)[keep]                 # mas/yr
    rv = df["radial_velocity"].to_numpy(float)[keep]
    rv[~np.isfinite(rv)] = 0.0                                # km/s
    gm = g[keep]
    sid = df["source_id"].to_numpy(np.int64)[keep]

    # J2016.0 → J2000.0 严格传播（含视差/视向速度的历元变换）
    # 注意 ERFA 约定: 自行 rad/yr 且赤经方向为 dRA/dt（Gaia 的 pmra 是
    # μα* 含 cosδ, 必须除以 cosδ）; 视差 arcsec; 视向速度 km/s
    mas2rad = np.pi / (180.0 * 3600.0 * 1000.0)
    px_as = px / 1000.0
    cosd = np.cos(dec)
    pmra_dradt = np.divide(pmra, cosd, out=np.zeros_like(pmra),
                           where=np.abs(cosd) > 1e-10)  # mas/yr, dRA/dt
    r2, d2, pr2, pd2, px2, rv2 = erfa.pmsafe(
        ra, dec, pmra_dradt * mas2rad, pmdec * mas2rad, px_as, rv,
        JD_GAIA_J2016, 0.0, JD_J2000, 0.0)

    rec = np.zeros(n1, dtype=REC_DTYPE)
    rec["source_id"] = sid
    rec["ra_j2000"] = r2
    rec["dec_j2000"] = d2
    rec["pmra_j2000"] = pr2
    rec["pmdec_j2000"] = pd2
    rec["parallax_as"] = px_as
    rec["rv_kms"] = rv
    rec["g_mag"] = gm
    rec["ra_j2016_deg"] = np.rad2deg(ra)
    rec["dec_j2016_deg"] = np.rad2deg(dec)

    header = struct.pack(
        HEADER_FMT, b"STDCAT01", 1, n1, float(args.gmax),
        JD_GAIA_J2016, JD_J2000, *([0] * 10))

    with open(bin_path, "wb") as f:
        f.write(header)
        f.write(rec.tobytes())
    np.savez(npz_path, **{name: rec[name] for name in REC_DTYPE.names})
    print(f"已生成 {bin_path} ({80 + 80 * n1} 字节)")
    print(f"已生成 {npz_path}")

    # 审计：极端恒星必须在表（此前流程曾出现极端星丢失）
    mu = np.hypot(pmra, pmdec)
    top = np.argsort(mu)[-5:][::-1]
    print("\n最大自行 5 星 (mas/yr, 应含 Barnard ~10393):")
    for i in top:
        print(f"  source_id={sid[i]:>19}  mu={mu[i]:9.1f}  pi={px[i]:7.2f} mas"
              f"  G={gm[i]:5.2f}  ra={np.rad2deg(ra[i]):8.3f} dec={np.rad2deg(dec[i]):8.3f}")
    topp = np.argsort(px)[-5:][::-1]
    print("最大视差 5 星 (mas, 应含 Proxima ~768):")
    for i in topp:
        print(f"  source_id={sid[i]:>19}  mu={mu[i]:9.1f}  pi={px[i]:7.2f} mas"
              f"  G={gm[i]:5.2f}  ra={np.rad2deg(ra[i]):8.3f} dec={np.rad2deg(dec[i]):8.3f}")

    meta = dict(build_utc=datetime.now(timezone.utc).isoformat(),
                version=__version__, input=args.input, g_max=args.gmax,
                n_stars=n1, jd_epoch_source=JD_GAIA_J2016,
                jd_epoch_propagated=JD_J2000)
    print("\n" + json.dumps(meta, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
