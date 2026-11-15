# -*- coding: utf-8 -*-
"""
validate_c_vs_python.py —— C 实现 (c/apparent_catalog.c) vs Python 参考实现
(apparent_catalog.py, ERFA) 逐星比对。

前提: 已运行 c/build/.../test_apparent_catalog 生成 c_sample_output.csv
      （默认观测条件: 2026-11-15 12:00 UTC, 550 km LEO 示例状态）。

用法:
  uv run --no-project --with numpy --with pyerfa \
      python -I validate_c_vs_python.py std_catalog_G9p0.npz c_sample_output.csv \
      --utc "2026-11-15 12:00:00" \
      --pos 6928.137 0 0 --vel 0 -1.0065 7.5219

版本: v1.0.261115-1
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import apparent_catalog as ac

__version__ = "v1.0.261115-1"

D2R = np.pi / 180.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("csv")
    ap.add_argument("--utc", default="2026-11-15 12:00:00")
    ap.add_argument("--pos", nargs=3, type=float, default=[6928.137, 0.0, 0.0])
    ap.add_argument("--vel", nargs=3, type=float, default=[0.0, -1.0065, 7.5219])
    args = ap.parse_args()

    cat = ac.load_catalog(args.npz)
    out = ac.apparent_catalog(cat, args.utc, args.pos, args.vel,
                              bodies=("SUN", "EARTH"), include_breakdown=True)

    by_index = {}
    with open(args.csv, newline="") as f:
        for row in csv.DictReader(f):
            by_index[int(row["index"])] = row

    n_cmp = 0
    d_max = {"ra": 0.0, "dec": 0.0, "total": 0.0, "pmpx": 0.0,
             "defl": 0.0, "abbr": 0.0}
    for i, row in sorted(by_index.items()):
        if cat["source_id"][i] != int(row["source_id"]):
            raise SystemExit(f"第 {i} 行 source_id 不一致: "
                             f"{cat['source_id'][i]} vs {row['source_id']}")
        dra = abs(out["ra_app_deg"][i] - float(row["ra_app_deg"])) * 3600.0
        dra *= np.cos(out["dec_app_deg"][i] * D2R)
        ddec = abs(out["dec_app_deg"][i] - float(row["dec_app_deg"])) * 3600.0
        d_max["ra"] = max(d_max["ra"], dra)
        d_max["dec"] = max(d_max["dec"], ddec)
        d_max["total"] = max(d_max["total"],
                             abs(out["total_as"][i] - float(row["total_as"])))
        d_max["pmpx"] = max(d_max["pmpx"],
                            abs(out["pm_parallax_as"][i] - float(row["pm_px_as"])))
        d_max["defl"] = max(d_max["defl"],
                            abs(out["deflection_mas"][i] - float(row["defl_mas"])))
        d_max["abbr"] = max(d_max["abbr"],
                            abs(out["aberration_as"][i] - float(row["abbr_as"])))
        n_cmp += 1

    print(f"比对 {n_cmp} 星 (C vs Python/ERFA):")
    for k, v in d_max.items():
        unit = "mas" if k == "defl" else "arcsec"
        scale = 1000.0 if k == "defl" else 1.0
        print(f"  max|Δ{k:>5}| = {v * scale:10.4f} {unit}")
    tol = 1e-3  # arcsec
    ok = max(d_max["ra"], d_max["dec"]) < tol
    print("结论:", "PASS (一致)" if ok else "FAIL — 差异超 1 mas，需排查！")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
