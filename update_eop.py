# -*- coding: utf-8 -*-
"""
update_eop.py —— IERS EOP C04 接入流水线：下载 → 更新 → 编译 → 测试 → Markdown 报告

数据源（按序尝试）:
  1. https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now   (IERS EOP 产品中心)
  2. https://hpiers.obspm.fr/iers/eop/eopc04/archive_11_12_2023/eopc04_IAU2000.62-now
     （历史存档, 仅应急, 数据较旧）
  也可 --input 指定本地已下载文件跳过网络。

用法（在本项目根目录）:
  uv run --no-project --with numpy --with pyerfa \
      python -I update_eop.py [--skip-download] [--skip-build] [--skip-test]

产物:
  eop/eopc04_IAU2000.62-now   EOP 数据（原子替换, 旧版备份 .prev）
  EOP_report.md               全流程报告

版本: v1.1.261115-4
"""
import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

import numpy as np
import erfa

__version__ = "v1.1.261115-4"

ROOT = os.path.dirname(os.path.abspath(__file__))
EOP_DIR = os.path.join(ROOT, "eop")
EOP_FILE = os.path.join(EOP_DIR, "eopc04_IAU2000.62-now")

SOURCES = [
    ("IERS EOP PC (hpiers, 主源)",
     "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now"),
    ("IERS EOP PC (hpiers, 2023-11 存档, 应急)",
     "https://hpiers.obspm.fr/iers/eop/eopc04/archive_11_12_2023/eopc04_IAU2000.62-now"),
]

# 快速产品源（精度优于长期产品外推, 详见报告对比节; 可达时自动启用）
RAPID_SOURCES = [
    ("USNO IERS RS/PC finals.daily (快速, ~1 日滞后)",
     "https://maia.usno.navy.mil/ser7/finals.daily"),
]
RAPID_FILE = os.path.join(EOP_DIR, "finals.daily")

R = {}   # 报告数据收集


def now_utc():
    return datetime.now(timezone.utc)


def mjd_now():
    return 40587.5 + time.time() / 86400.0        # Unix epoch 的 MJD = 40587.5


def mjd2date(mjd):
    y, m, d, f = erfa.jd2cal(2400000.5, mjd)
    return f"{y:04d}-{m:02d}-{d:02d}"


def log(msg):
    """控制台输出（目录名含零宽字符等, GBK 控制台下自动降级）。"""
    try:
        print(msg, flush=True)
    except UnicodeEncodeError:
        enc = sys.stdout.encoding or "gbk"
        print(msg.encode(enc, errors="replace").decode(enc), flush=True)


# ----------------------------------------------------------------------
# 1b) 快速产品 (finals.daily): 解析 + 与 C04 合并
# ----------------------------------------------------------------------
def parse_rapid(path):
    """解析 USNO finals.daily / finals2000A.daily（空格分隔, 误差列穿插）。

    列序: MJD x xerr y yerr UT1-UTC ut1err LOD loderr [dX dxerr dY dyerr]
    （1980 版末列为 dpsi/deps, 本项目视方向与 c2t06a 均不消费, 忽略。）
    返回 [(mjd, x_as, y_as, dut1_s, lod_s), ...] 按 mjd 升序。
    """
    rows = []
    with open(path, "r", encoding="ascii", errors="replace") as f:
        for line in f:
            t = line.split()
            if len(t) < 9:
                continue
            try:
                mjd = float(t[0])
                x = float(t[1])
                y = float(t[3])
                dut1 = float(t[5])
                lod = float(t[7])
            except ValueError:
                continue
            if not (30000.0 < mjd < 90000.0):
                continue
            rows.append((mjd, x, y, dut1, lod))
    rows.sort(key=lambda r: r[0])
    return rows


def merge_rapid_into_c04(c04_path, rapid_rows):
    """C04 历史 + 快速产品末端 → 写回 C04 同名同格式文件（C 库零改动）。

    只取快速产品中 MJD > C04 末端的行（C04 重叠段为最终值, 保留不动）。
    dX/dY 置 0（视方向不做章动修正、c2t06a 不消费该列, 无影响）。
    返回 (n_appended, mjd_first_rapid, mjd_last_rapid)。
    """
    rows = parse_eop(c04_path)
    last = rows[-1][1]
    add = [r for r in rapid_rows if r[0] > last]
    if not add:
        return 0, None, None
    with open(c04_path, "a", encoding="ascii") as f:
        for mjd, x, y, dut1, lod in add:
            y_, m_, d_, _fd = erfa.jd2cal(2400000.5, mjd)
            f.write(f"{y_:4d} {m_:2d} {d_:2d} {mjd:7.1f} {x:11.6f} {y:11.6f} "
                    f"{dut1:12.7f} {lod:11.7f} "
                    f"{'0.000000':>11s} {'0.000000':>11s} "
                    f"{'0.000000':>11s} {'0.000000':>11s} "
                    f"{'0.0000000':>12s} {'0.0000000':>12s} "
                    f"{'0.000000':>11s} {'0.000000':>11s}\n")
    return len(add), add[0][0], add[-1][0]


def rapid_vs_c04_diff(c04_rows, rapid_rows, n_days=30):
    """重叠段 |快速−最终| 统计（回答“快速产品精度”的实测口径）。"""
    import bisect
    m = [r[1] for r in c04_rows]
    dx, dy, dd = [], [], []
    n0 = max(0, len(m) - n_days)
    for mjd, x, y, dut1, lod in rapid_rows:
        k = bisect.bisect_left(m, mjd)
        if k < len(m) and m[k] == mjd and k >= n0:
            dx.append(abs(x - c04_rows[k][2]))
            dy.append(abs(y - c04_rows[k][3]))
            dd.append(abs(dut1 - c04_rows[k][4]))
    if not dx:
        return None
    return {"n": len(dx), "pm_max_mas": max(max(dx), max(dy)) * 1000.0,
            "dut1_max_ms": max(dd) * 1000.0}


def step_rapid(args):
    """下载/定位快速产品并与 C04 合并。返回报告信息 dict。"""
    log("\n===== [1b] 快速产品 (finals.daily) =====")
    info = {"mode": "none"}
    rapid_rows = None

    src = None
    if args.input_rapid and os.path.exists(args.input_rapid):
        shutil.copyfile(args.input_rapid, RAPID_FILE)
        src = f"本地文件 {args.input_rapid}"
    else:
        tmp = RAPID_FILE + ".tmp"
        for name, url in RAPID_SOURCES:
            try:
                http_download(url, tmp, timeout=180)
                os.replace(tmp, RAPID_FILE)
                src = f"{name} :: {url}"
                break
            except Exception as e:                    # noqa: BLE001
                log(f"  快速源失败: {name}: {e}")
    if src is None:
        log("  快速产品源不可达 → 退化为纯 C04 长期产品 + 外推（当前模式）")
        log("  手动启用: 下载 https://maia.usno.navy.mil/ser7/finals.daily"
            " 放到 eop/ 后加参数 --input-rapid eop/finals.daily")
        return info

    rapid_rows = parse_rapid(RAPID_FILE)
    if len(rapid_rows) < 10:
        log(f"  快速产品解析异常（{len(rapid_rows)} 行）→ 忽略")
        return info
    log(f"  快速产品: {len(rapid_rows)} 行, "
        f"覆盖 MJD {rapid_rows[0][0]:.0f}~{rapid_rows[-1][0]:.0f}")

    c04_rows = parse_eop(EOP_FILE)
    diff = rapid_vs_c04_diff(c04_rows, rapid_rows)
    n_app, m0, m1 = merge_rapid_into_c04(EOP_FILE, rapid_rows)
    info = {"mode": "merged", "source": src, "n_rapid": len(rapid_rows),
            "mjd_rapid_last": rapid_rows[-1][0],
            "date_rapid_last": mjd2date(rapid_rows[-1][0]),
            "n_appended": n_app, "diff": diff}
    if n_app:
        log(f"  已合并 {n_app} 行快速数据到 C04 之后（至 "
            f"{mjd2date(m1)}），滞后降至 {mjd_now() - m1:.1f} 天")
    else:
        log("  快速产品未超出 C04 末端, 无需合并")
    if diff:
        log(f"  重叠段实测 |快速−最终|: 极移 max {diff['pm_max_mas']:.2f} mas, "
            f"DUT1 max {diff['dut1_max_ms']:.3f} ms（{diff['n']} 天样本）")
    return info


# ----------------------------------------------------------------------
# 1) 下载 + 更新
# ----------------------------------------------------------------------
def parse_eop(path):
    rows = []
    with open(path, "r", encoding="ascii", errors="replace") as f:
        for line in f:
            t = line.split()
            if len(t) < 10:
                continue
            try:
                iy = int(t[0])
                vals = [float(x) for x in t[3:10]]
            except ValueError:
                continue
            if not (1962 <= iy <= 2100 and 30000 < vals[0] < 90000):
                continue
            rows.append((iy, *vals))
    if len(rows) < 2:
        raise ValueError(f"解析失败: 仅 {len(rows)} 行")
    return rows


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def http_download(url, dst, timeout=900):
    """流式下载（hpiers 较慢, ~22 KB/s, 3.7 MB 约需 3 分钟）。"""
    import urllib.request
    log(f"  下载 {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "star-catalog-eop-updater/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dst, "wb") as f:
        total = 0
        t0 = time.time()
        while True:
            chunk = resp.read(1 << 16)
            if not chunk:
                break
            f.write(chunk)
            total += len(chunk)
        log(f"  完成: {total} 字节, {time.time() - t0:.1f} s")
    return total


def step_download(args):
    log("\n===== [1/4] 下载 IERS EOP C04 =====")
    if args.input:
        src_url = f"本地文件 {args.input}"
        tmp = EOP_FILE + ".tmp"
        shutil.copyfile(args.input, tmp)
    else:
        tmp = EOP_FILE + ".tmp"
        src_url = None
        for name, url in SOURCES:
            try:
                http_download(url, tmp)
                src_url = f"{name} :: {url}"
                break
            except Exception as e:                       # noqa: BLE001
                log(f"  失败: {name}: {e}")
        if src_url is None:
            log("  所有数据源均失败！保留现有数据继续（如有）。")
            R["download"] = {"ok": False, "source": "全部失败"}
            if not os.path.exists(EOP_FILE):
                return False
            src_url = "保留旧文件（本次下载失败）"
    try:
        rows = parse_eop(tmp)
    except Exception as e:                               # noqa: BLE001
        log(f"  下载内容校验失败: {e}")
        R["download"] = {"ok": False, "source": src_url, "error": str(e)}
        os.remove(tmp)
        return False

    mjd_first, mjd_last = rows[0][1], rows[-1][1]
    lag = mjd_now() - mjd_last
    fresh = lag <= 50.0
    if not fresh:
        log(f"  警告: 数据末端滞后 {lag:.1f} 天 (>50 天)")

    # 原子替换 + 备份
    n_prev = None
    if os.path.exists(EOP_FILE):
        try:
            n_prev = len(parse_eop(EOP_FILE))
        except Exception:                                # noqa: BLE001
            n_prev = None
        shutil.copyfile(EOP_FILE, EOP_FILE + ".prev")
    os.replace(tmp, EOP_FILE)

    last = rows[-1]
    R["download"] = {
        "ok": True, "source": src_url,
        "size": os.path.getsize(EOP_FILE), "sha256": sha256_of(EOP_FILE),
        "n": len(rows), "n_prev": n_prev,
        "mjd_first": mjd_first, "mjd_last": mjd_last,
        "date_first": mjd2date(mjd_first), "date_last": mjd2date(mjd_last),
        "lag_days": lag, "fresh": fresh,
        "latest": {"xp_as": last[2], "yp_as": last[3], "dut1_s": last[4],
                   "lod_s": last[5], "dx_as": last[6], "dy_as": last[7]},
        "time": now_utc().strftime("%Y-%m-%d %H:%M:%S UTC"),
    }
    d = R["download"]
    log(f"  更新完成: {d['n']} 天记录, 覆盖 {d['date_first']} ~ {d['date_last']}"
        f"（滞后 {d['lag_days']:.1f} 天）")
    log(f"  最新 EOP: x={d['latest']['xp_as']:.6f}\" y={d['latest']['yp_as']:.6f}\" "
        f"DUT1={d['latest']['dut1_s']:.7f}s dX={d['latest']['dx_as']:.6f}\" "
        f"dY={d['latest']['dy_as']:.6f}\"")
    return True


# ----------------------------------------------------------------------
# 2) 编译
# ----------------------------------------------------------------------
def find_cmake():
    p = shutil.which("cmake")
    return p


def find_vs():
    """vswhere 动态查找（不硬编码 VS 版本），失败则回退目录遍历。"""
    vswhere = os.path.join(r"C:\Program Files (x86)",
                           "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if os.path.exists(vswhere):
        try:
            out = subprocess.run(
                [vswhere, "-latest", "-products", "*",
                 "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                 "-property", "installationPath"],
                capture_output=True, timeout=30).stdout.decode(
                    "utf-8", errors="replace").strip()
            if out:
                return out
        except Exception:                                # noqa: BLE001
            pass
    for pf in (r"C:\Program Files", r"C:\Program Files (x86)"):
        for ver in (18, 17, 19, 16):
            for ed in ("BuildTools", "Community", "Professional", "Enterprise"):
                p = os.path.join(pf, "Microsoft Visual Studio", str(ver), ed)
                if os.path.isdir(p):
                    return p
    return None


def step_build():
    log("\n===== [2/4] 编译 C 库与测试 =====")
    cmake = find_cmake()
    vs = find_vs()
    cdir = os.path.join(ROOT, "c")
    bdir = os.path.join(cdir, "build")
    log(f"  cmake: {cmake or '未找到'}")
    log(f"  VS  : {vs or '未找到'}")

    if not cmake:
        R["build"] = {"ok": False, "error": "cmake 不在 PATH"}
        return False

    ok = False
    toolsets = ["v145", None]          # 本机 VS2026 BuildTools v145, 其它机器回退默认
    cfg_out = ""
    for ts in toolsets:
        cmd = [cmake, "-S", cdir, "-B", bdir]
        if ts:
            cmd += ["-T", ts]
        r = subprocess.run(cmd, capture_output=True, timeout=600)
        cfg_out = (r.stdout + r.stderr).decode("utf-8", errors="replace")[-800:]
        if r.returncode == 0:
            ok = True
            break
        log(f"  配置失败 (toolset={ts or 'default'}), 尝试下一个")
    if not ok:
        log("  cmake 配置失败:\n" + cfg_out)
        R["build"] = {"ok": False, "error": cfg_out}
        return False

    r = subprocess.run([cmake, "--build", bdir, "--config", "Release"],
                       capture_output=True, timeout=1800)
    out = (r.stdout + r.stderr).decode("utf-8", errors="replace")
    n_warn = out.count("warning C")
    n_err = out.count("error C")
    ok = r.returncode == 0 and n_err == 0
    if not ok:
        log("  构建失败详情(末 500 字符):\n" + out[-500:])
    exes = [f for f in ("test_apparent_catalog.exe", "test_eop.exe")
            if os.path.exists(os.path.join(bdir, "Release", f))]
    log(f"  构建{'成功' if ok else '失败'}: returncode={r.returncode} "
        f"warning={n_warn} error={n_err}, 产物: {', '.join(exes)}")
    R["build"] = {"ok": ok, "warnings": n_warn, "errors": n_err,
                  "targets": exes, "vs": os.path.basename(vs) if vs else "?",
                  "time": now_utc().strftime("%Y-%m-%d %H:%M:%S UTC")}
    return ok


# ----------------------------------------------------------------------
# 3) 测试（C 测试 + Python 交叉验证）
# ----------------------------------------------------------------------
def run_cmd(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=1800,
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    out = (r.stdout + r.stderr).decode("utf-8", errors="replace")
    return r.returncode, out


# ----------------------------------------------------------------------
# 3b) 长期产品外推精度留出回测（不使用 IERS 快速产品）
# ----------------------------------------------------------------------
def backtest_extrapolation(horizons=(1, 5, 10, 20, 31), n_starts=6):
    """用 C04 自身历史做留出回测。

    对每个预测跨度 h: 取 6 个截断起点 t0（相隔 30 天），只用 t0 及之前的数据，
    以 [t0-30, t0] 线性趋势外推极移、以 -LOD 外推 DUT1，与 t0+k (k=1..h) 真值
    比对；同时给出"端点保持（不外推）"对照。折合地固旋转误差 = 矢量合成
    (极移误差 + DUT1×ω 换算, 1 ms ↔ 0.0150″)。
    """
    rows = parse_eop(EOP_FILE)
    m = np.array([r[1] for r in rows])
    xp = np.array([r[2] for r in rows])
    yp = np.array([r[3] for r in rows])
    dut = np.array([r[4] for r in rows])
    lod = np.array([r[5] for r in rows])
    out = []
    for h in horizons:
        e_ext = []      # (dut_ms, pm_as) 外推
        e_hold = []     # (dut_ms, pm_as) 端点保持
        for j in range(n_starts):
            i0 = len(m) - 1 - h - 30 * j
            if i0 - 30 < 0:
                break
            j0 = max(0, i0 - 30)
            x = m[j0:i0 + 1] - m[i0]
            bx, ax_ = np.polyfit(x, xp[j0:i0 + 1], 1)
            by, ay_ = np.polyfit(x, yp[j0:i0 + 1], 1)
            for k in range(1, h + 1):
                ik = i0 + k
                e_ext.append((abs((dut[i0] - lod[i0] * k) - dut[ik]) * 1e3,
                              np.hypot(ax_ + bx * k - xp[ik], ay_ + by * k - yp[ik])))
                e_hold.append((abs(dut[i0] - dut[ik]) * 1e3,
                               np.hypot(xp[i0] - xp[ik], yp[i0] - yp[ik])))
        ee = np.array(e_ext)
        eh = np.array(e_hold)
        rot = np.sqrt(ee[:, 1] ** 2 + (0.015 * ee[:, 0]) ** 2)   # 折合旋转误差, ″
        out.append({
            "h": h,
            "dut_ext_max_ms": float(ee[:, 0].max()),
            "dut_hold_max_ms": float(eh[:, 0].max()),
            "pm_ext_max_as": float(ee[:, 1].max()),
            "pm_hold_max_as": float(eh[:, 1].max()),
            "rot_ext_max_as": float(rot.max()),
            "n": len(e_ext),
        })
    return out


def cross_validate_rotation():
    """C 导出的 ITRS→GCRS 旋转样本 vs Python(erfa.c2t06a) 逐项比对。"""
    sample = os.path.join(ROOT, "c", "eop_rotation_sample.csv")
    rows = []
    with open(sample, newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    dr_max = dv_max = 0.0
    for row in rows:
        jd = float(row["jd_utc"])
        d1 = int(jd)
        d2 = jd - d1
        iy, im, idd, fd = erfa.jd2cal(d1, d2)
        # 与 C 端 ac_utc_to_jd 相同的组装方式
        jd_c = erfa.cal2jd(iy, im, idd)
        jd_utc = jd_c[0] + jd_c[1] + fd
        dat = erfa.dat(iy, im, idd, fd)
        jd_tt = jd_utc + (dat + 32.184) / 86400.0
        # EOP 插值
        mjd = jd_utc - 2400000.5
        rows_e = parse_eop(EOP_FILE)
        m = np.array([r[1] for r in rows_e])
        xp_a = np.array([r[2] for r in rows_e])
        yp_a = np.array([r[3] for r in rows_e])
        dut_a = np.array([r[4] for r in rows_e])
        i = int(np.clip(np.searchsorted(m, mjd), 1, len(m) - 1))
        t = (mjd - m[i - 1]) / (m[i] - m[i - 1])
        xp = xp_a[i - 1] + t * (xp_a[i] - xp_a[i - 1])
        yp = yp_a[i - 1] + t * (yp_a[i] - yp_a[i - 1])
        dut1 = dut_a[i - 1] + t * (dut_a[i] - dut_a[i - 1])
        ut1 = jd_utc + dut1 / 86400.0
        rc2t = erfa.c2t06a(jd_tt, 0.0, ut1, 0.0,
                           xp * erfa.DAS2R, yp * erfa.DAS2R)
        r_i = np.array([7000.0, 1000.0, -1500.0])
        v_i = np.array([-1.1, 6.9, 0.3])
        w = np.array([0.0, 0.0, 7.292115146706979e-5])
        r_g = rc2t.T @ r_i
        v_g = rc2t.T @ (v_i + np.cross(w, r_i))
        dr_max = max(dr_max, max(abs(r_g[k] - float(row[f"rg{'xyz'[k]}"])) for k in range(3)))
        dv_max = max(dv_max, max(abs(v_g[k] - float(row[f"vg{'xyz'[k]}"])) for k in range(3)))
    return len(rows), dr_max, dv_max


def step_test():
    log("\n===== [3/4] 测试 =====")
    cdir = os.path.join(ROOT, "c")
    bdir = os.path.join(cdir, "build", "Release")
    eop_arg = os.path.relpath(EOP_FILE, cdir)
    cat_arg = os.path.relpath(os.path.join(ROOT, "std_catalog_G9p0.bin"), cdir)

    # C: 星表链路回归
    cmd = [os.path.join(bdir, "test_apparent_catalog.exe"), cat_arg,
           "2026", "11", "15", "12.0", "6928.137", "0", "0",
           "0", "-1.0065", "7.5219", "200"]
    rc, out = run_cmd(cmd, cdir)
    reg_ok = rc == 0
    log(f"  test_apparent_catalog: {'PASS' if reg_ok else 'FAIL'}")
    R["t_regression"] = {"ok": reg_ok}
    if not reg_ok:
        log(out[-600:])

    # C: EOP 测试
    have_cat = os.path.exists(os.path.join(ROOT, "std_catalog_G9p0.bin"))
    cmd = [os.path.join(bdir, "test_eop.exe"), eop_arg] + \
          ([cat_arg] if have_cat else [])
    rc, out = run_cmd(cmd, cdir)
    eop_ok = rc == 0
    n_pass = out.count("[PASS]")
    n_fail = out.count("[FAIL]")
    log(f"  test_eop: {'PASS' if eop_ok else 'FAIL'} ({n_pass} 项通过, {n_fail} 项失败)")
    R["t_eop"] = {"ok": eop_ok, "pass": n_pass, "fail": n_fail,
                  "detail": [ln.strip() for ln in out.splitlines()
                             if ln.strip().startswith("[")]}
    if not eop_ok:
        log(out[-600:])

    # Python: 旋转交叉验证
    try:
        n, dr, dv = cross_validate_rotation()
        ok = dr < 1e-3 and dv < 1e-6      # 容差 1 m / 1 mm/s, 实际预期 ~1e-5 km
        log(f"  ITRS→GCRS 旋转 C vs Python: {n} 组, "
            f"max|Δr|={dr:.3e} km, max|Δv|={dv:.3e} km/s -> "
            f"{'PASS' if ok else 'FAIL'}")
        R["t_crossval"] = {"ok": ok, "n": n, "dr_km": dr, "dv_kms": dv}
    except Exception as e:                                # noqa: BLE001
        log(f"  旋转交叉验证异常: {e}")
        R["t_crossval"] = {"ok": False, "error": str(e)}

    # 长期产品外推精度回测
    try:
        bt = backtest_extrapolation()
        log("  长期产品(C04)外推留出回测: 跨度/外推DUT1(max ms)/外推极移(max \")/"
            "折合旋转(max \")")
        for b in bt:
            log(f"    {b['h']:>3} 天: {b['dut_ext_max_ms']:.2f} / "
                f"{b['pm_ext_max_as']:.4f} / {b['rot_ext_max_as']:.4f}"
                f"（端点保持对照 {b['dut_hold_max_ms']:.2f} ms, "
                f"{b['pm_hold_max_as']:.4f}\"）")
        R["t_backtest"] = {"ok": True, "rows": bt}
    except Exception as e:                                # noqa: BLE001
        log(f"  外推回测异常: {e}")
        R["t_backtest"] = {"ok": False, "error": str(e)}

    all_ok = reg_ok and eop_ok and R["t_crossval"].get("ok", False)
    R["test"] = {"ok": all_ok,
                 "time": now_utc().strftime("%Y-%m-%d %H:%M:%S UTC")}
    return all_ok


# ----------------------------------------------------------------------
# 4) 报告
# ----------------------------------------------------------------------
def step_report(path):
    log(f"\n===== [4/4] 生成报告 {os.path.basename(path)} =====")
    d = R.get("download", {})
    b = R.get("build", {})
    t = R.get("test", {})
    lines = []
    A = lines.append
    A("# IERS EOP C04 接入流水线报告")
    A("")
    A(f"- 生成时间：{now_utc().strftime('%Y-%m-%d %H:%M:%S UTC')}")
    A(f"- 流程版本：update_eop.py {__version__} / 库 {__version__}")
    A(f"- 运行主机：{os.environ.get('COMPUTERNAME', '?')} "
      f"({sys.platform})，Python {sys.version.split()[0]}，ERFA {erfa.__version__}")
    A("")
    A("## 一、数据下载与更新")
    if d.get("ok"):
        lag = d["lag_days"]
        A(f"- 数据源：{d['source']}")
        A(f"- 文件：`eop/eopc04_IAU2000.62-now`（{d['size']:,} 字节，"
          f"SHA256 前 16 位 `{d['sha256']}`）")
        A(f"- 记录数：{d['n']:,} 天"
          + (f"（上版 {d['n_prev']:,} 天，新增 {d['n'] - d['n_prev']:,} 天）"
             if d.get("n_prev") else ""))
        A(f"- 覆盖范围：{d['date_first']} ~ {d['date_last']}"
          f"（MJD {d['mjd_first']:.0f}–{d['mjd_last']:.0f}）")
        A(f"- 新鲜度：距今滞后 **{lag:.1f} 天** "
          f"{'✅ 正常（C04 约月更）' if d['fresh'] else '⚠️ 偏旧，检查数据源'}")
        lt = d["latest"]
        A("")
        A("| 最新历元参数 | 数值 |")
        A("|---|---|")
        A(f"| 极移 x | {lt['xp_as']:.6f}″ |")
        A(f"| 极移 y | {lt['yp_as']:.6f}″ |")
        A(f"| UT1−UTC | {lt['dut1_s']:.7f} s |")
        A(f"| LOD | {lt['lod_s']:.7f} s |")
        A(f"| 天极偏差 dX | {lt['dx_as']:.6f}″ |")
        A(f"| 天极偏差 dY | {lt['dy_as']:.6f}″ |")
        if lag > 1.0:
            hold = lag * abs(lt['lod_s']) * 1000.0
            rp = R.get("rapid", {})
            A("")
            if rp.get("mode") == "merged":
                A(f"> **策略**：长期产品 C04 为历史基底 + **快速产品 finals.daily**"
                  f"补充末端（至 {rp.get('date_rapid_last', '?')}），当前有效滞后 "
                  f"{mjd_now() - rp.get('mjd_rapid_last', mjd_now()):.1f} 天；"
                  f"仅超出快速产品末端的部分才外推（DUT1 按 −LOD·Δ、极移按末端"
                  f"30 天线性趋势）。")
            else:
                A(f"> **策略**：当前**快速产品源不可达，退化为纯长期产品 C04 + 外推**"
                  f"（DUT1 按 −LOD·Δ、极移按末端 30 天线性趋势）。数据滞后 "
                  f"{lag:.1f} 天，若仅端点保持 DUT1 误差约 {hold:.1f} ms"
                  f"（↔ {hold * 1e-3 * 7.292e-5 * 206265:.2f}″ 地固旋转）；"
                  f"外推后的实际残差见下方回测表。快速产品典型精度优于外推 1~2 个"
                  f"数量级，网络可达或手动提供 finals.daily（--input-rapid）后自动"
                  f"升级为合并模式。")
    else:
        A(f"- ❌ 下载/更新失败：{d.get('error', d.get('source', '未知'))}")
    rp = R.get("rapid", {})
    if rp.get("mode") == "merged":
        A("")
        A("### 快速产品合并详情")
        A(f"- 来源：{rp['source']}")
        A(f"- 快速序列：{rp['n_rapid']} 行，末端 {rp['date_rapid_last']}")
        A(f"- 本次并入 C04 之后：{rp['n_appended']} 行（重叠段保留 C04 最终值）")
        d = rp.get("diff")
        if d:
            A(f"- 重叠段实测 |快速−最终|（{d['n']} 天）：极移 max "
              f"{d['pm_max_mas']:.2f} mas，DUT1 max {d['dut1_max_ms']:.3f} ms"
              f" —— 即合并模式下末端数据的实际精度")
    A("")
    A("## 二、编译")
    if b.get("ok") is True:
        A(f"- ✅ 成功（Visual Studio：{b.get('vs', '?')}，工具集优先 v145）")
        A(f"- 编译告警 {b.get('warnings', 0)} 处，错误 {b.get('errors', 0)} 处")
        A(f"- 产物：{', '.join(b.get('targets', []))}")
    elif b.get("ok") is False:
        A(f"- ❌ 失败：`{str(b.get('error', b.get('targets', '')))[:400]}`")
    else:
        A("- ⏭️ 本轮跳过（--skip-build）")
    A("")
    A("## 三、测试")
    if t:
        reg = R.get("t_regression", {})
        A(f"- C 星表链路回归（test_apparent_catalog）："
          f"{'✅ PASS' if reg.get('ok') else '❌ FAIL'}")
        te = R.get("t_eop", {})
        if te:
            A(f"- C EOP 测试（test_eop）："
              f"{'✅ PASS' if te.get('ok') else '❌ FAIL'}"
              f"（{te.get('pass', 0)} 项通过 / {te.get('fail', 0)} 项失败）")
            for ln in te.get("detail", []):
                mark = "✅" if "[PASS]" in ln else ("❌" if "[FAIL]" in ln else "·")
                text = ln
                for tag in ("[PASS] ", "[FAIL] ", "[INFO] "):
                    if text.startswith(tag):
                        text = text[len(tag):]
                        break
                A(f"  - {mark} {text}")
        cv = R.get("t_crossval", {})
        if cv:
            if cv.get("ok"):
                A(f"- ITRS→GCRS 旋转 C vs Python(erfa.c2t06a)："
                  f"✅ PASS（{cv['n']} 组样本，max|Δr|={cv['dr_km']:.2e} km，"
                  f"max|Δv|={cv['dv_kms']:.2e} km/s）")
            else:
                A(f"- ITRS→GCRS 旋转交叉验证：❌ FAIL（{cv}）")
        bt = R.get("t_backtest", {})
        if bt.get("ok"):
            A("")
            A("- 长期产品（C04）外推留出回测（6 个截断起点×逐日评估，"
              "外推法 = DUT1 −LOD·Δ + 极移末端 30 天线性趋势）：")
            A("")
            A("| 预测跨度 | DUT1 外推 max (ms) | DUT1 端点保持 max (ms) | "
              "极移外推 max (″) | 极移保持 max (″) | 折合旋转误差 max (″) |")
            A("|---|---|---|---|---|---|")
            for row in bt["rows"]:
                A(f"| {row['h']} 天 | {row['dut_ext_max_ms']:.2f} | "
                  f"{row['dut_hold_max_ms']:.2f} | {row['pm_ext_max_as']:.4f} | "
                  f"{row['pm_hold_max_as']:.4f} | {row['rot_ext_max_as']:.4f} |")
        A("")
        A(f"**测试总结：{'✅ 全部通过' if t.get('ok') else '❌ 存在失败项'}**")
    else:
        A("- ⏭️ 本轮跳过（--skip-test）")
    A("")
    A("## 四、结论")
    overall = (d.get("ok") and (b.get("ok") is not False)
               and (t.get("ok", False) if t else True))
    A(f"- 总体状态：**{'✅ PASS — 数据、构建、测试全部正常' if overall else '❌ FAIL — 见上文分项'}**")
    A("- 本流水线可重复执行：下载(带镜像回退)→原子更新→重编译→全量测试→本报告")
    A("- EOP 用途：`ac_itrs_to_gcrs()` / `ac_apparent_catalog_itrs()`（ECEF 状态"
      "直接输入），DUT1/极移来自 C04 插值")
    A("")
    A("## 更新记录（最新在前）")
    stamp = now_utc().strftime("%Y-%m-%d %H:%M UTC")
    entry = f"- **{stamp}**：{d.get('date_last', '?')} 数据"
    if d.get("ok"):
        entry += f"（{d['n']:,} 天，滞后 {d['lag_days']:.1f} 天）；"
    entry += "编译" + ("✅" if b.get("ok") else "❌") + "，测试"
    entry += "✅" if t.get("ok") else "❌"
    old = []
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            # 仅保留形如 "- **2026-11-15 05:12 UTC**：" 的历史条目行
            old = [ln for ln in f
                   if ln.startswith("- **2") and "-**：" not in ln
                   and ln.count("UTC**：") == 1
                   and not ln.startswith("- **" + stamp)]
    A(entry)
    lines.extend(old[:20])
    A("")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    log(f"报告已写入: {path}（{len(lines)} 行）")
    return True


def main():
    ap = argparse.ArgumentParser(description="IERS EOP 下载/更新/编译/测试/报告")
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--skip-test", action="store_true")
    ap.add_argument("--input", default=None, help="本地 EOP C04 文件(跳过网络下载)")
    ap.add_argument("--input-rapid", default=None,
                    help="本地快速产品文件(finals.daily), 如网络不可达时手动提供")
    ap.add_argument("--report", default=os.path.join(ROOT, "EOP_report.md"))
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    log(f"update_eop {__version__} — IERS EOP C04 流水线")
    os.makedirs(EOP_DIR, exist_ok=True)

    ok = True
    if args.skip_download:
        log("\n===== [1/4] 下载（跳过，使用现有数据）=====")
        if os.path.exists(EOP_FILE):
            rows = parse_eop(EOP_FILE)
            last = rows[-1]
            R["download"] = {
                "ok": True, "source": "现有文件（跳过下载）",
                "size": os.path.getsize(EOP_FILE), "sha256": sha256_of(EOP_FILE),
                "n": len(rows), "n_prev": None,
                "mjd_first": rows[0][1], "mjd_last": last[1],
                "date_first": mjd2date(rows[0][1]), "date_last": mjd2date(last[1]),
                "lag_days": mjd_now() - last[1],
                "fresh": (mjd_now() - last[1]) <= 50.0,
                "latest": {"xp_as": last[2], "yp_as": last[3],
                           "dut1_s": last[4], "lod_s": last[5],
                           "dx_as": last[6], "dy_as": last[7]},
                "time": now_utc().strftime("%Y-%m-%d %H:%M:%S UTC")}
        else:
            ok = False
    else:
        ok = step_download(args)
        if ok:
            R["rapid"] = step_rapid(args)
    if not ok:
        log("数据不可用，终止。")
        step_report(args.report)
        return 1

    if not args.skip_build:
        ok = step_build()
    if not args.skip_test:
        ok = step_test() and ok

    step_report(args.report)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
