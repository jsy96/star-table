# 卫星在轨标准星表与视方向计算（星敏感器姿态确定支持）

从 Gaia DR3 生成的**标准星表**（定长二进制，C/C++ 直接 `fread`）+ 一个
**C/C++ 函数**：输入观测时刻（UTC）与卫星状态矢量（GCRS 位置/速度），
输出该时刻该位置的整表恒星**视方向**（自行传播、周年视差、太阳/地球引力
光线偏折、光行差全部改正）。算法链路与 IAU SOFA `atciq` 完全一致
（逐微角秒校验通过），但不施加 bias-岁差-章动旋转——该共模旋转由姿态
解算吸收，逐星物理改正才需要逐星施加。

## 更新记录（最新在前）

- **2026-11-15 v1.0.261115-1**：初始交付。
  - `build_std_catalog.py`：Gaia DR3（J2016.0）→ 标准星表 `.bin`(C) + `.npz`(Python)，
    含 `pmsafe` 严格历元传播 J2016.0→J2000.0；G≤9 共 177,426 星；
    Barnard(10.39″/yr)、Kapteyn、Groombridge 1830、比邻星(π=768 mas) 等极端星确认在表。
  - C/C++ 库 `c/apparent_catalog.{h,c}`（零第三方依赖，内置 ERFA 子集）。
  - 验证：Python 链路 vs SOFA `apci13+atciq` 最大差 **0.000 μas**；
    C vs Python 逐星交叉验证（2000 星）**全部 <0.0001″ PASS**；
    LEO 示例统计与理论预测吻合（自行+视差中位数 0.1776″ vs 理论 0.1783″，
    最大 112.90″ = Barnard 星 10.87 年自行）。
  - 修复记录：ERFA 视差单位约定（arcsec 非 mas）、ERFA 赤经自行约定
    （dRA/dt，需除 cosδ）——两处约定错误曾在统计中表现为 12%~10³ 倍偏差。

## 文件结构

```
louyifan/
├─ gaia_dr3_13.csv            # 源数据: Gaia DR3 摘要 7,369,632 星 (G≤13, J2016.0)
├─ build_std_catalog.py       # 标准星表生成器 (uv 运行)
├─ apparent_catalog.py        # Python 参考实现 + 自检 + ITRS 转换支持
├─ validate_c_vs_python.py    # C ↔ Python 逐星交叉验证
├─ std_catalog_G9p0.bin       # 生成的标准星表 (C/C++, 80B 头 + 80B/星)
├─ std_catalog_G9p0.npz       # 生成的标准星表 (Python)
└─ c/
   ├─ apparent_catalog.h      # C/C++ 公开 API (extern "C")
   ├─ apparent_catalog.c      # 实现
   ├─ test_apparent_catalog.c # 示例/测试程序 (含统计与 CSV 导出)
   ├─ CMakeLists.txt          # MSVC /utf-8, 工具集 v145
   └─ third_party/erfa/       # ERFA 子集 (LGPL-3.0, 原样收录)
```

## 快速开始

### 1. 生成标准星表（只需一次，或换星等截断时重跑）

```bash
cd louyifan
uv run --no-project --with numpy --with pandas --with pyerfa \
    python -I build_std_catalog.py --input gaia_dr3_13.csv --gmax 9.0
```

### 2. 编译 C 库（MSVC v145，经 vswhere 自动定位）

```bash
cd c
cmake -S . -B build -T v145
cmake --build build --config Release
# 测试 + 示例输出（LEO 550 km 示例状态）
./build/Release/test_apparent_catalog.exe ../std_catalog_G9p0.bin \
    2026 11 15 12.0  6928.137 0 0  0 -1.0065 7.5219  2000
```

### 3. C/C++ 调用

```c
#include "apparent_catalog.h"

ac_catalog_t cat;
if (ac_load_catalog("std_catalog_G9p0.bin", &cat) != AC_OK) { /* ... */ }

double jd = ac_utc_to_jd(2026, 11, 15, 12.0);        /* UTC → JD */
double r[3] = {6928.137, 0, 0};                      /* GCRS, km  */
double v[3] = {0, -1.0065, 7.5219};                  /* GCRS, km/s */

/* 整表视方向（ra/dec 可换 unit_gcrs、brk 按需传 NULL） */
static double ra[N], dec[N];
static ac_breakdown_t brk[N];
ac_apparent_catalog(&cat, jd, r, v, AC_DEFL_SUN_EARTH, ra, dec, NULL, brk);

ac_free_catalog(&cat);
```

性能（Release, 单线程）：整表 177,426 星约 67–132 ms，≈0.4–0.7 μs/星。

### 4. Python 参考实现（含自检与 ITRS 支持）

```bash
uv run --no-project --with numpy --with pyerfa \
    python -I apparent_catalog.py --catalog std_catalog_G9p0.npz
```

```python
import apparent_catalog as ac
cat = ac.load_catalog("std_catalog_G9p0.npz")
out = ac.apparent_catalog(cat, "2026-11-15 12:00:00",
                          r_sat=[6928.137, 0, 0], v_sat=[0, -1.0065, 7.5219],
                          include_breakdown=True)   # 各项改正量
```

## 标准星表格式

`.bin`（小端）＝ 80 字节文件头 + N×80 字节记录：

| 头字段 | 类型 | 说明 |
|---|---|---|
| magic | char[8] | `"STDCAT01"` |
| format_version | int32 | 1 |
| n_stars | int32 | 星数 |
| g_max | double | 生成时 G 截断 |
| jd_epoch_source | double | 2457389.0（Gaia J2016.0） |
| jd_epoch_propagated | double | 2451545.0（J2000.0） |
| reserved | int32[10] | 0 |

记录（与 `ac_star_record_t` 一一对应，无填充）：
`source_id`(i8), `ra_j2000`/`dec_j2000`(rad), `pmra_j2000`(**dRA/dt**, rad/yr,
不含 cosδ——ERFA 约定), `pmdec_j2000`(rad/yr), `parallax_as`(**arcsec**——ERFA 约定),
`rv_kms`(km/s, 缺测=0), `g_mag`, `ra_j2016_deg`/`dec_j2016_deg`(原始 Gaia 值,
误差分析基准)。

`.npz` 字段同名（`ra_j2000`…`dec_j2016_deg`），供 `apparent_catalog.py` 使用。

## 物理模型与精度

链路（每星）：`eraPmpx`（自行+视差，含 Roemer 光行时近似）→ `eraLd`
（太阳+地球引力偏折，dlim 与 SOFA `ldsun` 一致）→ `eraAb`（光行差，
观测者质心系总速度＝地球公转+卫星轨道）。
星历/时间：`eraEpv00`（1900–2100，速度精度 42 mm/s → 光行差误差 ~29 μas）、
`eraDat` 闰秒表（截至 2017-01-01 的 37 s，此后按 37 s 外推）、UTC→TT。
GCRS≈BCRS（差 <0.1 μas）；TT≈TDB（差 <0.1 μas）。

改正量量级（G≤9 星表，2026-11-15 LEO 550 km 示例，全表统计）：

| 项 | 均值 | 中位数 | P99 | 最大 |
|---|---|---|---|---|
| 自行+视差 | 0.3616″ | 0.1776″ | 2.9532″ | 112.90″(Barnard) |
| 引力偏折 | 6.48 mas | 4.43 mas | 36.3 mas | 2168 mas |
| 光行差 | 18.16″ | 19.79″ | 22.17″ | 22.17″ |
| 总改正 | 18.21″ | 19.81″ | 22.77″ | 126.11″ |

## 验证记录（最新在前）

- 2026-11-15：C ↔ Python 交叉验证 2000 星：max|Δra·cosδ|=max|Δdec| <0.0001″，PASS；
  Python 链路 vs SOFA `apci13+atciq`（地心、仅太阳偏折，bpn 逆旋转后比对）：
  max|Δ|=0.000 μas；C 测试程序全表统计与 Gaia 星表直接统计的理论预测一致
  （自行+视差 median 0.1776″ vs 0.1783″，P99 2.953″ vs 2.953″，max 112.90″ vs 112.98″）。

## 注意事项与已知限制

1. **输入坐标系**：C 版接受 GCRS（≈J2000 惯性系）状态矢量；ECEF/ITRS 数据请先转换
   （Python 版 `itrs_to_gcrs` 可用，精确转换需 DUT1）。
2. **时间范围**：1900–2100（`eraEpv00` 适用域）；新增闰秒后需更新 ERFA `dat.c` 表。
3. **不做岁差章动**：输出为 GCRS 视方向；如需 CIRS 视位置自行乘 bpn 矩阵。
4. **负视差星**：`pmsafe` 对视差 ≤0.1 mas 的 ~2000 星做了距离保护（ERFA 标准行为，
   μas 级影响）。
5. **视向速度缺测**置 0（丢失的透视加速度 μas 级）。
6. **许可**：`c/third_party/erfa/` 为 LGPL-3.0（ERFA，源自 IAU SOFA），原样收录并保留
   许可文本；对外发布产品时需遵守 LGPL 条款。
7. 星表二进制为小端格式（x86/ARM-LE 通用；大端平台需转换）。
